//
// Created by khara on 10/23/25.
//

#pragma once

#include <utility>

#include "config_reader/ConfigReader.h"
#include "config_reader/InputReader.h"
#include "Checkpointing.h"

using mycustom::ConfigReader;
using mylib::InputReader;

namespace fracture {
    enum ProblemCaseId : unsigned int {
        APS_2D = 0,
        APS_3D = 1,
        DYN_BRANCHING_2D = 3,
        KALTHOFF_WINKLER = 4,
        BEAM_2D = 100,
        MMS_2D = 101,
        NOTCHED_SHEAR = 102,
        STOKES_MMS_2D = 301,
    };
    enum CouplingFormulation : unsigned int {
        FULLY_COUPLED = 0,
        STAGGERED_COUPLED = 1,
    };
    enum BoundEnforcementType : unsigned int {
        SINGLE_SIDED_CONSTRAINT = 0,
        DOUBLE_SIDED_CONSTRAINT = 1,
    };
    enum class Norm : unsigned int {
        L0 = 0,
        L1 = 1,
        L2 = 2,
        H1 = 3,
        INFTY = 4,
    };
    enum ApplyChangeAt : int {
        NONE = 0,
        BEGINNING_OF_TIMESTEP = 1,
        BEGINNING_OF_STAG_ITERS = 2,
        BEGINNING_OF_PG_ITERS = 3,
        BEGINNING_OF_AL_ITERS = 4,
        END_OF_TIMESTEP = -1,
        END_OF_STAG_ITERS = -2,
        END_OF_PG_ITERS = -3,
        END_OF_AL_ITERS = -4,
    };
    enum PGCprevInterpolationType : unsigned int {
        BY_NODAL_VALUES = 0,
        BY_QUADRATURE_VALUES = 1,
    };
    enum class MacaulayRegType : unsigned int {
        NONE = 0, // exact piecewise (nonsmooth)
        SQRT = 1, // 0.5*(x - sqrt(x^2+eps^2))  (C^1 everywhere)
        SOFTPLUS = 2, // -eps*log(1+exp(-x/eps))     (C^∞)
    };
    // struct IEMethod { // Inequality Enforcement Method Name
    //     static constexpr std::string_view PENALTY = "penalty";
    //     static constexpr std::string_view PG = "pg";
    // };
    // struct FESpaceNames {
    //     static constexpr std::string_view H1 = "h1";
    //     static constexpr std::string_view L2 = "l2";
    // };
    namespace FESpaceNames {
        inline constexpr std::string_view H1 = "h1";
        inline constexpr std::string_view L2 = "l2";
    };
    namespace IEMethod { // Inequality Enforcement Method Name
        // NONE uses HISTFUNC integrators but does not apply the history functional
        inline constexpr std::string_view NO_ENFORCEMENT = "none"; // no enforcement of inequality
        inline constexpr std::string_view PENALTY = "pen"; // penalty
        inline constexpr std::string_view PG = "pg"; // proximal Galerkin
        inline constexpr std::string_view HISTFUNC = "hf"; // history functional
        inline constexpr std::string_view AL = "al"; // augmented Lagrangian
    }
    namespace PetscSolverPrefix {
        inline constexpr std::string_view ELASTICITY = "el_";
        inline constexpr std::string_view PHASEFIELD = "pf_";
    }
    namespace TimeDependenceType {
        inline constexpr std::string_view QUASISTATIC = "qs";
        inline constexpr std::string_view DYNAMIC = "dy";

    }
    namespace SurfaceEnergyRegularizationType {
        inline constexpr std::string_view AT1 = "at1";
        inline constexpr std::string_view AT2 = "at2";
        inline constexpr std::string_view WU = "wu";
    }

    static const char *NormToString(Norm n) {
        switch (n) {
        case Norm::L0: return "L0";
        case Norm::L1: return "L1";
        case Norm::L2: return "L2";
        case Norm::H1: return "H1";
        case Norm::INFTY: return "INFTY";
        }
        return "UNKNOWN";
    }

    static bool CheckValidNorm(unsigned int a) {
        bool is_valid = false;
        switch (static_cast<Norm>(a)) {
        case Norm::L0:
        case Norm::L1:
        case Norm::L2:
        case Norm::H1:
        case Norm::INFTY:
            is_valid = true;
        }
        return is_valid;
    }

    // structs below
    struct FractureData {
        double density = 1;
        double a1 = 1.0;
        double a2 = 0.25;
        double a3 = 1.0;

        double l = 0.02;
        double G = 1.;
        double Gc = 1.;
        double cw = 2.;
        double epsilon = 1e-5;

        double lambda = 1.;
        double mu = 1.;

        bool use_dirichlet_energy = false;

        bool use_tc_split = false; // tension-compression split
        bool use_voldev_2d_trace = false;
        bool use_spectral_split = false;

        bool at1_pf_fix_corner_node = false; // the corner node must be set up in the PCase class

        std::string at_reg_type = std::string(SurfaceEnergyRegularizationType::AT2);

        bool use_at1() const {
            return at_reg_type == SurfaceEnergyRegularizationType::AT1;
        }
        bool use_at2() const {
            return at_reg_type == SurfaceEnergyRegularizationType::AT2;
        }
        bool use_wu() const {
            return at_reg_type == SurfaceEnergyRegularizationType::WU;
        }

        void AdjustParameters() {
            if (use_at1()) {
                cw = 4./3.;
            } else if (use_at2()) {
                cw = 1.0;
            } else if (use_wu()) {
                cw = M_PI / 2.;
            }
        }

        void Print() const {
            if (mfem::Mpi::WorldRank() == 0) {
                std::cout << "FractureData: " << std::endl;
                std::cout << "a1: " << a1 << std::endl;
                std::cout << "a2: " << a2 << std::endl;
                std::cout << "a3: " << a3 << std::endl;
                std::cout << "l: " << l << std::endl;
                std::cout << "G: " << G << std::endl;
                std::cout << "Gc: " << Gc << std::endl;
                std::cout << "cw: " << cw << std::endl;
                std::cout << "epsilon: " << epsilon << std::endl;
                std::cout << "lambda: " << lambda << std::endl;
                std::cout << "mu: " << mu << std::endl;
            }
        }

        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading FractureData\n"; }
            reader.ReadValue("fracture_data.density", density);
            reader.ReadValue("fracture_data.a1", a1);
            reader.ReadValue("fracture_data.a2", a2);
            reader.ReadValue("fracture_data.a3", a3);
            reader.ReadValue("fracture_data.l", l);
            reader.ReadValue("fracture_data.lambda", lambda);
            reader.ReadValue("fracture_data.mu", mu);
            reader.ReadValue("fracture_data.Gc", Gc);
            reader.ReadValue("fracture_data.epsilon", epsilon);
            reader.ReadValue("fracture_data.use_dirichlet_energy", use_dirichlet_energy);
            reader.ReadValue("fracture_data.use_tc_split", use_tc_split);
            reader.ReadValue("fracture_data.use_spectral_split", use_spectral_split);
            reader.ReadValue("fracture_data.use_voldev_2d_trace", use_voldev_2d_trace);
            reader.ReadValue("fracture_data.at_reg_type", at_reg_type);
            reader.ReadValue("fracture_data.at1_pf_fix_corner_node", at1_pf_fix_corner_node);

            // now adjust the parameters based on the config inputs
            AdjustParameters();
        }
    };
    struct MeshConfig {
        bool use_external_mesh = false;
        std::string mesh_file = "aps_tri_unstr.msh";
        // below parameters are for when use_external_mesh = false
        int nsd = 2;
        std::vector<double> lengths = {1., 1.};
        std::vector<int> nels = {8, 8};
        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading MeshInfo\n"; }
            reader.ReadValue("mesh_config.use_external_mesh", use_external_mesh);
            if (use_external_mesh) {
                bool mesh_file_specified = reader.ReadValue("mesh_config.mesh_file", mesh_file);
                if (!mesh_file_specified) {
                    MFEM_ABORT("Mesh file not specified in the config file.\n")
                }
            } else {
                reader.ReadValue("mesh_config.nsd", nsd);
                reader.ReadArray("mesh_config.lengths", lengths);
                reader.ReadArray("mesh_config.nels", nels);
            }
        }
    };
    struct MethodConfig {
        CouplingFormulation formulation = STAGGERED_COUPLED;
        std::string ie_method = std::string(IEMethod::PG);
        bool no_enforcement = false; // experimental, and only works with HF
        bool use_two_sided_constraint = false;
        bool use_snes = true;

        BoundEnforcementType bound_type = SINGLE_SIDED_CONSTRAINT;

        bool use_ie_none() const {
            return ie_method == IEMethod::NO_ENFORCEMENT;
        }
        bool use_ie_penalty() const {
            return ie_method == IEMethod::PENALTY;
        }
        bool use_ie_pg() const {
            return ie_method == IEMethod::PG;
        }
        bool use_ie_hf() const {
            return ie_method == IEMethod::HISTFUNC;
        }
        bool use_ie_al() const {
            return ie_method == IEMethod::AL;
        }
        bool remove_ie() const {
            return no_enforcement;
        }
        bool use_double_bound() const {
            return use_two_sided_constraint;
        }

        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading MethodInfo\n"; }
            reader.ReadValue("method_config.ie_method", ie_method);
            if (use_ie_none()) {
                if (!mfem::Mpi::WorldRank()) {
                    std::cout << "ie_method=" << ie_method << " selected. This will internally set:\n";
                    std::cout << "   * ie_method = " << IEMethod::HISTFUNC << "\n";
                    std::cout << "   * no_enforcement = " << true << "\n";
                }
                ie_method = IEMethod::HISTFUNC;
                no_enforcement = true;
            }
            // reader.ReadValue("method_config.no_enforcement", no_enforcement);
            reader.ReadValue("method_config.use_snes", use_snes);
            reader.ReadValue("method_config.use_two_sided_constraint", use_two_sided_constraint);
            bound_type = SINGLE_SIDED_CONSTRAINT;
            if (use_two_sided_constraint) {
                bound_type = DOUBLE_SIDED_CONSTRAINT;
            }
        }
    };
    struct PGConfig {
        int max_it = 200;
        double stol = 1e-4; // solution tolerance (on the phase variable)
        double psi_stol = 1e-2; // solution tolerance (on the latent variable)
        double alpha_0 = 1.;
        double alpha_min = 0.5;
        double alpha_max = 1024.0;
        Norm norm_type = Norm::INFTY;
        std::string norm_type_string;

        bool add_coercive_mod = true;
        double coercive_mod_coeff = 1e-3;

        bool use_history_function = false;
        bool use_nodal_quadrature = false;
        int cprev_interpolation_type = BY_QUADRATURE_VALUES;
        bool cap_psi_at_qp = false;

        bool use_bddc_solver = false;
        bool use_block_nlf = false;

        bool use_custom_snes_monitor = false;

        int bound_preserving_strategy = 0; // 0: "by grad_rstar", 1: "by pull back"

        bool use_cprev_nodal_values() const {
            return not(use_nodal_quadrature) and (cprev_interpolation_type == BY_NODAL_VALUES);
        }
        bool use_cprev_qfunc_values() const {
            return not(use_nodal_quadrature) and (cprev_interpolation_type == BY_QUADRATURE_VALUES);
        }

        bool preserve_bounds_by_grad_rstar() const {
            return bound_preserving_strategy == 0;
        }
        bool preserve_bounds_by_pull_back() const {
            return bound_preserving_strategy == 1;
        }

        bool terminate_using_l2_norm() const    { return norm_type == Norm::L2; }
        bool terminate_using_linf_norm() const  { return norm_type == Norm::INFTY; }

        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading PGConfig\n"; }
            reader.ReadValue("pg_config.max_it", max_it);
            reader.ReadValue("pg_config.stol", stol);
            reader.ReadValue("pg_config.psi_stol", psi_stol);
            reader.ReadValue("pg_config.alpha_0", alpha_0);
            reader.ReadValue("pg_config.alpha_min", alpha_min);
            reader.ReadValue("pg_config.alpha_max", alpha_max);
            reader.ReadValue("pg_config.add_coercive_mod", add_coercive_mod);
            reader.ReadValue("pg_config.coercive_mod_coeff", coercive_mod_coeff);
            reader.ReadValue("pg_config.use_history_function", use_history_function);
            reader.ReadValue("pg_config.use_nodal_quadrature", use_nodal_quadrature);
            if (use_nodal_quadrature) {
                cprev_interpolation_type = BY_NODAL_VALUES; // only type available
            } else {
                reader.ReadValue("pg_config.cprev_interpolation_type", cprev_interpolation_type);
            }
            reader.ReadValue("pg_config.cap_psi_at_qp", cap_psi_at_qp);
            reader.ReadValue("pg_config.use_block_nlf", use_block_nlf);
            reader.ReadValue("pg_config.use_custom_snes_monitor", use_custom_snes_monitor);
            reader.ReadValue("pg_config.bound_preserving_strategy", bound_preserving_strategy);
            {
                unsigned int tmp = static_cast<unsigned int>(Norm::INFTY);
                reader.ReadValue("pg_config.norm_type", tmp);

                bool is_valid = CheckValidNorm(tmp);
                MFEM_VERIFY(is_valid, "Invalid pg_config.norm_type = " << tmp);

                norm_type = static_cast<Norm>(tmp);
                norm_type_string = NormToString(norm_type);
            }
        }
    };
    struct HFConfig {
        bool use_snes_vi = true;

        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading HFConfig\n"; }
            reader.ReadValue("hf_config.use_snes_vi", use_snes_vi);
        }
    };
    struct ALConfig {
        int max_it = 200;
        double stol = 1e-4; // solution tolerance (on the phase variable)
        double rho_initial = 1e3;
        double rho_min = 1e2;
        double rho_max = 1e5;

        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading ALConfig\n"; }
            reader.ReadValue("al_config.max_it", max_it);
            reader.ReadValue("al_config.stol", stol);
            reader.ReadValue("al_config.rho_initial", rho_initial);
            reader.ReadValue("al_config.rho_min", rho_min);
            reader.ReadValue("al_config.rho_max", rho_max);
        }
    };
    struct PenaltyConfig {
        double coefficient = 0.5;
        double tol_ir = 1e-2; // see page 8 of https://www.sciencedirect.com/science/article/pii/S0045782520305387

        bool use_nl_form = true;
        MacaulayRegType macaulay_reg_type = MacaulayRegType::NONE;
        double macaulay_reg_eps = 1e-12;

        double calc_gamma(double Gc, double l) const {
            return (Gc / (l)) * (1. / (tol_ir * tol_ir) - 1.);
        }

        double calc_gamma(const FractureData &fdata) const {
            return calc_gamma(fdata.Gc, fdata.l);
        }

        void ReadFromFile(InputReader &conf) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading PenaltyConfig\n"; }
            conf.ReadValue("penalty_config.coefficient", coefficient);
            conf.ReadValue("penalty_config.tol_ir", tol_ir);
            conf.ReadValue("penalty_config.use_nl_form", use_nl_form);
            unsigned int temp_reg_type;
            if (conf.ReadValue("penalty_config.macaulay_reg_type",temp_reg_type)) {
                if (temp_reg_type == 0 or temp_reg_type == 1 or temp_reg_type == 2) {
                    macaulay_reg_type = static_cast<MacaulayRegType>(temp_reg_type);
                } else {
                    if (Mpi::Root()) {
                        mfem::err << "macaulay_reg_type=" << temp_reg_type << " not available. Defaulting to NONE.\n";
                    }
                }
                conf.ReadValue("penalty_config.macaulay_reg_eps", macaulay_reg_eps);
            }
        }
    };
    struct StaggeredIterationConfig {
        Norm norm_type = Norm::INFTY;
        std::string norm_type_string;
        int max_it = 1000;
        double etol = 1e-3; // solution tolerance (on the phase variable)

        bool use_l0_norm() const    { return norm_type == Norm::L0; }
        bool use_l1_norm() const    { return norm_type == Norm::L1; }
        bool use_l2_norm() const    { return norm_type == Norm::L2; }
        bool use_h1_norm() const    { return norm_type == Norm::H1; }
        bool use_linf_norm() const  { return norm_type == Norm::INFTY; }

        void ReadFromFile(InputReader &conf) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading StaggeredIterationConfig\n"; }
            conf.ReadValue("stag_config.max_it", max_it);
            conf.ReadValue("stag_config.etol", etol);
            {
                unsigned int tmp = static_cast<unsigned int>(Norm::INFTY);
                conf.ReadValue("stag_config.norm_type", tmp);

                bool is_valid = CheckValidNorm(tmp);
                MFEM_VERIFY(is_valid, "Invalid stag_config.norm_type = " << tmp);

                norm_type = static_cast<Norm>(tmp);
                norm_type_string = NormToString(norm_type);
            }
        }
    };
    struct RefinementConfig {
        int num_initial_refine = 0;
        double reftol = 1e-3;
        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading RefinementConfig\n"; }
        }
    };
    struct TimeMarchingConfig {
        struct GenAlpha {
            bool use_explicit = false;
            bool use_newmark = false;
            bool use_hht = false;
            bool use_wbz = false;
            double rho_inf = 1.0;
            double am = 1.; // a_{n+am} = (1-am) * a_{n} + am * a_{n+1}
            double af = 1.; // t_{n+af} = (1-af) * t_{n} + af * t_{n+1}
            double gamma = 0.5;
            double beta = 0.25;
            void CalcParams() {
                if (use_explicit) {
                    am = (2. - rho_inf) / (1. + rho_inf);
                    beta = (5. - 3. * rho_inf) / (1. + rho_inf) / (1. + rho_inf) / (2. - rho_inf);
                    gamma = 0.5 + am;
                    af = 0.;
                } else {
                    am = (2. - rho_inf) / (1. + rho_inf);
                    af = 1. / (1. + rho_inf);
                    if (use_hht) {
                        am = 1;
                        af = 2. * rho_inf / (1. + rho_inf);
                    } else if (use_wbz) {
                        am = 2. / (1. + rho_inf);
                        af = 1.;
                    } else if (use_newmark) {
                        am = 1.;
                        af = 1.;
                    }
                    gamma = 0.5 + am - af;
                    beta = 0.25 * std::pow(1 + am - af, 2.);
                }
            }
            void PrintGenAlphaParams() const {
                if (mfem::Mpi::WorldRank()) { return; }
                mfem::out << "Printing gen-alpha parameters\n";
                mfem::out << "am=" << am << ", af=" << af << ", beta=" << beta << ", gamma=" << gamma << std::endl;
            }
        };
        GenAlpha gen_alpha;
        std::string time_dependence_type = std::string(TimeDependenceType::QUASISTATIC);
        bool compute_initial_acceleration = false;
        double t0 = 0.;
        double t_max = 1.;
        double dt = 0.1;
        double dt_min = 1e-10;
        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading TimeMarchingConfig\n"; }
            reader.ReadValue("time_marching.t0", t0);
            reader.ReadValue("time_marching.t_max", t_max);
            reader.ReadValue("time_marching.dt", dt);
            reader.ReadValue("time_marching.dt_min", dt_min);
            reader.ReadValue("time_marching.gen_alpha.use_explicit", gen_alpha.use_explicit);
            reader.ReadValue("time_marching.gen_alpha.use_newmark", gen_alpha.use_newmark);
            reader.ReadValue("time_marching.gen_alpha.use_hht", gen_alpha.use_hht);
            reader.ReadValue("time_marching.gen_alpha.use_wbz", gen_alpha.use_wbz);
            reader.ReadValue("time_marching.gen_alpha.rho_inf", gen_alpha.rho_inf);
            reader.ReadValue("time_marching.time_dependence_type", time_dependence_type);
            reader.ReadValue("time_marching.compute_initial_acceleration", compute_initial_acceleration);
            int selected =
                    static_cast<int>(gen_alpha.use_newmark) +
                    static_cast<int>(gen_alpha.use_hht) +
                    static_cast<int>(gen_alpha.use_wbz);
            if (selected > 1)
            {
                MFEM_ABORT("Exactly one of {use_newmark, use_hht, use_wbz} must be true. "
                           "Got selected=" << selected);
            }
            gen_alpha.CalcParams();
            gen_alpha.PrintGenAlphaParams();
        }
        bool is_quasi_static() const {
            return time_dependence_type == TimeDependenceType::QUASISTATIC;
        }
        bool is_dynamic() const {
            return time_dependence_type == TimeDependenceType::DYNAMIC;
        }
    };
    struct FESpaceConfig {
        std::string space = std::string(FESpaceNames::H1);
        int order = 1;
    };

    struct Experimental {
        // for experimental changes that are subject to deprecation later
        bool pf_use_jac_assy_2 = false;
        bool skip_elasticity_solve = false;
        bool elasticity_use_linear_solver = false;
        bool skip_phasefield_solve = false;
        bool elasticity_use_acceleration_based_solve = false;
        bool use_el_nlform = true;
        bool use_pf_nlform = true;
        bool calc_c_violation = true;
        double c_violation_tol = 1e-12;
        bool phase_pull_back = false;
        bool pg_phase_pull_back = false;
        bool check_matrix_symmetry_el = false;
        bool check_matrix_symmetry_pf = false;
        bool check_for_failed_linear_solve = false; //TODO part of debugging
        bool print_element_matrix_pf = false;  //TODO part of debugging. CAUTION: this is algorithmically set, not read from config file
        //TODO make this graceful once the dust settles (e.g., move it our of experiments)
        // 0: reset at the beginning of time step
        // 1: reset at the beginning of staggered iters
        int reset_psi_strategy = ApplyChangeAt::NONE;
        bool reset_psi_at_time_step() const {
            return reset_psi_strategy == ApplyChangeAt::BEGINNING_OF_TIMESTEP;
        }
        bool reset_psi_at_stag_step() const {
            return reset_psi_strategy == ApplyChangeAt::BEGINNING_OF_STAG_ITERS;
        }
        void ReadFromFile(mylib::InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading Experimental inputs\n"; }
            reader.ReadValue("experimental.pf_use_jac_assy_2", pf_use_jac_assy_2);
            reader.ReadValue("experimental.skip_elasticity_solve", skip_elasticity_solve);
            reader.ReadValue("experimental.skip_phasefield_solve", skip_phasefield_solve);
            reader.ReadValue("experimental.elasticity_use_linear_solver", elasticity_use_linear_solver);
            reader.ReadValue("experimental.elasticity_use_acceleration_based_solve", elasticity_use_acceleration_based_solve);
            reader.ReadValue("experimental.use_el_nlform", use_el_nlform);
            reader.ReadValue("experimental.use_pf_nlform", use_pf_nlform);
            reader.ReadValue("experimental.calc_c_violation", calc_c_violation);
            reader.ReadValue("experimental.c_violation_tol", c_violation_tol);
            reader.ReadValue("experimental.check_matrix_symmetry_el", check_matrix_symmetry_el);
            reader.ReadValue("experimental.check_matrix_symmetry_pf", check_matrix_symmetry_pf);
            reader.ReadValue("experimental.reset_psi_strategy", reset_psi_strategy);
            reader.ReadValue("experimental.check_for_failed_linear_solve", check_for_failed_linear_solve);
            reader.ReadValue("experimental.phase_pull_back", phase_pull_back);
            reader.ReadValue("experimental.pg_phase_pull_back", pg_phase_pull_back);
        }
    };

    struct ProblemCaseConfig {
        ProblemCaseId pcase_id = APS_2D;

        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading ProblemCase inputs\n"; }
            unsigned int temp = APS_2D;
            reader.ReadValue("pcase", temp);
            pcase_id = static_cast<ProblemCaseId>(temp);
        }
    };

    struct NotchedShearInputs {
        bool fix_uy_on_sides = true;
        bool fix_uy_at_notch_bdr = false;
        bool use_pre_crack_ic = false;
        bool use_pre_crack_histfunc = false;
        double pre_crack_thickness = 1.875e-3;
        double pre_crack_maxval = 1.;
        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading NotchedShearInputs\n"; }
            reader.ReadValue("nshear_inputs.fix_uy_on_sides", fix_uy_on_sides);
            reader.ReadValue("nshear_inputs.fix_uy_at_notch_bdr", fix_uy_at_notch_bdr);
            reader.ReadValue("nshear_inputs.use_pre_crack_ic", use_pre_crack_ic);
            reader.ReadValue("nshear_inputs.use_pre_crack_histfunc", use_pre_crack_histfunc);
            reader.ReadValue("nshear_inputs.pre_crack_thickness", pre_crack_thickness);
            reader.ReadValue("nshear_inputs.pre_crack_maxval", pre_crack_maxval);
        }
    };

    struct DynamicBranchingInputs {
        double traction_mag = 1e6;
        bool use_pre_crack_histfunc = false;
        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading DynamicBranchingInputs\n"; }
            reader.ReadValue("dy_branching.use_pre_crack_histfunc", use_pre_crack_histfunc);
            reader.ReadValue("dy_branching.sigma_top", traction_mag);
        }
    };

    struct KalthoffWinklerInputs {
        double t0 = 1e-6;
        double v0 = 16.5;
        bool symmetry_bc_fix_uy = false;
        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading KalthoffWinklerInputs\n"; }
            reader.ReadValue("kw_inputs.t0", t0);
            reader.ReadValue("kw_inputs.v0", v0);
            reader.ReadValue("kw_inputs.symmetry_bc_fix_uy", symmetry_bc_fix_uy);
        }
    };

    struct StokesMMS2DInputs {
        double rho = 1.0;
        double Re = 1.0;
        double mu = 1.0;
        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading StokesMMS2DInputs\n"; }
            reader.ReadValue("stokesmms_inputs.Re", Re);
            mu = 1./Re;
        }
    };

    struct StreamingInputs {
        bool if_stream = false;
        std::string host = std::string("localhost");
        int socket = 19916;
        int frequency = 2;
        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading StreamingInputs\n"; }
            reader.ReadValue("streaming_inputs.if_stream", if_stream);
            reader.ReadValue("streaming_inputs.host", host);
            reader.ReadValue("streaming_inputs.socket", socket);
            reader.ReadValue("streaming_inputs.frequency", frequency);
        }
    };

    struct RunContext {
        // THIS STRUCT IS ONLY FOR MONITORING / PRINTING PURPOSES
        // NOT FOR USING IN CALCULATIONS
        // struct for keeping track of running variables
        // such as time step, staggered iterations, PG iterations
        // and other related variables
        int time_step = -1;
        int stag_it = -1;
        int pg_it = -1;
        int al_it = -1;
        double time = 0.;
        double alpha_cur = 0.;

        double u_norm_l2 = 0.;
        double stag_dc_l2;
        double stag_dc_inf;
        double stag_dc_h1;
        double stag_dc_h1_semi;
        double pg_dc_l2;
        double pg_dc_inf;
        double pg_dpsi_l2;
        double pg_dpsi_inf;
        double al_dc_l2;
        double al_dc_inf;
        double al_dlam_l2;
        double al_dlam_inf;
    };

    struct CheckpointingInputs {
        // incoming
        bool restart_from_chkpt = false;
        std::string chkpt_root = "./checkpoint";

        // these will be calculated automatically
        char in_slot = 'A';
        std::string chkpt_dir_in;

        // outgoing
        int chkpt_write_frequency = 200;

        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading CheckpointingInputs\n"; }
            reader.ReadValue("checkpointing_inputs.restart_from_chkpt", restart_from_chkpt);
            reader.ReadValue("checkpointing_inputs.chkpt_root", chkpt_root);
            reader.ReadValue("checkpointing_inputs.chkpt_dir_in", chkpt_dir_in);
            reader.ReadValue("checkpointing_inputs.chkpt_write_frequency", chkpt_write_frequency);

            if (restart_from_chkpt) {
                if (chkpt_dir_in.empty()) {
                    in_slot = Checkpointing::ChooseValidSlot(chkpt_root, MPI_COMM_WORLD);
                    chkpt_dir_in = chkpt_root + "/" + std::string(1, in_slot);
                    if (!Mpi::WorldRank()) {
                        std::cout << "Restarting from checkpoint!\n";
                        std::cout << "      * Determined checkpoint slot = " << in_slot << "!\n";
                        std::cout << "      * Determined checkpoint path = " << chkpt_dir_in << "!\n";
                    }
                } else {
                    if (!Mpi::WorldRank()) {
                        std::cout << "Restarting from checkpoint!\n";
                        std::cout << "      * Given path = " << chkpt_dir_in << "!\n";
                    }
                }
            }
        }
    };

    class InputData {
        mycustom::ConfigReader conf; // not in use currently
        libconfig::Config config;
        mylib::InputReader reader;

    public:
        RunContext run_ctx;

        FESpaceConfig fes_config_primal;
        FESpaceConfig fes_config_latent;
        int file_write_freq = 2;
        int el_vdim = 1;

        CheckpointingInputs checkpointing_inputs;
        StreamingInputs streaming_inputs;
        MeshConfig mesh_config;
        TimeMarchingConfig time_marching;
        MethodConfig method_config;
        PGConfig pg_config;
        HFConfig hf_config;
        ALConfig al_config;
        PenaltyConfig penalty_config;
        StaggeredIterationConfig stag_config;
        RefinementConfig ref_config;
        FractureData fracture_data;
        Experimental experimental;
        ProblemCaseConfig pcase_config;
        NotchedShearInputs nshear_inputs;
        DynamicBranchingInputs dynamic_branching_inputs;
        KalthoffWinklerInputs kw_inputs;

        StokesMMS2DInputs stokes_mms2d_inputs;

        InputData() : conf(true), reader("", &config) {
        }

        void ReadFromFile(const std::string& cfg_filename = std::string("config.txt"), const std::string& petsc_filename = std::string("frac.petsc")) {
            conf.load(cfg_filename);
            reader.ReadConfigFile(cfg_filename);

            reader.ReadValue("el_vdim", el_vdim);
            reader.ReadValue("fe_primal.space", fes_config_primal.space);
            reader.ReadValue("fe_primal.order", fes_config_primal.order);
            reader.ReadValue("fe_latent.space", fes_config_latent.space);
            reader.ReadValue("fe_latent.order", fes_config_latent.order);
            reader.ReadValue("file_write_freq", file_write_freq);
            checkpointing_inputs.ReadFromFile(reader);
            streaming_inputs.ReadFromFile(reader);
            mesh_config.ReadFromFile(reader);
            time_marching.ReadFromFile(reader);
            pcase_config.ReadFromFile(reader);
            fracture_data.ReadFromFile(reader);
            method_config.ReadFromFile(reader);
            pg_config.ReadFromFile(reader);
            {
                if (petsc_filename == std::string("fract.petsc.bddc")) {
                    pg_config.use_bddc_solver = true;
                }
            }
            hf_config.ReadFromFile(reader);
            al_config.ReadFromFile(reader);
            penalty_config.ReadFromFile(reader);
            stag_config.ReadFromFile(reader);
            ref_config.ReadFromFile(reader);
            experimental.ReadFromFile(reader);
            // reader.PrintConfig("dump.cfg");
            nshear_inputs.ReadFromFile(reader);
            dynamic_branching_inputs.ReadFromFile(reader);
            kw_inputs.ReadFromFile(reader);
            stokes_mms2d_inputs.ReadFromFile(reader);
        }
    };
}
