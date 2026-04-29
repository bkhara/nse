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

namespace nse {
    enum ProblemCaseId : unsigned int {
        APS_2D = 0,
        APS_3D = 1,
        DYN_BRANCHING_2D = 3,
        KALTHOFF_WINKLER = 4,
        BEAM_2D = 100,
        MMS_2D = 101,
        NOTCHED_SHEAR = 102,
        STOKES_MMS_2D = 301,
        FPC_2D = 302,
    };
    enum CouplingFormulation : unsigned int {
        FULLY_COUPLED = 0, // both stabilized and unstabilized
        UNCOUPLED = 1, // projection method
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
    namespace NSStabilizationMethod { // Inequality Enforcement Method Name
        // NONE uses HISTFUNC integrators but does not apply the history functional
        inline constexpr std::string_view UNSTABILIZED = "none"; // plain NSE without any modification
        inline constexpr std::string_view VMS_STABILIZED = "vms"; // vms for velocity (and / or pressure)
        inline constexpr std::string_view SUPG_STABILIZED = "supg"; // supg
        inline constexpr std::string_view PSPG_STABILIZED = "pspg"; // pspg
        inline constexpr std::string_view SUPG_PSPG_STABILIZED = "sups"; // sups
    }
    namespace PetscSolverPrefix {
        inline constexpr std::string_view NS_COUPLED = "nse_";
        inline constexpr std::string_view NS_UNCOUPLED_MOMENTUM = "mom_";
        inline constexpr std::string_view NS_UNCOUPLED_PPE = "ppe_";
        inline constexpr std::string_view NS_UNCOUPLED_VUE = "vue_";
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
    namespace TimeMarchingScheme {
        inline constexpr std::string_view BDF2 = "bdf2";
        inline constexpr std::string_view CN = "cn";
    }

    enum class ProjectionScheme {
        ChorinFirstOrder,
        IncPressureBDF2
    };

    struct ProjectionConfig {
        ProjectionScheme scheme = ProjectionScheme::IncPressureBDF2;

        bool IsBDF2() const {
            return scheme == ProjectionScheme::IncPressureBDF2;
        }

        bool IsClassicChorin() const {
            return scheme == ProjectionScheme::ChorinFirstOrder;
        }
    };

    struct VMSConfig {
        double Ci = 36.0;
    };

    struct SUPSConfig {
        bool use_supg = true;
        bool use_pspg = true;
        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading SUPSConfig\n"; }
            reader.ReadValue("sups_config.use_supg", use_supg);
            reader.ReadValue("sups_config.use_pspg", use_pspg);
        }
    };

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
        CouplingFormulation coupling_form = FULLY_COUPLED;
        std::string stab_scheme = std::string(NSStabilizationMethod::SUPG_PSPG_STABILIZED);

        bool is_coupled() const {
            return coupling_form == FULLY_COUPLED;
        }
        bool is_uncoupled() const {
            return coupling_form == UNCOUPLED;
        }

        bool use_stab_none() const {
            return stab_scheme == NSStabilizationMethod::UNSTABILIZED;
        }
        bool use_stab_vms() const {
            return stab_scheme == NSStabilizationMethod::VMS_STABILIZED;
        }
        bool use_stab_sups() const {
            return stab_scheme == NSStabilizationMethod::SUPG_PSPG_STABILIZED;
        }

        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading MethodInfo\n"; }
            reader.ReadValue("method_config.stab_scheme", stab_scheme);
            {
                auto tmp = static_cast<int>(coupling_form);
                reader.ReadValue("method_config.coupling_form", tmp);
                coupling_form = static_cast<CouplingFormulation>(tmp);
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

        std::string time_dependence_type = std::string(TimeDependenceType::QUASISTATIC);
        bool compute_initial_acceleration = false;
        std::string marching_scheme = std::string(TimeMarchingScheme::BDF2);
        double t0 = 0.;
        double t_max = 1.;
        double dt = 0.1;
        double dt_min = 1e-10;
        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading TimeMarchingConfig\n"; }
            reader.ReadValue("time_marching.marching_scheme", marching_scheme);
            reader.ReadValue("time_marching.t0", t0);
            reader.ReadValue("time_marching.t_max", t_max);
            reader.ReadValue("time_marching.dt", dt);
            reader.ReadValue("time_marching.dt_min", dt_min);
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

    struct FlowPropertiesInputs {
        double Re = 1.0;
        double nu = 1.0;
        bool disable_convection = false;
        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading FlowProperties\n"; }
            reader.ReadValue("Re", Re);
            reader.ReadValue("disable_convection", disable_convection);
            nu = 1. / Re;
        }
    };

    struct MMS2DInputs {
        int mms_type = 0;
        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading MMS2DInputs\n"; }
            reader.ReadValue("mms2d_inputs.mms_type", mms_type);
        }
    };

    struct FPC2DInputs {
        struct InitialNoise {
            double num_initial_steps = 20;
            double eps = 1e-2;
        };
        double H = 5.0; // channel width
        bool cylinder_flip_sign_for_force_calc = true;
        InitialNoise initial_noise;
        void ReadFromFile(InputReader &reader) {
            if (!mfem::Mpi::WorldRank()) { mfem::out << "Reading FPC2DInputs\n"; }
            reader.ReadValue("fpc2d_inputs.H", H);
            reader.ReadValue("fpc2d_inputs.cylinder_flip_sign_for_force_calc", cylinder_flip_sign_for_force_calc);
            reader.ReadValue("fpc2d_inputs.initial_noise.num_initial_steps", initial_noise.num_initial_steps);
            reader.ReadValue("fpc2d_inputs.initial_noise.eps", initial_noise.eps);
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
        int vel_vdim = 2;

        ProblemCaseConfig pcase_config;
        FlowPropertiesInputs flow_properties;

        CheckpointingInputs checkpointing_inputs;
        StreamingInputs streaming_inputs;
        MeshConfig mesh_config;
        TimeMarchingConfig time_marching;

        MethodConfig method_config;
        ProjectionConfig projection_config;
        VMSConfig vms_config;
        SUPSConfig sups_config;

        RefinementConfig ref_config;
        StaggeredIterationConfig stag_config;
        Experimental experimental;

        MMS2DInputs mms2d_inputs;
        FPC2DInputs fpc2d_inputs;

        InputData() : conf(true), reader("", &config) {
        }

        void ReadFromFile(const std::string& cfg_filename = std::string("config.txt"), const std::string& petsc_filename = std::string("frac.petsc")) {
            conf.load(cfg_filename);
            reader.ReadConfigFile(cfg_filename);

            reader.ReadValue("vel_vdim", vel_vdim);
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
            flow_properties.ReadFromFile(reader);
            mms2d_inputs.ReadFromFile(reader);
            fpc2d_inputs.ReadFromFile(reader);
            method_config.ReadFromFile(reader);
            sups_config.ReadFromFile(reader);
        }
    };
}
