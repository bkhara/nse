//
// Created by Biswajit Khara on 4/18/26.
//
#pragma once

#include "ProblemCase.h"
#include "Integrators.h"

namespace nse {
    class TimeLevelFields;

    class MMS2DStokesExactVelocity : public VectorCoefficient {
    private:
        double omega_;

    public:
        explicit MMS2DStokesExactVelocity(double omega)
            : VectorCoefficient(2), omega_(omega) {
        }

        void Eval(Vector &u, ElementTransformation &T,
                  const IntegrationPoint &ip) override {
            Vector x(2);
            T.Transform(ip, x);
            const double X = x[0];
            const double Y = x[1];

            const double t = GetTime();
            const double ct = std::cos(omega_ * t);

            const double sx = std::sin(M_PI * X);
            const double cx = std::cos(M_PI * X);
            const double sy = std::sin(M_PI * Y);
            const double cy = std::cos(M_PI * Y);

            u.SetSize(2);
            u[0] = M_PI * sx * cy * ct;
            u[1] = -M_PI * cx * sy * ct;
        }
    };

    class MMS2DStokesExactPressure : public Coefficient {
    private:
        double omega_;

    public:
        explicit MMS2DStokesExactPressure(double omega)
            : omega_(omega) {
        }

        double Eval(ElementTransformation &T,
                    const IntegrationPoint &ip) override {
            Vector x(2);
            T.Transform(ip, x);
            const double X = x[0];
            const double Y = x[1];

            const double t = GetTime();
            const double ct = std::cos(omega_ * t);

            const double sx = std::sin(M_PI * X);
            const double sy = std::sin(M_PI * Y);

            return sx * sy * ct;
        }
    };

    class MMS2DStokesRHS : public VectorCoefficient {
    private:
        double nu_;
        double omega_;
        bool disable_convection_;

    public:
        MMS2DStokesRHS(double omega,
                       double nu,
                       bool disable_convection = false)
            : VectorCoefficient(2),
              omega_(omega),
              nu_(nu),
              disable_convection_(disable_convection) {
        }

        void Eval(Vector& f, ElementTransformation& T,
                  const IntegrationPoint& ip) override {
            Vector x(2);
            T.Transform(ip, x);

            const double X = x[0];
            const double Y = x[1];

            const double t = GetTime();

            const double ct = std::cos(omega_ * t);
            const double st = std::sin(omega_ * t);

            const double sx = std::sin(M_PI * X);
            const double cx = std::cos(M_PI * X);
            const double sy = std::sin(M_PI * Y);
            const double cy = std::cos(M_PI * Y);

            f.SetSize(2);

            // ------------------------------------------------------------
            // PDE:
            //
            // u_t + div(u \otimes u) - nu Delta u + grad p = f
            //
            // If disable_convection_ == true, this reduces to the Stokes MMS:
            //
            // u_t - nu Delta u + grad p = f
            // ------------------------------------------------------------

            // Stokes part:
            //
            // f = u_t - nu Delta u + grad p
            f[0] =
                -omega_ * M_PI * sx * cy * st
                + 2.0 * nu_ * M_PI * M_PI * M_PI * sx * cy * ct
                + M_PI * cx * sy * ct;

            f[1] =
                +omega_ * M_PI * cx * sy * st
                - 2.0 * nu_ * M_PI * M_PI * M_PI * cx * sy * ct
                + M_PI * sx * cy * ct;

            // NSE convection part:
            //
            // div(u \otimes u) = (u . grad)u + u div(u)
            //
            // For this MMS, div(u) = 0, so:
            //
            // div(u \otimes u) = (u . grad)u
            //
            // conv_1 = pi^3 sin(pi x) cos(pi x) cos^2(omega t)
            // conv_2 = pi^3 sin(pi y) cos(pi y) cos^2(omega t)
            if (!disable_convection_) {
                const double ct2 = ct * ct;
                const double pi3 = M_PI * M_PI * M_PI;

                f[0] += pi3 * sx * cx * ct2;
                f[1] += pi3 * sy * cy * ct2;
            }
        }
    };

    class PCase_Stokes_MMS_2D : public ProblemCase {
        double omega;
        Array<int> bdr_attr_u;
        Array<int> bdr_attr_p;
        Array<int> traction_bdr_attr;
        ConstantCoefficient one;

        MMS2DStokesExactVelocity exact_velocity;
        MMS2DStokesExactPressure exact_pressure;

        int local_corner_vertex = -1;

        std::string err_file_name;

        enum {
            LEFT = 1,
            RIGHT = 2,
            BTM = 3,
            TOP = 4,
        };

        public:
        PCase_Stokes_MMS_2D(InputData& idata, FEMachinery& fem, TimeLevelFields& tlf)
            : ProblemCase(idata, fem, tlf),
        omega(M_PI),
        exact_velocity(omega), exact_pressure(omega) {

            // instantiate forcing
            forcing_rhs = new MMS2DStokesRHS(omega, idata.flow_properties.nu, idata.flow_properties.disable_convection);

            // check parameters before proceeding
            {
                mfem::real_t h_min, h_max, kappa_min, kappa_max;
                fem.mesh->GetCharacteristics(h_min, h_max, kappa_min, kappa_max);
                double dt = idata.time_marching.dt;
                double VR = 2125; // Rayleigh wave speed in m/s
                double dt_acceptable = h_min / VR;
                if (!Mpi::WorldRank()) {
                    mfem::out << "h_min = " << h_min << ", h_max = " << h_max << "\n";
                    mfem::out << "kappa_min = " << kappa_min << ", kappa_max = " << kappa_max << "\n";
                    mfem::out << "Current dt=" << dt << "\n";
                }
            }

            int nsteps = (idata.time_marching.t_max - idata.time_marching.t0) / idata.time_marching.dt;
            err_file_name = std::string("error_Nt_") + std::to_string(nsteps) + ".txt";
            {
                std::cout << "fname=" << err_file_name << "\n";
                std::ofstream f(err_file_name.c_str());
                f << "t,ex,ey,ep\n";
            }
        }

        void ObtainBoundaryDOFs() override
        {
            // this case does not have outlet
            has_outlet_bc = false;
            outlet_marker.SetSize(0);

            // Dirichlet BCs
            bdr_attr_u.SetSize(fem.mesh->bdr_attributes.Max());
            bdr_attr_u = 0;
            bdr_attr_u[LEFT-1]  = 1;
            bdr_attr_u[RIGHT-1] = 1;
            bdr_attr_u[BTM-1]   = 1;
            bdr_attr_u[TOP-1]   = 1;

            fem.fespace_primal_u->GetEssentialTrueDofs(bdr_attr_u, ess_tdof_list_u);

            ess_tdof_list_p.DeleteAll();

            local_corner_vertex = -1;
            const double x_lft = 0.0;
            const double y_btm = 0.0;

            for (int i = 0; i < fem.mesh->GetNV(); i++)
            {
                double *coord = fem.mesh->GetVertex(i);
                if (std::abs(coord[0] - x_lft) < 1e-9 &&
                    std::abs(coord[1] - y_btm) < 1e-9)
                {
                    local_corner_vertex = i;
                    break;
                }
            }

            if (local_corner_vertex != -1)
            {
                mfem::Array<int> local_dofs;
                fem.fespace_p->GetVertexDofs(local_corner_vertex, local_dofs);

                for (int i = 0; i < local_dofs.Size(); i++)
                {
                    const int tdof = fem.fespace_p->GetLocalTDofNumber(local_dofs[i]);
                    if (tdof >= 0)
                    {
                        ess_tdof_list_p.Append(tdof);
                    }
                }
            }

            ess_tdof_list_p.Sort();
            // remove duplicates if needed in your MFEM version

            std::cout << "[PCASE] Rank = " << Mpi::WorldRank()
                      << " ess_tdof_list_p size = " << ess_tdof_list_p.Size() << "\n";
            ess_tdof_list_p.Print(std::cout);
        }

        void SetTime(double const t) override {
            ProblemCase::SetTime(t);
            ProblemCase::SetTime(t);
            exact_velocity.SetTime(t);
            exact_pressure.SetTime(t);
            forcing_rhs->SetTime(t);
        }

        void SetIC(NSEGridFields& fgf) override {
            exact_velocity.SetTime(0.);
            fgf.u.ProjectCoefficient(exact_velocity);
        }

        void ApplyBC(NSEGridFields &fgf) override {
            fgf.u.ProjectBdrCoefficient(exact_velocity, bdr_attr_u);
            // fgf.p.ProjectBdrCoefficient(exact_pressure, bdr_attr_p);
        }

        void SetAnalyticalSolution(NSEGridFields& fgf) override {
            fgf.u.ProjectCoefficient(exact_velocity);
            fgf.p.ProjectCoefficient(exact_pressure);
        }

        void PostStep(const double t, const double dt) override {
            SetTime(t);
            SetAnalyticalSolution(tlf.exact);

            // Update stored error fields
            tlf.error = tlf.current;
            tlf.error -= tlf.exact;

            ParGridFunction &u_h = tlf.current.u;
            VectorCoefficient &exact = exact_velocity;

            const int vdim = u_h.FESpace()->GetVDim();
            std::vector<real_t> l2_comp(vdim, 0.0);

            for (int c = 0; c < vdim; c++) {
                mfem::Vector e(vdim);
                e = 0.0;
                e(c) = 1.0;

                mfem::VectorConstantCoefficient v_weight(e);

                l2_comp[c] = u_h.ComputeLpError(2.0, exact,
                                                /*weight=*/nullptr,
                                                /*v_weight=*/&v_weight);
            }

            real_t l2_err = 0.0;
            for (int c = 0; c < vdim; c++) {
                l2_err += l2_comp[c] * l2_comp[c];
            }
            l2_err = std::sqrt(l2_err);

            ParGridFunction &p_h = tlf.current.p;
            Coefficient &exact_p = exact_pressure;
            double err_p = p_h.ComputeL2Error(exact_p);

            if (!Mpi::WorldRank()) {
                mfem::out << "||u-u_exact||_L2 = " << l2_err << ", ||p-p_exact||_L2 = " << err_p << "\n";
                mfem::out << "componentwise: "
                        << l2_comp[0] << ", " << l2_comp[1] << "\n";

                std::ofstream f(err_file_name.c_str(), std::ios::app);
                f << t << "," << l2_comp[0] << "," << l2_comp[1] << "," << err_p << "\n";
            }
        }

        void RegisterParaviewFields(ParaViewDataCollection& pvdc, ParaViewDataCollection& pvdc_q) override {
            pvdc.RegisterField("ex_u", &tlf.exact.u);
            pvdc.RegisterField("ex_p", &tlf.exact.p);
            pvdc.RegisterField("err_u", &tlf.error.u);
            pvdc.RegisterField("err_p", &tlf.error.p);
        }
    };
}