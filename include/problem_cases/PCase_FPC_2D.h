//
// Created by bkhara on 4/24/26.
//

#pragma once

#include "ProblemCase.h"

namespace fracture {
    class PCase_FPC_2D : public ProblemCase {

        enum {
            INLET = 1,
            OUTLET = 2,
            BOTTOM_WALL = 3,
            TOP_WALL = 4,
            CYLINDER = 7,
        };

        public:
        PCase_FPC_2D(InputData& idata, FEMachinery& fem, TimeLevelFields& tlf)
            : ProblemCase(idata, fem, tlf) {
            
            // check parameters before proceeding
            {
                mfem::real_t h_min, h_max, kappa_min, kappa_max;
                fem.mesh->GetCharacteristics(h_min, h_max, kappa_min, kappa_max);
                double dt = idata.time_marching.dt;
                if (!Mpi::WorldRank()) {
                    mfem::out << "h_min = " << h_min << ", h_max = " << h_max << "\n";
                    mfem::out << "kappa_min = " << kappa_min << ", kappa_max = " << kappa_max << "\n";
                    mfem::out << "Current dt=" << dt << "\n";
                }
            }
        }

        void ObtainElasticityBoundaryDOFs() override {
            Array<int> bdr_attr;
            Array<int> temp;

            ess_tdof_list_u.SetSize(0);
            bdr_attr.SetSize(fem.mesh->bdr_attributes.Max());

            auto AddEssentialComponent = [&](int attr, int comp) {
                bdr_attr = 0;
                bdr_attr[attr - 1] = 1;

                temp.SetSize(0);
                fem.fespace_primal_u->GetEssentialTrueDofs(bdr_attr, temp, comp);
                ess_tdof_list_u.Append(temp);
            };

            auto AddEssentialVector = [&](int attr) {
                for (int d = 0; d < vdim; d++) {
                    AddEssentialComponent(attr, d);
                }
            };

            // Inlet: prescribed velocity, e.g. u = (u_in(y), 0)
            AddEssentialVector(INLET);

            // Channel walls: no-slip
            AddEssentialVector(BOTTOM_WALL);
            AddEssentialVector(TOP_WALL);

            // Cylinder surface: no-slip
            AddEssentialVector(CYLINDER);

            // Outlet: usually NOT essential for velocity.
            // Pressure/outflow condition is handled elsewhere.

            ess_tdof_list_u.Sort();
            ess_tdof_list_u.Unique();

            // OUTLET Marker
            outlet_marker.SetSize(fem.mesh->bdr_attributes.Max());
            outlet_marker[OUTLET - 1] = 1;
        }

        void ApplyElasticityBC(FractureGridFields& fgf) override {
            Array<int> bdr_attr;
            bdr_attr.SetSize(fem.mesh->bdr_attributes.Max());

            std::vector<Coefficient*> coeffs(vdim, nullptr);

            ConstantCoefficient zero(0.0);

            // ------------------------------------------------------------------
            // Inlet velocity
            //
            // Standard channel inflow:
            //
            //     u_x(y) = 4 Umax (y - ymin)(ymax - y) / (ymax - ymin)^2
            //     u_y(y) = 0
            //
            // Replace ymin/ymax/Umax by your input-data values if available.
            // ------------------------------------------------------------------

            const double ymin = -2.5;
            const double ymax = 2.5; // common benchmark channel height
            const double Umax = 1.0; // example value; choose based on desired Re

            FunctionCoefficient inlet_ux(
                [=](const Vector& x) {
                    const double y = x(1);
                    const double H = ymax - ymin;
                    return 4.0 * Umax * (y - ymin) * (ymax - y) / (H * H);
                });

            {
                bdr_attr = 0;
                bdr_attr[INLET - 1] = 1;

                coeffs[0] = &inlet_ux;
                coeffs[1] = &zero;

                fgf.u.ProjectBdrCoefficient(coeffs.data(), bdr_attr);
            }

            // ------------------------------------------------------------------
            // No-slip walls: u = (0,0)
            // ------------------------------------------------------------------
            {
                bdr_attr = 0;
                bdr_attr[BOTTOM_WALL - 1] = 1;
                bdr_attr[TOP_WALL - 1] = 1;

                coeffs[0] = &zero;
                coeffs[1] = &zero;

                fgf.u.ProjectBdrCoefficient(coeffs.data(), bdr_attr);
            }

            // ------------------------------------------------------------------
            // No-slip cylinder: u = (0,0)
            // ------------------------------------------------------------------
            {
                bdr_attr = 0;
                bdr_attr[CYLINDER - 1] = 1;

                coeffs[0] = &zero;
                coeffs[1] = &zero;

                fgf.u.ProjectBdrCoefficient(coeffs.data(), bdr_attr);
            }

            // ------------------------------------------------------------------
            // Outlet:
            // Do not project velocity here for a standard open/outflow boundary.
            // Typically handled by natural traction / do-nothing condition, or by
            // fixing pressure elsewhere.
            // ------------------------------------------------------------------
        }

        void SetTime(double const t) override {
            ProblemCase::SetTime(t);
        }

        void SetElasticityIC(FractureGridFields& fgf) override {
        }

        void SetPhaseFieldIC(FractureGridFields& fgf) override {
        }

        void AddElasticityRHSIntegrators(ParLinearForm &b) override {
        }

        void PostStep(const double t, const double dt) override {
        }

        void UpdateQuadratureFunctions() override {
        }

        void RegisterParaviewFields(ParaViewDataCollection& pvdc, ParaViewDataCollection& pvdc_q) override {
        }
    };

}