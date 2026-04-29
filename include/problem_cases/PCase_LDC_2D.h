//
// Lid-driven cavity, 2D
//

#pragma once

#include "ProblemCase.h"

namespace nse {

    class PCase_LDC_2D : public ProblemCase {
        enum {
            LEFT  = 1,
            RIGHT = 2,
            BTM   = 3,
            TOP   = 4,
        };

        mfem::Array<int> bdr_attr_u;
        mfem::Array<int> bdr_attr_top;
        mfem::Array<int> bdr_attr_walls;

        int local_corner_vertex = -1;

        // Standard nondimensional lid speed.
        // Replace by an input-data value if you later add one, e.g.
        // idata.ldc2d_inputs.lid_velocity.
        double lid_velocity = 1.0;

    public:
        PCase_LDC_2D(InputData& idata, FEMachinery& fem, TimeLevelFields& tlf)
            : ProblemCase(idata, fem, tlf) {

            MFEM_VERIFY(vdim == 2, "PCase_LDC_2D requires a 2D velocity space.");

            // Print a small mesh/time-step summary, consistent with the other cases.
            {
                mfem::real_t h_min, h_max, kappa_min, kappa_max;
                fem.mesh->GetCharacteristics(h_min, h_max, kappa_min, kappa_max);
                const double dt = idata.time_marching.dt;

                if (!Mpi::WorldRank()) {
                    mfem::out << "[PCase_LDC_2D] h_min = " << h_min
                              << ", h_max = " << h_max << "\n";
                    mfem::out << "[PCase_LDC_2D] kappa_min = " << kappa_min
                              << ", kappa_max = " << kappa_max << "\n";
                    mfem::out << "[PCase_LDC_2D] Current dt = " << dt << "\n";
                    mfem::out << "[PCase_LDC_2D] Lid velocity = "
                              << lid_velocity << "\n";
                }
            }
        }

        void ObtainBoundaryDOFs() override {
            // This case is a closed cavity: no outlet/open boundary.
            has_outlet_bc = false;
            outlet_marker.SetSize(0);

            const int nbdr = fem.mesh->bdr_attributes.Max();

            // Velocity Dirichlet boundary: all cavity walls are essential.
            bdr_attr_u.SetSize(nbdr);
            bdr_attr_u = 0;
            bdr_attr_u[LEFT  - 1] = 1;
            bdr_attr_u[RIGHT - 1] = 1;
            bdr_attr_u[BTM   - 1] = 1;
            bdr_attr_u[TOP   - 1] = 1;

            fem.fespace_primal_u->GetEssentialTrueDofs(bdr_attr_u,
                                                       ess_tdof_list_u);
            ess_tdof_list_u.Sort();
            ess_tdof_list_u.Unique();

            // Separate markers used by ApplyBC().
            bdr_attr_top.SetSize(nbdr);
            bdr_attr_top = 0;
            bdr_attr_top[TOP - 1] = 1;

            bdr_attr_walls.SetSize(nbdr);
            bdr_attr_walls = 0;
            bdr_attr_walls[LEFT  - 1] = 1;
            bdr_attr_walls[RIGHT - 1] = 1;
            bdr_attr_walls[BTM   - 1] = 1;

            // Pressure pin: left-bottom corner, same style as the MMS case.
            // We first infer the local lower-left mesh coordinates. This avoids
            // hard-coding (0,0), while still selecting the geometric lower-left
            // corner on the rank that owns that vertex.
            ess_tdof_list_p.DeleteAll();
            local_corner_vertex = -1;

            double x_min = std::numeric_limits<double>::infinity();
            double y_min = std::numeric_limits<double>::infinity();

            for (int i = 0; i < fem.mesh->GetNV(); i++) {
                const double *coord = fem.mesh->GetVertex(i);
                x_min = std::min(x_min, coord[0]);
                y_min = std::min(y_min, coord[1]);
            }

            const double tol = 1e-9;
            for (int i = 0; i < fem.mesh->GetNV(); i++) {
                const double *coord = fem.mesh->GetVertex(i);
                if (std::abs(coord[0] - x_min) < tol &&
                    std::abs(coord[1] - y_min) < tol) {
                    local_corner_vertex = i;
                    break;
                }
            }

            if (local_corner_vertex != -1) {
                mfem::Array<int> local_dofs;
                fem.fespace_p->GetVertexDofs(local_corner_vertex, local_dofs);

                for (int i = 0; i < local_dofs.Size(); i++) {
                    const int tdof =
                        fem.fespace_p->GetLocalTDofNumber(local_dofs[i]);

                    if (tdof >= 0) {
                        ess_tdof_list_p.Append(tdof);
                    }
                }
            }

            ess_tdof_list_p.Sort();
            ess_tdof_list_p.Unique();

            if (!Mpi::WorldRank()) {
                mfem::out << "[PCase_LDC_2D] ess_tdof_list_u size = "
                          << ess_tdof_list_u.Size() << "\n";
            }

            mfem::out << "[PCase_LDC_2D] Rank = " << Mpi::WorldRank()
                      << " ess_tdof_list_p size = "
                      << ess_tdof_list_p.Size() << "\n";
            ess_tdof_list_p.Print(mfem::out);
        }

        void SetTime(const double t) override {
            ProblemCase::SetTime(t);
        }

        void SetIC(NSEGridFields& fgf) override {
            mfem::Vector zero_vec(vdim);
            zero_vec = 0.0;
            mfem::VectorConstantCoefficient zero_u(zero_vec);
            mfem::ConstantCoefficient zero_p(0.0);

            fgf.u.ProjectCoefficient(zero_u);
            fgf.p.ProjectCoefficient(zero_p);

            // Make sure the initial field also satisfies the lid Dirichlet data.
            ApplyBC(fgf);
        }

        void ApplyBC(NSEGridFields& fgf) override {
            mfem::ConstantCoefficient zero(0.0);
            mfem::ConstantCoefficient lid_u(lid_velocity);

            std::vector<mfem::Coefficient*> coeffs(vdim, nullptr);

            // No-slip on left, right, and bottom walls.
            coeffs[0] = &zero;
            coeffs[1] = &zero;
            fgf.u.ProjectBdrCoefficient(coeffs.data(), bdr_attr_walls);

            // Moving lid on the top boundary: u = (U_lid, 0).
            // This is applied after the side walls. Therefore, if a nodal DOF is
            // shared at a top corner, the lid value wins. This matches the usual
            // discontinuous cavity boundary data convention.
            coeffs[0] = &lid_u;
            coeffs[1] = &zero;
            fgf.u.ProjectBdrCoefficient(coeffs.data(), bdr_attr_top);
        }

        void RegisterParaviewFields(ParaViewDataCollection& pvdc,
                                    ParaViewDataCollection& pvdc_q) override {
            // Nothing special to register for this benchmark.
            // Keep this override so the case has the same shape as the others.
        }
    };
}
