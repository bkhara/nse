//
// Created by bkhara on 4/24/26.
//

#pragma once

#include "ProblemCase.h"
#include "Utils.h"

namespace nse {
    struct DragLiftForces {
        double drag = 0.0;
        double lift = 0.0;
    };

    class PCase_FPC_2D : public ProblemCase {

        enum {
            INLET = 1,
            OUTLET = 2,
            BOTTOM_WALL = 3,
            TOP_WALL = 4,
            CYLINDER = 7,
        };

        std::string forcefilename = "forces.txt";

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

            if (Mpi::Root()) {
                std::ofstream file(forcefilename.c_str());
                file << "t,Cd,Cl\n";
                file.close();
            }
        }

        void ObtainBoundaryDOFs() override {
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
            has_outlet_bc = true;
            outlet_marker.SetSize(fem.mesh->bdr_attributes.Max());
            outlet_marker[OUTLET - 1] = 1;
        }

        void ApplyBC(NSEGridFields& fgf) override {
            const double t = tlf.GetTime();
            const double dt = tlf.GetTimeStep();

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

            const double channel_width = idata.fpc2d_inputs.H;
            const double ymin = -channel_width / 2.0;
            const double ymax = channel_width / 2.0; // common benchmark channel height
            const double Umax = 1.0; // example value; choose based on desired Re

            FunctionCoefficient inlet_ux(
                [=](const Vector& x) {
                    const double y = x(1);
                    const double H = ymax - ymin;
                    return 4.0 * Umax * (y - ymin) * (ymax - y) / (H * H);
                });

            // Perturb only during startup.
            const double perturb_time = idata.fpc2d_inputs.initial_noise.num_initial_steps * dt;
            const double eps = idata.fpc2d_inputs.initial_noise.eps;
            FunctionCoefficient inlet_uy_perturbed(
                [=](const Vector &x) {
                    if (t >= perturb_time) {
                        return 0.0;
                    }

                    const double y = x(1);
                    const double yc = 0.5 * (ymin + ymax);
                    const double halfH = 0.5 * (ymax - ymin);

                    const double eta = (y - yc) / halfH; // eta in [-1, 1]

                    // Antisymmetric about channel centerline.
                    // Zero at centerline and at top/bottom walls.
                    return eps * Umax * eta * (1.0 - eta * eta);
                });

            {
                bdr_attr = 0;
                bdr_attr[INLET - 1] = 1;

                coeffs[0] = &inlet_ux;
                coeffs[1] = &inlet_uy_perturbed; // = 0 for t > perturb_time

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

        void PostStep(const double t, const double dt) override {
            mfem::Array<int> cylinder_marker(fem.mesh->bdr_attributes.Max());
            cylinder_marker = 0;
            cylinder_marker[CYLINDER - 1] = 1;

            DragLiftForces dl = ComputeDragLiftOnBoundary(cylinder_marker,
                                                    2,
                                                    fem.ordering,
                                                    -1,
                                                    idata.fpc2d_inputs.cylinder_flip_sign_for_force_calc);

            double rho = 1.0;
            double Uref = 1.0; // the max at mid-channel
            double D = 1.0; // diameter of the cylinder
            double reference_force = 0.5 * rho * Uref * Uref * D;
            double Cd = dl.drag / reference_force;
            double Cl = dl.lift / reference_force;
            if (Mpi::Root()) {
                std::ofstream file(forcefilename.c_str(), std::ios_base::app);
                file << t << "," << Cd << "," << Cl << "\n";
                file.close();
            }
            if (Mpi::Root()) {
                std::ofstream file(forcefilename.c_str(), std::ios_base::app);
                std::cout << "t=" << t << ", Cd=" << Cd << ", Cl=" << Cl << "\n";
                file.close();
            }
        }

        DragLiftForces ComputeDragLiftOnBoundary(const mfem::Array<int>& bdr_marker,
                                           const int vdim,
                                           const mfem::Ordering::Type ordering,
                                           const int quad_order = -1,
                                           const bool flip_sign = false) {
            MFEM_VERIFY(vdim == 2, "This drag/lift helper is written for 2D.");

            mfem::ParMesh &pmesh = *fem.mesh;
            const mfem::ParFiniteElementSpace* ufes = fem.fespace_primal_u;
            const mfem::ParFiniteElementSpace* pfes = fem.fespace_p;

            MFEM_VERIFY(ufes != nullptr, "u_gf must have a valid ParFiniteElementSpace.");
            MFEM_VERIFY(pfes != nullptr, "p_gf must have a valid ParFiniteElementSpace.");

            double local_force[2] = {0.0, 0.0};

            mfem::Array<int> vdofs_u;
            mfem::Array<int> dofs_p;

            mfem::Vector u_el;
            mfem::Vector p_el;

            mfem::Vector Nu;
            mfem::Vector Np;

            mfem::DenseMatrix dNu;
            mfem::DenseMatrix grad_u(vdim, vdim);

            mfem::Vector nor(vdim);
            mfem::Vector traction(vdim);

            const double nu = idata.flow_properties.nu;

            for (int be = 0; be < pmesh.GetNBE(); be++) {
                const int attr = pmesh.GetBdrAttribute(be);

                if (attr <= 0 || attr > bdr_marker.Size()) {
                    continue;
                }

                if (bdr_marker[attr - 1] == 0) {
                    continue;
                }

                mfem::FaceElementTransformations* Tr =
                    pmesh.GetBdrFaceTransformations(be);

                MFEM_VERIFY(Tr != nullptr, "Null boundary FaceElementTransformations.");
                MFEM_VERIFY(Tr->Elem1 != nullptr, "Null Elem1 transformation.");
                MFEM_VERIFY(Tr->Face != nullptr, "Null Face transformation.");

                const int elem = Tr->Elem1No;

                const mfem::FiniteElement* el_u = ufes->GetFE(elem);
                const mfem::FiniteElement* el_p = pfes->GetFE(elem);

                const int dof_u = el_u->GetDof();
                const int dof_p = el_p->GetDof();
                const int dim = Tr->Elem1->GetSpaceDim();

                MFEM_VERIFY(dim == 2, "This drag/lift helper is written for 2D.");
                MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");

                ufes->GetElementVDofs(elem, vdofs_u);
                pfes->GetElementDofs(elem, dofs_p);

                tlf.current.u.GetSubVector(vdofs_u, u_el);
                tlf.current.p.GetSubVector(dofs_p, p_el);

                Nu.SetSize(dof_u);
                Np.SetSize(dof_p);
                dNu.SetSize(dof_u, dim);

                const int order =
                    (quad_order > 0)
                        ? quad_order
                        : 2 * std::max(el_u->GetOrder(), el_p->GetOrder());

                const mfem::IntegrationRule& ir =
                    mfem::IntRules.Get(Tr->GetGeometryType(), order);

                for (int iq = 0; iq < ir.GetNPoints(); iq++) {
                    const mfem::IntegrationPoint& ip_face = ir.IntPoint(iq);

                    Tr->Face->SetIntPoint(&ip_face);

                    mfem::IntegrationPoint ip_el;
                    Tr->Loc1.Transform(ip_face, ip_el);
                    Tr->Elem1->SetIntPoint(&ip_el);

                    el_u->CalcShape(ip_el, Nu);
                    el_u->CalcPhysDShape(*Tr->Elem1, dNu);
                    el_p->CalcShape(ip_el, Np);

                    double p = 0.0;
                    for (int a = 0; a < dof_p; a++) {
                        p += p_el(a) * Np(a);
                    }

                    grad_u = 0.0;

                    for (int a = 0; a < dof_u; a++) {
                        for (int c = 0; c < vdim; c++) {
                            const int ia = VDofIndex(dof_u, vdim, a, c, ordering);
                            const double ua_c = u_el(ia);

                            for (int j = 0; j < dim; j++) {
                                grad_u(c, j) += ua_c * dNu(a, j);
                            }
                        }
                    }

                    // Scaled outward normal:
                    //
                    //     nor = n |J_face|
                    //
                    // Therefore ds is already contained in nor.
                    // Only multiply by the quadrature weight.
                    mfem::CalcOrtho(Tr->Face->Jacobian(), nor);

                    traction = 0.0;

                    // traction = sigma n
                    //
                    // sigma = -p I + nu (grad u + grad u^T)
                    //
                    // Because nor is scaled, this is really sigma * nor.
                    for (int i = 0; i < dim; i++) {
                        traction(i) += -p * nor(i);

                        for (int j = 0; j < dim; j++) {
                            traction(i) +=
                                nu * (grad_u(i, j) + grad_u(j, i)) * nor(j);
                        }
                    }

                    const double w = ip_face.weight;

                    local_force[0] += traction(0) * w;
                    local_force[1] += traction(1) * w;
                }
            }

            double global_force[2] = {0.0, 0.0};

            MPI_Allreduce(local_force, global_force, 2, MPI_DOUBLE, MPI_SUM,
                          pmesh.GetComm());

            DragLiftForces dl;

            dl.drag = global_force[0];
            dl.lift = global_force[1];

            if (flip_sign) {
                dl.drag *= -1.0;
                dl.lift *= -1.0;
            }

            return dl;
        }
    };

}
