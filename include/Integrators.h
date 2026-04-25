//
// Created by khara on 11/17/25.
//

#pragma once

#include "mfem.hpp"
#include "InputData.h"
#include "TimeLevelFields.h"
#include "Utils.h"

using namespace mfem;
using namespace nse;

namespace nse {
    class NSEIntegratorBase : public mfem::BlockNonlinearFormIntegrator {
    protected:
        const InputData &idata;
        const TimeLevelFields &tlf;

        mfem::VectorCoefficient *f_coeff;

        const int vdim;
        const mfem::Ordering::Type ordering;

        static void EvalVectorAtIP(const mfem::Vector &elfun,
                                   const int ndof,
                                   const int vdim,
                                   const mfem::Ordering::Type ordering,
                                   const mfem::Vector &N,
                                   mfem::Vector &u) {
            u.SetSize(vdim);
            u = 0.0;

            for (int c = 0; c < vdim; c++) {
                for (int a = 0; a < ndof; a++) {
                    u(c) += elfun(VDofIndex(ndof, vdim, a, c, ordering)) * N(a);
                }
            }
        }

        static void EvalVectorGradAtIP(const mfem::Vector &elfun,
                                       const int ndof,
                                       const int vdim,
                                       const mfem::Ordering::Type ordering,
                                       const mfem::DenseMatrix &dN,
                                       mfem::DenseMatrix &grad_u) {
            const int dim = dN.Width();
            grad_u.SetSize(vdim, dim);
            grad_u = 0.0;

            for (int c = 0; c < vdim; c++) {
                for (int a = 0; a < ndof; a++) {
                    const double ua = elfun(VDofIndex(ndof, vdim, a, c, ordering));
                    for (int j = 0; j < dim; j++) {
                        grad_u(c, j) += ua * dN(a, j);
                    }
                }
            }
        }

        static double EvalScalarAtIP(const mfem::Vector &elfun,
                                     const mfem::Vector &N) {
            double val = 0.0;
            for (int a = 0; a < N.Size(); a++) {
                val += elfun(a) * N(a);
            }
            return val;
        }

        static double Divergence(const mfem::DenseMatrix &grad_u) {
            const int vdim = grad_u.Height();
            double div_u = 0.0;
            for (int c = 0; c < vdim; c++) { div_u += grad_u(c, c); }
            return div_u;
        }

        // const mfem::IntegrationRule *GetRule(const mfem::FiniteElement &el_u,
        //                                      const mfem::FiniteElement &el_p) const {
        //     if (IntRule) { return IntRule; }
        //
        //     return &mfem::IntRules.Get(el_u.GetGeomType(),
        //                                std::max(2 * el_u.GetOrder(),
        //                                         2 * el_p.GetOrder()));
        // }

    public:
        NSEIntegratorBase(const InputData &idata,
                             const TimeLevelFields &tlf,
                             const int vdim,
                             const mfem::Ordering::Type ordering,
                             mfem::VectorCoefficient *f_coeff = nullptr)
            : idata(idata),
              tlf(tlf),
              f_coeff(f_coeff),
              vdim(vdim),
              ordering(ordering) {
        }
    };

    class NSEBlockIntegBDF2 : public NSEIntegratorBase {
    public:
        NSEBlockIntegBDF2(const InputData& idata,
                             const TimeLevelFields& tlf,
                             const int vdim,
                             const mfem::Ordering::Type ordering,
                             mfem::VectorCoefficient* f_coeff = nullptr)
            : NSEIntegratorBase(idata, tlf, vdim, ordering, f_coeff) {
        }

        void AssembleElementVector(const mfem::Array<const mfem::FiniteElement*>& el,
                                   mfem::ElementTransformation& T,
                                   const mfem::Array<const mfem::Vector*>& elfun,
                                   const mfem::Array<mfem::Vector*>& elvec) override {
            const mfem::FiniteElement& el_u = *el[0];
            const mfem::FiniteElement& el_p = *el[1];

            const int ndof_u = el_u.GetDof();
            const int ndof_p = el_p.GetDof();
            const int dim = T.GetSpaceDim();

            MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");

            elvec[0]->SetSize(vdim * ndof_u);
            *elvec[0] = 0.0;

            elvec[1]->SetSize(ndof_p);
            *elvec[1] = 0.0;

            const mfem::Vector& eu = *elfun[0];
            const mfem::Vector& ep = *elfun[1];

            mfem::Vector Nu(ndof_u), Np(ndof_p);
            mfem::DenseMatrix dNu(ndof_u, dim);

            mfem::Vector u_n(vdim), u_nm1(vdim), u_nm2(vdim), f_n(vdim);
            mfem::DenseMatrix grad_u_n(vdim, dim);

            const int e = T.ElementNo;
            const mfem::IntegrationRule* ir = &tlf.femach.qspace->GetIntRule(e);

            const double ctime = tlf.GetTime(); // t_{n+1}
            const double dt = tlf.GetTimeStep();

            const double nu = idata.flow_properties.nu;

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint& ip = ir->IntPoint(iq);
                T.SetIntPoint(&ip);

                el_u.CalcShape(ip, Nu);
                el_u.CalcPhysDShape(T, dNu);
                el_p.CalcShape(ip, Np);

                const double wdet = ip.weight * T.Weight();

                EvalVectorAtIP(eu, ndof_u, vdim, ordering, Nu, u_n);
                EvalVectorGradAtIP(eu, ndof_u, vdim, ordering, dNu, grad_u_n);

                const double p_n = EvalScalarAtIP(ep, Np);

                tlf.prev_1.u.GetVectorValue(T, ip, u_nm1);
                tlf.prev_2.u.GetVectorValue(T, ip, u_nm2);

                f_n.SetSize(vdim);
                f_n = 0.0;
                if (f_coeff) {
                    f_coeff->SetTime(ctime);
                    f_coeff->Eval(f_n, T, ip);
                }

                double div_u_n = 0.0;
                for (int c = 0; c < vdim; c++) {
                    div_u_n += grad_u_n(c, c);
                }

                // ------------------------------------------------------------
                // Momentum residual block
                //
                // Conservative/divergence-form NSE:
                //
                // (3u^{n+1} - 4u^n + u^{n-1})/(2dt)
                // + div(u^{n+1} \otimes u^{n+1})
                // - nu Delta u^{n+1}
                // + grad p^{n+1}
                // - f^{n+1} = 0
                //
                // Weak form:
                //
                // (3/(2dt) u^{n+1}, v)
                // - ((4u^n - u^{n-1})/(2dt), v)
                // - (u^{n+1} \otimes u^{n+1}, grad v)
                // + nu (grad u^{n+1}, grad v)
                // - (p^{n+1}, div v)
                // - (f^{n+1}, v)
                //
                // Boundary flux from convection is not included here:
                //
                // + <(u \otimes u)n, v>_{\partial\Omega}
                //
                // This term vanishes on strongly imposed no-slip/inflow Dirichlet
                // boundaries since the test function is zero there. For open/outflow
                // boundaries, add a separate boundary integrator if needed.
                // ------------------------------------------------------------
                for (int a = 0; a < ndof_u; a++) {
                    for (int c = 0; c < vdim; c++) {
                        const int ia = VDofIndex(ndof_u, vdim, a, c, ordering);

                        // BDF2 mass contribution from u^{n+1}
                        (*elvec[0])(ia) +=
                            (3.0 / (2.0 * dt)) * u_n(c) * Nu(a) * wdet;

                        // BDF2 history contribution
                        (*elvec[0])(ia) +=
                            -(4.0 * u_nm1(c) - u_nm2(c)) / (2.0 * dt) * Nu(a) * wdet;

                        // Conservative convection contribution:
                        //
                        // - (u \otimes u, grad v)
                        //
                        // Component form:
                        //
                        // - sum_j u_c u_j d_j v_c
                        //
                        if (not idata.flow_properties.disable_convection) {
                            for (int j = 0; j < dim; j++) {
                                (*elvec[0])(ia) +=
                                    -u_n(c) * u_n(j) * dNu(a, j) * wdet;
                            }
                        }

                        // Diffusion contribution
                        for (int j = 0; j < dim; j++) {
                            (*elvec[0])(ia) +=
                                nu * grad_u_n(c, j) * dNu(a, j) * wdet;
                        }

                        // Pressure contribution
                        (*elvec[0])(ia) += -p_n * dNu(a, c) * wdet;

                        // Forcing contribution
                        (*elvec[0])(ia) += -f_n(c) * Nu(a) * wdet;
                    }
                }

                // ------------------------------------------------------------
                // Continuity residual block
                //
                // (q, div u^{n+1})
                // ------------------------------------------------------------
                for (int a = 0; a < ndof_p; a++) {
                    (*elvec[1])(a) += Np(a) * div_u_n * wdet;
                }
            }
        }

        void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement*>& el,
                                 mfem::ElementTransformation& T,
                                 const mfem::Array<const mfem::Vector*>& elfun,
                                 const mfem::Array2D<mfem::DenseMatrix*>& elmat) override {
            const mfem::FiniteElement& el_u = *el[0];
            const mfem::FiniteElement& el_p = *el[1];

            const int ndof_u = el_u.GetDof();
            const int ndof_p = el_p.GetDof();
            const int dim = T.GetSpaceDim();

            MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");

            elmat(0, 0)->SetSize(vdim * ndof_u, vdim * ndof_u);
            *elmat(0, 0) = 0.0;

            elmat(0, 1)->SetSize(vdim * ndof_u, ndof_p);
            *elmat(0, 1) = 0.0;

            elmat(1, 0)->SetSize(ndof_p, vdim * ndof_u);
            *elmat(1, 0) = 0.0;

            elmat(1, 1)->SetSize(ndof_p, ndof_p);
            *elmat(1, 1) = 0.0;

            mfem::DenseMatrix& Auu = *elmat(0, 0);
            mfem::DenseMatrix& Aup = *elmat(0, 1);
            mfem::DenseMatrix& Apu = *elmat(1, 0);

            const mfem::Vector& eu = *elfun[0];

            mfem::Vector Nu(ndof_u), Np(ndof_p);
            mfem::DenseMatrix dNu(ndof_u, dim);

            mfem::Vector u_n(vdim);
            mfem::DenseMatrix grad_u_np1(vdim, dim);

            const int e = T.ElementNo;
            const mfem::IntegrationRule* ir = &tlf.femach.qspace->GetIntRule(e);

            const double dt = tlf.GetTimeStep();
            const double nu = idata.flow_properties.nu;

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint& ip = ir->IntPoint(iq);
                T.SetIntPoint(&ip);

                el_u.CalcShape(ip, Nu);
                el_u.CalcPhysDShape(T, dNu);
                el_p.CalcShape(ip, Np);

                const double wdet = ip.weight * T.Weight();

                EvalVectorAtIP(eu, ndof_u, vdim, ordering, Nu, u_n);
                EvalVectorGradAtIP(eu, ndof_u, vdim, ordering, dNu, grad_u_np1);

                // ------------------------------------------------------------
                // Auu block
                //
                // Linear contribution:
                //
                // d/du [
                //   3/(2dt) (u, v) + nu (grad u, grad v)
                // ]
                //
                // Conservative convection contribution:
                //
                // C(u; v) = - (u \otimes u, grad v)
                //
                // Newton derivative:
                //
                // dC(u)[du; v]
                // =
                // - (du \otimes u, grad v)
                // - (u \otimes du, grad v)
                //
                // Component form:
                //
                // C_c = - sum_j u_c u_j d_j v_c
                //
                // dC_c / d(du_k)
                // =
                // - delta_{ck} du_k sum_j u_j d_j v_c
                // - u_c du_k d_k v_c
                // ------------------------------------------------------------
                for (int a = 0; a < ndof_u; a++) {
                    for (int b = 0; b < ndof_u; b++) {
                        const double mass_ab =
                            (3.0 / (2.0 * dt)) * Nu(a) * Nu(b) * wdet;

                        double diff_ab = 0.0;
                        for (int j = 0; j < dim; j++) {
                            diff_ab += nu * dNu(a, j) * dNu(b, j) * wdet;
                        }

                        // Mass + diffusion: diagonal in velocity components
                        for (int c = 0; c < vdim; c++) {
                            const int ia = VDofIndex(ndof_u, vdim, a, c, ordering);
                            const int ib = VDofIndex(ndof_u, vdim, b, c, ordering);

                            Auu(ia, ib) += mass_ab + diff_ab;
                        }

                        if (not idata.flow_properties.disable_convection) {
                            // Compute u . grad(phi_a), where phi_a is the test shape.
                            //
                            // This appears in:
                            //
                            // - (du_c u_j, d_j v_c)
                            //
                            double u_dot_grad_Na = 0.0;
                            for (int j = 0; j < dim; j++) {
                                u_dot_grad_Na += u_n(j) * dNu(a, j);
                            }

                            for (int c = 0; c < vdim; c++) {
                                const int ia = VDofIndex(ndof_u, vdim, a, c, ordering);

                                for (int k = 0; k < vdim; k++) {
                                    const int ib = VDofIndex(ndof_u, vdim, b, k, ordering);

                                    double conv_jac = 0.0;

                                    // Variation in the first term:
                                    //
                                    // - (du \otimes u, grad v)
                                    //
                                    // Row component c, column component k:
                                    //
                                    // - delta_{ck} phi_b sum_j u_j d_j phi_a
                                    //
                                    if (c == k) {
                                        conv_jac += -Nu(b) * u_dot_grad_Na;
                                    }

                                    // Variation in the second term:
                                    //
                                    // - (u \otimes du, grad v)
                                    //
                                    // Row component c, column component k:
                                    //
                                    // - u_c phi_b d_k phi_a
                                    //
                                    conv_jac += -u_n(c) * Nu(b) * dNu(a, k);

                                    Auu(ia, ib) += conv_jac * wdet;
                                }
                            }
                        }
                    }
                }

                // ------------------------------------------------------------
                // Aup block
                //
                // d/dp [ -(p, div v) ]
                // ------------------------------------------------------------
                for (int a = 0; a < ndof_u; a++) {
                    for (int c = 0; c < vdim; c++) {
                        const int ia = VDofIndex(ndof_u, vdim, a, c, ordering);

                        for (int b = 0; b < ndof_p; b++) {
                            Aup(ia, b) += -Np(b) * dNu(a, c) * wdet;
                        }
                    }
                }

                // ------------------------------------------------------------
                // Apu block
                //
                // d/du [ (q, div u) ]
                // ------------------------------------------------------------
                for (int a = 0; a < ndof_p; a++) {
                    for (int b = 0; b < ndof_u; b++) {
                        for (int c = 0; c < vdim; c++) {
                            const int ib = VDofIndex(ndof_u, vdim, b, c, ordering);

                            Apu(a, ib) += Np(a) * dNu(b, c) * wdet;
                        }
                    }
                }

                // A11 block remains zero.
            }
        }
    };

    class NSEBlockIntegBDF2OutletConvectiveFlux : public NSEIntegratorBase {
    public:
        NSEBlockIntegBDF2OutletConvectiveFlux(const InputData& idata,
                                                 const TimeLevelFields& tlf,
                                                 const int vdim,
                                                 const mfem::Ordering::Type ordering,
                                                 mfem::VectorCoefficient* f_coeff = nullptr)
            : NSEIntegratorBase(idata, tlf, vdim, ordering, f_coeff) {
        }

        void AssembleFaceVector(const mfem::Array<const mfem::FiniteElement*>& el1,
                                const mfem::Array<const mfem::FiniteElement*>& el2,
                                mfem::FaceElementTransformations& Tr,
                                const mfem::Array<const mfem::Vector*>& elfun,
                                const mfem::Array<mfem::Vector*>& elvec) override {
            // Do not verify el2 here. This integrator is added through
            // AddBdrFaceIntegrator(...), so MFEM already restricts it to
            // boundary faces. In ParBlockNonlinearForm, el2 may still be
            // structurally non-empty.

            const mfem::FiniteElement& el_u = *el1[0];
            const mfem::FiniteElement& el_p = *el1[1];

            const int ndof_u = el_u.GetDof();
            const int ndof_p = el_p.GetDof();
            const int dim = Tr.Elem1->GetSpaceDim();

            MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");

            elvec[0]->SetSize(vdim * ndof_u);
            *elvec[0] = 0.0;

            // Pressure residual gets no outlet convective-flux contribution,
            // but it must still have the correct block size.
            elvec[1]->SetSize(ndof_p);
            *elvec[1] = 0.0;

            if (idata.flow_properties.disable_convection) {
                return;
            }

            const mfem::Vector& eu = *elfun[0];

            mfem::Vector Nu(ndof_u);
            mfem::Vector u(vdim);
            mfem::Vector nor(dim);

            const int order = 2 * el_u.GetOrder();
            const mfem::IntegrationRule* ir =
                &mfem::IntRules.Get(Tr.GetGeometryType(), order);

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint& ip_face = ir->IntPoint(iq);

                Tr.Face->SetIntPoint(&ip_face);

                // Map face integration point to element-1 integration point.
                mfem::IntegrationPoint ip_el;
                Tr.Loc1.Transform(ip_face, ip_el);
                Tr.Elem1->SetIntPoint(&ip_el);

                el_u.CalcShape(ip_el, Nu);

                EvalVectorAtIP(eu, ndof_u, vdim, ordering, Nu, u);

                // CalcOrtho gives the scaled physical normal:
                //
                //     nor = n |J_face|
                //
                // Therefore the surface measure is already contained in nor.
                // Only multiply by the quadrature weight.
                mfem::CalcOrtho(Tr.Face->Jacobian(), nor);

                double u_dot_n = 0.0;
                for (int j = 0; j < dim; j++) {
                    u_dot_n += u(j) * nor(j);
                }

                const double w = ip_face.weight;

                // Boundary flux:
                //
                // <(u \otimes u)n, v>
                // =
                // <u (u . n), v>
                for (int a = 0; a < ndof_u; a++) {
                    for (int c = 0; c < vdim; c++) {
                        const int ia = VDofIndex(ndof_u, vdim, a, c, ordering);

                        (*elvec[0])(ia) += Nu(a) * u(c) * u_dot_n * w;
                    }
                }
            }
        }

        void AssembleFaceGrad(const mfem::Array<const mfem::FiniteElement*>& el1,
                              const mfem::Array<const mfem::FiniteElement*>& el2,
                              mfem::FaceElementTransformations& Tr,
                              const mfem::Array<const mfem::Vector*>& elfun,
                              const mfem::Array2D<mfem::DenseMatrix*>& elmat) override {
            // Do not verify el2 here. This is called through AddBdrFaceIntegrator.

            const mfem::FiniteElement& el_u = *el1[0];
            const mfem::FiniteElement& el_p = *el1[1];

            const int dof_u = el_u.GetDof();
            const int dof_p = el_p.GetDof();
            const int dim = Tr.Elem1->GetSpaceDim();

            MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");

            elmat(0, 0)->SetSize(vdim * dof_u, vdim * dof_u);
            *elmat(0, 0) = 0.0;

            // These blocks get no nonzero contribution from the outlet
            // convective flux, but they must still be correctly sized.
            elmat(0, 1)->SetSize(vdim * dof_u, dof_p);
            *elmat(0, 1) = 0.0;

            elmat(1, 0)->SetSize(dof_p, vdim * dof_u);
            *elmat(1, 0) = 0.0;

            elmat(1, 1)->SetSize(dof_p, dof_p);
            *elmat(1, 1) = 0.0;

            if (idata.flow_properties.disable_convection) {
                return;
            }

            mfem::DenseMatrix& Auu = *elmat(0, 0);

            const mfem::Vector& eu = *elfun[0];

            mfem::Vector Nu(dof_u);
            mfem::Vector u(vdim);
            mfem::Vector nor(dim);

            const int order = 2 * el_u.GetOrder();
            const mfem::IntegrationRule* ir =
                &mfem::IntRules.Get(Tr.GetGeometryType(), order);

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint& ip_face = ir->IntPoint(iq);

                Tr.Face->SetIntPoint(&ip_face);

                mfem::IntegrationPoint ip_el;
                Tr.Loc1.Transform(ip_face, ip_el);
                Tr.Elem1->SetIntPoint(&ip_el);

                el_u.CalcShape(ip_el, Nu);

                EvalVectorAtIP(eu, dof_u, vdim, ordering, Nu, u);

                mfem::CalcOrtho(Tr.Face->Jacobian(), nor);

                double u_dot_n = 0.0;
                for (int j = 0; j < dim; j++) {
                    u_dot_n += u(j) * nor(j);
                }

                const double w = ip_face.weight;

                // Residual:
                //
                // R_c = u_c (u . n)
                //
                // Jacobian:
                //
                // dR_c / du_k = delta_ck (u . n) + u_c n_k
                //
                // Since nor is already scaled by |J_face|, this is really:
                //
                // dR_c / du_k = delta_ck (u . nor) + u_c nor_k
                for (int a = 0; a < dof_u; a++) {
                    for (int b = 0; b < dof_u; b++) {
                        for (int c = 0; c < vdim; c++) {
                            const int ia = VDofIndex(dof_u, vdim, a, c, ordering);

                            for (int k = 0; k < vdim; k++) {
                                const int ib = VDofIndex(dof_u, vdim, b, k, ordering);

                                double jac_ck = 0.0;

                                if (c == k) {
                                    jac_ck += u_dot_n;
                                }

                                jac_ck += u(c) * nor(k);

                                Auu(ia, ib) += Nu(a) * Nu(b) * jac_ck * w;
                            }
                        }
                    }
                }
            }
        }
    };
}
