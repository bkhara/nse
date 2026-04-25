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

    class NSEBlockIntegBDF2VMSConservative : public NSEIntegratorBase {
    private:
        static double Dot(const mfem::Vector &a, const mfem::Vector &b) {
            double v = 0.0;
            for (int i = 0; i < a.Size(); i++) {
                v += a(i) * b(i);
            }
            return v;
        }

        static double RowDot(const mfem::DenseMatrix &A,
                             const int row,
                             const mfem::Vector &x) {
            double v = 0.0;
            for (int j = 0; j < x.Size(); j++) {
                v += A(row, j) * x(j);
            }
            return v;
        }

        static double GradShapeDotU(const mfem::DenseMatrix &dN,
                                    const int a,
                                    const mfem::Vector &u) {
            double v = 0.0;
            for (int j = 0; j < u.Size(); j++) {
                v += u(j) * dN(a, j);
            }
            return v;
        }

        static double ElementLength(mfem::ElementTransformation &T,
                                    const mfem::FiniteElement &el) {
            const int dim = T.GetSpaceDim();

            // This is only a simple local length-scale estimate.
            // If you already have a preferred h_K, replace this.
            const double meas = std::abs(T.Weight());

            if (dim == 1) {
                return meas;
            } else if (dim == 2) {
                return std::sqrt(meas);
            } else {
                return std::cbrt(meas);
            }
        }

        void ComputeTau(const mfem::Vector &u,
                        const double h,
                        const double dt,
                        const double nu,
                        double &tauM,
                        double &tauC) const {
            const double unorm = std::sqrt(Dot(u, u));

            // Standard-ish residual-based VMS/SUPG scaling:
            //
            // tauM = [ (2/dt)^2 + (2|u|/h)^2 + (C_I nu / h^2)^2 ]^{-1/2}
            //
            // C_I is method/order dependent. 4 is a reasonable mild starting value.
            const double CI = 4.0;

            const double t_time = 2.0 / dt;
            const double t_conv = 2.0 * unorm / h;
            const double t_diff = CI * nu / (h * h);

            tauM = 1.0 / std::sqrt(t_time * t_time
                                   + t_conv * t_conv
                                   + t_diff * t_diff
                                   + 1.0e-30);

            // Continuity/grad-div stabilization scale.
            //
            // This can be strong. If the method feels too over-constrained,
            // multiply this by 0.1 or expose a user parameter.
            tauC = h * h / (tauM + 1.0e-30);
        }

    public:
        NSEBlockIntegBDF2VMSConservative(const InputData &idata,
                                         const TimeLevelFields &tlf,
                                         const int vdim,
                                         const mfem::Ordering::Type ordering,
                                         mfem::VectorCoefficient *f_coeff = nullptr)
            : NSEIntegratorBase(idata, tlf, vdim, ordering, f_coeff) {
        }

        void AssembleElementVector(const mfem::Array<const mfem::FiniteElement *> &el,
                                   mfem::ElementTransformation &T,
                                   const mfem::Array<const mfem::Vector *> &elfun,
                                   const mfem::Array<mfem::Vector *> &elvec) override {
            const mfem::FiniteElement &el_u = *el[0];
            const mfem::FiniteElement &el_p = *el[1];

            const int ndof_u = el_u.GetDof();
            const int ndof_p = el_p.GetDof();
            const int dim = T.GetSpaceDim();

            MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");

            elvec[0]->SetSize(vdim * ndof_u);
            *elvec[0] = 0.0;

            elvec[1]->SetSize(ndof_p);
            *elvec[1] = 0.0;

            const mfem::Vector &eu = *elfun[0];
            const mfem::Vector &ep = *elfun[1];

            mfem::Vector Nu(ndof_u), Np(ndof_p);
            mfem::DenseMatrix dNu(ndof_u, dim);
            mfem::DenseMatrix dNp(ndof_p, dim);

            mfem::Vector u(vdim), u_nm1(vdim), u_nm2(vdim), f(vdim);
            mfem::DenseMatrix grad_u(vdim, dim);

            mfem::Vector RM(vdim);

            const int e = T.ElementNo;
            const mfem::IntegrationRule *ir = &tlf.femach.qspace->GetIntRule(e);

            const double ctime = tlf.GetTime();
            const double dt = tlf.GetTimeStep();
            const double nu = idata.flow_properties.nu;
            const double alpha = 3.0 / (2.0 * dt);

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint &ip = ir->IntPoint(iq);
                T.SetIntPoint(&ip);

                el_u.CalcShape(ip, Nu);
                el_u.CalcPhysDShape(T, dNu);

                el_p.CalcShape(ip, Np);
                el_p.CalcPhysDShape(T, dNp);

                const double wdet = ip.weight * T.Weight();

                EvalVectorAtIP(eu, ndof_u, vdim, ordering, Nu, u);
                EvalVectorGradAtIP(eu, ndof_u, vdim, ordering, dNu, grad_u);

                tlf.prev_1.u.GetVectorValue(T, ip, u_nm1);
                tlf.prev_2.u.GetVectorValue(T, ip, u_nm2);

                f.SetSize(vdim);
                f = 0.0;
                if (f_coeff) {
                    f_coeff->SetTime(ctime);
                    f_coeff->Eval(f, T, ip);
                }

                const double div_u = Divergence(grad_u);

                const double h = ElementLength(T, el_u);

                double tauM = 0.0;
                double tauC = 0.0;
                ComputeTau(u, h, dt, nu, tauM, tauC);

                // ------------------------------------------------------------
                // Conservative strong momentum residual for stabilization:
                //
                // R_M =
                //   (3u^{n+1} - 4u^n + u^{n-1})/(2dt)
                // + div(u^{n+1} tensor u^{n+1})
                // + grad p^{n+1}
                // - f^{n+1}
                //
                // Expanded:
                //
                // div(u tensor u)_c =
                //   (u . grad) u_c + u_c div(u)
                //
                // We intentionally omit -nu Delta u here to avoid second
                // derivatives of the velocity basis. The Galerkin weak form
                // already contains the viscous term.
                // ------------------------------------------------------------
                RM.SetSize(vdim);
                RM = 0.0;

                for (int c = 0; c < vdim; c++) {
                    // BDF2 part
                    RM(c) += alpha * u(c);
                    RM(c) += -(4.0 * u_nm1(c) - u_nm2(c)) / (2.0 * dt);

                    // Conservative convection:
                    //
                    // div(u tensor u)_c = (u . grad) u_c + u_c div(u)
                    if (!idata.flow_properties.disable_convection) {
                        RM(c) += RowDot(grad_u, c, u);
                        RM(c) += u(c) * div_u;
                    }

                    // grad p
                    for (int b = 0; b < ndof_p; b++) {
                        RM(c) += ep(b) * dNp(b, c);
                    }

                    // forcing
                    RM(c) += -f(c);
                }

                // ------------------------------------------------------------
                // Velocity residual contribution:
                //
                // SUPG:
                //   tauM (R_M, u . grad v)
                //
                // LSIC / grad-div:
                //   tauC (div u, div v)
                // ------------------------------------------------------------
                for (int a = 0; a < ndof_u; a++) {
                    const double u_dot_grad_Na = GradShapeDotU(dNu, a, u);

                    for (int c = 0; c < vdim; c++) {
                        const int ia = VDofIndex(ndof_u, vdim, a, c, ordering);

                        // SUPG contribution
                        (*elvec[0])(ia) +=
                                tauM * u_dot_grad_Na * RM(c) * wdet;

                        // LSIC / grad-div contribution
                        (*elvec[0])(ia) +=
                                tauC * dNu(a, c) * div_u * wdet;
                    }
                }

                // ------------------------------------------------------------
                // Pressure residual contribution:
                //
                // PSPG:
                //   tauM (R_M, grad q)
                // ------------------------------------------------------------
                for (int a = 0; a < ndof_p; a++) {
                    for (int c = 0; c < vdim; c++) {
                        (*elvec[1])(a) +=
                                tauM * dNp(a, c) * RM(c) * wdet;
                    }
                }
            }
        }

        void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement *> &el,
                                 mfem::ElementTransformation &T,
                                 const mfem::Array<const mfem::Vector *> &elfun,
                                 const mfem::Array2D<mfem::DenseMatrix *> &elmat) override {
            const mfem::FiniteElement &el_u = *el[0];
            const mfem::FiniteElement &el_p = *el[1];

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

            mfem::DenseMatrix &Auu = *elmat(0, 0);
            mfem::DenseMatrix &Aup = *elmat(0, 1);
            mfem::DenseMatrix &Apu = *elmat(1, 0);
            mfem::DenseMatrix &App = *elmat(1, 1);

            const mfem::Vector &eu = *elfun[0];
            const mfem::Vector &ep = *elfun[1];

            mfem::Vector Nu(ndof_u), Np(ndof_p);
            mfem::DenseMatrix dNu(ndof_u, dim);
            mfem::DenseMatrix dNp(ndof_p, dim);

            mfem::Vector u(vdim), u_nm1(vdim), u_nm2(vdim), f(vdim);
            mfem::DenseMatrix grad_u(vdim, dim);

            mfem::Vector RM(vdim);

            const int e = T.ElementNo;
            const mfem::IntegrationRule *ir = &tlf.femach.qspace->GetIntRule(e);

            const double ctime = tlf.GetTime();
            const double dt = tlf.GetTimeStep();
            const double nu = idata.flow_properties.nu;
            const double alpha = 3.0 / (2.0 * dt);

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint &ip = ir->IntPoint(iq);
                T.SetIntPoint(&ip);

                el_u.CalcShape(ip, Nu);
                el_u.CalcPhysDShape(T, dNu);

                el_p.CalcShape(ip, Np);
                el_p.CalcPhysDShape(T, dNp);

                const double wdet = ip.weight * T.Weight();

                EvalVectorAtIP(eu, ndof_u, vdim, ordering, Nu, u);
                EvalVectorGradAtIP(eu, ndof_u, vdim, ordering, dNu, grad_u);

                tlf.prev_1.u.GetVectorValue(T, ip, u_nm1);
                tlf.prev_2.u.GetVectorValue(T, ip, u_nm2);

                f.SetSize(vdim);
                f = 0.0;
                if (f_coeff) {
                    f_coeff->SetTime(ctime);
                    f_coeff->Eval(f, T, ip);
                }

                const double div_u = Divergence(grad_u);

                const double h = ElementLength(T, el_u);

                double tauM = 0.0;
                double tauC = 0.0;
                ComputeTau(u, h, dt, nu, tauM, tauC);

                // Conservative strong momentum residual.
                RM.SetSize(vdim);
                RM = 0.0;

                for (int c = 0; c < vdim; c++) {
                    RM(c) += alpha * u(c);
                    RM(c) += -(4.0 * u_nm1(c) - u_nm2(c)) / (2.0 * dt);

                    if (!idata.flow_properties.disable_convection) {
                        RM(c) += RowDot(grad_u, c, u);
                        RM(c) += u(c) * div_u;
                    }

                    for (int b = 0; b < ndof_p; b++) {
                        RM(c) += ep(b) * dNp(b, c);
                    }

                    RM(c) += -f(c);
                }

                // ------------------------------------------------------------
                // Jacobian of conservative R_M.
                //
                // R_c =
                //   alpha u_c
                // + (u . grad) u_c
                // + u_c div(u)
                // + grad_c p
                // - history - f
                //
                // Variation wrt velocity component k:
                //
                // dR_c/du_k =
                //   alpha delta_ck phi_b
                // + phi_b grad_k u_c
                // + delta_ck (u . grad phi_b)
                // + delta_ck phi_b div(u)
                // + u_c d_k phi_b
                //
                // where d_k phi_b contributes to div(delta u).
                //
                // tauM and tauC are lagged here. That is a quasi-Newton
                // stabilization Jacobian. It is usually much simpler and robust.
                // ------------------------------------------------------------

                // ------------------------------------------------------------
                // Auu and Aup from SUPG velocity equation:
                //
                // tauM (R_M, u . grad v)
                //
                // Linearized wrt u:
                //   tauM (dR_M, u . grad v)
                // + tauM (R_M, du . grad v)
                //
                // Linearized wrt p:
                //   tauM (grad dp, u . grad v)
                //
                // Plus LSIC:
                //   tauC (div u, div v)
                // ------------------------------------------------------------
                for (int a = 0; a < ndof_u; a++) {
                    const double u_dot_grad_Na = GradShapeDotU(dNu, a, u);

                    for (int b = 0; b < ndof_u; b++) {
                        const double u_dot_grad_Nb = GradShapeDotU(dNu, b, u);

                        for (int c = 0; c < vdim; c++) {
                            const int ia = VDofIndex(ndof_u, vdim, a, c, ordering);

                            for (int k = 0; k < vdim; k++) {
                                const int ib = VDofIndex(ndof_u, vdim, b, k, ordering);

                                double dRM_du = 0.0;

                                // Time derivative: alpha u_c
                                if (c == k) {
                                    dRM_du += alpha * Nu(b);
                                }

                                if (!idata.flow_properties.disable_convection) {
                                    // d[(u . grad) u_c]
                                    //
                                    // = delta u . grad u_c
                                    // + u . grad delta u_c

                                    // delta u_k = phi_b in component k:
                                    // phi_b * partial_k u_c
                                    dRM_du += Nu(b) * grad_u(c, k);

                                    // u . grad(delta u_c)
                                    if (c == k) {
                                        dRM_du += u_dot_grad_Nb;
                                    }

                                    // d[u_c div(u)]
                                    //
                                    // = delta u_c div(u)
                                    // + u_c div(delta u)

                                    // delta u_c div(u)
                                    if (c == k) {
                                        dRM_du += Nu(b) * div_u;
                                    }

                                    // u_c div(delta u)
                                    // div(delta u) = partial_k phi_b
                                    dRM_du += u(c) * dNu(b, k);
                                }

                                // tauM (dR_c, u . grad v_c)
                                Auu(ia, ib) +=
                                        tauM * u_dot_grad_Na * dRM_du * wdet;

                                // tauM (R_c, delta u . grad v_c)
                                //
                                // This is the derivative of the streamline test
                                // direction u . grad(v_c).
                                if (!idata.flow_properties.disable_convection) {
                                    Auu(ia, ib) +=
                                            tauM * RM(c) * Nu(b) * dNu(a, k) * wdet;
                                }

                                // LSIC / grad-div:
                                //
                                // tauC (div delta u, div v)
                                Auu(ia, ib) +=
                                        tauC * dNu(b, k) * dNu(a, c) * wdet;
                            }
                        }
                    }

                    // Aup from SUPG:
                    //
                    // tauM (grad dp, u . grad v)
                    for (int c = 0; c < vdim; c++) {
                        const int ia = VDofIndex(ndof_u, vdim, a, c, ordering);

                        for (int b = 0; b < ndof_p; b++) {
                            Aup(ia, b) +=
                                    tauM * u_dot_grad_Na * dNp(b, c) * wdet;
                        }
                    }
                }

                // ------------------------------------------------------------
                // Apu and App from PSPG pressure equation:
                //
                // tauM (R_M, grad q)
                //
                // Linearized wrt u:
                //   tauM (dR_M, grad q)
                //
                // Linearized wrt p:
                //   tauM (grad dp, grad q)
                // ------------------------------------------------------------
                for (int a = 0; a < ndof_p; a++) {
                    for (int b = 0; b < ndof_u; b++) {
                        const double u_dot_grad_Nb = GradShapeDotU(dNu, b, u);

                        for (int c = 0; c < vdim; c++) {
                            for (int k = 0; k < vdim; k++) {
                                const int ib = VDofIndex(ndof_u, vdim, b, k, ordering);

                                double dRM_du = 0.0;

                                // Time derivative
                                if (c == k) {
                                    dRM_du += alpha * Nu(b);
                                }

                                if (!idata.flow_properties.disable_convection) {
                                    // d[(u . grad) u_c]
                                    dRM_du += Nu(b) * grad_u(c, k);

                                    if (c == k) {
                                        dRM_du += u_dot_grad_Nb;
                                    }

                                    // d[u_c div(u)]
                                    if (c == k) {
                                        dRM_du += Nu(b) * div_u;
                                    }

                                    dRM_du += u(c) * dNu(b, k);
                                }

                                Apu(a, ib) +=
                                        tauM * dNp(a, c) * dRM_du * wdet;
                            }
                        }
                    }

                    // App from PSPG:
                    //
                    // tauM (grad dp, grad q)
                    for (int b = 0; b < ndof_p; b++) {
                        for (int c = 0; c < vdim; c++) {
                            App(a, b) +=
                                    tauM * dNp(a, c) * dNp(b, c) * wdet;
                        }
                    }
                }
            }
        }
    };
}
