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
    inline void EvalVectorAtIP(const mfem::Vector& elfun,
                               const int ndof,
                               const int vdim,
                               const mfem::Ordering::Type ordering,
                               const mfem::Vector& N,
                               mfem::Vector& u) {
        u.SetSize(vdim);
        u = 0.0;

        for (int c = 0; c < vdim; c++) {
            for (int a = 0; a < ndof; a++) {
                u(c) += elfun(VDofIndex(ndof, vdim, a, c, ordering)) * N(a);
            }
        }
    }

    inline void EvalVectorGradAtIP(const mfem::Vector& elfun,
                                   const int ndof,
                                   const int vdim,
                                   const mfem::Ordering::Type ordering,
                                   const mfem::DenseMatrix& dN,
                                   mfem::DenseMatrix& grad_u) {
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

    inline double EvalScalarAtIP(const mfem::Vector& elfun,
                                 const mfem::Vector& N) {
        double val = 0.0;
        for (int a = 0; a < N.Size(); a++) {
            val += elfun(a) * N(a);
        }
        return val;
    }

    // -------------------------------------------------------------------
    // Shared strong-form Laplacian helpers.
    //
    // MFEM's FiniteElement::CalcPhysLaplacian(T, Laplacian) evaluates the
    // physical-space Laplacian (trace of the Hessian) of every scalar
    // shape function of `el` at the current integration point, and is
    // implemented for the standard H1 element families (segment/quad/hex,
    // triangle/tet) in 1D/2D/3D -- i.e. exactly the elements this code
    // uses (see FEMachinery.cpp, H1_FECollection). These are the single,
    // shared implementation of the "strong viscous term" second-derivative
    // evaluation: every integrator that needs -nu*Delta(u) should call
    // these instead of defining its own per-class private copy.
    // -------------------------------------------------------------------

    // Delta(u) for the vector field u = sum_a elfun_a N_a(x), given the
    // local dof vector `elfun` (size ndof*vdim, laid out per `ordering`)
    // and the per-shape-function Laplacians `lap_N` (dof-sized, from
    // el.CalcPhysLaplacian(T, lap_N)). Mirrors EvalVectorAtIP / EvalVectorGradAtIP.
    inline void EvalVectorLaplacianAtIP(const mfem::Vector& elfun,
                                        const int ndof,
                                        const int vdim,
                                        const mfem::Ordering::Type ordering,
                                        const mfem::Vector& lap_N,
                                        mfem::Vector& lap_u) {
        lap_u.SetSize(vdim);
        lap_u = 0.0;

        for (int c = 0; c < vdim; c++) {
            for (int a = 0; a < ndof; a++) {
                lap_u(c) += elfun(VDofIndex(ndof, vdim, a, c, ordering)) * lap_N(a);
            }
        }
    }

    // Delta(u) for a vector-valued GridFunction `gf` (e.g. tlf.current.u)
    // at the current integration point. For RHS-style integrators that
    // know the field values but don't have an elfun of gf's own space --
    // note `gf` may live on a different FE space than the integrator's own
    // `el` (e.g. the PPE RHS integrates over the pressure space but needs
    // Delta of the velocity field), so this looks up gf's own element and
    // ordering rather than trusting the caller's.
    inline void EvalGridFunctionVectorLaplacianAtIP(const mfem::GridFunction& gf,
                                                    mfem::ElementTransformation& T,
                                                    const int vdim,
                                                    mfem::Vector& lap_u) {
        lap_u.SetSize(vdim);
        lap_u = 0.0;

        const mfem::FiniteElement* gf_el = gf.FESpace()->GetFE(T.ElementNo);

        mfem::Vector lap_N(gf_el->GetDof());
        gf_el->CalcPhysLaplacian(T, lap_N);

        mfem::Array<int> vdofs;
        gf.FESpace()->GetElementVDofs(T.ElementNo, vdofs);

        mfem::Vector gf_elfun;
        gf.GetSubVector(vdofs, gf_elfun);

        EvalVectorLaplacianAtIP(gf_elfun, gf_el->GetDof(), vdim,
                                 gf.FESpace()->GetOrdering(), lap_N, lap_u);
    }

    // Shared helper: evaluate p_hat = p0*p^n + p1*p^{n-1} and its gradient
    // at the current integration point, using the coefficients from
    // GetProjectionCoefficients(). This is "p_h*" in the projection-method
    // algorithm (Require: extrapolated pressure p_h*), and every integrator
    // that needs it (momentum predictor, VUE RHS, PPE RHS) must use this
    // same helper so that all three steps of a given time step agree on
    // p_h* -- currently p1 == 0 for both defined schemes, so p_hat == p^n
    // numerically, but the three integrators should not each independently
    // decide how to combine prev_1/prev_2 pressure.
    inline void EvalPressureHatAndGradAtIP(const TimeLevelFields& tlf,
                                           mfem::ElementTransformation& T,
                                           const mfem::IntegrationPoint& ip,
                                           const ProjectionCoefficients& pc,
                                           const int vdim,
                                           double& p_hat,
                                           mfem::Vector& grad_p_hat) {
        grad_p_hat.SetSize(vdim);
        grad_p_hat = 0.0;

        const double p_n = tlf.prev_1.p.GetValue(T, ip);
        const double p_nm1 = tlf.prev_2.p.GetValue(T, ip);

        p_hat = pc.p0 * p_n + pc.p1 * p_nm1;

        mfem::Vector grad_p_n(vdim);
        mfem::Vector grad_p_nm1(vdim);

        tlf.prev_1.p.GetGradient(T, grad_p_n);
        tlf.prev_2.p.GetGradient(T, grad_p_nm1);

        for (int c = 0; c < vdim; c++) {
            grad_p_hat(c) = pc.p0 * grad_p_n(c) + pc.p1 * grad_p_nm1(c);
        }
    }

    inline void EvalConvection(const mfem::Vector& u,
                               const mfem::DenseMatrix& grad_u,
                               mfem::Vector& conv) {
        const int dim = u.Size();

        conv.SetSize(dim);
        conv = 0.0;

        // conv_i = u_j d_j u_i
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                conv(i) += u(j) * grad_u(i, j);
            }
        }
    }

    inline double Divergence(const mfem::DenseMatrix& grad_u) {
        const int vdim = grad_u.Height();
        double div_u = 0.0;
        for (int c = 0; c < vdim; c++) { div_u += grad_u(c, c); }
        return div_u;
    }

    static void CalcTau(mfem::ElementTransformation& T,
                                  const mfem::Vector& u,
                                  const double nu,
                                  const double Ci_f,
                                  const double dt,
                                  double& tauM,
                                  double& tauC) {
        const int dim = T.GetSpaceDim();

        MFEM_VERIFY(T.GetDimension() == dim,
                    "CalcTauMFEMWithNu assumes dim == space dim.");
        MFEM_VERIFY(u.Size() == dim,
                    "Velocity vector has wrong dimension.");

        const mfem::DenseMatrix& J = T.Jacobian();

        mfem::DenseMatrix ksiX(dim, dim);
        mfem::CalcInverse(J, ksiX);
        ksiX *= 2.0;

        mfem::DenseMatrix Ge(dim, dim);
        Ge = 0.0;

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                for (int k = 0; k < dim; k++) {
                    Ge(i, j) += ksiX(k, i) * ksiX(k, j);
                }
            }
        }

        double u_Gu = 0.0;
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                u_Gu += u(i) * Ge(i, j) * u(j);
            }
        }

        double G_G = 0.0;
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                G_G += Ci_f * nu * nu * Ge(i, j) * Ge(i, j);
            }
        }

        const double denom_squared =
            (4.0 / (dt * dt)) + u_Gu + G_G;

        tauM = 1.0 / std::sqrt(denom_squared + 1.0e-30);

        mfem::Vector ge(dim);
        ge = 0.0;

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                ge(i) += ksiX(j, i);
            }
        }

        double g_g = 0.0;
        for (int i = 0; i < dim; i++) {
            g_g += ge(i) * ge(i);
        }

        tauC = 1.0 / (tauM * g_g + 1.0e-30);
    }

    class NSEIntegratorBase : public mfem::BlockNonlinearFormIntegrator {
    protected:
        const InputData &idata;
        const TimeLevelFields &tlf;

        mfem::VectorCoefficient *f_coeff;

        const int vdim;
        const mfem::Ordering::Type ordering;

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

    class NSEBlockIntegCrankNicolson : public NSEIntegratorBase {
    private:
        const double theta;

    public:
        NSEBlockIntegCrankNicolson(const InputData& idata,
                                   const TimeLevelFields& tlf,
                                   const int vdim,
                                   const mfem::Ordering::Type ordering,
                                   mfem::VectorCoefficient* f_coeff = nullptr)
            : NSEIntegratorBase(idata, tlf, vdim, ordering, f_coeff),
              theta(0.5) {
            MFEM_VERIFY(theta >= 0.0 && theta <= 1.0,
                        "Crank-Nicolson theta must be in [0, 1].");
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

            mfem::Vector u_np1(vdim), u_n(vdim);
            mfem::Vector conv_np1(vdim), conv_n(vdim);
            mfem::Vector f_np1(vdim), f_n(vdim);

            mfem::DenseMatrix grad_u_np1(vdim, dim);
            mfem::DenseMatrix grad_u_n(vdim, dim);

            const int e = T.ElementNo;
            const mfem::IntegrationRule* ir = &tlf.femach.qspace->GetIntRule(e);

            const double ctime = tlf.GetTime(); // t_{n+1}
            const double dt = tlf.GetTimeStep();
            const double ptime = ctime - dt; // t_n
            const double nu = idata.flow_properties.nu;

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint& ip = ir->IntPoint(iq);
                T.SetIntPoint(&ip);

                el_u.CalcShape(ip, Nu);
                el_u.CalcPhysDShape(T, dNu);
                el_p.CalcShape(ip, Np);

                const double wdet = ip.weight * T.Weight();

                // Current Newton unknowns: u^{n+1}, p^{n+1}
                EvalVectorAtIP(eu, ndof_u, vdim, ordering, Nu, u_np1);
                EvalVectorGradAtIP(eu, ndof_u, vdim, ordering, dNu, grad_u_np1);
                const double p_np1 = EvalScalarAtIP(ep, Np);

                // Previous accepted fields: u^n, p^n
                tlf.prev_1.u.GetVectorValue(T, ip, u_n);
                tlf.prev_1.u.GetVectorGradient(T, grad_u_n);
                const double p_n = tlf.prev_1.p.GetValue(T, ip);

                // Forcing at t_{n+1} and t_n
                f_np1.SetSize(vdim);
                f_np1 = 0.0;

                f_n.SetSize(vdim);
                f_n = 0.0;

                if (f_coeff) {
                    f_coeff->SetTime(ctime);
                    f_coeff->Eval(f_np1, T, ip);

                    f_coeff->SetTime(ptime);
                    f_coeff->Eval(f_n, T, ip);
                }

                const double div_u_np1 = Divergence(grad_u_np1);

                conv_np1.SetSize(vdim);
                conv_np1 = 0.0;

                conv_n.SetSize(vdim);
                conv_n = 0.0;

                if (!idata.flow_properties.disable_convection) {
                    EvalConvection(u_np1, grad_u_np1, conv_np1);
                    EvalConvection(u_n, grad_u_n, conv_n);
                }

                const double p_theta =
                    theta * p_np1 + (1.0 - theta) * p_n;

                // ------------------------------------------------------------
                // Momentum residual:
                //
                // ((u^{n+1} - u^n)/dt, v)
                // + theta     ((u^{n+1}.grad)u^{n+1}, v)
                // + (1-theta) ((u^n.grad)u^n, v)
                // + nu (theta grad u^{n+1}
                //       + (1-theta) grad u^n, grad v)
                // - (theta p^{n+1} + (1-theta) p^n, div v)
                // - (theta f^{n+1} + (1-theta) f^n, v)
                //
                // Continuity:
                //
                // (q, div u^{n+1})
                // ------------------------------------------------------------
                for (int a = 0; a < ndof_u; a++) {
                    for (int c = 0; c < vdim; c++) {
                        const int ia = VDofIndex(ndof_u, vdim, a, c, ordering);

                        // Time derivative
                        (*elvec[0])(ia) +=
                            ((u_np1(c) - u_n(c)) / dt) * Nu(a) * wdet;

                        // Convective form
                        if (!idata.flow_properties.disable_convection) {
                            const double conv_theta =
                                theta * conv_np1(c) + (1.0 - theta) * conv_n(c);

                            (*elvec[0])(ia) += Nu(a) * conv_theta * wdet;
                        }

                        // Diffusion
                        for (int j = 0; j < dim; j++) {
                            const double grad_theta =
                                theta * grad_u_np1(c, j)
                                + (1.0 - theta) * grad_u_n(c, j);

                            (*elvec[0])(ia) +=
                                nu * grad_theta * dNu(a, j) * wdet;
                        }

                        // Pressure
                        (*elvec[0])(ia) +=
                            -p_np1 * dNu(a, c) * wdet;

                        // Forcing
                        const double f_theta =
                            theta * f_np1(c) + (1.0 - theta) * f_n(c);

                        (*elvec[0])(ia) +=
                            -f_theta * Nu(a) * wdet;
                    }
                }

                // Continuity residual
                for (int a = 0; a < ndof_p; a++) {
                    (*elvec[1])(a) += Np(a) * div_u_np1 * wdet;
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

            mfem::Vector u_np1(vdim);
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

                EvalVectorAtIP(eu, ndof_u, vdim, ordering, Nu, u_np1);
                EvalVectorGradAtIP(eu, ndof_u, vdim, ordering, dNu, grad_u_np1);

                // ------------------------------------------------------------
                // Auu block
                //
                // d/du [
                //   (1/dt) (u, v)
                // + theta ((u.grad)u, v)
                // + theta nu (grad u, grad v)
                // ]
                //
                // Convective form:
                //
                // C_c(u) = (u . grad) u_c
                //
                // dC_c/du_k =
                //   phi_b * grad_k u_c
                // + delta_ck * (u . grad phi_b)
                // ------------------------------------------------------------
                for (int a = 0; a < ndof_u; a++) {
                    for (int b = 0; b < ndof_u; b++) {
                        const double mass_ab =
                            (1.0 / dt) * Nu(a) * Nu(b) * wdet;

                        double diff_ab = 0.0;
                        for (int j = 0; j < dim; j++) {
                            diff_ab +=
                                theta * nu * dNu(a, j) * dNu(b, j) * wdet;
                        }

                        // Mass + diffusion: diagonal in velocity components
                        for (int c = 0; c < vdim; c++) {
                            const int ia = VDofIndex(ndof_u, vdim, a, c, ordering);
                            const int ib = VDofIndex(ndof_u, vdim, b, c, ordering);

                            Auu(ia, ib) += mass_ab + diff_ab;
                        }

                        if (!idata.flow_properties.disable_convection) {
                            double u_dot_grad_Nb = 0.0;
                            for (int j = 0; j < dim; j++) {
                                u_dot_grad_Nb += u_np1(j) * dNu(b, j);
                            }

                            for (int c = 0; c < vdim; c++) {
                                const int ia = VDofIndex(ndof_u, vdim, a, c, ordering);

                                for (int k = 0; k < vdim; k++) {
                                    const int ib = VDofIndex(ndof_u, vdim, b, k, ordering);

                                    double conv_jac = 0.0;

                                    // delta u . grad u_c
                                    conv_jac += Nu(b) * grad_u_np1(c, k);

                                    // u . grad(delta u_c)
                                    if (c == k) {
                                        conv_jac += u_dot_grad_Nb;
                                    }

                                    Auu(ia, ib) +=
                                        theta * Nu(a) * conv_jac * wdet;
                                }
                            }
                        }
                    }
                }

                // ------------------------------------------------------------
                // Aup block
                //
                // d/dp [ -(theta p^{n+1}, div v) ]
                // ------------------------------------------------------------
                for (int a = 0; a < ndof_u; a++) {
                    for (int c = 0; c < vdim; c++) {
                        const int ia = VDofIndex(ndof_u, vdim, a, c, ordering);

                        for (int b = 0; b < ndof_p; b++) {
                            Aup(ia, b) +=
                                - Np(b) * dNu(a, c) * wdet;
                        }
                    }
                }

                // ------------------------------------------------------------
                // Apu block
                //
                // d/du [ (q, div u^{n+1}) ]
                // ------------------------------------------------------------
                for (int a = 0; a < ndof_p; a++) {
                    for (int b = 0; b < ndof_u; b++) {
                        for (int c = 0; c < vdim; c++) {
                            const int ib = VDofIndex(ndof_u, vdim, b, c, ordering);

                            Apu(a, ib) += Np(a) * dNu(b, c) * wdet;
                        }
                    }
                }

                // App block remains zero.
            }
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

    class NSEBlockIntegBDF2SUPGConservativeAddOn
        : public NSEIntegratorBase {
    public:
        NSEBlockIntegBDF2SUPGConservativeAddOn(
            const InputData& idata,
            const TimeLevelFields& tlf,
            const int vdim,
            const mfem::Ordering::Type ordering,
            mfem::VectorCoefficient* f_coeff = nullptr)
            : NSEIntegratorBase(idata, tlf, vdim, ordering, f_coeff) {
        }

        void AssembleElementVector(
            const mfem::Array<const mfem::FiniteElement*>& el,
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
            mfem::DenseMatrix dNp(ndof_p, dim);

            mfem::Vector u(vdim), u_nm1(vdim), u_nm2(vdim), f(vdim);
            mfem::Vector conv(vdim), RM(vdim);
            mfem::DenseMatrix grad_u(vdim, dim);

            const int e = T.ElementNo;
            const mfem::IntegrationRule* ir = &tlf.femach.qspace->GetIntRule(e);

            const double ctime = tlf.GetTime();
            const double dt = tlf.GetTimeStep();
            const double nu = idata.flow_properties.nu;
            const double alpha = 3.0 / (2.0 * dt);

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint& ip = ir->IntPoint(iq);
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

                conv.SetSize(vdim);
                conv = 0.0;
                if (!idata.flow_properties.disable_convection) {
                    EvalConvection(u, grad_u, conv);
                }

                double tauM = 0.0;
                double tauC_unused = 0.0;
                CalcTau(T, u, nu, idata.vms_config.Ci, dt, tauM, tauC_unused);

                // Conservative strong momentum residual for SUPG:
                //
                // R_M =
                //   alpha u^{n+1}
                // - (4u^n - u^{n-1})/(2dt)
                // + div(u^{n+1} tensor u^{n+1})
                // + grad p^{n+1}
                // - f^{n+1}
                //
                // with
                //
                // div(u tensor u)_c = (u . grad) u_c + u_c div(u).
                //
                // We omit -nu Delta u to avoid requiring second derivatives.
                RM.SetSize(vdim);
                RM = 0.0;

                for (int c = 0; c < vdim; c++) {
                    RM(c) += alpha * u(c);
                    RM(c) += -(4.0 * u_nm1(c) - u_nm2(c)) / (2.0 * dt);

                    if (!idata.flow_properties.disable_convection) {
                        RM(c) += conv(c);
                        RM(c) += u(c) * div_u;
                    }

                    for (int b = 0; b < ndof_p; b++) {
                        RM(c) += ep(b) * dNp(b, c);
                    }

                    RM(c) += -f(c);
                }

                // SUPG contribution to the momentum equation:
                //
                // tauM (R_M, u . grad v)
                for (int a = 0; a < ndof_u; a++) {
                    double u_dot_grad_Na = 0.0;
                    for (int j = 0; j < dim; j++) {
                        u_dot_grad_Na += u(j) * dNu(a, j);
                    }

                    for (int c = 0; c < vdim; c++) {
                        const int ia = VDofIndex(ndof_u, vdim, a, c, ordering);

                        (*elvec[0])(ia) +=
                            tauM * u_dot_grad_Na * RM(c) * wdet;
                    }
                }

                // elvec[1] intentionally remains zero.
            }
        }

        void AssembleElementGrad(
            const mfem::Array<const mfem::FiniteElement*>& el,
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

            const mfem::Vector& eu = *elfun[0];
            const mfem::Vector& ep = *elfun[1];

            mfem::Vector Nu(ndof_u), Np(ndof_p);
            mfem::DenseMatrix dNu(ndof_u, dim);
            mfem::DenseMatrix dNp(ndof_p, dim);

            mfem::Vector u(vdim), u_nm1(vdim), u_nm2(vdim), f(vdim);
            mfem::Vector conv(vdim), RM(vdim);
            mfem::DenseMatrix grad_u(vdim, dim);

            const int e = T.ElementNo;
            const mfem::IntegrationRule* ir = &tlf.femach.qspace->GetIntRule(e);

            const double ctime = tlf.GetTime();
            const double dt = tlf.GetTimeStep();
            const double nu = idata.flow_properties.nu;
            const double alpha = 3.0 / (2.0 * dt);

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint& ip = ir->IntPoint(iq);
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

                conv.SetSize(vdim);
                conv = 0.0;
                if (!idata.flow_properties.disable_convection) {
                    EvalConvection(u, grad_u, conv);
                }

                double tauM = 0.0;
                double tauC_unused = 0.0;
                CalcTau(T, u, nu, idata.vms_config.Ci, dt, tauM, tauC_unused);

                RM.SetSize(vdim);
                RM = 0.0;

                for (int c = 0; c < vdim; c++) {
                    RM(c) += alpha * u(c);
                    RM(c) += -(4.0 * u_nm1(c) - u_nm2(c)) / (2.0 * dt);

                    if (!idata.flow_properties.disable_convection) {
                        RM(c) += conv(c);
                        RM(c) += u(c) * div_u;
                    }

                    for (int b = 0; b < ndof_p; b++) {
                        RM(c) += ep(b) * dNp(b, c);
                    }

                    RM(c) += -f(c);
                }

                for (int a = 0; a < ndof_u; a++) {
                    double u_dot_grad_Na = 0.0;
                    for (int j = 0; j < dim; j++) {
                        u_dot_grad_Na += u(j) * dNu(a, j);
                    }

                    for (int b = 0; b < ndof_u; b++) {
                        double u_dot_grad_Nb = 0.0;
                        for (int j = 0; j < dim; j++) {
                            u_dot_grad_Nb += u(j) * dNu(b, j);
                        }

                        for (int c = 0; c < vdim; c++) {
                            const int ia = VDofIndex(ndof_u, vdim, a, c, ordering);

                            for (int k = 0; k < vdim; k++) {
                                const int ib = VDofIndex(ndof_u, vdim, b, k, ordering);

                                double dRM_du = 0.0;

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

                                // tauM (dR_c, u . grad v_c)
                                Auu(ia, ib) +=
                                    tauM * u_dot_grad_Na * dRM_du * wdet;

                                // derivative of u . grad(v_c)
                                if (!idata.flow_properties.disable_convection) {
                                    Auu(ia, ib) +=
                                        tauM * RM(c) * Nu(b) * dNu(a, k) * wdet;
                                }
                            }
                        }
                    }

                    // pressure derivative:
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

                // Apu and App intentionally remain zero for SUPG-only integrator.
            }
        }
    };

    class NSEBlockIntegBDF2PSPGConservativeAddOn
        : public NSEIntegratorBase {
    public:
        NSEBlockIntegBDF2PSPGConservativeAddOn(
            const InputData& idata,
            const TimeLevelFields& tlf,
            const int vdim,
            const mfem::Ordering::Type ordering,
            mfem::VectorCoefficient* f_coeff = nullptr)
            : NSEIntegratorBase(idata, tlf, vdim, ordering, f_coeff) {
        }

        void AssembleElementVector(
            const mfem::Array<const mfem::FiniteElement*>& el,
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
            mfem::DenseMatrix dNp(ndof_p, dim);

            mfem::Vector u(vdim), u_nm1(vdim), u_nm2(vdim), f(vdim);
            mfem::Vector conv(vdim), RM(vdim);
            mfem::DenseMatrix grad_u(vdim, dim);

            const int e = T.ElementNo;
            const mfem::IntegrationRule* ir = &tlf.femach.qspace->GetIntRule(e);

            const double ctime = tlf.GetTime();
            const double dt = tlf.GetTimeStep();
            const double nu = idata.flow_properties.nu;
            const double alpha = 3.0 / (2.0 * dt);

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint& ip = ir->IntPoint(iq);
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

                conv.SetSize(vdim);
                conv = 0.0;
                if (!idata.flow_properties.disable_convection) {
                    EvalConvection(u, grad_u, conv);
                }

                double tauM = 0.0;
                double tauC_unused = 0.0;
                CalcTau(T, u, nu, idata.vms_config.Ci, dt, tauM, tauC_unused);

                RM.SetSize(vdim);
                RM = 0.0;

                for (int c = 0; c < vdim; c++) {
                    RM(c) += alpha * u(c);
                    RM(c) += -(4.0 * u_nm1(c) - u_nm2(c)) / (2.0 * dt);

                    if (!idata.flow_properties.disable_convection) {
                        RM(c) += conv(c);
                        RM(c) += u(c) * div_u;
                    }

                    for (int b = 0; b < ndof_p; b++) {
                        RM(c) += ep(b) * dNp(b, c);
                    }

                    RM(c) += -f(c);
                }

                // PSPG contribution to the pressure/continuity equation:
                //
                // tauM (R_M, grad q)
                for (int a = 0; a < ndof_p; a++) {
                    for (int c = 0; c < vdim; c++) {
                        (*elvec[1])(a) +=
                            tauM * dNp(a, c) * RM(c) * wdet;
                    }
                }

                // elvec[0] intentionally remains zero.
            }
        }

        void AssembleElementGrad(
            const mfem::Array<const mfem::FiniteElement*>& el,
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

            mfem::DenseMatrix& Apu = *elmat(1, 0);
            mfem::DenseMatrix& App = *elmat(1, 1);

            const mfem::Vector& eu = *elfun[0];
            const mfem::Vector& ep = *elfun[1];

            mfem::Vector Nu(ndof_u), Np(ndof_p);
            mfem::DenseMatrix dNu(ndof_u, dim);
            mfem::DenseMatrix dNp(ndof_p, dim);

            mfem::Vector u(vdim), u_nm1(vdim), u_nm2(vdim), f(vdim);
            mfem::Vector conv(vdim), RM(vdim);
            mfem::DenseMatrix grad_u(vdim, dim);

            const int e = T.ElementNo;
            const mfem::IntegrationRule* ir = &tlf.femach.qspace->GetIntRule(e);

            const double ctime = tlf.GetTime();
            const double dt = tlf.GetTimeStep();
            const double nu = idata.flow_properties.nu;
            const double alpha = 3.0 / (2.0 * dt);

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint& ip = ir->IntPoint(iq);
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

                conv.SetSize(vdim);
                conv = 0.0;
                if (!idata.flow_properties.disable_convection) {
                    EvalConvection(u, grad_u, conv);
                }

                double tauM = 0.0;
                double tauC_unused = 0.0;
                CalcTau(T, u, nu, idata.vms_config.Ci, dt, tauM, tauC_unused);

                RM.SetSize(vdim);
                RM = 0.0;

                for (int c = 0; c < vdim; c++) {
                    RM(c) += alpha * u(c);
                    RM(c) += -(4.0 * u_nm1(c) - u_nm2(c)) / (2.0 * dt);

                    if (!idata.flow_properties.disable_convection) {
                        RM(c) += conv(c);
                        RM(c) += u(c) * div_u;
                    }

                    for (int b = 0; b < ndof_p; b++) {
                        RM(c) += ep(b) * dNp(b, c);
                    }

                    RM(c) += -f(c);
                }

                // Apu from PSPG:
                //
                // tauM (dR_M/du, grad q)
                for (int a = 0; a < ndof_p; a++) {
                    for (int b = 0; b < ndof_u; b++) {
                        double u_dot_grad_Nb = 0.0;
                        for (int j = 0; j < dim; j++) {
                            u_dot_grad_Nb += u(j) * dNu(b, j);
                        }

                        for (int c = 0; c < vdim; c++) {
                            for (int k = 0; k < vdim; k++) {
                                const int ib = VDofIndex(ndof_u, vdim, b, k, ordering);

                                double dRM_du = 0.0;

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

                // Auu and Aup intentionally remain zero for PSPG-only integrator.
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
                // ComputeTau(u, h, dt, nu, tauM, tauC);
                CalcTau(T, u, nu, idata.vms_config.Ci, dt, tauM, tauC);


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
                // ComputeTau(u, h, dt, nu, tauM, tauC);
                CalcTau(T, u, nu, idata.vms_config.Ci, dt, tauM, tauC);

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

    // class NSEProjectionVelocityPredictorInteg : public mfem::NonlinearFormIntegrator {
    // private:
    //     const InputData& idata;
    //     const TimeLevelFields& tlf;
    //     const ProjectionScheme scheme;
    //     const int vdim;
    //     const mfem::Ordering::Type ordering;
    //     mfem::VectorCoefficient* f_coeff = nullptr;
    //
    // public:
    //     NSEProjectionVelocityPredictorInteg(
    //         const InputData& idata,
    //         const TimeLevelFields& tlf,
    //         const ProjectionScheme scheme,
    //         const int vdim,
    //         const mfem::Ordering::Type ordering,
    //         mfem::VectorCoefficient* f_coeff = nullptr)
    //         : idata(idata),
    //           tlf(tlf),
    //           scheme(scheme),
    //           vdim(vdim),
    //           ordering(ordering),
    //           f_coeff(f_coeff) {
    //     }
    //
    //     void AssembleElementVector(const mfem::FiniteElement& el,
    //                                mfem::ElementTransformation& T,
    //                                const mfem::Vector& elfun,
    //                                mfem::Vector& elvec) override {
    //         const int ndof = el.GetDof();
    //         const int dim = T.GetSpaceDim();
    //
    //         MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");
    //
    //         elvec.SetSize(vdim * ndof);
    //         elvec = 0.0;
    //
    //         mfem::Vector N(ndof);
    //         mfem::DenseMatrix dN(ndof, dim);
    //
    //         mfem::Vector u(vdim);
    //         mfem::DenseMatrix grad_u(vdim, dim);
    //
    //         mfem::Vector u_n(vdim), u_nm1(vdim), u_hist(vdim), f(vdim);
    //
    //         const int e = T.ElementNo;
    //         const mfem::IntegrationRule* ir = &tlf.femach.qspace->GetIntRule(e);
    //
    //         const double dt = tlf.GetTimeStep();
    //         const double ctime = tlf.GetTime();
    //         const double nu = idata.flow_properties.nu;
    //
    //         const ProjectionCoefficients pc = GetProjectionCoefficients(scheme, dt);
    //
    //         for (int iq = 0; iq < ir->GetNPoints(); iq++) {
    //             const mfem::IntegrationPoint& ip = ir->IntPoint(iq);
    //             T.SetIntPoint(&ip);
    //
    //             el.CalcShape(ip, N);
    //             el.CalcPhysDShape(T, dN);
    //
    //             const double wdet = ip.weight * T.Weight();
    //
    //             EvalVectorAtIP(elfun, ndof, vdim, ordering, N, u);
    //             EvalVectorGradAtIP(elfun, ndof, vdim, ordering, dN, grad_u);
    //
    //             tlf.prev_1.u.GetVectorValue(T, ip, u_n);
    //             tlf.prev_2.u.GetVectorValue(T, ip, u_nm1);
    //
    //             u_hist.SetSize(vdim);
    //             for (int c = 0; c < vdim; c++) {
    //                 u_hist(c) = pc.beta0 * u_n(c) + pc.beta1 * u_nm1(c);
    //             }
    //
    //             double p_hat = 0.0;
    //
    //             // Real p_hat evaluation.
    //             {
    //                 double p_n = 0.0;
    //                 double p_nm1 = 0.0;
    //
    //                 p_n = tlf.prev_1.p.GetValue(T, ip);
    //                 p_nm1 = tlf.prev_2.p.GetValue(T, ip);
    //
    //                 p_hat = pc.p0 * p_n + pc.p1 * p_nm1;
    //             }
    //
    //             f.SetSize(vdim);
    //             f = 0.0;
    //             if (f_coeff) {
    //                 f_coeff->SetTime(ctime);
    //                 f_coeff->Eval(f, T, ip);
    //             }
    //
    //             for (int a = 0; a < ndof; a++) {
    //                 for (int c = 0; c < vdim; c++) {
    //                     const int ia = VDofIndex(ndof, vdim, a, c, ordering);
    //
    //                     // alpha (u_star, v)
    //                     elvec(ia) += pc.alpha * u(c) * N(a) * wdet;
    //
    //                     // - (u_star tensor u_star, grad v)
    //                     if (!idata.flow_properties.disable_convection) {
    //                         for (int j = 0; j < dim; j++) {
    //                             elvec(ia) += -u(c) * u(j) * dN(a, j) * wdet;
    //                         }
    //                     }
    //
    //                     // nu (grad u_star, grad v)
    //                     for (int j = 0; j < dim; j++) {
    //                         elvec(ia) += nu * grad_u(c, j) * dN(a, j) * wdet;
    //                     }
    //
    //                     // - (p_hat, div v)
    //                     elvec(ia) += -p_hat * dN(a, c) * wdet;
    //
    //                     // - (u_history, v)
    //                     elvec(ia) += -u_hist(c) * N(a) * wdet;
    //
    //                     // - (f, v)
    //                     elvec(ia) += -f(c) * N(a) * wdet;
    //                 }
    //             }
    //         }
    //     }
    //
    //     void AssembleElementGrad(const mfem::FiniteElement& el,
    //                              mfem::ElementTransformation& T,
    //                              const mfem::Vector& elfun,
    //                              mfem::DenseMatrix& elmat) override {
    //         const int ndof = el.GetDof();
    //         const int dim = T.GetSpaceDim();
    //
    //         MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");
    //
    //         elmat.SetSize(vdim * ndof, vdim * ndof);
    //         elmat = 0.0;
    //
    //         mfem::Vector N(ndof);
    //         mfem::DenseMatrix dN(ndof, dim);
    //
    //         mfem::Vector u(vdim);
    //         mfem::DenseMatrix grad_u(vdim, dim);
    //
    //         const int e = T.ElementNo;
    //         const mfem::IntegrationRule* ir = &tlf.femach.qspace->GetIntRule(e);
    //
    //         const double dt = tlf.GetTimeStep();
    //         const double nu = idata.flow_properties.nu;
    //
    //         const ProjectionCoefficients pc = GetProjectionCoefficients(scheme, dt);
    //
    //         for (int iq = 0; iq < ir->GetNPoints(); iq++) {
    //             const mfem::IntegrationPoint& ip = ir->IntPoint(iq);
    //             T.SetIntPoint(&ip);
    //
    //             el.CalcShape(ip, N);
    //             el.CalcPhysDShape(T, dN);
    //
    //             const double wdet = ip.weight * T.Weight();
    //
    //             EvalVectorAtIP(elfun, ndof, vdim, ordering, N, u);
    //             EvalVectorGradAtIP(elfun, ndof, vdim, ordering, dN, grad_u);
    //
    //             for (int a = 0; a < ndof; a++) {
    //                 for (int b = 0; b < ndof; b++) {
    //                     const double mass_ab = pc.alpha * N(a) * N(b) * wdet;
    //
    //                     double diff_ab = 0.0;
    //                     for (int j = 0; j < dim; j++) {
    //                         diff_ab += nu * dN(a, j) * dN(b, j) * wdet;
    //                     }
    //
    //                     for (int c = 0; c < vdim; c++) {
    //                         const int ia = VDofIndex(ndof, vdim, a, c, ordering);
    //                         const int ib = VDofIndex(ndof, vdim, b, c, ordering);
    //
    //                         elmat(ia, ib) += mass_ab + diff_ab;
    //                     }
    //
    //                     if (!idata.flow_properties.disable_convection) {
    //                         double u_dot_grad_Na = 0.0;
    //                         for (int j = 0; j < dim; j++) {
    //                             u_dot_grad_Na += u(j) * dN(a, j);
    //                         }
    //
    //                         for (int c = 0; c < vdim; c++) {
    //                             const int ia = VDofIndex(ndof, vdim, a, c, ordering);
    //
    //                             for (int k = 0; k < vdim; k++) {
    //                                 const int ib =
    //                                     VDofIndex(ndof, vdim, b, k, ordering);
    //
    //                                 double conv_jac = 0.0;
    //
    //                                 // - (du tensor u, grad v)
    //                                 if (c == k) {
    //                                     conv_jac += -N(b) * u_dot_grad_Na;
    //                                 }
    //
    //                                 // - (u tensor du, grad v)
    //                                 conv_jac += -u(c) * N(b) * dN(a, k);
    //
    //                                 elmat(ia, ib) += conv_jac * wdet;
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // };
    //
    // class NSEProjectionOutletConvectiveFluxInteg : public mfem::NonlinearFormIntegrator {
    // private:
    //     const InputData& idata;
    //     const int vdim;
    //     const mfem::Ordering::Type ordering;
    //
    // public:
    //     NSEProjectionOutletConvectiveFluxInteg(
    //         const InputData& idata,
    //         const int vdim,
    //         const mfem::Ordering::Type ordering)
    //         : idata(idata), vdim(vdim), ordering(ordering) {
    //     }
    //
    //     void AssembleFaceVector(const mfem::FiniteElement& el1,
    //                             const mfem::FiniteElement& el2,
    //                             mfem::FaceElementTransformations& Tr,
    //                             const mfem::Vector& elfun,
    //                             mfem::Vector& elvec) override {
    //         const int ndof = el1.GetDof();
    //         const int dim = Tr.Elem1->GetSpaceDim();
    //
    //         MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");
    //
    //         elvec.SetSize(vdim * ndof);
    //         elvec = 0.0;
    //
    //         if (idata.flow_properties.disable_convection) {
    //             return;
    //         }
    //
    //         mfem::Vector N(ndof);
    //         mfem::Vector u(vdim);
    //         mfem::Vector nor(dim);
    //
    //         const int order = 2 * el1.GetOrder();
    //         const mfem::IntegrationRule* ir =
    //             &mfem::IntRules.Get(Tr.GetGeometryType(), order);
    //
    //         for (int iq = 0; iq < ir->GetNPoints(); iq++) {
    //             const mfem::IntegrationPoint& ip_face = ir->IntPoint(iq);
    //
    //             Tr.Face->SetIntPoint(&ip_face);
    //
    //             mfem::IntegrationPoint ip_el;
    //             Tr.Loc1.Transform(ip_face, ip_el);
    //             Tr.Elem1->SetIntPoint(&ip_el);
    //
    //             el1.CalcShape(ip_el, N);
    //             EvalVectorAtIP(elfun, ndof, vdim, ordering, N, u);
    //
    //             mfem::CalcOrtho(Tr.Face->Jacobian(), nor);
    //
    //             double u_dot_n = 0.0;
    //             for (int j = 0; j < dim; j++) {
    //                 u_dot_n += u(j) * nor(j);
    //             }
    //
    //             const double w = ip_face.weight;
    //
    //             for (int a = 0; a < ndof; a++) {
    //                 for (int c = 0; c < vdim; c++) {
    //                     const int ia = VDofIndex(ndof, vdim, a, c, ordering);
    //
    //                     // < u (u . n), v >
    //                     elvec(ia) += N(a) * u(c) * u_dot_n * w;
    //                 }
    //             }
    //         }
    //     }
    //
    //     void AssembleFaceGrad(const mfem::FiniteElement& el1,
    //                           const mfem::FiniteElement& el2,
    //                           mfem::FaceElementTransformations& Tr,
    //                           const mfem::Vector& elfun,
    //                           mfem::DenseMatrix& elmat) override {
    //         const int ndof = el1.GetDof();
    //         const int dim = Tr.Elem1->GetSpaceDim();
    //
    //         MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");
    //
    //         elmat.SetSize(vdim * ndof, vdim * ndof);
    //         elmat = 0.0;
    //
    //         if (idata.flow_properties.disable_convection) {
    //             return;
    //         }
    //
    //         mfem::Vector N(ndof);
    //         mfem::Vector u(vdim);
    //         mfem::Vector nor(dim);
    //
    //         const int order = 2 * el1.GetOrder();
    //         const mfem::IntegrationRule* ir =
    //             &mfem::IntRules.Get(Tr.GetGeometryType(), order);
    //
    //         for (int iq = 0; iq < ir->GetNPoints(); iq++) {
    //             const mfem::IntegrationPoint& ip_face = ir->IntPoint(iq);
    //
    //             Tr.Face->SetIntPoint(&ip_face);
    //
    //             mfem::IntegrationPoint ip_el;
    //             Tr.Loc1.Transform(ip_face, ip_el);
    //             Tr.Elem1->SetIntPoint(&ip_el);
    //
    //             el1.CalcShape(ip_el, N);
    //             EvalVectorAtIP(elfun, ndof, vdim, ordering, N, u);
    //
    //             mfem::CalcOrtho(Tr.Face->Jacobian(), nor);
    //
    //             double u_dot_n = 0.0;
    //             for (int j = 0; j < dim; j++) {
    //                 u_dot_n += u(j) * nor(j);
    //             }
    //
    //             const double w = ip_face.weight;
    //
    //             for (int a = 0; a < ndof; a++) {
    //                 for (int b = 0; b < ndof; b++) {
    //                     for (int c = 0; c < vdim; c++) {
    //                         const int ia = VDofIndex(ndof, vdim, a, c, ordering);
    //
    //                         for (int k = 0; k < vdim; k++) {
    //                             const int ib =
    //                                 VDofIndex(ndof, vdim, b, k, ordering);
    //
    //                             double jac_ck = 0.0;
    //
    //                             if (c == k) {
    //                                 jac_ck += u_dot_n;
    //                             }
    //
    //                             jac_ck += u(c) * nor(k);
    //
    //                             elmat(ia, ib) += N(a) * N(b) * jac_ck * w;
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // };
    //
    // class NSEPressureIncrementRHSInteg : public mfem::LinearFormIntegrator {
    // private:
    //     const TimeLevelFields& tlf;
    //     const mfem::ParGridFunction& u_star;
    //     const double alpha;
    //
    // public:
    //     NSEPressureIncrementRHSInteg(const TimeLevelFields& tlf,
    //                                  const mfem::ParGridFunction& u_star,
    //                                  const double alpha)
    //         : tlf(tlf), u_star(u_star), alpha(alpha) {
    //     }
    //
    //     void AssembleRHSElementVect(const mfem::FiniteElement& el,
    //                                 mfem::ElementTransformation& T,
    //                                 mfem::Vector& elvec) override {
    //         const int ndof = el.GetDof();
    //         const int dim = T.GetSpaceDim();
    //
    //         elvec.SetSize(ndof);
    //         elvec = 0.0;
    //
    //         mfem::Vector N(ndof);
    //         mfem::DenseMatrix grad_u;
    //
    //         const int e = T.ElementNo;
    //         const mfem::IntegrationRule* ir = &tlf.femach.qspace->GetIntRule(e);
    //
    //         for (int iq = 0; iq < ir->GetNPoints(); iq++) {
    //             const mfem::IntegrationPoint& ip = ir->IntPoint(iq);
    //             T.SetIntPoint(&ip);
    //
    //             el.CalcShape(ip, N);
    //
    //             const double wdet = ip.weight * T.Weight();
    //
    //             u_star.GetVectorGradient(T, grad_u);
    //
    //             double div_u = 0.0;
    //             for (int c = 0; c < dim; c++) {
    //                 div_u += grad_u(c, c);
    //             }
    //
    //             for (int a = 0; a < ndof; a++) {
    //                 elvec(a) += -alpha * N(a) * div_u * wdet;
    //             }
    //         }
    //     }
    // };
    //
    // class NSEVelocityCorrectionRHSInteg : public mfem::LinearFormIntegrator {
    // private:
    //     const TimeLevelFields& tlf;
    //     const mfem::ParGridFunction& u_star;
    //     const mfem::ParGridFunction& phi;
    //     const double inv_alpha;
    //     const int vdim;
    //     const mfem::Ordering::Type ordering;
    //
    // public:
    //     NSEVelocityCorrectionRHSInteg(
    //         const TimeLevelFields& tlf,
    //         const mfem::ParGridFunction& u_star,
    //         const mfem::ParGridFunction& phi,
    //         const double alpha,
    //         const int vdim,
    //         const mfem::Ordering::Type ordering)
    //         : tlf(tlf),
    //           u_star(u_star),
    //           phi(phi),
    //           inv_alpha(1.0 / alpha),
    //           vdim(vdim),
    //           ordering(ordering) {
    //     }
    //
    //     void AssembleRHSElementVect(const mfem::FiniteElement& el,
    //                                 mfem::ElementTransformation& T,
    //                                 mfem::Vector& elvec) override {
    //         const int ndof = el.GetDof();
    //         const int dim = T.GetSpaceDim();
    //
    //         MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");
    //
    //         elvec.SetSize(vdim * ndof);
    //         elvec = 0.0;
    //
    //         mfem::Vector N(ndof);
    //         mfem::Vector u(vdim);
    //         mfem::Vector grad_phi(dim);
    //
    //         const int e = T.ElementNo;
    //         const mfem::IntegrationRule* ir = &tlf.femach.qspace->GetIntRule(e);
    //
    //         for (int iq = 0; iq < ir->GetNPoints(); iq++) {
    //             const mfem::IntegrationPoint& ip = ir->IntPoint(iq);
    //             T.SetIntPoint(&ip);
    //
    //             el.CalcShape(ip, N);
    //
    //             const double wdet = ip.weight * T.Weight();
    //
    //             u_star.GetVectorValue(T, ip, u);
    //             phi.GetGradient(T, grad_phi);
    //
    //             for (int a = 0; a < ndof; a++) {
    //                 for (int c = 0; c < vdim; c++) {
    //                     const int ia = VDofIndex(ndof, vdim, a, c, ordering);
    //
    //                     elvec(ia) +=
    //                         N(a) * (u(c) - inv_alpha * grad_phi(c)) * wdet;
    //                 }
    //             }
    //         }
    //     }
    // };

    class NSEProjVMSIntegMomentumConvForm : public mfem::NonlinearFormIntegrator {
    private:
        const InputData& idata;
        const TimeLevelFields& tlf;
        const int vdim;
        const mfem::Ordering::Type ordering;
        mfem::VectorCoefficient* f_coeff = nullptr;

    public:
        NSEProjVMSIntegMomentumConvForm(
            const InputData& idata,
            const TimeLevelFields& tlf,
            const int vdim,
            const mfem::Ordering::Type ordering,
            mfem::VectorCoefficient* f_coeff = nullptr)
            : idata(idata),
              tlf(tlf),
              vdim(vdim),
              ordering(ordering),
              f_coeff(f_coeff) {
        }

    private:

        void EvalPressureHatAndGradAtIP(mfem::ElementTransformation& T,
                                        const mfem::IntegrationPoint& ip,
                                        const ProjectionCoefficients& pc,
                                        double& p_hat,
                                        mfem::Vector& grad_p_hat) const {
            // Delegates to the shared free function (InputData.h /
            // Integrators.h) so that the momentum predictor, VUE RHS, and
            // PPE RHS integrators (conv-form and div-form alike) all
            // compute p_h* identically.
            nse::EvalPressureHatAndGradAtIP(tlf, T, ip, pc, vdim, p_hat, grad_p_hat);
        }

    public:
        void AssembleElementVector(const mfem::FiniteElement& el,
                                   mfem::ElementTransformation& T,
                                   const mfem::Vector& elfun,
                                   mfem::Vector& elvec) override {
            const int ndof = el.GetDof();
            const int dim = T.GetSpaceDim();

            MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");

            elvec.SetSize(vdim * ndof);
            elvec = 0.0;

            mfem::Vector N(ndof);
            mfem::DenseMatrix dN(ndof, dim);

            mfem::Vector u(vdim);
            mfem::DenseMatrix grad_u(vdim, dim);

            mfem::Vector u_n(vdim), u_nm1(vdim), u_hist(vdim), f(vdim);

            mfem::Vector conv(vdim);
            mfem::Vector grad_p_hat(vdim);
            mfem::Vector lap_u(vdim);
            mfem::Vector lap_N(ndof);
            mfem::Vector ns_res(vdim);

            const int e = T.ElementNo;
            const mfem::IntegrationRule* ir = &tlf.femach.qspace->GetIntRule(e);

            const double dt = tlf.GetTimeStep();
            const double ctime = tlf.GetTime();
            const double nu = idata.flow_properties.nu;

            const ProjectionCoefficients pc = GetProjectionCoefficients(idata.projection_config.scheme);

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint& ip = ir->IntPoint(iq);
                T.SetIntPoint(&ip);

                el.CalcShape(ip, N);
                el.CalcPhysDShape(T, dN);
                el.CalcPhysLaplacian(T, lap_N);

                const double wdet = ip.weight * T.Weight();

                EvalVectorAtIP(elfun, ndof, vdim, ordering, N, u);
                EvalVectorGradAtIP(elfun, ndof, vdim, ordering, dN, grad_u);
                EvalVectorLaplacianAtIP(elfun, ndof, vdim, ordering, lap_N, lap_u);

                tlf.prev_1.u.GetVectorValue(T, ip, u_n);
                tlf.prev_2.u.GetVectorValue(T, ip, u_nm1);

                u_hist.SetSize(vdim);
                for (int c = 0; c < vdim; c++) {
                    u_hist(c) = pc.beta1 * u_n(c) + pc.beta2 * u_nm1(c);
                }

                double p_hat = 0.0;
                EvalPressureHatAndGradAtIP(T, ip, pc, p_hat, grad_p_hat);

                f.SetSize(vdim);
                f = 0.0;
                if (f_coeff) {
                    f_coeff->SetTime(ctime);
                    f_coeff->Eval(f, T, ip);
                }

                conv.SetSize(vdim);
                conv = 0.0;

                for (int i = 0; i < vdim; i++) {
                    for (int j = 0; j < dim; j++) {
                        conv(i) += u(j) * grad_u(i, j);
                    }
                }

                double tauM, tauC;
                CalcTau(T, u, nu, idata.vms_config.Ci, dt, tauM, tauC);

                // -----------------------------------------------------------------
                // Strong residual used by the VMS fine scale:
                //
                // old:
                // ns_res(i) =
                //   BDF(u)_i + convection_i + gradp_i
                //   - nu * laplacian(u_i) - forcing_i
                // -----------------------------------------------------------------
                ns_res.SetSize(vdim);
                for (int i = 0; i < vdim; i++) {
                    ns_res(i) =
                        (pc.beta0 * u(i) + u_hist(i)) / dt
                        + conv(i)
                        + grad_p_hat(i)
                        - nu * lap_u(i)
                        - f(i);
                }

                for (int a = 0; a < ndof; a++) {
                    double crossTermVelocityPart = 0.0;
                    for (int i = 0; i < dim; i++) {
                        crossTermVelocityPart += dN(a, i) * u(i);
                    }

                    double crossTermFineScalePart = 0.0;
                    for (int i = 0; i < dim; i++) {
                        crossTermFineScalePart += dN(a, i) * tauM * ns_res(i);
                    }

                    for (int c = 0; c < vdim; c++) {
                        const int ia = VDofIndex(ndof, vdim, a, c, ordering);

                        // ---------------------------------------------------------
                        // Normal Galerkin Navier--Stokes terms.
                        // This matches your old normal_NS contribution:
                        //
                        //   N(a) * BDF(u)
                        // + N(a) * convection
                        // + N(a) * grad_p_hat(c)
                        // + nu * grad_u : grad_v
                        // - N(a) * f
                        // ---------------------------------------------------------

                        elvec(ia) += N(a) * ((pc.beta0 * u(c) + u_hist(c)) / dt) * wdet;

                        elvec(ia) += N(a) * conv(c) * wdet;

                        elvec(ia) += N(a) * grad_p_hat(c) * wdet;

                        for (int j = 0; j < dim; j++) {
                            elvec(ia) += nu * dN(a, j) * grad_u(c, j) * wdet;
                        }

                        elvec(ia) += -N(a) * f(c) * wdet;

                        // ---------------------------------------------------------
                        // VMS cross term 1:
                        //
                        //   + tauM * (u . grad w) * Res_c
                        // ---------------------------------------------------------

                        elvec(ia) +=
                            tauM * crossTermVelocityPart * ns_res(c) * wdet;

                        // ---------------------------------------------------------
                        // VMS cross term 2:
                        //
                        //   - w * grad(u_c)_j * tauM * Res_j
                        // ---------------------------------------------------------

                        for (int j = 0; j < dim; j++) {
                            elvec(ia) +=
                                -N(a) * grad_u(c, j) * tauM * ns_res(j) * wdet;
                        }

                        // ---------------------------------------------------------
                        // Reynolds-stress / fine-scale term:
                        //
                        //   - grad(w)_j * tauM Res_c * tauM Res_j
                        // ---------------------------------------------------------

                        for (int j = 0; j < dim; j++) {
                            elvec(ia) +=
                                -dN(a, j)
                                * (tauM * ns_res(c))
                                * (tauM * ns_res(j))
                                * wdet;
                        }
                    }
                }
            }
        }

        void AssembleElementGrad(const mfem::FiniteElement& el,
                                 mfem::ElementTransformation& T,
                                 const mfem::Vector& elfun,
                                 mfem::DenseMatrix& elmat) override {
            const int ndof = el.GetDof();
            const int dim = T.GetSpaceDim();

            MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");

            elmat.SetSize(vdim * ndof, vdim * ndof);
            elmat = 0.0;

            mfem::Vector N(ndof);
            mfem::DenseMatrix dN(ndof, dim);

            mfem::Vector u(vdim);
            mfem::DenseMatrix grad_u(vdim, dim);

            mfem::Vector u_n(vdim), u_nm1(vdim), u_hist(vdim), f(vdim);

            mfem::Vector conv(vdim);
            mfem::Vector grad_p_hat(vdim);
            mfem::Vector lap_u(vdim);
            mfem::Vector lap_N(ndof);
            mfem::Vector ns_res(vdim);

            const int e = T.ElementNo;
            const mfem::IntegrationRule* ir = &tlf.femach.qspace->GetIntRule(e);

            const double dt = tlf.GetTimeStep();
            const double ctime = tlf.GetTime();
            const double nu = idata.flow_properties.nu;

            const ProjectionCoefficients pc = GetProjectionCoefficients(idata.projection_config.scheme);
            const double alpha = pc.beta0 / dt;

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint& ip = ir->IntPoint(iq);
                T.SetIntPoint(&ip);

                el.CalcShape(ip, N);
                el.CalcPhysDShape(T, dN);
                el.CalcPhysLaplacian(T, lap_N);

                const double wdet = ip.weight * T.Weight();

                EvalVectorAtIP(elfun, ndof, vdim, ordering, N, u);
                EvalVectorGradAtIP(elfun, ndof, vdim, ordering, dN, grad_u);
                EvalVectorLaplacianAtIP(elfun, ndof, vdim, ordering, lap_N, lap_u);

                tlf.prev_1.u.GetVectorValue(T, ip, u_n);
                tlf.prev_2.u.GetVectorValue(T, ip, u_nm1);

                u_hist.SetSize(vdim);
                for (int c = 0; c < vdim; c++) {
                    u_hist(c) = pc.beta1 * u_n(c) + pc.beta2 * u_nm1(c);
                }

                double p_hat = 0.0;
                EvalPressureHatAndGradAtIP(T, ip, pc, p_hat, grad_p_hat);

                f.SetSize(vdim);
                f = 0.0;
                if (f_coeff) {
                    f_coeff->SetTime(ctime);
                    f_coeff->Eval(f, T, ip);
                }

                conv.SetSize(vdim);
                conv = 0.0;

                for (int i = 0; i < vdim; i++) {
                    for (int j = 0; j < dim; j++) {
                        conv(i) += u(j) * grad_u(i, j);
                    }
                }

                double tauM, tauC;
                CalcTau(T, u, nu, idata.vms_config.Ci, dt, tauM, tauC);

                ns_res.SetSize(vdim);
                for (int i = 0; i < vdim; i++) {
                    ns_res(i) =
                        (pc.beta0 * u(i) + u_hist(i)) / dt
                        + conv(i)
                        + grad_p_hat(i)
                        - nu * lap_u(i)
                        - f(i);
                }

                for (int a = 0; a < ndof; a++) {
                    double crossTermVelocityPart = 0.0;
                    for (int i = 0; i < dim; i++) {
                        crossTermVelocityPart += dN(a, i) * u(i);
                    }

                    double crossTermFineScalePart = 0.0;
                    for (int i = 0; i < dim; i++) {
                        crossTermFineScalePart += dN(a, i) * tauM * ns_res(i);
                    }

                    for (int b = 0; b < ndof; b++) {
                        double conv_b = 0.0;
                        for (int i = 0; i < dim; i++) {
                            conv_b += dN(b, i) * u(i);
                        }

                        double gradNa_dot_gradNb = 0.0;
                        for (int j = 0; j < dim; j++) {
                            gradNa_dot_gradNb += dN(a, j) * dN(b, j);
                        }

                        // old:
                        //
                        // diff_J += nu * fe.d2N(b,j,j)
                        //
                        // Here lap_N(b) = sum_j d^2 N_b / dx_j^2.
                        const double diff_J = nu * lap_N(b);

                        // ---------------------------------------------------------
                        // 1. Normal Navier--Stokes Jacobian.
                        // ---------------------------------------------------------

                        for (int c = 0; c < vdim; c++) {
                            const int ia = VDofIndex(ndof, vdim, a, c, ordering);
                            const int ib = VDofIndex(ndof, vdim, b, c, ordering);

                            elmat(ia, ib) += alpha * N(a) * N(b) * wdet;

                            // diffusion
                            elmat(ia, ib) += nu * gradNa_dot_gradNb * wdet;
                        }

                        for (int c = 0; c < vdim; c++) {
                            const int ia = VDofIndex(ndof, vdim, a, c, ordering);

                            for (int k = 0; k < vdim; k++) {
                                const int ib = VDofIndex(ndof, vdim, b, k, ordering);

                                double val = 0.0;

                                // old diagonal term:
                                //
                                //   N(a) * conv_b
                                //
                                // from (w, u . grad(delta u_c))
                                if (c == k) {
                                    val += N(a) * conv_b;
                                }

                                // old off-diagonal term:
                                //
                                //   N(a) * grad_u(c,k) * N(b)
                                //
                                // from (w, delta u . grad u_c)
                                val += N(a) * grad_u(c, k) * N(b);

                                elmat(ia, ib) += val * wdet;
                            }
                        }

                        // ---------------------------------------------------------
                        // 2. Cross term 1 Jacobian.
                        //
                        // Residual contribution:
                        //
                        //   + tauM * (u . grad w) * Res_c
                        //
                        // old code block:
                        //
                        //   crossTermVelocityPart * tauM *
                        //       (alpha N_b + conv_b - diff_J)
                        //
                        //   crossTermVelocityPart * tauM *
                        //       grad_u(c,k) * N_b
                        //
                        //   tauM * dN(a,k) * ns_res(c) * N_b
                        // ---------------------------------------------------------

                        for (int c = 0; c < vdim; c++) {
                            const int ia = VDofIndex(ndof, vdim, a, c, ordering);

                            for (int k = 0; k < vdim; k++) {
                                const int ib = VDofIndex(ndof, vdim, b, k, ordering);

                                double val = 0.0;

                                if (c == k) {
                                    val +=
                                        crossTermVelocityPart * tauM *
                                        (alpha * N(b) + conv_b - diff_J);
                                }

                                val +=
                                    crossTermVelocityPart * tauM *
                                    grad_u(c, k) * N(b);

                                val +=
                                    tauM * dN(a, k) *
                                    ns_res(c) * N(b);

                                elmat(ia, ib) += val * wdet;
                            }
                        }

                        // ---------------------------------------------------------
                        // 3. Cross term 2 Jacobian.
                        //
                        // Residual contribution:
                        //
                        //   - N(a) * grad_u(c,j) * tauM * Res_j
                        //
                        // This follows your old signs.
                        // ---------------------------------------------------------

                        for (int c = 0; c < vdim; c++) {
                            const int ia = VDofIndex(ndof, vdim, a, c, ordering);

                            for (int k = 0; k < vdim; k++) {
                                const int ib = VDofIndex(ndof, vdim, b, k, ordering);

                                double val = 0.0;

                                // old:
                                //
                                // - grad_u(c,j) * N(a) * tauM * N(b) * grad_u(j,k)
                                for (int j = 0; j < dim; j++) {
                                    val +=
                                        -grad_u(c, j) *
                                        N(a) *
                                        tauM *
                                        N(b) *
                                        grad_u(j, k);
                                }

                                // old:
                                //
                                // - grad_u(c,k) * N(a) * tauM *
                                //      (alpha N_b + conv_b - diff_J)
                                val +=
                                    -grad_u(c, k) *
                                    N(a) *
                                    tauM *
                                    (alpha * N(b) + conv_b - diff_J);

                                // old diagonal contribution:
                                //
                                // - N(a) * tauM * ns_res(k) * dN(b,k)
                                //
                                // In the old code this is added to the diagonal
                                // equation/component c through column component c.
                                if (c == k) {
                                    for (int j = 0; j < dim; j++) {
                                        val +=
                                            -N(a) *
                                            tauM *
                                            ns_res(j) *
                                            dN(b, j);
                                    }
                                }

                                elmat(ia, ib) += val * wdet;
                            }
                        }

                        // ---------------------------------------------------------
                        // 4. Reynolds-stress Jacobian.
                        //
                        // Residual contribution:
                        //
                        //   - dN(a,j) * tauM Res_c * tauM Res_j
                        //
                        // This follows the four groups in your old code.
                        // ---------------------------------------------------------

                        for (int c = 0; c < vdim; c++) {
                            const int ia = VDofIndex(ndof, vdim, a, c, ordering);

                            for (int k = 0; k < vdim; k++) {
                                const int ib = VDofIndex(ndof, vdim, b, k, ordering);

                                double val = 0.0;

                                // old first group:
                                //
                                // -crossTermFineScalePart * tauM *
                                //      (alpha N_b + conv_b - diff_J)
                                //
                                // diagonal with respect to component c
                                if (c == k) {
                                    val +=
                                        -crossTermFineScalePart *
                                        tauM *
                                        (alpha * N(b) + conv_b - diff_J);
                                }

                                // old second group:
                                //
                                // -crossTermFineScalePart * tauM *
                                //      grad_u(c,k) * N_b
                                val +=
                                    -crossTermFineScalePart *
                                    tauM *
                                    grad_u(c, k) *
                                    N(b);

                                // old third group:
                                //
                                // -tauM * ns_res(c) * dN(a,j) *
                                //      tauM * N(b) * grad_u(j,k)
                                for (int j = 0; j < dim; j++) {
                                    val +=
                                        -tauM *
                                        ns_res(c) *
                                        dN(a, j) *
                                        tauM *
                                        N(b) *
                                        grad_u(j, k);
                                }

                                // old fourth group:
                                //
                                // -tauM * ns_res(c) * dN(a,k) * tauM *
                                //      (alpha N_b + conv_b - diff_J)
                                val +=
                                    -tauM *
                                    ns_res(c) *
                                    dN(a, k) *
                                    tauM *
                                    (alpha * N(b) + conv_b - diff_J);

                                elmat(ia, ib) += val * wdet;
                            }
                        }
                    }
                }
            }
        }
    };

    class NSEProjVMSIntegVUERHSConvForm : public mfem::LinearFormIntegrator {
    private:
        const InputData& idata;
        const TimeLevelFields& tlf;
        const int vdim;
        const mfem::Ordering::Type ordering;
        mfem::VectorCoefficient* f_coeff = nullptr;

    public:
        NSEProjVMSIntegVUERHSConvForm(
            const InputData& idata,
            const TimeLevelFields& tlf,
            const int vdim,
            const mfem::Ordering::Type ordering,
            mfem::VectorCoefficient* f_coeff = nullptr)
            : idata(idata),
              tlf(tlf),
              vdim(vdim),
              ordering(ordering),
              f_coeff(f_coeff) {
        }

        void AssembleRHSElementVect(const mfem::FiniteElement& el,
                                    mfem::ElementTransformation& T,
                                    mfem::Vector& elvec) override {
            const int ndof = el.GetDof();
            const int dim = T.GetSpaceDim();

            MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");

            elvec.SetSize(vdim * ndof);
            elvec = 0.0;

            mfem::Vector N(ndof);

            mfem::Vector u(vdim);
            mfem::Vector u_n(vdim), u_nm1(vdim), u_hist(vdim);
            mfem::Vector f(vdim);

            mfem::DenseMatrix grad_u(vdim, dim);

            mfem::Vector conv(vdim);
            mfem::Vector lap_u(vdim);
            mfem::Vector grad_p_curr(vdim);
            mfem::Vector grad_p_hat(vdim);
            mfem::Vector ns_res(vdim);

            const int e = T.ElementNo;
            const mfem::IntegrationRule* ir = &tlf.femach.qspace->GetIntRule(e);

            const double dt = tlf.GetTimeStep();
            const double ctime = tlf.GetTime();
            const double nu = idata.flow_properties.nu;

            const ProjectionCoefficients pc = GetProjectionCoefficients(idata.projection_config.scheme);

            const double sigma = pc.beta0 / dt;
            const double sigmainv = 1.0 / sigma;

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint& ip = ir->IntPoint(iq);
                T.SetIntPoint(&ip);

                el.CalcShape(ip, N);

                const double wdet = ip.weight * T.Weight();

                // Current velocity, old velocities.
                tlf.current.u.GetVectorValue(T, ip, u);
                tlf.prev_1.u.GetVectorValue(T, ip, u_n);
                tlf.prev_2.u.GetVectorValue(T, ip, u_nm1);

                // Current velocity gradient.
                tlf.current.u.GetVectorGradient(T, grad_u);

                // Current pressure gradient, and p_h* = p0*p^n + p1*p^{n-1}
                // (same helper the momentum predictor uses, so all three
                // projection steps agree on p_h*).
                tlf.current.p.GetGradient(T, grad_p_curr);

                double p_hat = 0.0;
                EvalPressureHatAndGradAtIP(tlf, T, ip, pc, vdim, p_hat, grad_p_hat);

                EvalGridFunctionVectorLaplacianAtIP(tlf.current.u, T, vdim, lap_u);

                u_hist.SetSize(vdim);
                for (int c = 0; c < vdim; c++) {
                    u_hist(c) = pc.beta1 * u_n(c) + pc.beta2 * u_nm1(c);
                }

                f.SetSize(vdim);
                f = 0.0;
                if (f_coeff) {
                    f_coeff->SetTime(ctime);
                    f_coeff->Eval(f, T, ip);
                }

                conv.SetSize(vdim);
                conv = 0.0;

                for (int i = 0; i < vdim; i++) {
                    for (int j = 0; j < dim; j++) {
                        conv(i) += u(j) * grad_u(i, j);
                    }
                }

                double tauM, tauC;
                CalcTau(T, u, nu, idata.vms_config.Ci, dt, tauM, tauC);

                ns_res.SetSize(vdim);
                for (int i = 0; i < vdim; i++) {
                    ns_res(i) =
                        (pc.beta0 * u(i) + u_hist(i)) / dt
                        + conv(i)
                        - nu * lap_u(i)
                        + grad_p_hat(i)
                        - f(i);
                }

                for (int a = 0; a < ndof; a++) {
                    for (int c = 0; c < vdim; c++) {
                        const int ia = VDofIndex(ndof, vdim, a, c, ordering);

                        elvec(ia) +=
                            N(a) * (u(c) - tauM * ns_res(c)) * wdet;

                        elvec(ia) +=
                            -sigmainv * N(a) *
                            (grad_p_curr(c) - grad_p_hat(c)) *
                            wdet;
                    }
                }
            }
        }
    };

    class NSEProjVMSIntegPPERHSConvForm : public mfem::LinearFormIntegrator {
    private:
        const InputData& idata;
        const TimeLevelFields& tlf;
        const int vdim;
        mfem::VectorCoefficient* f_coeff = nullptr;

    public:
        NSEProjVMSIntegPPERHSConvForm(
            const InputData& idata,
            const TimeLevelFields& tlf,
            const int vdim,
            mfem::VectorCoefficient* f_coeff = nullptr)
            : idata(idata),
              tlf(tlf),
              vdim(vdim),
              f_coeff(f_coeff) {
        }

        void AssembleRHSElementVect(const mfem::FiniteElement& el,
                                    mfem::ElementTransformation& T,
                                    mfem::Vector& elvec) override {
            const int ndof = el.GetDof();
            const int dim = T.GetSpaceDim();

            MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");

            elvec.SetSize(ndof);
            elvec = 0.0;

            mfem::Vector N(ndof);
            mfem::DenseMatrix dN(ndof, dim);

            mfem::Vector u(vdim);
            mfem::Vector u_n(vdim), u_nm1(vdim), u_hist(vdim);
            mfem::Vector f(vdim);

            mfem::DenseMatrix grad_u(vdim, dim);

            mfem::Vector conv(vdim);
            mfem::Vector lap_u(vdim);
            mfem::Vector grad_p_hat(vdim);
            mfem::Vector ns_res(vdim);

            const int e = T.ElementNo;
            const mfem::IntegrationRule* ir = &tlf.femach.qspace->GetIntRule(e);

            const double dt = tlf.GetTimeStep();
            const double ctime = tlf.GetTime();
            const double nu = idata.flow_properties.nu;

            const ProjectionCoefficients pc = GetProjectionCoefficients(idata.projection_config.scheme);

            const double sigma = pc.beta0 / dt;

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint& ip = ir->IntPoint(iq);
                T.SetIntPoint(&ip);

                el.CalcShape(ip, N);
                el.CalcPhysDShape(T, dN);

                const double wdet = ip.weight * T.Weight();

                // Current velocity and old velocities.
                tlf.current.u.GetVectorValue(T, ip, u);
                tlf.prev_1.u.GetVectorValue(T, ip, u_n);
                tlf.prev_2.u.GetVectorValue(T, ip, u_nm1);

                // Current velocity gradient.
                tlf.current.u.GetVectorGradient(T, grad_u);

                // p_h* = p0*p^n + p1*p^{n-1} (same helper the momentum
                // predictor and VUE RHS use, so all three projection steps
                // agree on p_h*).
                double p_hat = 0.0;
                EvalPressureHatAndGradAtIP(tlf, T, ip, pc, vdim, p_hat, grad_p_hat);

                EvalGridFunctionVectorLaplacianAtIP(tlf.current.u, T, vdim, lap_u);

                u_hist.SetSize(vdim);
                for (int c = 0; c < vdim; c++) {
                    u_hist(c) = pc.beta1 * u_n(c) + pc.beta2 * u_nm1(c);
                }

                f.SetSize(vdim);
                f = 0.0;
                if (f_coeff) {
                    f_coeff->SetTime(ctime);
                    f_coeff->Eval(f, T, ip);
                }

                conv.SetSize(vdim);
                conv = 0.0;

                for (int i = 0; i < vdim; i++) {
                    for (int j = 0; j < dim; j++) {
                        conv(i) += u(j) * grad_u(i, j);
                    }
                }

                double tauM, tauC;
                CalcTau(T, u, nu, idata.vms_config.Ci, dt, tauM, tauC);

                ns_res.SetSize(vdim);
                for (int i = 0; i < vdim; i++) {
                    ns_res(i) =
                        (pc.beta0 * u(i) + u_hist(i)) / dt
                        + conv(i)
                        - nu * lap_u(i)
                        + grad_p_hat(i)
                        - f(i);
                }

                double div_u = 0.0;
                for (int k = 0; k < dim; k++) {
                    div_u += grad_u(k, k);
                }

                for (int a = 0; a < ndof; a++) {
                    // -sigma * N(a) * div(u)
                    elvec(a) += -sigma * N(a) * div_u * wdet;

                    for (int k = 0; k < dim; k++) {
                        // sigma * grad(q)_k * (-tauM * NS_residual_k)
                        elvec(a) +=
                            sigma * dN(a, k) *
                            (-tauM * ns_res(k)) *
                            wdet;

                        // grad(q)_k * grad(p_h*)_k
                        elvec(a) +=
                            dN(a, k) *
                            grad_p_hat(k) *
                            wdet;
                    }
                }
            }
        }
    };

    // =====================================================================
    // Divergence (conservative) form of the projection-method integrators.
    //
    // These mirror NSEProjVMSIntegMomentumConvForm / NSEProjVMSIntegVUERHSConvForm /
    // NSEProjVMSIntegPPERHSConvForm, but build the nonlinear convective term
    // using the divergence form
    //
    //     div(u tensor u)_c = d_j(u_c u_j) = (u . grad) u_c + u_c (div u)
    //
    // wherever the coarse-scale strong residual Res is evaluated (i.e. in
    // the VMS closure terms), and use the corresponding integrated-by-parts
    // Galerkin convective term
    //
    //     -(u tensor u, grad v)
    //
    // in place of the direct/advective (u . grad u, v) term used by the
    // *ConvForm classes. This implements the VMS-stabilized projection
    // scheme:
    //
    //   Step 1 (momentum predictor, nonlinear):
    //     sigma(u,v) - (u tensor u, grad v) + nu(grad u, grad v) + (grad p*, v)
    //     + sum_K (u tensor tauM*Res, grad v)_K
    //     + sum_K (tauM*Res tensor u, grad v)_K
    //     - sum_K (tauM^2 Res tensor Res, grad v)_K
    //     - (g, v) = 0
    //
    //   Step 2 (pressure Poisson, linear):
    //     (grad p_hat, grad q) = (grad p*, grad q) - sigma (div u, q)
    //       - sigma sum_K (tauM*Res, grad q)_K
    //
    //   Step 3 (velocity projection, linear):
    //     (u_hat, w) = (u, w) - sum_K (tauM*Res, w)_K
    //       - (1/sigma) (grad(p_hat - p*), w)
    //
    // with Res = sigma*u + (history)/dt + div(u tensor u) - nu*Delta(u)
    //            + grad(p*) - g.
    //
    // NOTE: switching the coarse-scale Galerkin convective term to the
    // integrated-by-parts (divergence) form introduces a boundary flux
    //
    //     <n . [u tensor u - nu grad u], v>_{Gamma_N}
    //
    // on any boundary where velocity is not strongly imposed (open/outflow
    // boundaries). That term is NOT included here, mirroring how
    // NSEBlockIntegBDF2OutletConvectiveFlux is a *separate* boundary face
    // integrator for the coupled solver. If NSSolverUncoupled is used with
    // an open/outflow boundary, add an analogous face integrator (via
    // AddBdrFaceIntegrator on the momentum predictor's ParNonlinearForm)
    // before relying on these divergence-form integrators.
    // =====================================================================

    class NSEProjVMSIntegMomentumDivForm : public mfem::NonlinearFormIntegrator {
    private:
        const InputData& idata;
        const TimeLevelFields& tlf;
        const int vdim;
        const mfem::Ordering::Type ordering;
        mfem::VectorCoefficient* f_coeff = nullptr;

    public:
        NSEProjVMSIntegMomentumDivForm(
            const InputData& idata,
            const TimeLevelFields& tlf,
            const int vdim,
            const mfem::Ordering::Type ordering,
            mfem::VectorCoefficient* f_coeff = nullptr)
            : idata(idata),
              tlf(tlf),
              vdim(vdim),
              ordering(ordering),
              f_coeff(f_coeff) {
        }

    private:
        void EvalPressureHatAndGradAtIP(mfem::ElementTransformation& T,
                                        const mfem::IntegrationPoint& ip,
                                        const ProjectionCoefficients& pc,
                                        double& p_hat,
                                        mfem::Vector& grad_p_hat) const {
            // Delegates to the shared free function so that the momentum
            // predictor, VUE RHS, and PPE RHS integrators all compute
            // p_h* identically.
            nse::EvalPressureHatAndGradAtIP(tlf, T, ip, pc, vdim, p_hat, grad_p_hat);
        }

    public:
        void AssembleElementVector(const mfem::FiniteElement& el,
                                   mfem::ElementTransformation& T,
                                   const mfem::Vector& elfun,
                                   mfem::Vector& elvec) override {
            const int ndof = el.GetDof();
            const int dim = T.GetSpaceDim();

            MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");

            elvec.SetSize(vdim * ndof);
            elvec = 0.0;

            mfem::Vector N(ndof);
            mfem::DenseMatrix dN(ndof, dim);

            mfem::Vector u(vdim);
            mfem::DenseMatrix grad_u(vdim, dim);

            mfem::Vector u_n(vdim), u_nm1(vdim), u_hist(vdim), f(vdim);

            mfem::Vector conv(vdim);
            mfem::Vector grad_p_hat(vdim);
            mfem::Vector lap_u(vdim);
            mfem::Vector lap_N(ndof);
            mfem::Vector ns_res(vdim);

            const int e = T.ElementNo;
            const mfem::IntegrationRule* ir = &tlf.femach.qspace->GetIntRule(e);

            const double dt = tlf.GetTimeStep();
            const double ctime = tlf.GetTime();
            const double nu = idata.flow_properties.nu;

            const ProjectionCoefficients pc = GetProjectionCoefficients(idata.projection_config.scheme);

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint& ip = ir->IntPoint(iq);
                T.SetIntPoint(&ip);

                el.CalcShape(ip, N);
                el.CalcPhysDShape(T, dN);
                el.CalcPhysLaplacian(T, lap_N);

                const double wdet = ip.weight * T.Weight();

                EvalVectorAtIP(elfun, ndof, vdim, ordering, N, u);
                EvalVectorGradAtIP(elfun, ndof, vdim, ordering, dN, grad_u);
                EvalVectorLaplacianAtIP(elfun, ndof, vdim, ordering, lap_N, lap_u);

                tlf.prev_1.u.GetVectorValue(T, ip, u_n);
                tlf.prev_2.u.GetVectorValue(T, ip, u_nm1);

                u_hist.SetSize(vdim);
                for (int c = 0; c < vdim; c++) {
                    u_hist(c) = pc.beta1 * u_n(c) + pc.beta2 * u_nm1(c);
                }

                double p_hat = 0.0;
                EvalPressureHatAndGradAtIP(T, ip, pc, p_hat, grad_p_hat);

                f.SetSize(vdim);
                f = 0.0;
                if (f_coeff) {
                    f_coeff->SetTime(ctime);
                    f_coeff->Eval(f, T, ip);
                }

                // -------------------------------------------------------
                // Divergence-form convection:
                //
                //   div(u tensor u)_c = (u . grad) u_c + u_c (div u)
                // -------------------------------------------------------
                const double div_u = Divergence(grad_u);

                conv.SetSize(vdim);
                conv = 0.0;

                for (int i = 0; i < vdim; i++) {
                    for (int j = 0; j < dim; j++) {
                        conv(i) += u(j) * grad_u(i, j);
                    }
                    conv(i) += u(i) * div_u;
                }

                double tauM, tauC;
                CalcTau(T, u, nu, idata.vms_config.Ci, dt, tauM, tauC);

                // -----------------------------------------------------------------
                // Coarse-scale (strong) momentum residual, divergence form:
                //
                // Res_i =
                //   sigma*u_i + (history)/dt
                //   + div(u tensor u)_i
                //   + grad(p_hat)_i
                //   - nu * laplacian(u_i)
                //   - f_i
                // -----------------------------------------------------------------
                ns_res.SetSize(vdim);
                for (int i = 0; i < vdim; i++) {
                    ns_res(i) =
                        (pc.beta0 * u(i) + u_hist(i)) / dt
                        + conv(i)
                        + grad_p_hat(i)
                        - nu * lap_u(i)
                        - f(i);
                }

                for (int a = 0; a < ndof; a++) {
                    double crossTermVelocityPart = 0.0;
                    for (int i = 0; i < dim; i++) {
                        crossTermVelocityPart += dN(a, i) * u(i);
                    }

                    double crossTermFineScalePart = 0.0;
                    for (int i = 0; i < dim; i++) {
                        crossTermFineScalePart += dN(a, i) * tauM * ns_res(i);
                    }

                    for (int c = 0; c < vdim; c++) {
                        const int ia = VDofIndex(ndof, vdim, a, c, ordering);

                        // ---------------------------------------------------------
                        // Galerkin terms:
                        //
                        //   sigma*(u,v)
                        // - (u tensor u, grad v)     [divergence/by-parts convection]
                        // + nu (grad u, grad v)
                        // + (grad p_hat, v)
                        // - (f, v)
                        // ---------------------------------------------------------

                        elvec(ia) += N(a) * ((pc.beta0 * u(c) + u_hist(c)) / dt) * wdet;

                        for (int j = 0; j < dim; j++) {
                            elvec(ia) += -u(c) * u(j) * dN(a, j) * wdet;
                        }

                        elvec(ia) += N(a) * grad_p_hat(c) * wdet;

                        for (int j = 0; j < dim; j++) {
                            elvec(ia) += nu * dN(a, j) * grad_u(c, j) * wdet;
                        }

                        elvec(ia) += -N(a) * f(c) * wdet;

                        // ---------------------------------------------------------
                        // VMS cross term 1:  (u tensor tauM*Res, grad v)
                        // ---------------------------------------------------------

                        elvec(ia) +=
                            tauM * crossTermVelocityPart * ns_res(c) * wdet;

                        // ---------------------------------------------------------
                        // VMS cross term 2:  (tauM*Res tensor u, grad v)
                        // ---------------------------------------------------------

                        elvec(ia) +=
                            crossTermFineScalePart * u(c) * wdet;

                        // ---------------------------------------------------------
                        // Reynolds-stress / fine-scale term:
                        //   -(tauM^2 Res tensor Res, grad v)
                        // ---------------------------------------------------------

                        for (int j = 0; j < dim; j++) {
                            elvec(ia) +=
                                -dN(a, j)
                                * (tauM * ns_res(c))
                                * (tauM * ns_res(j))
                                * wdet;
                        }
                    }
                }
            }
        }

        void AssembleElementGrad(const mfem::FiniteElement& el,
                                 mfem::ElementTransformation& T,
                                 const mfem::Vector& elfun,
                                 mfem::DenseMatrix& elmat) override {
            const int ndof = el.GetDof();
            const int dim = T.GetSpaceDim();

            MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");

            elmat.SetSize(vdim * ndof, vdim * ndof);
            elmat = 0.0;

            mfem::Vector N(ndof);
            mfem::DenseMatrix dN(ndof, dim);

            mfem::Vector u(vdim);
            mfem::DenseMatrix grad_u(vdim, dim);

            mfem::Vector u_n(vdim), u_nm1(vdim), u_hist(vdim), f(vdim);

            mfem::Vector conv(vdim);
            mfem::Vector grad_p_hat(vdim);
            mfem::Vector lap_u(vdim);
            mfem::Vector lap_N(ndof);
            mfem::Vector ns_res(vdim);

            mfem::DenseMatrix dResdu(vdim, vdim);
            mfem::Vector dResdu_dot_dNa(vdim);

            const int e = T.ElementNo;
            const mfem::IntegrationRule* ir = &tlf.femach.qspace->GetIntRule(e);

            const double dt = tlf.GetTimeStep();
            const double ctime = tlf.GetTime();
            const double nu = idata.flow_properties.nu;

            const ProjectionCoefficients pc = GetProjectionCoefficients(idata.projection_config.scheme);
            const double alpha = pc.beta0 / dt;

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint& ip = ir->IntPoint(iq);
                T.SetIntPoint(&ip);

                el.CalcShape(ip, N);
                el.CalcPhysDShape(T, dN);
                el.CalcPhysLaplacian(T, lap_N);

                const double wdet = ip.weight * T.Weight();

                EvalVectorAtIP(elfun, ndof, vdim, ordering, N, u);
                EvalVectorGradAtIP(elfun, ndof, vdim, ordering, dN, grad_u);
                EvalVectorLaplacianAtIP(elfun, ndof, vdim, ordering, lap_N, lap_u);

                tlf.prev_1.u.GetVectorValue(T, ip, u_n);
                tlf.prev_2.u.GetVectorValue(T, ip, u_nm1);

                u_hist.SetSize(vdim);
                for (int c = 0; c < vdim; c++) {
                    u_hist(c) = pc.beta1 * u_n(c) + pc.beta2 * u_nm1(c);
                }

                double p_hat = 0.0;
                EvalPressureHatAndGradAtIP(T, ip, pc, p_hat, grad_p_hat);

                f.SetSize(vdim);
                f = 0.0;
                if (f_coeff) {
                    f_coeff->SetTime(ctime);
                    f_coeff->Eval(f, T, ip);
                }

                const double div_u = Divergence(grad_u);

                conv.SetSize(vdim);
                conv = 0.0;

                for (int i = 0; i < vdim; i++) {
                    for (int j = 0; j < dim; j++) {
                        conv(i) += u(j) * grad_u(i, j);
                    }
                    conv(i) += u(i) * div_u;
                }

                double tauM, tauC;
                CalcTau(T, u, nu, idata.vms_config.Ci, dt, tauM, tauC);

                ns_res.SetSize(vdim);
                for (int i = 0; i < vdim; i++) {
                    ns_res(i) =
                        (pc.beta0 * u(i) + u_hist(i)) / dt
                        + conv(i)
                        + grad_p_hat(i)
                        - nu * lap_u(i)
                        - f(i);
                }

                for (int a = 0; a < ndof; a++) {
                    double crossTermVelocityPart = 0.0;
                    for (int i = 0; i < dim; i++) {
                        crossTermVelocityPart += dN(a, i) * u(i);
                    }

                    double crossTermFineScalePart = 0.0;
                    for (int i = 0; i < dim; i++) {
                        crossTermFineScalePart += dN(a, i) * tauM * ns_res(i);
                    }

                    for (int b = 0; b < ndof; b++) {
                        double conv_b = 0.0; // u . grad(N_b)
                        for (int i = 0; i < dim; i++) {
                            conv_b += dN(b, i) * u(i);
                        }

                        double gradNa_dot_gradNb = 0.0;
                        for (int j = 0; j < dim; j++) {
                            gradNa_dot_gradNb += dN(a, j) * dN(b, j);
                        }

                        // old:
                        //
                        // diff_J += nu * fe.d2N(b,j,j)
                        const double diff_J = nu * lap_N(b);

                        // -----------------------------------------------------
                        // dResdu(c,k): derivative of the divergence-form
                        // strong residual component c w.r.t. trial velocity
                        // component k at trial dof b.
                        //
                        // Res_c = sigma*u_c + div(u tensor u)_c - nu*lap(u_c) + ...
                        //
                        // d[div(u tensor u)_c]/du_k =
                        //     N_b * grad_u(c,k)
                        //   + delta_ck * (u . grad N_b)
                        //   + delta_ck * N_b * div(u)
                        //   + u_c * dN(b,k)
                        // -----------------------------------------------------
                        dResdu = 0.0;
                        for (int c = 0; c < vdim; c++) {
                            for (int k = 0; k < vdim; k++) {
                                double val = 0.0;
                                if (c == k) {
                                    val += alpha * N(b);
                                    val -= diff_J;
                                }

                                val += N(b) * grad_u(c, k);
                                if (c == k) {
                                    val += conv_b;
                                    val += N(b) * div_u;
                                }
                                val += u(c) * dN(b, k);

                                dResdu(c, k) = val;
                            }
                        }

                        // ---------------------------------------------------------
                        // 1. Direct Galerkin Jacobian: mass + diffusion + convection
                        //    (convection now in by-parts/divergence form).
                        // ---------------------------------------------------------
                        for (int c = 0; c < vdim; c++) {
                            const int ia = VDofIndex(ndof, vdim, a, c, ordering);

                            for (int k = 0; k < vdim; k++) {
                                const int ib = VDofIndex(ndof, vdim, b, k, ordering);

                                double val = 0.0;

                                if (c == k) {
                                    val += alpha * N(a) * N(b);
                                    val += nu * gradNa_dot_gradNb;
                                }

                                if (c == k) {
                                    val += -N(b) * crossTermVelocityPart;
                                }
                                val += -u(c) * N(b) * dN(a, k);

                                elmat(ia, ib) += val * wdet;
                            }
                        }

                        // dResdu_dot_dNa(k) = sum_i dN(a,i) * dResdu(i,k)
                        dResdu_dot_dNa = 0.0;
                        for (int k = 0; k < vdim; k++) {
                            double s = 0.0;
                            for (int i = 0; i < dim; i++) {
                                s += dN(a, i) * dResdu(i, k);
                            }
                            dResdu_dot_dNa(k) = s;
                        }

                        // ---------------------------------------------------------
                        // 2. Cross term 1 Jacobian:  tauM * crossTermVelocityPart * Res_c
                        // ---------------------------------------------------------
                        for (int c = 0; c < vdim; c++) {
                            const int ia = VDofIndex(ndof, vdim, a, c, ordering);

                            for (int k = 0; k < vdim; k++) {
                                const int ib = VDofIndex(ndof, vdim, b, k, ordering);

                                double val = 0.0;

                                val += tauM * dN(a, k) * N(b) * ns_res(c);
                                val += tauM * crossTermVelocityPart * dResdu(c, k);

                                elmat(ia, ib) += val * wdet;
                            }
                        }

                        // ---------------------------------------------------------
                        // 3. Cross term 2 Jacobian:  crossTermFineScalePart * u_c
                        // ---------------------------------------------------------
                        for (int c = 0; c < vdim; c++) {
                            const int ia = VDofIndex(ndof, vdim, a, c, ordering);

                            for (int k = 0; k < vdim; k++) {
                                const int ib = VDofIndex(ndof, vdim, b, k, ordering);

                                double val = tauM * dResdu_dot_dNa(k) * u(c);
                                if (c == k) {
                                    val += crossTermFineScalePart * N(b);
                                }

                                elmat(ia, ib) += val * wdet;
                            }
                        }

                        // ---------------------------------------------------------
                        // 4. Reynolds-stress Jacobian:
                        //      -tauM^2 * Res_c * Res_j * dN(a,j)   (summed over j)
                        //    = -tauM * dResdu(c,k) * crossTermFineScalePart
                        //      - tauM^2 * Res_c * dResdu_dot_dNa(k)
                        // ---------------------------------------------------------
                        for (int c = 0; c < vdim; c++) {
                            const int ia = VDofIndex(ndof, vdim, a, c, ordering);

                            for (int k = 0; k < vdim; k++) {
                                const int ib = VDofIndex(ndof, vdim, b, k, ordering);

                                double val = 0.0;
                                val += -tauM * dResdu(c, k) * crossTermFineScalePart;
                                val += -tauM * tauM * ns_res(c) * dResdu_dot_dNa(k);

                                elmat(ia, ib) += val * wdet;
                            }
                        }
                    }
                }
            }
        }
    };

    // =====================================================================
    // Outflow (Neumann) boundary flux for the divergence-form projection
    // momentum predictor. Intended outflow condition: n . grad u = 0.
    //
    // The div-form predictor assembles, as VOLUME terms only,
    //     - (u tensor u, grad v)  +  nu (grad u, grad v)  +  (grad p*, v)
    // (note grad p* is kept in volume form, NOT integrated by parts, so the
    // pressure contributes no boundary term here -- unlike the coupled solver
    // where pressure enters as -(p, div v)). Integrating the volume terms by
    // parts back to strong form exposes the natural boundary term
    //
    //     nu <grad u . n, v>  -  <(u.n) u, v>   on Gamma_N.
    //
    // With NO face integrator, the implied natural BC would therefore be the
    // spurious  nu grad u . n = (u.n) u  -- an artifact of writing convection
    // conservatively. Adding back the convective flux
    //
    //     + <(u.n) u, v>_{Gamma_N}
    //
    // cancels that artifact and leaves the intended natural condition
    //
    //     nu grad u . n = 0   <=>   n . grad u = 0.
    //
    // So for n.grad u = 0 the CONVECTIVE flux alone is what is required (this
    // is the default). Mirrors NSEBlockIntegBDF2OutletConvectiveFlux, which
    // likewise adds only the convective flux for the coupled solver.
    //
    // include_viscous_flux is provided only for completeness: setting it true
    // ALSO adds -nu <grad u . n, v>, which cancels the viscous natural term
    // entirely and yields a *free* ("pure do-nothing") outflow with NO
    // constraint on n.grad u -- this is NOT n.grad u = 0. Leave it false to
    // impose the intended condition. Added as a separate boundary-face
    // integrator restricted to the outlet marker (v = 0 on Dirichlet
    // boundaries makes it vanish there anyway).
    // =====================================================================

    class NSEProjVMSIntegMomentumOutflowFlux : public mfem::NonlinearFormIntegrator {
    private:
        const InputData& idata;
        const TimeLevelFields& tlf;
        const int vdim;
        const mfem::Ordering::Type ordering;
        const bool include_viscous_flux;

    public:
        NSEProjVMSIntegMomentumOutflowFlux(
            const InputData& idata,
            const TimeLevelFields& tlf,
            const int vdim,
            const mfem::Ordering::Type ordering,
            const bool include_viscous_flux = false)
            : idata(idata),
              tlf(tlf),
              vdim(vdim),
              ordering(ordering),
              include_viscous_flux(include_viscous_flux) {
        }

        void AssembleFaceVector(const mfem::FiniteElement& el1,
                                const mfem::FiniteElement& el2,
                                mfem::FaceElementTransformations& Tr,
                                const mfem::Vector& elfun,
                                mfem::Vector& elvec) override {
            // el2 is a dummy here: this integrator is attached via
            // AddBdrFaceIntegrator(...), so MFEM restricts it to boundary
            // faces and element 1 is the (interior) owner.
            const int ndof = el1.GetDof();
            const int dim = Tr.Elem1->GetSpaceDim();

            MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");

            elvec.SetSize(vdim * ndof);
            elvec = 0.0;

            const double nu = idata.flow_properties.nu;

            mfem::Vector N(ndof);
            mfem::DenseMatrix dN(ndof, dim);
            mfem::Vector u(vdim);
            mfem::DenseMatrix grad_u(vdim, dim);
            mfem::Vector nor(dim);

            const int order = 2 * el1.GetOrder();
            const mfem::IntegrationRule* ir =
                &mfem::IntRules.Get(Tr.GetGeometryType(), order);

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint& ip_face = ir->IntPoint(iq);

                Tr.Face->SetIntPoint(&ip_face);

                mfem::IntegrationPoint ip_el;
                Tr.Loc1.Transform(ip_face, ip_el);
                Tr.Elem1->SetIntPoint(&ip_el);

                el1.CalcShape(ip_el, N);
                EvalVectorAtIP(elfun, ndof, vdim, ordering, N, u);

                if (include_viscous_flux) {
                    el1.CalcPhysDShape(*Tr.Elem1, dN);
                    EvalVectorGradAtIP(elfun, ndof, vdim, ordering, dN, grad_u);
                }

                // nor = n |J_face|: the surface measure is already contained
                // in nor, so only multiply by the quadrature weight.
                mfem::CalcOrtho(Tr.Face->Jacobian(), nor);

                double u_dot_n = 0.0;
                for (int j = 0; j < dim; j++) {
                    u_dot_n += u(j) * nor(j);
                }

                const double w = ip_face.weight;

                for (int a = 0; a < ndof; a++) {
                    for (int c = 0; c < vdim; c++) {
                        const int ia = VDofIndex(ndof, vdim, a, c, ordering);

                        // + < (u.n) u_c , v >
                        elvec(ia) += N(a) * u(c) * u_dot_n * w;

                        // - nu < (grad u . n)_c , v >
                        if (include_viscous_flux) {
                            double gradu_dot_n = 0.0;
                            for (int j = 0; j < dim; j++) {
                                gradu_dot_n += grad_u(c, j) * nor(j);
                            }
                            elvec(ia) += -nu * N(a) * gradu_dot_n * w;
                        }
                    }
                }
            }
        }

        void AssembleFaceGrad(const mfem::FiniteElement& el1,
                              const mfem::FiniteElement& el2,
                              mfem::FaceElementTransformations& Tr,
                              const mfem::Vector& elfun,
                              mfem::DenseMatrix& elmat) override {
            const int ndof = el1.GetDof();
            const int dim = Tr.Elem1->GetSpaceDim();

            MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");

            elmat.SetSize(vdim * ndof, vdim * ndof);
            elmat = 0.0;

            const double nu = idata.flow_properties.nu;

            mfem::Vector N(ndof);
            mfem::DenseMatrix dN(ndof, dim);
            mfem::Vector u(vdim);
            mfem::Vector nor(dim);

            const int order = 2 * el1.GetOrder();
            const mfem::IntegrationRule* ir =
                &mfem::IntRules.Get(Tr.GetGeometryType(), order);

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint& ip_face = ir->IntPoint(iq);

                Tr.Face->SetIntPoint(&ip_face);

                mfem::IntegrationPoint ip_el;
                Tr.Loc1.Transform(ip_face, ip_el);
                Tr.Elem1->SetIntPoint(&ip_el);

                el1.CalcShape(ip_el, N);
                EvalVectorAtIP(elfun, ndof, vdim, ordering, N, u);

                if (include_viscous_flux) {
                    el1.CalcPhysDShape(*Tr.Elem1, dN);
                }

                mfem::CalcOrtho(Tr.Face->Jacobian(), nor);

                double u_dot_n = 0.0;
                for (int j = 0; j < dim; j++) {
                    u_dot_n += u(j) * nor(j);
                }

                const double w = ip_face.weight;

                for (int a = 0; a < ndof; a++) {
                    for (int b = 0; b < ndof; b++) {
                        // d/du_k of the viscous flux uses grad(N_b).n
                        double gradNb_dot_n = 0.0;
                        if (include_viscous_flux) {
                            for (int j = 0; j < dim; j++) {
                                gradNb_dot_n += dN(b, j) * nor(j);
                            }
                        }

                        for (int c = 0; c < vdim; c++) {
                            const int ia = VDofIndex(ndof, vdim, a, c, ordering);

                            for (int k = 0; k < vdim; k++) {
                                const int ib = VDofIndex(ndof, vdim, b, k, ordering);

                                double val = 0.0;

                                // d/du_k [ (u.n) u_c ]
                                //   = delta_ck (u.n) + u_c n_k    (times N_b)
                                double jac_ck = u(c) * nor(k);
                                if (c == k) {
                                    jac_ck += u_dot_n;
                                }
                                val += N(a) * N(b) * jac_ck;

                                // d/du_k [ -nu (grad u . n)_c ]
                                //   = -nu delta_ck (grad N_b . n)
                                if (include_viscous_flux && c == k) {
                                    val += -nu * N(a) * gradNb_dot_n;
                                }

                                elmat(ia, ib) += val * w;
                            }
                        }
                    }
                }
            }
        }
    };

    class NSEProjVMSIntegVUERHSDivForm : public mfem::LinearFormIntegrator {
    private:
        const InputData& idata;
        const TimeLevelFields& tlf;
        const int vdim;
        const mfem::Ordering::Type ordering;
        mfem::VectorCoefficient* f_coeff = nullptr;

    public:
        NSEProjVMSIntegVUERHSDivForm(
            const InputData& idata,
            const TimeLevelFields& tlf,
            const int vdim,
            const mfem::Ordering::Type ordering,
            mfem::VectorCoefficient* f_coeff = nullptr)
            : idata(idata),
              tlf(tlf),
              vdim(vdim),
              ordering(ordering),
              f_coeff(f_coeff) {
        }

        void AssembleRHSElementVect(const mfem::FiniteElement& el,
                                    mfem::ElementTransformation& T,
                                    mfem::Vector& elvec) override {
            const int ndof = el.GetDof();
            const int dim = T.GetSpaceDim();

            MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");

            elvec.SetSize(vdim * ndof);
            elvec = 0.0;

            mfem::Vector N(ndof);

            mfem::Vector u(vdim);
            mfem::Vector u_n(vdim), u_nm1(vdim), u_hist(vdim);
            mfem::Vector f(vdim);

            mfem::DenseMatrix grad_u(vdim, dim);

            mfem::Vector conv(vdim);
            mfem::Vector lap_u(vdim);
            mfem::Vector grad_p_curr(vdim);
            mfem::Vector grad_p_hat(vdim);
            mfem::Vector ns_res(vdim);

            const int e = T.ElementNo;
            const mfem::IntegrationRule* ir = &tlf.femach.qspace->GetIntRule(e);

            const double dt = tlf.GetTimeStep();
            const double ctime = tlf.GetTime();
            const double nu = idata.flow_properties.nu;

            const ProjectionCoefficients pc = GetProjectionCoefficients(idata.projection_config.scheme);

            const double sigma = pc.beta0 / dt;
            const double sigmainv = 1.0 / sigma;

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint& ip = ir->IntPoint(iq);
                T.SetIntPoint(&ip);

                el.CalcShape(ip, N);

                const double wdet = ip.weight * T.Weight();

                // Current velocity, old velocities.
                tlf.current.u.GetVectorValue(T, ip, u);
                tlf.prev_1.u.GetVectorValue(T, ip, u_n);
                tlf.prev_2.u.GetVectorValue(T, ip, u_nm1);

                // Current velocity gradient.
                tlf.current.u.GetVectorGradient(T, grad_u);

                // Current pressure gradient, and p_h* = p0*p^n + p1*p^{n-1}
                // (same helper the momentum predictor uses, so all three
                // projection steps agree on p_h*).
                tlf.current.p.GetGradient(T, grad_p_curr);

                double p_hat = 0.0;
                EvalPressureHatAndGradAtIP(tlf, T, ip, pc, vdim, p_hat, grad_p_hat);

                EvalGridFunctionVectorLaplacianAtIP(tlf.current.u, T, vdim, lap_u);

                u_hist.SetSize(vdim);
                for (int c = 0; c < vdim; c++) {
                    u_hist(c) = pc.beta1 * u_n(c) + pc.beta2 * u_nm1(c);
                }

                f.SetSize(vdim);
                f = 0.0;
                if (f_coeff) {
                    f_coeff->SetTime(ctime);
                    f_coeff->Eval(f, T, ip);
                }

                // Divergence-form convection: div(u tensor u)_i = (u.grad)u_i + u_i*div(u)
                const double div_u = Divergence(grad_u);

                conv.SetSize(vdim);
                conv = 0.0;

                for (int i = 0; i < vdim; i++) {
                    for (int j = 0; j < dim; j++) {
                        conv(i) += u(j) * grad_u(i, j);
                    }
                    conv(i) += u(i) * div_u;
                }

                double tauM, tauC;
                CalcTau(T, u, nu, idata.vms_config.Ci, dt, tauM, tauC);

                ns_res.SetSize(vdim);
                for (int i = 0; i < vdim; i++) {
                    ns_res(i) =
                        (pc.beta0 * u(i) + u_hist(i)) / dt
                        + conv(i)
                        - nu * lap_u(i)
                        + grad_p_hat(i)
                        - f(i);
                }

                for (int a = 0; a < ndof; a++) {
                    for (int c = 0; c < vdim; c++) {
                        const int ia = VDofIndex(ndof, vdim, a, c, ordering);

                        elvec(ia) +=
                            N(a) * (u(c) - tauM * ns_res(c)) * wdet;

                        elvec(ia) +=
                            -sigmainv * N(a) *
                            (grad_p_curr(c) - grad_p_hat(c)) *
                            wdet;
                    }
                }
            }
        }
    };

    class NSEProjVMSIntegPPERHSDivForm : public mfem::LinearFormIntegrator {
    private:
        const InputData& idata;
        const TimeLevelFields& tlf;
        const int vdim;
        mfem::VectorCoefficient* f_coeff = nullptr;

    public:
        NSEProjVMSIntegPPERHSDivForm(
            const InputData& idata,
            const TimeLevelFields& tlf,
            const int vdim,
            mfem::VectorCoefficient* f_coeff = nullptr)
            : idata(idata),
              tlf(tlf),
              vdim(vdim),
              f_coeff(f_coeff) {
        }

        void AssembleRHSElementVect(const mfem::FiniteElement& el,
                                    mfem::ElementTransformation& T,
                                    mfem::Vector& elvec) override {
            const int ndof = el.GetDof();
            const int dim = T.GetSpaceDim();

            MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");

            elvec.SetSize(ndof);
            elvec = 0.0;

            mfem::Vector N(ndof);
            mfem::DenseMatrix dN(ndof, dim);

            mfem::Vector u(vdim);
            mfem::Vector u_n(vdim), u_nm1(vdim), u_hist(vdim);
            mfem::Vector f(vdim);

            mfem::DenseMatrix grad_u(vdim, dim);

            mfem::Vector conv(vdim);
            mfem::Vector lap_u(vdim);
            mfem::Vector grad_p_hat(vdim);
            mfem::Vector ns_res(vdim);

            const int e = T.ElementNo;
            const mfem::IntegrationRule* ir = &tlf.femach.qspace->GetIntRule(e);

            const double dt = tlf.GetTimeStep();
            const double ctime = tlf.GetTime();
            const double nu = idata.flow_properties.nu;

            const ProjectionCoefficients pc = GetProjectionCoefficients(idata.projection_config.scheme);

            const double sigma = pc.beta0 / dt;

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint& ip = ir->IntPoint(iq);
                T.SetIntPoint(&ip);

                el.CalcShape(ip, N);
                el.CalcPhysDShape(T, dN);

                const double wdet = ip.weight * T.Weight();

                // Current velocity and old velocities.
                tlf.current.u.GetVectorValue(T, ip, u);
                tlf.prev_1.u.GetVectorValue(T, ip, u_n);
                tlf.prev_2.u.GetVectorValue(T, ip, u_nm1);

                // Current velocity gradient.
                tlf.current.u.GetVectorGradient(T, grad_u);

                // p_h* = p0*p^n + p1*p^{n-1} (same helper the momentum
                // predictor and VUE RHS use, so all three projection steps
                // agree on p_h*).
                double p_hat = 0.0;
                EvalPressureHatAndGradAtIP(tlf, T, ip, pc, vdim, p_hat, grad_p_hat);

                EvalGridFunctionVectorLaplacianAtIP(tlf.current.u, T, vdim, lap_u);

                u_hist.SetSize(vdim);
                for (int c = 0; c < vdim; c++) {
                    u_hist(c) = pc.beta1 * u_n(c) + pc.beta2 * u_nm1(c);
                }

                f.SetSize(vdim);
                f = 0.0;
                if (f_coeff) {
                    f_coeff->SetTime(ctime);
                    f_coeff->Eval(f, T, ip);
                }

                const double div_u = Divergence(grad_u);

                conv.SetSize(vdim);
                conv = 0.0;

                for (int i = 0; i < vdim; i++) {
                    for (int j = 0; j < dim; j++) {
                        conv(i) += u(j) * grad_u(i, j);
                    }
                    conv(i) += u(i) * div_u;
                }

                double tauM, tauC;
                CalcTau(T, u, nu, idata.vms_config.Ci, dt, tauM, tauC);

                ns_res.SetSize(vdim);
                for (int i = 0; i < vdim; i++) {
                    ns_res(i) =
                        (pc.beta0 * u(i) + u_hist(i)) / dt
                        + conv(i)
                        - nu * lap_u(i)
                        + grad_p_hat(i)
                        - f(i);
                }

                for (int a = 0; a < ndof; a++) {
                    // -sigma * N(a) * div(u)
                    elvec(a) += -sigma * N(a) * div_u * wdet;

                    for (int k = 0; k < dim; k++) {
                        // sigma * grad(q)_k * (-tauM * NS_residual_k)
                        elvec(a) +=
                            sigma * dN(a, k) *
                            (-tauM * ns_res(k)) *
                            wdet;

                        // grad(q)_k * grad(p_h*)_k
                        elvec(a) +=
                            dN(a, k) *
                            grad_p_hat(k) *
                            wdet;
                    }
                }
            }
        }
    };
}
