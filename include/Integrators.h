//
// Created by khara on 11/17/25.
//

#pragma once

#include "mfem.hpp"
#include "InputData.h"
#include "TimeLevelFields.h"
#include "Utils.h"

using namespace mfem;
using namespace fracture;

namespace fracture {
    // this is a struct for holding various values (at a Gauss point)
    // that will be used for calculating the integrands
    // Furthermore, these values are at one time-level
    struct FracGPValues {
        int nsd;
        int el_vdim;

        Vector u;
        Vector v;
        Vector a;
        double c = 0.;
        double psi = 0.;
        double lambda = 0.;

        DenseMatrix du;
        DenseMatrix dv;
        DenseMatrix da;
        Vector dc;
        Vector dpsi;

        double e_density = 0.0;

        FracGPValues(int el_vdim, int nsd)
            : nsd(nsd), el_vdim(el_vdim),
              u(el_vdim), v(el_vdim), a(el_vdim),
              du(el_vdim, nsd), dv(el_vdim, nsd), da(el_vdim, nsd),
              dc(nsd), dpsi(nsd) {
        }
    };


    class FractureIntegratorBase {
    public:
        FractureIntegratorBase(const InputData& idata, const TimeLevelFields& tlf)
            : idata(idata),
              tlf(tlf),
              fdata(idata.fracture_data),
              nsd(tlf.current.u.ParFESpace()->GetMesh()->Dimension()),
              el_vdim(tlf.current.u.ParFESpace()->GetVDim()),
              vcurr(el_vdim, nsd),
              vprev_1(el_vdim, nsd),
              vprev_pg(el_vdim, nsd) {
        }

    protected:
        const InputData& idata;
        const TimeLevelFields& tlf;
        const FractureData &fdata;
        const int nsd;
        const int el_vdim;

        FracGPValues vcurr;
        FracGPValues vprev_1;
        FracGPValues vprev_pg;

        void ComputeGPValues(ElementTransformation &T, const IntegrationPoint &ip, const FractureGridFields &fgf, FracGPValues &vals) const {
            fgf.u.GetVectorValue(T, ip, vals.u);
            fgf.v.GetVectorValue(T, ip, vals.v);
            fgf.a.GetVectorValue(T, ip, vals.a);
            vals.c = fgf.c.GetValue(T, ip);
            vals.psi = fgf.psi.GetValue(T, ip);
            vals.lambda = fgf.lambda.GetValue(T, ip);

            fgf.u.GetVectorGradient(T, vals.du);
            fgf.v.GetVectorGradient(T, vals.dv);
            fgf.a.GetVectorGradient(T, vals.da);
            fgf.c.GetGradient(T, vals.dc);
            fgf.psi.GetGradient(T, vals.dpsi);

            if (vals.el_vdim == 1) {
                vals.e_density = 0;
                for (int i = 0; i < vals.nsd; i++) {
                    vals.e_density += vals.du(0,i) * vals.du(0,i);
                }
            } else if (vals.el_vdim > 1) {
                DenseMatrix epsilon(vals.el_vdim, vals.nsd);
                for (int i = 0; i < vals.el_vdim; i++) {
                    for (int j = 0; j < vals.nsd; j++) {
                        epsilon(i,j) = 0.5 * (vals.du(i, j) + vals.du(j, i));
                    }
                }

                double trace = 0.;
                for (int i = 0; i < vals.el_vdim; i++) {
                    trace += epsilon(i, i);
                }
                double eps_eps = 0.;
                for (int i = 0; i < vals.el_vdim; i++) {
                    for (int j = 0; j < vals.el_vdim; j++) {
                        eps_eps += epsilon(i, j) * epsilon(i, j);
                    }
                }
                vals.e_density = 0.5 * fdata.lambda * trace * trace + fdata.mu * eps_eps;
            }
        }

        static inline double g_AT2(double d, double kappa) {
            const double a = 1.0 - d;
            return a * a + kappa;
        }

        // Build B matrix (3 x 2*dof) for 2D engineering shear Voigt
        static inline void BuildB(const int dof, const DenseMatrix& dN, DenseMatrix& B) {
            B.SetSize(3, 2 * dof);
            B = 0.0;
            for (int a = 0; a < dof; a++) {
                const double dN_dx = dN(a, 0);
                const double dN_dy = dN(a, 1);

                const int ax = 0 * dof + a;
                const int ay = 1 * dof + a;

                B(0, ax) = dN_dx; // exx
                B(1, ay) = dN_dy; // eyy
                B(2, ax) = dN_dy; // gxy
                B(2, ay) = dN_dx; // gxy
            }
        }

        static double CalcElementSize(const IntegrationRule *ir, ElementTransformation &T, const int nsd) {
            double meas = 0.0;
            for (int i = 0; i < ir->GetNPoints(); i++) {
                const IntegrationPoint& ip = ir->IntPoint(i);
                T.SetIntPoint(&ip);
                meas += ip.weight * T.Weight();
            }
            double h = std::pow(meas, 1.0 / ((double) nsd));
            return h;
        }

        static double EvalAtQuad(const Vector& elfun, int dof, const Vector& N) {
            MFEM_VERIFY(N.Size() == dof, "EvalAtQuad: N has wrong size.");

            double result = 0.0;
            for (int a = 0; a < dof; a++) {
                result += elfun(a) * N(a);
            }
            return result;
        }

        static void EvalGradAtQuad(const Vector& elfun, int dof, const DenseMatrix& dN, Vector& grad) {
            const int nsd = dN.Width();
            MFEM_VERIFY(dN.Height() == dof, "EvalGradAtQuad: dN has wrong height.");

            grad.SetSize(nsd);
            grad = 0.;
            for (int k = 0; k < nsd; k++) {
                for (int a = 0; a < dof; a++) {
                    grad(k) += elfun(a) * dN(a, k);
                }
            }
        }

        static void EvalAtQuadMDOF(const Vector &elfun,
                               int dof,
                               int vdim,
                               Ordering::Type ordering,
                               const Vector &N, // size: dof
                               Vector &val_q) // size: vdim (output)
        {
            MFEM_VERIFY(N.Size() == dof, "EvalAtQuadMDOF: N has wrong size.");
            MFEM_VERIFY(elfun.Size() == vdim * dof,
                        "EvalAtQuadMDOF: elfun size must be vdim*dof.");

            val_q.SetSize(vdim);
            val_q = 0.0;

            auto local_idx = [=](int a, int comp) {
                if (ordering == Ordering::byNODES) {
                    // [u_x(node0), u_y(node0), u_x(node1), u_y(node1), ...]
                    return vdim * a + comp;
                } else // Ordering::byVDIM
                {
                    // [u_x(node0..dof-1), u_y(node0..dof-1), ...]
                    return dof * comp + a;
                }
            };

            for (int a = 0; a < dof; a++) {
                const double Na = N(a);
                for (int comp = 0; comp < vdim; comp++) {
                    const int idx = local_idx(a, comp);
                    val_q(comp) += elfun(idx) * Na;
                }
            }
        }

        static void EvalGradAtQuadMDOF(const Vector &elfun,
                                   int dof,
                                   int vdim,
                                   Ordering::Type ordering,
                                   const DenseMatrix &dN, // size: dof x nsd
                                   DenseMatrix &grad_q) // size: vdim x nsd (output)
        {
            const int nsd = dN.Width();
            MFEM_VERIFY(dN.Height() == dof, "EvalGradAtQuadMDOF: dN has wrong height.");
            MFEM_VERIFY(elfun.Size() == vdim * dof,
                        "EvalGradAtQuadMDOF: elfun size must be vdim*dof.");

            grad_q.SetSize(vdim, nsd);
            grad_q = 0.0;

            auto local_idx = [=](int a, int comp) {
                if (ordering == Ordering::byNODES) {
                    return vdim * a + comp;
                } else // Ordering::byVDIM
                {
                    return dof * comp + a;
                }
            };

            for (int a = 0; a < dof; a++) {
                for (int comp = 0; comp < vdim; comp++) {
                    const int idx = local_idx(a, comp);
                    const double ua = elfun(idx);
                    for (int j = 0; j < nsd; j++) {
                        // grad_q(comp, j) += u_a^comp * dN_a,j
                        grad_q(comp, j) += ua * dN(a, j);
                    }
                }
            }
        }



        void ComputeStressMatrix(const double lambda, const double mu, const DenseMatrix &strain, DenseMatrix &stress) const {
            double trace = 0.;
            for (int i = 0; i < nsd; i++) {
                trace += strain(i,i);
            }

            stress.SetSize(nsd, nsd);
            stress = 0.;
            for (int i = 0; i < stress.Height(); i++) {
                for (int j = 0; j < stress.Width(); j++) {
                    stress(i, j) = mu * (strain(i, j) + strain(j, i));
                    if (i == j) {
                        stress(i, j) += lambda * trace;
                    }
                }
            }
        }

        void ComputeStrainMatrix(const DenseMatrix &grad_u, DenseMatrix &strain) const {
            strain.SetSize(nsd, nsd);
            strain = 0.;
            for (int i = 0; i < strain.Height(); i++) {
                for (int j = 0; j < strain.Width(); j++) {
                    strain(i, j) = 0.5 * (grad_u(i, j) + grad_u(j, i));
                }
            }
        }

        void ComputeBMatrix(const FiniteElement &el, const ElementTransformation &T, const DenseMatrix &dN, /*out*/  DenseMatrix &B) const {
            int dof = el.GetDof();
            const int ord = tlf.current.u.ParFESpace()->GetOrdering();
            // Build B at this quadrature point
            // eps11 = du1/dx, eps22 = du2/dy, gamma12 = du1/dy + du2/dx
            B = 0.0;
            for (int a = 0; a < dof; a++) {
                const double dNa_dx = dN(a, 0);
                const double dNa_dy = dN(a, 1);

                int ux, uy; // column indices for u1_a and u2_a

                if (ord == Ordering::byNODES) {
                    // [u1_1, u2_1, u1_2, u2_2, ..., u1_dof, u2_dof]
                    ux = el_vdim * a; // u1 at node a
                    uy = el_vdim * a + 1; // u2 at node a
                } else // Ordering::byVDIM
                {
                    // [u1_1, ..., u1_dof, u2_1, ..., u2_dof]
                    ux = a; // u1 at node a
                    uy = dof + a; // u2 at node a
                }

                // eps11 = du1/dx
                B(0, ux) += dNa_dx;

                // eps22 = du2/dy
                B(1, uy) += dNa_dy;

                // gamma12 = du1/dy + du2/dx
                B(2, ux) += dNa_dy; // du1/dy
                B(2, uy) += dNa_dx; // du2/dx
            }
        }

    };

    class StokesIntegratorBase : public mfem::BlockNonlinearFormIntegrator {
    protected:
        const InputData &idata;
        const TimeLevelFields &tlf;

        double rho;
        double mu;

        mfem::VectorCoefficient *f_coeff;

        const int vdim;
        const mfem::Ordering::Type ordering;

        static inline int VDofIndex(const int dof,
                                    const int vdim,
                                    const int a,
                                    const int comp,
                                    const mfem::Ordering::Type ordering) {
            return (ordering == mfem::Ordering::byNODES)
                       ? a * vdim + comp
                       : comp * dof + a;
        }

        static void EvalVectorAtIP(const mfem::Vector &elfun,
                                   const int dof,
                                   const int vdim,
                                   const mfem::Ordering::Type ordering,
                                   const mfem::Vector &N,
                                   mfem::Vector &u) {
            u.SetSize(vdim);
            u = 0.0;

            for (int c = 0; c < vdim; c++) {
                for (int a = 0; a < dof; a++) {
                    u(c) += elfun(VDofIndex(dof, vdim, a, c, ordering)) * N(a);
                }
            }
        }

        static void EvalVectorGradAtIP(const mfem::Vector &elfun,
                                       const int dof,
                                       const int vdim,
                                       const mfem::Ordering::Type ordering,
                                       const mfem::DenseMatrix &dN,
                                       mfem::DenseMatrix &grad_u) {
            const int dim = dN.Width();
            grad_u.SetSize(vdim, dim);
            grad_u = 0.0;

            for (int c = 0; c < vdim; c++) {
                for (int a = 0; a < dof; a++) {
                    const double ua = elfun(VDofIndex(dof, vdim, a, c, ordering));
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
        StokesIntegratorBase(const InputData &idata,
                             const TimeLevelFields &tlf,
                             const int vdim,
                             const mfem::Ordering::Type ordering,
                             mfem::VectorCoefficient *f_coeff = nullptr)
            : idata(idata),
              tlf(tlf),
              f_coeff(f_coeff),
              vdim(vdim),
              ordering(ordering), rho(idata.stokes_mms2d_inputs.rho), mu(idata.stokes_mms2d_inputs.mu) {
        }
    };

    class StokesBlockIntegBDF2 : public StokesIntegratorBase {
    public:
        StokesBlockIntegBDF2(const InputData &idata,
                             const TimeLevelFields &tlf,
                             const int vdim,
                             const mfem::Ordering::Type ordering,
                             mfem::VectorCoefficient *f_coeff = nullptr)
            : StokesIntegratorBase(idata, tlf, vdim, ordering, f_coeff) {
        }

        void AssembleElementVector(const mfem::Array<const mfem::FiniteElement *> &el,
                                   mfem::ElementTransformation &T,
                                   const mfem::Array<const mfem::Vector *> &elfun,
                                   const mfem::Array<mfem::Vector *> &elvec) override {
            const mfem::FiniteElement &el_u = *el[0];
            const mfem::FiniteElement &el_p = *el[1];

            const int dof_u = el_u.GetDof();
            const int dof_p = el_p.GetDof();
            const int dim = T.GetSpaceDim();

            MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");

            elvec[0]->SetSize(vdim * dof_u);
            *elvec[0] = 0.0;

            elvec[1]->SetSize(dof_p);
            *elvec[1] = 0.0;

            const mfem::Vector &eu = *elfun[0];
            const mfem::Vector &ep = *elfun[1];

            mfem::Vector Nu(dof_u), Np(dof_p);
            mfem::DenseMatrix dNu(dof_u, dim);

            mfem::Vector u_np1(vdim), u_n(vdim), u_nm1(vdim), f_np1(vdim);
            mfem::DenseMatrix grad_u_np1(vdim, dim);

            const int e = T.ElementNo;
            const mfem::IntegrationRule *ir = &tlf.femach.qspace->GetIntRule(e);

            const double ctime = tlf.GetTime(); // t_{n+1}
            const double dt = tlf.GetTimeStep();

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint &ip = ir->IntPoint(iq);
                T.SetIntPoint(&ip);

                el_u.CalcShape(ip, Nu);
                el_u.CalcPhysDShape(T, dNu);
                el_p.CalcShape(ip, Np);

                const double wdet = ip.weight * T.Weight();

                // current solution values at quadrature point
                EvalVectorAtIP(eu, dof_u, vdim, ordering, Nu, u_np1);
                EvalVectorGradAtIP(eu, dof_u, vdim, ordering, dNu, grad_u_np1);
                const double p_np1 = EvalScalarAtIP(ep, Np);

                // previous time levels
                tlf.prev_1.u.GetVectorValue(T, ip, u_n);
                tlf.prev_2.u.GetVectorValue(T, ip, u_nm1);

                // forcing at current time
                f_np1.SetSize(vdim);
                f_np1 = 0.0;
                if (f_coeff) {
                    f_coeff->SetTime(ctime);
                    f_coeff->Eval(f_np1, T, ip);
                }

                // div(u^{n+1})
                double div_u_np1 = 0.0;
                for (int c = 0; c < vdim; c++) {
                    div_u_np1 += grad_u_np1(c, c);
                }

                // ------------------------------------------------------------
                // Momentum residual block
                //
                // rho*(3 u^{n+1} - 4 u^n + u^{n-1})/(2dt)
                // - mu*Delta u^{n+1} + grad p^{n+1} - f^{n+1} = 0
                //
                // Weak form:
                // (rho*3/(2dt) u^{n+1}, v)
                // - (rho*(4u^n-u^{n-1})/(2dt), v)
                // + (mu grad u^{n+1}, grad v)
                // - (p^{n+1}, div v)
                // - (f^{n+1}, v)
                // ------------------------------------------------------------
                for (int a = 0; a < dof_u; a++) {
                    for (int c = 0; c < vdim; c++) {
                        const int ia = VDofIndex(dof_u, vdim, a, c, ordering);

                        // BDF2 mass contribution from u^{n+1}
                        (*elvec[0])(ia) += rho * (3.0 / (2.0 * dt)) * u_np1(c) * Nu(a) * wdet;

                        // BDF2 history contribution
                        (*elvec[0])(ia) += -rho * (4.0 * u_n(c) - u_nm1(c)) / (2.0 * dt) * Nu(a) * wdet;

                        // diffusion contribution
                        for (int j = 0; j < dim; j++) {
                            (*elvec[0])(ia) += mu * grad_u_np1(c, j) * dNu(a, j) * wdet;
                        }

                        // pressure contribution
                        (*elvec[0])(ia) += -p_np1 * dNu(a, c) * wdet;

                        // forcing contribution
                        (*elvec[0])(ia) += -f_np1(c) * Nu(a) * wdet;
                    }
                }

                // ------------------------------------------------------------
                // Continuity residual block
                //
                // (q, div u^{n+1})
                // ------------------------------------------------------------
                for (int a = 0; a < dof_p; a++) {
                    (*elvec[1])(a) += Np(a) * div_u_np1 * wdet;
                }
            }
        }

        void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement *> &el,
                                 mfem::ElementTransformation &T,
                                 const mfem::Array<const mfem::Vector *> &elfun,
                                 const mfem::Array2D<mfem::DenseMatrix *> &elmat) override {
            const mfem::FiniteElement &el_u = *el[0];
            const mfem::FiniteElement &el_p = *el[1];

            const int dof_u = el_u.GetDof();
            const int dof_p = el_p.GetDof();
            const int dim = T.GetSpaceDim();

            MFEM_VERIFY(vdim == dim, "Assuming vdim == dim.");

            elmat(0, 0)->SetSize(vdim * dof_u, vdim * dof_u);
            *elmat(0, 0) = 0.0;

            elmat(0, 1)->SetSize(vdim * dof_u, dof_p);
            *elmat(0, 1) = 0.0;

            elmat(1, 0)->SetSize(dof_p, vdim * dof_u);
            *elmat(1, 0) = 0.0;

            elmat(1, 1)->SetSize(dof_p, dof_p);
            *elmat(1, 1) = 0.0;

            mfem::DenseMatrix &Auu = *elmat(0, 0);
            mfem::DenseMatrix &Aup = *elmat(0, 1);
            mfem::DenseMatrix &Apu = *elmat(1, 0);

            mfem::Vector Nu(dof_u), Np(dof_p);
            mfem::DenseMatrix dNu(dof_u, dim);

            const int e = T.ElementNo;
            const mfem::IntegrationRule *ir = &tlf.femach.qspace->GetIntRule(e);

            const double dt = tlf.GetTimeStep();

            for (int iq = 0; iq < ir->GetNPoints(); iq++) {
                const mfem::IntegrationPoint &ip = ir->IntPoint(iq);
                T.SetIntPoint(&ip);

                el_u.CalcShape(ip, Nu);
                el_u.CalcPhysDShape(T, dNu);
                el_p.CalcShape(ip, Np);

                const double wdet = ip.weight * T.Weight();

                // ------------------------------------------------------------
                // Auu block
                //
                // d/du [
                //   rho*3/(2dt) (u, v) + mu (grad u, grad v)
                // ]
                // ------------------------------------------------------------
                for (int a = 0; a < dof_u; a++) {
                    for (int b = 0; b < dof_u; b++) {
                        const double mass_ab =
                                rho * (3.0 / (2.0 * dt)) * Nu(a) * Nu(b) * wdet;

                        double diff_ab = 0.0;
                        for (int j = 0; j < dim; j++) {
                            diff_ab += mu * dNu(a, j) * dNu(b, j) * wdet;
                        }

                        for (int c = 0; c < vdim; c++) {
                            const int ia = VDofIndex(dof_u, vdim, a, c, ordering);
                            const int ib = VDofIndex(dof_u, vdim, b, c, ordering);
                            Auu(ia, ib) += mass_ab + diff_ab;
                        }
                    }
                }

                // ------------------------------------------------------------
                // Aup block
                //
                // d/dp [ -(p, div v) ]
                // ------------------------------------------------------------
                for (int a = 0; a < dof_u; a++) {
                    for (int c = 0; c < vdim; c++) {
                        const int ia = VDofIndex(dof_u, vdim, a, c, ordering);
                        for (int b = 0; b < dof_p; b++) {
                            Aup(ia, b) += -Np(b) * dNu(a, c) * wdet;
                        }
                    }
                }

                // ------------------------------------------------------------
                // Apu block
                //
                // d/du [ (q, div u) ]
                // ------------------------------------------------------------
                for (int a = 0; a < dof_p; a++) {
                    for (int b = 0; b < dof_u; b++) {
                        for (int c = 0; c < vdim; c++) {
                            const int ib = VDofIndex(dof_u, vdim, b, c, ordering);
                            Apu(a, ib) += Np(a) * dNu(b, c) * wdet;
                        }
                    }
                }

                // A11 block stays zero for standard unstabilized Stokes
            }
        }
    };
}
