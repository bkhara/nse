//
// Created by Biswajit Khara on 10/11/25.
//

#pragma once

#include "mfem.hpp"

using namespace mfem;

namespace fracture {
    struct CombinedGridFunctions {
        ParGridFunction &u;
        ParGridFunction &c;
        ParGridFunction &psi;

        CombinedGridFunctions(const CombinedGridFunctions &ogf) : u(ogf.u), c(ogf.c), psi(ogf.psi) {
        }

        CombinedGridFunctions(ParGridFunction &u, ParGridFunction &c, ParGridFunction &psi) : u(u), c(c),
            psi(psi) {
        }
    };

    struct StaggeredIterationError {
        real_t l2 = 0.0; // ||dc||_L2
        real_t linfty = 0.0; // ||dc||_inf
        real_t h1 = 0.0; // ||dc||_H1
        real_t h1_semi = 0.0; // ||grad dc||_L2
    };

    inline double gradR_star_single_bound(double c_prev, double psi,
                                      double max_exp_arg = 700.0) {
        // 700 ~ log(DBL_MAX) for IEEE double (approx 709), conservative to stay finite
        if (psi >  max_exp_arg) psi =  max_exp_arg;
        if (psi < -max_exp_arg) psi = -max_exp_arg; // optional; exp(-700) underflows nicely anyway
        return c_prev + std::exp(psi);
    }
    inline double gradR_star_double_bound(double c_prev, double psi) {
        if (psi >= 0.0) {
            const double z = std::exp(-psi);
            return (z * c_prev + 1.0) / (z + 1.0);
        } else {
            const double z = std::exp(psi);
            return (c_prev + z) / (1.0 + z);
        }
    }
    inline double hessR_star_single_bound(double c_prev, double psi,
                                      double max_exp_arg = 700.0) {
        // 700 ~ log(DBL_MAX) for IEEE double (approx 709), conservative to stay finite
        if (psi >  max_exp_arg) psi =  max_exp_arg;
        if (psi < -max_exp_arg) psi = -max_exp_arg; // optional; exp(-700) underflows nicely anyway
        return std::exp(psi);
    }
    inline double hessR_star_double_bound(double c_prev, double psi) {
        if (psi >= 0.0) {
            const double z = std::exp(-psi);
            return z * (1.0 - c_prev) / (z + 1.0) / (z + 1.0);
        } else {
            const double z = std::exp(psi);
            return z * (1.0 - c_prev) / (1.0 + z) / (1.0 + z);
        }
    }

    inline double ParGridFunctionMaxValue(const ParGridFunction& u) {
        // Work on true dofs so constrained/shared dofs are handled consistently
        Vector ut;
        u.GetTrueDofs(ut);

        // local max
        double lmax = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < ut.Size(); i++) {
            lmax = std::max(lmax, ut[i]);
        }

        // global max across MPI ranks
        double gmax = lmax;
        MPI_Allreduce(&lmax, &gmax, 1, MPI_DOUBLE, MPI_MAX, u.ParFESpace()->GetComm());
        return gmax;
    }

    inline double ParGridFunctionMinValue(const ParGridFunction& u) {
        Vector ut;
        u.GetTrueDofs(ut);

        double lmin = +std::numeric_limits<double>::infinity();
        for (int i = 0; i < ut.Size(); i++) {
            lmin = std::min(lmin, ut[i]);
        }

        double gmin = lmin;
        MPI_Allreduce(&lmin, &gmin, 1, MPI_DOUBLE, MPI_MIN, u.ParFESpace()->GetComm());
        return gmin;
    }



    inline void ComputePhaseIncrements(TimeLevelFields &tlf) {
        // MFEM_VERIFY(tlf.current.c.ParFESpace() == tlf.prev_1.c.ParFESpace(),
        //             "current.c and prev_1.c must live on the same ParFiniteElementSpace.");

        const double *c = tlf.current.c.GetData();
        const double *cprev = tlf.prev_1.c.GetData();
        double *dcm = tlf.dc_minus.GetData();
        double *dcm_c = tlf.dc_minus_cmltv.GetData();

        const int n = tlf.current.c.Size();
        for (int i = 0; i < n; i++) {
            const double dc = c[i] - cprev[i];
            // dcp[i] = std::max(dc, 0.0); // <dc>_+
            // dcm[i] = std::min(dc, 0.0); // <dc>_-  (signed)
            dcm[i] = std::max(-dc, 0.0); // negative part magnitude
            dcm_c[i] += std::max(-dc, 0.0);
        }
    }

    inline void PhasePullBack(ParGridFunction &cgf) {
        double *c = cgf.GetData();
        const int n = cgf.Size();
        for (int i = 0; i < n; i++) {
            const double cv = c[i];
            c[i] = std::max(cv, 0.0); // negative part magnitude
        }
    }

    inline void PhaseBPByGradRStar(const TimeLevelFields &tlf, const InputData &idata) {
        const int n = tlf.current.c.Size();

        double *c = tlf.current.c.GetData();
        const double *cprev = tlf.prev_1.c.GetData();
        const double *psicurr = tlf.current.psi.GetData();

        for (int i = 0; i < n; ++i) {
            const double cpi = cprev[i];
            const double psii = psicurr[i];

            double cval_updated;
            if (idata.method_config.use_double_bound()) {
                cval_updated = gradR_star_double_bound(cpi, psii);
            } else {
                cval_updated = gradR_star_single_bound(cpi, psii);
            }
            c[i] = cval_updated;
        }
    }

    inline void PhaseBPByNodalClamp(const TimeLevelFields &tlf, const InputData &idata) {
        const int n = tlf.current.c.Size();

        double *c = tlf.current.c.GetData();
        const double *cprev = tlf.prev_1.c.GetData();

        double bound_high = 1.0;

        for (int i = 0; i < n; ++i) {
            const double ci = c[i];
            const double cpi = cprev[i];

            double bound_low = cpi;
            c[i] = std::clamp(ci, bound_low, bound_high);
        }
    }

    inline void PGNodalBoundPreservingPhase(const TimeLevelFields &tlf, const InputData &idata) {
        if (idata.pg_config.preserve_bounds_by_grad_rstar()) {
            PhaseBPByGradRStar(tlf, idata);
        } else {
            PhaseBPByNodalClamp(tlf, idata);
        }
    }

    // Enforce irreversibility for c (0=intact, 1=broken):  c_n >= c_{n-1}
    // Projected AL update (pointwise / DoF-wise):
    //   lambda_n = max(0, lambda_{n-1} + rho * (c_{n-1} - c_n))
    void ALUpdateLambda(TimeLevelFields &tlf,
                        const InputData &idata,
                        const double rho_penalty,
                        const double tol = 0.0) // optional (e.g. 1e-12) to ignore tiny "healing"
    {
        MFEM_VERIFY(rho_penalty > 0.0, "rho_penalty must be positive.");

        mfem::ParGridFunction& c_k = tlf.current.c;
        mfem::ParGridFunction& lam_k = tlf.current.lambda;

        const mfem::ParGridFunction& c_prev = tlf.prev_1.c;
        const mfem::ParGridFunction& lam_km1 = tlf.prev_pg.lambda;

        MFEM_VERIFY(c_k.Size() == c_prev.Size(), "c size mismatch (current vs prev_1).");
        MFEM_VERIFY(lam_k.Size() == lam_km1.Size(), "lambda size mismatch (current vs prev_pg).");
        MFEM_VERIFY(c_k.Size() == lam_k.Size(), "c/lambda size mismatch.");

        const int ndof = c_k.Size();
        for (int i = 0; i < ndof; i++) {
            // violation > 0 means "healing attempted": c_{n-1} > c_n
            double bound_lo = c_prev(i);
            double viol_lo = bound_lo - c_k(i);
            if (tol > 0.0 && std::abs(viol_lo) < tol) { viol_lo = 0.0; }
            double trial_lo = std::max(0.0, lam_km1(i) + rho_penalty * viol_lo);

            double trial_hi = 0.;
            if (idata.method_config.use_double_bound()) {
                double bound_hi = 1.0;
                double viol_hi = bound_hi - c_k(i);
                if (tol > 0.0 && std::abs(viol_hi) < tol) { viol_hi = 0.0; }
                trial_hi = std::min(0.0, lam_km1(i) + rho_penalty * viol_hi);
            }

            lam_k(i) = trial_lo + trial_hi;
        }
    }

    inline void SaveSolution(int cycle, ParMesh &mesh, ParGridFunction &u_gf, ParGridFunction &c_gf, ParGridFunction &psi_gf) {
        // Give the output a name; ParaView will look for the .pvtu file
        // and automatically load the per-rank .vtu pieces.
        ParaViewDataCollection pvdc("pv", &mesh);

        // Optional: put all outputs in a subdirectory
        // pvdc.SetPrefixPath("ParaView");

        // Register the vector solution; ParaView will see it as a vector field
        pvdc.RegisterField("u_h", &u_gf);
        pvdc.RegisterField("c_h", &c_gf);
        pvdc.RegisterField("psi_h", &psi_gf);

        // Optional: annotate with a cycle/time if you like (steady-state here)
        pvdc.SetCycle(cycle);
        pvdc.SetTime(cycle);

        // Save to disk: writes e1_poisson_pv.pvtu + e1_poisson_pv_####.vtu (one per rank)
        pvdc.Save();
    }

    // Update/latch violation flags based on c_curr < c_prev
    static inline int UpdateCViolationLatch(const mfem::ParGridFunction& c_prev,
                                            const mfem::ParGridFunction& c_curr,
                                            mfem::ParGridFunction& c_violation,
                                            const double tol = 0.0) {
        MFEM_VERIFY(c_prev.ParFESpace() && c_curr.ParFESpace() && c_violation.ParFESpace(),
                    "ParGridFunctions must have valid ParFESpaces.");
        MFEM_VERIFY(c_prev.ParFESpace() == c_curr.ParFESpace() &&
                    c_prev.ParFESpace() == c_violation.ParFESpace(),
                    "c_prev, c_curr, c_violation must live on the same ParFESpace.");

        mfem::Vector cprev_td, ccurr_td, viol_td;
        c_prev.GetTrueDofs(cprev_td);
        c_curr.GetTrueDofs(ccurr_td);
        c_violation.GetTrueDofs(viol_td);

        MFEM_VERIFY(cprev_td.Size() == ccurr_td.Size() && ccurr_td.Size() == viol_td.Size(),
                    "True dof sizes do not match.");

        int newly_violated_local = 0;

        const int n = viol_td.Size();
        for (int i = 0; i < n; i++) {
            // latch: if already 1, keep it 1
            if (viol_td[i] >= 0.5) {
                viol_td[i] = 1.0;
                continue;
            }

            // set to 1 if irreversibility violated
            if (ccurr_td[i] < cprev_td[i] - tol) {
                viol_td[i] = 1.0;
                newly_violated_local++;
            }
            else {
                viol_td[i] = 0.0;
            }
        }

        c_violation.SetFromTrueDofs(viol_td);

        // return total newly violated (global)
        int newly_violated_global = newly_violated_local;
        MPI_Comm comm = c_curr.ParFESpace()->GetComm();
        MPI_Allreduce(&newly_violated_local, &newly_violated_global, 1, MPI_INT, MPI_SUM, comm);
        return newly_violated_global;
    }

    static inline int UpdateCViolation(TimeLevelFields& tlf,
                                       const double tol = 0.0) {
        return UpdateCViolationLatch(tlf.prev_1.c, tlf.current.c, tlf.c_violation, tol);
    }

    inline bool SameIntegrationRule(const IntegrationRule& a,
                                    const IntegrationRule& b,
                                    double tol = 1e-14) {
        if (a.GetNPoints() != b.GetNPoints()) { return false; }

        for (int i = 0; i < a.GetNPoints(); i++) {
            const IntegrationPoint& pa = a.IntPoint(i);
            const IntegrationPoint& pb = b.IntPoint(i);

            if (std::abs(pa.x - pb.x) > tol) return false;
            if (std::abs(pa.y - pb.y) > tol) return false;
            if (std::abs(pa.z - pb.z) > tol) return false;
            if (std::abs(pa.weight - pb.weight) > tol) return false;
        }
        return true;
    }

    static const IntegrationRule& GetTensorGLLRule(int p, int dim) {
        MFEM_VERIFY(dim == 2 || dim == 3, "Tensor-product GLL intended for QUAD/HEX only.");
        const int np = p + 1;

        static std::map<std::pair<int, int>, IntegrationRule> cache;
        auto key = std::make_pair(np, dim);
        auto it = cache.find(key);
        if (it != cache.end()) { return it->second; }

        IntegrationRule ir1d(np);
        QuadratureFunctions1D qf1d;
        qf1d.GaussLobatto(np, &ir1d); // fills points+weights in 1D

        IntegrationRule ir_td;
        if (dim == 2) { ir_td = IntegrationRule(ir1d, ir1d); }
        else if (dim == 3) { ir_td = IntegrationRule(ir1d, ir1d, ir1d); }

        auto res = cache.emplace(key, ir_td);
        return res.first->second;
    }

    static const IntegrationRule&
    GetNodalLumpedRule(const FiniteElement& el,
                       int geom, // Geometry::Type
                       int dim,
                       int mass_qorder = -1) {
        MFEM_VERIFY(geom == Geometry::SQUARE || geom == Geometry::CUBE ||
                    geom == Geometry::TRIANGLE || geom == Geometry::TETRAHEDRON,
                    "Supported geoms: QUAD/HEX/TRI/TET.");

        const int dof = el.GetDof();
        const int p = el.GetOrder();

        if (mass_qorder < 0) { mass_qorder = 2 * p; }

        using Key = std::tuple<int, int, int, int, int>;
        static std::map<Key, IntegrationRule> cache;

        Key key{geom, dim, p, mass_qorder, dof};
        auto it = cache.find(key);
        if (it != cache.end()) { return it->second; }

        // 1) Nodal points
        const IntegrationRule& nodes = el.GetNodes();
        MFEM_VERIFY(nodes.GetNPoints() == dof,
                    "Expected nodal element where #nodes == dof.");

        // 2) Choose a "mass" rule on the reference element
        const IntegrationRule* ir_mass_ptr = &IntRules.Get(geom, mass_qorder);
        const IntegrationRule& ir_mass = *ir_mass_ptr;

        // 3) Compute nodal weights w(a)
        // ----------------------------------------------------------
        Vector w(dof);
        w = 0.0;

        // Compute reference measure by integrating 1: volume_hat = Σ w_q
        double volume_hat = 0.0;
        for (int iq = 0; iq < ir_mass.GetNPoints(); ++iq) {
            volume_hat += ir_mass.IntPoint(iq).weight;
        }

        const bool is_p1_tri = (geom == Geometry::TRIANGLE && p == 1 && dof == 3);
        const bool is_p1_tet = (geom == Geometry::TETRAHEDRON && p == 1 && dof == 4);

        if (is_p1_tri || is_p1_tet) {
            // Equal weights at vertex nodes
            const double wi = volume_hat / double(dof);
            for (int a = 0; a < dof; ++a) { w(a) = wi; }
        }
        else {
            // General path: w(a) = ∫ N_a dξ
            Vector N(dof);
            for (int iq = 0; iq < ir_mass.GetNPoints(); ++iq) {
                const IntegrationPoint& ip = ir_mass.IntPoint(iq);
                el.CalcShape(ip, N);

                const double wq = ip.weight;
                for (int a = 0; a < dof; ++a) { w(a) += N(a) * wq; }
            }
        }

        // 4) Build nodal rule: points=nodes, weights=w(a)
        // ----------------------------------------------------------
        IntegrationRule ir_nodal(dof);
        for (int a = 0; a < dof; ++a) {
            IntegrationPoint ipa = nodes.IntPoint(a);
            ipa.weight = w(a);
            ir_nodal.IntPoint(a) = ipa;
        }

        auto res = cache.emplace(key, ir_nodal);
        return res.first->second;
    }

    /// Print all IntegrationPoints (x,y,z) and weights in an IntegrationRule.
    /// Optionally also print detJxW by evaluating T.Weight() at each point.
    inline void PrintIntegrationRule(const mfem::IntegrationRule& ir,
                                     const char* label = nullptr,
                                     mfem::ElementTransformation* T = nullptr,
                                     int max_points = -1,
                                     int precision = 16) {
        const int nq = ir.GetNPoints();
        const int nprint = (max_points > 0) ? std::min(nq, max_points) : nq;

        // Column formatting (keep these consistent for header + rows)
        const int WIDX = 4;
        const int WNUM = 20;
        const char* SEP1 = ": ";
        const char* SEP2 = "  ";

        if (label) { mfem::out << "\n=== " << label << " ===\n"; }
        mfem::out << "Nq = " << nq << "\n";

        mfem::out << std::setprecision(precision);

        // Header (printed using the same widths/separators as the data)
        auto PrintHeader = [&](bool with_detJxW) {
            mfem::out
                << std::setw(WIDX) << "idx" << SEP1
                << std::setw(WNUM) << "x" << SEP2
                << std::setw(WNUM) << "y" << SEP2
                << std::setw(WNUM) << "z" << SEP2
                << std::setw(WNUM) << "w_ref";

            if (with_detJxW) {
                mfem::out << SEP2 << std::setw(WNUM) << "detJxW";
            }
            mfem::out << "\n";
        };

        PrintHeader(T != nullptr);

        for (int iq = 0; iq < nprint; iq++) {
            const mfem::IntegrationPoint& ip = ir.IntPoint(iq);
            const double z = (T && T->GetDimension() == 3) ? ip.z : 0.0;

            mfem::out
                << std::setw(WIDX) << iq << SEP1
                << std::setw(WNUM) << ip.x << SEP2
                << std::setw(WNUM) << ip.y << SEP2
                << std::setw(WNUM) << z << SEP2
                << std::setw(WNUM) << ip.weight;

            if (T) {
                T->SetIntPoint(&ip);
                const double detJxW = ip.weight * T->Weight();
                mfem::out << SEP2 << std::setw(WNUM) << detJxW;
            }

            mfem::out << "\n";
        }

        if (nprint < nq) {
            mfem::out << "... (printed " << nprint << " of " << nq << " points)\n";
        }

        mfem::out << std::flush;
    }


    static void PrintIntegrationRule(const IntegrationRule &ir,
                                     std::ostream &out = std::cout,
                                     const char *name = "IR") {
        out << name << ": NPoints = " << ir.GetNPoints() << "\n";
        for (int i = 0; i < ir.GetNPoints(); i++) {
            const IntegrationPoint &ip = ir.IntPoint(i);
            // ip has x,y,z always (z may be 0 for 2D)
            out << "  i=" << i
                    << "  (x,y,z)=(" << ip.x << ", " << ip.y << ", " << ip.z << ")"
                    << "  w=" << ip.weight << "\n";
        }
    }

    static void PrintReferenceRuleAndNodes(const FiniteElement& el,
                                           int geom, int qorder,
                                           std::ostream& os = std::cout) {
        // 1) Print a built-in reference integration rule (on the reference element)
        const IntegrationRule& ir = IntRules.Get(geom, qorder);
        PrintIntegrationRule(ir, os, "Reference IntRules.Get(geom,qorder)");

        // 2) Print the nodal points for this element (also reference coordinates)
        // const IntegrationRule& nodes = el.GetNodes();
        const IntegrationRule& nodes = GetNodalLumpedRule(el, geom, el.GetDim());
        PrintIntegrationRule(nodes, os, "el.GetNodes() (nodal points)");
    }

    inline void PrintElementWiseIntegrationRule(ParMesh* pmesh,
                                                const ParFiniteElementSpace* fes,
                                                int qorder = -1,
                                                std::ostream& os = std::cout) {
        MFEM_VERIFY(pmesh, "pmesh is null.");
        MFEM_VERIFY(fes, "fes is null.");
        MFEM_VERIFY(fes->GetMesh() == pmesh,
                    "fes must be defined on the same mesh as pmesh.");

        Vector x; // physical coordinates (size will become dim)
        os << std::setprecision(16);

        for (int e = 0; e < pmesh->GetNE(); ++e) {
            ElementTransformation* T = pmesh->GetElementTransformation(e);
            MFEM_VERIFY(T, "Null ElementTransformation.");

            const FiniteElement* el = fes->GetFE(e);
            MFEM_VERIFY(el, "Null FiniteElement from fes->GetFE(e).");

            const int geom = el->GetGeomType();
            const int p = el->GetOrder();
            const int qo = (qorder >= 0) ? qorder : (2 * p);

            const IntegrationRule& ir = IntRules.Get(geom, qo);

            os << "Element e=" << e
                << "  geom=" << geom
                << "  p=" << p
                << "  qorder=" << qo
                << "  Nq=" << ir.GetNPoints()
                << "\n";

            // for (int iq = 0; iq < ir.GetNPoints(); ++iq) {
            //     const IntegrationPoint& ip = ir.IntPoint(iq);
            //     // SetIntPoint updates internal Jacobian etc. for Weight()
            //     T->SetIntPoint(&ip);
            //     // physical point
            //     T->Transform(ip, x);
            // }
            PrintReferenceRuleAndNodes(*el, geom, qorder, os);
        }
    }

    /// Check if a square DenseMatrix is (numerically) diagonal.
    /// Prints a short report to mfem::out.
    /// Returns true if the largest off-diagonal entry is below the tolerance.
    inline bool CheckDiagonal(const mfem::DenseMatrix& A,
                              double rel_tol = 1e-12,
                              double abs_tol = 1e-14,
                              const char* label = nullptr) {
        const int n = A.Height();
        MFEM_VERIFY(A.Width() == n, "CheckDiagonal expects a square matrix.");

        double max_diag = 0.0;
        double max_off = 0.0;
        int imax = -1, jmax = -1;

        for (int i = 0; i < n; i++) {
            max_diag = std::max(max_diag, std::abs(A(i, i)));
            for (int j = 0; j < n; j++) {
                if (i == j) { continue; }
                const double aij = std::abs(A(i, j));
                if (aij > max_off) {
                    max_off = aij;
                    imax = i;
                    jmax = j;
                }
            }
        }

        const double thresh = std::max(abs_tol, rel_tol * max_diag);
        const bool is_diag = (max_off <= thresh);

        if (label) { mfem::out << label << ": "; }
        mfem::out << "Diagonal check: max|diag|=" << max_diag
            << ", max|off|=" << max_off
            << ", threshold=" << thresh
            << "  => " << (is_diag ? "DIAGONAL-ish" : "NOT diagonal")
            << "\n";

        if (!is_diag && imax >= 0) {
            mfem::out << "Largest off-diagonal at (" << imax << "," << jmax
                << ") = " << A(imax, jmax) << "\n";
        }

        return is_diag;
    }

    inline void StreamReopen(ParMesh &pmesh, ParGridFunction &pgf, int step, int every_n, socketstream &sock) {
        if (step % every_n != 0) { return; }

        sock.precision(16);
        sock << "parallel " << pmesh.GetNRanks() << " " << pmesh.GetMyRank() << "\n";
        sock << "solution\n";
        sock << pmesh << pgf << std::flush;
        if (step == 0) {
            sock << "keys Rj" << std::endl;
        }
    }


}

#ifdef MFEM_USE_PETSC
#include <petscmat.h>
namespace fracture {
    static void CheckSymmetryByDotTests(const HypreParMatrix &A,
                                    MPI_Comm comm,
                                    int ntests = 3,
                                    double rtol = 1e-10) {
        const int n = A.Width(); // true-dof size on this rank

        Vector x(n), y(n), Ax(n), Ay(n);

        double worst_rel = 0.0;

        for (int t = 0; t < ntests; t++) {
            // simple deterministic-ish fill (avoid RNG hassles)
            for (int i = 0; i < n; i++) {
                x[i] = std::sin(0.1 * (i + 1) + 1.37 * (t + 1));
                y[i] = std::cos(0.07 * (i + 1) + 2.11 * (t + 1));
            }

            A.Mult(x, Ax);
            A.Mult(y, Ay);

            // global dot products
            const double xTAy = InnerProduct(comm, x, Ay);
            const double yTAx = InnerProduct(comm, y, Ax);

            const double num = std::fabs(xTAy - yTAx);
            const double den = std::fabs(xTAy) + std::fabs(yTAx) + 1e-30;
            const double rel = num / den;

            worst_rel = std::max(worst_rel, rel);
        }

        int rank;
        MPI_Comm_rank(comm, &rank);
        if (rank == 0) {
            mfem::out << "Jacobian symmetry dot-test worst_rel = " << worst_rel
                    << " (rtol=" << rtol << ")\n";
            if (worst_rel > rtol) {
                mfem::out << "WARNING: Jacobian appears non-symmetric (dot-test).\n";
            }
        }
    }

    static void CheckSymmetryViaPetscAij(const mfem::HypreParMatrix &H,
                                     double rtol = 1e-10)
    {
        // Build a *temporary* PETSc AIJ matrix from the Hypre matrix
        mfem::PetscParMatrix A_petsc(&H, mfem::Operator::PETSC_MATAIJ); // default is MATAIJ :contentReference[oaicite:2]{index=2}
        Mat A = (Mat)A_petsc; // PetscParMatrix -> Mat :contentReference[oaicite:3]{index=3}

        Mat AT = nullptr, D = nullptr;
        PetscReal nA = 0.0, nD = 0.0;

        MatTranspose(A, MAT_INITIAL_MATRIX, &AT);
        MatDuplicate(A, MAT_COPY_VALUES, &D);
        MatAXPY(D, -1.0, AT, DIFFERENT_NONZERO_PATTERN); // safe

        MatNorm(A, NORM_FROBENIUS, &nA);
        MatNorm(D, NORM_FROBENIUS, &nD);

        const PetscReal rel = nD / (nA + PETSC_SMALL);
        PetscPrintf(PetscObjectComm((PetscObject)A),
                    "Symmetry check: ||A-A^T||_F = %.3e, ||A||_F = %.3e, rel = %.3e\n",
                    (double)nD, (double)nA, (double)rel);

        if (rel > rtol)
        {
            PetscPrintf(PetscObjectComm((PetscObject)A),
                        "WARNING: Jacobian appears NON-symmetric (rtol=%.1e)\n", rtol);
        }

        MatDestroy(&D);
        MatDestroy(&AT);
    }
}
#endif