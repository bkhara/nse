//
// Created by Biswajit Khara on 2/11/26.
//

#include "TimeLevelFields.h"

namespace fracture {
    FractureGridFields::FractureGridFields(const FEMachinery &fem): u(fem.fespace_primal_u),
                                                                           v(fem.fespace_primal_u),
                                                                           a(fem.fespace_primal_u),
                                                                           c(fem.fespace_primal_c),
                                                                           p(fem.fespace_p),
                                                                           psi(fem.fespace_latent),
                                                                           lambda(fem.fespace_latent),
                                                                           Hq(fem.qspace, 1),
                                                                           ePlus_q(fem.qspace, 1),
                                                                           c_q(fem.qspace, 1) {
    }

    FractureGridFields & FractureGridFields::operator+=(const FractureGridFields &rhs) {
        MFEM_VERIFY(u.Size() == rhs.u.Size(), "size mismatch in u");
        u += rhs.u;
        v += rhs.v;
        a += rhs.a;
        c += rhs.c;
        p += rhs.p;
        psi += rhs.psi;
        lambda += rhs.lambda;
        Hq += rhs.Hq;
        c_q += rhs.c_q;
        return *this;
    }

    FractureGridFields & FractureGridFields::operator-=(const FractureGridFields &rhs) {
        MFEM_VERIFY(u.Size() == rhs.u.Size(), "size mismatch in u");
        u -= rhs.u;
        v -= rhs.v;
        a -= rhs.a;
        c -= rhs.c;
        p -= rhs.p;
        psi -= rhs.psi;
        lambda -= rhs.lambda;
        Hq -= rhs.Hq;
        c_q -= rhs.c_q;
        return *this;
    }

    void FractureGridFields::CopyFrom(const FractureGridFields &src) {
        u = src.u; // copies true dof vector data
        v = src.v;
        a = src.a;
        c = src.c;
        p = src.p;
        psi = src.psi;
        lambda = src.lambda;

        Hq = src.Hq; // copies quad data
        ePlus_q = src.ePlus_q;
        c_q = src.c_q;
    }

    FractureGridFields operator+(FractureGridFields lhs, const FractureGridFields &rhs) {
        lhs += rhs;
        return lhs;
    }

    FractureGridFields operator-(FractureGridFields lhs, const FractureGridFields &rhs) {
        lhs -= rhs;
        return lhs;
    }

    TimeLevelFields::TimeLevelFields(FEMachinery &fem): femach(fem),
                                                               time(fem.idata.time_marching.t0),
                                                               dt(fem.idata.time_marching.dt),
                                                               current(fem),
                                                               prev_1(fem),
                                                               prev_2(fem),
                                                               prev_stag(fem),
                                                               prev_pg(fem),
                                                               exact(fem),
                                                               error(fem),
                                                               dc_minus(fem.fespace_primal_c),
                                                               dc_minus_cmltv(fem.fespace_primal_c),
                                                               c_violation(fem.fespace_primal_c),
                                                               residual_u(fem.fespace_primal_u),
                                                               residual_c(fem.fespace_primal_c),
                                                               energy_p(fem.fespace_cell_vdim_1),
                                                               energy_m(fem.fespace_cell_vdim_1),
                                                               stress_p(fem.fespace_cell_vdim_3),
                                                               stress_m(fem.fespace_cell_vdim_3),
                                                               Hplus(fem.fespace_cell_vdim_1) {
        Initialize();
    }

    void TimeLevelFields::Initialize() {
        {
            dc_minus = 0.;
            dc_minus_cmltv = 0.;
            c_violation = 0.;
            residual_u = 0.;
            residual_c = 0.;
            energy_p = 0.;
            energy_m = 0.;
            stress_p = 0.;
            stress_m = 0.;
            Hplus = 0.;
        }
        {
            current.u = 0;
            current.v = 0;
            current.a = 0;
            current.c = 0;
            current.p = 0;
            current.psi = PSI_MIN_CLAMP;
            current.lambda = 0.;
            current.Hq = 0;
            current.ePlus_q = 0;
            current.c_q = 0;
        }
        UpdatePGIterates();
        UpdateStaggeredIterates();
        UpdateTimeStepIterates();
    }

    void TimeLevelFields::UpdateStaggeredIterates() {
        // prev_stag = current;
        prev_stag.CopyFrom(current);
    }

    void TimeLevelFields::RollBackFromStaggeredIterates() {
        // current = prev_stag;
        current.CopyFrom(prev_stag);
    }

    void TimeLevelFields::UpdatePGIterates() {
        // prev_pg = current;
        prev_pg.CopyFrom(current);
    }

    void TimeLevelFields::RollBackFromPGIterates() {
        // current = prev_pg;
        current.CopyFrom(prev_pg);
    }

    void TimeLevelFields::UpdateTimeStepIterates() {
        // prev_2 = prev_1;
        // prev_1 = current;
        prev_2.CopyFrom(prev_1);
        prev_1.CopyFrom(current);
    }

    void TimeLevelFields::RollBackFromTimeStepIterates() {
        // current = prev_1;;
        current.CopyFrom(prev_1);
    }

    void TimeLevelFields::ReInitializePGLatentVariable() {
        current.psi = PSI_MIN_CLAMP;
        prev_pg.psi = PSI_MIN_CLAMP;
        prev_1.psi = PSI_MIN_CLAMP;
        prev_2.psi = PSI_MIN_CLAMP;
        prev_stag.psi = PSI_MIN_CLAMP;
    }

    void TimeLevelFields::ResetPGLatentVariableByClamping() {
        const int n = current.psi.Size();
        double *psi = current.psi.GetData();

        for (int i = 0; i < n; ++i) {
            psi[i] = std::clamp(psi[i], PSI_MIN_CLAMP, PSI_MAX_CLAMP);
        }
    }

    void TimeLevelFields::SetTime(const double t) {
        time = t;
    }

    void TimeLevelFields::SetTimeStep(const double k) {
        dt = k;
    }

    double TimeLevelFields::GetTime() const {
        return time;
    }

    double TimeLevelFields::GetTimeStep() const {
        return dt;
    }
}