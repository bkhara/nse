//
// Created by Biswajit Khara on 2/11/26.
//

#include "TimeLevelFields.h"

namespace nse {
    NSEGridFields::NSEGridFields(const FEMachinery &fem): u(fem.fespace_primal_u),
                                                                           p(fem.fespace_p),
                                                                           Hq(fem.qspace, 1) {
    }

    NSEGridFields & NSEGridFields::operator+=(const NSEGridFields &rhs) {
        MFEM_VERIFY(u.Size() == rhs.u.Size(), "size mismatch in u");
        u += rhs.u;
        p += rhs.p;
        Hq += rhs.Hq;
        return *this;
    }

    NSEGridFields & NSEGridFields::operator-=(const NSEGridFields &rhs) {
        MFEM_VERIFY(u.Size() == rhs.u.Size(), "size mismatch in u");
        u -= rhs.u;
        p -= rhs.p;
        Hq -= rhs.Hq;
        return *this;
    }

    void NSEGridFields::CopyFrom(const NSEGridFields &src) {
        u = src.u; // copies true dof vector data
        p = src.p;

        Hq = src.Hq; // copies quad data
    }

    NSEGridFields operator+(NSEGridFields lhs, const NSEGridFields &rhs) {
        lhs += rhs;
        return lhs;
    }

    NSEGridFields operator-(NSEGridFields lhs, const NSEGridFields &rhs) {
        lhs -= rhs;
        return lhs;
    }

    TimeLevelFields::TimeLevelFields(FEMachinery &fem): femach(fem),
                                                               time(fem.idata.time_marching.t0),
                                                               dt(fem.idata.time_marching.dt),
                                                               current(fem),
                                                               prev_1(fem),
                                                               prev_2(fem),
                                                               exact(fem),
                                                               error(fem) {
        Initialize();
    }

    void TimeLevelFields::Initialize() {
        {
            current.u = 0;
            current.p = 0;
            current.Hq = 0;
        }
        UpdateTimeStepIterates();
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