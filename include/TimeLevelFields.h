//
// Created by khara on 11/14/25.
//

#pragma once
#include "mfem.hpp"
#include "FEMachinery.h"

using namespace mfem;

namespace nse {
    struct NSEGridFields {
        ParGridFunction u;
        ParGridFunction p;
        QuadratureFunction Hq;

        explicit NSEGridFields(const FEMachinery &fem);

        // In-place addition
        NSEGridFields& operator+=(const NSEGridFields& rhs);

        // In-place subtraction
        NSEGridFields& operator-=(const NSEGridFields& rhs);

        // Out-of-place addition
        friend NSEGridFields operator+(NSEGridFields lhs,
                                            const NSEGridFields& rhs);

        // Out-of-place subtraction
        friend NSEGridFields operator-(NSEGridFields lhs,
                                            const NSEGridFields& rhs);

        void CopyFrom(const NSEGridFields &src);
    };

    class TimeLevelFields {
        double time;
        double dt;

    public:
        FEMachinery &femach;
        NSEGridFields current;
        NSEGridFields prev_1;
        NSEGridFields prev_2;
        NSEGridFields exact;
        NSEGridFields error;

        const double PSI_MIN_CLAMP = -36.0;
        const double PSI_MAX_CLAMP =  36.0;

        explicit TimeLevelFields(FEMachinery& fem);

        void Initialize();

        void UpdateTimeStepIterates();

        void RollBackFromTimeStepIterates();

        void SetTime(const double t);

        void SetTimeStep(const double k);

        double GetTime() const;

        double GetTimeStep() const;
    };
}
