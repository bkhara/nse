//
// Created by khara on 11/14/25.
//

#pragma once
#include "mfem.hpp"
#include "FEMachinery.h"

using namespace mfem;

namespace fracture {
    struct FractureGridFields {
        ParGridFunction u;
        ParGridFunction v;
        ParGridFunction a;
        ParGridFunction c;
        ParGridFunction p;
        ParGridFunction psi;
        ParGridFunction lambda;
        QuadratureFunction Hq;
        QuadratureFunction ePlus_q;
        QuadratureFunction c_q;

        explicit FractureGridFields(const FEMachinery &fem);

        // In-place addition
        FractureGridFields& operator+=(const FractureGridFields& rhs);

        // In-place subtraction
        FractureGridFields& operator-=(const FractureGridFields& rhs);

        // Out-of-place addition
        friend FractureGridFields operator+(FractureGridFields lhs,
                                            const FractureGridFields& rhs);

        // Out-of-place subtraction
        friend FractureGridFields operator-(FractureGridFields lhs,
                                            const FractureGridFields& rhs);

        void CopyFrom(const FractureGridFields &src);
    };

    class TimeLevelFields {
        double time;
        double dt;

    public:
        FEMachinery &femach;
        FractureGridFields current;
        FractureGridFields prev_1;
        FractureGridFields prev_2;
        FractureGridFields prev_stag;
        FractureGridFields prev_pg;
        FractureGridFields exact;
        FractureGridFields error;
        ParGridFunction dc_minus; // the magnitude
        ParGridFunction dc_minus_cmltv;
        ParGridFunction c_violation;
        ParGridFunction residual_u;
        ParGridFunction residual_c;
        ParGridFunction energy_p;
        ParGridFunction energy_m;
        ParGridFunction stress_p;
        ParGridFunction stress_m;
        ParGridFunction Hplus;

        const double PSI_MIN_CLAMP = -36.0;
        const double PSI_MAX_CLAMP =  36.0;

        explicit TimeLevelFields(FEMachinery& fem);

        void Initialize();

        void UpdateStaggeredIterates();

        void RollBackFromStaggeredIterates();

        void UpdatePGIterates();

        void RollBackFromPGIterates();

        void UpdateTimeStepIterates();

        void RollBackFromTimeStepIterates();

        void ReInitializePGLatentVariable();

        void ResetPGLatentVariableByClamping();

        void SetTime(const double t);

        void SetTimeStep(const double k);

        double GetTime() const;

        double GetTimeStep() const;
    };
}
