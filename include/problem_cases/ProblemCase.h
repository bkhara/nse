//
// Created by khara on 11/22/25.
//

#pragma once
#include "mfem.hpp"
#include "FEMachinery.h"

using namespace mfem;

namespace nse {
    class TimeLevelFields;
    class ProblemCase {
    protected:
        double time;
        double dt;
    public:
        int vdim;
        const InputData& idata;
        const FEMachinery& fem;
        TimeLevelFields &tlf;
        Array<int> ess_tdof_list_u;
        Array<int> ess_tdof_list_c;
        Array<int> ess_tdof_list_p;
        Array<int> ess_tdof_list_u_loaded_X;
        Array<int> ess_tdof_list_u_loaded_Y;

        bool has_outlet_bc = false;
        Array<int> outlet_marker;

        VectorCoefficient *exact_velocity = nullptr;
        Coefficient *exact_pressure = nullptr;
        VectorCoefficient *forcing_rhs = nullptr;

        ProblemCase(InputData& idata, FEMachinery& fem, TimeLevelFields &tlf)
            : idata(idata),
              fem(fem),
              tlf(tlf),
              vdim(fem.fespace_primal_u->GetVDim()) {
        }

        virtual ~ProblemCase() {
            delete forcing_rhs;
        }

        // Optional exact fields (MMS / benchmarks)
        virtual const VectorCoefficient* ExactU() const { return nullptr; }
        virtual const VectorCoefficient* ExactV() const { return nullptr; }
        virtual const VectorCoefficient* ExactA() const { return nullptr; }

        virtual const VectorCoefficient* ExactTraction() const { return nullptr; }

        /** This function fills the essential BC nodes in @p ess_tdof_list_u
         */
        virtual void ObtainBoundaryDOFs() {}

        /** Project the analytical solution onto this specified gridfield
         * The derived class needs to have an analytical solution member function
         * @param fgf
         */
        virtual void SetAnalyticalSolution(NSEGridFields& fgf) {}

        /** Project the initial conditions onto this specified gridfield
         *
         * @param fgf
         */
        virtual void SetIC(NSEGridFields& fgf) {}

        /** Project the Dirichlet BCs onto
         * the given gridfields. The particular gridfunction where the BCs
         * will be filled might differ from case to case
         * @param fgf
         */
        virtual void ApplyBC(NSEGridFields &fgf) {}

        /** This is a helper function to set the time
         * in various FunctionCoefficients owned by an instance of ProblemCase
         * @param t time
         */
        virtual void SetTime(double const t) { time = t; }

        /** This is a helper function to set the time step (=dt)
         * in various FunctionCoefficients and/or integrators owned by an instance of ProblemCase
         * @param k time
         */
        virtual void SetTimeStep(double const k) { dt = k; }

        /** Function to define various post processing tasks
         * at the end of a time step
         * @param t time
         * @param k stepsize
         */
        virtual void PostStep(const double t, const double k) {}

        virtual void RegisterParaviewFields(ParaViewDataCollection& pvdc, ParaViewDataCollection& npvdc_q) {}
    };
}