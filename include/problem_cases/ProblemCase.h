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

        /** Project the velocity Dirichlet BCs onto the given gridfields
         * (inflow, no-slip walls, etc.)
         * @param fgf
         */
        virtual void ApplyVelocityBC(NSEGridFields &fgf) {}

        /** Project the pressure Dirichlet BCs onto the given gridfields
         * (e.g. an outlet pin for the projection method, or a single-point
         * pin for a fully enclosed domain).
         * @param fgf
         */
        virtual void ApplyPressureBC(NSEGridFields &fgf) {}

        /** Convenience helper that applies both velocity and pressure BCs
         * together. The coupled solver, which solves for (u, p)
         * simultaneously, calls this once per step. The uncoupled
         * (projection-method) solver instead calls ApplyVelocityBC() and
         * ApplyPressureBC() individually, at the point in its step sequence
         * where each is actually needed (velocity BC before the tentative
         * velocity solve, pressure BC before the pressure Poisson solve).
         * Not virtual: derived classes should override the two functions
         * above, not this one.
         * @param fgf
         */
        void ApplyBC(NSEGridFields &fgf) {
            ApplyVelocityBC(fgf);
            ApplyPressureBC(fgf);
        }

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