//
// Created by khara on 11/22/25.
//

#pragma once
#include "mfem.hpp"
#include "FEMachinery.h"

using namespace mfem;

namespace fracture {
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
        Array<int> outlet_marker;

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
        virtual void ObtainElasticityBoundaryDOFs() {}

        /** This function fills the essential BC nodes in @p ess_tdof_list_c
         */
        virtual void ObtainPhaseFieldBoundaryDOFs() {}

        /** Project the analytical solution onto this specified gridfield
         * The derived class needs to have an analytical solution member function
         * @param fgf
         */
        virtual void SetAnalyticalSolution(FractureGridFields& fgf) {}

        /** Project the initial conditions onto this specified gridfield
         *
         * @param fgf
         */
        virtual void SetElasticityIC(FractureGridFields& fgf) {}

        /** Project the initial conditions onto this specified gridfield (for c and psi)
         *
         * @param fgf
         */
        virtual void SetPhaseFieldIC(FractureGridFields& fgf) {}

        /** Project the Dirichlet BCs onto
         * the given gridfields. The particular gridfunction where the BCs
         * will be filled might differ from case to case
         * @param fgf
         */
        virtual void ApplyElasticityBC(FractureGridFields &fgf) {}

        /** Project the Dirichlet BCs onto
         * the given gridfields.
         * @param fgf
         */
        virtual void ApplyPhaseFieldBC(FractureGridFields &fgf) {}

        /** Add the integrators for solving for the first step acceleration.
         * @param a bilinear form (lhs)
         * @param b linear form (rhs)
         */
        virtual void AddElasticityIntegratorsZeroStep(ParBilinearForm &a, ParLinearForm &b) {}

        /** Add the integrators for this case.
         * The integrators can change based on the case, the time marching algorith
         * or other method related considerations
         * @param a bilinear form (lhs)
         * @param b linear form (rhs)
         */
        virtual void AddElasticityIntegrators(ParBilinearForm &a, ParLinearForm &b) {}
        virtual void AddElasticityNonlinearIntegrators(ParNonlinearForm &nlf) {}

        /** Add the RHS integrators for this case.
         * The integrators can change based on the case, the time marching algorith
         * or other method related considerations
         * @param a bilinear form (lhs)
         * @param b linear form (rhs)
         */
        virtual void AddElasticityRHSIntegrators(ParLinearForm &b) {}

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

        virtual double CalculateStaggeredIterationError() {
            return 0.;
        }

        virtual void RegisterParaviewFields(ParaViewDataCollection& pvdc, ParaViewDataCollection& npvdc_q) {}

        virtual void CalcElasticLoading(Operator &Mop, ParNonlinearForm & nlfK, ParLinearForm &b) {}
        virtual void UpdateQuadratureFunctions() {}
    };
}