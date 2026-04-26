//
// Created by Biswajit Khara on 4/25/26.
//

#pragma once
#include <petscsnes.h>
#include "mfem.hpp"
#include "NSOperators.h"

namespace nse {
    class NSSolver {
    protected:
        int myrank;

        InputData &idata;
        const FEMachinery &fem;
        TimeLevelFields &tlf;

        ProblemCase *pcase;

    public:
        NSSolver(InputData &idata, FEMachinery &fem, TimeLevelFields &tlf, ProblemCase *pcase)
            : myrank(Mpi::WorldRank()),
              idata(idata), fem(fem), tlf(tlf), pcase(pcase) {
            if (myrank == 0) {
                mfem::out << "NSSolver constructor\n";
            }

            pcase->ObtainBoundaryDOFs();
        }

        virtual ~NSSolver() = default;

        virtual void SolveStep(const double t, const double dt) {
            MFEM_ABORT("SolverNSBase::SolveStep: base class is being instantiated. Use a derived class.");
        }
    };

    class NSSolverCoupled : public NSSolver {

        NSEBlockOperator *nse_block_op;
        PetscNonlinearSolver *petsc_nonlinear_solver = nullptr;
        SNES snes;

    public:
        NSSolverCoupled(InputData &idata, FEMachinery &fem, TimeLevelFields &tlf, ProblemCase *pcase)
            : NSSolver(idata, fem, tlf, pcase) {
            if (myrank == 0) {
                mfem::out << "SolverNSCoupled constructor\n";
            }

            nse_block_op = new NSEBlockOperator(fem.fespace_primal_u->GetTrueVSize() + fem.fespace_p->GetTrueVSize(),idata, fem, tlf, pcase);
            petsc_nonlinear_solver = new PetscNonlinearSolver(fem.fespace_primal_u->GetComm(), *nse_block_op,
                                                                  std::string(nse::PetscSolverPrefix::ELASTICITY));
            petsc_nonlinear_solver->iterative_mode = true;
            snes = static_cast<petsc::SNES>(*petsc_nonlinear_solver);
        }

        ~NSSolverCoupled() override {
            delete petsc_nonlinear_solver;
            delete nse_block_op;
        }

        void SolveStep(const double t, const double dt) override {
            pcase->SetTimeStep(dt);
            pcase->SetTime(t);

            pcase->ApplyBC(tlf.current);
            SolveStepNonlinear();

            pcase->SetTime(t); // set t again because it might have been reset by the stage calls
        }

        void SolveStepNonlinear() const {
            PetscInt snes_its;
            PetscReal snes_residual;
            SNESConvergedReason snes_converged_reason = SNES_DIVERGED_FUNCTION_DOMAIN;
            // just using a negative value to initialize
            {
                Vector zerovec;
                BlockVector xvec(nse_block_op->t_offsets);
                tlf.current.u.GetTrueDofs(xvec.GetBlock(0));
                tlf.current.p.GetTrueDofs(xvec.GetBlock(1));

                petsc_nonlinear_solver->Mult(zerovec, xvec);

                tlf.current.u.SetFromTrueDofs(xvec.GetBlock(0));
                tlf.current.p.SetFromTrueDofs(xvec.GetBlock(1));
                //PhaseField_ProjectBoundaryValues(tlf.current.c);
            }
            SNESGetIterationNumber(snes, &snes_its);
            SNESGetConvergedReason(snes, &snes_converged_reason);
            SNESGetFunctionNorm(snes, &snes_residual);
        }
    };
}