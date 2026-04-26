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

    class NSSolverUncoupled : public NSSolver {
    private:
        int myrank;

        InputData& idata;
        FEMachinery& fem;
        TimeLevelFields& tlf;
        ProblemCase* pcase;

        ProjectionScheme scheme;

        NSEProjectionVelocityPredictorOperator* vel_pred_op = nullptr;
        mfem::PetscNonlinearSolver* vel_pred_solver = nullptr;

        ParBilinearForm A_ppe;
        ParBilinearForm A_vue;

        ParLinearForm B_ppe;
        ParLinearForm B_vue;

        mfem::ParGridFunction phi;

    public:
        NSSolverUncoupled(InputData& idata,
                           FEMachinery& fem,
                           TimeLevelFields& tlf,
                           ProblemCase* pcase)
            : NSSolver(idata, fem, tlf, pcase),
              myrank(mfem::Mpi::WorldRank()),
              idata(idata),
              fem(fem),
              tlf(tlf),
              pcase(pcase),
              scheme(idata.projection_config.scheme),
              A_ppe(fem.fespace_p), B_ppe(fem.fespace_p),
              A_vue(fem.fespace_primal_u), B_vue(fem.fespace_primal_u),
              phi(fem.fespace_p) {
            if (myrank == 0) {
                mfem::out << "NSSolverUncoupled constructor\n";
            }

            // velocity prediction setup
            vel_pred_op =
                new NSEProjectionVelocityPredictorOperator(
                    idata, fem, tlf, pcase, scheme);

            vel_pred_solver =
                new mfem::PetscNonlinearSolver(
                    fem.fespace_primal_u->GetComm(),
                    *vel_pred_op,
                    std::string("mo_"));
            vel_pred_solver->iterative_mode = true;

            // set up operators for PPE
            A_ppe.AddDomainIntegrator(new mfem::DiffusionIntegrator);
            A_ppe.Assemble();
            A_ppe.Finalize();

            B_ppe.AddDomainIntegrator(new NSEProjPPERHSInteg(idata, tlf, scheme, fem.el_vdim, pcase->forcing_rhs));

            // set up operators for VUE
            A_vue.AddDomainIntegrator(new mfem::VectorMassIntegrator);
            A_vue.Assemble();
            A_vue.Finalize();

            B_vue.AddDomainIntegrator(
                new NSEProjVUERHSInteg(idata, tlf, scheme, fem.el_vdim, fem.ordering, pcase->forcing_rhs)
            );
        }

        ~NSSolverUncoupled() override {
            delete vel_pred_solver;
            delete vel_pred_op;
        }

        void SolveStep(const double t, const double dt) override {
            pcase->SetTimeStep(dt);
            pcase->SetTime(t);

            const ProjectionCoefficients pc =
                GetProjectionCoefficients(scheme, dt);

            // ------------------------------------------------------------
            // Step 1: tentative velocity solve
            // ------------------------------------------------------------
            pcase->ApplyBC(tlf.current);

            {
                mfem::Vector rhs; //(fem.fespace_primal_u->GetTrueVSize());
                // rhs = 0.0;

                mfem::Vector x;
                tlf.current.u.GetTrueDofs(x);

                vel_pred_solver->Mult(rhs, x);

                tlf.current.u.SetFromTrueDofs(x);
            }

            // ------------------------------------------------------------
            // Step 2: pressure increment Poisson solve
            // ------------------------------------------------------------
            SolvePressureIncrement(pc.alpha);

            // ------------------------------------------------------------
            // Step 3: velocity correction
            // ------------------------------------------------------------
            CorrectVelocity(pc.alpha);

            // ------------------------------------------------------------
            // Step 4: pressure update
            // ------------------------------------------------------------
            // UpdatePressure();

            tlf.current.u = tlf.current.u_corrected;
        }

    private:
        void SolvePressureIncrement(const double alpha) {
            B_ppe.Assemble();

            mfem::HypreParMatrix A;
            mfem::Vector X, B;

            A_ppe.FormLinearSystem(pcase->ess_tdof_list_p, tlf.current.p, B_ppe, A, X, B);

            mfem::HypreBoomerAMG amg(A);
            mfem::CGSolver cg(A.GetComm());

            cg.SetRelTol(1.0e-12);
            cg.SetAbsTol(1.0e-14);
            cg.SetMaxIter(500);
            cg.SetPrintLevel(0);
            cg.SetPreconditioner(amg);
            cg.SetOperator(A);

            cg.Mult(B, X);

            A_ppe.RecoverFEMSolution(X, B_ppe, tlf.current.p);
        }

        void CorrectVelocity(const double alpha) {
            Array<int> empty;

            B_vue.Assemble();

            mfem::HypreParMatrix A;
            mfem::Vector X, B;

            A_vue.FormLinearSystem(
                empty,
                tlf.current.u,
                B_vue,
                A,
                X,
                B);

            mfem::CGSolver cg(A.GetComm());
            mfem::HypreSmoother jacobi(A, mfem::HypreSmoother::Jacobi);

            cg.SetRelTol(1.0e-12);
            cg.SetAbsTol(1.0e-14);
            cg.SetMaxIter(200);
            cg.SetPrintLevel(0);
            cg.SetPreconditioner(jacobi);
            cg.SetOperator(A);

            cg.Mult(B, X);

            A_vue.RecoverFEMSolution(X, B_vue, tlf.current.u_corrected);
        }
    };
}