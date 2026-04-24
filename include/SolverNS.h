//
// Created by Biswajit Khara on 4/18/26.
//

#pragma once
#include <petscsnes.h>

#include "mfem.hpp"
#include "TimeLevelFields.h"
#include "Integrators.h"
#include "problem_cases/ProblemCase.h"
#include "problem_cases/PCase_Stokes_MMS_2D.h"

using namespace mfem;

namespace fracture {
    class NSEBlockOperator : public mfem::Operator {
    public:
        InputData &idata;
        FEMachinery &femach;
        TimeLevelFields &tlf;
        ProblemCase *pcase;

        mfem::Array<int> t_offsets;
        mfem::Array<int> ess_tdof_list_u, ess_tdof_list_p;
        Array<int> ess_tdof_list_all;

        mfem::ParBlockNonlinearForm *bnlf = nullptr;
        mutable mfem::PetscParMatrix *matJ = nullptr;

        NSEBlockOperator(int operator_size,
                            InputData &idata,
                            FEMachinery &femach,
                            TimeLevelFields &tlf,
                            ProblemCase *pcase)
            : mfem::Operator(operator_size),
              idata(idata), femach(femach), tlf(tlf), pcase(pcase), t_offsets(3) {
            t_offsets[0] = 0;
            t_offsets[1] = femach.fespace_primal_u->GetTrueVSize();
            t_offsets[2] = femach.fespace_p->GetTrueVSize();
            t_offsets.PartialSum();

            // set up the tdof lists
            {
                ess_tdof_list_u = pcase->ess_tdof_list_u;
                ess_tdof_list_p = pcase->ess_tdof_list_p;

                // Append velocity DOFs (indices remain [0, t_offsets[1]))
                ess_tdof_list_all.Append(ess_tdof_list_u);

                // Append pressure DOFs (indices shifted by t_offsets[1])
                int u_offset = t_offsets[1];
                for (int i = 0; i < ess_tdof_list_p.Size(); i++) {
                    ess_tdof_list_all.Append(ess_tdof_list_p[i] + u_offset);
                }
            }

            bnlf = new mfem::ParBlockNonlinearForm(femach.fespace_block_up);
            bnlf->AddDomainIntegrator(new StokesBlockIntegBDF2(idata, tlf, 2, femach.ordering, pcase->forcing_rhs));
        }

        ~NSEBlockOperator() override {
            delete bnlf;
            delete matJ;
        }

        void Mult(const mfem::Vector &U, mfem::Vector &B) const override {
            // mfem::BlockVector bU(const_cast<double *>(U.GetData()), t_offsets);

            mfem::BlockVector bB(B.GetData(), t_offsets);
            bB = 0.0;

            // tlf.current.u.SetFromTrueDofs(bU.GetBlock(0));
            // tlf.current.p.SetFromTrueDofs(bU.GetBlock(1));

            bnlf->Mult(U, B);

            bB.GetBlock(0).SetSubVector(ess_tdof_list_u, 0.0);
            bB.GetBlock(1).SetSubVector(ess_tdof_list_p, 0.0);
        }

        // mfem::Operator &GetGradient(const mfem::Vector &U) const override {
        //     delete matJ;
        //     matJ = nullptr;
        //
        //     // mfem::BlockVector bU(const_cast<double *>(U.GetData()), t_offsets);
        //     // tlf.current.u.SetFromTrueDofs(bU.GetBlock(0));
        //     // tlf.current.p.SetFromTrueDofs(bU.GetBlock(1));
        //
        //     mfem::Operator &J = bnlf->GetGradient(U);
        //     const Operator::Type OPTYPE = Operator::PETSC_MATAIJ;
        //     matJ = new PetscParMatrix(femach.fespace_primal_u->GetComm(), &J, OPTYPE);
        //     MFEM_VERIFY(matJ, "Expected HypreParMatrix from bnlf->GetGradient(U)");
        //     {
        //         auto *Je = matJ->EliminateRowsCols(pcase->ess_tdof_list_u);
        //         delete Je;
        //     }
        //     {
        //         auto *Je = matJ->EliminateRowsCols(pcase->ess_tdof_list_p);
        //         delete Je;
        //     }
        //     return *matJ;
        // }

        mfem::Operator &GetGradient(const mfem::Vector &U) const override {
            // 1. Clean up the previous Jacobian to prevent memory leaks
            delete matJ;
            matJ = nullptr;

            // bnlf->GetGradient(U) returns a BlockOperator containing HypreParMatrix blocks
            mfem::BlockOperator &blockJ = static_cast<mfem::BlockOperator &>(bnlf->GetGradient(U));

            // 3. Prepare blocks for merging using the non-deprecated (const) API
            mfem::Array2D<const mfem::HypreParMatrix *> blocks_array(2, 2);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    // Retrieve the block and cast to const HypreParMatrix*
                    mfem::Operator *op = &blockJ.GetBlock(i, j);
                    blocks_array(i, j) = dynamic_cast<const mfem::HypreParMatrix *>(op);
                }
            }

            // 4. Flatten the 2x2 block structure into a single monolithic HypreParMatrix
            mfem::HypreParMatrix *hypreJ = mfem::HypreParMatrixFromBlocks(blocks_array);
            matJ = new mfem::PetscParMatrix(hypreJ, mfem::Operator::PETSC_MATAIJ);
            auto *Je = matJ->EliminateRowsCols(ess_tdof_list_all);

            delete hypreJ;
            delete Je; // Je contains the eliminated entries; usually not needed here

            return *matJ;
        }
    };
    class SolverNS {
        int myrank;

        InputData &idata;
        const FEMachinery &fem;
        TimeLevelFields &tlf;

        ProblemCase *pcase;

        NSEBlockOperator *nse_block_op;

        Vector B, X;

        PetscNonlinearSolver *petsc_nonlinear_solver;
        SNES snes;

    public:
        SolverNS(InputData &idata, FEMachinery &fem, TimeLevelFields &tlf, ProblemCase *pcase)
            : myrank(Mpi::WorldRank()),
              idata(idata), fem(fem), tlf(tlf), pcase(pcase) {
            if (myrank == 0) {
                mfem::out << "NSSolver constructor\n";
            }

            pcase->ObtainElasticityBoundaryDOFs();

            nse_block_op = new NSEBlockOperator(fem.fespace_primal_u->GetTrueVSize() + fem.fespace_p->GetTrueVSize(),idata, fem, tlf, pcase);
            petsc_nonlinear_solver = new PetscNonlinearSolver(fem.fespace_primal_u->GetComm(), *nse_block_op,
                                                                  std::string(fracture::PetscSolverPrefix::ELASTICITY));
            petsc_nonlinear_solver->iterative_mode = true;
            snes = static_cast<petsc::SNES>(*petsc_nonlinear_solver);
        }

        ~SolverNS() {
            delete petsc_nonlinear_solver;
            delete nse_block_op;
        }

        void SolveStep(const double t, const double dt) {
            pcase->SetTimeStep(dt);
            pcase->SetTime(t);

            pcase->ApplyElasticityBC(tlf.current);
            SolveStepNonlinear();

            pcase->SetTime(t); // set t again because it might have been reset by the stage calls
            // pcase->ApplyElasticityBC(tlf.current); // just to make sure
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