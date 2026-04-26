//
// Created by Biswajit Khara on 4/25/26.
//

#pragma once
#include "mfem.hpp"
#include "TimeLevelFields.h"
#include "Integrators.h"
#include "problem_cases/ProblemCase.h"

namespace nse {
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
            bnlf->AddDomainIntegrator(new NSEBlockIntegBDF2(idata, tlf, 2, femach.ordering, pcase->forcing_rhs));
            if (pcase->has_outlet_bc) {
                bnlf->AddBdrFaceIntegrator(
                    new NSEBlockIntegBDF2OutletConvectiveFlux(idata, tlf, 2, femach.ordering, nullptr),
                    pcase->outlet_marker
                    );
            }
            if (idata.method_config.use_stab_vms()) {
                bnlf->AddDomainIntegrator(new NSEBlockIntegBDF2VMSConservative(idata, tlf, 2, femach.ordering, pcase->forcing_rhs));
            }
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
                    const auto *hpm = dynamic_cast<const mfem::HypreParMatrix *>(&blockJ.GetBlock(i, j));
                    MFEM_VERIFY(hpm != nullptr,
                                "Expected HypreParMatrix block in NSEBlockOperator::GetGradient.");
                    blocks_array(i, j) = hpm;
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

    class NSEProjectionVelocityPredictorOperator : public mfem::Operator {
    private:
        InputData& idata;
        FEMachinery& femach;
        TimeLevelFields& tlf;
        ProblemCase* pcase;
        ProjectionScheme scheme;

        mfem::ParNonlinearForm* nlf = nullptr;
        mutable mfem::PetscParMatrix* matJ = nullptr;

    public:
        mfem::Array<int> ess_tdof_list_u;

        NSEProjectionVelocityPredictorOperator(
            InputData& idata,
            FEMachinery& femach,
            TimeLevelFields& tlf,
            ProblemCase* pcase,
            ProjectionScheme scheme)
            : mfem::Operator(femach.fespace_primal_u->GetTrueVSize()),
              idata(idata),
              femach(femach),
              tlf(tlf),
              pcase(pcase),
              scheme(scheme) {
            ess_tdof_list_u = pcase->ess_tdof_list_u;

            nlf = new mfem::ParNonlinearForm(femach.fespace_primal_u);

            nlf->AddDomainIntegrator(
                new NSEProjMomentumVMSInteg(
                    idata,
                    tlf,
                    scheme,
                    femach.el_vdim,
                    femach.ordering,
                    pcase->forcing_rhs));
        }

        ~NSEProjectionVelocityPredictorOperator() override {
            delete nlf;
            delete matJ;
        }

        void Mult(const mfem::Vector& U, mfem::Vector& B) const override {
            nlf->Mult(U, B);

            // Essential velocity rows.
            B.SetSubVector(ess_tdof_list_u, 0.0);
        }

        mfem::Operator& GetGradient(const mfem::Vector& U) const override {
            delete matJ;
            matJ = nullptr;

            mfem::Operator& Jop = nlf->GetGradient(U);

            auto* hypreJ = dynamic_cast<mfem::HypreParMatrix*>(&Jop);
            MFEM_VERIFY(hypreJ != nullptr,
                        "Expected HypreParMatrix from velocity predictor Jacobian.");

            matJ = new mfem::PetscParMatrix(hypreJ, mfem::Operator::PETSC_MATAIJ);

            auto* Je = matJ->EliminateRowsCols(ess_tdof_list_u);
            delete Je;

            return *matJ;
        }
    };
}