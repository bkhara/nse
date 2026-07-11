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
                            FEMachinery &fem,
                            TimeLevelFields &tlf,
                            ProblemCase *pcase)
            : mfem::Operator(operator_size),
              idata(idata), femach(fem), tlf(tlf), pcase(pcase), t_offsets(3) {
            t_offsets[0] = 0;
            t_offsets[1] = fem.fespace_primal_u->GetTrueVSize();
            t_offsets[2] = fem.fespace_p->GetTrueVSize();
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

            bnlf = new mfem::ParBlockNonlinearForm(fem.fespace_block_up);
            if (idata.time_marching.marching_scheme == TimeMarchingScheme::BDF2) {
                // NSEBlockIntegBDF2 (and its outlet-flux/VMS/SUPG/PSPG companions)
                // are all written in the divergence (conservative) convection
                // form. No convective- or skew-symmetric-form BDF2 integrator
                // exists yet, so anything else must abort rather than silently
                // running with the wrong discrete form.
                if (!idata.method_config.is_divergence_form()) {
                    MFEM_ABORT("NSEBlockOperator: BDF2 time marching is only "
                        "implemented with the divergence (conservative) "
                        "convection form (ConvectionForms::DIV_FORM). The "
                        "convective and skew-symmetric forms are not yet "
                        "implemented for BDF2.");
                }

                bnlf->AddDomainIntegrator(new NSEBlockIntegBDF2(idata, tlf, fem.vel_vdim, fem.ordering, pcase->forcing_rhs));
                if (pcase->has_outlet_bc) {
                    bnlf->AddBdrFaceIntegrator(
                        new NSEBlockIntegBDF2OutletConvectiveFlux(idata, tlf, fem.vel_vdim, fem.ordering, nullptr),
                        pcase->outlet_marker
                        );
                }
                if (idata.method_config.use_stab_vms()) {
                    bnlf->AddDomainIntegrator(new NSEBlockIntegBDF2VMSConservative(idata, tlf, fem.vel_vdim, fem.ordering, pcase->forcing_rhs));
                } else if (idata.method_config.use_stab_sups()) {
                    if (idata.method_config.sups.use_supg) {
                        bnlf->AddDomainIntegrator(new NSEBlockIntegBDF2SUPGConservativeAddOn(idata, tlf, fem.vel_vdim, fem.ordering, pcase->forcing_rhs));
                    }
                    if (idata.method_config.sups.use_pspg) {
                        bnlf->AddDomainIntegrator(new NSEBlockIntegBDF2PSPGConservativeAddOn(idata, tlf, fem.vel_vdim, fem.ordering, pcase->forcing_rhs));
                    }
                }
            } else if (idata.time_marching.marching_scheme == TimeMarchingScheme::CN) {
                // NSEBlockIntegCrankNicolson is written in the direct
                // convective form. No divergence- or skew-symmetric-form CN
                // integrator exists yet.
                if (!idata.method_config.is_convective_form()) {
                    MFEM_ABORT("NSEBlockOperator: Crank-Nicolson time marching "
                        "is only implemented with the convective convection "
                        "form (ConvectionForms::CONV_FORM). The divergence and "
                        "skew-symmetric forms are not yet implemented for CN.");
                }

                bnlf->AddDomainIntegrator(new NSEBlockIntegCrankNicolson(idata, tlf, fem.vel_vdim,fem.ordering, pcase->forcing_rhs));
            } else {
                MFEM_ABORT("Unknown TimeMarchingScheme. Must be either bdf2 or cn");
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

        mfem::ParNonlinearForm* nlf = nullptr;
        mutable mfem::PetscParMatrix* matJ = nullptr;

    public:
        mfem::Array<int> ess_tdof_list_u;

        NSEProjectionVelocityPredictorOperator(
            InputData& idata,
            FEMachinery& femach,
            TimeLevelFields& tlf,
            ProblemCase* pcase)
            : mfem::Operator(femach.fespace_primal_u->GetTrueVSize()),
              idata(idata),
              femach(femach),
              tlf(tlf),
              pcase(pcase) {
            ess_tdof_list_u = pcase->ess_tdof_list_u;

            nlf = new mfem::ParNonlinearForm(femach.fespace_primal_u);

            if (idata.method_config.is_convective_form()) {
                nlf->AddDomainIntegrator(
                    new NSEProjVMSIntegMomentumConvForm(idata, tlf, femach.vel_vdim, femach.ordering, pcase->forcing_rhs));
            } else if (idata.method_config.is_divergence_form()) {
                nlf->AddDomainIntegrator(
                    new NSEProjVMSIntegMomentumDivForm(idata, tlf, femach.vel_vdim, femach.ordering, pcase->forcing_rhs));

                // Intended projection outflow BC: n . grad u = 0. Writing
                // convection conservatively, -(u tensor u, grad v), otherwise
                // imposes the spurious natural condition nu grad u.n = (u.n)u.
                // Adding the convective flux <(u.n)u, v> on the outlet cancels
                // that artifact and restores n . grad u = 0. (Only convective
                // flux; the viscous flux is left as the natural h = 0 term.)
                if (pcase->has_outlet_bc) {
                    nlf->AddBdrFaceIntegrator(
                        new NSEProjVMSIntegMomentumOutflowFlux(
                            idata, tlf, femach.vel_vdim, femach.ordering,
                            /*include_viscous_flux=*/false),
                        pcase->outlet_marker);
                }
            } else {
                MFEM_ABORT("NSEProjectionVelocityPredictorOperator: "
                    "skew-symmetric convection form is not yet implemented "
                    "for the projection momentum predictor.");
            }
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