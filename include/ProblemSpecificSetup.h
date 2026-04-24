//
// Created by Biswajit Khara on 11/25/25.
//

#pragma once
#include "mfem.hpp"

#include "problem_cases/PCase_Stokes_MMS_2D.h"
#include "problem_cases/PCase_FPC_2D.h"

namespace fracture {
    inline void SelectProblemCase(InputData &idata, FEMachinery &fem, TimeLevelFields &tlf, ProblemCase *&pcase) {
        if (idata.pcase_config.pcase_id == fracture::STOKES_MMS_2D) {
            pcase = new PCase_Stokes_MMS_2D(idata, fem, tlf);
        }
        if (idata.pcase_config.pcase_id == fracture::FPC_2D) {
            pcase = new PCase_FPC_2D(idata, fem, tlf);
        }
        if (not(pcase)) {
            MFEM_ABORT("EquationCase object not initialized");
        }
    }
}
