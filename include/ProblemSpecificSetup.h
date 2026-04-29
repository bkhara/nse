//
// Created by Biswajit Khara on 11/25/25.
//

#pragma once
#include "mfem.hpp"

#include "problem_cases/PCase_MMS_2D.h"
#include "problem_cases/PCase_FPC_2D.h"
#include "problem_cases/PCase_LDC_2D.h"

namespace nse {
    inline void SelectProblemCase(InputData &idata, FEMachinery &fem, TimeLevelFields &tlf, ProblemCase *&pcase) {
        if (idata.pcase_config.pcase_id == nse::NS_MMS_2D) {
            pcase = new PCase_MMS_2D(idata, fem, tlf);
        }
        if (idata.pcase_config.pcase_id == nse::LDC_2D) {
            pcase = new PCase_LDC_2D(idata, fem, tlf);
        }
        if (idata.pcase_config.pcase_id == nse::FPC_2D) {
            pcase = new PCase_FPC_2D(idata, fem, tlf);
        }
        if (not(pcase)) {
            MFEM_ABORT("EquationCase object not initialized");
        }
    }
}
