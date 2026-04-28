//
// Created by khara on 11/14/25.
//

#pragma once
#include "mfem.hpp"
#include "InputData.h"

using namespace mfem;

namespace nse {
    class FEMachinery {
    public:
        InputData &idata;
        ParMesh* mesh = nullptr;
        FiniteElementCollection *fec_primal = nullptr;
        FiniteElementCollection *fec_latent = nullptr;
        FiniteElementCollection *fec_cell_data = nullptr;
        ParFiniteElementSpace *fespace_primal_u = nullptr;
        ParFiniteElementSpace *fespace_primal_c = nullptr;
        ParFiniteElementSpace *fespace_p = nullptr;
        ParFiniteElementSpace *fespace_latent = nullptr;
        ParFiniteElementSpace *fespace_cell_vdim_1 = nullptr;
        ParFiniteElementSpace *fespace_cell_vdim_2 = nullptr;
        ParFiniteElementSpace *fespace_cell_vdim_3 = nullptr;
        QuadratureSpace *qspace = nullptr;

        Array<ParFiniteElementSpace *> fespace_block_cpsi;
        Array<ParFiniteElementSpace *> fespace_block_up;

        Ordering::Type ordering = Ordering::byVDIM;

        int vel_vdim = 2;

        explicit FEMachinery(InputData& idata);

        ~FEMachinery();

        void SetUpMesh();

        void SetUpFESpaces();

        void SaveMesh() const;
    };
}
