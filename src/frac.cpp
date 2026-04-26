//
// Created by Biswajit Khara on 10/9/25.
//

#include "mfem.hpp"
#include <fstream>
#include <iostream>

#ifndef MFEM_USE_PETSC
#error This example requires that MFEM is built with MFEM_USE_PETSC=YES
#endif

#include "Checkpointing.h"
#include "InputData.h"
#include "FEMachinery.h"
#include "TimeLevelFields.h"
#include "ParallelTimer.h"
#include "Utils.h"
#include "problem_cases/ProblemCase.h"
#include "ProblemSpecificSetup.h"
#include "NSSolver.h"
#include "TimeStepper.h"

using namespace mfem;
using std::endl;
using std::flush;
using std::max;
using std::min;

int main(int argc, char *argv[]) {
    // MPI and Hypre initializations
    Mpi::Init(&argc, &argv);
    Hypre::Init();

    // parse arguments
    const char *petscrc_file = "ns.petsc";
    const char *config_file = "config.txt";

    OptionsParser args(argc, argv);
    args.AddOption(&petscrc_file, "-petscopts", "--petscopts", "PetscOptions file to use.");
    args.AddOption(&config_file, "-config", "--config", "Simulation configuration input file.");
    args.Parse();
    args.ParseCheck();

    MFEMInitializePetsc(&argc, &argv, petscrc_file, nullptr);

    {
        nse::ParallelTimer t_solve(MPI_COMM_WORLD, "solve");
        nse::ParallelTimer t_total(MPI_COMM_WORLD, "total");

        t_total.Start();
        {
            nse::InputData idata;
            idata.ReadFromFile(config_file, petscrc_file);

            FEMachinery fem(idata);
            TimeLevelFields tlf(fem);

            ProblemCase *pcase = nullptr; // polymorphic
            SelectProblemCase(idata, fem, tlf, pcase);

            SolverNS *ns_solver;
            if (idata.method_config.is_coupled()) {
                ns_solver = new SolverNSCoupled(idata, fem, tlf, pcase);
            } else if (idata.method_config.is_uncoupled()) {

            }

            TimeStepper time_stepper(idata, fem, tlf, pcase, ns_solver);

            t_solve.Start();
            time_stepper.SetInitialConditions();
            time_stepper.RegisterParaviewFields();
            time_stepper.WriteInitialConditions();
            time_stepper.MarchInTime();
            t_solve.Stop();

            delete pcase;
        }
        t_total.Stop();

        t_solve.Print();
        t_total.Print();
    }
    // We finalize PETSc
    // Hypre::Finalize();
    MFEMFinalizePetsc();

    return 0;
}
