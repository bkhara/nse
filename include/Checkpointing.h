//
// Created by khara on 2/10/26.
//

#pragma once

#include "mfem.hpp"
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <sys/stat.h>   // mkdir
#include <sys/types.h>

using namespace mfem;

namespace nse {
    class InputData;
    class FEMachinery;
    class TimeLevelFields;
    struct RestartMeta {
        int step = 0;
        double t = 0.0;
        double dt = 0.0;
        int comm_size = 1;
    };

    class Checkpointing {
    public:
        // mkdir that only rank 0 calls; everyone waits.
        static void MakeCheckpointDir(const std::string &dir, MPI_Comm comm);

        static void MkdirP(const std::string &path, MPI_Comm comm);

        static void RemoveAllInDir(const std::string &path, MPI_Comm comm);

        static std::string RankSuffix(MPI_Comm comm, int width = 6);

        static std::string StepDir(const std::string &root, int step);

        static void WriteRestartMeta(const std::string &fname,
                                     const RestartMeta &m, MPI_Comm comm);

        static RestartMeta ReadRestartMeta(const std::string &fname,
                                           MPI_Comm comm);

        // ---------------------- A/B rotation helpers ----------------------

        static char SlotForStep(int step);

        static char SlotForStep(int step, int chk_every);

        static char OtherSlot(char s);

        static void WriteCurrentSlot(const std::string &root, char slot, MPI_Comm comm);

        static char ReadCurrentSlot(const std::string &root, MPI_Comm comm);

        static bool HasDoneMarker(const std::string &dir, MPI_Comm comm);

        static void WriteDoneMarker(const std::string &dir, MPI_Comm comm);

        static void WriteCurrentCycle(const std::string &dir, int cycle, MPI_Comm comm);
        static int ReadCurrentCycle(const std::string& dir);

        // --------- MAIN CHECKPOINT WRITE ---------

        static void WriteParallelMesh(const std::string &dir, mfem::ParMesh &pmesh, MPI_Comm comm);

        static mfem::ParMesh *ReadParallelMesh(const std::string &dir, MPI_Comm comm);


        static void WriteCheckpointPlain(InputData &idata,
                                    FEMachinery &fem,
                                    TimeLevelFields &tlf,
                                    const RestartMeta &meta,
                                    MPI_Comm comm);

        // --------- MAIN CHECKPOINT READ ---------

        static char ChooseValidSlot(const std::string &root, MPI_Comm comm);

        static void ReadCheckpointMeshPlain(InputData &idata,
                                       ParMesh *&pmesh,
                                       MPI_Comm comm);

        static void ReadCheckpointPlain(InputData &idata,
                                   TimeLevelFields &tlf,
                                   RestartMeta &meta_out,
                                   MPI_Comm comm);

        static void WriteCheckpointParaView(InputData &idata, FEMachinery &fem, TimeLevelFields &tlf, const RestartMeta &meta, MPI_Comm comm);
        static void ReadCheckpointParaView(InputData &idata, FEMachinery& fem, TimeLevelFields &tlf, RestartMeta &meta_out, MPI_Comm comm);

        static void ReadCheckpointMeshVisIt(InputData &idata,
                                            ParMesh*& mesh,
                                            MPI_Comm comm);

        static void ReadCheckpointVisIt(InputData &idata,
                                        std::string &INPUT_DIR,
                                        FEMachinery &fem,
                                        TimeLevelFields &tlf,
                                        RestartMeta &meta_out,
                                        MPI_Comm comm);

        static void WriteCheckpointVisIt(InputData &idata, FEMachinery &fem,
                                         TimeLevelFields &tlf,
                                         const RestartMeta &meta,
                                         MPI_Comm comm);

        static void WriteCheckpoint(InputData &idata,
                             FEMachinery &fem,
                             TimeLevelFields &tlf,
                             const RestartMeta &meta,
                             MPI_Comm comm);

        static void ReadCheckpointMesh(InputData &idata,
                                ParMesh *&pmesh,
                                MPI_Comm comm);

        static void ReadCheckpointData(InputData &idata,
                            FEMachinery &fem,
                            TimeLevelFields &tlf,
                            RestartMeta &meta_out,
                            MPI_Comm comm);

    };
}
