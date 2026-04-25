//
// Created by Biswajit Khara on 2/11/26.
//

#include "InputData.h"
#include "FEMachinery.h"
#include "TimeLevelFields.h"
#include "Checkpointing.h"

namespace nse {
    void Checkpointing::MkdirP(const std::string &path, MPI_Comm comm) {
        int rank;
        MPI_Comm_rank(comm, &rank);

        if (rank == 0) {
            std::error_code ec;
            std::filesystem::create_directories(path, ec); // mkdir -p
            MFEM_VERIFY(!ec, "create_directories failed for " << path << ": " << ec.message());
        }

        MPI_Barrier(comm);
    }

    void Checkpointing::RemoveAllInDir(const std::string &path, MPI_Comm comm) {
        int rank;
        MPI_Comm_rank(comm, &rank);
        if (rank == 0) {
            std::error_code ec;
            std::filesystem::remove_all(path, ec); // ok if missing
            std::filesystem::create_directories(path, ec);
            MFEM_VERIFY(!ec, "remove_all/create_directories failed for " << path << ": " << ec.message());
        }
        MPI_Barrier(comm);
    }

    std::string Checkpointing::RankSuffix(MPI_Comm comm, int width) {
        int rank;
        MPI_Comm_rank(comm, &rank);
        std::ostringstream ss;
        ss << "." << std::setw(width) << std::setfill('0') << rank;
        return ss.str();
    }

    void Checkpointing::WriteRestartMeta(const std::string &fname, const RestartMeta &m, MPI_Comm comm) {
        int rank;
        MPI_Comm_rank(comm, &rank);
        if (rank != 0) { return; }

        std::ofstream os(fname);
        MFEM_VERIFY(os, "Cannot open meta file: " << fname);

        os << "step " << m.step << "\n";
        os << std::setprecision(17);
        os << "t " << m.t << "\n";
        os << "dt " << m.dt << "\n";
        os << "comm_size " << m.comm_size << "\n";
    }

    RestartMeta Checkpointing::ReadRestartMeta(const std::string &fname, MPI_Comm comm) {
        RestartMeta m;
        int rank;
        MPI_Comm_rank(comm, &rank);

        if (rank == 0) {
            std::ifstream is(fname);
            MFEM_VERIFY(is, "Cannot open meta file: " << fname);
            std::string key;
            while (is >> key) {
                if (key == "step") { is >> m.step; }
                else if (key == "t") { is >> m.t; }
                else if (key == "dt") { is >> m.dt; }
                else if (key == "comm_size") { is >> m.comm_size; }
                else {
                    std::string dummy;
                    is >> dummy; // tolerate extras
                }
            }
        }

        int ibuf[2] = {m.step, m.comm_size};
        MPI_Bcast(ibuf, 2, MPI_INT, 0, comm);
        m.step = ibuf[0];
        m.comm_size = ibuf[1];

        double dbuf[2] = {m.t, m.dt};
        MPI_Bcast(dbuf, 2, MPI_DOUBLE, 0, comm);
        m.t  = dbuf[0];
        m.dt = dbuf[1];

        return m;
    }

    char Checkpointing::SlotForStep(int step, int chk_every) {
        const int chk_id = step / chk_every; // 0,1,2,3,... at each checkpoint
        return (chk_id % 2 == 0) ? 'A' : 'B';
    }

    char Checkpointing::OtherSlot(char s) { return (s == 'A') ? 'B' : 'A'; }

    void Checkpointing::WriteCurrentSlot(const std::string &root, char slot, MPI_Comm comm) {
        int rank;
        MPI_Comm_rank(comm, &rank);
        if (rank != 0) { return; }
        std::ofstream os(root + "/CURRENT");
        MFEM_VERIFY(os, "Cannot write CURRENT file");
        os << slot << "\n";
    }

    char Checkpointing::ReadCurrentSlot(const std::string &root, MPI_Comm comm) {
        int rank;
        MPI_Comm_rank(comm, &rank);
        char slot = 'A';
        if (rank == 0) {
            std::ifstream is(root + "/CURRENT");
            if (is) { is >> slot; }
            // if missing, default to 'A' (caller should validate DONE anyway)
        }
        MPI_Bcast(&slot, 1, MPI_CHAR, 0, comm);
        return slot;
    }

    bool Checkpointing::HasDoneMarker(const std::string &dir, MPI_Comm comm) {
        int rank;
        MPI_Comm_rank(comm, &rank);
        int ok = 0;
        if (rank == 0) {
            ok = std::filesystem::exists(dir + "/DONE") ? 1 : 0;
        }
        MPI_Bcast(&ok, 1, MPI_INT, 0, comm);
        return ok == 1;
    }

    void Checkpointing::WriteDoneMarker(const std::string &dir, MPI_Comm comm) {
        int rank;
        MPI_Comm_rank(comm, &rank);
        if (rank != 0) { return; }
        std::ofstream os(dir + "/DONE");
        MFEM_VERIFY(os, "Cannot write DONE marker");
        os << "ok\n";
    }

    void Checkpointing::WriteCurrentCycle(const std::string &dir, int cycle, MPI_Comm comm) {
        int rank;
        MPI_Comm_rank(comm, &rank);
        if (rank != 0) { return; }
        std::ofstream os(dir + "/CYCLE");
        MFEM_VERIFY(os, "Cannot write CYCLE file");
        os << cycle << "\n";
    }

    int Checkpointing::ReadCurrentCycle(const std::string& dir) {
        std::string filename = dir + "/CYCLE";
        std::ifstream in(filename);
        if (!in) { throw std::runtime_error("Could not open: " + filename); }

        int cycle;
        if (!(in >> cycle)) { throw std::runtime_error("Failed to read: " + filename); }
        return cycle;
    }

    void Checkpointing::WriteParallelMesh(const std::string &dir, mfem::ParMesh &pmesh, MPI_Comm comm) {
        int rank; MPI_Comm_rank(comm, &rank);

        const std::string rs = RankSuffix(comm);

        std::ostringstream fn;
        fn << dir << "/pmesh" << rs;

        std::ofstream os(fn.str());
        MFEM_VERIFY(os, "Cannot open parallel mesh file for write: " << fn.str());

        pmesh.ParPrint(os);  // <-- THIS is the key for true parallel mesh format
    }

    mfem::ParMesh * Checkpointing::ReadParallelMesh(const std::string &dir, MPI_Comm comm) {
        int rank; MPI_Comm_rank(comm, &rank);

        const std::string rs = RankSuffix(comm);

        std::ostringstream fn;
        fn << dir << "/pmesh" << rs;

        std::ifstream is(fn.str());
        MFEM_VERIFY(is, "Cannot open parallel mesh file for read: " << fn.str());

        return new mfem::ParMesh(comm, is);  // expects "communication_groups" etc.
    }

    /**
     * @brief Manually loads a QuadratureFunction from an MFEM-formatted file.
     * This bypasses the header and forces a raw vector load to avoid silent failures.
     */
    void LoadQuadratureFunc(std::istream &is, mfem::QuadratureFunction &qf) {
        std::string line;
        bool found_data = false;

        // 1. Advance the stream to the end of the header
        // MFEM headers for QFields end after the "VDim: X" line
        while (std::getline(is, line)) {
            if (line.find("VDim:") != std::string::npos) {
                found_data = true;
                break;
            }
        }

        if (!found_data) {
            MFEM_ABORT("LoadQuadratureSafe: Could not find end of header (VDim line)!");
        }

        // 2. The stream head is now at the newline after "VDim: X".
        // We skip the next line if it's empty to get to the raw numbers.
        is >> std::ws;

        // 3. Perform a 'Force Load'.
        // We use the Vector parent class Load to avoid the 'int s; in >> s' logic
        // that fails on text. We tell it exactly how many doubles to expect.
        qf.Vector::Load(is, qf.Size());

        // 4. Verification
        if (is.fail()) {
            MFEM_ABORT("LoadQuadratureSafe: Stream failed during data extraction. "
                "Check if file size matches QuadratureFunction size (" << qf.Size() << ").");
        }
    }

    char Checkpointing::ChooseValidSlot(const std::string &root, MPI_Comm comm) {
        // Prefer CURRENT; if invalid, fall back to the other slot if it has DONE.
        char cur = ReadCurrentSlot(root, comm);
        std::string curdir = root + "/" + std::string(1, cur);
        if (HasDoneMarker(curdir, comm)) { return cur; }

        char other = OtherSlot(cur);
        std::string otherdir = root + "/" + std::string(1, other);
        MFEM_VERIFY(HasDoneMarker(otherdir, comm),
                    "No valid checkpoint found: neither slot has DONE");
        return other;
    }

    void Checkpointing::ReadCheckpointVisIt(InputData &idata,
                                            std::string &INPUT_DIR,
                                            FEMachinery &fem,
                                            TimeLevelFields &tlf,
                                            RestartMeta &meta_out,
                                            MPI_Comm comm) {
        const std::string dir = INPUT_DIR;

        // meta first (so caller can set step/t/dt)
        meta_out = ReadRestartMeta(dir + "/restart.meta", comm);

        int cycle = ReadCurrentCycle(dir);

        // ---- Load mesh + fields ----
        // Use a stable collection name (NOT the path). Keep it the same in writer and reader.
        const char *dc_name = "chkpt";

        mfem::VisItDataCollection dc(comm, dc_name, fem.mesh);
        dc.SetPrefixPath(dir);
        dc.SetCycle(cycle);
        dc.Load(cycle);

        std::string visit_folder = dir + "/" + std::string(dc_name) + "_" + std::to_string(cycle).insert(0, 6 - std::to_string(cycle).length(), '0');

        // Now (re)construct fields on the new spaces if needed, then pull data by name:
        // Easiest: register pointers to your already-constructed objects and call dc.Load().
        // But since we've already loaded, we can fetch GridFunction pointers from the dc:
        auto &g = tlf.current;

        auto load_GF =[&](const std::string &name, ParGridFunction &gf) {
            mfem::ParGridFunction *v = dc.GetParField(name);
            if (v) {
                gf = *v;
                if (Mpi::Root()) {
                    std::cout << "Read from VisIt checkpoint data | "  << visit_folder << ": " << name << std::endl;
                }
            } else {
                if (Mpi::Root()) {
                    MFEM_VERIFY(v, "Missing field in checkpoint: " + name);
                    std::cout << "ParGridFunction " << name << ": checkpoint file loading ERROR (File not found or permission denied)" << std::endl;
                }
            }
        };

        auto load_Qf = [&](const std::string &name, QuadratureFunction &qf) {
            std::string q_file = visit_folder + "/" + name + RankSuffix(comm);
            std::ifstream q_is(q_file);
            // Check if file exists and is accessible
            if (q_is.is_open()) {
                LoadQuadratureFunc(q_is, qf);
                if (Mpi::Root()) {
                    std::cout << "Read from VisIt checkpoint data (custom loader) | " << visit_folder << ": " << name << std::endl;
                }
            } else {
                if (Mpi::Root()) {
                    std::cout << "QFunction " << name << ": checkpoint file loading ERROR (File not found or permission denied)" << std::endl;
                }
            }
        };
        // quadrature functions -- uses my custom LoadQuadrature
        load_Qf("Hq", g.Hq);

        // pargridfunctions -- uses the VisItDataCollection::GetParField
        load_GF("u", g.u);
        load_GF("p", g.p);

        MPI_Barrier(comm);
    }

    void Checkpointing::WriteCheckpointVisIt(InputData &idata, FEMachinery &fem,
                                             TimeLevelFields &tlf,
                                             const RestartMeta &meta,
                                             MPI_Comm comm) {
        std::string root = idata.checkpointing_inputs.chkpt_root;
        MkdirP(root, comm);

        const char slot = SlotForStep(meta.step, idata.checkpointing_inputs.chkpt_write_frequency);
        const std::string dir = root + "/" + std::string(1, slot);

        const std::string rs = RankSuffix(comm);
        (void) rs;

        // Clean slot directory before writing new checkpoint into it
        RemoveAllInDir(dir, comm);

        // write the mesh in plain text format - this works more reliably
        WriteParallelMesh(dir, *fem.mesh, MPI_COMM_WORLD);

        // write the data files
        {
            const char *dc_name = "chkpt";
            mfem::VisItDataCollection dc(comm, dc_name, fem.mesh);
            dc.SetPrefixPath(dir);
            dc.SetCycle(meta.step);     // e.g. 42 -> chkpt_000042
            dc.SetTime(meta.t);         // optional, but nice

            auto &g = tlf.current;

            dc.RegisterField("u", &g.u);
            dc.RegisterField("p", &g.p);

            dc.Save(); // writes mesh + fields in a restart-loadable format
        }

        // metadata (rank 0)
        WriteRestartMeta(dir + "/restart.meta", meta, comm);

        MPI_Barrier(comm);
        WriteDoneMarker(dir, comm);
        WriteCurrentCycle(dir, meta.step, comm);
        WriteCurrentSlot(root, slot, comm);
        MPI_Barrier(comm);
    }

    void Checkpointing::WriteCheckpoint(InputData &idata, FEMachinery &fem, TimeLevelFields &tlf,
        const RestartMeta &meta, MPI_Comm comm) {
        WriteCheckpointVisIt(idata, fem, tlf, meta, comm);
        // WriteCheckpointPlain(idata, fem, tlf, meta, comm);
    }

    void Checkpointing::ReadCheckpointMesh(InputData &idata, ParMesh *&pmesh, MPI_Comm comm) {
        // ReadCheckpointMeshVisIt(idata, pmesh, comm);
        int size_now;
        MPI_Comm_size(comm, &size_now);
        RestartMeta m = ReadRestartMeta(idata.checkpointing_inputs.chkpt_dir_in + "/restart.meta", comm);
        MFEM_VERIFY(
            m.comm_size == size_now,
            "Restart was written with " << m.comm_size
            << " MPI ranks, but current run uses " << size_now
            << ". Restart with a matching number of ranks.");
        pmesh = ReadParallelMesh(idata.checkpointing_inputs.chkpt_dir_in, comm);

        // quick check
        // {
        //     H1_FECollection fec(idata.fes_config_primal.order, pmesh->Dimension());
        //     ParFiniteElementSpace fes(pmesh, &fec);
        //
        //     ParGridFunction temp(&fes);
        //     temp = 0.;
        //
        //     ParaViewDataCollection pvdc("loadedMesh", pmesh);
        //     pvdc.RegisterField("temp", &temp);
        //     pvdc.SetCycle(0);
        //     pvdc.SetTime(0);
        //     pvdc.Save();
        // }
    }

    void Checkpointing::ReadCheckpointData(InputData &idata, FEMachinery &fem, TimeLevelFields &tlf, RestartMeta &meta_out,
                                           MPI_Comm comm) {
        ReadCheckpointVisIt(idata, idata.checkpointing_inputs.chkpt_dir_in, fem, tlf, meta_out, comm);
        // ReadCheckpointPlain(idata, tlf, meta_out, comm);
        {
            ParaViewDataCollection pvdc("LoadedData", fem.mesh);
            pvdc.RegisterField("loaded_u", &tlf.current.u);
            pvdc.RegisterField("loaded_p", &tlf.current.p);
            pvdc.SetCycle(meta_out.step);
            pvdc.SetTime(meta_out.t);
            pvdc.Save();
            ParaViewDataCollection pvdc_q("LoadedData_q", fem.mesh);
            pvdc_q.SetCycle(meta_out.step);
            pvdc_q.SetTime(meta_out.t);
            pvdc_q.Save();
        }
    }
}
