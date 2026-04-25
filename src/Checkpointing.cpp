//
// Created by Biswajit Khara on 2/11/26.
//

#include "InputData.h"
#include "FEMachinery.h"
#include "TimeLevelFields.h"
#include "Checkpointing.h"

namespace fracture {
    void Checkpointing::MakeCheckpointDir(const std::string &dir, MPI_Comm comm) {
        int rank;
        MPI_Comm_rank(comm, &rank);

        if (rank == 0) {
            // mkdir -p behavior (one level). If you want nested, add your own helper.
            // If it exists, that's fine.
            ::mkdir(dir.c_str(), 0775);
        }
        MPI_Barrier(comm);
    }

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

    std::string Checkpointing::StepDir(const std::string &root, int step) {
        std::ostringstream ss;
        ss << root << "/chkpt_" << std::setw(8) << std::setfill('0') << step;
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

    char Checkpointing::SlotForStep(int step) {
        return (step % 2 == 0) ? 'A' : 'B';
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

    void Checkpointing::WriteCheckpointPlain(InputData &idata, FEMachinery &fem,
                                        TimeLevelFields &tlf, const RestartMeta &meta, MPI_Comm comm) {
        std::string root = idata.checkpointing_inputs.chkpt_root;
        MkdirP(root, comm);

        const char slot = SlotForStep(meta.step, idata.checkpointing_inputs.chkpt_write_frequency);
        const std::string dir = root + "/" + std::string(1, slot);

        const std::string rs = RankSuffix(comm);

        // Clean slot directory before writing new checkpoint into it
        RemoveAllInDir(dir, comm);

        // 1) mesh (rank-local)
        {
            std::ofstream os(dir + "/mesh" + rs);
            MFEM_VERIFY(os, "Cannot open mesh file for write");
            fem.mesh->Print(os);
        }
        WriteParallelMesh(dir, *fem.mesh, comm);

        // 2) ParGridFunctions (rank-local)
        auto save_pgf = [&](const std::string &name, const ParGridFunction &gf) {
            std::ofstream os(dir + "/" + name + rs);
            MFEM_VERIFY(os, "Cannot open " << name << " for write");
            gf.Save(os);
        };
        auto save_Qf = [&](const std::string &name, const QuadratureFunction &qf) {
            std::ofstream os(dir + "/" + name + rs);
            MFEM_VERIFY(os, "Cannot open " << name << " for write");
            qf.Save(os);
        };
        save_pgf("u", tlf.current.u);
        save_pgf("c", tlf.current.c);
        if (idata.time_marching.is_dynamic()) {
            save_pgf("v", tlf.current.v);
            save_pgf("a", tlf.current.a);
        }

        // 4) metadata (rank 0)
        WriteRestartMeta(dir + "/restart.meta", meta, comm);

        MPI_Barrier(comm);
        WriteDoneMarker(dir, comm);
        WriteCurrentSlot(root, slot, comm);
        MPI_Barrier(comm);
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

     void Checkpointing::ReadCheckpointPlain(InputData &idata, TimeLevelFields &tlf, RestartMeta &meta_out, MPI_Comm comm) {
        const std::string dir = idata.checkpointing_inputs.chkpt_root + "/" + std::string(1, idata.checkpointing_inputs.in_slot);
        const std::string rs = RankSuffix(comm);

        // meta first (so caller can set step/t/dt)
        meta_out = ReadRestartMeta(dir + "/restart.meta", comm);

        // IMPORTANT:
        // After reading the mesh, you must ensure fes_u, fes_c, and qspace are
        // defined on *this* mesh in the same way they were during the run.
        // If your code currently builds fes_u/fes_c earlier, refactor so you can rebuild them here.

        // load PGFs
        auto load_pgf = [&](const std::string &name, ParGridFunction &gf) {
            std::stringstream fn;
            fn << dir + "/" + name + rs;
            std::ifstream is(fn.str());
            if (!Mpi::WorldRank()) {
                std::cout << "Reading " << name << " from " << name << ".*" << std::endl;
            }
            MFEM_VERIFY(is, "Cannot open " << name << " for read");
            gf.Load(is);
        };
        auto load_Qf = [&](const std::string &name, QuadratureFunction &qf) {
            std::stringstream fn;
            fn << dir + "/" + name + rs;
            std::ifstream is(fn.str());
            if (!Mpi::WorldRank()) {
                std::cout << "Reading " << name << " from " << name << ".*" << std::endl;
            }
            MFEM_VERIFY(is, "Cannot open " << name << " for read");
            qf.Load(is);
        };

        load_pgf("u", tlf.current.u);
        load_pgf("c", tlf.current.c);
        if (idata.time_marching.is_dynamic()) {
            load_pgf("v", tlf.current.v);
            load_pgf("a", tlf.current.a);
        }

        MPI_Barrier(comm);
    }

    void Checkpointing::WriteCheckpointParaView(InputData &idata, FEMachinery &fem,
                                            TimeLevelFields &tlf, const RestartMeta &meta,
                                            MPI_Comm comm)
    {
        std::string root = idata.checkpointing_inputs.chkpt_root;
        MkdirP(root, comm);

        const char slot = SlotForStep(meta.step, idata.checkpointing_inputs.chkpt_write_frequency);
        const std::string dir = root + "/" + std::string(1, slot);

        const std::string rs = RankSuffix(comm);
        (void)rs;

        // Clean slot directory before writing new checkpoint into it
        RemoveAllInDir(dir, comm);

        // WRITE THE MESH
        {
            WriteParallelMesh(dir, *fem.mesh, comm);
        }

        // read the data files
        const std::string dc_name = "paraview";
        {
            // IMPORTANT:
            // The "collection name" here is the DataCollection name (it affects the .pvdc filename).
            // To match your reader, we use `root` as the name and write into `dir` via PrefixPath.
            mfem::ParaViewDataCollection dc(dc_name, fem.mesh);
            dc.SetPrefixPath(dir);

            // Optional, keep consistent if you care about ParaView output fidelity
            // dc.SetLevelsOfDetail(1);
            // dc.SetHighOrderOutput(true);

            auto &g = tlf.current;

            // ParGridFunctions
            dc.RegisterField("u",      &g.u);
            dc.RegisterField("v",      &g.v);
            dc.RegisterField("a",      &g.a);
            dc.RegisterField("c",      &g.c);
            dc.RegisterField("psi",    &g.psi);
            dc.RegisterField("lambda", &g.lambda);

            // QuadratureFunctions
            dc.RegisterQField("Hq",      &g.Hq);
            dc.RegisterQField("ePlus_q", &g.ePlus_q);
            dc.RegisterQField("c_q",     &g.c_q);

            dc.Save();
        }

        // 4) metadata (rank 0)
        WriteRestartMeta(dir + "/restart.meta", meta, comm);

        MPI_Barrier(comm);
        WriteDoneMarker(dir, comm);
        WriteCurrentSlot(root, slot, comm);
        MPI_Barrier(comm);
    }

    void Checkpointing::ReadCheckpointParaView(InputData &idata, FEMachinery& fem, TimeLevelFields &tlf, RestartMeta &meta_out, MPI_Comm comm) {
        const std::string dir = idata.checkpointing_inputs.chkpt_root + "/" + std::string(1, idata.checkpointing_inputs.in_slot);
        const std::string rs = RankSuffix(comm);

        // meta first (so caller can set step/t/dt)
        meta_out = ReadRestartMeta(dir + "/restart.meta", comm);

        const std::string dc_name = "paraview";

        {
            mfem::ParaViewDataCollection dc(dc_name, fem.mesh);
            dc.SetPrefixPath(dir);

            // Register fields that exist inside tlf (these must be allocated on the right spaces already)
            // Assuming: tlf.current is your FractureGridFields
            auto &g = tlf.current;

            dc.RegisterField("u", &g.u);
            dc.RegisterField("v", &g.v);
            dc.RegisterField("a", &g.a);
            dc.RegisterField("c", &g.c);
            dc.RegisterField("psi", &g.psi);
            dc.RegisterField("lambda", &g.lambda);

            dc.RegisterQField("Hq", &g.Hq);
            dc.RegisterQField("ePlus_q", &g.ePlus_q);
            dc.RegisterQField("c_q", &g.c_q);

            // Actually read from disk into the registered objects
            dc.Load();
        }

        // Optional: a couple sanity checks (cheap and very useful)
        MFEM_VERIFY(tlf.current.u.Size() > 0 || fem.mesh->GetNE() == 0,
                    "ReadCheckpoint: loaded u has Size()==0 (mesh/space mismatch or wrong directory/name).");
        MFEM_VERIFY(tlf.current.Hq.Size() > 0 || fem.mesh->GetNE() == 0,
                    "ReadCheckpoint: loaded Hq has Size()==0 (quad-space mismatch or wrong directory/name).");

        MPI_Barrier(comm);
    }

    void Checkpointing::ReadCheckpointMeshVisIt(InputData &idata,
                                            ParMesh*& mesh,
                                            MPI_Comm comm) {
        const std::string root = idata.checkpointing_inputs.chkpt_root;
        const std::string dir = root + "/" + std::string(1, idata.checkpointing_inputs.in_slot);

        // ---- Load mesh + fields ----
        // Use a stable collection name (NOT the path). Keep it the same in writer and reader.
        const char *dc_name = "chkpt";

        mfem::VisItDataCollection dc(comm, dc_name);
        dc.SetPrefixPath(dir);
        dc.SetCycle(10);
        dc.Load(10); // <- actually implemented (loads mesh + any fields in the checkpoint)

        // // Replace fem.mesh with the loaded mesh (careful with ownership!)
        // auto *loaded_pm = dynamic_cast<mfem::ParMesh *>(dc.GetMesh());
        // MFEM_VERIFY(loaded_pm, "Checkpoint mesh is not a ParMesh.");
        //
        // mesh = new mfem::ParMesh(*loaded_pm);

        mfem::Mesh *m = dc.GetMesh();
        MFEM_VERIFY(m, "VisItDataCollection returned null mesh.");
        if (auto *pm = dynamic_cast<mfem::ParMesh *>(m)) {
            std::cout << "dynamic route worked\n";
            // Parallel checkpoint: clone directly
            mesh = new mfem::ParMesh(*pm);
        } else {
            std::cout << "dynamic route didn't work\n";
            // Serial checkpoint: partition into ParMesh
            // (This requires that m is a real mfem::Mesh.)
            auto *sm = dynamic_cast<mfem::Mesh *>(m);
            MFEM_VERIFY(sm, "Loaded mesh is neither ParMesh nor Mesh.");

            mesh = new mfem::ParMesh(comm, *sm);
        }

        MPI_Barrier(comm);
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

        // pargridfunctions -- uses the VisItDataCollection::GetParField
        load_GF("u", g.u);
        load_GF("v", g.v);
        load_GF("a", g.a);
        load_GF("c", g.c);
        load_GF("psi", g.psi);
        load_GF("lambda", g.lambda);

        // quadrature functions -- uses my custom LoadQuadrature
        load_Qf("Hq", g.Hq);
        load_Qf("ePlus_q", g.ePlus_q);
        load_Qf("c_q", g.c_q);

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
            dc.RegisterField("v", &g.v);
            dc.RegisterField("a", &g.a);
            dc.RegisterField("c", &g.c);
            dc.RegisterField("psi", &g.psi);
            dc.RegisterField("lambda", &g.lambda);

            dc.RegisterQField("Hq", &g.Hq);
            dc.RegisterQField("ePlus_q", &g.ePlus_q);
            dc.RegisterQField("c_q", &g.c_q);

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
            pvdc.RegisterField("loaded_v", &tlf.current.v);
            pvdc.RegisterField("loaded_a", &tlf.current.a);
            pvdc.RegisterField("loaded_c", &tlf.current.c);
            pvdc.RegisterField("loaded_psi", &tlf.current.psi);
            pvdc.SetCycle(meta_out.step);
            pvdc.SetTime(meta_out.t);
            pvdc.Save();
            ParaViewDataCollection pvdc_q("LoadedData_q", fem.mesh);
            pvdc_q.RegisterQField("loaded_Hq", &tlf.current.Hq);
            pvdc_q.RegisterQField("loaded_ePlus_q", &tlf.current.ePlus_q);
            pvdc_q.RegisterQField("loaded_c_q", &tlf.current.c_q);
            pvdc_q.SetCycle(meta_out.step);
            pvdc_q.SetTime(meta_out.t);
            pvdc_q.Save();
        }
    }
}
