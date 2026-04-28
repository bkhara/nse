//
// Created by Biswajit Khara on 11/25/25.
//

#pragma once

#include "mfem.hpp"
#include "InputData.h"
#include "FEMachinery.h"
#include "TimeLevelFields.h"
#include "ParallelTimer.h"

using namespace mfem;
using std::endl;

namespace nse {
    class TimeStepper {
        int myrank;
        InputData &idata;
        FEMachinery &fem;
        TimeLevelFields &tlf;
        ProblemCase *pcase;

        NSSolver *ns_solver = nullptr;

        ParaViewDataCollection pvdc;
        ParaViewDataCollection pvdc_q;

        ParallelTimer t_ns;
        ParallelTimer t_pv_writer;

        socketstream vis;

        std::string stag_history_filename = "stag_history.txt";

        RestartMeta restart_meta;

    public:
        TimeStepper(InputData &idata, FEMachinery &fem, TimeLevelFields &tlf, ProblemCase *pcase, NSSolver *ns_solver)
            : myrank(Mpi::WorldRank()), idata(idata), fem(fem), tlf(tlf), pcase(pcase), ns_solver(ns_solver),
              pvdc("pv", fem.mesh), pvdc_q("pvq", fem.mesh),
              t_ns(MPI_COMM_WORLD, "solve_elasticity"),
              t_pv_writer(MPI_COMM_WORLD, "pv_writer") {
            if (!myrank) {
                std::ofstream stag_history_file(stag_history_filename.c_str());
                stag_history_file << "t,stag_it,dc_l2,dc_h1,dc_h1semi,dc_infty\n";
                stag_history_file.close();
            }

            // std::stringstream oss;
            // oss << "rank_" << std::setfill('0') << std::setw(2) << myrank;
            // std::ofstream ofs(oss.str().c_str());
            // fracture::PrintElementWiseIntegrationRule(fem.mesh, fem.fespace_primal_c, 2*1, ofs);
            if (idata.streaming_inputs.if_stream) {
                vis.open(idata.streaming_inputs.host.c_str(), idata.streaming_inputs.socket);
            }
        }

        void MarchInTime() {
            double stag_ctol = idata.stag_config.etol;
            int stag_max_it = idata.stag_config.max_it;
            int staggered_iters_total = 0;

            int n = 0;
            double t = 0.0;
            double dt = idata.time_marching.dt;
            if (idata.checkpointing_inputs.restart_from_chkpt) {
                n = restart_meta.step;
                t = restart_meta.t;
                dt = restart_meta.dt;
                if (!Mpi::WorldRank()) {
                    std::cout << "Starting values: dt=" << dt << ", step_0=" << n << ", t_0=" << t << "\n";
                }
            }

            bool is_last_step = false;
            while (t < idata.time_marching.t_max) {
                double time_remaining = (idata.time_marching.t_max - t);
                if (time_remaining < dt or (time_remaining - dt) < idata.time_marching.dt_min) {
                    dt = time_remaining;
                    is_last_step = true;
                    if (!myrank) mfem::out << "dt = " << dt << ". This will be the last step\n";
                }

                if (n == 0) {
                    idata.projection_config.scheme = ProjectionScheme::ChorinFirstOrder;
                } else {
                    idata.projection_config.scheme = ProjectionScheme::IncPressureBDF2;
                }

                n += 1;
                t += dt;
                tlf.SetTime(t);
                tlf.SetTimeStep(dt);
                idata.run_ctx.time_step = n;
                idata.run_ctx.time = t;


                // staggered iteration (or alternating solve iteration)
                double energy_curr = 0;
                double energy_prev = 0;
                double energy_0 = 0;
                double staggered_iters_this_step = 0;
                double stage_timing_ns = 0;

                {
                    t_ns.Start();
                    // el_solver.SolveStep(t, dt);
                    ns_solver->SolveStep(t, dt);
                    t_ns.Stop();

                    stage_timing_ns += t_ns.GetLastTimeSeconds();
                }

                // post processing tasks
                pcase->PostStep(t, dt);
                
                // upate time fields
                tlf.UpdateTimeStepIterates();

                // write paraview files
                t_pv_writer.Start();
                if (n % idata.file_write_freq == 0 or is_last_step) {
                    pvdc.SetCycle(n);
                    pvdc.SetTime(t);
                    pvdc.Save();
                    pvdc_q.SetCycle(n);
                    pvdc_q.SetTime(t);
                    pvdc_q.Save();
                }
                t_pv_writer.Stop();
                if (n % idata.checkpointing_inputs.chkpt_write_frequency == 0 or is_last_step) {
                    restart_meta.comm_size = Mpi::WorldSize();
                    restart_meta.step = n;
                    restart_meta.t = t;
                    restart_meta.dt = dt;
                    nse::Checkpointing::WriteCheckpoint(idata, fem, tlf, restart_meta, fem.mesh->GetComm());
                }
                if (idata.streaming_inputs.if_stream) {
                    nse::StreamingInputs si = idata.streaming_inputs;
                    StreamReopen(*(fem.mesh), tlf.current.p, n, si.frequency, vis);
                }

                // print final info
                if (!myrank) {
                    // std::cout << "[STEP-time::NS](total,sec)=" << stage_timing_ns << ", [STAG-time::NS](avg,sec)=" << stage_timing_ns/staggered_iters_this_step << "\n";
                    std::cout << "Step=" << n << ", t=" << t << "\n";
                    std::cout << "\n";
                }
            }

            t_ns.Print();
            t_pv_writer.Print();
        }

        // all fields are usually set to zero when tlf is initialized
        // but here we do case-by-case initialization, if it exists
        void SetInitialConditions() {
            if (idata.checkpointing_inputs.restart_from_chkpt) {
                nse::Checkpointing::ReadCheckpointData(idata, fem, tlf, restart_meta, fem.mesh->GetComm());
                if (!Mpi::WorldRank()) {
                    std::cout << "Reading from checkpoint file successful!\n";
                }
            }
            else {
                pcase->SetIC(tlf.current);
            }
            tlf.UpdateTimeStepIterates();
        }

        // this is optional
        void WriteInitialConditions() {
            if (idata.checkpointing_inputs.restart_from_chkpt) { return; }
            pvdc.SetCycle(0);
            pvdc.SetTime(0);
            pvdc.Save();
            pvdc_q.SetCycle(0);
            pvdc_q.SetTime(0);
            pvdc_q.Save();
        }

        // set the fields to be written in paraview readable files
        void RegisterParaviewFields() {
            // this is required (needed for all paraview writing)
            pvdc.RegisterField("u_h", &tlf.current.u);
            pvdc.RegisterField("p_h", &tlf.current.p);
            // pvdc.RegisterField("dc_minus", &tlf.dc_minus);
            // pvdc.RegisterField("dc_minus_cmltv", &tlf.dc_minus_cmltv);

            // if (fem.el_vdim > 1) {
            //     pvdc.RegisterField("energy_p", &tlf.energy_p);
            //     pvdc.RegisterField("energy_m", &tlf.energy_m);
            //     pvdc.RegisterField("sigma_p", &tlf.stress_p);
            //     pvdc.RegisterField("sigma_m", &tlf.stress_m);
            //     pvdc.RegisterField("Hplus", &tlf.Hplus);
            // }
            // if (fem.el_vdim > 1 && idata.time_marching.is_dynamic()) {
            //     pvdc.RegisterField("v_h", &tlf.current.v);
            //     pvdc.RegisterField("a_h", &tlf.current.a);
            // }
            //
            // if (idata.method_config.use_ie_pg())
            //     pvdc.RegisterField("psi_h", &tlf.current.psi);
            //
            // if (idata.method_config.use_ie_al())
            //     pvdc.RegisterField("lam_h", &tlf.current.lambda);
            //
            // // pvdc.RegisterField("c_viol", &tlf.c_violation);
            // // pvdc.RegisterField("residual_u", &tlf.residual_u);
            // // pvdc.RegisterField("residual_c", &tlf.residual_c);
            //
            //
            // pvdc_q.RegisterQField("Hq", &tlf.current.Hq);
            // pvdc_q.RegisterQField("cq", &tlf.current.c_q);
            //
            //
            // case specific fields
            pcase->RegisterParaviewFields(pvdc, pvdc_q);
        }
    };
}