//
// Created by Biswajit Khara on 10/22/25.
//

#pragma once

#include <mpi.h>
#include <string>
#include <iostream>
#include <iomanip>

namespace nse {
    class ParallelTimer {
    private:
        MPI_Comm comm_;
        int rank_ = 0, size_ = 1;
        std::string name_;

        mutable double t0_ = 0.0;
        mutable bool   running_ = false;
        mutable double total_ = 0.0;
        mutable double last_time_ = 0.0;
        int segments_ = 0;

    public:
        struct Stats {
            double min = 0.0, max = 0.0, sum = 0.0, avg = 0.0;
            int    nprocs = 1;
        };

        explicit ParallelTimer(MPI_Comm comm = MPI_COMM_WORLD, std::string name = {})
            : comm_(comm), name_(std::move(name)) { MPI_Comm_rank(comm_, &rank_); MPI_Comm_size(comm_, &size_); }

        // Start/stop accumulate into total_
        void Start() {
            if (!running_) { t0_ = MPI_Wtime(); running_ = true; }
        }
        void Stop() {
            double end_time = MPI_Wtime();
            if (running_) {
                last_time_ = end_time - t0_;
                total_ += last_time_;
                running_ = false;
                ++segments_;
            }
        }

        // Reset the accumulated time (doesn't change name/comm)
        void Reset() { total_ = 0.0; running_ = false; segments_ = 0; }

        // Local (per-rank) elapsed totals; if running, includes the open segment up to "now".
        double LocalElapsed() const {
            return running_ ? (total_ + (MPI_Wtime() - t0_)) : total_;
        }

        int Segments() const { return segments_; }

        // Optional barriered section timing helper:
        //   double dt = timer.TimeThis([]{ work(); }, /*barrier=*/true);
        template <typename F>
        double TimeThis(F&& f, bool barrier = true) {
            if (barrier) MPI_Barrier(comm_);
            double s = MPI_Wtime();
            std::forward<F>(f)();
            if (barrier) MPI_Barrier(comm_);
            double dt = MPI_Wtime() - s;
            total_ += dt;
            ++segments_;
            return dt;
        }

        // Reduce the provided value (defaults to the timer's LocalElapsed) to get min/avg/max/sum.
        Stats Reduce(bool barrier = true,
                     double value = std::numeric_limits<double>::quiet_NaN()) const {
            if (barrier) MPI_Barrier(comm_);
            const double v = std::isnan(value) ? LocalElapsed() : value;

            Stats st;
            MPI_Allreduce(&v, &st.min, 1, MPI_DOUBLE, MPI_MIN, comm_);
            MPI_Allreduce(&v, &st.max, 1, MPI_DOUBLE, MPI_MAX, comm_);
            MPI_Allreduce(&v, &st.sum, 1, MPI_DOUBLE, MPI_SUM, comm_);
            st.nprocs = size_;
            st.avg = st.sum / st.nprocs;
            return st;
        }

        // Print on rank 0.
        void Print(const Stats& st, const std::string& label = "") const {
            if (rank_ == 0) {
                const std::string& nm = label.empty() ? name_ : label;
                std::cout << std::fixed << std::setprecision(6);
                if (!nm.empty()) std::cout << "[" << nm << "] ";
                std::cout << "n="   << st.nprocs << ", "
                          << "min=" << st.min << " sec, "
                          << "max=" << st.max << " sec, "
                          << "sum=" << st.sum << " sec, "
                          << "avg=" << st.avg << " sec"
                          << '\n';
            }
        }

        void Print(const std::string& label = "") const {
            auto st = this->Reduce(true);
            Print(st, label);
        }

        double GetLastTimeSeconds() const {
            return last_time_;
        }

        // Small RAII scope helper: starts on construction, stops on destruction, accumulates into this timer.
        class Scope {
        public:
            explicit Scope(ParallelTimer& T, bool barrier = false) : T_(T), barrier_(barrier) {
                if (barrier_) MPI_Barrier(T_.comm_); T_.Start();
            }
            ~Scope() {
                T_.Stop(); if (barrier_) MPI_Barrier(T_.comm_);
            }
        private:
            ParallelTimer& T_;
            bool barrier_;
        };
    };
}