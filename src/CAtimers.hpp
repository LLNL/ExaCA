// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_TIMERS_HPP
#define EXACA_TIMERS_HPP

#include "mpi.h"
#include <iostream>
#include <sstream>

// Print timing info
struct Timers {

    int id;
    double init_time = 0.0, nucl_time = 0.0, create_sv_time = 0.0, capture_time = 0.0, ghost_time = 0.0;
    double start_nucl_time, start_create_sv_time, start_capture_time, start_ghost_time;
    double start_init_time, start_run_time, start_out_time;
    double run_time, out_time;
    double layer_time_1, layer_time_2;

    double init_max_time, init_min_time, out_max_time, out_min_time = 0.0;
    double nucl_max_time, nucl_min_time, create_sv_min_time, create_sv_max_time, capture_max_time, capture_min_time,
        ghost_max_time, ghost_min_time = 0.0;

    Timers(const int mpi_id)
        : id(mpi_id) {}

    void startInit() { start_init_time = MPI_Wtime(); }
    void endInit() {
        init_time = MPI_Wtime() - start_init_time;
        if (id == 0)
            std::cout << "Data initialized: Time spent: " << init_time << " s" << std::endl;
    }

    void startRun() { start_run_time = MPI_Wtime(); }

    void endRun() { run_time = MPI_Wtime() - start_run_time; }

    void startOutput() { start_out_time = MPI_Wtime(); }

    void endOutput() { out_time = MPI_Wtime() - start_out_time; }

    void startNucleation() { start_nucl_time = MPI_Wtime(); }
    void endNucleation() { nucl_time += MPI_Wtime() - start_nucl_time; }

    void startSV() { start_create_sv_time = MPI_Wtime(); }
    void endSV() { create_sv_time += MPI_Wtime() - start_create_sv_time; }

    void startCapture() { start_capture_time = MPI_Wtime(); }
    void endCapture() { capture_time += MPI_Wtime() - start_capture_time; }

    void startComm() { start_ghost_time = MPI_Wtime(); }
    void endComm() { ghost_time += MPI_Wtime() - start_ghost_time; }

    void startLayer() { layer_time_1 = MPI_Wtime(); }
    void endLayer(const int layernumber) {
        double layer_time_2 = MPI_Wtime();
        if (id == 0)
            std::cout << "Time for layer number " << layernumber << " was " << layer_time_2 - layer_time_1
                      << " s, starting layer " << layernumber + 1 << std::endl;
    }
    void endLayer() {
        double layer_time_2 = MPI_Wtime();
        if (id == 0)
            std::cout << "Time for final layer was " << layer_time_2 - layer_time_1 << " s" << std::endl;
    }

    double getTotal() { return init_time + run_time + out_time; }

    auto printLog() {
        // This assumes reduceMPI() has already been called.
        std::stringstream log;
        log << "   \"Timing\": {" << std::endl;
        log << "       \"Runtime\": " << getTotal() << "," << std::endl;
        log << "       \"InitRunOutputBreakdown\": [" << init_time << "," << run_time << "," << out_time << "],"
            << std::endl;
        log << "       \"MaxMinInitTime\": [" << init_max_time << "," << init_min_time << "]," << std::endl;
        log << "       \"MaxMinNucleationTime\": [" << nucl_max_time << "," << nucl_min_time << "]," << std::endl;
        log << "       \"MaxMinSteeringVectorCreationTime\": [" << create_sv_max_time << "," << create_sv_min_time
            << "]," << std::endl;
        log << "       \"MaxMinCellCaptureTime\": [" << capture_max_time << "," << capture_min_time << "],"
            << std::endl;
        log << "       \"MaxMinGhostExchangeTime\": [" << ghost_max_time << "," << ghost_min_time << "]," << std::endl;
        log << "       \"MaxMinOutputTime\": [" << out_max_time << "," << out_min_time << "]" << std::endl;
        log << "   }" << std::endl;
        return log.str();
    }

    void reduceMPI() {
        // Reduce all times across MPI ranks
        MPI_Allreduce(&init_time, &init_max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&init_time, &init_min_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&nucl_time, &nucl_max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&nucl_time, &nucl_min_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&create_sv_time, &create_sv_max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&create_sv_time, &create_sv_min_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&capture_time, &capture_max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&capture_time, &capture_min_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&ghost_time, &ghost_max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&ghost_time, &ghost_min_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&out_time, &out_max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&out_time, &out_min_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    }

    void printFinal(const int np, const int cycle) {

        if (id != 0)
            return;

        std::cout << "===================================================================================" << std::endl;
        std::cout << "Having run with = " << np << " processors" << std::endl;
        std::cout << "Output written at cycle = " << cycle << std::endl;
        std::cout << "Total time = " << init_time + run_time + out_time << std::endl;
        std::cout << "Time spent initializing data = " << init_time << " s" << std::endl;
        std::cout << "Time spent performing CA calculations = " << run_time << " s" << std::endl;
        std::cout << "Time spent collecting and printing output data = " << out_time << " s\n" << std::endl;

        std::cout << "Max/min rank time initializing data  = " << init_max_time << " / " << init_min_time << " s"
                  << std::endl;
        std::cout << "Max/min rank time in CA nucleation   = " << nucl_max_time << " / " << nucl_min_time << " s"
                  << std::endl;
        std::cout << "Max/min rank time in CA steering vector creation = " << create_sv_max_time << " / "
                  << create_sv_min_time << " s" << std::endl;
        std::cout << "Max/min rank time in CA cell capture = " << capture_max_time << " / " << capture_min_time << " s"
                  << std::endl;
        std::cout << "Max/min rank time in CA ghosting     = " << ghost_max_time << " / " << ghost_min_time << " s"
                  << std::endl;
        std::cout << "Max/min rank time exporting data     = " << out_max_time << " / " << out_min_time << " s\n"
                  << std::endl;

        std::cout << "===================================================================================" << std::endl;
    }
};

#endif
