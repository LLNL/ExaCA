// Copyright Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_TIMERS_HPP
#define EXACA_TIMERS_HPP

#include "mpi.h"
#include <iostream>
#include <sstream>

class Timer {
    double _time = 0.0;
    double start_time = 0.0;
    double max_time, min_time = 0.0;
    int num_calls = 0;

  public:
    void start() { start_time = MPI_Wtime(); }
    void stop() {
        _time += MPI_Wtime() - start_time;
        num_calls++;
    }
    void reset() { _time = 0.0; }
    auto time() { return _time; }
    auto numCalls() { return num_calls; }
    auto minTime() { return max_time; }
    auto maxTime() { return min_time; }

    void reduceMPI() {
        MPI_Allreduce(&_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    }

    auto print(std::string description) {
        std::stringstream out;
        out << "Time spent " << description << " = " << _time << " s" << std::endl;
        return out.str();
    }

    auto printMinMax(std::string description) {
        std::stringstream out;
        out << "Max/min rank time " << description << " = " << max_time << " / " << min_time << " s" << std::endl;
        return out.str();
    }
};

// Print timing info
struct Timers {

    int id;
    Timer init, run, output;
    Timer nucl, create_sv, capture, ghost;
    Timer layer;
    Timer heat_transfer;

    Timers(const int mpi_id)
        : id(mpi_id)
        , init()
        , run()
        , output()
        , nucl()
        , create_sv()
        , capture()
        , ghost()
        , layer()
        , heat_transfer() {}

    void startInit() { init.start(); }
    void stopInit() {
        init.stop();
        if (id == 0)
            std::cout << "Data initialized: Time spent: " << init.time() << " s" << std::endl;
    }

    void startRun() { run.start(); }

    void stopRun() { run.stop(); }

    void startOutput() { output.start(); }

    void stopOutput() { output.stop(); }

    void startNucleation() { nucl.start(); }
    void stopNucleation() { nucl.stop(); }

    void startSV() { create_sv.start(); }
    void stopSV() { create_sv.stop(); }

    void startCapture() { capture.start(); }
    void stopCapture() { capture.stop(); }

    void startComm() { ghost.start(); }
    void stopComm() { ghost.stop(); }

    void startLayer() { layer.start(); }
    void stopLayer(const int layernumber) {
        layer.stop();
        if (id == 0)
            std::cout << "Time for layer number " << layernumber << " was " << layer.time() << " s, starting layer "
                      << layernumber + 1 << std::endl;
        layer.reset();
    }
    void stopLayer() {
        layer.stop();
        if (id == 0)
            std::cout << "Time for final layer was " << layer.time() << " s" << std::endl;
    }

    void startHeatTransfer() { heat_transfer.start(); }
    void stopHeatTransfer() { heat_transfer.stop(); }

    double getTotal() { return init.time() + run.time() + output.time(); }

    auto printLog() {
        // This assumes reduceMPI() has already been called.
        std::stringstream log;
        log << "   \"Timing\": {" << std::endl;
        log << "       \"Runtime\": " << getTotal() << "," << std::endl;
        log << "       \"InitRunOutputBreakdown\": [" << init.time() << "," << run.time() << "," << output.time()
            << "]," << std::endl;
        log << "       \"MaxMinInitTime\": [" << init.maxTime() << "," << init.minTime() << "]," << std::endl;
        log << "       \"MaxMinNucleationTime\": [" << nucl.maxTime() << "," << nucl.minTime() << "]," << std::endl;
        log << "       \"MaxMinSteeringVectorCreationTime\": [" << create_sv.maxTime() << "," << create_sv.minTime()
            << "]," << std::endl;
        log << "       \"MaxMinCellCaptureTime\": [" << capture.maxTime() << "," << capture.minTime() << "],"
            << std::endl;
        log << "       \"MaxMinGhostExchangeTime\": [" << ghost.maxTime() << "," << ghost.minTime() << "],"
            << std::endl;
        log << "       \"MaxMinOutputTime\": [" << output.maxTime() << "," << output.minTime() << "]" << std::endl;
        log << "   }" << std::endl;
        return log.str();
    }

    void reduceMPI() {
        // Reduce all times across MPI ranks
        init.reduceMPI();
        nucl.reduceMPI();
        create_sv.reduceMPI();
        capture.reduceMPI();
        ghost.reduceMPI();
        output.reduceMPI();

        if (heat_transfer.time() > 0)
            heat_transfer.reduceMPI();
    }

    void printFinal(const int np, const int cycle) {

        if (id != 0)
            return;

        std::cout << "===================================================================================" << std::endl;
        std::cout << "Having run with = " << np << " processors" << std::endl;
        std::cout << "Output written at cycle = " << cycle << std::endl;
        std::cout << "Total time = " << getTotal() << std::endl;
        if (heat_transfer.time() > 0)
            std::cout << heat_transfer.print("performing heat transfer simulation");
        std::cout << init.print("initializing data");
        std::cout << run.print("performing CA calculations");
        std::cout << output.print("collecting and printing output data");

        std::cout << init.printMinMax("initializing data");
        std::cout << nucl.printMinMax("in CA nucleation");
        std::cout << create_sv.printMinMax("in CA steering vector creation");
        std::cout << capture.printMinMax("in CA cell capture");
        std::cout << ghost.printMinMax("in CA cell communication");
        std::cout << output.printMinMax("exporting data");

        std::cout << "===================================================================================" << std::endl;
    }
};

#endif
