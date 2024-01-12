// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT
#include "CAprint.hpp"

#include "mpi.h"

// Print timing info to console
void printExaCATiming(const int np, const double init_time, const double run_time, const double out_time,
                      const int cycle, const double init_max_time, const double init_min_time,
                      const double nucl_max_time, const double nucl_min_time, const double create_sv_min_time,
                      const double create_sv_max_time, const double capture_max_time, const double capture_min_time,
                      const double ghost_max_time, const double ghost_min_time, const double out_max_time,
                      const double out_min_time) {

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
