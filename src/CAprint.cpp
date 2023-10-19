// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT
#include "CAprint.hpp"

#include "mpi.h"

// Print timing info to console
void PrintExaCATiming(int np, double InitTime, double RunTime, double OutTime, int cycle, double InitMaxTime,
                      double InitMinTime, double NuclMaxTime, double NuclMinTime, double CreateSVMinTime,
                      double CreateSVMaxTime, double CaptureMaxTime, double CaptureMinTime, double GhostMaxTime,
                      double GhostMinTime, double OutMaxTime, double OutMinTime) {

    std::cout << "===================================================================================" << std::endl;
    std::cout << "Having run with = " << np << " processors" << std::endl;
    std::cout << "Output written at cycle = " << cycle << std::endl;
    std::cout << "Total time = " << InitTime + RunTime + OutTime << std::endl;
    std::cout << "Time spent initializing data = " << InitTime << " s" << std::endl;
    std::cout << "Time spent performing CA calculations = " << RunTime << " s" << std::endl;
    std::cout << "Time spent collecting and printing output data = " << OutTime << " s\n" << std::endl;

    std::cout << "Max/min rank time initializing data  = " << InitMaxTime << " / " << InitMinTime << " s" << std::endl;
    std::cout << "Max/min rank time in CA nucleation   = " << NuclMaxTime << " / " << NuclMinTime << " s" << std::endl;
    std::cout << "Max/min rank time in CA steering vector creation = " << CreateSVMaxTime << " / " << CreateSVMinTime
              << " s" << std::endl;
    std::cout << "Max/min rank time in CA cell capture = " << CaptureMaxTime << " / " << CaptureMinTime << " s"
              << std::endl;
    std::cout << "Max/min rank time in CA ghosting     = " << GhostMaxTime << " / " << GhostMinTime << " s"
              << std::endl;
    std::cout << "Max/min rank time exporting data     = " << OutMaxTime << " / " << OutMinTime << " s\n" << std::endl;

    std::cout << "===================================================================================" << std::endl;
}
