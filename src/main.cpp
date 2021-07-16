// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAinitialize.hpp"
#include "runCA.hpp"
#include "runCAremelt.hpp"
#include <Kokkos_Core.hpp>

#include "mpi.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char *argv[]) {
    // Initialize MPI
    int id, np;
    MPI_Init(&argc, &argv);
    // Initialize Kokkos
    Kokkos::initialize(argc, argv);
    {

        // Get number of processes
        MPI_Comm_size(MPI_COMM_WORLD, &np);
        // Get individual process ID
        MPI_Comm_rank(MPI_COMM_WORLD, &id);

        if (id == 0)
            Kokkos::DefaultExecutionSpace::print_configuration(std::cout);
        if (id == 0)
            std::cout << "Number of MPI ranks = " << np << std::endl;

        if (argc < 2) {
            throw std::runtime_error("Error: Must provide path to input file on the command line.");
        }
        else {
            std::string SimulationType;
            // Determine simulation type from the input file:
            // If "C" - constrained solidification - no temperature data file nor remelting needed
            // If "R" - reduced temperature format - temperature data needed, no remelting
            // If "RM" - full time-temperature history - needs remelting
            std::string InputFile = argv[1];
            size_t backslash = InputFile.find_last_of("/");
            std::string FilePath = InputFile.substr(0, backslash);
            std::ifstream InputData;
            InputData.open(InputFile);
            if (id == 0)
                std::cout << "Input file " << InputFile << " opened" << std::endl;
            skipLines(InputData);
            SimulationType = parseInput(InputData, "Problem type");
            InputData.close();
            if (SimulationType == "RM")
                RunProgram_Remelt(id, np, SimulationType, InputFile);
            else
                RunProgram_Reduced(id, np, SimulationType, InputFile);
        }
    }
    // Finalize Kokkos
    Kokkos::finalize();
    // Finalize MPI
    MPI_Finalize();
    return 0;
}
