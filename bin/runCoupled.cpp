// Copyright 2021-2024 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "ExaCA.hpp"
#include "Finch_Core.hpp"

#include <Kokkos_Core.hpp>

#include "mpi.h"

#include <cctype>
#include <stdexcept>
#include <string>

void runCoupled(int id, int np, Finch::Inputs finch_inputs, Inputs exaca_inputs) {
    // Set up Finch
    using exec_space = Kokkos::DefaultExecutionSpace;
    using memory_space = exec_space::memory_space;

    // Create timers
    Timers timers(id);
    timers.startHeatTransfer();

    // initialize a moving beam
    Finch::MovingBeam beam(finch_inputs.source.scan_path_file);

    // Define boundary condition details.
    std::array<std::string, 6> bc_types = {"adiabatic", "adiabatic", "adiabatic",
                                           "adiabatic", "adiabatic", "adiabatic"};

    // create the global mesh
    Finch::Grid<memory_space> finch_grid(MPI_COMM_WORLD, finch_inputs.space.cell_size,
                                         finch_inputs.space.global_low_corner, finch_inputs.space.global_high_corner,
                                         finch_inputs.space.ranks_per_dim, bc_types,
                                         finch_inputs.space.initial_temperature);

    // Run Finch heat transfer
    Finch::Layer app(finch_inputs, finch_grid);
    auto fd = Finch::createSolver(finch_inputs, finch_grid);
    app.run(exec_space(), finch_inputs, finch_grid, beam, fd);

    // Start of ExaCA init time.
    timers.stopHeatTransfer();
    timers.startInit();

    // Setup local and global grids, decomposing domain (needed to construct temperature)
    Grid exaca_grid(id, np, exaca_inputs.domain, finch_inputs.space.cell_size, finch_inputs.space.global_low_corner,
                    finch_inputs.space.global_high_corner);
    // Temperature fields characterized by data in this structure
    Temperature<memory_space> temperature(id, np, exaca_grid, exaca_inputs.temperature, app.getSolidificationData(),
                                          exaca_inputs.print.store_solidification_start);

    // Now run ExaCA
    runExaCA(id, np, exaca_inputs, timers, exaca_grid, temperature);
}

int main(int argc, char *argv[]) {
    // Initialize Kokkos & MPI
    int id, np;
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);
    {
        // Get number of processes
        MPI_Comm_size(MPI_COMM_WORLD, &np);
        // Get individual process ID
        MPI_Comm_rank(MPI_COMM_WORLD, &id);

        if (id == 0) {
            std::cout << "ExaCA version: " << version() << " \nExaCA commit:  " << gitCommitHash()
                      << "\nKokkos version: " << kokkosVersion() << std::endl;
            Kokkos::DefaultExecutionSpace().print_configuration(std::cout);
            std::cout << "Number of MPI ranks = " << np << std::endl;
        }
        if (argc < 3) {
            throw std::runtime_error("Error: Must provide path to Finch and ExaCA input files on the command line.");
        }
        else {
            std::string finch_input_file = argv[1];
            std::string exaca_input_file = argv[2];

            // Setup Finch simulation
            Finch::Inputs finch_inputs(MPI_COMM_WORLD, finch_input_file);

            // Setup ExaCA simulation details
            Inputs inputs(id, exaca_input_file);

            runCoupled(id, np, finch_inputs, inputs);
        }
    }

    // Finalize Kokkos & MPI
    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}
