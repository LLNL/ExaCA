// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
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

void runCoupled(int id, int np, Finch::Inputs inputs, Inputs exaca_inputs) {
    // Set up Finch
    using exec_space = Kokkos::DefaultExecutionSpace;
    using memory_space = exec_space::memory_space;

    // Create timers
    Timers timers(id);
    timers.startHeatTransfer();

    // initialize a moving beam
    Finch::MovingBeam beam(inputs.source.scan_path_file);

    // Define boundary condition details.
    std::array<std::string, 6> bc_types = {"adiabatic", "adiabatic", "adiabatic",
                                           "adiabatic", "adiabatic", "adiabatic"};

    // create the global mesh
    Finch::Grid<memory_space> grid(MPI_COMM_WORLD, inputs.space.cell_size, inputs.space.global_low_corner,
                                   inputs.space.global_high_corner, inputs.space.ranks_per_dim, bc_types,
                                   inputs.space.initial_temperature);

    // Run Finch heat transfer
    Finch::SingleLayer app(inputs, grid);
    app.run(exec_space(), inputs, grid, beam, fd);

    // Start of ExaCA init time.
    timers.endHeatTransfer();
    timers.startInit();

    // Setup local and global grids, decomposing domain (needed to construct temperature)
    Grid grid(inputs.simulation_type, id, np, inputs.domain.number_of_layers, inputs.domain, inputs.temperature);
    // Temperature fields characterized by data in this structure
    Temperature<memory_space> temperature(grid, app.getSolidificationData(), inputs.print.store_solidification_start);

    // Now run ExaCA
    runExaCA(id, np, exaca_inputs, timers, grid, temperature);
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
        if (argc < 2) {
            throw std::runtime_error("Error: Must provide path to input file on the command line.");
        }
        else {
            std::string input_file = argv[1];

            // Setup Finch simulation
            Finch::Inputs finch_inputs(MPI_COMM_WORLD, argc, argv);

            // Setup ExaCA simulation details
            Inputs inputs(id, input_file);

            runCoupled(id, np, finch_inputs, inputs);
        }
    }

    // Finalize Kokkos & MPI
    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}
