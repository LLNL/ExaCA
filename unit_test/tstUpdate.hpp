// Copyright Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>
#include <nlohmann/json.hpp>

#include "CAcelldata.hpp"
#include "CAinputs.hpp"
#include "CAinterface.hpp"
#include "CAparsefiles.hpp"
#include "CAtypes.hpp"
#include "ExaCA.hpp"

#include <gtest/gtest.h>

#include "mpi.h"

#include <fstream>
#include <string>
#include <vector>

namespace Test {

//---------------------------------------------------------------------------//
// full_simulations
//---------------------------------------------------------------------------//
void testSmallDirS() {

    using memory_space = TEST_MEMSPACE;
    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    std::string input_file = "Inp_SmallDirSolidification.json";
    Inputs inputs(id, input_file);
    Timers timers(id);

    // Setup local and global grids, decomposing domain (needed to construct temperature)
    Grid grid(inputs.simulation_type, id, np, inputs.domain.number_of_layers, inputs.domain, inputs.substrate,
              inputs.temperature);
    // Temperature fields characterized by data in this structure
    Temperature<memory_space> temperature(grid, inputs.temperature, inputs.print);

    // Run SmallDirS problem and check volume fraction of nucleated grains with 1% tolerance of expected value (to
    // account for the non-deterministic nature of the cell capture)
    runExaCA(id, np, inputs, timers, grid, temperature);

    // MPI barrier to ensure that log file has been written
    MPI_Barrier(MPI_COMM_WORLD);
    std::string log_file = "TestProblemSmallDirS.json";
    std::ifstream log_data_stream(log_file);
    nlohmann::json log_data = nlohmann::json::parse(log_data_stream);
    float vol_fraction_nucleated = log_data["Nucleation"]["VolFractionNucleated"];
    EXPECT_NEAR(vol_fraction_nucleated, 0.1882, 0.0100);
}

void testSmallEquiaxedGrain() {

    using memory_space = TEST_MEMSPACE;
    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    std::string input_file = "Inp_SmallEquiaxedGrain.json";
    Inputs inputs(id, input_file);
    Timers timers(id);

    // Setup local and global grids, decomposing domain (needed to construct temperature)
    Grid grid(inputs.simulation_type, id, np, inputs.domain.number_of_layers, inputs.domain, inputs.substrate,
              inputs.temperature);
    // Temperature fields characterized by data in this structure
    Temperature<memory_space> temperature(grid, inputs.temperature, inputs.print);

    // Run Small equiaxed grain problem and check time step at which the grain reaches the domain edge
    runExaCA(id, np, inputs, timers, grid, temperature);

    // MPI barrier to ensure that log file has been written
    MPI_Barrier(MPI_COMM_WORLD);
    std::string log_file = "TestProblemSmallEquiaxedGrain.json";
    std::ifstream log_data_stream(log_file);
    nlohmann::json log_data = nlohmann::json::parse(log_data_stream);
    int time_step_of_output = log_data["TimeStepOfOutput"];
    // FIXME: Output time step is usually 4820, but may be 1 time step larger - need to investigate this possible race
    // condition. Finishes 1 time step earlier now after update to active cell init for no remelting case
    EXPECT_NEAR(time_step_of_output, 4820, 1);
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, full_simulations) {
    testSmallDirS();
    testSmallEquiaxedGrain();
}
} // end namespace Test
