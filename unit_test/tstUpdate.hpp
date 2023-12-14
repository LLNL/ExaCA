// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>
#include <nlohmann/json.hpp>

#include "CAcelldata.hpp"
#include "CAinitialize.hpp"
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

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    std::string InputFile = "Inp_SmallDirSolidification.json";

    // Run SmallDirS problem and check volume fraction of nucleated grains with 1% tolerance of expected value (to
    // account for the non-deterministic nature of the cell capture)
    RunProgram_Reduced(id, np, InputFile);

    // MPI barrier to ensure that log file has been written
    MPI_Barrier(MPI_COMM_WORLD);
    std::string LogFile = "TestProblemSmallDirS.json";
    std::ifstream LogDataStream(LogFile);
    nlohmann::json logdata = nlohmann::json::parse(LogDataStream);
    float VolFractionNucleated = logdata["Nucleation"]["VolFractionNucleated"];
    EXPECT_NEAR(VolFractionNucleated, 0.1784, 0.0100);
}

void testSmallEquiaxedGrain() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    std::string InputFile = "Inp_SmallEquiaxedGrain.json";

    // Run Small equiaxed grain problem and check time step at which the grain reaches the domain edge
    RunProgram_Reduced(id, np, InputFile);

    // MPI barrier to ensure that log file has been written
    MPI_Barrier(MPI_COMM_WORLD);
    std::string LogFile = "TestProblemSmallEquiaxedGrain.json";
    std::ifstream LogDataStream(LogFile);
    nlohmann::json logdata = nlohmann::json::parse(LogDataStream);
    int TimeStepOfOutput = logdata["TimeStepOfOutput"];
    // FIXME: Output time step is usually 4820, but may be 4821 - need to investigate this possible race condition
    EXPECT_NEAR(TimeStepOfOutput, 4820, 1);
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, full_simulations) {
    testSmallDirS();
    testSmallEquiaxedGrain();
}
} // end namespace Test
