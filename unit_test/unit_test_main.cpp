// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#include <mpi.h>

#include <unit_test_cmd_line.hpp>

std::string ExaCAUnitTestingInput::command_line_arg = "";

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);

    // Only for input case with optional commandline input.
    if (argc > 1)
        ExaCAUnitTestingInput::command_line_arg = argv[1];

    // Run.
    int return_val = RUN_ALL_TESTS();

    Kokkos::finalize();
    MPI_Finalize();
    return return_val;
}
