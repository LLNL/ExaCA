// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#include <mpi.h>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    int return_val = RUN_ALL_TESTS();
    Kokkos::finalize();
    MPI_Finalize();
    return return_val;
}
