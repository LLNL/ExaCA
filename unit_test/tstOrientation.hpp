// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAconfig.hpp"
#include "CAorientation.hpp"
#include "CAparsefiles.hpp"

#include <gtest/gtest.h>

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
// orientation init tests
//---------------------------------------------------------------------------//
void testOrientationInit_Vectors() {

    using memory_space = TEST_MEMSPACE;

    std::string grain_orientation_file = checkFileInstalled("GrainOrientationVectors.csv", 0);

    // Initialize grain orientations - unit vector form only
    Orientation<memory_space> orientation(grain_orientation_file, false);

    // Check results for first 2 orientations (first 18 values in the file)
    EXPECT_EQ(orientation.n_grain_orientations, 10000);

    std::vector<float> expected_grain_unit_vector = {0.848294,  0.493303,  0.19248,  -0.522525, 0.720911,  0.455253,
                                                     0.0858167, -0.486765, 0.869308, 0.685431,  0.188182,  0.7034,
                                                     -0.468504, 0.85348,   0.228203, -0.557394, -0.485963, 0.673166};
    auto grain_unit_vector_host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orientation.grain_unit_vector);
    for (int n = 0; n < 18; n++) {
        EXPECT_FLOAT_EQ(grain_unit_vector_host(n), expected_grain_unit_vector[n]);
    }
}

void testOrientationInit_Angles() {

    using memory_space = TEST_MEMSPACE;

    std::string grain_orientation_file = checkFileInstalled("GrainOrientationVectors.csv", 0);

    // Initialize grain orientations - unit vector form and data from GrainOrientationEulerAnglesBungeZXZ.csv should be
    // read
    Orientation<memory_space> orientation(grain_orientation_file, true);

    // Check results
    EXPECT_EQ(orientation.n_grain_orientations, 10000);

    // Check first two orientations (first 6 values in the file)
    std::vector<float> expected_euler_angles = {9.99854, 29.62172, 22.91854, 311.08350, 47.68814, 72.02547};
    for (int n = 0; n < 6; n++) {
        EXPECT_FLOAT_EQ(orientation.grain_bunge_euler_host(n), expected_euler_angles[n]);
    }
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, orientation) {
    testOrientationInit_Vectors();
    testOrientationInit_Angles();
}
} // end namespace Test
