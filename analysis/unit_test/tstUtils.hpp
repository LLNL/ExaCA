// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "GAutils.hpp"

#include <gtest/gtest.h>

#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
// GAutils tests without Kokkos
//---------------------------------------------------------------------------//
// Test a function to take the unique grain IDs from a vector of grain IDs with repeats
void testgetUniqueGrains() {

    std::vector<int> GrainIDVector = {1, 1, 1, 2, 2, 3, 5, 5, 2, 3, 10};
    int NumberOfGrains;
    std::vector<int> UniqueGrainIDVector = getUniqueGrains(GrainIDVector, NumberOfGrains);
    EXPECT_EQ(NumberOfGrains, 5);
    std::vector<int> UniqueGrainIDVector_Expected = {1, 2, 3, 5, 10};
    for (int n = 0; n < NumberOfGrains; n++)
        EXPECT_EQ(UniqueGrainIDVector[n], UniqueGrainIDVector_Expected[n]);
}

// Test a function that obtains the volume total for each of the unique grain ID values appears in a larger list of
// grain ID values with repeats
void testgetGrainSizes() {

    std::vector<int> GrainIDVector = {1, 1, 1, 2, 2, 3, 5, 5, 2, 3, 10};
    std::vector<int> UniqueGrainIDVector = {1, 2, 3, 5, 10};
    int NumberOfGrains = 5;
    double deltax = 2 * pow(10, -6);
    std::vector<float> GrainSizeVector =
        getGrainSizes(GrainIDVector, UniqueGrainIDVector, NumberOfGrains, deltax, "volume");
    std::vector<float> GrainSizeVector_Expected = {3, 3, 2, 2, 1};
    for (int n = 0; n < NumberOfGrains; n++) {
        // Convert expected sizes to volumes
        float GrainSizeExpected = GrainSizeVector_Expected[n] * convertToMicrons(deltax, "volume");
        EXPECT_FLOAT_EQ(GrainSizeVector[n], GrainSizeExpected);
    }
}

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, grain_util_tests) {
    testgetUniqueGrains();
    testgetGrainSizes();
}
} // end namespace Test
