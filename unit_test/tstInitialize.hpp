// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAconfig.hpp"
#include "CAinitialize.hpp"
#include "CAinterfacialresponse.hpp"
#include "CAparsefiles.hpp"
#include "CAprint.hpp"

#include <gtest/gtest.h>

#include "mpi.h"

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
// file_read_tests
//---------------------------------------------------------------------------//
void testReadWrite(bool PrintReadBinary) {

    // Make lists of some int and float data
    int IntData[5] = {-2, 0, 2, 4, 6};
    float FloatData[5] = {-1.0, 0.0, 1.0, 2.0, 3.0};

    // Write data as binary to be used as input
    std::ofstream TestIntData;
    std::ofstream TestFloatData;
    if (PrintReadBinary) {
        TestIntData.open("TestIntData.txt", std::ios::out | std::ios::binary);
        TestFloatData.open("TestFloatData.txt", std::ios::out | std::ios::binary);
    }
    else {
        TestIntData.open("TestIntData.txt");
        TestFloatData.open("TestFloatData.txt");
    }
    for (int n = 0; n < 5; n++) {
        // Write to files
        WriteData(TestIntData, IntData[n], PrintReadBinary, true);
        WriteData(TestFloatData, FloatData[n], PrintReadBinary, true);
    }
    TestIntData.close();
    TestFloatData.close();

    // Read data and convert back to ints and floats, compare to original values
    std::ifstream TestIntDataRead;
    TestIntDataRead.open("TestIntData.txt");
    std::ifstream TestFloatDataRead;
    TestFloatDataRead.open("TestFloatData.txt");
    // For reading ASCII data, obtain the lines from the files first, then parse the string stream at the spaces
    if (PrintReadBinary) {
        for (int n = 0; n < 5; n++) {
            int IntToCompare = ReadBinaryData<int>(TestIntDataRead, true);
            float FloatToCompare = ReadBinaryData<float>(TestFloatDataRead, true);
            // Compare to expected values
            EXPECT_EQ(IntToCompare, IntData[n]);
            EXPECT_FLOAT_EQ(FloatToCompare, FloatData[n]);
        }
    }
    else {
        std::string intline, floatline;
        getline(TestIntDataRead, intline);
        getline(TestFloatDataRead, floatline);
        std::istringstream intss(intline);
        std::istringstream floatss(floatline);
        for (int n = 0; n < 5; n++) {
            // Get values from string stream
            int IntToCompare = ParseASCIIData<float>(intss);
            float FloatToCompare = ParseASCIIData<float>(floatss);
            // Compare to expected values
            EXPECT_EQ(IntToCompare, IntData[n]);
            EXPECT_FLOAT_EQ(FloatToCompare, FloatData[n]);
        }
    }
}

// TODO: Relocate to future tstOrientation.hpp
void testOrientationInit_Vectors() {

    using memory_space = TEST_MEMSPACE;

    int id;
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    int ValsPerLine = 9;
    int NGrainOrientations = 0;
    std::string GrainOrientationFile = checkFileInstalled("GrainOrientationVectors.csv", id);

    // View for storing orientation data
    Kokkos::View<float *, memory_space> GrainOrientationData(
        Kokkos::ViewAllocateWithoutInitializing("GrainOrientationData"), 0);

    // Call OrientationInit - without optional final argument
    OrientationInit(id, NGrainOrientations, GrainOrientationData, GrainOrientationFile);

    // Copy orientation data back to the host
    auto GrainOrientationData_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainOrientationData);

    // Check results
    EXPECT_EQ(NGrainOrientations, 10000);

    std::vector<float> ExpectedGrainOrientations = {0.848294,  0.493303,  0.19248,  -0.522525, 0.720911,  0.455253,
                                                    0.0858167, -0.486765, 0.869308, 0.685431,  0.188182,  0.7034,
                                                    -0.468504, 0.85348,   0.228203, -0.557394, -0.485963, 0.673166};
    for (int n = 0; n < 2 * ValsPerLine; n++) {
        EXPECT_FLOAT_EQ(GrainOrientationData_Host(n), ExpectedGrainOrientations[n]);
    }
}

void testOrientationInit_Angles() {

    using memory_space = TEST_MEMSPACE;

    int id;
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    int ValsPerLine = 3;
    int NGrainOrientations = 0;
    std::string GrainOrientationFile = checkFileInstalled("GrainOrientationEulerAnglesBungeZXZ.csv", id);

    // View for storing orientation data
    Kokkos::View<float *, memory_space> GrainOrientationData(
        Kokkos::ViewAllocateWithoutInitializing("GrainOrientationData"), 0);

    // Call OrientationInit - with optional final argument
    OrientationInit(id, NGrainOrientations, GrainOrientationData, GrainOrientationFile, ValsPerLine);

    // Copy orientation data back to the host
    auto GrainOrientationData_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainOrientationData);

    // Check results
    EXPECT_EQ(NGrainOrientations, 10000);

    // Check first two orientations
    std::vector<float> ExpectedGrainOrientations = {9.99854, 29.62172, 22.91854, 311.08350, 47.68814, 72.02547};
    for (int n = 0; n < 2 * ValsPerLine; n++) {
        EXPECT_FLOAT_EQ(GrainOrientationData_Host(n), ExpectedGrainOrientations[n]);
    }
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, fileread_test) {
    // test functions for reading and writing data as binary (true) and ASCII (false)
    testReadWrite(true);
    testReadWrite(false);
}
TEST(TEST_CATEGORY, orientation_init_tests) {
    testOrientationInit_Vectors();
    testOrientationInit_Angles();
}
} // end namespace Test
