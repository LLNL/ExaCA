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

//---------------------------------------------------------------------------//
// activedomainsizecalc
//---------------------------------------------------------------------------//
void calcz_layer_bottom() {

    int LayerHeight = 10;
    int NumberOfLayers = 10;
    double deltax = 1 * pow(10, -6);
    double ZMin = -0.5 * pow(10, -6);
    double *ZMinLayer = new double[NumberOfLayers];
    for (int layernumber = 0; layernumber < NumberOfLayers; layernumber++) {
        // Set ZMinLayer for each layer to be offset by LayerHeight cells from the previous one (lets solution for
        // both problem types by the same)
        ZMinLayer[layernumber] = ZMin + layernumber * LayerHeight * deltax;
        // Call function for each layernumber, and for simulation types "S" and "R"
        int z_layer_bottom_S = calc_z_layer_bottom("S", LayerHeight, layernumber, ZMinLayer, ZMin, deltax);
        EXPECT_EQ(z_layer_bottom_S, LayerHeight * layernumber);
        int z_layer_bottom_R = calc_z_layer_bottom("R", LayerHeight, layernumber, ZMinLayer, ZMin, deltax);
        EXPECT_EQ(z_layer_bottom_R, LayerHeight * layernumber);
    }
}

void calcz_layer_top() {

    // A separate function is now used for ZBound_High calculation
    int SpotRadius = 100;
    int LayerHeight = 10;
    double ZMin = 0.5 * pow(10, -6);
    double deltax = 1.0 * pow(10, -6);
    int nz = 101;
    int NumberOfLayers = 10;
    double *ZMaxLayer = new double[NumberOfLayers];
    for (int layernumber = 0; layernumber < NumberOfLayers; layernumber++) {
        // Set ZMaxLayer for each layer to be offset by LayerHeight cells from the previous one, with layer 0 having a
        // ZMax value of ZMin + SpotRadius (lets solution for both problem types be the same)
        ZMaxLayer[layernumber] = ZMin + SpotRadius * deltax + layernumber * LayerHeight * deltax;
        // Call function for each layernumber, and for simulation types "S" and "R"
        int z_layer_top_S = calc_z_layer_top("S", SpotRadius, LayerHeight, layernumber, ZMin, deltax, nz, ZMaxLayer);
        EXPECT_EQ(z_layer_top_S, SpotRadius + LayerHeight * layernumber);
        int z_layer_top_R = calc_z_layer_top("R", SpotRadius, LayerHeight, layernumber, ZMin, deltax, nz, ZMaxLayer);
        EXPECT_EQ(z_layer_top_R, SpotRadius + LayerHeight * layernumber);
        // For simulation type C, should be independent of layernumber
        int z_layer_top_C = calc_z_layer_top("C", SpotRadius, LayerHeight, layernumber, ZMin, deltax, nz, ZMaxLayer);
        EXPECT_EQ(z_layer_top_C, nz - 1);
    }
}

void testcalc_nz_layer() {

    int id = 0;
    int z_layer_bottom = 5;
    int NumberOfLayers = 10;
    for (int layernumber = 0; layernumber < NumberOfLayers; layernumber++) {
        int z_layer_top = 6 + layernumber;
        int nz_layer = calc_nz_layer(z_layer_bottom, z_layer_top, id, layernumber);
        EXPECT_EQ(nz_layer, 2 + layernumber);
    }
}

void testcalcLayerDomainSize() {

    int nx = 5;
    int ny_local = 4;
    int nz_layer = 10;
    int DomainSize = calcLayerDomainSize(nx, ny_local, nz_layer);
    EXPECT_EQ(DomainSize, 10 * 5 * 4);
}

//---------------------------------------------------------------------------//
// bounds_init_test
//---------------------------------------------------------------------------//
void testFindXYZBounds(bool TestBinaryInputRead) {

    using memory_space = TEST_MEMSPACE;

    // Write fake OpenFOAM data - temperature data should be of type double
    double deltax = 1 * pow(10, -6);
    std::string TestFilename = "TestData";
    if (TestBinaryInputRead)
        TestFilename = TestFilename + ".catemp";
    else
        TestFilename = TestFilename + ".txt";
    std::ofstream TestData;
    if (TestBinaryInputRead)
        TestData.open(TestFilename, std::ios::out | std::ios::binary);
    else {
        TestData.open(TestFilename);
        TestData << "x, y, z, tm, tl, cr" << std::endl;
    }
    // only x,y,z data should be read, tm, tl, cr should not affect result
    for (int k = 0; k < 3; k++) {
        for (int j = 0; j < 3; j++) {
            for (int i = 0; i < 4; i++) {
                if (TestBinaryInputRead) {
                    WriteData(TestData, static_cast<double>(i * deltax), TestBinaryInputRead);
                    WriteData(TestData, static_cast<double>(j * deltax), TestBinaryInputRead);
                    WriteData(TestData, static_cast<double>(k * deltax), TestBinaryInputRead);
                    WriteData(TestData, static_cast<double>(-1.0), TestBinaryInputRead);
                    WriteData(TestData, static_cast<double>(-1.0), TestBinaryInputRead);
                    WriteData(TestData, static_cast<double>(-1.0), TestBinaryInputRead);
                }
                else
                    TestData << i * deltax << "," << j * deltax << "," << k * deltax << "," << static_cast<double>(-1.0)
                             << "," << static_cast<double>(-1.0) << "," << static_cast<double>(-1.0) << std::endl;
            }
        }
    }
    TestData.close();

    // Empty inputs struct with default values - manually set non-default substrateInputs values
    Inputs<memory_space> inputs;
    inputs.SimulationType = "R";
    inputs.temperatureInputs.TempFilesInSeries = 1;
    inputs.temperatureInputs.temp_paths.push_back(TestFilename);
    int LayerHeight = 2;
    int NumberOfLayers = 2;
    // Values to be calculated in FindXYZBounds
    int nx, ny, nz;
    double XMin, XMax, YMin, YMax, ZMin, ZMax;
    double *ZMinLayer = new double[NumberOfLayers];
    double *ZMaxLayer = new double[NumberOfLayers];

    FindXYZBounds(0, deltax, nx, ny, nz, XMin, XMax, YMin, YMax, ZMin, ZMax, ZMinLayer, ZMaxLayer, NumberOfLayers,
                  LayerHeight, inputs);

    EXPECT_DOUBLE_EQ(XMin, 0.0);
    EXPECT_DOUBLE_EQ(YMin, 0.0);
    EXPECT_DOUBLE_EQ(ZMin, 0.0);
    EXPECT_DOUBLE_EQ(XMax, 3 * deltax);
    EXPECT_DOUBLE_EQ(YMax, 2 * deltax);
    // ZMax is equal to the largest Z coordinate in the file, offset by LayerHeight cells in the build direction due
    // to the second layer
    EXPECT_DOUBLE_EQ(ZMax, 4 * deltax);
    // Bounds for each individual layer - 2nd layer offset by LayerHeight cells from the first
    EXPECT_DOUBLE_EQ(ZMinLayer[0], 0.0);
    EXPECT_DOUBLE_EQ(ZMaxLayer[0], 2 * deltax);
    EXPECT_DOUBLE_EQ(ZMinLayer[1], 2 * deltax);
    EXPECT_DOUBLE_EQ(ZMaxLayer[1], 4 * deltax);
    // Size of overall domain
    EXPECT_EQ(nx, 4);
    EXPECT_EQ(ny, 3);
    EXPECT_EQ(nz, 5);
}

void testOrientationInit_Vectors() {

    int id;
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    int ValsPerLine = 9;
    int NGrainOrientations = 0;
    std::string GrainOrientationFile = checkFileInstalled("GrainOrientationVectors.csv", id);

    // View for storing orientation data
    ViewF GrainOrientationData(Kokkos::ViewAllocateWithoutInitializing("GrainOrientationData"), 0);

    // Call OrientationInit - without optional final argument
    OrientationInit(id, NGrainOrientations, GrainOrientationData, GrainOrientationFile);

    // Copy orientation data back to the host
    ViewF_H GrainOrientationData_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainOrientationData);

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

    int id;
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    int ValsPerLine = 3;
    int NGrainOrientations = 0;
    std::string GrainOrientationFile = checkFileInstalled("GrainOrientationEulerAnglesBungeZXZ.csv", id);

    // View for storing orientation data
    ViewF GrainOrientationData(Kokkos::ViewAllocateWithoutInitializing("GrainOrientationData"), 0);

    // Call OrientationInit - with optional final argument
    OrientationInit(id, NGrainOrientations, GrainOrientationData, GrainOrientationFile, ValsPerLine);

    // Copy orientation data back to the host
    ViewF_H GrainOrientationData_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainOrientationData);

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
TEST(TEST_CATEGORY, activedomainsizecalc) {
    calcz_layer_bottom();
    calcz_layer_top();
    testcalc_nz_layer();
    testcalcLayerDomainSize();
}
TEST(TEST_CATEGORY, bounds_init_test) {
    // reading temperature files to obtain xyz bounds, using binary/non-binary format
    testFindXYZBounds(false);
    testFindXYZBounds(true);
}
TEST(TEST_CATEGORY, orientation_init_tests) {
    testOrientationInit_Vectors();
    testOrientationInit_Angles();
}
} // end namespace Test
