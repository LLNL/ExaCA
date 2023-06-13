// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAfunctions.hpp"
#include "CAinitialize.hpp"
#include "CAparsefiles.hpp"
#include "CAprint.hpp"
#include "CAtypes.hpp"

#include <gtest/gtest.h>

#include "mpi.h"

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
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

void testReadTemperatureData(int NumberOfLayers, bool LayerwiseTempRead, bool TestBinaryInputRead) {

    // Create test data
    double deltax = 1 * pow(10, -6);
    double HT_deltax = 1 * pow(10, -6);
    int HTtoCAratio;
    // Domain size is a 3 by 12 region
    int nx = 3;
    int ny = 12;
    // Write fake OpenFOAM data - only rank 0. Temperature data should be of type double
    // Write two files, one or both of which should be read
    std::string TestTempFileName1 = "TestData1";
    if (TestBinaryInputRead)
        TestTempFileName1 = TestTempFileName1 + ".catemp";
    else
        TestTempFileName1 = TestTempFileName1 + ".txt";
    std::ofstream TestDataFile1;
    if (TestBinaryInputRead)
        TestDataFile1.open(TestTempFileName1, std::ios::out | std::ios::binary);
    else {
        TestDataFile1.open(TestTempFileName1);
        TestDataFile1 << "x, y, z, tm, tl, cr" << std::endl;
    }
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            if (TestBinaryInputRead) {
                WriteData(TestDataFile1, static_cast<double>(i * deltax), TestBinaryInputRead);
                WriteData(TestDataFile1, static_cast<double>(j * deltax), TestBinaryInputRead);
                WriteData(TestDataFile1, static_cast<double>(0.0), TestBinaryInputRead);
                WriteData(TestDataFile1, static_cast<double>(i * j), TestBinaryInputRead);
                WriteData(TestDataFile1, static_cast<double>(i * j + i), TestBinaryInputRead);
                WriteData(TestDataFile1, static_cast<double>(i * j + j), TestBinaryInputRead);
            }
            else
                TestDataFile1 << i * deltax << "," << j * deltax << "," << 0.0 << "," << static_cast<double>(i * j)
                              << "," << static_cast<double>(i * j + i) << "," << static_cast<double>(i * j + j)
                              << std::endl;
        }
    }
    TestDataFile1.close();
    std::string TestTempFileName2 = "TestData2";
    if (TestBinaryInputRead)
        TestTempFileName2 = TestTempFileName2 + ".catemp";
    else
        TestTempFileName2 = TestTempFileName2 + ".txt";
    std::ofstream TestDataFile2;
    if (TestBinaryInputRead)
        TestDataFile2.open(TestTempFileName2, std::ios::out | std::ios::binary);
    else {
        TestDataFile2.open(TestTempFileName2);
        TestDataFile2 << "x, y, z, tm, tl, cr" << std::endl;
    }
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            if (TestBinaryInputRead) {
                WriteData(TestDataFile2, static_cast<double>(i * deltax), TestBinaryInputRead);
                WriteData(TestDataFile2, static_cast<double>(j * deltax), TestBinaryInputRead);
                WriteData(TestDataFile2, static_cast<double>(deltax), TestBinaryInputRead);
                WriteData(TestDataFile2, static_cast<double>(i * j), TestBinaryInputRead);
                WriteData(TestDataFile2, static_cast<double>(i * j + i), TestBinaryInputRead);
                WriteData(TestDataFile2, static_cast<double>(i * j + j), TestBinaryInputRead);
            }
            else
                TestDataFile2 << i * deltax << "," << j * deltax << "," << deltax << "," << static_cast<double>(i * j)
                              << "," << static_cast<double>(i * j + i) << "," << static_cast<double>(i * j + j)
                              << std::endl;
        }
    }
    TestDataFile2.close();

    // Test each 12 by 3 subdomain
    for (int id = 0; id < 4; id++) {
        int MyYSlices = 3;
        int MyYOffset = 3 * id; // each col is separated from the others by 3 cells
        double YMin = 0.0;
        int layernumber = 0;
        std::vector<std::string> temp_paths(2);
        temp_paths[0] = TestTempFileName1;
        temp_paths[1] = TestTempFileName2;
        int TempFilesInSeries = 2;
        int *FirstValue = new int[NumberOfLayers];
        int *LastValue = new int[NumberOfLayers];

        // Read in data to "RawData"
        ViewD_H RawTemperatureData(Kokkos::ViewAllocateWithoutInitializing("RawTemperatureData"), 9);

        ReadTemperatureData(id, deltax, HT_deltax, HTtoCAratio, MyYSlices, MyYOffset, YMin, temp_paths, NumberOfLayers,
                            TempFilesInSeries, FirstValue, LastValue, LayerwiseTempRead, layernumber,
                            RawTemperatureData);

        // Check the results.
        // Does each rank have the right number of temperature data points? Each rank should have six (x,y,z,tm,tl,cr)
        // for each of the 9 cells in the subdomain
        // If both files were read, twice as many temperature data points per file should be present
        int NumberOfTemperatureDataPoints = RawTemperatureData.extent(0);
        int NumTempPointsMultiplier;
        if (LayerwiseTempRead)
            NumTempPointsMultiplier = 1;
        else
            NumTempPointsMultiplier = std::min(NumberOfLayers, TempFilesInSeries);
        EXPECT_EQ(NumberOfTemperatureDataPoints, 54 * NumTempPointsMultiplier);
        // Ratio of HT cell size and CA cell size should be 1
        EXPECT_EQ(HTtoCAratio, 1);
        int NumberOfCellsPerRank = 9;
        // Does each rank have the right temperature data values?
        for (int layercounter = 0; layercounter < NumTempPointsMultiplier; layercounter++) {
            for (int n = 0; n < NumberOfCellsPerRank; n++) {
                double ExpectedValues_ThisDataPoint[6];
                // Location on local grid
                int CARow = n % 3;
                int CACol = n / 3;
                // X Coordinate
                ExpectedValues_ThisDataPoint[0] = CARow * deltax;
                // Y Coordinate
                ExpectedValues_ThisDataPoint[1] = (CACol + 3 * id) * deltax;
                // Z Coordinate
                ExpectedValues_ThisDataPoint[2] = deltax * layercounter;
                int XInt = ExpectedValues_ThisDataPoint[0] / deltax;
                int YInt = ExpectedValues_ThisDataPoint[1] / deltax;
                // Melting time
                ExpectedValues_ThisDataPoint[3] = XInt * YInt;
                // Liquidus time
                ExpectedValues_ThisDataPoint[4] = XInt * YInt + XInt;
                // Cooling rate
                ExpectedValues_ThisDataPoint[5] = XInt * YInt + YInt;
                for (int nn = 0; nn < 6; nn++) {
                    EXPECT_DOUBLE_EQ(ExpectedValues_ThisDataPoint[nn],
                                     RawTemperatureData(NumberOfCellsPerRank * 6 * layercounter + 6 * n + nn));
                }
            }
        }
    }
}

void testgetTempCoords() {

    // cell size and simulation lower bounds
    double deltax = 0.5;
    double XMin = -5;
    double YMin = 0.0;

    // test - each rank initializes data for layer id, of 4 total layers
    for (int layernumber = 0; layernumber < 4; layernumber++) {
        double *ZMinLayer = new double[4];
        for (int n = 0; n < 4; n++) {
            ZMinLayer[n] = 2 * n;
        }
        int LayerHeight = 10;

        // Fill "RawData" with example temperature field data - 2 data point (x,y,z,tm,tl,cr)
        ViewD_H RawTemperatureData(Kokkos::ViewAllocateWithoutInitializing("RawTemperatureData"), 12);

        std::vector<double> RawData = {0.0, 0.0, 3.0, 0.000050, 0.000055, 2000.0,
                                       5.0, 2.0, 2.0, 0.000150, 0.000155, 3500.0};
        for (int n = 0; n < 12; n++)
            RawTemperatureData(n) = RawData[n];

        // Read and check RawData values against expected results from getTempCoordX,Y,Z,TM,TL,CR functions
        // First data point (i = 0 through i = 5 values)
        int i = 0;
        int XInt = getTempCoordX(i, XMin, deltax, RawTemperatureData);
        int YInt = getTempCoordY(i, YMin, deltax, RawTemperatureData);
        int ZInt = getTempCoordZ(i, deltax, RawTemperatureData, LayerHeight, layernumber, ZMinLayer);
        double TMelting = getTempCoordTM(i, RawTemperatureData);
        double TLiquidus = getTempCoordTL(i, RawTemperatureData);
        double CoolingRate = getTempCoordCR(i, RawTemperatureData);
        EXPECT_EQ(XInt, 10);
        EXPECT_EQ(YInt, 0);
        EXPECT_EQ(ZInt, 6 + 6 * layernumber); // different for each rank, since LayerCounter = rank id
        EXPECT_DOUBLE_EQ(TMelting, 0.000050);
        EXPECT_DOUBLE_EQ(TLiquidus, 0.000055);
        EXPECT_DOUBLE_EQ(CoolingRate, 2000.0);
        // Second data point (i = 6 through i = 11 values)
        i = 6;
        XInt = getTempCoordX(i, XMin, deltax, RawTemperatureData);
        YInt = getTempCoordY(i, YMin, deltax, RawTemperatureData);
        ZInt = getTempCoordZ(i, deltax, RawTemperatureData, LayerHeight, layernumber, ZMinLayer);
        TMelting = getTempCoordTM(i, RawTemperatureData);
        TLiquidus = getTempCoordTL(i, RawTemperatureData);
        CoolingRate = getTempCoordCR(i, RawTemperatureData);
        EXPECT_EQ(XInt, 20);
        EXPECT_EQ(YInt, 4);
        EXPECT_EQ(ZInt, 4 + 6 * layernumber); // different for each rank, since LayerCounter = rank id
        EXPECT_DOUBLE_EQ(TMelting, 0.000150);
        EXPECT_DOUBLE_EQ(TLiquidus, 0.000155);
        EXPECT_DOUBLE_EQ(CoolingRate, 3500.0);
    }
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, orientation_init_tests) {
    testOrientationInit_Vectors();
    testOrientationInit_Angles();
}
TEST(TEST_CATEGORY, temperature_init_test) {
    // Multiple permutations of inputs: NumberOfLayers, LayerwiseTempRead, TestBinaryInputRead
    // reading temperature data is performed in the same manner with and without remelting
    std::vector<int> NumberOfLayers_vals = {1, 2, 2, 2};
    std::vector<bool> LayerwiseTempRead_vals = {false, true, false, true};
    std::vector<bool> TestBinaryInputRead_vals = {false, false, false, true};
    int num_vals = TestBinaryInputRead_vals.size();
    for (int test_count = 0; test_count < num_vals; test_count++) {
        testReadTemperatureData(NumberOfLayers_vals[test_count], LayerwiseTempRead_vals[test_count],
                                TestBinaryInputRead_vals[test_count]);
    }
    testgetTempCoords();
}
} // end namespace Test
