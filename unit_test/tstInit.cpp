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

void testInterfacialResponse() {

    // Test that the interfacial response can be read for the new file format
    std::vector<std::string> material_file_names = {"Inconel625.json", "Inconel625_Quadratic.json", "SS316.json"};
    double deltax = 0.5;
    double deltat = 1.0;
    for (auto file_name : material_file_names) {
        std::cout << "Reading " << file_name << std::endl;
        InterfacialResponseFunction irf(0, file_name, deltat, deltax);

        // Check that fitting parameters were correctly initialized and normalized
        // Fitting parameters should've been normalized by deltat / deltax, i.e. twice as large as the numbers in the
        // file
        double ATest, BTest, CTest, DTest, FreezingRangeTest, ExpectedV;
        double LocU = 11.0;
        if (file_name == "Inconel625.json") {
            ATest = -0.00000010302;
            BTest = 0.00010533;
            CTest = 0.0022196;
            DTest = 0;
            FreezingRangeTest = 210;
            ExpectedV = (deltat / deltax) * (ATest * pow(LocU, 3.0) + BTest * pow(LocU, 2.0) + CTest * LocU + DTest);
        }
        else if (file_name == "Inconel625_Quadratic.json") {
            ATest = 0.000072879;
            BTest = 0.004939;
            CTest = -0.047024;
            FreezingRangeTest = 210;
            ExpectedV = (deltat / deltax) * (ATest * pow(LocU, 2.0) + BTest * LocU + CTest);
        }
        else if (file_name == "SS316.json") {
            ATest = 0.000007325;
            BTest = 3.12;
            CTest = 0;
            FreezingRangeTest = 26.5;
            ExpectedV = (deltat / deltax) * (ATest * pow(LocU, BTest) + CTest);
        }
        else {
            throw std::runtime_error("File not set up for testing.");
        }
        // For all IRFs, A and C should be normalized by deltat/deltax (i.e., 2)
        // For the power law IRF (SS316), B is dimensionless and should not be normalized unlike the other IRFs where
        // all coefficients are normalized
        EXPECT_DOUBLE_EQ(irf.A, ATest * 2);
        if (file_name == "SS316.json")
            EXPECT_DOUBLE_EQ(irf.B, BTest);
        else
            EXPECT_DOUBLE_EQ(irf.B, BTest * 2);
        EXPECT_DOUBLE_EQ(irf.C, CTest * 2);
        if (file_name == "Inconel625.json") {
            EXPECT_DOUBLE_EQ(irf.D, DTest * 2);
        }
        EXPECT_DOUBLE_EQ(irf.FreezingRange, FreezingRangeTest);
        double ComputedV = irf.compute(LocU);
        EXPECT_DOUBLE_EQ(ComputedV, ExpectedV);
    }
}
//---------------------------------------------------------------------------//
// activedomainsizecalc
//---------------------------------------------------------------------------//
void testcalcZBound_Low() {

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
        int ZBound_Low_S = calcZBound_Low("S", LayerHeight, layernumber, ZMinLayer, ZMin, deltax);
        EXPECT_EQ(ZBound_Low_S, LayerHeight * layernumber);
        int ZBound_Low_R = calcZBound_Low("R", LayerHeight, layernumber, ZMinLayer, ZMin, deltax);
        EXPECT_EQ(ZBound_Low_R, LayerHeight * layernumber);
    }
}

void testcalcZBound_High() {

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
        int ZBound_Max_S = calcZBound_High("S", SpotRadius, LayerHeight, layernumber, ZMin, deltax, nz, ZMaxLayer);
        EXPECT_EQ(ZBound_Max_S, SpotRadius + LayerHeight * layernumber);
        int ZBound_Max_R = calcZBound_High("R", SpotRadius, LayerHeight, layernumber, ZMin, deltax, nz, ZMaxLayer);
        EXPECT_EQ(ZBound_Max_R, SpotRadius + LayerHeight * layernumber);
        // For simulation type C, should be independent of layernumber
        int ZBound_Max_C = calcZBound_High("C", SpotRadius, LayerHeight, layernumber, ZMin, deltax, nz, ZMaxLayer);
        EXPECT_EQ(ZBound_Max_C, nz - 1);
    }
}

void testcalcnzActive() {

    int id = 0;
    int ZBound_Low = 5;
    int NumberOfLayers = 10;
    for (int layernumber = 0; layernumber < NumberOfLayers; layernumber++) {
        int ZBound_High = 6 + layernumber;
        int nzActive = calcnzActive(ZBound_Low, ZBound_High, id, layernumber);
        EXPECT_EQ(nzActive, 2 + layernumber);
    }
}

void testcalcLocalActiveDomainSize() {

    int nx = 5;
    int MyYSlices = 4;
    int nzActive = 10;
    int LocalActiveDomainSize = calcLocalActiveDomainSize(nx, MyYSlices, nzActive);
    EXPECT_EQ(LocalActiveDomainSize, 10 * 5 * 4);
}

//---------------------------------------------------------------------------//
// temp_init_test
//---------------------------------------------------------------------------//
void testFindXYZBounds(bool TestBinaryInputRead) {

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

    int TempFilesInSeries = 1;
    std::vector<std::string> temp_paths(TempFilesInSeries);
    temp_paths[0] = TestFilename;
    int LayerHeight = 2;
    int NumberOfLayers = 2;
    // Values to be calculated in FindXYZBounds
    int nx, ny, nz;
    double XMin, XMax, YMin, YMax, ZMin, ZMax;
    double *ZMinLayer = new double[NumberOfLayers];
    double *ZMaxLayer = new double[NumberOfLayers];

    FindXYZBounds("R", 0, deltax, nx, ny, nz, temp_paths, XMin, XMax, YMin, YMax, ZMin, ZMax, LayerHeight,
                  NumberOfLayers, TempFilesInSeries, ZMinLayer, ZMaxLayer, 0);

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

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, fileread_test) {
    // test functions for reading and writing data as binary (true) and ASCII (false)
    testReadWrite(true);
    testReadWrite(false);
    testInterfacialResponse();
}
TEST(TEST_CATEGORY, activedomainsizecalc) {
    testcalcZBound_Low();
    testcalcZBound_High();
    testcalcnzActive();
    testcalcLocalActiveDomainSize();
}
TEST(TEST_CATEGORY, temperature_init_test) {
    // reading temperature files to obtain xyz bounds, using binary/non-binary format
    testFindXYZBounds(false);
    testFindXYZBounds(true);
}
} // end namespace Test
