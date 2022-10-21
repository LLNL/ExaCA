// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAinitialize.hpp"
#include "CAinterfacialresponse.hpp"

#include <gtest/gtest.h>

#include "mpi.h"

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
// file_read_tests
//---------------------------------------------------------------------------//
void testInputReadFromFile() {

    int id = 0;
    // Three input files - one of each type
    // Inp_DirSolidification.txt and Inp_SpotMelt.txt were installed from the examples directory
    // Since no temperature files exist in the repo, and there is no ability to write temperature files to a different
    // directory ( would need examples/Temperatures) using the C++11 standard, dummy input files are written and parsed
    // to test an example problem that uses temperature data from a file.
    std::vector<std::string> InputFilenames = {"Inp_DirSolidification.txt", "Inp_SpotMelt.txt",
                                               "Inp_TemperatureTest_Old.txt", "Inp_TemperatureTest.txt"};
    std::vector<std::string> TemperatureFNames = {"1DummyTemperature.txt", "2DummyTemperature.txt"};

    // Write dummy input files for using read temperature data (InputFilenames[2] and [3])
    for (int n = 2; n < 4; n++) {
        std::ofstream TestDataFile;
        TestDataFile.open(InputFilenames[n]);
        // Write required inputs to all files
        TestDataFile << "Test problem data set" << std::endl;
        TestDataFile << "*****" << std::endl;
        TestDataFile << "Problem type: R" << std::endl;
        TestDataFile << "Decomposition strategy: 1" << std::endl;
        TestDataFile << "Material: Inconel625" << std::endl;
        TestDataFile << "Cell size: 1" << std::endl;
        TestDataFile << "Heterogeneous nucleation density: 10" << std::endl;
        TestDataFile << "Mean nucleation undercooling: 5" << std::endl;
        TestDataFile << "Standard deviation of nucleation undercooling: 0.5" << std::endl;
        TestDataFile << "Path to output: ExaCA" << std::endl;
        TestDataFile << "Output file base name: Test" << std::endl;
        TestDataFile << "File of grain orientations: GrainOrientationVectors.csv" << std::endl;
        TestDataFile << "Print file of grain misorientation values: Y" << std::endl;
        TestDataFile << "Print file of final undercooling values: Y" << std::endl;
        TestDataFile << "Print default RVE output: Y" << std::endl;
        TestDataFile << "Print file of all ExaCA data: N" << std::endl;
        TestDataFile << "Time step: 1.5" << std::endl;
        TestDataFile << "Density of powder surface sites active: 1000" << std::endl;
        if (n == 2) {
            // Deprecated temperature input lines
            TestDataFile << "Path to temperature file(s): ./" << std::endl;
            TestDataFile << "Temperature filename(s): DummyTemperature.txt" << std::endl;
            TestDataFile << "Number of temperature files: 2" << std::endl;
            TestDataFile << "Number of layers: 2" << std::endl;
            TestDataFile << "Offset between layers: 1" << std::endl;
            TestDataFile << "Heat transport data mesh size: 12" << std::endl;
        }
        else {
            // New temperature input lines
            TestDataFile << "Path to and name of temperature field assembly instructions: TInstructions.txt"
                         << std::endl;
            // Write test temperature instructions file - don't give HT_deltax, let value default to deltax
            std::ofstream TestTField;
            TestTField.open("TInstructions.txt");
            TestTField << "Number of layers: 2" << std::endl;
            TestTField << "Offset between layers: 1" << std::endl;
            TestTField << "*****" << std::endl;
            TestTField << ".//" << TemperatureFNames[0] << std::endl;
            TestTField << ".//" << TemperatureFNames[1] << std::endl;
            TestTField.close();
        }
        TestDataFile << "Extra set of wall cells in lateral domain directions: N" << std::endl;
        TestDataFile << "Random seed for grains and nuclei generation: 2.0" << std::endl;
        TestDataFile << "Substrate filename: DummySubstrate.txt" << std::endl;
        TestDataFile.close();
    }

    // Create test temperature files "1DummyTemperature.txt", "2DummyTemperature.txt"
    std::ofstream TestTemp1, TestTemp2;
    TestTemp1.open(TemperatureFNames[0]);
    TestTemp2.open(TemperatureFNames[1]);
    TestTemp1 << "X" << std::endl;
    TestTemp2 << "X" << std::endl;
    TestTemp1.close();
    TestTemp2.close();

    // Create test substrate file "DummySubstrate.txt"
    std::ofstream TestSub;
    TestSub.open("DummySubstrate.txt");
    TestSub << "X" << std::endl;
    TestSub.close();

    // Read and parse each input file
    for (auto FileName : InputFilenames) {
        int TempFilesInSeries, NumberOfLayers, LayerHeight, nx, ny, nz, PrintDebug, NSpotsX, NSpotsY, SpotOffset,
            SpotRadius, TimeSeriesInc, RVESize;
        float SubstrateGrainSpacing;
        double deltax, NMax, dTN, dTsigma, HT_deltax, deltat, G, R, FractSurfaceSitesActive, RNGSeed, PowderDensity;
        bool RemeltingYN, PrintMisorientation, PrintFinalUndercoolingVals, PrintFullOutput, PrintTimeSeries,
            UseSubstrateFile, PrintIdleTimeSeriesFrames, PrintDefaultRVE = false, BaseplateThroughPowder;
        std::string SimulationType, OutputFile, GrainOrientationFile, temppath, tempfile, SubstrateFileName,
            PathToOutput, MaterialFileName;
        std::vector<std::string> temp_paths;
        InputReadFromFile(id, FileName, SimulationType, deltax, NMax, dTN, dTsigma, OutputFile, GrainOrientationFile,
                          TempFilesInSeries, temp_paths, HT_deltax, RemeltingYN, deltat, NumberOfLayers, LayerHeight,
                          MaterialFileName, SubstrateFileName, SubstrateGrainSpacing, UseSubstrateFile, G, R, nx, ny,
                          nz, FractSurfaceSitesActive, PathToOutput, PrintDebug, PrintMisorientation,
                          PrintFinalUndercoolingVals, PrintFullOutput, NSpotsX, NSpotsY, SpotOffset, SpotRadius,
                          PrintTimeSeries, TimeSeriesInc, PrintIdleTimeSeriesFrames, PrintDefaultRVE, RNGSeed,
                          BaseplateThroughPowder, PowderDensity, RVESize);
        InterfacialResponseFunction irf(MaterialFileName, deltat, deltax);

        // Check the results
        // The existence of the specified orientation, substrate, and temperature filenames was already checked within
        // InputReadFromFile
        // These should be the same for all 3 test problems
        EXPECT_DOUBLE_EQ(deltax, 1.0 * pow(10, -6));
        EXPECT_DOUBLE_EQ(NMax, 1.0 * pow(10, 13));
        EXPECT_DOUBLE_EQ(dTN, 5.0);
        EXPECT_DOUBLE_EQ(dTsigma, 0.5);
        EXPECT_EQ(PrintDebug, 0);
        EXPECT_DOUBLE_EQ(irf.A, -0.00000010302 * deltat / deltax);
        EXPECT_DOUBLE_EQ(irf.B, 0.00010533 * deltat / deltax);
        EXPECT_DOUBLE_EQ(irf.C, 0.0022196 * deltat / deltax);
        EXPECT_DOUBLE_EQ(irf.D, 0);
        EXPECT_DOUBLE_EQ(irf.FreezingRange, 210);

        // These are different for all 3 test problems
        if (FileName == "Inp_DirSolidification.txt") {
            EXPECT_TRUE(PrintTimeSeries);
            EXPECT_EQ(TimeSeriesInc, 5250);
            EXPECT_FALSE(PrintIdleTimeSeriesFrames);
            EXPECT_DOUBLE_EQ(G, 500000.0);
            EXPECT_DOUBLE_EQ(R, 300000.0);
            EXPECT_DOUBLE_EQ(deltat, pow(10, -6) / 15.0); // based on N, deltax, G, and R from input file
            EXPECT_EQ(nx, 200);
            EXPECT_EQ(ny, 200);
            EXPECT_EQ(nz, 200);
            EXPECT_DOUBLE_EQ(FractSurfaceSitesActive, 0.08);
            EXPECT_TRUE(OutputFile == "TestProblemDirS");
            EXPECT_TRUE(PrintMisorientation);
            EXPECT_FALSE(PrintFinalUndercoolingVals);
            EXPECT_TRUE(PrintFullOutput);
            EXPECT_DOUBLE_EQ(RNGSeed, 0.0);
        }
        else if (FileName == "Inp_SpotMelt.txt") {
            EXPECT_TRUE(PrintTimeSeries);
            EXPECT_EQ(TimeSeriesInc, 37500);
            EXPECT_TRUE(PrintIdleTimeSeriesFrames);
            EXPECT_DOUBLE_EQ(G, 500000.0);
            EXPECT_DOUBLE_EQ(R, 300000.0);
            EXPECT_DOUBLE_EQ(deltat, pow(10, -6) / 15.0); // based on N, deltax, G, and R from input file
            EXPECT_EQ(NSpotsX, 3);
            EXPECT_EQ(NSpotsY, 2);
            EXPECT_EQ(SpotOffset, 100);
            EXPECT_EQ(SpotRadius, 75);
            EXPECT_EQ(NumberOfLayers, 2);
            EXPECT_EQ(LayerHeight, 20);
            EXPECT_FALSE(UseSubstrateFile);
            EXPECT_FALSE(BaseplateThroughPowder);
            EXPECT_FLOAT_EQ(SubstrateGrainSpacing, 25.0);
            EXPECT_TRUE(OutputFile == "TestProblemSpot");
            EXPECT_FALSE(PrintFinalUndercoolingVals);
            EXPECT_TRUE(PrintFullOutput);
            EXPECT_DOUBLE_EQ(RNGSeed, 0.0);
        }
        else if ((FileName == "Inp_TemperatureTest_Old.txt") || (FileName == "Inp_TemperatureTest_Old.txt")) {
            EXPECT_DOUBLE_EQ(deltat, 1.5 * pow(10, -6));
            EXPECT_EQ(TempFilesInSeries, 2);
            EXPECT_EQ(NumberOfLayers, 2);
            EXPECT_EQ(LayerHeight, 1);
            EXPECT_TRUE(UseSubstrateFile);
            EXPECT_FALSE(BaseplateThroughPowder);
            EXPECT_DOUBLE_EQ(PowderDensity, 0.001);
            if (FileName == "Inp_TemperatureTest_Old.txt")
                EXPECT_DOUBLE_EQ(HT_deltax, 12.0 * pow(10, -6));
            else
                EXPECT_DOUBLE_EQ(HT_deltax, deltax);
            EXPECT_TRUE(OutputFile == "Test");
            EXPECT_TRUE(temp_paths[0] == ".//1DummyTemperature.txt");
            EXPECT_TRUE(temp_paths[1] == ".//2DummyTemperature.txt");
            EXPECT_TRUE(PrintMisorientation);
            EXPECT_TRUE(PrintFinalUndercoolingVals);
            EXPECT_TRUE(PrintDefaultRVE);
            EXPECT_FALSE(PrintFullOutput);
            EXPECT_DOUBLE_EQ(RNGSeed, 2.0);
        }
    }
}
//---------------------------------------------------------------------------//
// activedomainsizecalc
//---------------------------------------------------------------------------//
void testcalcZBound_Low_Remelt() {

    int LayerHeight = 10;
    int NumberOfLayers = 10;
    double deltax = 1 * pow(10, -6);
    float ZMin = -0.5 * pow(10, -6);
    float *ZMinLayer = new float[NumberOfLayers];
    for (int layernumber = 0; layernumber < NumberOfLayers; layernumber++) {
        // Set ZMinLayer for each layer to be offset by LayerHeight cells from the previous one (lets solution for
        // both problem types by the same)
        ZMinLayer[layernumber] = ZMin + layernumber * LayerHeight * deltax;
        // Call function for each layernumber, and for simulation types "S" and "R"
        int ZBound_Low_S = calcZBound_Low_Remelt("S", LayerHeight, layernumber, ZMinLayer, ZMin, deltax);
        EXPECT_EQ(ZBound_Low_S, LayerHeight * layernumber);
        int ZBound_Low_R = calcZBound_Low_Remelt("R", LayerHeight, layernumber, ZMinLayer, ZMin, deltax);
        EXPECT_EQ(ZBound_Low_R, LayerHeight * layernumber);
    }
}

void testcalcZBound_High_Remelt() {

    // A separate function is now used for ZBound_High calculation without remelting
    int SpotRadius = 100;
    int LayerHeight = 10;
    float ZMin = 0.5 * pow(10, -6);
    double deltax = 1.0 * pow(10, -6);
    int nz = 101;
    int NumberOfLayers = 10;
    float *ZMaxLayer = new float[NumberOfLayers];
    for (int layernumber = 0; layernumber < NumberOfLayers; layernumber++) {
        // Set ZMaxLayer for each layer to be offset by LayerHeight cells from the previous one, with layer 0 having a
        // ZMax value of ZMin + SpotRadius (lets solution for both problem types be the same)
        ZMaxLayer[layernumber] = ZMin + SpotRadius * deltax + layernumber * LayerHeight * deltax;
        // Call function for each layernumber, and for simulation types "S" and "R"
        int ZBound_Max_S =
            calcZBound_High_Remelt("S", SpotRadius, LayerHeight, layernumber, ZMin, deltax, nz, ZMaxLayer);
        EXPECT_EQ(ZBound_Max_S, SpotRadius + LayerHeight * layernumber);
        int ZBound_Max_R =
            calcZBound_High_Remelt("R", SpotRadius, LayerHeight, layernumber, ZMin, deltax, nz, ZMaxLayer);
        EXPECT_EQ(ZBound_Max_R, SpotRadius + LayerHeight * layernumber);
        // For simulation type C, should be independent of layernumber
        int ZBound_Max_C =
            calcZBound_High_Remelt("C", SpotRadius, LayerHeight, layernumber, ZMin, deltax, nz, ZMaxLayer);
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
void testReadTemperatureData() {

    // Create test data
    double deltax = 1 * pow(10, -6);
    double HT_deltax = 1 * pow(10, -6);
    int HTtoCAratio;
    // Domain size is a 3 by 12 region
    int nx = 3;
    int ny = 12;
    // Write fake OpenFOAM data - only rank 0.
    std::string TestTempFileName = "TestData.txt";
    std::ofstream TestDataFile;
    TestDataFile.open(TestTempFileName);
    TestDataFile << "x, y, z, tm, tl, cr" << std::endl;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            TestDataFile << i * deltax << "," << j * deltax << "," << 0.0 << "," << (float)(i * j) << ","
                         << (float)(i * j + i) << "," << (float)(i * j + j) << std::endl;
        }
    }
    TestDataFile.close();

    // Test each 12 by 3 subdomain
    for (int id = 0; id < 4; id++) {
        int MyYSlices = 3;
        int MyYOffset = 3 * id; // each col is separated from the others by 3 cells
        float YMin = 0.0;
        std::vector<std::string> temp_paths(1);
        temp_paths[0] = TestTempFileName;
        int NumberOfLayers = 1;
        int TempFilesInSeries = 1;
        int *FirstValue = new int[NumberOfLayers];
        int *LastValue = new int[NumberOfLayers];
        unsigned int NumberOfTemperatureDataPoints = 0;

        // Read in data to "RawData"
        std::vector<double> RawData(12);

        float *ZMinLayer = new float[NumberOfLayers];
        float *ZMaxLayer = new float[NumberOfLayers];
        ZMinLayer[0] = 0.0;
        ZMaxLayer[0] = 0.0;
        ReadTemperatureData(id, deltax, HT_deltax, HTtoCAratio, MyYSlices, MyYOffset, YMin, temp_paths, NumberOfLayers,
                            TempFilesInSeries, NumberOfTemperatureDataPoints, RawData, FirstValue, LastValue);

        // Check the results.
        // Does each rank have the right number of temperature data points? Each rank should have six (x,y,z,tm,tl,cr)
        // for each of the 12 cells in the subdomain
        EXPECT_EQ(NumberOfTemperatureDataPoints, 54);
        // Ratio of HT cell size and CA cell size should be 1
        EXPECT_EQ(HTtoCAratio, 1);
        int NumberOfCellsPerRank = 9;
        // Does each rank have the right temperature data values?
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
            ExpectedValues_ThisDataPoint[2] = 0.0;
            int XInt = ExpectedValues_ThisDataPoint[0] / deltax;
            int YInt = ExpectedValues_ThisDataPoint[1] / deltax;
            // Melting time
            ExpectedValues_ThisDataPoint[3] = XInt * YInt;
            // Liquidus time
            ExpectedValues_ThisDataPoint[4] = XInt * YInt + XInt;
            // Cooling rate
            ExpectedValues_ThisDataPoint[5] = XInt * YInt + YInt;
            for (int nn = 0; nn < 6; nn++) {
                EXPECT_DOUBLE_EQ(ExpectedValues_ThisDataPoint[nn], RawData[6 * n + nn]);
            }
        }
    }
}

void testgetTempCoords() {

    // cell size and simulation lower bounds
    double deltax = 0.5;
    float XMin = -5;
    float YMin = 0.0;

    // test - each rank initializes data for layer id, of 4 total layers
    for (int layernumber = 0; layernumber < 4; layernumber++) {
        float *ZMinLayer = new float[4];
        for (int n = 0; n < 4; n++) {
            ZMinLayer[n] = 2 * n;
        }
        int LayerHeight = 10;

        // Fill "RawData" with example temperature field data - 2 data point (x,y,z,tm,tl,cr)
        std::vector<double> RawData = {0.0, 0.0, 3.0, 0.000050, 0.000055, 2000.0,
                                       5.0, 2.0, 2.0, 0.000150, 0.000155, 3500.0};

        // Read and check RawData values against expected results from getTempCoordX,Y,Z,TM,TL,CR functions
        // First data point (i = 0 through i = 5 values)
        int i = 0;
        int XInt = getTempCoordX(i, XMin, deltax, RawData);
        int YInt = getTempCoordY(i, YMin, deltax, RawData);
        int ZInt = getTempCoordZ(i, deltax, RawData, LayerHeight, layernumber, ZMinLayer);
        double TMelting = getTempCoordTM(i, RawData);
        double TLiquidus = getTempCoordTL(i, RawData);
        double CoolingRate = getTempCoordCR(i, RawData);
        EXPECT_EQ(XInt, 10);
        EXPECT_EQ(YInt, 0);
        EXPECT_EQ(ZInt, 6 + 6 * layernumber); // different for each rank, since LayerCounter = rank id
        EXPECT_DOUBLE_EQ(TMelting, 0.000050);
        EXPECT_DOUBLE_EQ(TLiquidus, 0.000055);
        EXPECT_DOUBLE_EQ(CoolingRate, 2000.0);
        // Second data point (i = 6 through i = 11 values)
        i = 6;
        XInt = getTempCoordX(i, XMin, deltax, RawData);
        YInt = getTempCoordY(i, YMin, deltax, RawData);
        ZInt = getTempCoordZ(i, deltax, RawData, LayerHeight, layernumber, ZMinLayer);
        TMelting = getTempCoordTM(i, RawData);
        TLiquidus = getTempCoordTL(i, RawData);
        CoolingRate = getTempCoordCR(i, RawData);
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
TEST(TEST_CATEGORY, fileread_test) { testInputReadFromFile(); }
TEST(TEST_CATEGORY, activedomainsizecalc) {
    testcalcZBound_Low_Remelt();
    testcalcZBound_High_Remelt();
    testcalcnzActive();
    testcalcLocalActiveDomainSize();
}
TEST(TEST_CATEGORY, temperature_init_test) {
    testReadTemperatureData(); // reading temperature data is performed in the same manner with and without remelting
    testgetTempCoords();
}
} // end namespace Test
