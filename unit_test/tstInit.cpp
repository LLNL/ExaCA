
#include <Kokkos_Core.hpp>

#include "CAinitialize.hpp"

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
        int DecompositionStrategy, TempFilesInSeries, NumberOfLayers, LayerHeight, nx, ny, nz, PrintDebug, NSpotsX,
            NSpotsY, SpotOffset, SpotRadius, TimeSeriesInc, RVESize;
        float SubstrateGrainSpacing;
        double AConst, BConst, CConst, DConst, FreezingRange, deltax, NMax, dTN, dTsigma, HT_deltax, deltat, G, R,
            FractSurfaceSitesActive, RNGSeed, PowderDensity;
        bool RemeltingYN, PrintMisorientation, PrintFinalUndercoolingVals, PrintFullOutput, PrintTimeSeries,
            UseSubstrateFile, PrintIdleTimeSeriesFrames, PrintDefaultRVE = false, BaseplateThroughPowder;
        std::string SimulationType, OutputFile, GrainOrientationFile, temppath, tempfile, SubstrateFileName,
            PathToOutput;
        std::vector<std::string> temp_paths;
        InputReadFromFile(
            id, FileName, SimulationType, DecompositionStrategy, AConst, BConst, CConst, DConst, FreezingRange, deltax,
            NMax, dTN, dTsigma, OutputFile, GrainOrientationFile, TempFilesInSeries, temp_paths, HT_deltax, RemeltingYN,
            deltat, NumberOfLayers, LayerHeight, SubstrateFileName, SubstrateGrainSpacing, UseSubstrateFile, G, R, nx,
            ny, nz, FractSurfaceSitesActive, PathToOutput, PrintDebug, PrintMisorientation, PrintFinalUndercoolingVals,
            PrintFullOutput, NSpotsX, NSpotsY, SpotOffset, SpotRadius, PrintTimeSeries, TimeSeriesInc,
            PrintIdleTimeSeriesFrames, PrintDefaultRVE, RNGSeed, BaseplateThroughPowder, PowderDensity, RVESize);

        // Check the results
        // The existence of the specified orientation, substrate, and temperature filenames was already checked within
        // InputReadFromFile
        // These should be the same for all 3 test problems
        EXPECT_EQ(DecompositionStrategy, 1);
        EXPECT_DOUBLE_EQ(AConst, -0.00000010302);
        EXPECT_DOUBLE_EQ(BConst, 0.00010533);
        EXPECT_DOUBLE_EQ(CConst, 0.0022196);
        EXPECT_DOUBLE_EQ(DConst, 0);
        EXPECT_DOUBLE_EQ(FreezingRange, 210);
        EXPECT_DOUBLE_EQ(deltax, 1.0 * pow(10, -6));
        EXPECT_DOUBLE_EQ(NMax, 1.0 * pow(10, 13));
        EXPECT_DOUBLE_EQ(dTN, 5.0);
        EXPECT_DOUBLE_EQ(dTsigma, 0.5);
        EXPECT_EQ(PrintDebug, 0);
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

    int MyXSlices = 5;
    int MyYSlices = 4;
    int nzActive = 10;
    int LocalActiveDomainSize = calcLocalActiveDomainSize(MyXSlices, MyYSlices, nzActive);
    EXPECT_EQ(LocalActiveDomainSize, 10 * 5 * 4);
}

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, fileread_test) { testInputReadFromFile(); }
TEST(TEST_CATEGORY, activedomainsizecalc) {
    testcalcnzActive();
    testcalcLocalActiveDomainSize();
}
} // end namespace Test
