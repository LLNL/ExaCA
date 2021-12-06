
#include "CAinitialize.hpp"

#include <gtest/gtest.h>

#include "mpi.h"

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
void testInputReadFromFile() {

    int id;
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Three input files - one of each type
    // Inp_DirSolidification.txt and Inp_SpotMelt.txt were installed from the examples directory
    // Since no temperature files exist in the repo, and there is no ability to write temperature files to a different
    // directory ( would need examples/Temperatures) using the C++11 standard, a dummy input file is written and parsed
    // to test an example problem that uses temperature data from a file.
    std::vector<std::string> InputFilenames = {"Inp_DirSolidification.txt", "Inp_SpotMelt.txt",
                                               "Inp_TemperatureTest.txt"};
    if (id == 0) {
        // Write dummy input file for using read temperature data
        std::ofstream TestDataFile;
        TestDataFile.open(InputFilenames[2]);
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
        TestDataFile << "File of grain orientations: GrainOrientationVectors_Robert.csv" << std::endl;
        TestDataFile << "Print file of grain misorientation values: Y" << std::endl;
        TestDataFile << "Print file of all ExaCA data: N" << std::endl;
        TestDataFile << "Time step: 1.5" << std::endl;
        TestDataFile << "Temperature filename(s): DummyTemperature.txt" << std::endl;
        TestDataFile << "Number of temperature files: 2" << std::endl;
        TestDataFile << "Number of layers: 2" << std::endl;
        TestDataFile << "Offset between layers: 1" << std::endl;
        TestDataFile << "Extra set of wall cells in lateral domain directions: N" << std::endl;
        TestDataFile << "Substrate filename: DummySubstrate.txt" << std::endl;
        TestDataFile << "Heat transport data mesh size: 12" << std::endl;
        TestDataFile << "Path to temperature file(s): ./" << std::endl;
        TestDataFile.close();

        // Create test temperature files "1DummyTemperature.txt", "2DummyTemperature.txt"
        std::ofstream TestTemp1, TestTemp2;
        TestTemp1.open("1DummyTemperature.txt");
        TestTemp2.open("2DummyTemperature.txt");
        TestTemp1 << "X" << std::endl;
        TestTemp2 << "X" << std::endl;
        TestTemp1.close();
        TestTemp2.close();

        // Create test substrate file "DummySubstrate.txt"
        std::ofstream TestSub;
        TestSub.open("DummySubstrate.txt");
        TestSub << "X" << std::endl;
        TestSub.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Read and parse each input file
    for (auto FileName : InputFilenames) {
        int DecompositionStrategy, TempFilesInSeries, NumberOfLayers, LayerHeight, nx, ny, nz, PrintDebug, NSpotsX,
            NSpotsY, SpotOffset, SpotRadius, TimeSeriesInc;
        float SubstrateGrainSpacing;
        double AConst, BConst, CConst, DConst, FreezingRange, deltax, NMax, dTN, dTsigma, HT_deltax, deltat, G, R,
            FractSurfaceSitesActive;
        bool RemeltingYN, ExtraWalls, PrintMisorientation, PrintFullOutput, PrintTimeSeries, UseSubstrateFile,
            PrintIdleTimeSeriesFrames;
        std::string SimulationType, OutputFile, GrainOrientationFile, temppath, tempfile, SubstrateFileName,
            PathToOutput;
        std::vector<std::string> temp_paths;
        InputReadFromFile(id, FileName, SimulationType, DecompositionStrategy, AConst, BConst, CConst, DConst,
                          FreezingRange, deltax, NMax, dTN, dTsigma, OutputFile, GrainOrientationFile, temppath,
                          tempfile, TempFilesInSeries, temp_paths, ExtraWalls, HT_deltax, RemeltingYN, deltat,
                          NumberOfLayers, LayerHeight, SubstrateFileName, SubstrateGrainSpacing, UseSubstrateFile, G, R,
                          nx, ny, nz, FractSurfaceSitesActive, PathToOutput, PrintDebug, PrintMisorientation,
                          PrintFullOutput, NSpotsX, NSpotsY, SpotOffset, SpotRadius, PrintTimeSeries, TimeSeriesInc,
                          PrintIdleTimeSeriesFrames);

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
            EXPECT_EQ(nx, 202);                           // +2 from input due to wall cells
            EXPECT_EQ(ny, 202);                           // +2 from input due to wall cells
            EXPECT_EQ(nz, 201);                           // +1 from input due to wall cells
            EXPECT_DOUBLE_EQ(FractSurfaceSitesActive, 0.08);
            EXPECT_TRUE(OutputFile == "TestProblemDirS");
            EXPECT_TRUE(PrintMisorientation);
            EXPECT_TRUE(PrintFullOutput);
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
            EXPECT_FLOAT_EQ(SubstrateGrainSpacing, 25.0);
            EXPECT_TRUE(OutputFile == "TestProblemSpot");
            EXPECT_TRUE(PrintFullOutput);
        }
        else if (FileName == "Inp_TemperatureTest.txt") {
            EXPECT_DOUBLE_EQ(deltat, 1.5 * pow(10, -6));
            EXPECT_EQ(TempFilesInSeries, 2);
            EXPECT_EQ(NumberOfLayers, 2);
            EXPECT_EQ(LayerHeight, 1);
            EXPECT_FALSE(ExtraWalls);
            EXPECT_TRUE(UseSubstrateFile);
            EXPECT_DOUBLE_EQ(HT_deltax, 12.0 * pow(10, -6));
            EXPECT_TRUE(OutputFile == "Test");
            EXPECT_TRUE(PrintMisorientation);
            EXPECT_FALSE(PrintFullOutput);
        }
    }
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, fileread_test) { testInputReadFromFile(); }

} // end namespace Test
