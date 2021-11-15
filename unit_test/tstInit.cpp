
#include <Kokkos_Core.hpp>

#include "CAinitialize.hpp"

#include <gtest/gtest.h>

#include "mpi.h"

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
void testReadTemperatureData() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Create test data - amount of data depends on the number of ranks
    int MyXOffset, MyYOffset;
    double deltax = 1 * pow(10, -6);
    double HT_deltax = 1 * pow(10, -6);
    int MyXSlices = 6;
    int MyYSlices = 6;
    // Two "rows" of different MPI rank subdomains (Y direction)
    // Varied numbers of "columns" depending on total number of MPI ranks
    int ProcRow = id % 2;
    int ProcCol = id / 2;
    MyXOffset = 6 * ProcCol + ProcCol;
    if (ProcRow == 0)
        MyYOffset = 0;
    else
        MyYOffset = 6;
    int nx = 50;
    int ny = 50;
    float XMin = 0.0;
    float YMin = 0.0;
    std::string TestTempFileName = "TestData.txt";
    std::vector<std::string> temp_paths(1);
    temp_paths[0] = TestTempFileName;
    int NumberOfLayers = 1;
    int TempFilesInSeries = 1;
    unsigned int NumberOfTemperatureDataPoints = 0;
    std::vector<double> RawData(12);
    int *FirstValue = new int[NumberOfLayers];
    int *LastValue = new int[NumberOfLayers];
    // If np = 1 or 2, one "column" of MPI ranks (rank 0 in row 0, rank 1 in row 1)
    // If np = 3 or 4, two "columns" of MPI ranks (even ranks in row 0, odd ranks in row 1)
    int GlobalProcCols = (np + 1) / 2;
    int AllRanks_NDataPoints = 6 * GlobalProcCols;
    std::vector<std::vector<double>> TestData(AllRanks_NDataPoints, std::vector<double>(6));
    // Fill "TestData" with values to be read/stored on different MPI ranks
    // Temperature data will be printed at two Y values: 4 * deltax and 6 * deltax
    // Temperature data will be printed at X = 2 * deltax, 3 * deltax, 4 * deltax if np = 1 or 2
    // If np > 2 (multiple "columns" of MPI ranks), an additional 3 data points will be printed for each "column"
    // For example, if np = 4, there will be additional data at X = 5 * deltax, X = 6 * deltax, X = 7 * deltax
    // Even ranks should only read and store data at Y = 4 * deltax, and at X = (ProcCol + 2) * deltax through X =
    // (ProcCol + 5) * deltax Odd ranks should store data at Y = 4 * deltax and Y = 6 * deltax, and at the same X
    // coordinates as (id - 1)
    for (int i = 0; i < GlobalProcCols; i++) {
        for (int j = 0; j < 3; j++) {
            int Row1Loc = 3 * i + j;                    // This is the Y = 4 * deltax data
            int Row2Loc = 3 * GlobalProcCols + Row1Loc; // This is the Y = 6 * deltax data
            TestData[Row1Loc][0] = deltax * (2 + (3 * i) + j);
            TestData[Row1Loc][1] = 4 * deltax;
            TestData[Row1Loc][2] = 0;
            TestData[Row1Loc][3] = 10 * i + j;
            TestData[Row1Loc][4] = TestData[Row1Loc][3] + 2.5;
            TestData[Row1Loc][5] = i * i;
            TestData[Row2Loc][0] = deltax * (2 + (3 * i) + j);
            TestData[Row2Loc][1] = 6 * deltax;
            TestData[Row2Loc][2] = 0;
            TestData[Row2Loc][3] = 10 * i + j;
            TestData[Row2Loc][4] = TestData[Row2Loc][3] + 2.5;
            TestData[Row2Loc][5] = i * i;
        }
    }
    // Write fake OpenFOAM data - only rank 0.
    if (id == 0) {
        std::ofstream TestDataFile;
        TestDataFile.open(TestTempFileName);
        TestDataFile << "x, y, z, tm, tl, cr" << std::endl;
        for (int i = 0; i < AllRanks_NDataPoints; i++) {
            TestDataFile << TestData[i][0] << "," << TestData[i][1] << "," << TestData[i][2] << "," << TestData[i][3]
                         << "," << TestData[i][4] << "," << TestData[i][5] << std::endl;
        }
        TestDataFile.close();
    }
    // Wait for data to be printed before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    // Read in data to "RawData"
    ReadTemperatureData(deltax, HT_deltax, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nx, ny, XMin, YMin, temp_paths,
                        NumberOfLayers, TempFilesInSeries, NumberOfTemperatureDataPoints, RawData, FirstValue,
                        LastValue);
    // Check the results.
    // Does each rank have the right number of temperature data points?
    if (ProcRow == 0) {
        // Even ranks - these stored data at 3 X values, and 1 Y value (Y = 4 * deltax)
        // As each X,Y,Z coordinate has 6 associated data points (X, Y, Z, tm, ts, cr), this rank should have 3 * 6 = 18
        // temperature data points stored
        EXPECT_EQ(NumberOfTemperatureDataPoints, 18);
        // TestData values should match RawData
        int CompValues[3] = {0 + 3 * ProcCol, 1 + 3 * ProcCol, 2 + 3 * ProcCol};
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 6; j++) {
                EXPECT_DOUBLE_EQ(RawData[6 * i + j], TestData[CompValues[i]][j]);
            }
        }
    }
    else {
        // Odd ranks - these stored data at 3 X values, and 2 Y values (Y = 4 * deltax and Y = 6 * deltax)
        // As each X,Y,Z coordinate has 6 associated data points (X, Y, Z, tm, ts, cr), this rank should have 6 * 6 = 36
        // temperature data points stored
        EXPECT_EQ(NumberOfTemperatureDataPoints, 36);
        // TestData values should match RawData
        int CompValues[6] = {0 + 3 * ProcCol,
                             1 + 3 * ProcCol,
                             2 + 3 * ProcCol,
                             3 * GlobalProcCols + 0 + 3 * ProcCol,
                             3 * GlobalProcCols + 1 + 3 * ProcCol,
                             3 * GlobalProcCols + 2 + 3 * ProcCol};
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                EXPECT_DOUBLE_EQ(RawData[6 * i + j], TestData[CompValues[i]][j]);
            }
        }
    }
}

void testInputReadFromFile() {

    int id;
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Three input files - one of each type
    std::vector<std::string> InputFilenames = {"TestC.txt", "TestS.txt", "TestR.txt"};
    if (id == 0) {
        // Write dummy input file for directional solidification (problem type C, iteration 0), spot melt (problem type
        // S, iteration 1), and read temperature data (problem type R, iteration 2)
        for (int FileNumber = 0; FileNumber < 3; FileNumber++) {
            std::ofstream TestDataFile;
            TestDataFile.open(InputFilenames[FileNumber]);
            // Write required inputs to all files
            TestDataFile << "Test problem data set" << std::endl;
            TestDataFile << "*****" << std::endl;
            if (FileNumber == 0)
                TestDataFile << "Problem type: C" << std::endl;
            else if (FileNumber == 1)
                TestDataFile << "Problem type: S" << std::endl;
            else if (FileNumber == 2)
                TestDataFile << "Problem type: R" << std::endl;
            TestDataFile << "Decomposition strategy: 2" << std::endl;
            TestDataFile << "Material: Inconel625" << std::endl;
            TestDataFile << "Cell size: 1" << std::endl;
            TestDataFile << "Heterogeneous nucleation density: 6" << std::endl;
            TestDataFile << "Mean nucleation undercooling: 3" << std::endl;
            TestDataFile << "Standard deviation of nucleation undercooling: 0.11" << std::endl;
            TestDataFile << "Path to output: ExaCA" << std::endl;
            TestDataFile << "Output file base name: Test" << std::endl;
            TestDataFile << "File of grain orientations: GrainOrientationVectors_Robert.csv" << std::endl;
            TestDataFile << "Print file of grain misorientation values: Y" << std::endl;
            TestDataFile << "Print file of all ExaCA data: N" << std::endl;
            // Print optional file inputs for 2 of the 3 problems, but not the 3rd
            // to check for proper behavior in both situations
            if ((FileNumber == 0) || (FileNumber == 1)) {
                // Print optional inputs, and problem type specific inputs
                TestDataFile << "Debug check (reduced): Y" << std::endl;
                TestDataFile << "Debug check (extensive): N" << std::endl;
                TestDataFile << "Print intermediate output frames: Y" << std::endl;
                TestDataFile << "Increment to separate frames: 2" << std::endl;
                TestDataFile << "Intermediate output even if system is unchanged from previous state: N" << std::endl;
            }
            // Problem type specific required inputs/optional inputs
            if (FileNumber == 0) {
                TestDataFile << "Thermal gradient: 1000" << std::endl;
                TestDataFile << "Cooling rate: 100" << std::endl;
                TestDataFile << "Time step ratio (from steady-state Velocity): 25" << std::endl;
                TestDataFile << "Domain size in x: 10" << std::endl;
                TestDataFile << "Domain size in y: 15" << std::endl;
                TestDataFile << "Domain size in z: 20" << std::endl;
                TestDataFile << "Fraction surface sites active: 0.35" << std::endl;
            }
            else if (FileNumber == 1) {
                TestDataFile << "Thermal gradient: 1000" << std::endl;
                TestDataFile << "Cooling rate: 100" << std::endl;
                TestDataFile << "Time step ratio (from steady-state Velocity): 25" << std::endl;
                TestDataFile << "Num spots in x: 1" << std::endl;
                TestDataFile << "Num spots in y: 5" << std::endl;
                TestDataFile << "Offset between spot centers: 5" << std::endl;
                TestDataFile << "Radii of spots: 10" << std::endl;
                TestDataFile << "Number of layers: 2" << std::endl;
                TestDataFile << "Offset between layers: 1" << std::endl;
                TestDataFile << "Substrate grain spacing: 10" << std::endl;
            }
            else if (FileNumber == 2) {
                TestDataFile << "Time step: 1.5" << std::endl;
                TestDataFile << "Temperature filename(s): DummyTemperature.txt" << std::endl;
                TestDataFile << "Number of temperature files: 2" << std::endl;
                TestDataFile << "Number of layers: 2" << std::endl;
                TestDataFile << "Offset between layers: 1" << std::endl;
                TestDataFile << "Extra set of wall cells in lateral domain directions: N" << std::endl;
                TestDataFile << "Substrate filename: DummySubstrate.txt" << std::endl;
                TestDataFile << "Heat transport data mesh size: 12" << std::endl;
                TestDataFile << "Path to temperature file(s): ./" << std::endl;
            }
            TestDataFile.close();
        }
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
    for (int FileNumber = 0; FileNumber < 3; FileNumber++) {
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
        InputReadFromFile(id, InputFilenames[FileNumber], SimulationType, DecompositionStrategy, AConst, BConst, CConst,
                          DConst, FreezingRange, deltax, NMax, dTN, dTsigma, OutputFile, GrainOrientationFile, temppath,
                          tempfile, TempFilesInSeries, temp_paths, ExtraWalls, HT_deltax, RemeltingYN, deltat,
                          NumberOfLayers, LayerHeight, SubstrateFileName, SubstrateGrainSpacing, UseSubstrateFile, G, R,
                          nx, ny, nz, FractSurfaceSitesActive, PathToOutput, PrintDebug, PrintMisorientation,
                          PrintFullOutput, NSpotsX, NSpotsY, SpotOffset, SpotRadius, PrintTimeSeries, TimeSeriesInc,
                          PrintIdleTimeSeriesFrames);

        // Check the results - required inputs for all problems
        // The existence of the specified orientation, substrate, and temperature filenames was already checked within
        // InputReadFromFile
        EXPECT_EQ(DecompositionStrategy, 2);
        EXPECT_DOUBLE_EQ(AConst, -0.00000010302);
        EXPECT_DOUBLE_EQ(BConst, 0.00010533);
        EXPECT_DOUBLE_EQ(CConst, 0.0022196);
        EXPECT_DOUBLE_EQ(DConst, 0);
        EXPECT_DOUBLE_EQ(FreezingRange, 210);
        EXPECT_DOUBLE_EQ(deltax, 1.0 * pow(10, -6));
        EXPECT_DOUBLE_EQ(NMax, 6.0 * pow(10, 12));
        EXPECT_DOUBLE_EQ(dTN, 3.0);
        EXPECT_DOUBLE_EQ(dTsigma, 0.11);
        EXPECT_TRUE(OutputFile == "Test");
        EXPECT_TRUE(PrintMisorientation);
        EXPECT_FALSE(PrintFullOutput);
        // Input file specific results to check
        if (FileNumber == 0) {
            EXPECT_EQ(PrintDebug, 1);
            EXPECT_TRUE(PrintTimeSeries);
            EXPECT_EQ(TimeSeriesInc, 5);
            EXPECT_FALSE(PrintIdleTimeSeriesFrames);
            EXPECT_DOUBLE_EQ(G, 1000.0);
            EXPECT_DOUBLE_EQ(R, 100.0);
            EXPECT_DOUBLE_EQ(deltat, 0.4 * pow(10, -6)); // based on N, deltax, G, and R from input file
            EXPECT_EQ(nx, 12);                           // +2 from input due to wall cells
            EXPECT_EQ(ny, 17);                           // +2 from input due to wall cells
            EXPECT_EQ(nz, 21);                           // +1 from input due to wall cells
            EXPECT_DOUBLE_EQ(FractSurfaceSitesActive, 0.35);
        }
        else if (FileNumber == 1) {
            EXPECT_EQ(PrintDebug, 1);
            EXPECT_TRUE(PrintTimeSeries);
            EXPECT_EQ(TimeSeriesInc, 5);
            EXPECT_FALSE(PrintIdleTimeSeriesFrames);
            EXPECT_DOUBLE_EQ(G, 1000.0);
            EXPECT_DOUBLE_EQ(R, 100.0);
            EXPECT_DOUBLE_EQ(deltat, 0.4 * pow(10, -6)); // based on N, deltax, G, and R from input file
            EXPECT_EQ(NSpotsX, 1);
            EXPECT_EQ(NSpotsY, 5);
            EXPECT_EQ(SpotOffset, 5);
            EXPECT_EQ(SpotRadius, 10);
            EXPECT_EQ(NumberOfLayers, 2);
            EXPECT_EQ(LayerHeight, 1);
            EXPECT_FALSE(UseSubstrateFile);
            EXPECT_FLOAT_EQ(SubstrateGrainSpacing, 10.0);
        }
        else if (FileNumber == 2) {
            EXPECT_DOUBLE_EQ(deltat, 1.5 * pow(10, -6));
            EXPECT_EQ(TempFilesInSeries, 2);
            EXPECT_EQ(NumberOfLayers, 2);
            EXPECT_EQ(LayerHeight, 1);
            EXPECT_FALSE(ExtraWalls);
            EXPECT_TRUE(UseSubstrateFile);
            EXPECT_DOUBLE_EQ(HT_deltax, 12.0 * pow(10, -6));
        }
    }
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, temperature_init_test) { testReadTemperatureData(); }
TEST(TEST_CATEGORY, fileread_test) { testInputReadFromFile(); }

} // end namespace Test
