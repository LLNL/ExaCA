
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
            std::cout << TestData[i][0] << "," << TestData[i][1] << "," << TestData[i][2] << "," << TestData[i][3]
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

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, temperature_init_test) { testReadTemperatureData(); }

} // end namespace Test
