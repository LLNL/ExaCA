
#include <Kokkos_Core.hpp>

#include "CAinitialize.hpp"
#include "CAtypes.hpp"

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
    double FreezingRange = 210.0;
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
    int LayerHeight = 1;
    int TempFilesInSeries = 1;
    unsigned int NumberOfTemperatureDataPoints = 0;
    std::vector<double> RawData(12);
    int *FirstValue = new int[NumberOfLayers];
    int *LastValue = new int[NumberOfLayers];
    float *ZMinLayer = new float[NumberOfLayers];
    float *ZMaxLayer = new float[NumberOfLayers];
    ViewI_H MaxSolidificationEvents_H(Kokkos::ViewAllocateWithoutInitializing("MaxSolidificationEvents"),
                                      NumberOfLayers);
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
    ReadTemperatureData(id, false, MaxSolidificationEvents_H, MyXSlices, MyYSlices, MyXOffset, MyYOffset, deltax,
                        HT_deltax, nx, ny, XMin, YMin, temp_paths, FreezingRange, LayerHeight, NumberOfLayers,
                        TempFilesInSeries, NumberOfTemperatureDataPoints, ZMinLayer, ZMaxLayer, FirstValue, LastValue,
                        RawData);
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

void testSubstrateInit_FromGrainSpacing() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Create test data
    int nz = 5;
    int nzActive = 4;
    int nx = 3;
    int MyXSlices = 3;
    int MyXOffset = 0;
    // Each rank is assigned a different portion of the domain in Y
    int ny = 4;
    int MyYSlices = 2;
    int MyYOffset;
    if (id % 2 == 0)
        MyYOffset = 0;
    else
        MyYOffset = 2;
    double deltax = 1 * pow(10, -6);
    int LocalActiveDomainSize = MyXSlices * MyYSlices * nzActive;
    int LocalDomainSize = MyXSlices * MyYSlices * nz;
    // There are 48 total cells in this domain (nx * ny * nz)
    // Even ranks have 30, and odd ranks have the other 30
    double SubstrateGrainSpacing = 1.95; // This grain spacing ensures that there will be 3 substrate grains
    ViewI_H GrainID_H(Kokkos::ViewAllocateWithoutInitializing("GrainID"), LocalDomainSize);
    ViewI_H CritTimeStep_H(Kokkos::ViewAllocateWithoutInitializing("CritTimeStep"), LocalDomainSize);
    // Initialize GrainID to 0 everywhere
    // Initialize CritTimeStep to 0 everywhere to allow substrate generation for these cells
    for (int i = 0; i < MyXSlices * MyYSlices * nz; i++) {
        CritTimeStep_H(i) = 0;
        GrainID_H(i) = 0;
    }
    CritTimeStep_H(13) = 1; // Let 1 cell have a non-zero CritTimeStep - without remelting, this cell should end up
                            // being assigned a GrainID of 0, while the others have GrainIDs of 1, 2, or 3
    SubstrateInit_FromGrainSpacing(false, SubstrateGrainSpacing, nx, ny, nz, nzActive, MyXSlices, MyYSlices, MyXOffset,
                                   MyYOffset, LocalActiveDomainSize, id, np, deltax, GrainID_H, CritTimeStep_H);

    // Check the results
    if (id % 2 == 0) {
        // Expected GrainID values: Cells 0-5 are walls, should have Grain ID 0
        // Cell 13 has CritTimeStep > 1, should also have GrainID 0
        // Other cells should have GrainID 1, 2, or 3, depending on the closest grain center
        ViewI_H ExpectedGID(Kokkos::ViewAllocateWithoutInitializing("ExGrainID"), LocalDomainSize);
        for (int i = 0; i < 6; i++) {
            ExpectedGID(i) = 0;
        }
        for (int i = 6; i < 19; i++) {
            ExpectedGID(i) = 1;
        }
        ExpectedGID(13) = 0;
        ExpectedGID(19) = 2;
        ExpectedGID(20) = 1;
        ExpectedGID(21) = 2;
        ExpectedGID(22) = 1;
        ExpectedGID(23) = 2;
        for (int i = 24; i < 30; i++) {
            ExpectedGID(i) = 3;
        }
        for (int i = 0; i < LocalDomainSize; i++) {
            EXPECT_EQ(GrainID_H(i), ExpectedGID(i));
        }
    }
    else {
        // Expected GrainID values: Cells 0-5 are walls, should have Grain ID 0
        // Cell 13 has CritTimeStep > 1, should also have GrainID 0
        // Other cells should have GrainID 1, 2, or 3, depending on the closest grain center
        ViewI_H ExpectedGID(Kokkos::ViewAllocateWithoutInitializing("ExGrainID"), LocalDomainSize);
        for (int i = 0; i < 6; i++) {
            ExpectedGID(i) = 0;
        }
        ExpectedGID(6) = 1;
        for (int i = 7; i < 24; i++) {
            ExpectedGID(i) = 2;
        }
        ExpectedGID(8) = 1;
        ExpectedGID(10) = 1;
        ExpectedGID(13) = 0;
        for (int i = 24; i < 30; i++) {
            ExpectedGID(i) = 3;
        }
        for (int i = 0; i < LocalDomainSize; i++) {
            EXPECT_EQ(GrainID_H(i), ExpectedGID(i));
        }
    }
}

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, kokkos_init_test) { testReadTemperatureData(); testSubstrateInit_FromGrainSpacing();}
} // end namespace Test
