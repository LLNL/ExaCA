
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
    // Create test data
    double deltax = 1 * pow(10, -6);
    double HT_deltax = 1 * pow(10, -6);
    int HTtoCAratio;
    // Domain size in y is variable based on the number of ranks
    // Each MPI rank contains a 3 by 3 subdomain
    int nx = 6;
    int ny = 3 * (1 + np / 2);
    // Each rank gets 1/(number of ranks) of the overall domain
    int MyXSlices = 3;
    int MyYSlices = 3;
    int ProcRow, ProcCol, MyXOffset, MyYOffset;
    if (id % 2 == 0) {
        ProcRow = 0;   // even ranks
        MyXOffset = 0; // no offset in X: contains X = 0-3
    }
    else {
        ProcRow = 1;   // odd ranks
        MyXOffset = 3; // X = 3 offset, contains X = 4-6
    }
    ProcCol = id / 2;        // Ranks 0-1 in 0th col, 2-3 in 1st col, 4-5 in 2nd col, etc
    MyYOffset = 3 * ProcCol; // each col is separated from the others by 3 cells
    float XMin = 0.0;
    float YMin = 0.0;
    std::string TestTempFileName = "TestData.txt";
    std::vector<std::string> temp_paths(1);
    temp_paths[0] = TestTempFileName;
    int NumberOfLayers = 1;
    int TempFilesInSeries = 1;
    int *FirstValue = new int[NumberOfLayers];
    int *LastValue = new int[NumberOfLayers];
    unsigned int NumberOfTemperatureDataPoints = 0;
    // Write fake OpenFOAM data - only rank 0.
    if (id == 0) {
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
    }
    // Wait for data to be printed before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    // Read in data to "RawData"
    std::vector<double> RawData(12);

    float *ZMinLayer = new float[NumberOfLayers];
    float *ZMaxLayer = new float[NumberOfLayers];
    ZMinLayer[0] = 0.0;
    ZMaxLayer[0] = 0.0;
    ReadTemperatureData(id, deltax, HT_deltax, HTtoCAratio, MyXSlices, MyYSlices, MyXOffset, MyYOffset, XMin, YMin,
                        temp_paths, NumberOfLayers, TempFilesInSeries, NumberOfTemperatureDataPoints, RawData,
                        FirstValue, LastValue);

    // Check the results.
    // Does each rank have the right number of temperature data points? Each rank should have six (x,y,z,tm,tl,cr) for
    // each of the 9 cells in the subdomain
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
        ExpectedValues_ThisDataPoint[0] = (CARow + 3 * ProcRow) * deltax;
        // Y Coordinate
        ExpectedValues_ThisDataPoint[1] = (CACol + 3 * ProcCol) * deltax;
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

void testcalcnzActive() {

    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    int ZBound_Low = 5;
    int NumberOfLayers = 10;
    for (int layernumber = 0; layernumber < NumberOfLayers; layernumber++) {
        int ZBound_High = 6 + id + layernumber;
        int nzActive = calcnzActive(ZBound_Low, ZBound_High, id, layernumber);
        EXPECT_EQ(nzActive, 2 + id + layernumber);
    }
}

void testcalcLocalActiveDomainSize() {

    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    int MyXSlices = 5 + id;
    int MyYSlices = id;
    int nzActive = 10;
    int LocalActiveDomainSize = calcLocalActiveDomainSize(MyXSlices, MyYSlices, nzActive);
    EXPECT_EQ(LocalActiveDomainSize, 10 * id * id + 50 * id);
}

void testgetTempCoords() {

    int id, np;
    // Get individual process id
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    // cell size and simulation lower bounds
    double deltax = 0.5;
    float XMin = -5;
    float YMin = 0.0;

    // test - each rank initializes data for layer "id", of "np" total layers
    int LayerCounter = id;
    float *ZMinLayer = new float[np];
    for (int n = 0; n < np; n++) {
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
    int ZInt = getTempCoordZ(i, deltax, RawData, LayerHeight, LayerCounter, ZMinLayer);
    double TMelting = getTempCoordTM(i, RawData);
    double TLiquidus = getTempCoordTL(i, RawData);
    double CoolingRate = getTempCoordCR(i, RawData);
    EXPECT_EQ(XInt, 10);
    EXPECT_EQ(YInt, 0);
    EXPECT_EQ(ZInt, 6 + 6 * id); // different for each rank, since LayerCounter = rank id
    EXPECT_DOUBLE_EQ(TMelting, 0.000050);
    EXPECT_DOUBLE_EQ(TLiquidus, 0.000055);
    EXPECT_DOUBLE_EQ(CoolingRate, 2000.0);
    // Second data point (i = 6 through i = 11 values)
    i = 6;
    XInt = getTempCoordX(i, XMin, deltax, RawData);
    YInt = getTempCoordY(i, YMin, deltax, RawData);
    ZInt = getTempCoordZ(i, deltax, RawData, LayerHeight, LayerCounter, ZMinLayer);
    TMelting = getTempCoordTM(i, RawData);
    TLiquidus = getTempCoordTL(i, RawData);
    CoolingRate = getTempCoordCR(i, RawData);
    EXPECT_EQ(XInt, 20);
    EXPECT_EQ(YInt, 4);
    EXPECT_EQ(ZInt, 4 + 6 * id); // different for each rank, since LayerCounter = rank id
    EXPECT_DOUBLE_EQ(TMelting, 0.000150);
    EXPECT_DOUBLE_EQ(TLiquidus, 0.000155);
    EXPECT_DOUBLE_EQ(CoolingRate, 3500.0);
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, activedomainsizecalc) {
    testcalcnzActive();
    testcalcLocalActiveDomainSize();
}
TEST(TEST_CATEGORY, temperature_init_test) {
    testReadTemperatureData(); // reading temperature data is performed in the same manner with and without remelting
    testgetTempCoords();
}
} // end namespace Test
