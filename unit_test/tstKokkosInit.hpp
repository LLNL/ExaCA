#include <Kokkos_Core.hpp>

#include "CAfunctions.hpp"
#include "CAinitialize.hpp"
#include "CAparsefiles.hpp"
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
//---------------------------------------------------------------------------//
// activedomainsizecalc
//---------------------------------------------------------------------------//
void testcalcZBound_Low(bool RemeltingYN) {

    if (!(RemeltingYN)) {
        ViewF_H ZMinLayer_Host(Kokkos::ViewAllocateWithoutInitializing("ZMinLayer_Host"), 1);
        int ZBound_Low_C = calcZBound_Low(RemeltingYN, "C", 0, 0, ZMinLayer_Host, 0, 0);
        // Should have just set the lower bound to zero, other inputs don't matter
        EXPECT_EQ(ZBound_Low_C, 0);
    }
    else {
        int LayerHeight = 10;
        int NumberOfLayers = 10;
        double deltax = 1 * pow(10, -6);
        float ZMin = -0.5 * pow(10, -6);
        ViewF_H ZMinLayer_Host(Kokkos::ViewAllocateWithoutInitializing("ZMinLayer_Host"), NumberOfLayers);
        for (int layernumber = 0; layernumber < NumberOfLayers; layernumber++) {
            // Set ZMinLayer for each layer to be offset by LayerHeight cells from the previous one (lets solution for
            // both problem types by the same)
            ZMinLayer_Host(layernumber) = ZMin + layernumber * LayerHeight * deltax;
            // Call function for each layernumber, and for simulation types "S" and "R"
            int ZBound_Low_S = calcZBound_Low(RemeltingYN, "S", LayerHeight, layernumber, ZMinLayer_Host, ZMin, deltax);
            EXPECT_EQ(ZBound_Low_S, LayerHeight * layernumber);
            int ZBound_Low_R = calcZBound_Low(RemeltingYN, "R", LayerHeight, layernumber, ZMinLayer_Host, ZMin, deltax);
            EXPECT_EQ(ZBound_Low_R, LayerHeight * layernumber);
        }
    }
}

void testcalcZBound_High() {

    int SpotRadius = 100;
    int LayerHeight = 10;
    float ZMin = 0.5 * pow(10, -6);
    double deltax = 1.0 * pow(10, -6);
    int nz = 101;
    int NumberOfLayers = 10;
    ViewF_H ZMaxLayer_Host(Kokkos::ViewAllocateWithoutInitializing("ZMaxLayer_Host"), NumberOfLayers);
    for (int layernumber = 0; layernumber < NumberOfLayers; layernumber++) {
        // Set ZMaxLayer for each layer to be offset by LayerHeight cells from the previous one, with layer 0 having a
        // ZMax value of ZMin + SpotRadius (lets solution for both problem types be the same)
        ZMaxLayer_Host(layernumber) = ZMin + SpotRadius * deltax + layernumber * LayerHeight * deltax;
        // Call function for each layernumber, and for simulation types "S" and "R"
        int ZBound_Max_S = calcZBound_High("S", SpotRadius, LayerHeight, layernumber, ZMin, deltax, nz, ZMaxLayer_Host);
        EXPECT_EQ(ZBound_Max_S, SpotRadius + LayerHeight * layernumber);
        int ZBound_Max_R = calcZBound_High("R", SpotRadius, LayerHeight, layernumber, ZMin, deltax, nz, ZMaxLayer_Host);
        EXPECT_EQ(ZBound_Max_R, SpotRadius + LayerHeight * layernumber);
        // For simulation type C, should be independent of layernumber
        int ZBound_Max_C = calcZBound_High("C", SpotRadius, LayerHeight, layernumber, ZMin, deltax, nz, ZMaxLayer_Host);
        EXPECT_EQ(ZBound_Max_C, nz - 1);
    }
}

//---------------------------------------------------------------------------//
// temp_init_test
//---------------------------------------------------------------------------//
void testReadTemperatureData() {

    // Create test data
    double deltax = 1 * pow(10, -6);
    double HT_deltax = 1 * pow(10, -6);
    int HTtoCAratio;
    // Domain size is a 6 by 6 region
    int nx = 6;
    int ny = 6;
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

    // Test each 3 by 3 subdomain
    for (int id = 0; id < 4; id++) {
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
        ProcCol = id / 2;        // Ranks 0-1 in 0th col, 2-3 in 1st col
        MyYOffset = 3 * ProcCol; // each col is separated from the others by 3 cells
        float XMin = 0.0;
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
        ReadTemperatureData(id, deltax, HT_deltax, HTtoCAratio, MyXSlices, MyYSlices, MyXOffset, MyYOffset, XMin, YMin,
                            temp_paths, NumberOfLayers, TempFilesInSeries, NumberOfTemperatureDataPoints, RawData,
                            FirstValue, LastValue);

        // Check the results.
        // Does each rank have the right number of temperature data points? Each rank should have six (x,y,z,tm,tl,cr)
        // for each of the 9 cells in the subdomain
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
}

void testgetTempCoords() {

    // cell size and simulation lower bounds
    double deltax = 0.5;
    float XMin = -5;
    float YMin = 0.0;

    // test - each rank initializes data for layer id, of 4 total layers
    for (int layernumber = 0; layernumber < 4; layernumber++) {
        ViewF_H ZMinLayer_Host(Kokkos::ViewAllocateWithoutInitializing("ZMinLayer_Host"), 4);

        for (int n = 0; n < 4; n++) {
            ZMinLayer_Host(n) = 2 * n;
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
        int ZInt = getTempCoordZ(i, deltax, RawData, LayerHeight, layernumber, ZMinLayer_Host);
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
        ZInt = getTempCoordZ(i, deltax, RawData, LayerHeight, layernumber, ZMinLayer_Host);
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

void testInterpolateSparseData() {

    // Domain size
    int nx = 10;
    int ny = 10;
    int nz = 7;
    int HTtoCAratio = 3;

    // Views to hold data - cells start with values of -1
    ViewD3D_H TL_Host(Kokkos::ViewAllocateWithoutInitializing("TL"), nz, nx, ny);
    Kokkos::deep_copy(TL_Host, -1);
    ViewD3D_H CR_Host(Kokkos::ViewAllocateWithoutInitializing("CR"), nz, nx, ny);
    Kokkos::deep_copy(CR_Host, -1);

    // Seed TL data at coordinates of 0, 3, 6, and 9 for interpolation
    // Seeded TL data for a point (k,i,j) is equal to one more than the sum of i, j, and k

    // These values ensure that the temperature data in a volume bounded by 8 points will be interpolated
    TL_Host(0, 0, 0) = 1;
    TL_Host(3, 0, 0) = 4;
    TL_Host(0, 3, 0) = 4;
    TL_Host(0, 0, 3) = 4;
    TL_Host(3, 3, 0) = 7;
    TL_Host(3, 0, 3) = 7;
    TL_Host(0, 3, 3) = 7;
    TL_Host(3, 3, 3) = 10;

    // These values ensure that the temperature data bounded by 4 points in an XY plane will be interpolated
    TL_Host(0, 6, 0) = 7;
    TL_Host(0, 6, 3) = 10;

    // These values ensure that the temperature data bounded by 4 points in a YZ plane will be interpolated
    TL_Host(0, 3, 6) = 10;
    TL_Host(3, 3, 6) = 13;

    // This value ensures that temperature data bounded by 4 points in an XZ plane will be interpolated
    TL_Host(3, 6, 3) = 13;

    // This value ensures that the temperature data bounded by 2 points in a line in X will be interpolated
    TL_Host(0, 9, 3) = 13;

    // This value ensures that the temperature data bounded by 2 points in a line in Y will be interpolated
    TL_Host(3, 3, 9) = 16;

    // This value ensures that the temperature data bounded by 2 points in a line in Z will be interpolated
    TL_Host(6, 3, 3) = 13;

    // CR data: leave as -1 where no TL data exists, set to 1 where TL data does exist
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                if (TL_Host(k, i, j) != -1)
                    CR_Host(k, i, j) = 1;
            }
        }
    }

    // Copy data to device
    ViewD3D TL = Kokkos::create_mirror_view_and_copy(device_memory_space(), TL_Host);
    ViewD3D CR = Kokkos::create_mirror_view_and_copy(device_memory_space(), CR_Host);

    // Run interpolation
    InterpolateSparseData(TL, CR, nx, ny, nz, HTtoCAratio);

    // Expected values - start with -1s
    ViewD3D_H TL_Expected_Host("TL_Expected_Host", nz, nx, ny);
    Kokkos::deep_copy(TL_Expected_Host, -1);
    ViewD3D_H CR_Expected_Host("CR_Expected_Host", nz, nx, ny);
    Kokkos::deep_copy(CR_Expected_Host, -1);

    // These cells should have interpolated CR data (which will be all 1s, since no other CR values were given):
    // Bounded by 8 points
    for (int k = 0; k <= 3; k++) {
        for (int i = 0; i <= 3; i++) {
            for (int j = 0; j <= 3; j++) {
                CR_Expected_Host(k, i, j) = 1;
            }
        }
    }
    // Bounded by 4 points (XY plane)
    for (int i = 3; i <= 6; i++) {
        for (int j = 0; j <= 3; j++) {
            CR_Expected_Host(0, i, j) = 1;
        }
    }
    // Bounded by 4 points (YZ plane)
    for (int k = 0; k <= 3; k++) {
        for (int j = 3; j <= 6; j++) {
            CR_Expected_Host(k, 3, j) = 1;
        }
    }
    // Bounded by 4 points (XZ plane)
    for (int k = 0; k <= 3; k++) {
        for (int i = 3; i <= 6; i++) {
            CR_Expected_Host(k, i, 3) = 1;
        }
    }
    // Bounded by 2 points (X line)
    for (int i = 6; i <= 9; i++) {
        CR_Expected_Host(0, i, 3) = 1;
    }
    // Bounded by 2 points (Y line)
    for (int j = 6; j <= 9; j++) {
        CR_Expected_Host(3, 3, j) = 1;
    }
    // Bounded by 2 points (Z line)
    for (int k = 3; k <= 6; k++) {
        CR_Expected_Host(k, 3, 3) = 1;
    }

    // The same cells that have interpolated CR data should have interpolated TL data
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                if (CR_Expected_Host(k, i, j) == 1)
                    TL_Expected_Host(k, i, j) = k + i + j + 1;
            }
        }
    }

    // Copy device TL and CR back to host and check against expected values
    TL_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), TL);
    CR_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CR);
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                std::cout << "k = " << k << ", i = " << i << ", j = " << j << std::endl;
                EXPECT_DOUBLE_EQ(TL_Host(k, i, j), TL_Expected_Host(k, i, j));
                EXPECT_DOUBLE_EQ(CR_Host(k, i, j), CR_Expected_Host(k, i, j));
            }
        }
    }
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, orientation_init_tests) {
    testOrientationInit_Vectors();
    testOrientationInit_Angles();
}
TEST(TEST_CATEGORY, activedomainsizecalc) {
    // calcZBound_Low tested with and without remelting
    testcalcZBound_Low(true);
    testcalcZBound_Low(false);
    testcalcZBound_High();
}
TEST(TEST_CATEGORY, temperature_init_test) {
    testReadTemperatureData(); // reading temperature data is performed in the same manner with and without remelting
    testgetTempCoords();
    testInterpolateSparseData(); // used for no remelting cases when input data resolution and CA cell size don't match
}
} // end namespace Test
