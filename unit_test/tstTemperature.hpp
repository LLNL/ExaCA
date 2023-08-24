// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAfunctions.hpp"
#include "CAinitialize.hpp"
#include "CAparsefiles.hpp"
#include "CAtemperature.hpp"
#include "CAtypes.hpp"

#include <gtest/gtest.h>

#include "mpi.h"

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
// tests for Temperature struct
//---------------------------------------------------------------------------//
// Tests constructing temperature object in the case of reading data from a file(s)
void testReadTemperatureData(int NumberOfLayers, bool LayerwiseTempRead, bool TestBinaryInputRead) {

    using memory_space = TEST_MEMSPACE;

    // Create test data
    double deltax = 1 * pow(10, -6);
    double HT_deltax = 1 * pow(10, -6);
    int HTtoCAratio;
    // Domain size is a 3 by 12 by 3 region
    int nx = 3;
    int ny = 12;
    int nz = 3;
    int DomainSize = nx * ny * nz;
    // Write fake OpenFOAM data - only rank 0. Temperature data should be of type double
    // Write two files, one or both of which should be read
    std::string TestTempFileName1 = "TestData1";
    std::string TestTempFileName2 = "TestData2";
    if (TestBinaryInputRead) {
        TestTempFileName1 = TestTempFileName1 + ".catemp";
        TestTempFileName2 = TestTempFileName2 + ".catemp";
    }
    else {
        TestTempFileName1 = TestTempFileName1 + ".txt";
        TestTempFileName2 = TestTempFileName2 + ".txt";
    }
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
        int ny_local = 3;
        int y_offset = 3 * id; // each col is separated from the others by 3 cells
        double YMin = 0.0;
        std::vector<std::string> temp_paths(2);
        temp_paths[0] = TestTempFileName1;
        temp_paths[1] = TestTempFileName2;
        int TempFilesInSeries = 2;

        // Read in data to "RawTemperatureData"
        Temperature<memory_space> temperature(DomainSize, NumberOfLayers);
        temperature.readTemperatureData(id, deltax, HT_deltax, HTtoCAratio, y_offset, ny_local, YMin, temp_paths,
                                        NumberOfLayers, TempFilesInSeries, LayerwiseTempRead, 0);

        // Check the results.
        // Does each rank have the right number of temperature data points? Each rank should have six (x,y,z,tm,tl,cr)
        // for each of the 9 cells in the subdomain
        // If both files were read, twice as many temperature data points per file should be present
        int NumberOfTemperatureDataPoints = temperature.RawTemperatureData.extent(0);
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
                    EXPECT_DOUBLE_EQ(
                        ExpectedValues_ThisDataPoint[nn],
                        temperature.RawTemperatureData(NumberOfCellsPerRank * 6 * layercounter + 6 * n + nn));
                }
            }
        }
    }
}

// Test unidirectional solidification problem for either directional solidification or growth of a single grain seed,
// with thermal gradient G in the domain
void testInit_UnidirectionalGradient(std::string SimulationType, double G) {

    using memory_space = TEST_MEMSPACE;

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    int nx = 2;
    int ny_local = 5;
    int nz = 6; // (Front is at Z = 0 for directional growth, single grain seed at Z = 2 for singlegrain problem)
    int DomainSize = nx * ny_local * nz;
    int coord_z_Center = floorf(static_cast<float>(nz) / 2.0);
    double initUndercooling;
    if (SimulationType == "C")
        initUndercooling = 0.0;
    else
        initUndercooling = 2.0;

    // For problems with non-zero thermal gradient, 1 K difference between each cell and its neighbor in Z
    double deltax;
    if (G == 0)
        deltax = 1 * pow(10, -6);
    else
        deltax = 1.0 / G;
    double GNorm = G * deltax;
    // Cells cool at rate of 1 K per time step
    double R = 1000000;
    double deltat = 1 * pow(10, -6);
    double RNorm = R * deltat;

    // Temperature struct
    Temperature<memory_space> temperature(DomainSize, 1);
    if (G == 0)
        temperature.initialize(R, id, deltat, DomainSize, initUndercooling);
    else
        temperature.initialize(SimulationType, G, R, id, nx, ny_local, nz, deltax, deltat, DomainSize,
                               initUndercooling);

    // Copy temperature views back to host
    auto NumberOfSolidificationEvents_Host = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), temperature.NumberOfSolidificationEvents); // Copy orientation data back to the host
    auto SolidificationEventCounter_Host = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), temperature.SolidificationEventCounter); // Copy orientation data back to the host
    auto MaxSolidificationEvents_Host = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), temperature.MaxSolidificationEvents); // Copy orientation data back to the host
    auto UndercoolingCurrent_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), temperature.UndercoolingCurrent);
    auto LayerTimeTempHistory_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), temperature.LayerTimeTempHistory);

    // Check results
    int locationOfLiquidus;
    if (SimulationType == "C")
        locationOfLiquidus = 0;
    else
        locationOfLiquidus = coord_z_Center + round(initUndercooling / GNorm);
    EXPECT_EQ(MaxSolidificationEvents_Host(0), 1);
    for (int coord_z = 0; coord_z < nz; coord_z++) {
        for (int coord_x = 0; coord_x < nx; coord_x++) {
            for (int coord_y = 0; coord_y < ny_local; coord_y++) {
                int index = get1Dindex(coord_x, coord_y, coord_z, nx, ny_local);
                // Each cell solidifies once, and counter should start at 0, associated with the zeroth layer
                // MeltTimeStep should be -1 for all cells
                // Cells cool at 1 K per time step
                EXPECT_FLOAT_EQ(LayerTimeTempHistory_Host(index, 0, 0), -1.0);
                EXPECT_FLOAT_EQ(LayerTimeTempHistory_Host(index, 0, 2), RNorm);
                EXPECT_EQ(NumberOfSolidificationEvents_Host(index), 1);
                EXPECT_EQ(SolidificationEventCounter_Host(index), 0);
                // UndercoolingCurrent should be zero for cells if a positive G is given, or initUndercooling if being
                // initialized with a uniform undercooling field (all cells initially below liquidus)
                if (G == 0) {
                    EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(index), initUndercooling);
                    EXPECT_FLOAT_EQ(LayerTimeTempHistory_Host(index, 0, 1), -1);
                }
                else {
                    int distFromLiquidus = coord_z - locationOfLiquidus;
                    if (distFromLiquidus < 0) {
                        // Undercooled cell (liquidus time step already passed, set to -1)
                        // Cell is more undercooled from initUndercooling based on the magnitude of the distance from
                        // the liquidus and the thermal gradient
                        EXPECT_FLOAT_EQ(LayerTimeTempHistory_Host(index, 0, 1), -1);
                        EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(index),
                                        initUndercooling - distFromLiquidus * (G * deltax));
                    }
                    else {
                        // Cell has not yet reached the nonzero liquidus time yet (or reaches it at time = 0), either
                        // does not have an assigned undercooling or the assigned undercooling is zero
                        EXPECT_FLOAT_EQ(LayerTimeTempHistory_Host(index, 0, 1), distFromLiquidus * GNorm / RNorm);
                        EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(index), 0.0);
                    }
                }
            }
        }
    }
}

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, temperature) {
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
    // Test for directional and single grain problems, and with and without a thermal gradient for the single grain
    // problem
    testInit_UnidirectionalGradient("C", 1000000);
    testInit_UnidirectionalGradient("SingleGrain", 0);
    testInit_UnidirectionalGradient("SingleGrain", 1000000);
}
} // end namespace Test
