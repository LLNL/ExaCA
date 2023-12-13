// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAgrid.hpp"
#include "CAinputs.hpp"
#include "CAprint.hpp"

#include <gtest/gtest.h>

#include "mpi.h"

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
// activedomainsizecalc
//---------------------------------------------------------------------------//
void testCalcZLayerBottom() {

    // Default initialized inputs and grid structs
    Inputs inputs;
    int NumberOfLayers_temp = 10;
    Grid grid(NumberOfLayers_temp);

    // Manually set layer info
    grid.LayerHeight = 10;
    grid.NumberOfLayers = 10;
    grid.deltax = 1 * pow(10, -6);
    grid.ZMin = -0.5 * pow(10, -6);
    for (int layernumber = 0; layernumber < NumberOfLayers; layernumber++) {
        // Set ZMinLayer for each layer to be offset by LayerHeight cells from the previous one (lets solution for
        // both problem types by the same)
        grid.ZMinLayer(layernumber) = grid.ZMin + layernumber * grid.LayerHeight * grid.deltax;
        // Call function for each layernumber, and for simulation types "S" and "R"
        int z_layer_bottom_S = grid.calc_z_layer_bottom("S", layernumber);
        EXPECT_EQ(z_layer_bottom_S, grid.LayerHeight * layernumber);
        int z_layer_bottom_R = grid.calc_z_layer_bottom("R", layernumber);
        EXPECT_EQ(z_layer_bottom_R, grid.LayerHeight * layernumber);
    }
}

void testCalcZLayerTop() {

    // A separate function is now used for ZBound_High calculation
    // Default initialized inputs and grid structs with manually set values for tests
    Inputs inputs;
    inputs.domain.SpotRadius = 100;
    inputs.domain.LayerHeight = 10;
    int NumberOfLayers_temp = 10;
    Grid grid(NumberOfLayers_temp);
    grid.ZMin = 0.5 * pow(10, -6);
    grid.deltax = 1.0 * pow(10, -6);
    grid.nz = 101;
    for (int layernumber = 0; layernumber < grid.NumberOfLayers; layernumber++) {
        // Set ZMaxLayer for each layer to be offset by LayerHeight cells from the previous one, with layer 0 having a
        // ZMax value of ZMin + SpotRadius (lets solution for both problem types be the same)
        grid.ZMaxLayer(layernumber) =
            grid.ZMin + inputs.domain.SpotRadius * grid.deltax + layernumber * grid.LayerHeight * grid.deltax;
        // Call function for each layernumber, and for simulation types "S" and "R"
        int z_layer_top_S = grid.calc_z_layer_top("S", inputs.domain.SpotRadius, layernumber);
        EXPECT_EQ(z_layer_top_S, inputs.domain.SpotRadius + grid.LayerHeight * layernumber);
        int z_layer_top_R = calc_z_layer_top("R", inputs.domain.SpotRadius, layernumber);
        EXPECT_EQ(z_layer_top_R, inputs.domain.SpotRadius + grid.LayerHeight * layernumber);
        // For simulation type C, should be independent of layernumber
        int z_layer_top_C = calc_z_layer_top("C", inputs.domain.SpotRadius, layernumber);
        EXPECT_EQ(z_layer_top_C, grid.nz - 1);
    }
}

void testCalcnzLayer() {

    int id = 0;
    Inputs inputs;
    int NumberOfLayers_temp = 10;
    Grid grid(NumberOfLayers_temp);
    grid.z_layer_bottom = 5;
    for (int layernumber = 0; layernumber < grid.NumberOfLayers; layernumber++) {
        grid.z_layer_top = 6 + layernumber;
        int nz_layer = grid.calc_nz_layer(id, layernumber);
        EXPECT_EQ(nz_layer, 2 + layernumber);
    }
}

void testCalcDomainSize() {

    Grid grid();
    grid.nx = 5;
    grid.ny_local = 4;
    grid.nz_layer = 10;
    int DomainSize = calcDomainSize();
    EXPECT_EQ(DomainSize, 10 * 5 * 4);
}

//---------------------------------------------------------------------------//
// bounds_init_test
//---------------------------------------------------------------------------//
void testFindXYZBounds(bool TestBinaryInputRead) {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

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

    // Set up grid
    Inputs inputs;
    inputs.temperature.TempFilesInSeries = 1;
    inputs.temperature.temp_paths.push_back(TestFilename);
    inputs.LayerHeight = 2;
    inputs.NumberOfLayers = 2;

    // Fill grid struct
    Grid grid("R", id, np, 1);

    EXPECT_DOUBLE_EQ(grid.XMin, 0.0);
    EXPECT_DOUBLE_EQ(grid.YMin, 0.0);
    EXPECT_DOUBLE_EQ(grid.ZMin, 0.0);
    EXPECT_DOUBLE_EQ(grid.XMax, 3 * grid.deltax);
    EXPECT_DOUBLE_EQ(grid.YMax, 2 * grid.deltax);
    // ZMax is equal to the largest Z coordinate in the file, offset by LayerHeight cells in the build direction due
    // to the second layer
    EXPECT_DOUBLE_EQ(grid.ZMax, 4 * grid.deltax);
    // Bounds for each individual layer - 2nd layer offset by LayerHeight cells from the first
    EXPECT_DOUBLE_EQ(grid.ZMinLayer(0), 0.0);
    EXPECT_DOUBLE_EQ(grid.ZMaxLayer(0), 2 * grid.deltax);
    EXPECT_DOUBLE_EQ(grid.ZMinLayer(1), 2 * grid.deltax);
    EXPECT_DOUBLE_EQ(grid.ZMaxLayer(1), 4 * grid.deltax);
    // Size of overall domain
    EXPECT_EQ(grid.nx, 4);
    EXPECT_EQ(grid.ny, 3);
    EXPECT_EQ(grid.nz, 5);
}

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, activedomainsizecalc) {
    testCalcZLayerBottom();
    testCalcZLayerTop();
    testCalcnzLayer();
    testCalcDomainSize();
}
TEST(TEST_CATEGORY, bounds_init_test) {
    // reading temperature files to obtain xyz bounds, using binary/non-binary format
    testFindXYZBounds(false);
    testFindXYZBounds(true);
}
} // end namespace Test
