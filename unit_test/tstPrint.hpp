// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAgrid.hpp"
#include "CAprint.hpp"
#include "CAtypes.hpp"

#include "mpi.h"

#include <gtest/gtest.h>

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
void testPrintExaConstitDefaultRVE() {

    // Create test grid - set up so that the RVE is 5 cells in X, Y, and Z
    int NumberOfLayers_temp = 10;
    Grid grid(NumberOfLayers_temp);
    grid.nx = 10;
    grid.ny = 10;
    grid.ny_local = 10;
    grid.y_offset = 0;
    grid.nz = 10;
    grid.deltax = 0.0001; // in meters

    // default inputs struct - manually set non-default substrateInputs values
    Inputs inputs;
    inputs.print.PrintDefaultRVE = true;
    inputs.print.RVESize = 0.0005 / grid.deltax;
    // File name/path for test RVE output
    inputs.print.BaseFileName = "TestRVE";
    // Initialize printing struct from inputs
    Print print(grid, 1, inputs.print);

    // Check that inputs in print struct match the initialization from inputs
    EXPECT_TRUE(print._inputs.PrintDefaultRVE);
    EXPECT_DOUBLE_EQ(inputs.print.RVESize, print._inputs.RVESize);
    EXPECT_EQ(inputs.print.BaseFileName, print._inputs.BaseFileName);

    // Create test data
    ViewI3D_H GrainID_WholeDomain(Kokkos::ViewAllocateWithoutInitializing("GrainID_WholeDomain"), grid.nz, grid.nx,
                                  grid.ny);
    ViewS3D_H LayerID_WholeDomain(Kokkos::ViewAllocateWithoutInitializing("LayerID_WholeDomain"), grid.nz, grid.nx,
                                  grid.ny);
    for (int k = 0; k < grid.nz; k++) {
        for (int j = 0; j < grid.ny; j++) {
            for (int i = 0; i < grid.nx; i++) {
                LayerID_WholeDomain(k, i, j) = k;
                GrainID_WholeDomain(k, i, j) = i + j;
            }
        }
    }

    // Print RVE
    print.printExaConstitDefaultRVE(grid, LayerID_WholeDomain, GrainID_WholeDomain);

    // Check printed RVE
    std::ifstream GrainplotE;
    std::string ExpectedFilename = "TestRVE_ExaConstit.csv";
    GrainplotE.open(ExpectedFilename);
    std::string line;
    std::getline(GrainplotE, line);
    EXPECT_TRUE(line == "Coordinates are in CA units, 1 cell = 0.0001 m. Data is cell-centered. Origin at 3,3,4 , "
                        "domain size is 5 by 5 by 5 cells");
    std::getline(GrainplotE, line);
    EXPECT_TRUE(line == "X coord, Y coord, Z coord, Grain ID");
    for (int k = 4; k < 9; k++) {
        for (int j = 3; j < 8; j++) {
            for (int i = 3; i < 8; i++) {
                std::string ExpectedLine =
                    std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + "," + std::to_string(i + j);
                std::getline(GrainplotE, line);
                EXPECT_TRUE(line == ExpectedLine);
            }
        }
    }
    GrainplotE.close();
}

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, print_test) { testPrintExaConstitDefaultRVE(); }

} // end namespace Test
