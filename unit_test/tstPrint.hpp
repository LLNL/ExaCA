// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAprint.hpp"
#include "CAtypes.hpp"

#include <gtest/gtest.h>

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
void testPrintExaConstitDefaultRVE() {

    // Create test grid - set up so that the RVE is 5 cells in X, Y, and Z
    int nx = 10;
    int ny = 10;
    int ny_local = 10;
    int y_offset = 0;
    int nz = 10;
    int NumberOfLayers = 10;
    double deltax = 0.0001; // in meters

    // default inputs struct - manually set non-default substrateInputs values
    Inputs inputs;
    inputs.print.PrintDefaultRVE = true;
    inputs.print.RVESize = 0.0005 / deltax;
    // File name/path for test RVE output
    inputs.print.BaseFileName = "TestRVE";
    // Initialize printing struct from inputs
    Print print(nx, ny, nz, ny_local, y_offset, 1, inputs.print);

    // Check that inputs in print struct match the initialization from inputs
    EXPECT_TRUE(print._inputs.PrintDefaultRVE);
    EXPECT_DOUBLE_EQ(inputs.print.RVESize, print._inputs.RVESize);
    EXPECT_EQ(inputs.print.BaseFileName, print._inputs.BaseFileName);

    // Create test data
    ViewI3D_H GrainID_WholeDomain(Kokkos::ViewAllocateWithoutInitializing("GrainID_WholeDomain"), nz, nx, ny);
    ViewI3D_H LayerID_WholeDomain(Kokkos::ViewAllocateWithoutInitializing("LayerID_WholeDomain"), nz, nx, ny);
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                LayerID_WholeDomain(k, i, j) = k;
                GrainID_WholeDomain(k, i, j) = i + j;
            }
        }
    }

    // Print RVE
    print.printExaConstitDefaultRVE(nx, ny, nz, LayerID_WholeDomain, GrainID_WholeDomain, deltax, NumberOfLayers);

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
