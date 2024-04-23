// Copyright Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAgrid.hpp"
#include "CAinputs.hpp"
#include "CAprint.hpp"

#include "mpi.h"

#include <gtest/gtest.h>

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
void testPrintExaConstitDefaultRVE() {

    // Create test grid - set up so that the RVE is 5 cells in X, Y, and Z
    int number_of_layers_temp = 10;
    Grid grid(number_of_layers_temp);
    grid.nx = 10;
    grid.ny = 10;
    grid.ny_local = 10;
    grid.y_offset = 0;
    grid.nz = 10;
    grid.deltax = 0.0001; // in meters

    // default inputs struct - manually set non-default substrateInputs values
    Inputs inputs;
    inputs.print.print_default_rve = true;
    inputs.print.rve_size = Kokkos::round(0.0005 / grid.deltax);
    // File name/path for test RVE output
    inputs.print.base_filename = "TestRVE";
    // Initialize printing struct from inputs
    Print print(grid, 1, inputs.print);

    // Check that inputs in print struct match the initialization from inputs
    EXPECT_TRUE(print._inputs.print_default_rve);
    EXPECT_EQ(inputs.print.rve_size, print._inputs.rve_size);
    EXPECT_EQ(inputs.print.base_filename, print._inputs.base_filename);

    // Create test data
    Kokkos::View<int ***, Kokkos::HostSpace> grain_id_whole_domain(
        Kokkos::ViewAllocateWithoutInitializing("GrainID_WholeDomain"), grid.nz, grid.nx, grid.ny);
    Kokkos::View<short ***, Kokkos::HostSpace> layer_id_whole_domain(
        Kokkos::ViewAllocateWithoutInitializing("LayerID_WholeDomain"), grid.nz, grid.nx, grid.ny);
    for (int k = 0; k < grid.nz; k++) {
        for (int j = 0; j < grid.ny; j++) {
            for (int i = 0; i < grid.nx; i++) {
                layer_id_whole_domain(k, i, j) = k;
                grain_id_whole_domain(k, i, j) = i + j;
            }
        }
    }

    // Print RVE
    print.printExaConstitDefaultRVE(grid, layer_id_whole_domain, grain_id_whole_domain);

    // Check printed RVE
    std::ifstream grainplot_e;
    std::string expected_filename = "TestRVE_ExaConstit.csv";
    grainplot_e.open(expected_filename);
    std::string line;
    std::getline(grainplot_e, line);
    EXPECT_TRUE(line == "Coordinates are in CA units, 1 cell = 0.0001 m. Data is cell-centered. Origin at 3,3,4 , "
                        "domain size is 5 by 5 by 5 cells");
    std::getline(grainplot_e, line);
    EXPECT_TRUE(line == "X coord, Y coord, Z coord, Grain ID");
    for (int k = 4; k < 9; k++) {
        for (int j = 3; j < 8; j++) {
            for (int i = 3; i < 8; i++) {
                std::string ExpectedLine =
                    std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + "," + std::to_string(i + j);
                std::getline(grainplot_e, line);
                EXPECT_TRUE(line == ExpectedLine);
            }
        }
    }
    grainplot_e.close();
}

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, print_test) { testPrintExaConstitDefaultRVE(); }

} // end namespace Test
