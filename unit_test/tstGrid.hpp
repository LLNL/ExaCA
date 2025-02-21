// Copyright Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
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
    int number_of_layers_temp = 10;
    Grid grid(number_of_layers_temp);

    // Manually set layer info
    grid.layer_height = 10;
    grid.deltax = 1 * pow(10, -6);
    grid.z_min = -0.5 * pow(10, -6);

    // Test for simulation type "Spot"
    int z_layer_bottom = grid.calcZLayerBottom("Spot", 0);

    // Call function for each layernumber with simulation type and "FromFile"
    EXPECT_EQ(z_layer_bottom, 0);
    for (int layernumber = 0; layernumber < grid.number_of_layers; layernumber++) {
        // Set z_min_layer for each layer to be offset by layer_height cells from the previous one (lets solution for
        // both problem types by the same)
        grid.z_min_layer(layernumber) = grid.z_min + layernumber * grid.layer_height * grid.deltax;
        int z_layer_bottom = grid.calcZLayerBottom("FromFile", layernumber);
        EXPECT_EQ(z_layer_bottom, grid.layer_height * layernumber);
    }
}

void testCalcZLayerTop() {

    // Default initialized inputs and grid structs with manually set values for tests
    Inputs inputs;
    inputs.domain.layer_height = 10;
    int number_of_layers_temp = 10;
    Grid grid(number_of_layers_temp);
    grid.z_min = 0.5 * pow(10, -6);
    grid.deltax = 1.0 * pow(10, -6);
    for (int layernumber = 0; layernumber < number_of_layers_temp; layernumber++) {
        grid.z_max_layer(layernumber) = grid.z_min + grid.deltax * inputs.domain.layer_height;
    }
    grid.nz = 111;
    for (int layernumber = 0; layernumber < grid.number_of_layers; layernumber++) {
        // Set ZMaxLayer for each layer to be offset by layer_height cells from the previous one, with layer 0 having a
        // z_max_layer value of z_min + spot_radius (lets solution for both problem types be the same)
        // Call function for each layernumber for simulation types "FromFile"
        int z_layer_top_R = grid.calcZLayerTop("FromFile", layernumber);
        EXPECT_EQ(z_layer_top_R, (grid.z_max_layer(layernumber) - grid.z_min) / grid.deltax);
        // For simulation type Directional, should be independent of layernumber
        int z_layer_top_C = grid.calcZLayerTop("Directional", layernumber);
        EXPECT_EQ(z_layer_top_C, grid.nz - 1);
        int z_layer_top_S = grid.calcZLayerTop("Spot", layernumber);
        EXPECT_EQ(z_layer_top_S, grid.nz - 1);
    }
}

void testCalcNzLayer() {

    int id = 0;
    Inputs inputs;
    int number_of_layers_temp = 10;
    Grid grid(number_of_layers_temp);
    grid.z_layer_bottom = 5;
    for (int layernumber = 0; layernumber < grid.number_of_layers; layernumber++) {
        grid.z_layer_top = 6 + layernumber;
        int nz_layer = grid.calcNzLayer(id, layernumber);
        EXPECT_EQ(nz_layer, 2 + layernumber);
    }
}

void testCalcDomainSize() {

    Grid grid;
    grid.nx = 5;
    grid.ny_local = 4;
    grid.nz_layer = 10;
    int domain_size = grid.calcDomainSize();
    EXPECT_EQ(domain_size, 10 * 5 * 4);
}

//---------------------------------------------------------------------------//
// bounds_init_test
//---------------------------------------------------------------------------//
void testFindXYZBounds(bool test_binary_input_read, bool use_fixed_z) {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // Write fake OpenFOAM data - temperature data should be of type double
    const int nx = 4;
    const int ny = 8;
    const int nz = 3;
    double deltax = 1 * pow(10, -6);
    std::string test_filename = "TestData";
    if (test_binary_input_read)
        test_filename = test_filename + ".catemp";
    else
        test_filename = test_filename + ".txt";
    if (id == 0) {
        std::ofstream test_data_stream;
        if (test_binary_input_read)
            test_data_stream.open(test_filename, std::ios::out | std::ios::binary);
        else {
            test_data_stream.open(test_filename);
            test_data_stream << "x, y, z, tm, tl, cr" << std::endl;
        }
        // only x,y,z data should be read, tm, tl, cr should not affect result
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    if (test_binary_input_read) {
                        writeData(test_data_stream, static_cast<double>(i * deltax), test_binary_input_read);
                        writeData(test_data_stream, static_cast<double>(j * deltax), test_binary_input_read);
                        writeData(test_data_stream, static_cast<double>(k * deltax), test_binary_input_read);
                        writeData(test_data_stream, static_cast<double>(-1.0), test_binary_input_read);
                        writeData(test_data_stream, static_cast<double>(-1.0), test_binary_input_read);
                        writeData(test_data_stream, static_cast<double>(-1.0), test_binary_input_read);
                    }
                    else
                        test_data_stream << i * deltax << "," << j * deltax << "," << k * deltax << ","
                                         << static_cast<double>(-1.0) << "," << static_cast<double>(-1.0) << ","
                                         << static_cast<double>(-1.0) << std::endl;
                }
            }
        }
        test_data_stream.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Set up grid
    Inputs inputs;
    inputs.temperature.temp_files_in_series = 1;
    inputs.temperature.temp_paths.push_back(test_filename);
    inputs.domain.deltax = deltax;
    inputs.domain.layer_height = 2;
    inputs.domain.number_of_layers = 2;
    inputs.substrate.use_fixed_z_bounds = use_fixed_z;
    if (use_fixed_z) {
        // Should drop the z min value of the domain by one cell spacing, and increase the nz value by 1
        inputs.substrate.baseplate_bottom_z = -deltax;
    }

    // Fill grid struct
    Grid grid("FromFile", id, np, inputs.domain.number_of_layers, inputs.domain, inputs.substrate, inputs.temperature);

    EXPECT_DOUBLE_EQ(grid.x_min, 0.0);
    EXPECT_DOUBLE_EQ(grid.y_min, 0.0);
    if (use_fixed_z)
        EXPECT_DOUBLE_EQ(grid.z_min, -1.0 * pow(10, -6));
    else
        EXPECT_DOUBLE_EQ(grid.z_min, 0.0);
    EXPECT_DOUBLE_EQ(grid.x_max, (nx - 1) * deltax);
    EXPECT_DOUBLE_EQ(grid.y_max, (ny - 1) * deltax);
    // z_max is equal to the largest Z coordinate in the file, offset by layer_height cells in the build direction due
    // to the second layer
    EXPECT_DOUBLE_EQ(grid.z_max, 4 * grid.deltax);
    // Bounds for each individual layer - 2nd layer offset by layer_height cells from the first
    EXPECT_DOUBLE_EQ(grid.z_min_layer(0), 0.0);
    EXPECT_DOUBLE_EQ(grid.z_max_layer(0), (nz - 1) * deltax);
    EXPECT_DOUBLE_EQ(grid.z_min_layer(1), (nz - 1) * deltax);
    EXPECT_DOUBLE_EQ(grid.z_max_layer(1), (nz - 1 + inputs.domain.layer_height) * deltax);
    // Size of overall domain
    EXPECT_EQ(grid.nx, nx);
    EXPECT_EQ(grid.ny, ny);
    if (use_fixed_z)
        EXPECT_EQ(grid.nz, nz + inputs.domain.layer_height + 1);
    else
        EXPECT_EQ(grid.nz, nz + inputs.domain.layer_height);
}

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, activedomainsizecalc) {
    testCalcZLayerBottom();
    testCalcZLayerTop();
    testCalcNzLayer();
    testCalcDomainSize();
}
TEST(TEST_CATEGORY, bounds_init_test) {
    // reading temperature files to obtain xyz bounds, using binary/non-binary format and with/without the bottom Z of
    // the domain fixed via the BaseplateBottomZ input
    testFindXYZBounds(false, false);
    testFindXYZBounds(true, false);
    testFindXYZBounds(false, true);
}
} // end namespace Test
