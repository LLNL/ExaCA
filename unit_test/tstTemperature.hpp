// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAgrid.hpp"
#include "CAparsefiles.hpp"
#include "CAprint.hpp"
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
void testReadTemperatureData(int number_of_layers, bool layerwise_temp_read, bool test_binary_input_read) {

    using memory_space = TEST_MEMSPACE;

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // Create test data
    // Default grid, manually set values
    Grid grid(number_of_layers);
    grid.ny_local = 3;
    grid.y_offset = 3 * id; // each col is separated from the others by 3 cells
    grid.y_min = 0.0;
    grid.deltax = 1 * pow(10, -6);

    // Domain size is a 3 by 12 by 3 region
    grid.nx = 3;
    grid.ny = 12;
    grid.nz = 3;
    grid.domain_size = grid.nx * grid.ny * grid.nz;
    grid.domain_size_all_layers = grid.domain_size;
    grid.layer_range = std::make_pair(0, grid.domain_size);
    // Write fake OpenFOAM data - only rank 0. Temperature data should be of type double
    // Write two files, one or both of which should be read
    std::string test_temp_filename_1 = "TestData1";
    std::string test_temp_filename_2 = "TestData2";
    if (test_binary_input_read) {
        test_temp_filename_1 = test_temp_filename_1 + ".catemp";
        test_temp_filename_2 = test_temp_filename_2 + ".catemp";
    }
    else {
        test_temp_filename_1 = test_temp_filename_1 + ".txt";
        test_temp_filename_2 = test_temp_filename_2 + ".txt";
    }
    if (id == 0) {
        std::ofstream test_data_file_1;
        if (test_binary_input_read)
            test_data_file_1.open(test_temp_filename_1, std::ios::out | std::ios::binary);
        else {
            test_data_file_1.open(test_temp_filename_1);
            test_data_file_1 << "x, y, z, tm, tl, cr" << std::endl;
        }
        for (int j = 0; j < grid.ny; j++) {
            for (int i = 0; i < grid.nx; i++) {
                if (test_binary_input_read) {
                    writeData(test_data_file_1, static_cast<double>(i * grid.deltax), test_binary_input_read);
                    writeData(test_data_file_1, static_cast<double>(j * grid.deltax), test_binary_input_read);
                    writeData(test_data_file_1, static_cast<double>(0.0), test_binary_input_read);
                    writeData(test_data_file_1, static_cast<double>(i * j), test_binary_input_read);
                    writeData(test_data_file_1, static_cast<double>(i * j + i), test_binary_input_read);
                    writeData(test_data_file_1, static_cast<double>(i * j + j), test_binary_input_read);
                }
                else
                    test_data_file_1 << i * grid.deltax << "," << j * grid.deltax << "," << 0.0 << ","
                                     << static_cast<double>(i * j) << "," << static_cast<double>(i * j + i) << ","
                                     << static_cast<double>(i * j + j) << std::endl;
            }
        }
        test_data_file_1.close();

        std::ofstream test_data_file_2;
        if (test_binary_input_read)
            test_data_file_2.open(test_temp_filename_2, std::ios::out | std::ios::binary);
        else {
            test_data_file_2.open(test_temp_filename_2);
            test_data_file_2 << "x, y, z, tm, tl, cr" << std::endl;
        }
        for (int j = 0; j < grid.ny; j++) {
            for (int i = 0; i < grid.nx; i++) {
                if (test_binary_input_read) {
                    writeData(test_data_file_2, static_cast<double>(i * grid.deltax), test_binary_input_read);
                    writeData(test_data_file_2, static_cast<double>(j * grid.deltax), test_binary_input_read);
                    writeData(test_data_file_2, static_cast<double>(grid.deltax), test_binary_input_read);
                    writeData(test_data_file_2, static_cast<double>(i * j), test_binary_input_read);
                    writeData(test_data_file_2, static_cast<double>(i * j + i), test_binary_input_read);
                    writeData(test_data_file_2, static_cast<double>(i * j + j), test_binary_input_read);
                }
                else
                    test_data_file_2 << i * grid.deltax << "," << j * grid.deltax << "," << grid.deltax << ","
                                     << static_cast<double>(i * j) << "," << static_cast<double>(i * j + i) << ","
                                     << static_cast<double>(i * j + j) << std::endl;
            }
        }
        test_data_file_2.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Test each 12 by 3 subdomain (this test uses MPI - each rank has its own subdomain)
    // Default inputs struct - manually set non-default substrateInputs values
    Inputs inputs;
    inputs.temperature.temp_paths.push_back(test_temp_filename_1);
    inputs.temperature.temp_paths.push_back(test_temp_filename_2);
    inputs.temperature.temp_files_in_series = 2;
    inputs.temperature.layerwise_temp_read = layerwise_temp_read;

    // Ensure that constructor correctly initialized the local values of inputs
    Temperature<memory_space> temperature(grid, inputs.temperature);
    if (layerwise_temp_read)
        EXPECT_TRUE(temperature._inputs.layerwise_temp_read);
    else
        EXPECT_FALSE(temperature._inputs.layerwise_temp_read);
    EXPECT_EQ(inputs.temperature.temp_files_in_series, temperature._inputs.temp_files_in_series);
    EXPECT_TRUE(temperature._inputs.temp_paths[0] == inputs.temperature.temp_paths[0]);
    EXPECT_TRUE(temperature._inputs.temp_paths[1] == inputs.temperature.temp_paths[1]);

    // Read in data to "RawTemperatureData"
    temperature.readTemperatureData(id, grid, 0);
    // Check the results.
    // Does each rank have the right number of temperature data points? Each rank should have six (x,y,z,tm,tl,cr)
    // for each of the 9 cells in the subdomain
    // If both files were read, twice as many temperature data points per file should be present
    int number_of_temperature_data_points = temperature.raw_temperature_data.extent(0);
    int num_temp_points_multiplier;
    if (layerwise_temp_read)
        num_temp_points_multiplier = 1;
    else
        num_temp_points_multiplier = std::min(number_of_layers, inputs.temperature.temp_files_in_series);
    EXPECT_EQ(number_of_temperature_data_points, 54 * num_temp_points_multiplier);
    int number_of_cells_per_rank = 9;
    // Does each rank have the right temperature data values?
    for (int layercounter = 0; layercounter < num_temp_points_multiplier; layercounter++) {
        for (int n = 0; n < number_of_cells_per_rank; n++) {
            double expected_values_this_data_point[6];
            // Location on local grid
            int ca_row = n % 3;
            int ca_col = n / 3;
            // X Coordinate
            expected_values_this_data_point[0] = ca_row * grid.deltax;
            // Y Coordinate
            expected_values_this_data_point[1] = (ca_col + 3 * id) * grid.deltax;
            // Z Coordinate
            expected_values_this_data_point[2] = grid.deltax * layercounter;
            int x_int = expected_values_this_data_point[0] / grid.deltax;
            int y_int = expected_values_this_data_point[1] / grid.deltax;
            // Melting time
            expected_values_this_data_point[3] = x_int * y_int;
            // Liquidus time
            expected_values_this_data_point[4] = x_int * y_int + x_int;
            // Cooling rate
            expected_values_this_data_point[5] = x_int * y_int + y_int;
            for (int nn = 0; nn < 6; nn++) {
                EXPECT_DOUBLE_EQ(
                    expected_values_this_data_point[nn],
                    temperature.raw_temperature_data(number_of_cells_per_rank * 6 * layercounter + 6 * n + nn));
            }
        }
    }
}

// Test unidirectional solidification problem for either directional solidification or growth of a single grain seed,
// with thermal gradient G in the domain and initial undercooling init_undercooling at the solidification interface
void testInit_UnidirectionalGradient(const std::string simulation_type, const double G,
                                     const double init_undercooling) {

    using memory_space = TEST_MEMSPACE;

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // Default grid struct with manually set values for test
    Grid grid;
    grid.nx = 2;
    grid.ny_local = 5;
    grid.nz = 6; // (Front is at Z = 0 for directional growth, single grain seed at Z = 2 for singlegrain problem)
    grid.domain_size = grid.nx * grid.ny_local * grid.nz;
    grid.domain_size_all_layers = grid.domain_size;
    grid.number_of_layers = 1;
    grid.layer_range = std::make_pair(0, grid.domain_size);
    int coord_z_center = Kokkos::floorf(static_cast<float>(grid.nz) / 2.0);

    // default inputs struct - manually set non-default substrateInputs values
    Inputs inputs;
    inputs.simulation_type = simulation_type;
    inputs.temperature.G = G;
    inputs.temperature.init_undercooling = init_undercooling;

    // For problems with non-zero thermal gradient, 1 K difference between each cell and its neighbor in Z
    if (G == 0)
        grid.deltax = 1 * pow(10, -6);
    else
        grid.deltax = 1.0 / G;
    inputs.domain.deltax = grid.deltax;
    double G_norm = G * grid.deltax;
    // Cells cool at rate of 1 K per time step
    inputs.temperature.R = 1000000;
    double deltat = 1 * pow(10, -6);
    inputs.domain.deltat = deltat;
    double R_norm = inputs.temperature.R * deltat;

    // Temperature struct
    Temperature<memory_space> temperature(grid, inputs.temperature);
    // Test constructor initialization of _inputs
    // These should've been initialized with default values
    EXPECT_FALSE(temperature._inputs.layerwise_temp_read);
    EXPECT_EQ(temperature._inputs.temp_files_in_series, 0);
    // These should have assigned values
    EXPECT_DOUBLE_EQ(temperature._inputs.R, inputs.temperature.R);
    EXPECT_DOUBLE_EQ(temperature._inputs.G, G);
    EXPECT_DOUBLE_EQ(temperature._inputs.init_undercooling, inputs.temperature.init_undercooling);
    temperature.initialize(id, simulation_type, grid, deltat);

    // Copy temperature views back to host
    auto number_of_solidification_events_host = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), temperature.number_of_solidification_events); // Copy orientation data back to the host
    auto solidification_event_counter_host = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), temperature.solidification_event_counter); // Copy orientation data back to the host
    auto max_solidification_events_host = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), temperature.max_solidification_events); // Copy orientation data back to the host
    auto undercooling_current_host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), temperature.undercooling_current);
    auto layer_time_temp_history_host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), temperature.layer_time_temp_history);

    // Check results
    int location_init_undercooling, location_liquidus_isotherm;
    if (simulation_type == "C")
        location_init_undercooling = 0;
    else
        location_init_undercooling = coord_z_center;

    // If thermal gradient is 0, liquidus isotherm does not exist - initialize to nz to avoid divide by zero error and
    // ensure all cells are initialized as undercooled (i.e., at Z coordinates less than nz)
    if (G == 0)
        location_liquidus_isotherm = grid.nz;
    else
        location_liquidus_isotherm =
            location_init_undercooling + Kokkos::round(inputs.temperature.init_undercooling / (G * grid.deltax));

    EXPECT_EQ(max_solidification_events_host(0), 1);
    for (int coord_z = 0; coord_z < grid.nz; coord_z++) {
        for (int coord_x = 0; coord_x < grid.nx; coord_x++) {
            for (int coord_y = 0; coord_y < grid.ny_local; coord_y++) {
                int index = grid.get1DIndex(coord_x, coord_y, coord_z);
                // Each cell solidifies once, and counter should start at 0, associated with the zeroth layer
                // MeltTimeStep should be -1 for all cells
                // Cells cool at 1 K per time step
                EXPECT_FLOAT_EQ(layer_time_temp_history_host(index, 0, 0), -1.0);
                EXPECT_FLOAT_EQ(layer_time_temp_history_host(index, 0, 2), R_norm);
                EXPECT_EQ(number_of_solidification_events_host(index), 1);
                EXPECT_EQ(solidification_event_counter_host(index), 0);
                // undercooling_current should be zero for cells if a positive G is given, or init_undercooling if being
                // initialized with a uniform undercooling field (all cells initially below liquidus)
                if (G == 0) {
                    EXPECT_FLOAT_EQ(undercooling_current_host(index), inputs.temperature.init_undercooling);
                    EXPECT_FLOAT_EQ(layer_time_temp_history_host(index, 0, 1), -1);
                }
                else {
                    int dist_from_liquidus = coord_z - location_liquidus_isotherm;
                    if (dist_from_liquidus < 0) {
                        // Undercooled cell (liquidus time step already passed, set to -1)
                        // The undercooling at Z = locationOfLiquidus should have been set to init_undercooling
                        EXPECT_FLOAT_EQ(layer_time_temp_history_host(index, 0, 1), -1);
                        int dist_from_init_undercooling = coord_z - location_init_undercooling;
                        EXPECT_FLOAT_EQ(undercooling_current_host(index),
                                        inputs.temperature.init_undercooling -
                                            dist_from_init_undercooling * G * grid.deltax);
                    }
                    else {
                        // Cell has not yet reached the nonzero liquidus time yet (or reaches it at time = 0), either
                        // does not have an assigned undercooling or the assigned undercooling is zero
                        EXPECT_FLOAT_EQ(layer_time_temp_history_host(index, 0, 1),
                                        dist_from_liquidus * G_norm / R_norm);
                        EXPECT_FLOAT_EQ(undercooling_current_host(index), 0.0);
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
    // Multiple permutations of inputs: number_of_layers, layerwise_temp_read, test_binary_input_read
    // reading temperature data is performed in the same manner with and without remelting
    std::vector<int> number_of_layers_vals = {1, 2, 2, 2};
    std::vector<bool> layerwise_temp_read_vals = {false, true, false, true};
    std::vector<bool> test_binary_input_read_vals = {false, false, false, true};
    int num_vals = test_binary_input_read_vals.size();
    for (int test_count = 0; test_count < num_vals; test_count++) {
        testReadTemperatureData(number_of_layers_vals[test_count], layerwise_temp_read_vals[test_count],
                                test_binary_input_read_vals[test_count]);
    }
    // Test for directional and single grain problems, and with and without a thermal gradient for the single grain
    // problem and with/without an initial undercooling present at the initial solidification front
    testInit_UnidirectionalGradient("C", 1000000, 0);
    testInit_UnidirectionalGradient("C", 1000000, 2);
    testInit_UnidirectionalGradient("SingleGrain", 0, 2);
    testInit_UnidirectionalGradient("SingleGrain", 1000000, 2);
}
} // end namespace Test
