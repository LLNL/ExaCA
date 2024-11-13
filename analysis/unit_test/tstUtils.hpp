// Copyright Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAprint.hpp"
#include "GAutils.hpp"

#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
// GAutils tests without Kokkos
//---------------------------------------------------------------------------//
template <typename Print3DViewType>
void writeTestMicrostructureData(std::ofstream &output_fstream, const int nx, const int ny, const int nz,
                                 const std::string data_label, Print3DViewType view_data, const bool write_binary_vtk) {
    output_fstream << "LOOKUP_TABLE default" << std::endl;
    for (int coord_z = 0; coord_z < nz; coord_z++) {
        for (int coord_y = 0; coord_y < ny; coord_y++) {
            for (int coord_x = 0; coord_x < nx; coord_x++) {
                if (data_label == "int") {
                    int writeval = static_cast<int>(view_data(coord_z, coord_x, coord_y));
                    writeData(output_fstream, writeval, write_binary_vtk, true);
                }
                else if (data_label == "short") {
                    short writeval = static_cast<short>(view_data(coord_z, coord_x, coord_y));
                    writeData(output_fstream, writeval, write_binary_vtk, true);
                }
                else if (data_label == "float") {
                    float writeval = static_cast<float>(view_data(coord_z, coord_x, coord_y));
                    writeData(output_fstream, writeval, write_binary_vtk, true);
                }
            }
        }
        // Do not insert newline character if using binary writing, as this will break the binary data read by
        // adding a blank line
        if (!(write_binary_vtk))
            output_fstream << std::endl;
    }
}

void testInitializeData(const bool write_binary_vtk, const bool read_write_layer_id,
                        const bool read_write_unused_field) {

    // Fixed domain size
    const int nx = 2;
    const int ny = 4;
    const int nz = 3;
    const int num_points = nx * ny * nz;

    // Analysis data currently stored only on host
    Kokkos::View<float ***, Kokkos::HostSpace> dummy_id(Kokkos::ViewAllocateWithoutInitializing("dummy_id"), nz, nx,
                                                        ny);
    Kokkos::View<short ***, Kokkos::HostSpace> layer_id_write(Kokkos::ViewAllocateWithoutInitializing("layer_id_w"), nz,
                                                              nx, ny);
    Kokkos::View<int ***, Kokkos::HostSpace> grain_id_write(Kokkos::ViewAllocateWithoutInitializing("grain_id_w"), nz,
                                                            nx, ny);
    for (int coord_z = 0; coord_z < nz; coord_z++) {
        for (int coord_y = 0; coord_y < ny; coord_y++) {
            for (int coord_x = 0; coord_x < nx; coord_x++) {
                dummy_id(coord_z, coord_x, coord_y) = 0.7;
                layer_id_write(coord_z, coord_x, coord_y) = static_cast<short>(coord_z);
                grain_id_write(coord_z, coord_x, coord_y) = coord_x;
            }
        }
    }
    Kokkos::View<short ***, Kokkos::HostSpace> layer_id_read(Kokkos::ViewAllocateWithoutInitializing("layer_id_r"), nz,
                                                             nx, ny);
    Kokkos::View<int ***, Kokkos::HostSpace> grain_id_read(Kokkos::ViewAllocateWithoutInitializing("grain_id_r"), nz,
                                                           nx, ny);
    Kokkos::View<short ***, Kokkos::HostSpace> phase_id_read(Kokkos::ViewAllocateWithoutInitializing("phase_id_r"), nz,
                                                             nx, ny);

    // Write microstructure header data to test file
    std::string microstructure_file = "TestMicrostructure.vtk";
    std::ofstream test_microstructure;
    test_microstructure.open(microstructure_file);
    test_microstructure << "# vtk DataFile Version 3.0" << std::endl;
    test_microstructure << "vtk output" << std::endl;
    if (write_binary_vtk)
        test_microstructure << "BINARY" << std::endl;
    else
        test_microstructure << "ASCII" << std::endl;
    test_microstructure << "DATASET STRUCTURED_POINTS" << std::endl;
    test_microstructure << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;
    test_microstructure << "ORIGIN 0 0 0" << std::endl;
    test_microstructure << "SPACING 1e-06 1e-06 1e-06" << std::endl;
    test_microstructure << "POINT_DATA " << num_points << std::endl;

    // Write microstructure field data to test file - always grain_id, optionally other fields
    test_microstructure << "SCALARS GrainID int 1" << std::endl;
    writeTestMicrostructureData(test_microstructure, nx, ny, nz, "int", grain_id_write, write_binary_vtk);
    if (read_write_unused_field) {
        test_microstructure << "SCALARS DummyID float 1" << std::endl;
        writeTestMicrostructureData(test_microstructure, nx, ny, nz, "float", dummy_id, write_binary_vtk);
    }
    if (read_write_layer_id) {
        test_microstructure << "SCALARS LayerID short 1" << std::endl;
        writeTestMicrostructureData(test_microstructure, nx, ny, nz, "short", layer_id_write, write_binary_vtk);
    }
    test_microstructure.close();

    // Check that fields can be initialized correctly
    bool found_layer_id = false;
    initializeData(microstructure_file, nx, ny, nz, grain_id_read, layer_id_read, phase_id_read, found_layer_id);
    if (read_write_layer_id)
        EXPECT_TRUE(found_layer_id);
    else
        EXPECT_FALSE(found_layer_id);
    for (int coord_z = 0; coord_z < nz; coord_z++) {
        for (int coord_y = 0; coord_y < ny; coord_y++) {
            for (int coord_x = 0; coord_x < nx; coord_x++) {
                EXPECT_EQ(grain_id_read(coord_z, coord_x, coord_y), grain_id_write(coord_z, coord_x, coord_y));
                if (read_write_layer_id) {
                    EXPECT_EQ(layer_id_read(coord_z, coord_x, coord_y), layer_id_write(coord_z, coord_x, coord_y));
                }
                // For single phase w/o phase id data printed, should be initialized to zeros
                EXPECT_EQ(phase_id_read(coord_z, coord_x, coord_y), 0.0);
            }
        }
    }
}

TEST(TEST_CATEGORY, grain_util_tests) {
    // Test init for BINARY/ASCII data and 1) a case where layer_id is present, 2) a case where layer_id is not present,
    // and 3) a case where layer_id is present as well as a field that should be ignored by the analysis parsing
    std::vector<bool> layer_id_used = {true, false, true};
    std::vector<bool> dummy_id_used = {false, false, true};
    for (int test_num = 0; test_num < 3; test_num++) {
        testInitializeData(true, layer_id_used[test_num], dummy_id_used[test_num]);
        testInitializeData(false, layer_id_used[test_num], dummy_id_used[test_num]);
    }
}
} // end namespace Test
