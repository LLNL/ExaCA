// Copyright Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAcelldata.hpp"
#include "CAgrid.hpp"
#include "CAinputs.hpp"
#include "CAtemperature.hpp"
#include "CAtypes.hpp"
#include "ExaCA.hpp"

#include <gtest/gtest.h>

#include "mpi.h"

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
// communication tests
//---------------------------------------------------------------------------//
void testHaloUpdate() {

    using memory_space = TEST_MEMSPACE;

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // Initialize empty inputs/grid structs - set these manually for test
    Inputs inputs;
    Grid grid;

    // Get neighbor rank process IDs
    if (id == np - 1) {
        grid.neighbor_rank_north = MPI_PROC_NULL;
        grid.at_north_boundary = true;
    }
    else {
        grid.neighbor_rank_north = id + 1;
        grid.at_north_boundary = false;
    }
    if (id == 0) {
        grid.neighbor_rank_south = MPI_PROC_NULL;
        grid.at_south_boundary = true;
    }
    else {
        grid.neighbor_rank_south = id - 1;
        grid.at_south_boundary = false;
    }

    // Domain is a 4 by (3 * np) by 10 region
    // The top half of the domain is the current layer and the portion of it of interest of this test
    grid.nx = 4;
    grid.nz_layer = 5;
    grid.nz = 10;
    grid.z_layer_bottom = 5;
    // Domain is subdivided in Y, with ghost nodes between ranks
    // Each rank is size 3 in Y, plus ghost nodes if needed (i.e, not at problem boundary)
    // For example, if np = 4:
    // Rank 0: Y = 0, 1, 2, 3 (Y = 3 are ghost nodes)
    // Rank 1: Y = 2, 3, 4, 5, 6 (Y = 2, 6 are ghost nodes)
    // Rank 2: Y = 5, 6, 7, 8, 9 (Y = 5, 9 are ghost nodes)
    // Rank 3: Y = 8, 9, 10, 11 (Y = 8 are ghost nodes)
    grid.ny_local = 5; // assuming two sets of ghost nodes - removed if at boundary
    if (grid.at_south_boundary)
        grid.ny_local--;
    if (grid.at_north_boundary)
        grid.ny_local--;
    if (id == 0)
        grid.y_offset = 0;
    else
        grid.y_offset = 3 * (id - 1) + 2;
    grid.domain_size = grid.nx * grid.ny_local * grid.nz_layer;
    grid.domain_size_all_layers = grid.nx * grid.ny_local * grid.nz;

    // Initialize grain orientations
    std::string grain_orientation_file = checkFileInstalled("GrainOrientationVectors.csv", id);
    Orientation<memory_space> orientation(id, grain_orientation_file, false);

    // Initialize host views - set initial GrainID values to 0, all CellType values to liquid
    CellData<memory_space> celldata(grid, inputs.substrate);

    // Subviews are the portion of the domain of interest for the test (i.e., the current layer of the problem, cells
    // located at the top 5 Z coordinates)
    auto grain_id = celldata.getGrainIDSubview(grid);
    auto cell_type_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), celldata.cell_type);
    Kokkos::deep_copy(cell_type_host, Liquid);
    auto grain_id_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), grain_id);

    // Interface struct
    // Initial size large enough to hold all data
    int buf_size_initial_estimate = grid.nx * grid.nz_layer;
    Interface<memory_space> interface(id, grid.domain_size, 0.01, buf_size_initial_estimate);
    // Copy to host for initialization
    auto diagonal_length_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.diagonal_length);
    auto octahedron_center_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.octahedron_center);
    auto crit_diagonal_length_host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.crit_diagonal_length);
    // Initialize active domain views to 0
    Kokkos::deep_copy(diagonal_length_host, 0.0);
    Kokkos::deep_copy(octahedron_center_host, 0.0);
    Kokkos::deep_copy(crit_diagonal_length_host, 0.0);

    // Testing of loading of ghost nodes data, sending/receiving, unpacking, and calculations on ghost node data:
    // X = 2, Z = 1 is chosen for active cell placement on all ranks
    // Active cells will be located at Y = 1 and Y = ny_local-2 on each rank... these are located in the halo regions
    // and should be loaded into the send buffers
    int halo_locations_active_region[2];
    halo_locations_active_region[0] = 1 * grid.nx * grid.ny_local + 2 * grid.ny_local + 1;
    halo_locations_active_region[1] = 1 * grid.nx * grid.ny_local + 2 * grid.ny_local + grid.ny_local - 2;
    // Physical cell centers in Y are different for each location and each rank
    float oct_centers_y[2];
    oct_centers_y[0] = grid.y_offset + 1.5;
    oct_centers_y[1] = grid.y_offset + grid.ny_local - 1.5;
    for (int n = 0; n < 2; n++) {
        cell_type_host(halo_locations_active_region[n]) = Active;
        grain_id_host(halo_locations_active_region[n]) = 29;
        // Initialize with an initial diagonal length, octahedra centered at cell center (X = 2.5, Y = varied, Z = 1.5)
        diagonal_length_host(halo_locations_active_region[n]) = 0.01;
        octahedron_center_host(3 * halo_locations_active_region[n]) = 2.5;
        octahedron_center_host(3 * halo_locations_active_region[n] + 1) = oct_centers_y[n];
        octahedron_center_host(3 * halo_locations_active_region[n] + 2) = 1.5;
    }

    // Also testing an alternate situation where the ghost node data should NOT be unpacked, and the cells do not need
    // to be updated: X = 2, Z = 2 is chosen for active cell placement on all ranks Four active cells are located at Y =
    // 0, Y = 1, Y = ny_local - 2, and Y = ny_local - 1
    int halo_locations_alt_active_region[4];
    halo_locations_alt_active_region[0] = 2 * grid.nx * grid.ny_local + 2 * grid.ny_local + 0;
    halo_locations_alt_active_region[1] = 2 * grid.nx * grid.ny_local + 2 * grid.ny_local + 1;
    halo_locations_alt_active_region[2] = 2 * grid.nx * grid.ny_local + 2 * grid.ny_local + grid.ny_local - 2;
    halo_locations_alt_active_region[3] = 2 * grid.nx * grid.ny_local + 2 * grid.ny_local + grid.ny_local - 1;
    // Physical cell centers in Y are different for each location and each rank
    float oct_centers_y_alt[4];
    oct_centers_y_alt[0] = grid.y_offset + 0.5;
    oct_centers_y_alt[1] = grid.y_offset + 1.5;
    oct_centers_y_alt[2] = grid.y_offset + grid.ny_local - 1.5;
    oct_centers_y_alt[3] = grid.y_offset + grid.ny_local - 0.5;
    for (int n = 0; n < 4; n++) {
        cell_type_host(halo_locations_alt_active_region[n]) = Active;
        grain_id_host(halo_locations_alt_active_region[n]) = -id;
        // Initialize with an initial diagonal length, octahedra centered at cell center (X = 2.5, Y = varied, Z = 2.5)
        diagonal_length_host(halo_locations_alt_active_region[n]) = 0.01;
        octahedron_center_host(3 * halo_locations_alt_active_region[n]) = 2.5;
        octahedron_center_host(3 * halo_locations_alt_active_region[n] + 1) = oct_centers_y_alt[n];
        octahedron_center_host(3 * halo_locations_alt_active_region[n] + 2) = 2.5;
    }

    // Copy view data to the device
    celldata.cell_type = Kokkos::create_mirror_view_and_copy(memory_space(), cell_type_host);
    grain_id = Kokkos::create_mirror_view_and_copy(memory_space(), grain_id_host);
    interface.diagonal_length = Kokkos::create_mirror_view_and_copy(memory_space(), diagonal_length_host);
    interface.crit_diagonal_length = Kokkos::create_mirror_view_and_copy(memory_space(), crit_diagonal_length_host);
    interface.octahedron_center = Kokkos::create_mirror_view_and_copy(memory_space(), octahedron_center_host);

    // Fill send buffers
    Kokkos::parallel_for(
        "testloadGhostNodes", grid.domain_size, KOKKOS_LAMBDA(const int &index) {
            // 3D Coordinate of this cell on the "global" (all cells in the Z direction) grid
            if (celldata.cell_type(index) == Active) {
                int coord_z = grid.getCoordZ(index);
                int coord_x = grid.getCoordX(index);
                int coord_y = grid.getCoordY(index);
                int ghost_gid = grain_id(index);
                float ghost_docx = interface.octahedron_center(3 * index);
                float ghost_docy = interface.octahedron_center(3 * index + 1);
                float ghost_docz = interface.octahedron_center(3 * index + 2);
                float ghost_dl = interface.diagonal_length(index);
                interface.loadGhostNodes(ghost_gid, ghost_docx, ghost_docy, ghost_docz, ghost_dl, grid.ny_local,
                                         coord_x, coord_y, coord_z, grid.at_north_boundary, grid.at_south_boundary,
                                         orientation.n_grain_orientations);
            }
        });

    haloUpdate(0, 0, grid, celldata, interface, orientation);

    // Copy views to host to check values
    grain_id = celldata.getGrainIDSubview(grid);
    grain_id_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), grain_id);
    cell_type_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), celldata.cell_type);
    diagonal_length_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.diagonal_length);
    crit_diagonal_length_host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.crit_diagonal_length);
    octahedron_center_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.octahedron_center);

    // These cells should have new data based on received buffer information from neighboring ranks
    // Calculated critical diagonal lengths should match those expected based on the buffer values and the grain
    // orientation
    std::vector<float> crit_diagonal_length_host_expected{
        1.527831, 1.715195, 1.576758, 1.527831, 1.715195, 1.576758, 1.927272, 2.248777, 1.851513,
        1.927272, 2.248777, 1.851513, 1.914002, 2.266173, 2.291471, 1.914002, 2.266173, 2.291471,
        2.902006, 2.613426, 2.236490, 2.346452, 2.346452, 2.236490, 2.613426, 2.902006};

    int halo_locations_unpacked_active_region[2];
    float oct_centers_y_unpacked[2];
    halo_locations_unpacked_active_region[0] = 1 * grid.nx * grid.ny_local + 2 * grid.ny_local;
    halo_locations_unpacked_active_region[1] = 1 * grid.nx * grid.ny_local + 2 * grid.ny_local + grid.ny_local - 1;
    oct_centers_y_unpacked[0] = grid.y_offset + 0.5;
    oct_centers_y_unpacked[1] = grid.y_offset + grid.ny_local - 0.5;
    for (int n = 0; n < 2; n++) {
        if (((n == 0) && (!(grid.at_south_boundary))) || ((n == 1) && (!(grid.at_north_boundary)))) {
            EXPECT_EQ(cell_type_host(halo_locations_unpacked_active_region[n]), Active);
            EXPECT_EQ(grain_id_host(halo_locations_unpacked_active_region[n]), 29);
            EXPECT_FLOAT_EQ(diagonal_length_host(halo_locations_unpacked_active_region[n]), 0.01);
            EXPECT_FLOAT_EQ(octahedron_center_host(3 * halo_locations_unpacked_active_region[n]), 2.5);
            EXPECT_FLOAT_EQ(octahedron_center_host(3 * halo_locations_unpacked_active_region[n] + 1),
                            oct_centers_y_unpacked[n]);
            EXPECT_FLOAT_EQ(octahedron_center_host(3 * halo_locations_unpacked_active_region[n] + 2), 1.5);
            // Were critical diagonal lengths correctly calculated from the unloaded buffer data?
            for (int l = 0; l < 26; l++) {
                EXPECT_FLOAT_EQ(crit_diagonal_length_host(26 * halo_locations_unpacked_active_region[n] + l),
                                crit_diagonal_length_host_expected[l]);
            }
        }
    }
    // These cells should not have been modified, as their data was unchanged
    for (int n = 0; n < 4; n++) {
        EXPECT_EQ(cell_type_host(halo_locations_alt_active_region[n]), Active);
        EXPECT_EQ(grain_id_host(halo_locations_alt_active_region[n]), -id);
        EXPECT_FLOAT_EQ(diagonal_length_host(halo_locations_alt_active_region[n]), 0.01);
        EXPECT_FLOAT_EQ(octahedron_center_host(3 * halo_locations_alt_active_region[n]), 2.5);
        EXPECT_FLOAT_EQ(octahedron_center_host(3 * halo_locations_alt_active_region[n] + 1), oct_centers_y_alt[n]);
        EXPECT_FLOAT_EQ(octahedron_center_host(3 * halo_locations_alt_active_region[n] + 2), 2.5);
    }
}

void testResizeRefillBuffers() {
    // Tests load_ghost_nodes, refill_buffers, check_buffers, and resize_buffers during a potential buffer overload in
    // cell_capture
    using memory_space = TEST_MEMSPACE;

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // Initialize empty inputs/grid structs
    Inputs inputs;
    Grid grid;

    // Manually fill grid struct with test values
    // MPI rank locations relative to the global grid
    if (id == 0)
        grid.at_south_boundary = true;
    else
        grid.at_south_boundary = false;
    if (id == np - 1)
        grid.at_north_boundary = true;
    else
        grid.at_north_boundary = false;

    // 3 by 4 by 10 domain, only the top half of it is part of the current layer
    // Each rank has a portion of the domain subdivided in the Y direction
    grid.nx = 3;
    grid.ny_local = 4;
    grid.nz = 20;
    grid.z_layer_bottom = 5;
    grid.nz_layer = 10;
    grid.domain_size = grid.nx * grid.ny_local * grid.nz_layer;
    grid.domain_size_all_layers = grid.nx * grid.ny_local * grid.nz;
    int n_grain_orientations = 10000;

    // Allocate device views: entire domain on each rank
    // Default to wall cells (CellType(index) = 0) with GrainID of 0
    CellData<memory_space> celldata(grid, inputs.substrate);
    Kokkos::deep_copy(celldata.cell_type, Liquid);
    auto grain_id = celldata.getGrainIDSubview(grid);

    // Orientation struct
    // Initialize grain orientations
    std::string grain_orientation_file = checkFileInstalled("GrainOrientationVectors.csv", id);
    Orientation<memory_space> orientation(id, grain_orientation_file, false);

    // Interface struct - set buffer size to 1
    int buf_size_initial_estimate = 1;
    Interface<memory_space> interface(id, grid.domain_size, 0.01, buf_size_initial_estimate);

    // Start with 2 cells in the current layer active (one in each buffer), GrainID equal to the X coordinate, Diagonal
    // length equal to the Y coordinate, octahedron center at (x + 0.5, y + 0.5, z + 0.5)
    Kokkos::parallel_for(
        "InitDomainActiveCellsNorth", 1, KOKKOS_LAMBDA(const int &) {
            int coord_x = 1;
            int coord_y = grid.ny_local - 2;
            int coord_z = 0;
            int index = grid.get1DIndex(coord_x, coord_y, coord_z);
            celldata.cell_type(index) = Active;
            grain_id(index) = coord_x;
            int ghost_gid = coord_x;
            interface.octahedron_center(3 * index) = coord_x + 0.5;
            float ghost_docx = coord_x + 0.5;
            interface.octahedron_center(3 * index + 1) = coord_y + 0.5;
            float ghost_docy = coord_y + 0.5;
            interface.octahedron_center(3 * index + 2) = coord_z + 0.5;
            float ghost_docz = coord_z + 0.5;
            interface.diagonal_length(index) = static_cast<float>(coord_y);
            float ghost_dl = static_cast<float>(coord_y);
            // Load into appropriate buffers
            interface.loadGhostNodes(ghost_gid, ghost_docx, ghost_docy, ghost_docz, ghost_dl, grid.ny_local, coord_x,
                                     coord_y, coord_z, grid.at_north_boundary, grid.at_south_boundary,
                                     n_grain_orientations);
        });
    Kokkos::parallel_for(
        "InitDomainActiveCellsSouth", 1, KOKKOS_LAMBDA(const int &) {
            int coord_x = 1;
            int coord_y = 1;
            int coord_z = 0;
            int index = grid.get1DIndex(coord_x, coord_y, coord_z);
            celldata.cell_type(index) = Active;
            grain_id(index) = coord_x;
            int ghost_gid = coord_x;
            interface.octahedron_center(3 * index) = coord_x + 0.5;
            float ghost_docx = coord_x + 0.5;
            interface.octahedron_center(3 * index + 1) = coord_y + 0.5;
            float ghost_docy = coord_y + 0.5;
            interface.octahedron_center(3 * index + 2) = coord_z + 0.5;
            float ghost_docz = coord_z + 0.5;
            interface.diagonal_length(index) = static_cast<float>(coord_y);
            float ghost_dl = static_cast<float>(coord_y);
            // Load into appropriate buffers
            interface.loadGhostNodes(ghost_gid, ghost_docx, ghost_docy, ghost_docz, ghost_dl, grid.ny_local, coord_x,
                                     coord_y, coord_z, grid.at_north_boundary, grid.at_south_boundary,
                                     n_grain_orientations);
        });

    // Each rank will have "id % 4" cells of additional data to send to the south, and 1 cell of additional data to send
    // to the north. These cells will be marked as having failed to be loaded into the buffers, as the buffer sizes are
    // only 1
    int failed_buffer_cells_south = id % 4;
    int failed_buffer_cells_north = 1;
    Kokkos::parallel_for(
        "InitDomainFailedNorth", failed_buffer_cells_north, KOKKOS_LAMBDA(const int &k) {
            int coord_x = 2;
            int coord_y = grid.ny_local - 2;
            int coord_z = k;
            int index = grid.get1DIndex(coord_x, coord_y, coord_z);
            grain_id(index) = coord_x;
            int ghost_gid = coord_x;
            interface.octahedron_center(3 * index) = coord_x + 0.5;
            float ghost_docx = coord_x + 0.5;
            interface.octahedron_center(3 * index + 1) = coord_x + 0.5;
            float ghost_docy = coord_y + 0.5;
            interface.octahedron_center(3 * index + 2) = coord_y + 0.5;
            float ghost_docz = coord_z + 0.5;
            interface.diagonal_length(index) = static_cast<float>(coord_y);
            float ghost_dl = static_cast<float>(coord_y);
            // Attempt to load into appropriate buffers
            bool data_fits_in_buffer = interface.loadGhostNodes(
                ghost_gid, ghost_docx, ghost_docy, ghost_docz, ghost_dl, grid.ny_local, coord_x, coord_y, coord_z,
                grid.at_north_boundary, grid.at_south_boundary, n_grain_orientations);
            if (!(data_fits_in_buffer)) {
                // This cell's data did not fit in the buffer with current size buf_size - mark with temporary type
                celldata.cell_type(index) = ActiveFailedBufferLoad;
            }
        });

    Kokkos::parallel_for(
        "InitDomainFailedSouth", failed_buffer_cells_south, KOKKOS_LAMBDA(const int &k) {
            int coord_x = 2;
            int coord_y = 1;
            int coord_z = k;
            int index = grid.get1DIndex(coord_x, coord_y, coord_z);
            grain_id(index) = coord_x;
            int ghost_gid = coord_x;
            interface.octahedron_center(3 * index) = coord_x + 0.5;
            float ghost_docx = coord_x + 0.5;
            interface.octahedron_center(3 * index + 1) = coord_y + 0.5;
            float ghost_docy = coord_y + 0.5;
            interface.octahedron_center(3 * index + 2) = coord_z + 0.5;
            float ghost_docz = coord_z + 0.5;
            interface.diagonal_length(index) = static_cast<float>(coord_y);
            float ghost_dl = static_cast<float>(coord_y);
            // Attempt to load into appropriate buffers
            bool data_fits_in_buffer = interface.loadGhostNodes(
                ghost_gid, ghost_docx, ghost_docy, ghost_docz, ghost_dl, grid.ny_local, coord_x, coord_y, coord_z,
                grid.at_north_boundary, grid.at_south_boundary, n_grain_orientations);
            if (!(data_fits_in_buffer)) {
                // This cell's data did not fit in the buffer with current size buf_size - mark with temporary type
                celldata.cell_type(index) = ActiveFailedBufferLoad;
            }
        });

    // Attempt to resize buffers and load the remaining data
    checkBuffers(id, 0, grid, celldata, interface, orientation.n_grain_orientations);

    // If there was 1 rank, buffer size should still be 1, as no data was loaded
    // Otherwise, 25 cells should have been added to the buffer in
    // addition to the required capacity increase during the resize
    if (np == 1)
        EXPECT_EQ(interface.buf_size, 1);
    else {
        int capacity_inc = std::min(3, np - 1) + 25;
        EXPECT_EQ(interface.buf_size, 1 + capacity_inc);
    }

    // Check that the correct information was retained in the original buffer positions
    // Check that the new information was correctly entered into the buffers
    auto buffer_north_send_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.buffer_north_send);
    auto buffer_south_send_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.buffer_south_send);
    // Check buffers other than Rank 0's south and Rank np-1's north
    if (!(grid.at_south_boundary)) {
        // Data previously stored in buffer
        EXPECT_FLOAT_EQ(buffer_south_send_host(0, 0), 1.0); // RankX
        EXPECT_FLOAT_EQ(buffer_south_send_host(0, 1), 0.0); // RankZ
        EXPECT_FLOAT_EQ(buffer_south_send_host(0, 2), 1.0); // grain orientation
        EXPECT_FLOAT_EQ(buffer_south_send_host(0, 3), 1.0); // grain number
        EXPECT_FLOAT_EQ(buffer_south_send_host(0, 4), 1.5); // oct center X
        EXPECT_FLOAT_EQ(buffer_south_send_host(0, 5), 1.5); // oct center Y
        EXPECT_FLOAT_EQ(buffer_south_send_host(0, 6), 0.5); // oct center Z
        EXPECT_FLOAT_EQ(buffer_south_send_host(0, 7), 1.0); // diagonal length
        // Data that should've been added to expanded buffer
        EXPECT_FLOAT_EQ(buffer_south_send_host(1, 0), 2.0); // RankX
        EXPECT_FLOAT_EQ(buffer_south_send_host(1, 1), 0.0); // RankZ
        EXPECT_FLOAT_EQ(buffer_south_send_host(1, 2), 2.0); // grain orientation
        EXPECT_FLOAT_EQ(buffer_south_send_host(1, 3), 1.0); // grain number
        EXPECT_FLOAT_EQ(buffer_south_send_host(1, 4), 2.5); // oct center X
        EXPECT_FLOAT_EQ(buffer_south_send_host(1, 5), 1.5); // oct center Y
        EXPECT_FLOAT_EQ(buffer_south_send_host(1, 6), 0.5); // oct center Z
        EXPECT_FLOAT_EQ(buffer_south_send_host(1, 7), 1.0); // diagonal length
    }
    if (!(grid.at_north_boundary)) {
        // Data previously stored in buffer
        EXPECT_FLOAT_EQ(buffer_north_send_host(0, 0), 1.0);                                     // RankX
        EXPECT_FLOAT_EQ(buffer_north_send_host(0, 1), 0.0);                                     // RankZ
        EXPECT_FLOAT_EQ(buffer_north_send_host(0, 2), 1.0);                                     // grain orientation
        EXPECT_FLOAT_EQ(buffer_north_send_host(0, 3), 1.0);                                     // grain number
        EXPECT_FLOAT_EQ(buffer_north_send_host(0, 4), 1.5);                                     // oct center X
        EXPECT_FLOAT_EQ(buffer_north_send_host(0, 5), static_cast<float>(grid.ny_local - 1.5)); // oct center Y
        EXPECT_FLOAT_EQ(buffer_north_send_host(0, 6), 0.5);                                     // oct center Z
        EXPECT_FLOAT_EQ(buffer_north_send_host(0, 7), static_cast<float>(grid.ny_local - 2));   // diagonal length
    }
}

void testResizeBuffers() {

    using memory_space = TEST_MEMSPACE;

    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // Interface struct
    int domain_size = 50;
    int buf_size_initial_estimate = 50;
    // Init buffers to large size
    Interface<memory_space> interface(id, domain_size, 0.01, 50);

    // Fill buffers with test data
    Kokkos::parallel_for(
        "InitBuffers", buf_size_initial_estimate, KOKKOS_LAMBDA(const int &i) {
            for (int j = 0; j < interface.buf_components; j++) {
                interface.buffer_south_send(i, j) = i + j;
                interface.buffer_north_send(i, j) = i + j;
                interface.buffer_south_recv(i, j) = i;
                interface.buffer_north_recv(i, j) = i;
            }
        });

    // Reduce size back to the default of 25
    interface.buf_size = 25;
    interface.resizeBuffers(0, 0);

    // Copy buffers back to host
    auto buffer_north_send_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.buffer_north_send);
    auto buffer_south_send_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.buffer_south_send);
    auto buffer_north_recv_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.buffer_north_recv);
    auto buffer_south_recv_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.buffer_south_recv);

    // Check that original values in first 25 x 8 positions were preserved
    for (int i = 0; i < interface.buf_size; i++) {
        for (int j = 0; j < interface.buf_components; j++) {
            EXPECT_EQ(buffer_south_send_host(i, j), i + j);
            EXPECT_EQ(buffer_north_send_host(i, j), i + j);
            EXPECT_EQ(buffer_south_recv_host(i, j), i);
            EXPECT_EQ(buffer_north_recv_host(i, j), i);
        }
    }
}

//---------------------------------------------------------------------------//
// Steering vector and decentered octahedron tests
//---------------------------------------------------------------------------//
void testFillSteeringVector_Remelt() {

    using memory_space = TEST_MEMSPACE;

    // Grid struct with manually set values
    // Create views - each rank has 125 cells, 75 of which are part of the active region of the domain
    Grid grid;
    grid.nx = 5;
    grid.ny_local = 5;
    grid.y_offset = 5;
    grid.nz = 5;
    grid.z_layer_bottom = 2;
    grid.nz_layer = 3;
    grid.domain_size = grid.nx * grid.ny_local * grid.nz_layer;
    grid.domain_size_all_layers = grid.nx * grid.ny_local * grid.nz;

    // default inputs struct
    Inputs inputs;

    // Initialize cell/temperature structures
    CellData<memory_space> celldata(grid, inputs.substrate);
    auto grain_id = celldata.getGrainIDSubview(grid);
    Temperature<memory_space> temperature(grid, inputs.temperature);

    // Fill temperature structure
    Kokkos::parallel_for(
        "TempDataInit", grid.domain_size, KOKKOS_LAMBDA(const int &index) {
            grain_id(index) = 1;
            celldata.cell_type(index) = TempSolid;
            // Cell coordinates on this rank in X, Y, and Z (GlobalZ = relative to domain bottom)
            int coord_z = grid.getCoordZ(index);
            int coord_y = grid.getCoordY(index);
            // Cells at Z = 0 through Z = 2 are Solid, Z = 3 and 4 are TempSolid
            if (coord_z <= 2) {
                // Solid cells should have -1 assigned as their melt/crit time steps
                temperature.liquidus_time(index, 0, 0) = -1;
                temperature.liquidus_time(index, 0, 1) = -1;
                temperature.cooling_rate(index, 0) = 0;
            }
            else {
                // Cells "melt" at a time step corresponding to their Y location in the overall domain (depends on
                // y_offset of the rank)
                temperature.liquidus_time(index, 0, 0) = coord_y + grid.y_offset + 1;
                // Cells reach liquidus during cooling 2 time steps after melting
                temperature.liquidus_time(index, 0, 1) = temperature.liquidus_time(index, 0, 0) + 2;
            }
            temperature.cooling_rate(index, 0) = 0.2;
            temperature.number_of_solidification_events(index) = 1;
        });

    // Interface struct
    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    Interface<memory_space> interface(id, grid.domain_size, 0.01);

    int numcycles = 15;
    for (int cycle = 1; cycle <= numcycles; cycle++) {
        // Update cell types, local undercooling each time step, and fill the steering vector
        fillSteeringVector_Remelt(cycle, grid, celldata, temperature, interface);
    }

    // Copy data back to host to check steering vector construction results
    auto cell_type_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), celldata.cell_type);
    auto steering_vector_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.steering_vector);
    auto num_steer_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.num_steer);
    auto undercooling_current_host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), temperature.undercooling_current);
    auto liquidus_time_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), temperature.liquidus_time);

    // Check the modified cell_type and undercooling_current values on the host:
    // Check that the cells corresponding to outside of the "active" portion of the domain have unchanged values
    // Check that the cells corresponding to the "active" portion of the domain have potentially changed values
    // Z = 3: Rank 0, with all cells having GrainID = 0, should have Liquid cells everywhere, with undercooling
    // depending on the Y position Z = 3, Rank > 0, with GrainID > 0, should have TempSolid cells (if not enough time
    // steps to melt), Liquid cells (if enough time steps to melt but not reach the liquidus again), or FutureActive
    // cells (if below liquidus time step). If the cell is FutureActive type, it should also have a local undercooling
    // based on the Rank ID and the time step that it reached the liquidus relative to numcycles Z = 4, all ranks:
    // should have TempSolid cells (if not enough time steps to melt) or Liquid (if enough time steps to melt). Local
    // undercooling should be based on the Rank ID and the time step that it reached the liquidus relative to numcycles
    int future_active_cells = 0;
    for (int index = 0; index < grid.domain_size; index++) {
        int coord_z = grid.getCoordZ(index);
        if (coord_z <= 2)
            EXPECT_FLOAT_EQ(undercooling_current_host(index), 0.0);
        else {
            if (numcycles < liquidus_time_host(index, 0, 0)) {
                EXPECT_EQ(cell_type_host(index), TempSolid);
                EXPECT_FLOAT_EQ(undercooling_current_host(index), 0.0);
            }
            else if ((numcycles >= liquidus_time_host(index, 0, 0)) && (numcycles <= liquidus_time_host(index, 0, 1))) {
                EXPECT_EQ(cell_type_host(index), Liquid);
                EXPECT_FLOAT_EQ(undercooling_current_host(index), 0.0);
            }
            else {
                EXPECT_FLOAT_EQ(undercooling_current_host(index), (numcycles - liquidus_time_host(index, 0, 1)) * 0.2);
                if (coord_z == 4)
                    EXPECT_EQ(cell_type_host(index), Liquid);
                else {
                    EXPECT_EQ(cell_type_host(index), FutureActive);
                    future_active_cells++;
                }
            }
        }
    }
    // Check the steering vector values on the host
    EXPECT_EQ(future_active_cells, num_steer_host(0));
    for (int i = 0; i < future_active_cells; i++) {
        // This cell should correspond to a cell at GlobalZ = 3 (RankZ = 1), and some X and Y
        int lower_bound_cell_location = grid.nx * grid.ny_local - 1;
        int upper_bound_cell_location = 2 * grid.nx * grid.ny_local;
        EXPECT_GT(steering_vector_host(i), lower_bound_cell_location);
        EXPECT_LT(steering_vector_host(i), upper_bound_cell_location);
    }
}

void testCalcCritDiagonalLength() {
    using memory_space = TEST_MEMSPACE;
    using view_type = Kokkos::View<float *, memory_space>;

    // 2 "cells" in the domain
    int domain_size = 2;
    // First orientation is orientation 12 (starting indexing at 0) of GrainOrientationVectors.csv
    // Second orientation is orientation 28 (starting indexing at 0)
    std::vector<float> grain_unit_vector_v{0.52877,  0.651038, -0.544565, -0.572875, 0.747163,  0.336989,
                                           0.626272, 0.133778, 0.768041,  0.736425,  -0.530983, 0.419208,
                                           0.664512, 0.683971, -0.301012, -0.126894, 0.500241,  0.856538};
    // Load unit vectors into view
    view_type grain_unit_vector(Kokkos::ViewAllocateWithoutInitializing("grain_unit_vector"), 9 * domain_size);
    auto grain_unit_vector_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), grain_unit_vector);
    for (int i = 0; i < 9 * domain_size; i++) {
        grain_unit_vector_host(i) = grain_unit_vector_v[i];
    }
    // Copy octahedron center and grain unit vector data to the device
    grain_unit_vector = Kokkos::create_mirror_view_and_copy(memory_space(), grain_unit_vector_host);

    // Initialize interface struct
    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    Interface<memory_space> interface(id, domain_size, 0.01);

    // Load octahedron centers into test view
    view_type octahedron_center_test(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * domain_size);
    Kokkos::parallel_for(
        "TestInitCritDiagonalLength", 1, KOKKOS_LAMBDA(const int) {
            // Octahedron center for first cell (does not align with CA cell center, which is at 31.5, 3.5, 1.5)
            octahedron_center_test(0) = 30.611142;
            octahedron_center_test(1) = 3.523741;
            octahedron_center_test(2) = 0.636301;
            // Octahedron center for second cell (aligns with CA cell center at 110.5, 60.5, 0.5))
            octahedron_center_test(3) = 110.5;
            octahedron_center_test(4) = 60.5;
            octahedron_center_test(5) = 0.5;
        });

    std::vector<float> crit_diagonal_length_expected{
        3.1648562, 2.497145, 2.646351,  1.472129,  2.732952, 1.5504522, 3.2025092, 3.114523,  3.057610,
        2.062575,  0.177466, 1.7573346, 1.914978,  3.121724, 3.708569,  3.403329,  3.2783692, 1.9366806,
        3.639777,  3.351627, 4.378946,  3.3160222, 3.038192, 1.4418885, 3.2407162, 1.7094444, 1.527831,
        1.715195,  1.576758, 1.527831,  1.715195,  1.576758, 1.927272,  2.248777,  1.851513,  1.927272,
        2.248777,  1.851513, 1.914002,  2.266173,  2.291471, 1.914002,  2.266173,  2.291471,  2.902006,
        2.613426,  2.236490, 2.346452,  2.346452,  2.236490, 2.613426,  2.902006};

    // For first grain, calculate critical diagonal lengths using test octahedron centers
    interface.calcCritDiagonalLength(0, 31.5, 3.5, 1.5, octahedron_center_test(0), octahedron_center_test(1),
                                     octahedron_center_test(2), 0, grain_unit_vector);
    // For second grain, calculate critical diagonal lengths
    interface.calcCritDiagonalLength(1, 110.5, 60.5, 0.5, octahedron_center_test(3), octahedron_center_test(4),
                                     octahedron_center_test(5), 1, grain_unit_vector);

    // Copy calculated critical diagonal lengths to the host to check against expected values
    auto crit_diagonal_length_host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.crit_diagonal_length);

    // Check results
    for (int i = 0; i < 26 * domain_size; i++) {
        EXPECT_FLOAT_EQ(crit_diagonal_length_host(i), crit_diagonal_length_expected[i]);
    }
}

void testCreateNewOctahedron() {
    using memory_space = TEST_MEMSPACE;

    // Grid struct manually setting a 18-cell domain with origin at (10, 10, 0)
    Grid grid;
    grid.nx = 3;
    grid.ny_local = 2;
    grid.nz_layer = 3;
    grid.domain_size = grid.nx * grid.ny_local * grid.nz_layer;
    grid.y_offset = 10;

    // Create interface struct
    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    Interface<memory_space> interface(id, grid.domain_size, 0.01);

    // Octahedra now use the layer coordinates, not the coordinates of the multilayer domain
    for (int coord_z = 0; coord_z < grid.nz_layer; coord_z++) {
        for (int coord_x = 0; coord_x < grid.nx; coord_x++) {
            for (int coord_y = 0; coord_y < grid.ny_local; coord_y++) {
                int index = grid.get1DIndex(coord_x, coord_y, coord_z);
                interface.createNewOctahedron(index, coord_x, coord_y, grid.y_offset, coord_z);
            }
        }
    }

    // Copy back to host and check values
    auto diagonal_length_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.diagonal_length);
    auto octahedron_center_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.octahedron_center);
    for (int coord_z = 0; coord_z < grid.nz_layer; coord_z++) {
        for (int coord_x = 0; coord_x < grid.nx; coord_x++) {
            for (int coord_y = 0; coord_y < grid.ny_local; coord_y++) {
                int index = grid.get1DIndex(coord_x, coord_y, coord_z);
                EXPECT_FLOAT_EQ(diagonal_length_host(index), 0.01);
                EXPECT_FLOAT_EQ(octahedron_center_host(3 * index), coord_x + 0.5);
                EXPECT_FLOAT_EQ(octahedron_center_host(3 * index + 1), coord_y + grid.y_offset + 0.5);
                EXPECT_FLOAT_EQ(octahedron_center_host(3 * index + 2), coord_z + 0.5);
            }
        }
    }
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, communication) {
    testHaloUpdate();
    testResizeRefillBuffers();
    testResizeBuffers();
}
TEST(TEST_CATEGORY, cell_update_tests) {
    testFillSteeringVector_Remelt();
    testCalcCritDiagonalLength();
    testCreateNewOctahedron();
}
} // end namespace Test
