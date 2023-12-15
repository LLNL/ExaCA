// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
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
        grid.NeighborRank_North = MPI_PROC_NULL;
        grid.AtNorthBoundary = true;
    }
    else {
        grid.NeighborRank_North = id + 1;
        grid.AtNorthBoundary = false;
    }
    if (id == 0) {
        grid.NeighborRank_South = MPI_PROC_NULL;
        grid.AtSouthBoundary = true;
    }
    else {
        grid.NeighborRank_South = id - 1;
        grid.AtSouthBoundary = false;
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
    if (grid.AtSouthBoundary)
        grid.ny_local--;
    if (grid.AtNorthBoundary)
        grid.ny_local--;
    if (id == 0)
        grid.y_offset = 0;
    else
        grid.y_offset = 3 * (id - 1) + 2;
    grid.DomainSize = grid.nx * grid.ny_local * grid.nz_layer;
    grid.DomainSize_AllLayers = grid.nx * grid.ny_local * grid.nz;

    // Intialize grain orientations
    std::string GrainOrientationFile = checkFileInstalled("GrainOrientationVectors.csv", id);
    int NGrainOrientations = 10000; // Number of grain orientations considered in the simulation
    ViewF GrainUnitVector(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), 9 * NGrainOrientations);
    OrientationInit(id, NGrainOrientations, GrainUnitVector, GrainOrientationFile);

    // Initialize host views - set initial GrainID values to 0, all CellType values to liquid
    CellData<memory_space> cellData(grid.DomainSize_AllLayers, inputs.substrate);

    // Subviews are the portion of the domain of interst for the test (i.e., the current layer of the problem, cells
    // located at the top 5 Z coordinates)
    auto CellType = cellData.getCellTypeSubview();
    auto GrainID = cellData.getGrainIDSubview();
    auto CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
    Kokkos::deep_copy(CellType_Host, Liquid);
    auto GrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);

    // Interface struct
    // Initial size large enough to hold all data
    int BufSizeInitialEstimate = grid.nx * grid.nz_layer;
    Interface<memory_space> interface(grid.DomainSize, BufSizeInitialEstimate);
    // Copy to host for initialization
    auto DiagonalLength_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.DiagonalLength);
    auto DOCenter_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.DOCenter);
    auto CritDiagonalLength_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.CritDiagonalLength);
    // Initialize active domain views to 0
    Kokkos::deep_copy(DiagonalLength_Host, 0.0);
    Kokkos::deep_copy(DOCenter_Host, 0.0);
    Kokkos::deep_copy(CritDiagonalLength_Host, Liquid);

    // Testing of loading of ghost nodes data, sending/receiving, unpacking, and calculations on ghost node data:
    // X = 2, Z = 1 is chosen for active cell placement on all ranks
    // Active cells will be located at Y = 1 and Y = MyYSlices-2 on each rank... these are located in the halo regions
    // and should be loaded into the send buffers
    int HaloLocations_ActiveRegion[2];
    HaloLocations_ActiveRegion[0] = 1 * grid.nx * grid.ny_local + 2 * grid.ny_local + 1;
    HaloLocations_ActiveRegion[1] = 1 * grid.nx * grid.ny_local + 2 * grid.ny_local + grid.ny_local - 2;
    // Physical cell centers in Y are different for each location and each rank
    float OctCentersY[2];
    OctCentersY[0] = grid.y_offset + 1.5;
    OctCentersY[1] = grid.y_offset + grid.ny_local - 1.5;
    for (int n = 0; n < 2; n++) {
        CellType_Host(HaloLocations_ActiveRegion[n]) = Active;
        GrainID_Host(HaloLocations_ActiveRegion[n]) = 29;
        // Initialize with an initial diagonal length, octahedra centered at cell center (X = 2.5, Y = varied, Z = 1.5)
        DiagonalLength_Host(HaloLocations_ActiveRegion[n]) = 0.01;
        DOCenter_Host(3 * HaloLocations_ActiveRegion[n]) = 2.5;
        DOCenter_Host(3 * HaloLocations_ActiveRegion[n] + 1) = OctCentersY[n];
        DOCenter_Host(3 * HaloLocations_ActiveRegion[n] + 2) = 1.5;
    }

    // Also testing an alternate situation where the ghost node data should NOT be unpacked, and the cells do not need
    // to be updated: X = 2, Z = 2 is chosen for active cell placement on all ranks Four active cells are located at Y =
    // 0, Y = 1, Y = ny_local - 2, and Y = ny_local - 1
    int HaloLocations_Alt_ActiveRegion[4];
    HaloLocations_Alt_ActiveRegion[0] = 2 * grid.nx * grid.ny_local + 2 * grid.ny_local + 0;
    HaloLocations_Alt_ActiveRegion[1] = 2 * grid.nx * grid.ny_local + 2 * grid.ny_local + 1;
    HaloLocations_Alt_ActiveRegion[2] = 2 * grid.nx * grid.ny_local + 2 * grid.ny_local + grid.ny_local - 2;
    HaloLocations_Alt_ActiveRegion[3] = 2 * grid.nx * grid.ny_local + 2 * grid.ny_local + grid.ny_local - 1;
    // Physical cell centers in Y are different for each location and each rank
    float OctCentersY_Alt[4];
    OctCentersY_Alt[0] = grid.y_offset + 0.5;
    OctCentersY_Alt[1] = grid.y_offset + 1.5;
    OctCentersY_Alt[2] = grid.y_offset + grid.ny_local - 1.5;
    OctCentersY_Alt[3] = grid.y_offset + grid.ny_local - 0.5;
    for (int n = 0; n < 4; n++) {
        CellType_Host(HaloLocations_Alt_ActiveRegion[n]) = Active;
        GrainID_Host(HaloLocations_Alt_ActiveRegion[n]) = -id;
        // Initialize with an initial diagonal length, octahedra centered at cell center (X = 2.5, Y = varied, Z = 2.5)
        DiagonalLength_Host(HaloLocations_Alt_ActiveRegion[n]) = 0.01;
        DOCenter_Host(3 * HaloLocations_Alt_ActiveRegion[n]) = 2.5;
        DOCenter_Host(3 * HaloLocations_Alt_ActiveRegion[n] + 1) = OctCentersY_Alt[n];
        DOCenter_Host(3 * HaloLocations_Alt_ActiveRegion[n] + 2) = 2.5;
    }

    // Copy view data to the device
    CellType = Kokkos::create_mirror_view_and_copy(memory_space(), CellType_Host);
    GrainID = Kokkos::create_mirror_view_and_copy(memory_space(), GrainID_Host);
    interface.DiagonalLength = Kokkos::create_mirror_view_and_copy(memory_space(), DiagonalLength_Host);
    interface.CritDiagonalLength = Kokkos::create_mirror_view_and_copy(memory_space(), CritDiagonalLength_Host);
    interface.DOCenter = Kokkos::create_mirror_view_and_copy(memory_space(), DOCenter_Host);

    // Fill send buffers
    auto BufferNorthSend_local = interface.BufferNorthSend;
    auto BufferSouthSend_local = interface.BufferSouthSend;
    auto SendSizeNorth_local = interface.SendSizeNorth;
    auto SendSizeSouth_local = interface.SendSizeSouth;
    auto BufSize_local = interface.BufSize Kokkos::parallel_for(
        "testloadghostnodes", DomainSize, KOKKOS_LAMBDA(const int &index) {
            // 3D Coordinate of this cell on the "global" (all cells in the Z direction) grid
            if (CellType(index) == Active) {
                int coord_z = grid.getCoordZ(index);
                int coord_x = grid.getCoordX(index);
                int coord_y = grid.getCoordY(index);
                int GhostGID = GrainID(index);
                float GhostDOCX = interface.DOCenter(3 * index);
                float GhostDOCY = interface.DOCenter(3 * index + 1);
                float GhostDOCZ = interface.DOCenter(3 * index + 2);
                float GhostDL = DiagonalLength(index);
                interface.loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth_local,
                                         SendSizeSouth_local, grid.ny_local, coord_x, coord_y, coord_z,
                                         grid.AtNorthBoundary, grid.AtSouthBoundary, BufferSouthSend_local,
                                         BufferNorthSend_local, NGrainOrientations, BufSize_local);
            }
        });

    interface.halo_update(0, grid, cellData, NGrainOrientations, GrainUnitVector);

    // Copy CellType, GrainID, DiagonalLength, DOCenter, CritDiagonalLength views to host to check values
    CellType = cellData.getCellTypeSubview();
    GrainID = cellData.getGrainIDSubview();
    GrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);
    CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
    DiagonalLength_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.DiagonalLength);
    CritDiagonalLength_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.CritDiagonalLength);
    DOCenter_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.DOCenter);

    // These cells should have new data based on received buffer information from neighboring ranks
    // Calculated critical diagonal lengths should match those expected based on the buffer values and the grain
    // orientation
    std::vector<float> CritDiagonalLength_Expected{1.715195, 1.914002, 1.927272, 2.291471, 1.851513, 2.346452, 2.236490,
                                                   2.902006, 2.613426, 1.576758, 1.576758, 2.248777, 2.266173, 2.266173,
                                                   2.248777, 1.527831, 1.527831, 1.715195, 1.927272, 1.914002, 1.851513,
                                                   2.291471, 2.902006, 2.613426, 2.346452, 2.236490};
    int HaloLocations_Unpacked_ActiveRegion[2];
    float OctCentersY_Unpacked[2];
    HaloLocations_Unpacked_ActiveRegion[0] = 1 * grid.nx * grid.ny_local + 2 * grid.ny_local;
    HaloLocations_Unpacked_ActiveRegion[1] = 1 * grid.nx * grid.ny_local + 2 * grid.ny_local + grid.ny_local - 1;
    OctCentersY_Unpacked[0] = grid.y_offset + 0.5;
    OctCentersY_Unpacked[1] = grid.y_offset + grid.ny_local - 0.5;
    for (int n = 0; n < 2; n++) {
        if (((n == 0) && (!(grid.AtSouthBoundary))) || ((n == 1) && (!(grid.AtNorthBoundary)))) {
            EXPECT_EQ(CellType_Host(HaloLocations_Unpacked_ActiveRegion[n]), Active);
            EXPECT_EQ(GrainID_Host(HaloLocations_Unpacked_ActiveRegion[n]), 29);
            EXPECT_FLOAT_EQ(DiagonalLength_Host(HaloLocations_Unpacked_ActiveRegion[n]), 0.01);
            EXPECT_FLOAT_EQ(DOCenter_Host(3 * HaloLocations_Unpacked_ActiveRegion[n]), 2.5);
            EXPECT_FLOAT_EQ(DOCenter_Host(3 * HaloLocations_Unpacked_ActiveRegion[n] + 1), OctCentersY_Unpacked[n]);
            EXPECT_FLOAT_EQ(DOCenter_Host(3 * HaloLocations_Unpacked_ActiveRegion[n] + 2), 1.5);
            // Were critical diagonal lengths correctly calculated from the unloaded buffer data?
            for (int l = 0; l < 26; l++) {
                EXPECT_FLOAT_EQ(CritDiagonalLength_Host(26 * HaloLocations_Unpacked_ActiveRegion[n] + l),
                                CritDiagonalLength_Expected[l]);
            }
        }
    }
    // These cells should not have been modified, as their data was unchanged
    for (int n = 0; n < 4; n++) {
        EXPECT_EQ(CellType_Host(HaloLocations_Alt_ActiveRegion[n]), Active);
        EXPECT_EQ(GrainID_Host(HaloLocations_Alt_ActiveRegion[n]), -id);
        EXPECT_FLOAT_EQ(DiagonalLength_Host(HaloLocations_Alt_ActiveRegion[n]), 0.01);
        EXPECT_FLOAT_EQ(DOCenter_Host(3 * HaloLocations_Alt_ActiveRegion[n]), 2.5);
        EXPECT_FLOAT_EQ(DOCenter_Host(3 * HaloLocations_Alt_ActiveRegion[n] + 1), OctCentersY_Alt[n]);
        EXPECT_FLOAT_EQ(DOCenter_Host(3 * HaloLocations_Alt_ActiveRegion[n] + 2), 2.5);
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
        grid.AtSouthBoundary = true;
    else
        grid.AtSouthBoundary = false;
    if (id == np - 1)
        grid.AtNorthBoundary = true;
    else
        grid.AtNorthBoundary = false;

    // 3 by 4 by 10 domain, only the top half of it is part of the current layer
    // Each rank has a portion of the domain subdivided in the Y direction
    grid.nx = 3;
    grid.ny_local = 4;
    grid.nz = 20;
    grid.z_layer_bottom = 5;
    grid.nz_layer = 10;
    grid.DomainSize = nx * ny_local * nz_layer;
    grid.DomainSize_AllLayers = nx * ny_local * nz;
    int NGrainOrientations = 10000;

    // Allocate device views: entire domain on each rank
    // Default to wall cells (CellType(index) = 0) with GrainID of 0
    CellData<memory_space> cellData(grid.DomainSize_AllLayers, inputs.substrate);
    Kokkos::deep_copy(cellData.CellType_AllLayers, Liquid);
    auto CellType = cellData.getCellTypeSubview();
    auto GrainID = cellData.getGrainIDSubview();

    // Interface struct - set buffer size to 1
    int BufSizeInitialEstimate = 1;
    Interface<memory_space> interface(grid.DomainSize, BufSizeInitialEstimate);

    // Start with 2 cells in the current layer active (one in each buffer), GrainID equal to the X coordinate, Diagonal
    // length equal to the Y coordinate, octahedron center at (x + 0.5, y + 0.5, z + 0.5)
    auto BufferNorthSend_local = interface.BufferNorthSend;
    auto BufferSouthSend_local = interface.BufferSouthSend;
    auto SendSizeNorth_local = interface.SendSizeNorth;
    auto SendSizeSouth_local = interface.SendSizeSouth;
    auto BufSize_local = interface.BufSize Kokkos::parallel_for(
        "InitDomainActiveCellsNorth", 1, KOKKOS_LAMBDA(const int &) {
            int coord_x = 1;
            int coord_y = ny_local - 2;
            int coord_z = 0;
            int index = grid.get1Dindex(coord_x, coord_y, coord_z);
            CellType(index) = Active;
            GrainID(index) = coord_x;
            int GhostGID = coord_x;
            interface.DOCenter(3 * index) = coord_x + 0.5;
            float GhostDOCX = coord_x + 0.5;
            interface.DOCenter(3 * index + 1) = coord_y + 0.5;
            float GhostDOCY = coord_y + 0.5;
            interface.DOCenter(3 * index + 2) = coord_z + 0.5;
            float GhostDOCZ = coord_z + 0.5;
            interface.DiagonalLength(index) = static_cast<float>(coord_y);
            float GhostDL = static_cast<float>(coord_y);
            // Load into appropriate buffers
            interface.loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth_local,
                                     SendSizeSouth_local, grid.ny_local, coord_x, coord_y, coord_z,
                                     grid.AtNorthBoundary, grid.AtSouthBoundary, BufferSouthSend_local,
                                     BufferNorthSend_local, NGrainOrientations, BufSize_local);
        });
    Kokkos::parallel_for(
        "InitDomainActiveCellsSouth", 1, KOKKOS_LAMBDA(const int &) {
            int coord_x = 1;
            int coord_y = 1;
            int coord_z = 0;
            int index = grid.get1Dindex(coord_x, coord_y, coord_z);
            CellType(index) = Active;
            GrainID(index) = coord_x;
            int GhostGID = coord_x;
            interface.DOCenter(3 * index) = coord_x + 0.5;
            float GhostDOCX = coord_x + 0.5;
            interface.DOCenter(3 * index + 1) = coord_y + 0.5;
            float GhostDOCY = coord_y + 0.5;
            interface.DOCenter(3 * index + 2) = coord_z + 0.5;
            float GhostDOCZ = coord_z + 0.5;
            interface.DiagonalLength(index) = static_cast<float>(coord_y);
            float GhostDL = static_cast<float>(coord_y);
            // Load into appropriate buffers
            interface.loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth_local,
                                     SendSizeSouth_local, grid.ny_local, coord_x, coord_y, coord_z,
                                     grid.AtNorthBoundary, grid.AtSouthBoundary, BufferSouthSend_local,
                                     BufferNorthSend_local, NGrainOrientations, BufSize_local);
        });

    // Each rank will have "id % 4" cells of additional data to send to the south, and 1 cell of additional data to send
    // to the north. These cells will be marked as having failed to be loaded into the buffers, as the buffer sizes are
    // only 1
    int FailedBufferCellsSouth = id % 4;
    int FailedBufferCellsNorth = 1;
    Kokkos::parallel_for(
        "InitDomainFailedNorth", FailedBufferCellsNorth, KOKKOS_LAMBDA(const int &k) {
            int coord_x = 2;
            int coord_y = ny_local - 2;
            int coord_z = k;
            int index = grid.get1Dindex(coord_x, coord_y, coord_z);
            GrainID(index) = coord_x;
            int GhostGID = coord_x;
            interface.DOCenter(3 * index) = coord_x + 0.5;
            float GhostDOCX = coord_x + 0.5;
            interface.DOCenter(3 * index + 1) = coord_x + 0.5;
            float GhostDOCY = coord_y + 0.5;
            interface.DOCenter(3 * index + 2) = coord_y + 0.5;
            float GhostDOCZ = coord_z + 0.5;
            interface.DiagonalLength(index) = static_cast<float>(coord_y);
            float GhostDL = static_cast<float>(coord_y);
            // Attempt to load into appropriate buffers
            bool DataFitsInBuffer = interface.loadghostnodes(
                GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth_local, SendSizeSouth_local,
                grid.ny_local, coord_x, coord_y, coord_z, grid.AtNorthBoundary, grid.AtSouthBoundary,
                BufferSouthSend_local, BufferNorthSend_local, NGrainOrientations, BufSize_local);
            if (!(DataFitsInBuffer)) {
                // This cell's data did not fit in the buffer with current size BufSize - mark with temporary type
                CellType(index) = ActiveFailedBufferLoad;
            }
        });

    Kokkos::parallel_for(
        "InitDomainFailedSouth", FailedBufferCellsSouth, KOKKOS_LAMBDA(const int &k) {
            int coord_x = 2;
            int coord_y = 1;
            int coord_z = k;
            int index = grid.get1Dindex(coord_x, coord_y, coord_z);
            GrainID(index) = coord_x;
            int GhostGID = coord_x;
            interface.DOCenter(3 * index) = coord_x + 0.5;
            float GhostDOCX = coord_x + 0.5;
            interface.DOCenter(3 * index + 1) = coord_y + 0.5;
            float GhostDOCY = coord_y + 0.5;
            interface.DOCenter(3 * index + 2) = coord_z + 0.5;
            float GhostDOCZ = coord_z + 0.5;
            interface.DiagonalLength(index) = static_cast<float>(coord_y);
            float GhostDL = static_cast<float>(coord_y);
            // Attempt to load into appropriate buffers
            bool DataFitsInBuffer = interface.loadghostnodes(
                GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth_local, SendSizeSouth_local,
                grid.ny_local, coord_x, coord_y, coord_z, grid.AtNorthBoundary, grid.AtSouthBoundary,
                BufferSouthSend_local, BufferNorthSend_local, NGrainOrientations, BufSize_local);
            if (!(DataFitsInBuffer)) {
                // This cell's data did not fit in the buffer with current size BufSize - mark with temporary type
                CellType(index) = ActiveFailedBufferLoad;
            }
        });

    // Attempt to resize buffers and load the remaining data
    interface.check_buffers(id, grid, cellData, NGrainOrientations);

    // If there was 1 rank, buffer size should still be 1, as no data was loaded
    // Otherwise, 25 cells should have been added to the buffer in
    // addition to the required capacity increase during the resize
    if (np == 1)
        EXPECT_EQ(interface.BufSize, 1);
    else {
        int CapacityInc = std::min(3, np - 1) + 25;
        EXPECT_EQ(interface.BufSize, 1 + CapacityInc);
    }

    // Check that the correct information was retained in the original buffer positions
    // Check that the new information was correctly entered into the buffers
    auto BufferNorthSend_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.BufferNorthSend);
    auto BufferSouthSend_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.BufferSouthSend);
    // Check buffers other than Rank 0's south and Rank np-1's north
    if (!(grid.AtSouthBoundary)) {
        // Data previously stored in buffer
        EXPECT_FLOAT_EQ(BufferSouthSend_Host(0, 0), 1.0); // RankX
        EXPECT_FLOAT_EQ(BufferSouthSend_Host(0, 1), 0.0); // RankZ
        EXPECT_FLOAT_EQ(BufferSouthSend_Host(0, 2), 1.0); // grain orientation
        EXPECT_FLOAT_EQ(BufferSouthSend_Host(0, 3), 1.0); // grain number
        EXPECT_FLOAT_EQ(BufferSouthSend_Host(0, 4), 1.5); // oct center X
        EXPECT_FLOAT_EQ(BufferSouthSend_Host(0, 5), 1.5); // oct center Y
        EXPECT_FLOAT_EQ(BufferSouthSend_Host(0, 6), 0.5); // oct center Z
        EXPECT_FLOAT_EQ(BufferSouthSend_Host(0, 7), 1.0); // diagonal length
        // Data that should've been added to expanded buffer
        EXPECT_FLOAT_EQ(BufferSouthSend_Host(1, 0), 2.0); // RankX
        EXPECT_FLOAT_EQ(BufferSouthSend_Host(1, 1), 0.0); // RankZ
        EXPECT_FLOAT_EQ(BufferSouthSend_Host(1, 2), 2.0); // grain orientation
        EXPECT_FLOAT_EQ(BufferSouthSend_Host(1, 3), 1.0); // grain number
        EXPECT_FLOAT_EQ(BufferSouthSend_Host(1, 4), 2.5); // oct center X
        EXPECT_FLOAT_EQ(BufferSouthSend_Host(1, 5), 1.5); // oct center Y
        EXPECT_FLOAT_EQ(BufferSouthSend_Host(1, 6), 0.5); // oct center Z
        EXPECT_FLOAT_EQ(BufferSouthSend_Host(1, 7), 1.0); // diagonal length
    }
    if (!(grid.AtNorthBoundary)) {
        // Data previously stored in buffer
        EXPECT_FLOAT_EQ(BufferNorthSend_Host(0, 0), 1.0);                                // RankX
        EXPECT_FLOAT_EQ(BufferNorthSend_Host(0, 1), 0.0);                                // RankZ
        EXPECT_FLOAT_EQ(BufferNorthSend_Host(0, 2), 1.0);                                // grain orientation
        EXPECT_FLOAT_EQ(BufferNorthSend_Host(0, 3), 1.0);                                // grain number
        EXPECT_FLOAT_EQ(BufferNorthSend_Host(0, 4), 1.5);                                // oct center X
        EXPECT_FLOAT_EQ(BufferNorthSend_Host(0, 5), static_cast<float>(ny_local - 1.5)); // oct center Y
        EXPECT_FLOAT_EQ(BufferNorthSend_Host(0, 6), 0.5);                                // oct center Z
        EXPECT_FLOAT_EQ(BufferNorthSend_Host(0, 7), static_cast<float>(ny_local - 2));   // diagonal length
    }
}

void testResizeBuffers() {

    using memory_space = TEST_MEMSPACE;

    // Interface struct
    int DomainSize = 50;
    int BufSizeInitialEstimate = 50;
    // Init buffers to large size
    Interface<memory_space> interface(DomainSize, 50);

    // Fill buffers with test data
    Kokkos::parallel_for(
        "InitBuffers", BufSizeInitialEstimate, KOKKOS_LAMBDA(const int &i) {
            for (int j = 0; j < interface.BufComponents; j++) {
                interface.BufferSouthSend(i, j) = i + j;
                interface.BufferNorthSend(i, j) = i + j;
                interface.BufferSouthRecv(i, j) = i;
                interface.BufferNorthRecv(i, j) = i;
            }
        });

    // Reduce size back to the default of 25
    interface.BufSize = 25;
    interface.resize_buffers();

    // Copy buffers back to host
    auto BufferNorthSend_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.BufferNorthSend);
    auto BufferSouthSend_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.BufferSouthSend);
    auto BufferNorthRecv_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.BufferNorthRecv);
    auto BufferSouthRecv_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.BufferSouthRecv);

    // Check that original values in first 25 x 8 positions were preserved
    for (int i = 0; i < interface.BufSize; i++) {
        for (int j = 0; j < interface.BufComponents; j++) {
            EXPECT_EQ(BufferSouthSend_Host(i, j), i + j);
            EXPECT_EQ(BufferNorthSend_Host(i, j), i + j);
            EXPECT_EQ(BufferSouthRecv_Host(i, j), i);
            EXPECT_EQ(BufferNorthRecv_Host(i, j), i);
        }
    }
}

//---------------------------------------------------------------------------//
// Steering vector and decentered octahedron tests
//---------------------------------------------------------------------------//
void testFillSteeringVector_Remelt() {

    // Grid struct with manually set values
    // Create views - each rank has 125 cells, 75 of which are part of the active region of the domain
    Grid grid;
    grid.nx = 5;
    grid.ny_local = 5;
    grid.y_offset = 5;
    grid.nz = 5;
    grid.z_layer_bottom = 2;
    grid.nz_layer = 3;
    grid.DomainSize = nx * ny_local * nz_layer;
    grid.DomainSize_AllLayers = nx * ny_local * nz;

    // default inputs struct
    Inputs inputs;

    // Initialize cell/temperature structures
    CellData<memory_space> cellData(DomainSize_AllLayers, DomainSize, nx, ny_local, z_layer_bottom, inputs.substrate);
    auto CellType = cellData.getCellTypeSubview();
    auto GrainID = cellData.getGrainIDSubview();
    Temperature<memory_space> temperature(grid.DomainSize, 1, inputs.temperature);

    // Fill temperature structure
    Kokkos::parallel_for(
        "TempDataInit", DomainSize, KOKKOS_LAMBDA(const int &index) {
            GrainID(index) = 1;
            CellType(index) = TempSolid;
            // Cell coordinates on this rank in X, Y, and Z (GlobalZ = relative to domain bottom)
            int coord_z = getCoordZ(index, nx, ny_local);
            int coord_y = getCoordY(index, nx, ny_local);
            // Cells at Z = 0 through Z = 2 are Solid, Z = 3 and 4 are TempSolid
            if (coord_z <= 2) {
                // Solid cells should have -1 assigned as their melt/crit time steps
                temperature.LayerTimeTempHistory(index, 0, 0) = -1;
                temperature.LayerTimeTempHistory(index, 0, 1) = -1;
                temperature.LayerTimeTempHistory(index, 0, 2) = 0;
            }
            else {
                // Cells "melt" at a time step corresponding to their Y location in the overall domain (depends on
                // MyYOffset of the rank)
                temperature.LayerTimeTempHistory(index, 0, 0) = coord_y + y_offset + 1;
                // Cells reach liquidus during cooling 2 time steps after melting
                temperature.LayerTimeTempHistory(index, 0, 1) = temperature.LayerTimeTempHistory(index, 0, 0) + 2;
            }
            temperature.LayerTimeTempHistory(index, 0, 2) = 0.2;
            temperature.NumberOfSolidificationEvents(index) = 1;
        });

    // Interface struct
    Interface interface(grid.DomainSize);

    int numcycles = 15;
    for (int cycle = 1; cycle <= numcycles; cycle++) {
        // Update cell types, local undercooling each time step, and fill the steering vector
        std::cout << cycle << std::endl;
        interface.fill_steering_vector_remelt(cycle, grid, temperature, cellData);
    }

    // Copy CellType, SteeringVector, numSteer, UndercoolingCurrent, Buffers back to host to check steering vector
    // construction results
    CellType = cellData.getCellTypeSubview();
    Kokkos::View<int *, memory_space> CellType_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
    auto SteeringVector_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.SteeringVector);
    auto numSteer_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.numSteer);
    ViewF_H UndercoolingCurrent_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), temperature.UndercoolingCurrent);
    ViewF3D_H LayerTimeTempHistory_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), temperature.LayerTimeTempHistory);

    // Check the modified CellType and UndercoolingCurrent values on the host:
    // Check that the cells corresponding to outside of the "active" portion of the domain have unchanged values
    // Check that the cells corresponding to the "active" portion of the domain have potentially changed values
    // Z = 3: Rank 0, with all cells having GrainID = 0, should have Liquid cells everywhere, with undercooling
    // depending on the Y position Z = 3, Rank > 0, with GrainID > 0, should have TempSolid cells (if not enough time
    // steps to melt), Liquid cells (if enough time steps to melt but not reach the liquidus again), or FutureActive
    // cells (if below liquidus time step). If the cell is FutureActive type, it should also have a local undercooling
    // based on the Rank ID and the time step that it reached the liquidus relative to numcycles Z = 4, all ranks:
    // should have TempSolid cells (if not enough time steps to melt) or Liquid (if enough time steps to melt). Local
    // undercooling should be based on the Rank ID and the time step that it reached the liquidus relative to numcycles
    int FutureActiveCells = 0;
    for (int index = 0; index < grid.DomainSize; index++) {
        int coord_z = grid.getCoordZ(index);
        if (coord_z <= 2)
            EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(index), 0.0);
        else {
            if (numcycles < LayerTimeTempHistory_Host(index, 0, 0)) {
                EXPECT_EQ(CellType_Host(index), TempSolid);
                EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(index), 0.0);
            }
            else if ((numcycles >= LayerTimeTempHistory_Host(index, 0, 0)) &&
                     (numcycles <= LayerTimeTempHistory_Host(index, 0, 1))) {
                EXPECT_EQ(CellType_Host(index), Liquid);
                EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(index), 0.0);
            }
            else {
                EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(index),
                                (numcycles - LayerTimeTempHistory_Host(index, 0, 1)) * 0.2);
                if (coord_z == 4)
                    EXPECT_EQ(CellType_Host(index), Liquid);
                else {
                    EXPECT_EQ(CellType_Host(index), FutureActive);
                    FutureActiveCells++;
                }
            }
        }
    }
    // Check the steering vector values on the host
    EXPECT_EQ(FutureActiveCells, numSteer_Host(0));
    for (int i = 0; i < FutureActiveCells; i++) {
        // This cell should correspond to a cell at GlobalZ = 3 (RankZ = 1), and some X and Y
        int LowerBoundCellLocation = grid.nx * grid.ny_local - 1;
        int UpperBoundCellLocation = 2 * grid.nx * grid.ny_local;
        EXPECT_GT(SteeringVector_Host(i), LowerBoundCellLocation);
        EXPECT_LT(SteeringVector_Host(i), UpperBoundCellLocation);
    }
}

void testCalcCritDiagonalLength() {
    using memory_space = TEST_MEMSPACE;
    using view_type = Kokkos::View<float *, memory_space>;

    // 2 "cells" in the domain
    int DomainSize = 2;
    // First orientation is orientation 12 (starting indexing at 0) of GrainOrientationVectors.csv
    // Second orientation is orientation 28 (starting indexing at 0)
    std::vector<float> GrainUnitVectorV{0.52877,  0.651038, -0.544565, -0.572875, 0.747163,  0.336989,
                                        0.626272, 0.133778, 0.768041,  0.736425,  -0.530983, 0.419208,
                                        0.664512, 0.683971, -0.301012, -0.126894, 0.500241,  0.856538};
    // Load unit vectors into view
    view_type GrainUnitVector(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), 9 * DomainSize);
    auto GrainUnitVector_Host = Kokkos::create_mirror_view(Kokkos::HostSpace(), GrainUnitVector);
    for (int i = 0; i < 9 * DomainSize; i++) {
        GrainUnitVector_Host(i) = GrainUnitVectorV[i];
    }
    // Copy octahedron center and grain unit vector data to the device
    GrainUnitVector = Kokkos::create_mirror_view_and_copy(memory_space(), GrainUnitVector_Host);

    // Initialize interface struct
    Interface interface(DomainSize);

    // Load octahedron centers into test view
    view_type DOCenter_Test(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * DomainSize);
    Kokkos::parallel_for(
        "TestInitCritDiagonalLength", 1, KOKKOS_LAMBDA(const int) {
            // Octahedron center for first cell (does not align with CA cell center, which is at 31.5, 3.5, 1.5)
            DOCenter_Test(0) = 30.611142;
            DOCenter_Test(1) = 3.523741;
            DOCenter_Test(2) = 0.636301;
            // Octahedron center for second cell (aligns with CA cell center at 110.5, 60.5, 0.5))
            DOCenter_Test(3) = 110.5;
            DOCenter_Test(4) = 60.5;
            DOCenter_Test(5) = 0.5;
        });

    // Expected results of calculations
    std::vector<float> CritDiagonalLength_Expected{
        2.732952,  3.403329,  2.062575, 3.708569,  1.7573346, 3.038192, 4.378946,  1.7094444, 3.2407162,
        2.646351,  1.5504522, 3.114523, 3.121724,  3.2783692, 0.177466, 3.1648562, 1.472129,  2.497145,
        3.2025092, 1.914978,  3.057610, 1.9366806, 3.639777,  3.351627, 3.3160222, 1.4418885, 1.715195,
        1.914002,  1.927272,  2.291471, 1.851513,  2.346452,  2.236490, 2.902006,  2.613426,  1.576758,
        1.576758,  2.248777,  2.266173, 2.266173,  2.248777,  1.527831, 1.527831,  1.715195,  1.927272,
        1.914002,  1.851513,  2.291471, 2.902006,  2.613426,  2.346452, 2.236490};

    // For first grain, calculate critical diagonal lenghts using test octahedron centers
    calcCritDiagonalLength(0, 31.5, 3.5, 1.5, DOCenter_Test(0), DOCenter_Test(1), DOCenter_Test(2), interface.NeighborX,
                           interface.NeighborY, interface.NeighborZ, 0, GrainUnitVector, interface.CritDiagonalLength);
    // For second grain, calculate critical diagonal lenghts
    calcCritDiagonalLength(1, 110.5, 60.5, 0.5, DOCenter_Test(3), DOCenter_Test(4), DOCenter_Test(5),
                           interface.NeighborX, interface.NeighborY, interface.NeighborZ, 1, GrainUnitVector,
                           interface.CritDiagonalLength);

    // Copy calculated critical diagonal lengths to the host to check against expected values
    auto CritDiagonalLength_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.CritDiagonalLength);

    // Check results
    for (int i = 0; i < 26 * DomainSize; i++) {
        EXPECT_FLOAT_EQ(CritDiagonalLength_Host(i), CritDiagonalLength_Expected[i]);
    }
}

void testCreateNewOctahedron() {
    using memory_space = TEST_MEMSPACE;
    using view_type = Kokkos::View<float *, memory_space>;

    // Grid struct manually setting a 18-cell domain with origin at (10, 10, 0)
    Grid grid;
    grid.nx = 3;
    grid.ny_local = 2;
    grid.nz_layer = 3;
    grid.DomainSize = grid.nx * grid.ny_local * grid.nz_layer;
    grid.y_offset = 10;

    // Create interface struct
    Interface interface(grid.DomainSize);

    // Octahedra now use the layer coordinates, not the coordinates of the multilayer domain
    for (int coord_z = 0; coord_z < grid.nz_layer; coord_z++) {
        for (int coord_x = 0; coord_x < grid.nx; coord_x++) {
            for (int coord_y = 0; coord_y < grid.ny_local; coord_y++) {
                int index = grid.get1Dindex(coord_x, coord_y, coord_z);
                interface.create_new_octahedron(index, interface.DiagonalLength, interface.DOCenter, coord_x, coord_y,
                                                y_offset, coord_z);
            }
        }
    }

    // Copy back to host and check values
    auto DiagonalLength_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.DiagonalLength);
    auto DOCenter_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.DOCenter);
    for (int coord_z = 0; coord_z < grid.nz_layer; coord_z++) {
        for (int coord_x = 0; coord_x < grid.nx; coord_x++) {
            for (int coord_y = 0; coord_y < grid.ny_local; coord_y++) {
                int index = grid.get1Dindex(coord_x, coord_y, coord_z);
                EXPECT_FLOAT_EQ(DiagonalLength_Host(index), 0.01);
                EXPECT_FLOAT_EQ(DOCenter_Host(3 * index), coord_x + 0.5);
                EXPECT_FLOAT_EQ(DOCenter_Host(3 * index + 1), coord_y + y_offset + 0.5);
                EXPECT_FLOAT_EQ(DOCenter_Host(3 * index + 2), coord_z + 0.5);
            }
        }
    }
}

// TODO: Relocate to future tstOrientation.hpp
void testConvertGrainIDForBuffer() {
    using memory_space = TEST_MEMSPACE;
    using view_type = Kokkos::View<int *, memory_space>;

    // Create a list of integer grain ID values
    // Test positive and negative values
    int NGrainIDValues = 7;
    int TotalNGrainIDValues = 2 * NGrainIDValues;
    std::vector<int> GrainIDV(TotalNGrainIDValues);
    GrainIDV[0] = 1;
    GrainIDV[1] = 2;
    GrainIDV[2] = 3;
    GrainIDV[3] = 10000;
    GrainIDV[4] = 10001;
    GrainIDV[5] = 19999;
    GrainIDV[6] = 30132129;
    for (int n = NGrainIDValues; n < TotalNGrainIDValues; n++) {
        GrainIDV[n] = -GrainIDV[n - NGrainIDValues];
    }
    view_type GrainID(Kokkos::ViewAllocateWithoutInitializing("GrainID"), TotalNGrainIDValues);
    auto GrainID_Host = Kokkos::create_mirror_view(Kokkos::HostSpace(), GrainID);
    for (int n = 0; n < TotalNGrainIDValues; n++) {
        GrainID_Host(n) = GrainIDV[n];
    }
    GrainID = Kokkos::create_mirror_view_and_copy(memory_space(), GrainID_Host);

    int NGrainOrientations = 10000;

    // Check that these were converted to their components and back correctly
    view_type GrainID_Converted(Kokkos::ViewAllocateWithoutInitializing("GrainID_Converted"), TotalNGrainIDValues);

    Kokkos::parallel_for(
        "TestInitGrainIDs", TotalNGrainIDValues, KOKKOS_LAMBDA(const int &n) {
            int MyGrainOrientation = getGrainOrientation(GrainID(n), NGrainOrientations, false);
            int MyGrainNumber = getGrainNumber(GrainID(n), NGrainOrientations);
            GrainID_Converted[n] = getGrainID(NGrainOrientations, MyGrainOrientation, MyGrainNumber);
        });
    auto GrainID_Converted_Host = Kokkos::create_mirror_view(Kokkos::HostSpace(), GrainID_Converted);
    for (int n = 0; n < TotalNGrainIDValues; n++)
        EXPECT_EQ(GrainIDV[n], GrainID_Converted_Host(n));
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
    testConvertGrainIDForBuffer();
}
} // end namespace Test
