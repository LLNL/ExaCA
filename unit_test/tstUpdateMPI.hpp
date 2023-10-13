// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>
#include <nlohmann/json.hpp>

#include "CAcelldata.hpp"
#include "CAghostnodes.hpp"
#include "CAinitialize.hpp"
#include "CAinputs.hpp"
#include "CAparsefiles.hpp"
#include "CAtypes.hpp"
#include "ExaCA.hpp"

#include <gtest/gtest.h>

#include "mpi.h"

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
// grain_init_tests
//---------------------------------------------------------------------------//
void testGhostNodes1D() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Get neighbor rank process IDs
    int NeighborRank_North, NeighborRank_South;
    bool AtNorthBoundary, AtSouthBoundary;
    if (id == np - 1) {
        NeighborRank_North = MPI_PROC_NULL;
        AtNorthBoundary = true;
    }
    else {
        NeighborRank_North = id + 1;
        AtNorthBoundary = false;
    }
    if (id == 0) {
        NeighborRank_South = MPI_PROC_NULL;
        AtSouthBoundary = true;
    }
    else {
        NeighborRank_South = id - 1;
        AtSouthBoundary = false;
    }

    // Domain is a 4 by (3 * np) by 10 region
    // The top half of the domain is the current layer and the portion of it of interest of this test
    int nx = 4;
    int nz_layer = 5;
    int nz = 10;
    int z_layer_bottom = 5;
    // Domain is subdivided in Y, with ghost nodes between ranks
    // Each rank is size 3 in Y, plus ghost nodes if needed (i.e, not at problem boundary)
    // For example, if np = 4:
    // Rank 0: Y = 0, 1, 2, 3 (Y = 3 are ghost nodes)
    // Rank 1: Y = 2, 3, 4, 5, 6 (Y = 2, 6 are ghost nodes)
    // Rank 2: Y = 5, 6, 7, 8, 9 (Y = 5, 9 are ghost nodes)
    // Rank 3: Y = 8, 9, 10, 11 (Y = 8 are ghost nodes)
    int ny_local = 5; // assuming two sets of ghost nodes - removed if at boundary
    if (AtSouthBoundary)
        ny_local--;
    if (AtNorthBoundary)
        ny_local--;
    int y_offset;
    if (id == 0)
        y_offset = 0;
    else
        y_offset = 3 * (id - 1) + 2;
    int DomainSize = nx * ny_local * nz_layer;
    int DomainSize_AllLayers = nx * ny_local * nz;

    // Intialize grain orientations
    std::string GrainOrientationFile = checkFileInstalled("GrainOrientationVectors.csv", id);
    int NGrainOrientations = 10000; // Number of grain orientations considered in the simulation
    ViewF GrainUnitVector(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), 9 * NGrainOrientations);
    OrientationInit(id, NGrainOrientations, GrainUnitVector, GrainOrientationFile);

    // Initialize neighbor lists
    NList NeighborX, NeighborY, NeighborZ;
    NeighborListInit(NeighborX, NeighborY, NeighborZ);

    // Initialize empty inputs struct
    Inputs inputs;

    // Initialize host views - set initial GrainID values to 0, all CellType values to liquid
    CellData<device_memory_space> cellData(DomainSize_AllLayers, DomainSize, nx, ny_local, z_layer_bottom,
                                           inputs.substrate);
    // Subviews are the portion of the domain of interst for the test (i.e., the current layer of the problem, cells
    // located at the top 5 Z coordinates)
    auto CellType = cellData.getCellTypeSubview();
    auto GrainID = cellData.getGrainIDSubview();
    auto CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
    Kokkos::deep_copy(CellType_Host, Liquid);
    auto GrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);
    // Initialize active domain views to 0
    ViewF_H DiagonalLength_Host("DiagonalLength_Host", DomainSize);
    ViewF_H DOCenter_Host("DOCenter_Host", 3 * DomainSize);
    ViewF_H CritDiagonalLength_Host("CritDiagonalLength_Host", 26 * DomainSize);

    // Testing of loading of ghost nodes data, sending/receiving, unpacking, and calculations on ghost node data:
    // X = 2, Z = 1 is chosen for active cell placement on all ranks
    // Active cells will be located at Y = 1 and Y = MyYSlices-2 on each rank... these are located in the halo regions
    // and should be loaded into the send buffers
    int HaloLocations_ActiveRegion[2];
    HaloLocations_ActiveRegion[0] = 1 * nx * ny_local + 2 * ny_local + 1;
    HaloLocations_ActiveRegion[1] = 1 * nx * ny_local + 2 * ny_local + ny_local - 2;
    // Physical cell centers in Y are different for each location and each rank
    float OctCentersY[2];
    OctCentersY[0] = y_offset + 1.5;
    OctCentersY[1] = y_offset + ny_local - 1.5;
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
    HaloLocations_Alt_ActiveRegion[0] = 2 * nx * ny_local + 2 * ny_local + 0;
    HaloLocations_Alt_ActiveRegion[1] = 2 * nx * ny_local + 2 * ny_local + 1;
    HaloLocations_Alt_ActiveRegion[2] = 2 * nx * ny_local + 2 * ny_local + ny_local - 2;
    HaloLocations_Alt_ActiveRegion[3] = 2 * nx * ny_local + 2 * ny_local + ny_local - 1;
    // Physical cell centers in Y are different for each location and each rank
    float OctCentersY_Alt[4];
    OctCentersY_Alt[0] = y_offset + 0.5;
    OctCentersY_Alt[1] = y_offset + 1.5;
    OctCentersY_Alt[2] = y_offset + ny_local - 1.5;
    OctCentersY_Alt[3] = y_offset + ny_local - 0.5;
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
    CellType = Kokkos::create_mirror_view_and_copy(device_memory_space(), CellType_Host);
    GrainID = Kokkos::create_mirror_view_and_copy(device_memory_space(), GrainID_Host);
    ViewF DiagonalLength = Kokkos::create_mirror_view_and_copy(device_memory_space(), DiagonalLength_Host);
    ViewF CritDiagonalLength = Kokkos::create_mirror_view_and_copy(device_memory_space(), CritDiagonalLength_Host);
    ViewF DOCenter = Kokkos::create_mirror_view_and_copy(device_memory_space(), DOCenter_Host);

    // Buffer size large enough to hold all data
    int BufSize = nx * nz_layer;

    // Send/recv buffers for ghost node data should be initialized with -1s in the first index as placeholders for empty
    // positions in the buffer
    Buffer2D BufferSouthSend("BufferSouthSend", BufSize, 8);
    Buffer2D BufferNorthSend("BufferNorthSend", BufSize, 8);
    Kokkos::parallel_for(
        "BufferReset", BufSize, KOKKOS_LAMBDA(const int &i) {
            BufferNorthSend(i, 0) = -1.0;
            BufferSouthSend(i, 0) = -1.0;
        });
    Buffer2D BufferSouthRecv(Kokkos::ViewAllocateWithoutInitializing("BufferSouthRecv"), BufSize, 8);
    Buffer2D BufferNorthRecv(Kokkos::ViewAllocateWithoutInitializing("BufferNorthRecv"), BufSize, 8);
    // Init to zero
    ViewI SendSizeSouth("SendSizeSouth", 1);
    ViewI SendSizeNorth("SendSizeNorth", 1);

    // Fill send buffers
    Kokkos::parallel_for(
        "testloadghostnodes", DomainSize, KOKKOS_LAMBDA(const int &index) {
            // 3D Coordinate of this cell on the "global" (all cells in the Z direction) grid
            if (CellType(index) == Active) {
                int coord_z = getCoordZ(index, nx, ny_local);
                int coord_x = getCoordX(index, nx, ny_local);
                int coord_y = getCoordY(index, nx, ny_local);
                int GhostGID = GrainID(index);
                float GhostDOCX = DOCenter(3 * index);
                float GhostDOCY = DOCenter(3 * index + 1);
                float GhostDOCZ = DOCenter(3 * index + 2);
                float GhostDL = DiagonalLength(index);
                loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth, SendSizeSouth,
                               ny_local, coord_x, coord_y, coord_z, AtNorthBoundary, AtSouthBoundary, BufferSouthSend,
                               BufferNorthSend, NGrainOrientations, BufSize);
            }
        });

    GhostNodes1D(0, id, NeighborRank_North, NeighborRank_South, nx, ny_local, y_offset, NeighborX, NeighborY, NeighborZ,
                 cellData, DOCenter, GrainUnitVector, DiagonalLength, CritDiagonalLength, NGrainOrientations,
                 BufferNorthSend, BufferSouthSend, BufferNorthRecv, BufferSouthRecv, BufSize, SendSizeNorth,
                 SendSizeSouth);

    // Copy CellType, GrainID, DiagonalLength, DOCenter, CritDiagonalLength views to host to check values
    CellType = cellData.getCellTypeSubview();
    GrainID = cellData.getGrainIDSubview();
    GrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);
    CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
    DiagonalLength_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), DiagonalLength);
    CritDiagonalLength_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CritDiagonalLength);
    DOCenter_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), DOCenter);

    // These cells should have new data based on received buffer information from neighboring ranks
    // Calculated critical diagonal lengths should match those expected based on the buffer values and the grain
    // orientation
    std::vector<float> CritDiagonalLength_Expected{1.715195, 1.914002, 1.927272, 2.291471, 1.851513, 2.346452, 2.236490,
                                                   2.902006, 2.613426, 1.576758, 1.576758, 2.248777, 2.266173, 2.266173,
                                                   2.248777, 1.527831, 1.527831, 1.715195, 1.927272, 1.914002, 1.851513,
                                                   2.291471, 2.902006, 2.613426, 2.346452, 2.236490};
    int HaloLocations_Unpacked_ActiveRegion[2];
    float OctCentersY_Unpacked[2];
    HaloLocations_Unpacked_ActiveRegion[0] = 1 * nx * ny_local + 2 * ny_local;
    HaloLocations_Unpacked_ActiveRegion[1] = 1 * nx * ny_local + 2 * ny_local + ny_local - 1;
    OctCentersY_Unpacked[0] = y_offset + 0.5;
    OctCentersY_Unpacked[1] = y_offset + ny_local - 0.5;
    for (int n = 0; n < 2; n++) {
        if (((n == 0) && (!(AtSouthBoundary))) || ((n == 1) && (!(AtNorthBoundary)))) {
            EXPECT_EQ(CellType_Host(HaloLocations_Unpacked_ActiveRegion[n]), Active);
            EXPECT_EQ(GrainID_Host(HaloLocations_Unpacked_ActiveRegion[n]), 29);
            EXPECT_FLOAT_EQ(DiagonalLength_Host(HaloLocations_Unpacked_ActiveRegion[n]), 0.01);
            EXPECT_FLOAT_EQ(DOCenter_Host(3 * HaloLocations_Unpacked_ActiveRegion[n]), 2.5);
            EXPECT_FLOAT_EQ(DOCenter_Host(3 * HaloLocations_Unpacked_ActiveRegion[n] + 1), OctCentersY_Unpacked[n]);
            EXPECT_FLOAT_EQ(DOCenter_Host(3 * HaloLocations_Unpacked_ActiveRegion[n] + 2), 1.5);
            // Were critical diagonal lengths correctly calculated from the unloaded buffer data?
            for (int l = 0; l < 26; l++) {
                EXPECT_FLOAT_EQ(CritDiagonalLength(26 * HaloLocations_Unpacked_ActiveRegion[n] + l),
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
    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // MPI rank locations relative to the global grid
    bool AtNorthBoundary, AtSouthBoundary;
    if (id == 0)
        AtSouthBoundary = true;
    else
        AtSouthBoundary = false;
    if (id == np - 1)
        AtNorthBoundary = true;
    else
        AtNorthBoundary = false;

    // 3 by 4 by 10 domain, only the top half of it is part of the current layer
    // Each rank has a portion of the domain subdivided in the Y direction
    int nx = 3;
    int ny_local = 4;
    int nz = 20;
    int z_layer_bottom = 5;
    int nz_layer = 10;
    int DomainSize = nx * ny_local * nz_layer;
    int DomainSize_AllLayers = nx * ny_local * nz;
    int NGrainOrientations = 10000;

    // Initialize empty inputs struct
    Inputs inputs;
    // Allocate device views: entire domain on each rank
    // Default to wall cells (CellType(index) = 0) with GrainID of 0
    CellData<device_memory_space> cellData(DomainSize_AllLayers, DomainSize, nx, ny_local, z_layer_bottom,
                                           inputs.substrate);
    Kokkos::deep_copy(cellData.CellType_AllLayers, Liquid);
    auto CellType = cellData.getCellTypeSubview();
    auto GrainID = cellData.getGrainIDSubview();

    // Allocate device views: only the active layer on each rank
    ViewF DiagonalLength("DiagonalLength", DomainSize);
    ViewF DOCenter(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * DomainSize);

    // Create send/receive buffers with a buffer size too small to hold all of the active cell data
    int BufSize = 1;
    Buffer2D BufferSouthSend(Kokkos::ViewAllocateWithoutInitializing("BufferSouthSend"), BufSize, 8);
    Buffer2D BufferNorthSend(Kokkos::ViewAllocateWithoutInitializing("BufferNorthSend"), BufSize, 8);
    Buffer2D BufferSouthRecv(Kokkos::ViewAllocateWithoutInitializing("BufferSouthRecv"), BufSize, 8);
    Buffer2D BufferNorthRecv(Kokkos::ViewAllocateWithoutInitializing("BufferNorthRecv"), BufSize, 8);
    // Init counts to 0 on device
    ViewI SendSizeSouth("SendSizeSouth", 1);
    ViewI SendSizeNorth("SendSizeNorth", 1);
    // Allocate counts on host (no init needed, will be copied from device)
    ViewI_H SendSizeSouth_Host(Kokkos::ViewAllocateWithoutInitializing("SendSizeSouth_Host"), 1);
    ViewI_H SendSizeNorth_Host(Kokkos::ViewAllocateWithoutInitializing("SendSizeNorth_Host"), 1);

    // Start with 2 cells in the current layer active (one in each buffer), GrainID equal to the X coordinate, Diagonal
    // length equal to the Y coordinate, octahedron center at (x + 0.5, y + 0.5, z + 0.5)
    Kokkos::parallel_for(
        "InitDomainActiveCellsNorth", 1, KOKKOS_LAMBDA(const int &) {
            int coord_x = 1;
            int coord_y = ny_local - 2;
            int coord_z = 0;
            int index = get1Dindex(coord_x, coord_y, coord_z, nx, ny_local);
            CellType(index) = Active;
            GrainID(index) = coord_x;
            int GhostGID = coord_x;
            DOCenter(3 * index) = coord_x + 0.5;
            float GhostDOCX = coord_x + 0.5;
            DOCenter(3 * index + 1) = coord_y + 0.5;
            float GhostDOCY = coord_y + 0.5;
            DOCenter(3 * index + 2) = coord_z + 0.5;
            float GhostDOCZ = coord_z + 0.5;
            DiagonalLength(index) = static_cast<float>(coord_y);
            float GhostDL = static_cast<float>(coord_y);
            // Load into appropriate buffers
            loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth, SendSizeSouth, ny_local,
                           coord_x, coord_y, coord_z, AtNorthBoundary, AtSouthBoundary, BufferSouthSend,
                           BufferNorthSend, NGrainOrientations, BufSize);
        });
    Kokkos::parallel_for(
        "InitDomainActiveCellsSouth", 1, KOKKOS_LAMBDA(const int &) {
            int coord_x = 1;
            int coord_y = 1;
            int coord_z = 0;
            int index = get1Dindex(coord_x, coord_y, coord_z, nx, ny_local);
            CellType(index) = Active;
            GrainID(index) = coord_x;
            int GhostGID = coord_x;
            DOCenter(3 * index) = coord_x + 0.5;
            float GhostDOCX = coord_x + 0.5;
            DOCenter(3 * index + 1) = coord_y + 0.5;
            float GhostDOCY = coord_y + 0.5;
            DOCenter(3 * index + 2) = coord_z + 0.5;
            float GhostDOCZ = coord_z + 0.5;
            DiagonalLength(index) = static_cast<float>(coord_y);
            float GhostDL = static_cast<float>(coord_y);
            // Load into appropriate buffers
            loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth, SendSizeSouth, ny_local,
                           coord_x, coord_y, coord_z, AtNorthBoundary, AtSouthBoundary, BufferSouthSend,
                           BufferNorthSend, NGrainOrientations, BufSize);
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
            int index = get1Dindex(coord_x, coord_y, coord_z, nx, ny_local);
            GrainID(index) = coord_x;
            int GhostGID = coord_x;
            DOCenter(3 * index) = coord_x + 0.5;
            float GhostDOCX = coord_x + 0.5;
            DOCenter(3 * index + 1) = coord_x + 0.5;
            float GhostDOCY = coord_y + 0.5;
            DOCenter(3 * index + 2) = coord_y + 0.5;
            float GhostDOCZ = coord_z + 0.5;
            DiagonalLength(index) = static_cast<float>(coord_y);
            float GhostDL = static_cast<float>(coord_y);
            // Attempt to load into appropriate buffers
            bool DataFitsInBuffer =
                loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth, SendSizeSouth,
                               ny_local, coord_x, coord_y, coord_z, AtNorthBoundary, AtSouthBoundary, BufferSouthSend,
                               BufferNorthSend, NGrainOrientations, BufSize);
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
            int index = get1Dindex(coord_x, coord_y, coord_z, nx, ny_local);
            GrainID(index) = coord_x;
            int GhostGID = coord_x;
            DOCenter(3 * index) = coord_x + 0.5;
            float GhostDOCX = coord_x + 0.5;
            DOCenter(3 * index + 1) = coord_y + 0.5;
            float GhostDOCY = coord_y + 0.5;
            DOCenter(3 * index + 2) = coord_z + 0.5;
            float GhostDOCZ = coord_z + 0.5;
            DiagonalLength(index) = static_cast<float>(coord_y);
            float GhostDL = static_cast<float>(coord_y);
            // Attempt to load into appropriate buffers
            bool DataFitsInBuffer =
                loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth, SendSizeSouth,
                               ny_local, coord_x, coord_y, coord_z, AtNorthBoundary, AtSouthBoundary, BufferSouthSend,
                               BufferNorthSend, NGrainOrientations, BufSize);
            if (!(DataFitsInBuffer)) {
                // This cell's data did not fit in the buffer with current size BufSize - mark with temporary type
                CellType(index) = ActiveFailedBufferLoad;
            }
        });

    // Attempt to resize buffers and load the remaining data
    int OldBufSize = BufSize;
    BufSize = ResizeBuffers(BufferNorthSend, BufferSouthSend, BufferNorthRecv, BufferSouthRecv, SendSizeNorth,
                            SendSizeSouth, SendSizeNorth_Host, SendSizeSouth_Host, OldBufSize);
    if (OldBufSize != BufSize)
        RefillBuffers(nx, nz_layer, ny_local, cellData, BufferNorthSend, BufferSouthSend, SendSizeNorth, SendSizeSouth,
                      AtNorthBoundary, AtSouthBoundary, DOCenter, DiagonalLength, NGrainOrientations, BufSize);
    // If there was 1 rank, buffer size should still be 1, as no data was loaded
    // Otherwise, 25 cells should have been added to the buffer in
    // addition to the required capacity increase during the resize
    if (np == 1)
        EXPECT_EQ(BufSize, 1);
    else {
        int CapacityInc = std::min(3, np - 1) + 25;
        EXPECT_EQ(BufSize, 1 + CapacityInc);
    }

    // Check that the correct information was retained in the original buffer positions
    // Check that the new information was correctly entered into the buffers
    Buffer2D_H BufferNorthSend_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), BufferNorthSend);
    Buffer2D_H BufferSouthSend_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), BufferSouthSend);
    // Check buffers other than Rank 0's south and Rank np-1's north
    if (!(AtSouthBoundary)) {
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
    if (!(AtNorthBoundary)) {
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

void testResetBufferCapacity() {

    // Init buffers to large size
    int BufSize = 50;
    Buffer2D BufferSouthSend(Kokkos::ViewAllocateWithoutInitializing("BufferSouthSend"), BufSize, 8);
    Buffer2D BufferNorthSend(Kokkos::ViewAllocateWithoutInitializing("BufferNorthSend"), BufSize, 8);
    Buffer2D BufferSouthRecv(Kokkos::ViewAllocateWithoutInitializing("BufferSouthRecv"), BufSize, 8);
    Buffer2D BufferNorthRecv(Kokkos::ViewAllocateWithoutInitializing("BufferNorthRecv"), BufSize, 8);

    // Fill buffers with test data
    Kokkos::parallel_for(
        "InitBuffers", BufSize, KOKKOS_LAMBDA(const int &i) {
            for (int j = 0; j < 8; j++) {
                BufferSouthSend(i, j) = i + j;
                BufferNorthSend(i, j) = i + j;
                BufferSouthRecv(i, j) = i;
                BufferNorthRecv(i, j) = i;
            }
        });

    // Reduce size back to the default of 25
    BufSize = 25;
    ResetBufferCapacity(BufferNorthSend, BufferSouthSend, BufferNorthRecv, BufferSouthRecv, BufSize);

    // Copy buffers back to host
    Buffer2D_H BufferNorthSend_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), BufferNorthSend);
    Buffer2D_H BufferSouthSend_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), BufferSouthSend);
    Buffer2D_H BufferNorthRecv_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), BufferNorthRecv);
    Buffer2D_H BufferSouthRecv_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), BufferSouthRecv);

    // Check that original values in first 25 x 8 positions were preserved
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < 8; j++) {
            EXPECT_EQ(BufferSouthSend_Host(i, j), i + j);
            EXPECT_EQ(BufferNorthSend_Host(i, j), i + j);
            EXPECT_EQ(BufferSouthRecv_Host(i, j), i);
            EXPECT_EQ(BufferNorthRecv_Host(i, j), i);
        }
    }
}

//---------------------------------------------------------------------------//
// domain_calculations
//---------------------------------------------------------------------------//
void testcalcVolFractionNucleated() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Get neighbor rank process IDs
    bool AtNorthBoundary, AtSouthBoundary;
    if (id == np - 1)
        AtNorthBoundary = true;
    else
        AtNorthBoundary = false;
    if (id == 0)
        AtSouthBoundary = true;
    else
        AtSouthBoundary = false;

    // Simulation domain
    int nx = 3;
    int ny_local = 3;
    int nz = 3;
    int DomainSize = nx * ny_local * nz;
    // Let all cells except those at Z = 0 have undergone solidification
    // Let the cells at Z = 1 consist of positive grain IDs, and those at Z = 2 of negative grain IDs
    ViewI_H GrainID_Host(Kokkos::ViewAllocateWithoutInitializing("GrainID"), DomainSize);
    ViewS_H LayerID_Host(Kokkos::ViewAllocateWithoutInitializing("LayerID"), DomainSize);
    for (int coord_z = 0; coord_z < nz; coord_z++) {
        for (int coord_x = 0; coord_x < nx; coord_x++) {
            for (int coord_y = 0; coord_y < ny_local; coord_y++) {
                int index = get1Dindex(coord_x, coord_y, coord_z, nx, ny_local);
                if (coord_z == 0)
                    LayerID_Host(index) = -1;
                else
                    LayerID_Host(index) = 0;
                if (coord_z == 2)
                    GrainID_Host(index) = -1;
                else
                    GrainID_Host(index) = 1;
            }
        }
    }
    ViewI GrainID = Kokkos::create_mirror_view_and_copy(device_memory_space(), GrainID_Host);
    ViewS LayerID = Kokkos::create_mirror_view_and_copy(device_memory_space(), LayerID_Host);

    // Perform calculation and compare to expected value (half of the solidified portion of the domain should consist of
    // nucleated grains, regardless of the number of MPI ranks used)
    float VolFractionNucleated =
        calcVolFractionNucleated(id, nx, ny_local, DomainSize, LayerID, GrainID, AtNorthBoundary, AtSouthBoundary);
    EXPECT_FLOAT_EQ(VolFractionNucleated, 0.5);
}

//---------------------------------------------------------------------------//
// full_simulations
//---------------------------------------------------------------------------//
void testSmallDirS() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    std::string InputFile = "Inp_SmallDirSolidification.json";

    // Run SmallDirS problem and check volume fraction of nucleated grains with 1% tolerance of expected value (to
    // account for the non-deterministic nature of the cell capture)
    RunProgram_Reduced(id, np, InputFile);

    // MPI barrier to ensure that log file has been written
    MPI_Barrier(MPI_COMM_WORLD);
    std::string LogFile = "TestProblemSmallDirS.json";
    std::ifstream LogDataStream(LogFile);
    nlohmann::json logdata = nlohmann::json::parse(LogDataStream);
    float VolFractionNucleated = logdata["Nucleation"]["VolFractionNucleated"];
    EXPECT_NEAR(VolFractionNucleated, 0.1784, 0.0100);
}

void testSmallEquiaxedGrain() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    std::string InputFile = "Inp_SmallEquiaxedGrain.json";

    // Run Small equiaxed grain problem and check time step at which the grain reaches the domain edge
    RunProgram_Reduced(id, np, InputFile);

    // MPI barrier to ensure that log file has been written
    MPI_Barrier(MPI_COMM_WORLD);
    std::string LogFile = "TestProblemSmallEquiaxedGrain.json";
    std::ifstream LogDataStream(LogFile);
    nlohmann::json logdata = nlohmann::json::parse(LogDataStream);
    int TimeStepOfOutput = logdata["TimeStepOfOutput"];
    EXPECT_EQ(TimeStepOfOutput, 4820);
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, communication) {
    testGhostNodes1D();
    testResizeRefillBuffers();
    testResetBufferCapacity();
}
TEST(TEST_CATEGORY, domain_calculations) { testcalcVolFractionNucleated(); }
TEST(TEST_CATEGORY, full_simulations) {
    testSmallDirS();
    testSmallEquiaxedGrain();
}
} // end namespace Test
