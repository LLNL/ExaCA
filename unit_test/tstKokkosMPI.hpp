// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>
#include <nlohmann/json.hpp>

#include "CAfunctions.hpp"
#include "CAghostnodes.hpp"
#include "CAinitialize.hpp"
#include "CAparsefiles.hpp"
#include "CAtypes.hpp"
#include "CAupdate.hpp"
#include "runCA.hpp"

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
    // The top half of the domain is the active portion of it
    int nx = 4;
    int MyXSlices = nx;
    int nz = 10;
    int ZBound_Low = 5;
    int nzActive = 5;
    // Domain is subdivided in Y, with ghost nodes between ranks
    // Each rank is size 3 in Y, plus ghost nodes if needed (i.e, not at problem boundary)
    // For example, if np = 4:
    // Rank 0: Y = 0, 1, 2, 3 (Y = 3 are ghost nodes)
    // Rank 1: Y = 2, 3, 4, 5, 6 (Y = 2, 6 are ghost nodes)
    // Rank 2: Y = 5, 6, 7, 8, 9 (Y = 5, 9 are ghost nodes)
    // Rank 3: Y = 8, 9, 10, 11 (Y = 8 are ghost nodes)
    int MyYSlices = 5; // assuming two sets of ghost nodes - removed if at boundary
    if (AtSouthBoundary)
        MyYSlices--;
    if (AtNorthBoundary)
        MyYSlices--;
    int MyYOffset;
    if (id == 0)
        MyYOffset = 0;
    else
        MyYOffset = 3 * (id - 1) + 2;
    int LocalDomainSize = MyXSlices * MyYSlices * nz;
    int LocalActiveDomainSize = MyXSlices * MyYSlices * nzActive;

    // Intialize grain orientations
    std::string GrainOrientationFile = checkFileInstalled("GrainOrientationVectors.csv", id);
    int NGrainOrientations = 10000; // Number of grain orientations considered in the simulation
    ViewF GrainUnitVector(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), 9 * NGrainOrientations);
    OrientationInit(id, NGrainOrientations, GrainUnitVector, GrainOrientationFile);

    // Initialize neighbor lists
    NList NeighborX, NeighborY, NeighborZ;
    NeighborListInit(NeighborX, NeighborY, NeighborZ);

    // Initialize host views - set initial GrainID values to 0, all CellType values to liquid
    ViewI_H CellType_Host(Kokkos::ViewAllocateWithoutInitializing("CellType_Host"), LocalDomainSize);
    Kokkos::deep_copy(CellType_Host, Liquid);
    ViewI_H GrainID_Host("GrainID_Host", LocalDomainSize);
    // Initialize active domain views to 0
    ViewF_H DiagonalLength_Host("DiagonalLength_Host", LocalActiveDomainSize);
    ViewF_H DOCenter_Host("DOCenter_Host", 3 * LocalActiveDomainSize);
    ViewF_H CritDiagonalLength_Host("CritDiagonalLength_Host", 26 * LocalActiveDomainSize);

    // Testing of loading of ghost nodes data, sending/receiving, unpacking, and calculations on ghost node data:
    // X = 2, Z = 6 is chosen for active cell placement on all ranks
    // Active cells will be located at Y = 1 and Y = MyYSlices-2 on each rank... these are located in the halo regions
    // and should be loaded into the send buffers
    int HaloLocations[2], HaloLocations_ActiveRegion[2];
    HaloLocations[0] = 6 * MyXSlices * MyYSlices + 2 * MyYSlices + 1;
    HaloLocations[1] = 6 * MyXSlices * MyYSlices + 2 * MyYSlices + MyYSlices - 2;
    HaloLocations_ActiveRegion[0] = (6 - ZBound_Low) * MyXSlices * MyYSlices + 2 * MyYSlices + 1;
    HaloLocations_ActiveRegion[1] = (6 - ZBound_Low) * MyXSlices * MyYSlices + 2 * MyYSlices + MyYSlices - 2;
    // Physical cell centers in Y are different for each location and each rank
    float OctCentersY[2];
    OctCentersY[0] = MyYOffset + 1.5;
    OctCentersY[1] = MyYOffset + MyYSlices - 1.5;
    for (int n = 0; n < 2; n++) {
        CellType_Host(HaloLocations[n]) = Active;
        GrainID_Host(HaloLocations[n]) = 29;
        // Initialize with an initial diagonal length, octahedra centered at cell center (X = 2.5, Y = varied, Z = 7.5)
        DiagonalLength_Host(HaloLocations_ActiveRegion[n]) = 0.01;
        DOCenter_Host(3 * HaloLocations_ActiveRegion[n]) = 2.5;
        DOCenter_Host(3 * HaloLocations_ActiveRegion[n] + 1) = OctCentersY[n];
        DOCenter_Host(3 * HaloLocations_ActiveRegion[n] + 2) = 6.5;
    }

    // Also testing an alternate situation where the ghost node data should NOT be unpacked, and the cells do not need
    // to be updated: X = 2, Z = 7 is chosen for active cell placement on all ranks Four active cells are located at Y =
    // 0, Y = 1, Y = MyYSlices - 2, and Y = MyYSlices - 1
    int HaloLocations_Alt[4], HaloLocations_Alt_ActiveRegion[4];
    HaloLocations_Alt[0] = 7 * MyXSlices * MyYSlices + 2 * MyYSlices + 0;
    HaloLocations_Alt[1] = 7 * MyXSlices * MyYSlices + 2 * MyYSlices + 1;
    HaloLocations_Alt[2] = 7 * MyXSlices * MyYSlices + 2 * MyYSlices + MyYSlices - 2;
    HaloLocations_Alt[3] = 7 * MyXSlices * MyYSlices + 2 * MyYSlices + MyYSlices - 1;
    HaloLocations_Alt_ActiveRegion[0] = (7 - ZBound_Low) * MyXSlices * MyYSlices + 2 * MyYSlices + 0;
    HaloLocations_Alt_ActiveRegion[1] = (7 - ZBound_Low) * MyXSlices * MyYSlices + 2 * MyYSlices + 1;
    HaloLocations_Alt_ActiveRegion[2] = (7 - ZBound_Low) * MyXSlices * MyYSlices + 2 * MyYSlices + MyYSlices - 2;
    HaloLocations_Alt_ActiveRegion[3] = (7 - ZBound_Low) * MyXSlices * MyYSlices + 2 * MyYSlices + MyYSlices - 1;
    // Physical cell centers in Y are different for each location and each rank
    float OctCentersY_Alt[4];
    OctCentersY_Alt[0] = MyYOffset + 0.5;
    OctCentersY_Alt[1] = MyYOffset + 1.5;
    OctCentersY_Alt[2] = MyYOffset + MyYSlices - 1.5;
    OctCentersY_Alt[3] = MyYOffset + MyYSlices - 0.5;
    for (int n = 0; n < 4; n++) {
        CellType_Host(HaloLocations_Alt[n]) = Active;
        GrainID_Host(HaloLocations_Alt[n]) = -id;
        // Initialize with an initial diagonal length, octahedra centered at cell center (X = 2.5, Y = varied, Z = 7.5)
        DiagonalLength_Host(HaloLocations_Alt_ActiveRegion[n]) = 0.01;
        DOCenter_Host(3 * HaloLocations_Alt_ActiveRegion[n]) = 2.5;
        DOCenter_Host(3 * HaloLocations_Alt_ActiveRegion[n] + 1) = OctCentersY_Alt[n];
        DOCenter_Host(3 * HaloLocations_Alt_ActiveRegion[n] + 2) = 7.5;
    }

    // Copy view data to the device
    ViewI CellType = Kokkos::create_mirror_view_and_copy(device_memory_space(), CellType_Host);
    ViewI GrainID = Kokkos::create_mirror_view_and_copy(device_memory_space(), GrainID_Host);
    ViewF DiagonalLength = Kokkos::create_mirror_view_and_copy(device_memory_space(), DiagonalLength_Host);
    ViewF CritDiagonalLength = Kokkos::create_mirror_view_and_copy(device_memory_space(), CritDiagonalLength_Host);
    ViewF DOCenter = Kokkos::create_mirror_view_and_copy(device_memory_space(), DOCenter_Host);

    // Buffer size large enough to hold all data
    int BufSize = MyXSlices * nzActive;

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
        "testloadghostnodes", LocalActiveDomainSize, KOKKOS_LAMBDA(const int &D3D1ConvPosition) {
            // 3D Coordinate of this cell on the "global" (all cells in the Z direction) grid
            int RankZ = D3D1ConvPosition / (MyXSlices * MyYSlices);
            int Rem = D3D1ConvPosition % (MyXSlices * MyYSlices);
            int RankX = Rem / MyYSlices;
            int RankY = Rem % MyYSlices;
            int GlobalZ = RankZ + ZBound_Low;
            int GlobalD3D1ConvPosition = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
            if (CellType(GlobalD3D1ConvPosition) == Active) {
                int GhostGID = GrainID(GlobalD3D1ConvPosition);
                float GhostDOCX = DOCenter(3 * D3D1ConvPosition);
                float GhostDOCY = DOCenter(3 * D3D1ConvPosition + 1);
                float GhostDOCZ = DOCenter(3 * D3D1ConvPosition + 2);
                float GhostDL = DiagonalLength(D3D1ConvPosition);
                loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth, SendSizeSouth,
                               MyYSlices, RankX, RankY, RankZ, AtNorthBoundary, AtSouthBoundary, BufferSouthSend,
                               BufferNorthSend, NGrainOrientations, BufSize);
            }
        });

    GhostNodes1D(0, id, NeighborRank_North, NeighborRank_South, MyXSlices, MyYSlices, MyYOffset, NeighborX, NeighborY,
                 NeighborZ, CellType, DOCenter, GrainID, GrainUnitVector, DiagonalLength, CritDiagonalLength,
                 NGrainOrientations, BufferNorthSend, BufferSouthSend, BufferNorthRecv, BufferSouthRecv, BufSize,
                 ZBound_Low, SendSizeNorth, SendSizeSouth);

    // Copy CellType, GrainID, DiagonalLength, DOCenter, CritDiagonalLength views to host to check values
    CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
    GrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);
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
    int HaloLocations_Unpacked[2], HaloLocations_Unpacked_ActiveRegion[2];
    float OctCentersY_Unpacked[2];
    HaloLocations_Unpacked[0] = 6 * MyXSlices * MyYSlices + 2 * MyYSlices;
    HaloLocations_Unpacked[1] = 6 * MyXSlices * MyYSlices + 2 * MyYSlices + MyYSlices - 1;
    HaloLocations_Unpacked_ActiveRegion[0] = (6 - ZBound_Low) * MyXSlices * MyYSlices + 2 * MyYSlices;
    HaloLocations_Unpacked_ActiveRegion[1] = (6 - ZBound_Low) * MyXSlices * MyYSlices + 2 * MyYSlices + MyYSlices - 1;
    OctCentersY_Unpacked[0] = MyYOffset + 0.5;
    OctCentersY_Unpacked[1] = MyYOffset + MyYSlices - 0.5;
    for (int n = 0; n < 2; n++) {
        if (((n == 0) && (!(AtSouthBoundary))) || ((n == 1) && (!(AtNorthBoundary)))) {
            EXPECT_EQ(CellType_Host(HaloLocations_Unpacked[n]), Active);
            EXPECT_EQ(GrainID_Host(HaloLocations_Unpacked[n]), 29);
            EXPECT_FLOAT_EQ(DiagonalLength_Host(HaloLocations_Unpacked_ActiveRegion[n]), 0.01);
            EXPECT_FLOAT_EQ(DOCenter_Host(3 * HaloLocations_Unpacked_ActiveRegion[n]), 2.5);
            EXPECT_FLOAT_EQ(DOCenter_Host(3 * HaloLocations_Unpacked_ActiveRegion[n] + 1), OctCentersY_Unpacked[n]);
            EXPECT_FLOAT_EQ(DOCenter_Host(3 * HaloLocations_Unpacked_ActiveRegion[n] + 2), 6.5);
            // Were critical diagonal lengths correctly calculated from the unloaded buffer data?
            for (int l = 0; l < 26; l++) {
                EXPECT_FLOAT_EQ(CritDiagonalLength(26 * HaloLocations_Unpacked_ActiveRegion[n] + l),
                                CritDiagonalLength_Expected[l]);
            }
        }
    }
    // These cells should not have been modified, as their data was unchanged
    for (int n = 0; n < 4; n++) {
        EXPECT_EQ(CellType_Host(HaloLocations_Alt[n]), Active);
        EXPECT_EQ(GrainID_Host(HaloLocations_Alt[n]), -id);
        EXPECT_FLOAT_EQ(DiagonalLength_Host(HaloLocations_Alt_ActiveRegion[n]), 0.01);
        EXPECT_FLOAT_EQ(DOCenter_Host(3 * HaloLocations_Alt_ActiveRegion[n]), 2.5);
        EXPECT_FLOAT_EQ(DOCenter_Host(3 * HaloLocations_Alt_ActiveRegion[n] + 1), OctCentersY_Alt[n]);
        EXPECT_FLOAT_EQ(DOCenter_Host(3 * HaloLocations_Alt_ActiveRegion[n] + 2), 7.5);
    }
}

void testResizeRefillBuffers() {
    using memory_space = TEST_MEMSPACE;
    using view_type_float = Kokkos::View<float *, memory_space>;
    using view_type_int = Kokkos::View<int *, memory_space>;

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
    int MyYSlices = 4;
    int nz = 20;
    int ZBound_Low = 5;
    int nzActive = 10;
    int LocalActiveDomainSize = nx * MyYSlices * nzActive;
    int LocalDomainSize = nx * MyYSlices * nz;
    int NGrainOrientations = 10000;

    // Allocate device views: entire domain on each rank
    // Default to wall cells (CellType(index) = 0) with GrainID of 0
    view_type_int GrainID("GrainID", LocalDomainSize);
    view_type_int CellType("CellType", LocalDomainSize);
    // Allocate device views: only the active layer on each rank
    view_type_float DiagonalLength("DiagonalLength", LocalActiveDomainSize);
    view_type_float DOCenter(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * LocalActiveDomainSize);

    // Create send/receive buffers with a buffer size too small to hold all of the active cell data
    int BufSize = 1;
    Buffer2D BufferSouthSend(Kokkos::ViewAllocateWithoutInitializing("BufferSouthSend"), BufSize, 8);
    Buffer2D BufferNorthSend(Kokkos::ViewAllocateWithoutInitializing("BufferNorthSend"), BufSize, 8);
    Buffer2D BufferSouthRecv(Kokkos::ViewAllocateWithoutInitializing("BufferSouthRecv"), BufSize, 8);
    Buffer2D BufferNorthRecv(Kokkos::ViewAllocateWithoutInitializing("BufferNorthRecv"), BufSize, 8);
    // Init counts to 0 on device
    view_type_int SendSizeSouth("SendSizeSouth", 1);
    view_type_int SendSizeNorth("SendSizeNorth", 1);
    // Allocate counts on host (no init needed, will be copied from device)
    ViewI_H SendSizeSouth_Host(Kokkos::ViewAllocateWithoutInitializing("SendSizeSouth_Host"), 1);
    ViewI_H SendSizeNorth_Host(Kokkos::ViewAllocateWithoutInitializing("SendSizeNorth_Host"), 1);

    // Start with 2 cells in the current layer active (one in each buffer), GrainID equal to the X coordinate, Diagonal
    // length equal to the Y coordinate, octahedron center at (x + 0.5, y + 0.5, z + 0.5)
    Kokkos::parallel_for(
        "InitDomainActiveCellsNorth", 1, KOKKOS_LAMBDA(const int &) {
            int RankX = 1;
            int RankY = MyYSlices - 2;
            int RankZ = 0;
            int D3D1ConvPosition = RankZ * nx * MyYSlices + RankX * MyYSlices + RankY;
            int GlobalD3D1ConvPosition = (RankZ + ZBound_Low) * nx * MyYSlices + RankX * MyYSlices + RankY;
            CellType(GlobalD3D1ConvPosition) = Active;
            GrainID(GlobalD3D1ConvPosition) = RankX;
            int GhostGID = RankX;
            DOCenter(3 * D3D1ConvPosition) = RankX + 0.5;
            float GhostDOCX = RankX + 0.5;
            DOCenter(3 * D3D1ConvPosition + 1) = RankY + 0.5;
            float GhostDOCY = RankY + 0.5;
            DOCenter(3 * D3D1ConvPosition + 2) = RankZ + 0.5;
            float GhostDOCZ = RankZ + 0.5;
            DiagonalLength(D3D1ConvPosition) = static_cast<float>(RankY);
            float GhostDL = static_cast<float>(RankY);
            // Load into appropriate buffers
            bool DataFitsInBuffer =
                loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth, SendSizeSouth,
                               MyYSlices, RankX, RankY, RankZ, AtNorthBoundary, AtSouthBoundary, BufferSouthSend,
                               BufferNorthSend, NGrainOrientations, BufSize);
            if (!(DataFitsInBuffer))
                throw std::runtime_error("Error: North send buffer overflowed too soon in testResizeRefillBuffers");
        });
    Kokkos::parallel_for(
        "InitDomainActiveCellsSouth", 1, KOKKOS_LAMBDA(const int &) {
            int RankX = 1;
            int RankY = 1;
            int RankZ = 0;
            int D3D1ConvPosition = RankZ * nx * MyYSlices + RankX * MyYSlices + RankY;
            int GlobalD3D1ConvPosition = (RankZ + ZBound_Low) * nx * MyYSlices + RankX * MyYSlices + RankY;
            CellType(GlobalD3D1ConvPosition) = Active;
            GrainID(GlobalD3D1ConvPosition) = RankX;
            int GhostGID = RankX;
            DOCenter(3 * D3D1ConvPosition) = RankX + 0.5;
            float GhostDOCX = RankX + 0.5;
            DOCenter(3 * D3D1ConvPosition + 1) = RankY + 0.5;
            float GhostDOCY = RankY + 0.5;
            DOCenter(3 * D3D1ConvPosition + 2) = RankZ + 0.5;
            float GhostDOCZ = RankZ + 0.5;
            DiagonalLength(D3D1ConvPosition) = static_cast<float>(RankY);
            float GhostDL = static_cast<float>(RankY);
            // Load into appropriate buffers
            bool DataFitsInBuffer =
                loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth, SendSizeSouth,
                               MyYSlices, RankX, RankY, RankZ, AtNorthBoundary, AtSouthBoundary, BufferSouthSend,
                               BufferNorthSend, NGrainOrientations, BufSize);
            if (!(DataFitsInBuffer))
                throw std::runtime_error("Error: South send buffer overflowed too soon in testResizeRefillBuffers");
        });

    // Each rank will have "id % 4" cells of additional data to send to the south, and 1 cell of additional data to send
    // to the north. These cells will be marked as having failed to be loaded into the buffers, as the buffer sizes are
    // only 1
    int FailedBufferCellsSouth = id % 4;
    int FailedBufferCellsNorth = 1;
    Kokkos::parallel_for(
        "InitDomainFailedNorth", FailedBufferCellsNorth, KOKKOS_LAMBDA(const int &k) {
            int RankX = 2;
            int RankY = MyYSlices - 2;
            int RankZ = k;
            int D3D1ConvPosition = RankZ * nx * MyYSlices + RankX * MyYSlices + RankY;
            int GlobalD3D1ConvPosition = (RankZ + ZBound_Low) * nx * MyYSlices + RankX * MyYSlices + RankY;
            GrainID(GlobalD3D1ConvPosition) = RankX;
            int GhostGID = RankX;
            DOCenter(3 * D3D1ConvPosition) = RankX + 0.5;
            float GhostDOCX = RankX + 0.5;
            DOCenter(3 * D3D1ConvPosition + 1) = RankY + 0.5;
            float GhostDOCY = RankY + 0.5;
            DOCenter(3 * D3D1ConvPosition + 2) = RankZ + 0.5;
            float GhostDOCZ = RankZ + 0.5;
            DiagonalLength(D3D1ConvPosition) = static_cast<float>(RankY);
            float GhostDL = static_cast<float>(RankY);
            // Attempt to load into appropriate buffers
            bool DataFitsInBuffer =
                loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth, SendSizeSouth,
                               MyYSlices, RankX, RankY, RankZ, AtNorthBoundary, AtSouthBoundary, BufferSouthSend,
                               BufferNorthSend, NGrainOrientations, BufSize);
            if (!(DataFitsInBuffer)) {
                // This cell's data did not fit in the buffer with current size BufSize - mark with temporary type
                CellType(GlobalD3D1ConvPosition) = ActiveFailedBufferLoad;
            }
        });

    Kokkos::parallel_for(
        "InitDomainFailedSouth", FailedBufferCellsSouth, KOKKOS_LAMBDA(const int &k) {
            int RankX = 2;
            int RankY = 1;
            int RankZ = k;
            int D3D1ConvPosition = RankZ * nx * MyYSlices + RankX * MyYSlices + RankY;
            int GlobalD3D1ConvPosition = (RankZ + ZBound_Low) * nx * MyYSlices + RankX * MyYSlices + RankY;
            GrainID(GlobalD3D1ConvPosition) = RankX;
            int GhostGID = RankX;
            DOCenter(3 * D3D1ConvPosition) = RankX + 0.5;
            float GhostDOCX = RankX + 0.5;
            DOCenter(3 * D3D1ConvPosition + 1) = RankY + 0.5;
            float GhostDOCY = RankY + 0.5;
            DOCenter(3 * D3D1ConvPosition + 2) = RankZ + 0.5;
            float GhostDOCZ = RankZ + 0.5;
            DiagonalLength(D3D1ConvPosition) = static_cast<float>(RankY);
            float GhostDL = static_cast<float>(RankY);
            // Attempt to load into appropriate buffers
            bool DataFitsInBuffer =
                loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth, SendSizeSouth,
                               MyYSlices, RankX, RankY, RankZ, AtNorthBoundary, AtSouthBoundary, BufferSouthSend,
                               BufferNorthSend, NGrainOrientations, BufSize);
            if (!(DataFitsInBuffer)) {
                // This cell's data did not fit in the buffer with current size BufSize - mark with temporary type
                CellType(GlobalD3D1ConvPosition) = ActiveFailedBufferLoad;
            }
        });

    // Attempt to resize buffers and load the remaining data
    int OldBufSize = BufSize;
    BufSize = ResizeBuffers(BufferNorthSend, BufferSouthSend, BufferNorthRecv, BufferSouthRecv, SendSizeNorth,
                            SendSizeSouth, SendSizeNorth_Host, SendSizeSouth_Host, OldBufSize);
    if (OldBufSize != BufSize)
        RefillBuffers(nx, nzActive, MyYSlices, ZBound_Low, CellType, BufferNorthSend, BufferSouthSend, SendSizeNorth,
                      SendSizeSouth, AtNorthBoundary, AtSouthBoundary, GrainID, DOCenter, DiagonalLength,
                      NGrainOrientations, BufSize);
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
        EXPECT_FLOAT_EQ(BufferNorthSend_Host(0, 0), 1.0);                                 // RankX
        EXPECT_FLOAT_EQ(BufferNorthSend_Host(0, 1), 0.0);                                 // RankZ
        EXPECT_FLOAT_EQ(BufferNorthSend_Host(0, 2), 1.0);                                 // grain orientation
        EXPECT_FLOAT_EQ(BufferNorthSend_Host(0, 3), 1.0);                                 // grain number
        EXPECT_FLOAT_EQ(BufferNorthSend_Host(0, 4), 1.5);                                 // oct center X
        EXPECT_FLOAT_EQ(BufferNorthSend_Host(0, 5), static_cast<float>(MyYSlices - 1.5)); // oct center Y
        EXPECT_FLOAT_EQ(BufferNorthSend_Host(0, 6), 0.5);                                 // oct center Z
        EXPECT_FLOAT_EQ(BufferNorthSend_Host(0, 7), static_cast<float>(MyYSlices - 2));   // diagonal length
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
    int MyYSlices = 3;
    int nz = 3;
    int LocalDomainSize = nx * MyYSlices * nz;
    // Let all cells except those at Z = 0 have undergone solidification
    // Let the cells at Z = 1 consist of positive grain IDs, and those at Z = 2 of negative grain IDs
    ViewI_H GrainID_Host(Kokkos::ViewAllocateWithoutInitializing("GrainID"), LocalDomainSize);
    ViewI_H LayerID_Host(Kokkos::ViewAllocateWithoutInitializing("LayerID"), LocalDomainSize);
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                int D3D1ConvPosition = k * nx * MyYSlices + i * MyYSlices + j;
                if (k == 0)
                    LayerID_Host(D3D1ConvPosition) = -1;
                else
                    LayerID_Host(D3D1ConvPosition) = 0;
                if (k == 2)
                    GrainID_Host(D3D1ConvPosition) = -1;
                else
                    GrainID_Host(D3D1ConvPosition) = 1;
            }
        }
    }
    ViewI GrainID = Kokkos::create_mirror_view_and_copy(device_memory_space(), GrainID_Host);
    ViewI LayerID = Kokkos::create_mirror_view_and_copy(device_memory_space(), LayerID_Host);

    // Perform calculation and compare to expected value (half of the solidified portion of the domain should consist of
    // nucleated grains, regardless of the number of MPI ranks used)
    float VolFractionNucleated = calcVolFractionNucleated(id, nx, MyYSlices, LocalDomainSize, LayerID, GrainID,
                                                          AtNorthBoundary, AtSouthBoundary);
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
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, communication) {
    testGhostNodes1D();
    testResizeRefillBuffers();
    testResetBufferCapacity();
}
TEST(TEST_CATEGORY, domain_calculations) { testcalcVolFractionNucleated(); }
TEST(TEST_CATEGORY, full_simulations) { testSmallDirS(); }
} // end namespace Test
