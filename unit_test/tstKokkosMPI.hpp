// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAfunctions.hpp"
#include "CAghostnodes.hpp"
#include "CAinitialize.hpp"
#include "CAparsefiles.hpp"
#include "CAtypes.hpp"
#include "CAupdate.hpp"

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
    int MyXOffset = 0;
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

    // Buffer sizes
    int BufSizeX = MyXSlices;
    int BufSizeZ = nzActive;

    // Send/recv buffers for ghost node data should be initialized with zeros
    Buffer2D BufferSouthSend("BufferSouthSend", BufSizeX * BufSizeZ, 5);
    Buffer2D BufferNorthSend("BufferNorthSend", BufSizeX * BufSizeZ, 5);
    Buffer2D BufferSouthRecv("BufferSouthRecv", BufSizeX * BufSizeZ, 5);
    Buffer2D BufferNorthRecv("BufferNorthRecv", BufSizeX * BufSizeZ, 5);

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
                double GhostGID = static_cast<double>(GrainID(GlobalD3D1ConvPosition));
                double GhostDOCX = static_cast<double>(DOCenter(3 * D3D1ConvPosition));
                double GhostDOCY = static_cast<double>(DOCenter(3 * D3D1ConvPosition + 1));
                double GhostDOCZ = static_cast<double>(DOCenter(3 * D3D1ConvPosition + 2));
                double GhostDL = static_cast<double>(DiagonalLength(D3D1ConvPosition));
                loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, BufSizeX, MyYSlices, RankX, RankY,
                               RankZ, AtNorthBoundary, AtSouthBoundary, BufferSouthSend, BufferNorthSend);
            }
        });

    // Perform halo exchange in 1D
    GhostNodes1D(0, 0, NeighborRank_North, NeighborRank_South, MyXSlices, MyYSlices, MyXOffset, MyYOffset, NeighborX,
                 NeighborY, NeighborZ, CellType, DOCenter, GrainID, GrainUnitVector, DiagonalLength, CritDiagonalLength,
                 NGrainOrientations, BufferNorthSend, BufferSouthSend, BufferNorthRecv, BufferSouthRecv, BufSizeX, 0,
                 BufSizeZ, ZBound_Low);

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

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, communication) { testGhostNodes1D(); }
} // end namespace Test
