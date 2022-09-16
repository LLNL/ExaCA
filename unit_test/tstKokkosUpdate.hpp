// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAinitialize.hpp"
#include "CAtypes.hpp"
#include "CAupdate.hpp"

#include <gtest/gtest.h>

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
void testNucleation() {

    int SuccessfulNucEvents_ThisRank = 0; // nucleation event counter
    // Counters for nucleation events (successful or not) - host and device
    int NucleationCounter = 0;

    // Create views - each rank has 125 cells, 75 of which are part of the active region of the domain
    int nx = 5;
    int MyYSlices = 5;
    int nz = 5;
    int ZBound_Low = 2;
    int nzActive = 3;
    int LocalActiveDomainSize = nx * MyYSlices * nzActive;
    int LocalDomainSize = nx * MyYSlices * nz;

    // All cells have GrainID of 0, CellType of Liquid - with the exception of the locations where the nucleation events
    // are unable to occur
    ViewI_H GrainID_Host("GrainID_Host", LocalDomainSize);
    ViewI_H CellType_Host(Kokkos::ViewAllocateWithoutInitializing("CellType_Host"), LocalDomainSize);
    Kokkos::deep_copy(CellType_Host, Liquid);

    // Create test nucleation data - 10 possible events per MPI rank
    int PossibleNuclei_ThisRankThisLayer = 10;
    ViewI_H NucleationTimes_Host(Kokkos::ViewAllocateWithoutInitializing("NucleationTimes_Host"),
                                 PossibleNuclei_ThisRankThisLayer);
    ViewI_H NucleiLocation_Host(Kokkos::ViewAllocateWithoutInitializing("NucleiLocation_Host"),
                                PossibleNuclei_ThisRankThisLayer);
    ViewI_H NucleiGrainID_Host(Kokkos::ViewAllocateWithoutInitializing("NucleiGrainID_Host"),
                               PossibleNuclei_ThisRankThisLayer);
    for (int n = 0; n < PossibleNuclei_ThisRankThisLayer; n++) {
        // NucleationTimes values should be in the order in which the events occur - let these times depend on the
        // process id
        NucleationTimes_Host(n) = 4 + n;
        NucleiLocation_Host(n) = ZBound_Low * nx * MyYSlices + n;
        // Give these nucleation events grain IDs based on their order, starting with -1 and counting down
        NucleiGrainID_Host(n) = -(n + 1);
    }
    // Include the case where 2 potential nucleation events (3 and 4) happen on the same time step - both successful
    NucleationTimes_Host(3) = NucleationTimes_Host(4);

    // Include the case where a potential nucleation event (2) is unsuccessful
    int UnsuccessfulLocA = ZBound_Low * nx * MyYSlices + 2;
    CellType_Host(UnsuccessfulLocA) = Active;
    GrainID_Host(UnsuccessfulLocA) = 1;

    // Include the case where 2 potential nucleation events (5 and 6) happen on the same time step - the first
    // successful, the second unsuccessful
    NucleationTimes_Host(5) = NucleationTimes_Host(6);
    int UnsuccessfulLocC = ZBound_Low * nx * MyYSlices + 6;
    CellType_Host(UnsuccessfulLocC) = Active;
    GrainID_Host(UnsuccessfulLocC) = 2;

    // Include the case where 2 potential nucleation events (8 and 9) happen on the same time step - the first
    // unsuccessful, the second successful
    NucleationTimes_Host(8) = NucleationTimes_Host(9);
    int UnsuccessfulLocB = ZBound_Low * nx * MyYSlices + 8;
    CellType_Host(UnsuccessfulLocB) = Active;
    GrainID_Host(UnsuccessfulLocB) = 3;

    // Copy host views to device
    using memory_space = Kokkos::DefaultExecutionSpace::memory_space;
    ViewI CellType = Kokkos::create_mirror_view_and_copy(memory_space(), CellType_Host);
    ViewI GrainID = Kokkos::create_mirror_view_and_copy(memory_space(), GrainID_Host);
    ViewI NucleiLocation = Kokkos::create_mirror_view_and_copy(memory_space(), NucleiLocation_Host);
    ViewI NucleiGrainID = Kokkos::create_mirror_view_and_copy(memory_space(), NucleiGrainID_Host);

    // Steering Vector
    ViewI SteeringVector(Kokkos::ViewAllocateWithoutInitializing("SteeringVector"), LocalActiveDomainSize);
    // Initialize steering vector size to 0
    ViewI numSteer("SteeringVector", 1);

    // Take enough time steps such that every nucleation event has a chance to occur
    for (int cycle = 0; cycle <= (15); cycle++) {
        Nucleation(cycle, SuccessfulNucEvents_ThisRank, NucleationCounter, PossibleNuclei_ThisRankThisLayer,
                   NucleationTimes_Host, NucleiLocation, NucleiGrainID, CellType, GrainID, ZBound_Low, nx, MyYSlices,
                   SteeringVector, numSteer);
    }

    // Copy CellType, SteeringVector, numSteer, GrainID back to host to check nucleation results
    CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
    ViewI_H SteeringVector_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), SteeringVector);
    ViewI_H numSteer_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), numSteer);
    GrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);

    // Check that all 10 possible nucleation events were attempted
    EXPECT_EQ(NucleationCounter, 10);
    // Check that 7 of the 10 nucleation events were successful
    EXPECT_EQ(SuccessfulNucEvents_ThisRank, 7);
    EXPECT_EQ(numSteer_Host(0), 7);

    // Ensure that the 3 events that should not have occurred, did not occur
    // These cells should be untouched - active type, and same GrainID they was initialized with
    EXPECT_EQ(CellType_Host(ZBound_Low * nx * MyYSlices + 2), Active);
    EXPECT_EQ(GrainID_Host(ZBound_Low * nx * MyYSlices + 2), 1);
    EXPECT_EQ(CellType_Host(ZBound_Low * nx * MyYSlices + 6), Active);
    EXPECT_EQ(GrainID_Host(ZBound_Low * nx * MyYSlices + 6), 2);
    EXPECT_EQ(CellType_Host(ZBound_Low * nx * MyYSlices + 8), Active);
    EXPECT_EQ(GrainID_Host(ZBound_Low * nx * MyYSlices + 8), 3);

    // Check that the successful nucleation events occurred as expected
    std::vector<int> SuccessfulNuc_GrainIDs{-1, -2, -4, -5, -6, -8, -10};
    for (int n = 0; n < 7; n++) {
        // Cell's location within the active layer
        int CellLocation_ThisLayer = SteeringVector_Host(n);
        int RankZ = CellLocation_ThisLayer / (nx * MyYSlices);
        int Rem = CellLocation_ThisLayer % (nx * MyYSlices);
        int RankX = Rem / MyYSlices;
        int RankY = Rem % MyYSlices;
        int GlobalZ = RankZ + ZBound_Low;
        // Cell location within the overall domain
        int CellLocation_AllLayers = GlobalZ * nx * MyYSlices + RankX * MyYSlices + RankY;
        // Check that this cell type is marked as "FutureActive", and the GrainID matches the expected value
        EXPECT_EQ(CellType_Host(CellLocation_AllLayers), FutureActive);
        EXPECT_EQ(GrainID_Host(CellLocation_AllLayers), SuccessfulNuc_GrainIDs[n]);
    }
}

void testFillSteeringVector_Remelt() {

    // Create views - each rank has 125 cells, 75 of which are part of the active region of the domain
    int nx = 5;
    int MyYSlices = 5;
    int MyYOffset = 5;
    int nz = 5;
    int ZBound_Low = 2;
    int nzActive = 3;
    int LocalActiveDomainSize = nx * MyYSlices * nzActive;
    int LocalDomainSize = nx * MyYSlices * nz;

    // MPI rank locations relative to the global grid
    bool AtNorthBoundary = false;
    bool AtSouthBoundary = false;

    // Buffers for ghost node data (fixed size)
    int BufSizeX = nx;
    int BufSizeZ = nzActive;

    // Send/recv buffers for ghost node data - initialize with values of 1.0
    Buffer2D BufferSouthSend(Kokkos::ViewAllocateWithoutInitializing("BufferSouthSend"), BufSizeX * BufSizeZ, 5);
    Buffer2D BufferNorthSend(Kokkos::ViewAllocateWithoutInitializing("BufferNorthSend"), BufSizeX * BufSizeZ, 5);
    Kokkos::deep_copy(BufferSouthSend, 1.0);
    Kokkos::deep_copy(BufferNorthSend, 1.0);

    // Initialize neighbor lists
    NList NeighborX, NeighborY, NeighborZ;
    NeighborListInit(NeighborX, NeighborY, NeighborZ);

    ViewI_H GrainID_Host(Kokkos::ViewAllocateWithoutInitializing("GrainID_Host"), LocalDomainSize);
    ViewI_H CellType_Host(Kokkos::ViewAllocateWithoutInitializing("CellType_Host"), LocalDomainSize);
    ViewI_H MeltTimeStep_Host(Kokkos::ViewAllocateWithoutInitializing("MeltTimeStep_Host"), LocalDomainSize);
    ViewI_H CritTimeStep_Host(Kokkos::ViewAllocateWithoutInitializing("CritTimeStep_Host"), LocalDomainSize);
    ViewF_H UndercoolingChange_Host(Kokkos::ViewAllocateWithoutInitializing("UndercoolingChange_Host"),
                                    LocalDomainSize);
    ViewF_H UndercoolingCurrent_Host("UndercoolingCurrent_Host",
                                     LocalDomainSize); // initialize to 0, no initial undercooling

    for (int i = 0; i < LocalDomainSize; i++) {
        // Cell coordinates on this rank in X, Y, and Z (GlobalZ = relative to domain bottom)
        int GlobalZ = i / (nx * MyYSlices);
        int Rem = i % (nx * MyYSlices);
        int RankY = Rem % MyYSlices;
        // Let cells be assigned GrainIDs based on the rank ID
        // Cells have grain ID 1
        GrainID_Host(i) = 1;
        // Cells at Z = 0 through Z = 2 are Solid, Z = 3 and 4 are TempSolid
        if (GlobalZ <= 2)
            CellType_Host(i) = Solid;
        else
            CellType_Host(i) = TempSolid;
        // Cells "melt" at a time step corresponding to their Y location in the overall domain (depends on MyYOffset of
        // the rank)
        MeltTimeStep_Host(i) = RankY + MyYOffset + 1;
        // Cells reach liquidus during cooling 2 time steps after melting
        CritTimeStep_Host(i) = MeltTimeStep_Host(i) + 2;
        // Let the cooling rate of the cells from the liquidus depend on the rank ID
        UndercoolingChange_Host(i) = 0.2;
    }

    // Steering Vector
    ViewI SteeringVector(Kokkos::ViewAllocateWithoutInitializing("SteeringVector"), LocalActiveDomainSize);
    ViewI_H numSteer_Host(Kokkos::ViewAllocateWithoutInitializing("SteeringVectorSize"), 1);
    numSteer_Host(0) = 0;

    // Copy views to device for test
    ViewI numSteer = Kokkos::create_mirror_view_and_copy(device_memory_space(), numSteer_Host);
    ViewI GrainID = Kokkos::create_mirror_view_and_copy(device_memory_space(), GrainID_Host);
    ViewI CellType = Kokkos::create_mirror_view_and_copy(device_memory_space(), CellType_Host);
    ViewI MeltTimeStep = Kokkos::create_mirror_view_and_copy(device_memory_space(), MeltTimeStep_Host);
    ViewI CritTimeStep = Kokkos::create_mirror_view_and_copy(device_memory_space(), CritTimeStep_Host);
    ViewF UndercoolingChange = Kokkos::create_mirror_view_and_copy(device_memory_space(), UndercoolingChange_Host);
    ViewF UndercoolingCurrent = Kokkos::create_mirror_view_and_copy(device_memory_space(), UndercoolingCurrent_Host);

    int numcycles = 15;
    for (int cycle = 1; cycle <= numcycles; cycle++) {
        // Update cell types, local undercooling each time step, and fill the steering vector
        FillSteeringVector_Remelt(cycle, LocalActiveDomainSize, nx, MyYSlices, NeighborX, NeighborY, NeighborZ,
                                  CritTimeStep, UndercoolingCurrent, UndercoolingChange, CellType, GrainID, ZBound_Low,
                                  nzActive, SteeringVector, numSteer, numSteer_Host, MeltTimeStep, BufSizeX,
                                  AtNorthBoundary, AtSouthBoundary, BufferNorthSend, BufferSouthSend);
    }

    // Copy CellType, SteeringVector, numSteer, UndercoolingCurrent, Buffers back to host to check steering vector
    // construction results
    CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
    ViewI_H SteeringVector_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), SteeringVector);
    numSteer_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), numSteer);
    UndercoolingCurrent_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), UndercoolingCurrent);
    Buffer2D_H BufferSouthSend_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), BufferSouthSend);
    Buffer2D_H BufferNorthSend_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), BufferNorthSend);

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
    for (int i = 0; i < LocalDomainSize; i++) {
        // Cell coordinates on this rank in X, Y, and Z (GlobalZ = relative to domain bottom)
        int GlobalZ = i / (nx * MyYSlices);
        if (GlobalZ <= 2) {
            EXPECT_EQ(CellType_Host(i), Solid);
            EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(i), 0.0);
        }
        else {
            if (numcycles < MeltTimeStep_Host(i)) {
                EXPECT_EQ(CellType_Host(i), TempSolid);
                EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(i), 0.0);
            }
            else if ((numcycles >= MeltTimeStep_Host(i)) && (numcycles <= CritTimeStep_Host(i))) {
                EXPECT_EQ(CellType_Host(i), Liquid);
                EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(i), 0.0);
            }
            else {
                EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(i), (numcycles - CritTimeStep_Host(i)) * 0.2);
                if (GlobalZ == 4)
                    EXPECT_EQ(CellType_Host(i), Liquid);
                else {
                    EXPECT_EQ(CellType_Host(i), FutureActive);
                    FutureActiveCells++;
                }
            }
        }
    }
    // Check the steering vector values on the host
    EXPECT_EQ(FutureActiveCells, numSteer_Host(0));
    for (int i = 0; i < FutureActiveCells; i++) {
        // This cell should correspond to a cell at GlobalZ = 3 (RankZ = 1), and some X and Y
        int LowerBoundCellLocation = nx * MyYSlices - 1;
        int UpperBoundCellLocation = 2 * nx * MyYSlices;
        EXPECT_GT(SteeringVector_Host(i), LowerBoundCellLocation);
        EXPECT_LT(SteeringVector_Host(i), UpperBoundCellLocation);
    }

    // Check that the buffer values were either appropriately set to zeros (if the cell underwent melting) or remained
    // at 1.0
    for (int i = 0; i < nx * BufSizeZ; i++) {
        int RankZ = i / nx;
        int GlobalZ = RankZ + ZBound_Low;
        int RankX = i % nx;
        int NorthCellCoordinate = GlobalZ * nx * MyYSlices + RankX * MyYSlices + (MyYSlices - 1);
        if ((CellType_Host(NorthCellCoordinate) == TempSolid) || (CellType_Host(NorthCellCoordinate) == Solid)) {
            for (int j = 0; j < 5; j++) {
                EXPECT_EQ(BufferNorthSend_Host(i, j), 1.0);
            }
        }
        else {
            for (int j = 0; j < 5; j++) {
                EXPECT_EQ(BufferNorthSend_Host(i, j), 0.0);
            }
        }
        int SouthCellCoordinate = GlobalZ * nx * MyYSlices + RankX * MyYSlices + (MyYSlices - 1);
        if ((CellType_Host(SouthCellCoordinate) == TempSolid) || (CellType_Host(SouthCellCoordinate) == Solid)) {
            for (int j = 0; j < 5; j++) {
                EXPECT_EQ(BufferNorthSend_Host(i, j), 1.0);
            }
        }
        else {
            for (int j = 0; j < 5; j++) {
                EXPECT_EQ(BufferNorthSend_Host(i, j), 0.0);
            }
        }
    }
}

void testcalcCritDiagonalLength() {

    // 2 "cells" in the domain
    int LocalDomainSize = 2;
    // First orientation is orientation 12 (starting indexing at 0) of GrainOrientationVectors.csv
    // Second orientation is orientation 28 (starting indexing at 0)
    std::vector<float> GrainUnitVectorV{0.52877,  0.651038, -0.544565, -0.572875, 0.747163,  0.336989,
                                        0.626272, 0.133778, 0.768041,  0.736425,  -0.530983, 0.419208,
                                        0.664512, 0.683971, -0.301012, -0.126894, 0.500241,  0.856538};
    // Load unit vectors into view
    ViewF_H GrainUnitVector_Host(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector_Host"), 9 * LocalDomainSize);
    for (int i = 0; i < 9 * LocalDomainSize; i++) {
        GrainUnitVector_Host(i) = GrainUnitVectorV[i];
    }

    // Initialize neighbor lists
    NList NeighborX, NeighborY, NeighborZ;
    NeighborListInit(NeighborX, NeighborY, NeighborZ);

    // Load octahedron centers into view
    ViewF_H DOCenter_Host(Kokkos::ViewAllocateWithoutInitializing("DOCenter_Host"), 3 * LocalDomainSize);
    // Octahedron center for first cell (does not align with CA cell center, which is at 31.5, 3.5, 1.5)
    DOCenter_Host(0) = 30.611142;
    DOCenter_Host(1) = 3.523741;
    DOCenter_Host(2) = 0.636301;
    // Octahedron center for second cell (aligns with CA cell center at 110.5, 60.5, 0.5))
    DOCenter_Host(3) = 110.5;
    DOCenter_Host(4) = 60.5;
    DOCenter_Host(5) = 0.5;
    // Copy octahedron center and grain unit vector data to the device
    ViewF DOCenter = Kokkos::create_mirror_view_and_copy(device_memory_space(), DOCenter_Host);
    ViewF GrainUnitVector = Kokkos::create_mirror_view_and_copy(device_memory_space(), GrainUnitVector_Host);

    // View for storing calculated critical diagonal lengths
    ViewF CritDiagonalLength(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength"), 26 * LocalDomainSize);

    // Expected results of calculations
    std::vector<float> CritDiagonalLength_Expected{
        2.732952,  3.403329,  2.062575, 3.708569,  1.7573346, 3.038192, 4.378946,  1.7094444, 3.2407162,
        2.646351,  1.5504522, 3.114523, 3.121724,  3.2783692, 0.177466, 3.1648562, 1.472129,  2.497145,
        3.2025092, 1.914978,  3.057610, 1.9366806, 3.639777,  3.351627, 3.3160222, 1.4418885, 1.715195,
        1.914002,  1.927272,  2.291471, 1.851513,  2.346452,  2.236490, 2.902006,  2.613426,  1.576758,
        1.576758,  2.248777,  2.266173, 2.266173,  2.248777,  1.527831, 1.527831,  1.715195,  1.927272,
        1.914002,  1.851513,  2.291471, 2.902006,  2.613426,  2.346452, 2.236490};

    // For first grain, calculate critical diagonal lenghts
    calcCritDiagonalLength(0, 31.5, 3.5, 1.5, DOCenter(0), DOCenter(1), DOCenter(2), NeighborX, NeighborY, NeighborZ, 0,
                           GrainUnitVector, CritDiagonalLength);
    // For second grain, calculate critical diagonal lenghts
    calcCritDiagonalLength(1, 110.5, 60.5, 0.5, DOCenter(3), DOCenter(4), DOCenter(5), NeighborX, NeighborY, NeighborZ,
                           1, GrainUnitVector, CritDiagonalLength);

    // Copy calculated critical diagonal lengths to the host to check against expected values
    ViewF_H CritDiagonalLength_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CritDiagonalLength);

    // Check results
    for (int i = 0; i < 26 * LocalDomainSize; i++) {
        EXPECT_FLOAT_EQ(CritDiagonalLength_Host(i), CritDiagonalLength_Expected[i]);
    }
}

void testcreateNewOctahedron() {

    // Create a 18-cell domain with origin at (10, 10, 0)
    int nx = 3;
    int MyYSlices = 2;
    int nz = 3;
    int LocalDomainSize = nx * MyYSlices * nz;
    int MyYOffset = 10;
    int ZBound_Low = 0;

    // Create diagonal length and octahedron center data structures on device
    ViewF DiagonalLength(Kokkos::ViewAllocateWithoutInitializing("DiagonalLength"), LocalDomainSize);
    ViewF DOCenter(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * LocalDomainSize);

    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                int D3D1ConvPosition = k * nx * MyYSlices + i * MyYSlices + j;
                int GlobalX = i;
                int GlobalY = j + MyYOffset;
                int GlobalZ = k + ZBound_Low;
                createNewOctahedron(D3D1ConvPosition, DiagonalLength, DOCenter, GlobalX, GlobalY, GlobalZ);
            }
        }
    }

    // Copy back to host and check values
    ViewF_H DiagonalLength_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), DiagonalLength);
    ViewF_H DOCenter_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), DOCenter);
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                int D3D1ConvPosition = k * nx * MyYSlices + i * MyYSlices + j;
                EXPECT_FLOAT_EQ(DiagonalLength_Host(D3D1ConvPosition), 0.01);
                EXPECT_FLOAT_EQ(DOCenter_Host(3 * D3D1ConvPosition), i + 0.5);
                EXPECT_FLOAT_EQ(DOCenter_Host(3 * D3D1ConvPosition + 1), j + MyYOffset + 0.5);
                EXPECT_FLOAT_EQ(DOCenter_Host(3 * D3D1ConvPosition + 2), k + ZBound_Low + 0.5);
            }
        }
    }
}

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, cell_update_tests) {
    testNucleation();
    testFillSteeringVector_Remelt();
    testcalcCritDiagonalLength();
    testcreateNewOctahedron();
}

} // end namespace Test
