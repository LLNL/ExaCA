// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAcelldata.hpp"
#include "CAfunctions.hpp"
#include "CAinitialize.hpp"
#include "CAnucleation.hpp"
#include "CAtypes.hpp"
#include "CAupdate.hpp"

#include <gtest/gtest.h>

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
void testNucleation() {

    using memory_space = TEST_MEMSPACE;
    using view_int = Kokkos::View<int *, TEST_MEMSPACE>;
    using view_int_host = typename view_int::HostMirror;

    // Create views - 125 cells, 75 of which are part of the active layer of the domain (Z = 2-5)
    int nx = 5;
    int MyYSlices = 5;
    int nz = 5;
    int ZBound_Low = 2;
    int nzActive = 3;
    int LocalActiveDomainSize = nx * MyYSlices * nzActive;
    int LocalDomainSize = nx * MyYSlices * nz;

    // All cells have GrainID of 0, CellType of Liquid - with the exception of the locations where the nucleation events
    // are unable to occur
    CellData<TEST_MEMSPACE> cellData(LocalDomainSize, LocalActiveDomainSize, nx, MyYSlices, ZBound_Low);
    Kokkos::deep_copy(cellData.CellType_AllLayers, Liquid);
    auto CellType = cellData.getCellTypeSubview();
    ViewI_H CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
    ViewI GrainID = cellData.getGrainIDSubview();
    ViewI_H GrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);

    // Create test nucleation data - 10 possible events
    int PossibleNuclei = 10;
    Nucleation<memory_space> nucleation(PossibleNuclei, 1.0, 1.0);
    nucleation.PossibleNuclei = 10;
    view_int_host NucleiLocations_Host(Kokkos::ViewAllocateWithoutInitializing("NucleiLocations_Host"), PossibleNuclei);
    view_int_host NucleiGrainID_Host(Kokkos::ViewAllocateWithoutInitializing("NucleiGrainID_Host"), PossibleNuclei);
    for (int n = 0; n < PossibleNuclei; n++) {
        // NucleationTimes values should be in the order in which the events occur - start with setting them between 0
        // and 9 NucleiLocations are in order starting with 0, through 9 (locations relative to the bottom of the
        // layer)
        nucleation.NucleationTimes_Host(n) = n;
        NucleiLocations_Host(n) = n;
        // Give these nucleation events grain IDs based on their order, starting with -1 and counting down
        NucleiGrainID_Host(n) = -(n + 1);
    }
    // Include the case where 2 potential nucleation events (3 and 4) happen on the same time step - both successful
    // Let nucleation events 3 and 4 both occur on time step 4
    nucleation.NucleationTimes_Host(3) = nucleation.NucleationTimes_Host(4);

    // Include the case where a potential nucleation event (2) is unsuccessful (time step 2)
    int UnsuccessfulLocA = 2;
    CellType_Host(UnsuccessfulLocA) = Active;
    GrainID_Host(UnsuccessfulLocA) = 1;

    // Include the case where 2 potential nucleation events (5 and 6) happen on the same time step (time step 6) - the
    // first successful, the second unsuccessful
    nucleation.NucleationTimes_Host(5) = nucleation.NucleationTimes_Host(6);
    int UnsuccessfulLocC = 6;
    CellType_Host(UnsuccessfulLocC) = Active;
    GrainID_Host(UnsuccessfulLocC) = 2;

    // Include the case where 2 potential nucleation events (8 and 9) happen on the same time step (time step 8) - the
    // first unsuccessful, the second successful
    nucleation.NucleationTimes_Host(8) = nucleation.NucleationTimes_Host(9);
    int UnsuccessfulLocB = 8;
    CellType_Host(UnsuccessfulLocB) = Active;
    GrainID_Host(UnsuccessfulLocB) = 3;

    // Copy host views to device
    CellType = Kokkos::create_mirror_view_and_copy(TEST_MEMSPACE(), CellType_Host);
    GrainID = Kokkos::create_mirror_view_and_copy(TEST_MEMSPACE(), GrainID_Host);
    nucleation.NucleiLocations = Kokkos::create_mirror_view_and_copy(TEST_MEMSPACE(), NucleiLocations_Host);
    nucleation.NucleiGrainID = Kokkos::create_mirror_view_and_copy(TEST_MEMSPACE(), NucleiGrainID_Host);

    // Steering Vector
    view_int SteeringVector(Kokkos::ViewAllocateWithoutInitializing("SteeringVector"), LocalActiveDomainSize);
    // Initialize steering vector size to 0
    view_int numSteer("SteeringVector", 1);

    // Take enough time steps such that every nucleation event has a chance to occur
    for (int cycle = 0; cycle < 10; cycle++) {
        nucleation.nucleate_grain(cycle, cellData, ZBound_Low, nx, MyYSlices, SteeringVector, numSteer);
    }

    // Copy CellType, SteeringVector, numSteer, GrainID back to host to check nucleation results
    CellType = cellData.getCellTypeSubview();
    CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
    auto SteeringVector_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), SteeringVector);
    auto numSteer_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), numSteer);
    GrainID = cellData.getGrainIDSubview();
    GrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);

    // Check that all 10 possible nucleation events were attempted
    EXPECT_EQ(nucleation.NucleationCounter, 10);
    // Check that 7 of the 10 nucleation events were successful
    EXPECT_EQ(nucleation.SuccessfulNucleationCounter, 7);
    EXPECT_EQ(numSteer_Host(0), 7);

    // Ensure that the 3 events that should not have occurred, did not occur
    // These cells should be untouched - active type, and same GrainID they was initialized with
    EXPECT_EQ(CellType_Host(2), Active);
    EXPECT_EQ(GrainID_Host(2), 1);
    EXPECT_EQ(CellType_Host(6), Active);
    EXPECT_EQ(GrainID_Host(6), 2);
    EXPECT_EQ(CellType_Host(8), Active);
    EXPECT_EQ(GrainID_Host(8), 3);

    // Check that the successful nucleation events occurred as expected
    // For each cell location (relative to the layer bottom) that should be home to a successful nucleation event,
    // check that the CellType has been set to FutureActive and the GrainID matches the expected value Also check that
    // the adjusted cell coordinate (relative to the current layer bounds) appears somewhere within the steering vector
    std::vector<int> SuccessfulNuc_GrainIDs{-1, -2, -4, -5, -6, -8, -10};
    std::vector<int> SuccessfulNuc_CellLocations{0, 1, 3, 4, 5, 7, 9};
    for (int nevent = 0; nevent < 7; nevent++) {
        int CellLocation = SuccessfulNuc_CellLocations[nevent];
        EXPECT_EQ(CellType_Host(CellLocation), FutureActive);
        EXPECT_EQ(GrainID_Host(CellLocation), SuccessfulNuc_GrainIDs[nevent]);
        bool OnSteeringVector = false;
        for (int svloc = 0; svloc < 7; svloc++) {
            if (SteeringVector_Host(svloc) == CellLocation) {
                OnSteeringVector = true;
                break;
            }
        }
        EXPECT_TRUE(OnSteeringVector);
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

    // Initialize neighbor lists
    NList NeighborX, NeighborY, NeighborZ;
    NeighborListInit(NeighborX, NeighborY, NeighborZ);

    CellData<TEST_MEMSPACE> cellData(LocalDomainSize, LocalActiveDomainSize, nx, MyYSlices, ZBound_Low);
    ViewI_H MeltTimeStep_Host(Kokkos::ViewAllocateWithoutInitializing("MeltTimeStep_Host"), LocalDomainSize);
    ViewI_H CritTimeStep_Host(Kokkos::ViewAllocateWithoutInitializing("CritTimeStep_Host"), LocalDomainSize);
    ViewI_H SolidificationEventCounter_Host("SolidificationEventCounter_Host", LocalDomainSize); // init to 0
    ViewI_H NumberOfSolidificationEvents_Host(
        Kokkos::ViewAllocateWithoutInitializing("NumberOfSolidificationEvents_Host"), LocalDomainSize);
    ViewF_H UndercoolingChange_Host(Kokkos::ViewAllocateWithoutInitializing("UndercoolingChange_Host"),
                                    LocalDomainSize);
    ViewF_H UndercoolingCurrent_Host("UndercoolingCurrent_Host",
                                     LocalDomainSize); // initialize to 0, no initial undercooling
    ViewF3D_H LayerTimeTempHistory_Host("LayerTimeTempHistory_Host", LocalActiveDomainSize, 1,
                                        3); // initialize to 0, no initial undercooling
    ViewI CellType = cellData.getCellTypeSubview();
    ViewI_H CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
    ViewI GrainID = cellData.getGrainIDSubview();
    ViewI_H GrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);

    for (int D3D1ConvPosition = 0; D3D1ConvPosition < LocalActiveDomainSize; D3D1ConvPosition++) {
        GrainID_Host(D3D1ConvPosition) = 1;
        CellType_Host(D3D1ConvPosition) = TempSolid;
    }
    for (int GlobalD3D1ConvPosition = 0; GlobalD3D1ConvPosition < LocalDomainSize; GlobalD3D1ConvPosition++) {
        // Cell coordinates on this rank in X, Y, and Z (GlobalZ = relative to domain bottom)
        int GlobalZ = GlobalD3D1ConvPosition / (nx * MyYSlices);
        int Rem = GlobalD3D1ConvPosition % (nx * MyYSlices);
        int RankY = Rem % MyYSlices;
        // Cells at Z = 0 through Z = 2 are Solid, Z = 3 and 4 are TempSolid
        if (GlobalZ <= 2) {
            // Solid cells should have -1 assigned as their melt/crit time steps
            MeltTimeStep_Host(GlobalD3D1ConvPosition) = -1;
            CritTimeStep_Host(GlobalD3D1ConvPosition) = -1;
        }
        else {
            // Cells "melt" at a time step corresponding to their Y location in the overall domain (depends on MyYOffset
            // of the rank)
            MeltTimeStep_Host(GlobalD3D1ConvPosition) = RankY + MyYOffset + 1;
            // Cells reach liquidus during cooling 2 time steps after melting
            CritTimeStep_Host(GlobalD3D1ConvPosition) = MeltTimeStep_Host(GlobalD3D1ConvPosition) + 2;
        }
        UndercoolingChange_Host(GlobalD3D1ConvPosition) = 0.2;
    }
    for (int D3D1ConvPosition = 0; D3D1ConvPosition < LocalActiveDomainSize; D3D1ConvPosition++) {
        // Cell coordinates on this rank in X, Y, and Z (GlobalZ = relative to domain bottom)
        int RankZ = D3D1ConvPosition / (nx * MyYSlices);
        int Rem = D3D1ConvPosition % (nx * MyYSlices);
        int RankX = Rem / MyYSlices;
        int RankY = Rem % MyYSlices;
        int GlobalZ = RankZ + ZBound_Low;
        int GlobalD3D1ConvPosition = GlobalZ * nx * MyYSlices + RankX * MyYSlices + RankY;
        LayerTimeTempHistory_Host(D3D1ConvPosition, 0, 0) = MeltTimeStep_Host(GlobalD3D1ConvPosition);
        LayerTimeTempHistory_Host(D3D1ConvPosition, 0, 1) = CritTimeStep_Host(GlobalD3D1ConvPosition);
        LayerTimeTempHistory_Host(D3D1ConvPosition, 0, 2) = UndercoolingChange_Host(GlobalD3D1ConvPosition);
        NumberOfSolidificationEvents_Host(D3D1ConvPosition) = 1;
    }

    // Steering Vector
    ViewI SteeringVector(Kokkos::ViewAllocateWithoutInitializing("SteeringVector"), LocalActiveDomainSize);
    ViewI_H numSteer_Host(Kokkos::ViewAllocateWithoutInitializing("SteeringVectorSize"), 1);
    numSteer_Host(0) = 0;
    // Copy views to device for test
    ViewI numSteer = Kokkos::create_mirror_view_and_copy(device_memory_space(), numSteer_Host);
    GrainID = Kokkos::create_mirror_view_and_copy(device_memory_space(), GrainID_Host);
    CellType = Kokkos::create_mirror_view_and_copy(device_memory_space(), CellType_Host);
    ViewI MeltTimeStep = Kokkos::create_mirror_view_and_copy(device_memory_space(), MeltTimeStep_Host);
    ViewI CritTimeStep = Kokkos::create_mirror_view_and_copy(device_memory_space(), CritTimeStep_Host);
    ViewI SolidificationEventCounter =
        Kokkos::create_mirror_view_and_copy(device_memory_space(), SolidificationEventCounter_Host);
    ViewI NumberOfSolidificationEvents =
        Kokkos::create_mirror_view_and_copy(device_memory_space(), NumberOfSolidificationEvents_Host);
    ViewF UndercoolingChange = Kokkos::create_mirror_view_and_copy(device_memory_space(), UndercoolingChange_Host);
    ViewF UndercoolingCurrent = Kokkos::create_mirror_view_and_copy(device_memory_space(), UndercoolingCurrent_Host);
    ViewF3D LayerTimeTempHistory =
        Kokkos::create_mirror_view_and_copy(device_memory_space(), LayerTimeTempHistory_Host);

    int numcycles = 15;
    for (int cycle = 1; cycle <= numcycles; cycle++) {
        // Update cell types, local undercooling each time step, and fill the steering vector
        std::cout << cycle << std::endl;
        FillSteeringVector_Remelt(cycle, LocalActiveDomainSize, nx, MyYSlices, NeighborX, NeighborY, NeighborZ,
                                  CritTimeStep, UndercoolingCurrent, UndercoolingChange, cellData, ZBound_Low, nzActive,
                                  SteeringVector, numSteer, numSteer_Host, MeltTimeStep, SolidificationEventCounter,
                                  NumberOfSolidificationEvents, LayerTimeTempHistory);
    }

    // Copy CellType, SteeringVector, numSteer, UndercoolingCurrent, Buffers back to host to check steering vector
    // construction results
    CellType = cellData.getCellTypeSubview();
    CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
    ViewI_H SteeringVector_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), SteeringVector);
    numSteer_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), numSteer);
    UndercoolingCurrent_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), UndercoolingCurrent);

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
    for (int D3D1ConvPosition = 0; D3D1ConvPosition < LocalActiveDomainSize; D3D1ConvPosition++) {
        int RankZ = D3D1ConvPosition / (nx * MyYSlices);
        int GlobalZ = RankZ + ZBound_Low;
        int GlobalD3D1ConvPosition = D3D1ConvPosition + ZBound_Low * nx * MyYSlices;
        if (GlobalZ <= 2)
            EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(GlobalD3D1ConvPosition), 0.0);
        else {
            if (numcycles < MeltTimeStep_Host(GlobalD3D1ConvPosition)) {
                EXPECT_EQ(CellType_Host(D3D1ConvPosition), TempSolid);
                EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(GlobalD3D1ConvPosition), 0.0);
            }
            else if ((numcycles >= MeltTimeStep_Host(GlobalD3D1ConvPosition)) &&
                     (numcycles <= CritTimeStep_Host(GlobalD3D1ConvPosition))) {
                EXPECT_EQ(CellType_Host(D3D1ConvPosition), Liquid);
                EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(GlobalD3D1ConvPosition), 0.0);
            }
            else {
                EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(GlobalD3D1ConvPosition),
                                (numcycles - CritTimeStep_Host(GlobalD3D1ConvPosition)) * 0.2);
                if (GlobalZ == 4)
                    EXPECT_EQ(CellType_Host(D3D1ConvPosition), Liquid);
                else {
                    EXPECT_EQ(CellType_Host(D3D1ConvPosition), FutureActive);
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
}

void testcalcCritDiagonalLength() {
    using memory_space = TEST_MEMSPACE;
    using view_type = Kokkos::View<float *, memory_space>;

    // 2 "cells" in the domain
    int LocalDomainSize = 2;
    // First orientation is orientation 12 (starting indexing at 0) of GrainOrientationVectors.csv
    // Second orientation is orientation 28 (starting indexing at 0)
    std::vector<float> GrainUnitVectorV{0.52877,  0.651038, -0.544565, -0.572875, 0.747163,  0.336989,
                                        0.626272, 0.133778, 0.768041,  0.736425,  -0.530983, 0.419208,
                                        0.664512, 0.683971, -0.301012, -0.126894, 0.500241,  0.856538};
    // Load unit vectors into view
    view_type GrainUnitVector(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), 9 * LocalDomainSize);
    auto GrainUnitVector_Host = Kokkos::create_mirror_view(Kokkos::HostSpace(), GrainUnitVector);
    for (int i = 0; i < 9 * LocalDomainSize; i++) {
        GrainUnitVector_Host(i) = GrainUnitVectorV[i];
    }
    // Copy octahedron center and grain unit vector data to the device
    GrainUnitVector = Kokkos::create_mirror_view_and_copy(memory_space(), GrainUnitVector_Host);

    // Initialize neighbor lists
    NList NeighborX, NeighborY, NeighborZ;
    NeighborListInit(NeighborX, NeighborY, NeighborZ);

    // Load octahedron centers into view
    view_type DOCenter(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * LocalDomainSize);
    Kokkos::parallel_for(
        "TestInitCritDiagonalLength", 1, KOKKOS_LAMBDA(const int) {
            // Octahedron center for first cell (does not align with CA cell center, which is at 31.5, 3.5, 1.5)
            DOCenter(0) = 30.611142;
            DOCenter(1) = 3.523741;
            DOCenter(2) = 0.636301;
            // Octahedron center for second cell (aligns with CA cell center at 110.5, 60.5, 0.5))
            DOCenter(3) = 110.5;
            DOCenter(4) = 60.5;
            DOCenter(5) = 0.5;
        });

    // View for storing calculated critical diagonal lengths
    view_type CritDiagonalLength(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength"), 26 * LocalDomainSize);

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
    auto CritDiagonalLength_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CritDiagonalLength);

    // Check results
    for (int i = 0; i < 26 * LocalDomainSize; i++) {
        EXPECT_FLOAT_EQ(CritDiagonalLength_Host(i), CritDiagonalLength_Expected[i]);
    }
}

void testcreateNewOctahedron() {
    using memory_space = TEST_MEMSPACE;
    using view_type = Kokkos::View<float *, memory_space>;

    // Create a 18-cell domain with origin at (10, 10, 0)
    int nx = 3;
    int MyYSlices = 2;
    int nz = 3;
    int LocalDomainSize = nx * MyYSlices * nz;
    int MyYOffset = 10;
    int ZBound_Low = 0;

    // Create diagonal length and octahedron center data structures on device
    view_type DiagonalLength(Kokkos::ViewAllocateWithoutInitializing("DiagonalLength"), LocalDomainSize);
    view_type DOCenter(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * LocalDomainSize);

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
    auto DiagonalLength_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), DiagonalLength);
    auto DOCenter_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), DOCenter);
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
TEST(TEST_CATEGORY, cell_update_tests) {
    testNucleation();
    testFillSteeringVector_Remelt();
    testcalcCritDiagonalLength();
    testcreateNewOctahedron();
    testConvertGrainIDForBuffer();
}

} // end namespace Test
