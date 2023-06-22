// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAfunctions.hpp"
#include "CAinitialize.hpp"
#include "CAnucleation.hpp"
#include "CAparsefiles.hpp"
#include "CAtypes.hpp"

#include <gtest/gtest.h>

#include "mpi.h"

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
// grain_init_tests
//---------------------------------------------------------------------------//
void testSubstrateInit_ConstrainedGrowth() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // Create test data
    int nz = 2;
    int nzActive = 2;
    int nx = 1;
    // Domain size in Y depends on the number of ranks - each rank has 4 cells in Y
    // Each rank is assigned a different portion of the domain in Y
    int ny = 4 * np;
    int MyYSlices = 4;
    int MyYOffset = 4 * id;
    int LocalActiveDomainSize = nx * MyYSlices * nzActive;
    int LocalDomainSize = nx * MyYSlices * nz;

    double FractSurfaceSitesActive = 0.5; // Each rank will have 2 active cells each, on average
    double RNGSeed = 0.0;
    std::string GrainOrientationFile = checkFileInstalled("GrainOrientationVectors.csv", id);
    int NGrainOrientations = 10000; // Number of grain orientations considered in the simulation
    ViewF GrainUnitVector(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), 9 * NGrainOrientations);
    OrientationInit(id, NGrainOrientations, GrainUnitVector, GrainOrientationFile);

    // Initialize neighbor lists
    NList NeighborX, NeighborY, NeighborZ;
    NeighborListInit(NeighborX, NeighborY, NeighborZ);

    // Initialize views - set initial GrainID values to 0, all CellType values to liquid
    ViewI CellType(Kokkos::ViewAllocateWithoutInitializing("CellType"), LocalDomainSize);
    Kokkos::deep_copy(CellType, Liquid);
    ViewI GrainID("GrainID", LocalDomainSize);
    ViewF DiagonalLength(Kokkos::ViewAllocateWithoutInitializing("DiagonalLength"), LocalActiveDomainSize);
    ViewF DOCenter(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * LocalActiveDomainSize);
    ViewF CritDiagonalLength(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength"), 26 * LocalActiveDomainSize);

    SubstrateInit_ConstrainedGrowth(id, FractSurfaceSitesActive, MyYSlices, nx, ny, MyYOffset, NeighborX, NeighborY,
                                    NeighborZ, GrainUnitVector, NGrainOrientations, CellType, GrainID, DiagonalLength,
                                    DOCenter, CritDiagonalLength, RNGSeed);

    // Copy CellType, GrainID views to host to check values
    ViewI_H CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
    ViewI_H GrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);
    for (int i = 0; i < LocalDomainSize; i++) {
        if (i >= nx * MyYSlices) {
            // Not at bottom surface - should be liquid cells with GrainID still equal to 0
            EXPECT_EQ(GrainID_Host(i), 0);
            EXPECT_EQ(CellType_Host(i), Liquid);
        }
        else {
            // Check that active cells have GrainIDs > 0, and less than 2 * np + 1 (there are 2 * np different positive
            // GrainIDs used for epitaxial grain seeds)
            if (CellType_Host(i) == Active) {
                EXPECT_GT(GrainID_Host(i), 0);
                EXPECT_LT(GrainID_Host(i), 2 * np + 1);
            }
            else {
                // Liquid cells should still have GrainID = 0
                EXPECT_EQ(GrainID_Host(i), 0);
            }
        }
    }
}

void testBaseplateInit_FromGrainSpacing(bool PowderFirstLayer) {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Create test data
    int nz = 4;
    int nzActive = 3;
    double ZMinLayer[1];
    double ZMaxLayer[1];
    ZMinLayer[0] = 0;
    double ZMin = ZMinLayer[0];
    ZMaxLayer[0] = 2 * pow(10, -6);
    int nx = 3;
    // Each rank is assigned a different portion of the domain in Y
    int ny = 3 * np;
    int MyYSlices = 3;
    int MyYOffset = 3 * id;
    double deltax = 1 * pow(10, -6);
    int BaseplateSize;
    int LayerHeight = 2;
    if (PowderFirstLayer)
        BaseplateSize = nx * MyYSlices * (round((ZMaxLayer[0] - ZMin) / deltax) + 1 - LayerHeight);
    else
        BaseplateSize = nx * MyYSlices * (round((ZMaxLayer[0] - ZMin) / deltax) + 1);
    int LocalDomainSize = nx * MyYSlices * nz;
    // There are 36 * np total cells in this domain (nx * ny * nz)
    // Each rank has 36 cells - the bottom 27 cells are assigned baseplate Grain ID values, unless the powder layer of
    // height 2 is used, in which case only the bottom 9 cells should be assigned baseplate Grain ID values The top
    // cells (Z = 3) are outside the "active" portion of the domain and are not assigned Grain IDs with the rest of the
    // baseplate This grain spacing ensures that there will be 1 grain per number of MPI ranks present (larger when
    // powder layer is present as the baseplate will only have a third as many cells)
    double SubstrateGrainSpacing;
    if (PowderFirstLayer)
        SubstrateGrainSpacing = 2.08;
    else
        SubstrateGrainSpacing = 3.0;
    double RNGSeed = 0.0;
    int NextLayer_FirstEpitaxialGrainID;
    // Initialize GrainIDs to 0 on device
    ViewI GrainID("GrainID_Device", LocalDomainSize);

    BaseplateInit_FromGrainSpacing(SubstrateGrainSpacing, nx, ny, ZMin, ZMaxLayer, MyYSlices, MyYOffset, id, deltax,
                                   GrainID, RNGSeed, NextLayer_FirstEpitaxialGrainID, nz, false, PowderFirstLayer,
                                   LayerHeight);

    // Copy results back to host to check
    ViewI_H GrainID_H = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);

    // Check the results for baseplate grains - cells should have GrainIDs between 1 and np (inclusive) if they are part
    // of the active domain, or 0 (unassigned) if not part of the active domain
    // In the case where there is powder on the first layer, the baseplate size is "LayerHeight" fewer cells in the Z
    // direction
    if (PowderFirstLayer)
        EXPECT_EQ(BaseplateSize, nx * MyYSlices * (nzActive - LayerHeight));
    else
        EXPECT_EQ(BaseplateSize, nx * MyYSlices * nzActive);
    for (int i = 0; i < BaseplateSize; i++) {
        EXPECT_GT(GrainID_H(i), 0);
        EXPECT_LT(GrainID_H(i), np + 1);
    }
    for (int i = BaseplateSize; i < LocalDomainSize; i++) {
        EXPECT_EQ(GrainID_H(i), 0);
    }
    // Next unused GrainID should be the number of grains present in the baseplate plus 2 (since GrainID = 0 is not used
    // for any baseplate grains)
    EXPECT_EQ(NextLayer_FirstEpitaxialGrainID, np + 2);
}

void testPowderInit() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // Create test data
    int nz = 6;
    int layernumber = 1;
    // Z = 4 and Z = 5 should be the region seeded with powder layer Grain IDs
    double ZMaxLayer[2];
    ZMaxLayer[1] = 5.0 * pow(10, -6);
    double deltax = 1.0 * pow(10, -6);
    double ZMin = 0;
    int LayerHeight = 2;
    int nx = 1;
    // Each rank is assigned a different portion of the domain in Y
    int ny = 2 * np;
    int MyYSlices = 2;
    int MyYOffset = 2 * id;
    int LocalDomainSize = nx * MyYSlices * nz;

    // Initialize GrainIDs to 0 on device
    ViewI GrainID("GrainID_Device", LocalDomainSize);

    // Starting Grain ID for powder grains
    int NextLayer_FirstEpitaxialGrainID = 10;
    int PreviousLayer_FirstEpitaxialGrainID = NextLayer_FirstEpitaxialGrainID;

    // Seed used to shuffle powder layer grain IDs
    double RNGSeed = 0.0;

    PowderInit(layernumber, nx, ny, LayerHeight, ZMaxLayer, ZMin, deltax, MyYSlices, MyYOffset, id, GrainID, RNGSeed,
               NextLayer_FirstEpitaxialGrainID, 1.0);

    // Copy results back to host to check
    ViewI_H GrainID_H = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);

    // Check the results - were the right number of powder grain IDs used?
    int Expected_NextLayer_FirstEpitaxialGrainID = PreviousLayer_FirstEpitaxialGrainID + nx * ny * LayerHeight;
    EXPECT_EQ(NextLayer_FirstEpitaxialGrainID, Expected_NextLayer_FirstEpitaxialGrainID);

    // Check the results - powder grains should have unique Grain ID values larger than
    // PreviousLayer_FirstEpitaxialGrainID-1 and smaller than NextLayer_FirstEpitaxialGrainID. Other cells should still
    // have Grain IDs of 0
    int TopPowderLayer = nx * MyYSlices * (nz - LayerHeight);
    for (int i = 0; i < TopPowderLayer; i++) {
        EXPECT_EQ(GrainID_H(i), 0);
    }
    for (int i = TopPowderLayer; i < LocalDomainSize; i++) {
        EXPECT_GT(GrainID_H(i), PreviousLayer_FirstEpitaxialGrainID - 1);
        EXPECT_LT(GrainID_H(i), NextLayer_FirstEpitaxialGrainID);
    }
}
//---------------------------------------------------------------------------//
// cell_init_tests
//---------------------------------------------------------------------------//
void testCellTypeInit() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // Domain for each rank
    int nx = id + 5;
    int MyYSlices = id + 5;
    int nzActive = 5;
    int ZBound_Low = 2;
    int nz = nzActive + ZBound_Low;
    int LocalActiveDomainSize = nx * MyYSlices * nzActive;
    int LocalDomainSize = nx * MyYSlices * nz;

    // Temporary host view for initializing NumberOfSolidificationEvents
    ViewI_H NumberOfSolidificationEvents_Host(
        Kokkos::ViewAllocateWithoutInitializing("NumberOfSolidificationEvents_Host"), LocalActiveDomainSize);
    for (int RankZ = 0; RankZ < nzActive; RankZ++) {
        for (int RankX = 0; RankX < nx; RankX++) {
            for (int RankY = 0; RankY < MyYSlices; RankY++) {
                int D3D1ConvPosition = RankZ * nx * MyYSlices + RankX * MyYSlices + RankY;
                // Assign some of these a value of 0 (these will be solid cells), and others a positive value
                // (these will be tempsolid cells)
                if (RankX + RankY % 2 == 0)
                    NumberOfSolidificationEvents_Host(D3D1ConvPosition) = 0;
                else
                    NumberOfSolidificationEvents_Host(D3D1ConvPosition) = 1;
            }
        }
    }

    ViewI NumberOfSolidificationEvents =
        Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), NumberOfSolidificationEvents_Host);
    // Start with cell type of 0
    ViewI CellType("CellType", LocalDomainSize);

    // Initialize cell type values
    CellTypeInit(nx, MyYSlices, LocalActiveDomainSize, CellType, NumberOfSolidificationEvents, id, ZBound_Low);

    // Copy cell types back to host to check
    ViewI_H CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);

    for (int RankZ = ZBound_Low; RankZ < nzActive; RankZ++) {
        for (int RankX = 0; RankX < nx; RankX++) {
            for (int RankY = 0; RankY < MyYSlices; RankY++) {
                int D3D1ConvPosition = RankZ * nx * MyYSlices + RankX * MyYSlices + RankY;
                // Cells with no associated solidification events should be solid, others TempSolid
                if (RankX + RankY % 2 == 0)
                    EXPECT_EQ(CellType_Host(D3D1ConvPosition), Solid);
                else
                    EXPECT_EQ(CellType_Host(D3D1ConvPosition), TempSolid);
            }
        }
    }
}
//---------------------------------------------------------------------------//
// nuclei_init_tests
//---------------------------------------------------------------------------//
void testNucleiInit() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Create test data
    // Only Z = 1 through 4 is part of this layer
    int nzActive = 4;
    int ZBound_Low = 1;
    int nx = 4;
    // Each rank is assigned a different portion of the domain in Y
    int ny = 2 * np;
    int MyYSlices = 2;
    int MyYOffset = 2 * id;
    double deltax = 1;
    int LocalActiveDomainSize = nx * MyYSlices * nzActive;
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

    // There are 40 * np total cells in this domain (nx * ny * nz)
    // Each rank has 40 cells - the top 32 cells are part of the active layer and are candidates for nucleation
    // assignment
    double NMax = 0.125; // This nucleation density ensures there will be 4 potential nuclei per MPI rank present
                         // without remelting (each cell solidifies once)
    int MaxPotentialNuclei_PerPass = 4 * np;
    // A cell can solidify 1-3 times
    int MaxSolidificationEvents_Count = 3;
    double dTN = 1;
    double dTsigma = 0.0001;
    double RNGSeed = 0.0;

    // Initialize MaxSolidificationEvents to 3 for each layer. LayerTimeTempHistory and NumberOfSolidificationEvents are
    // initialized for each cell on the host and copied to the device
    ViewI_H MaxSolidificationEvents_Host(Kokkos::ViewAllocateWithoutInitializing("MaxSolidificationEvents_Host"), 2);
    ViewI_H NumberOfSolidificationEvents_Host(
        Kokkos::ViewAllocateWithoutInitializing("NumberOfSolidificationEvents_Host"), LocalActiveDomainSize);
    ViewF3D LayerTimeTempHistory_Host(Kokkos::ViewAllocateWithoutInitializing("LayerTimeTempHistory_Host"),
                                      LocalActiveDomainSize, MaxSolidificationEvents_Count, 3);
    MaxSolidificationEvents_Host(0) = MaxSolidificationEvents_Count;
    MaxSolidificationEvents_Host(1) = MaxSolidificationEvents_Count;
    // Cells solidify 1, 2, or 3 times, depending on their X coordinate
    for (int k = 0; k < nzActive; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                int D3D1ConvPosition = k * nx * MyYSlices + i * MyYSlices + j;
                if (i < nx / 2 - 1)
                    NumberOfSolidificationEvents_Host(D3D1ConvPosition) = 3;
                else if (i < nx / 2)
                    NumberOfSolidificationEvents_Host(D3D1ConvPosition) = 2;
                else
                    NumberOfSolidificationEvents_Host(D3D1ConvPosition) = 1;
            }
        }
    }
    for (int n = 0; n < MaxSolidificationEvents_Count; n++) {
        for (int RankZ = 0; RankZ < nzActive; RankZ++) {
            for (int RankX = 0; RankX < nx; RankX++) {
                for (int RankY = 0; RankY < MyYSlices; RankY++) {
                    int D3D1ConvPosition = RankZ * nx * MyYSlices + RankX * MyYSlices + RankY;
                    int GlobalZ = RankZ + ZBound_Low;
                    if (n < NumberOfSolidificationEvents_Host(D3D1ConvPosition)) {
                        LayerTimeTempHistory_Host(D3D1ConvPosition, n, 0) =
                            GlobalZ + RankY + MyYOffset +
                            (LocalActiveDomainSize * n); // melting time step depends on solidification event number
                        LayerTimeTempHistory_Host(D3D1ConvPosition, n, 1) =
                            GlobalZ + RankY + MyYOffset + 1 +
                            (LocalActiveDomainSize * n); // liquidus time stemp depends on solidification event number
                        LayerTimeTempHistory_Host(D3D1ConvPosition, n, 2) =
                            1.2; // ensures that a cell's nucleation time will be 1 time step after its CritTimeStep
                                 // value
                    }
                }
            }
        }
    }

    using memory_space = TEST_MEMSPACE;
    ViewI MaxSolidificationEvents = Kokkos::create_mirror_view_and_copy(memory_space(), MaxSolidificationEvents_Host);
    ViewI NumberOfSolidificationEvents =
        Kokkos::create_mirror_view_and_copy(memory_space(), NumberOfSolidificationEvents_Host);
    ViewF3D LayerTimeTempHistory = Kokkos::create_mirror_view_and_copy(memory_space(), LayerTimeTempHistory_Host);

    // Nucleation data structure, containing views of nuclei locations, time steps, and ids, and nucleation event
    // counters - initialized with an estimate on the number of nuclei in the layer Without knowing
    // PossibleNuclei_ThisRankThisLayer yet, initialize nucleation data structures to estimated sizes, resize inside of
    // NucleiInit when the number of nuclei per rank is known
    int EstimatedNuclei_ThisRankThisLayer = NMax * pow(deltax, 3) * LocalActiveDomainSize;
    Nucleation<memory_space> nucleation(
        EstimatedNuclei_ThisRankThisLayer, NMax, deltax,
        100); // NucleiGrainID should start at -101 - supply optional input arg to constructor

    // Fill in nucleation data structures, and assign nucleation undercooling values to potential nucleation events
    // Potential nucleation grains are only associated with liquid cells in layer 1 - they will be initialized for each
    // successive layer when layer 1 in complete
    nucleation.placeNuclei(MaxSolidificationEvents, NumberOfSolidificationEvents, LayerTimeTempHistory, RNGSeed, 1, nx,
                           ny, nzActive, dTN, dTsigma, MyYSlices, MyYOffset, ZBound_Low, id, AtNorthBoundary,
                           AtSouthBoundary);

    // Copy results back to host to check
    ViewI_H NucleiLocation_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), nucleation.NucleiLocations);
    ViewI_H NucleiGrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), nucleation.NucleiGrainID);

    // Was the nucleation counter initialized to zero?
    EXPECT_EQ(nucleation.NucleationCounter, 0);

    // Is the total number of nuclei in the system correct, based on the number of remelting events? Equal probability
    // of creating a nucleus each time a cell resolidifies
    int ExpectedNucleiPerRank = 100 + MaxSolidificationEvents_Count * MaxPotentialNuclei_PerPass;
    EXPECT_EQ(nucleation.Nuclei_WholeDomain, ExpectedNucleiPerRank);
    for (int n = 0; n < nucleation.PossibleNuclei; n++) {
        // Are the nuclei grain IDs negative numbers in the expected range based on the inputs?
        EXPECT_GT(NucleiGrainID_Host(n), -(100 + ExpectedNucleiPerRank * np + 1));
        EXPECT_LT(NucleiGrainID_Host(n), -100);
        // Are the correct undercooling values associated with the correct cell locations?
        // Cell location is a global position (relative to the bottom of the whole domain, not the layer)
        int GlobalCellLocation = NucleiLocation_Host(n);
        int GlobalZ = GlobalCellLocation / (nx * MyYSlices);
        int Rem = GlobalCellLocation % (nx * MyYSlices);
        int RankY = Rem % MyYSlices;
        // Expected nucleation time with remelting can be one of 3 possibilities, depending on the associated
        // solidification event
        int Expected_NucleationTimeNoRM = GlobalZ + RankY + MyYOffset + 2;
        int AssociatedSEvent = nucleation.NucleationTimes_Host(n) / LocalActiveDomainSize;
        int Expected_NucleationTimeRM = Expected_NucleationTimeNoRM + AssociatedSEvent * LocalActiveDomainSize;
        EXPECT_EQ(nucleation.NucleationTimes_Host(n), Expected_NucleationTimeRM);

        // Are the nucleation events in order of the time steps at which they may occur?
        if (n < nucleation.PossibleNuclei - 2) {
            EXPECT_LE(nucleation.NucleationTimes_Host(n), nucleation.NucleationTimes_Host(n + 1));
        }
    }
}

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, grain_init_tests) {
    testSubstrateInit_ConstrainedGrowth();
    // w/ and w/o space left for powder layer
    testBaseplateInit_FromGrainSpacing(true);
    testBaseplateInit_FromGrainSpacing(false);
    testPowderInit();
}
TEST(TEST_CATEGORY, cell_init_test) { testCellTypeInit(); }
TEST(TEST_CATEGORY, nuclei_init_test) { testNucleiInit(); }
} // end namespace Test
