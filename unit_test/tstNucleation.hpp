// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAcelldata.hpp"
#include "CAfunctions.hpp"
#include "CAinitialize.hpp"
#include "CAinputs.hpp"
#include "CAnucleation.hpp"
#include "CAparsefiles.hpp"
#include "CAtemperature.hpp"
#include "CAtypes.hpp"

#include <gtest/gtest.h>

#include "mpi.h"

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
// tests for Nucleation struct
//---------------------------------------------------------------------------//
// Tests Nucleation object construction and placeNuclei
void testNucleiInit() {

    using memory_space = TEST_MEMSPACE;

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Create test data
    // Only Z = 1 through 4 is part of this layer
    int nz_layer = 4;
    int z_layer_bottom = 1;
    int nx = 4;
    // Each rank is assigned a different portion of the domain in Y
    int ny = 2 * np;
    int ny_local = 2;
    int y_offset = 2 * id;
    double deltax = 1;
    int DomainSize = nx * ny_local * nz_layer;
    int NumberOfLayers = 2;
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
    int MaxPotentialNuclei_PerPass = 4 * np;
    // A cell can solidify 1-3 times
    int MaxSolidificationEvents_Count = 3;
    // default inputs struct with manually set nucleation parameters
    // This nucleation density ensures there will be 4 potential nuclei per MPI rank present
    // without remelting (each cell solidifies once)
    Inputs inputs;
    inputs.nucleation.dTN = 1;
    inputs.nucleation.dTsigma = 0.0001;
    inputs.nucleation.NMax = 0.125;

    // Allocate temperature data structures
    Temperature<memory_space> temperature(DomainSize, NumberOfLayers, inputs.temperature, 1);
    // Resize LayerTimeTempHistory with the known max number of solidification events
    Kokkos::resize(temperature.LayerTimeTempHistory, DomainSize, MaxSolidificationEvents_Count, 3);
    // Initialize MaxSolidificationEvents to 3 for each layer. LayerTimeTempHistory and NumberOfSolidificationEvents are
    // initialized for each cell on the host and copied to the device
    ViewI_H MaxSolidificationEvents_Host(Kokkos::ViewAllocateWithoutInitializing("MaxSolidificationEvents_Host"),
                                         NumberOfLayers);
    MaxSolidificationEvents_Host(0) = MaxSolidificationEvents_Count;
    MaxSolidificationEvents_Host(1) = MaxSolidificationEvents_Count;
    // Cells solidify 1, 2, or 3 times, depending on their X coordinate
    Kokkos::parallel_for(
        "NumSolidificationEventsInit", nz_layer, KOKKOS_LAMBDA(const int &coord_z) {
            for (int coord_x = 0; coord_x < nx; coord_x++) {
                for (int coord_y = 0; coord_y < ny_local; coord_y++) {
                    int index = get1Dindex(coord_x, coord_y, coord_z, nx, ny_local);
                    if (coord_x < nx / 2 - 1)
                        temperature.NumberOfSolidificationEvents(index) = 3;
                    else if (coord_x < nx / 2)
                        temperature.NumberOfSolidificationEvents(index) = 2;
                    else
                        temperature.NumberOfSolidificationEvents(index) = 1;
                }
            }
        });
    Kokkos::fence();
    Kokkos::parallel_for(
        "LayerTimeTempHistoryInit", MaxSolidificationEvents_Count, KOKKOS_LAMBDA(const int &n) {
            for (int coord_z = 0; coord_z < nz_layer; coord_z++) {
                for (int coord_x = 0; coord_x < nx; coord_x++) {
                    for (int coord_y = 0; coord_y < ny_local; coord_y++) {
                        int index = get1Dindex(coord_x, coord_y, coord_z, nx, ny_local);
                        int coord_z_AllLayers = coord_z + z_layer_bottom;
                        if (n < temperature.NumberOfSolidificationEvents(index)) {
                            // melting time step depends on solidification event number
                            temperature.LayerTimeTempHistory(index, n, 0) =
                                coord_z_AllLayers + coord_y + y_offset + (DomainSize * n);
                            // liquidus time stemp depends on solidification event number
                            temperature.LayerTimeTempHistory(index, n, 1) =
                                coord_z_AllLayers + coord_y + y_offset + 1 + (DomainSize * n);
                            // ensures that a cell's nucleation time will be 1 time step after its CritTimeStep value
                            temperature.LayerTimeTempHistory(index, n, 2) = 1.2;
                        }
                    }
                }
            }
        });
    Kokkos::fence();
    temperature.MaxSolidificationEvents =
        Kokkos::create_mirror_view_and_copy(memory_space(), MaxSolidificationEvents_Host);

    // Nucleation data structure, containing views of nuclei locations, time steps, and ids, and nucleation event
    // counters - initialized with an estimate on the number of nuclei in the layer Without knowing
    // PossibleNuclei_ThisRankThisLayer yet, initialize nucleation data structures to estimated sizes, resize inside of
    // NucleiInit when the number of nuclei per rank is known
    int EstimatedNuclei_ThisRankThisLayer = inputs.nucleation.NMax * pow(deltax, 3) * DomainSize;
    Nucleation<memory_space> nucleation(
        EstimatedNuclei_ThisRankThisLayer, deltax, inputs.nucleation,
        100); // NucleiGrainID should start at -101 - supply optional input arg to constructor
    // Ensure nucleation inputs in nucleation struct were correctly initialized
    EXPECT_DOUBLE_EQ(inputs.nucleation.NMax, nucleation._inputs.NMax);
    EXPECT_DOUBLE_EQ(inputs.nucleation.dTN, nucleation._inputs.dTN);
    EXPECT_DOUBLE_EQ(inputs.nucleation.dTsigma, nucleation._inputs.dTsigma);

    // Fill in nucleation data structures, and assign nucleation undercooling values to potential nucleation events
    // Potential nucleation grains are only associated with liquid cells in layer 1 - they will be initialized for each
    // successive layer when layer 1 in complete
    nucleation.placeNuclei(temperature, inputs.RNGSeed, 1, nx, ny, nz_layer, ny_local, y_offset, z_layer_bottom, id,
                           AtNorthBoundary, AtSouthBoundary);

    // Copy results back to host to check
    auto NucleiLocation_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), nucleation.NucleiLocations);
    auto NucleiGrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), nucleation.NucleiGrainID);

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
        // Cell location is a local position (relative to the bottom of the layer)
        int index = NucleiLocation_Host(n);
        int coord_z = getCoordZ(index, nx, ny_local);
        int coord_y = getCoordY(index, nx, ny_local);
        int coord_z_AllLayers = coord_z + z_layer_bottom;
        // Expected nucleation time with remelting can be one of 3 possibilities, depending on the associated
        // solidification event
        int Expected_NucleationTimeNoRM = coord_z_AllLayers + coord_y + y_offset + 2;
        int AssociatedSEvent = nucleation.NucleationTimes_Host(n) / DomainSize;
        int Expected_NucleationTimeRM = Expected_NucleationTimeNoRM + AssociatedSEvent * DomainSize;
        EXPECT_EQ(nucleation.NucleationTimes_Host(n), Expected_NucleationTimeRM);

        // Are the nucleation events in order of the time steps at which they may occur?
        if (n < nucleation.PossibleNuclei - 2) {
            EXPECT_LE(nucleation.NucleationTimes_Host(n), nucleation.NucleationTimes_Host(n + 1));
        }
    }
}

// Tests nucleate_grain member function
void testNucleateGrain() {

    using memory_space = TEST_MEMSPACE;
    using view_int = Kokkos::View<int *, TEST_MEMSPACE>;
    using view_int_host = typename view_int::HostMirror;

    // Create views - 125 cells, 75 of which are part of the active layer of the domain (Z = 2-5)
    int nx = 5;
    int ny_local = 5;
    int nz = 5;
    int z_layer_bottom = 2;
    int nz_layer = 3;
    int DomainSize = nx * ny_local * nz_layer;
    int DomainSize_AllLayers = nx * ny_local * nz;

    // Empty inputs struct
    Inputs inputs;

    // All cells have GrainID of 0, CellType of Liquid - with the exception of the locations where the nucleation events
    // are unable to occur
    CellData<memory_space> cellData(DomainSize_AllLayers, DomainSize, nx, ny_local, z_layer_bottom, inputs.substrate);
    Kokkos::deep_copy(cellData.CellType_AllLayers, Liquid);
    auto CellType = cellData.getCellTypeSubview();
    auto CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
    auto GrainID = cellData.getGrainIDSubview();
    auto GrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);

    // Create test nucleation data - 10 possible events
    int PossibleNuclei = 10;
    Nucleation<memory_space> nucleation(PossibleNuclei, 1.0, inputs.nucleation);
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
    view_int SteeringVector(Kokkos::ViewAllocateWithoutInitializing("SteeringVector"), DomainSize);
    // Initialize steering vector size to 0
    view_int numSteer("SteeringVector", 1);

    // Take enough time steps such that every nucleation event has a chance to occur
    for (int cycle = 0; cycle < 10; cycle++) {
        nucleation.nucleate_grain(cycle, cellData, z_layer_bottom, nx, ny_local, SteeringVector, numSteer);
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
        int index = SuccessfulNuc_CellLocations[nevent];
        EXPECT_EQ(CellType_Host(index), FutureActive);
        EXPECT_EQ(GrainID_Host(index), SuccessfulNuc_GrainIDs[nevent]);
        bool OnSteeringVector = false;
        for (int svloc = 0; svloc < 7; svloc++) {
            if (SteeringVector_Host(svloc) == index) {
                OnSteeringVector = true;
                break;
            }
        }
        EXPECT_TRUE(OnSteeringVector);
    }
}

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, nucleation) {
    testNucleiInit();
    testNucleateGrain();
}
} // end namespace Test
