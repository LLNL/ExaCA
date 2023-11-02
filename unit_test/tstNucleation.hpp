// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAcelldata.hpp"
#include "CAfunctions.hpp"
#include "CAgrid.hpp"
#include "CAinputs.hpp"
#include "CAinterface.hpp"
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
    // default inputs struct
    Inputs inputs;
    // manually set grid for test of 2 layer problem
    int NumberOfLayers_temp = 2;
    Grid grid(NumberOfLayers_temp);
    // Create test data
    // Only Z = 1 through 4 is part of this layer
    grid.nz_layer = 4;
    grid.z_layer_bottom = 1;
    grid.nx = 4;
    // Each rank is assigned a different portion of the domain in Y
    grid.ny = 2 * np;
    grid.ny_local = 2;
    grid.y_offset = 2 * id;
    grid.deltax = 1;
    grid.DomainSize = grid.nx * grid.ny_local * grid.nz_layer;
    grid.BottomOfCurrentLayer = grid.getBottomOfCurrentLayer();
    grid.TopOfCurrentLayer = grid.getTopOfCurrentLayer();
    grid.LayerRange = std::make_pair(grid.BottomOfCurrentLayer, grid.TopOfCurrentLayer);
    // MPI rank locations relative to the global grid
    if (id == 0)
        grid.AtSouthBoundary = true;
    else
        grid.AtSouthBoundary = false;
    if (id == np - 1)
        grid.AtNorthBoundary = true;
    else
        grid.AtNorthBoundary = false;

    // There are 40 * np total cells in this domain (nx * ny * nz)
    // Each rank has 40 cells - the top 32 cells are part of the active layer and are candidates for nucleation
    // assignment
    int MaxPotentialNuclei_PerPass = 4 * np;
    // A cell can solidify 1-3 times
    int MaxSolidificationEvents_Count = 3;
    // Manually set nucleation parameters
    // This nucleation density ensures there will be 4 potential nuclei per MPI rank present
    // without remelting (each cell solidifies once)
    inputs.nucleation.dTN = 1;
    inputs.nucleation.dTsigma = 0.0001;
    inputs.nucleation.NMax = 0.125;

    // Allocate temperature data structures
    Temperature<memory_space> temperature(grid.DomainSize, grid.NumberOfLayers, inputs.temperature, 1);
    // Resize LayerTimeTempHistory with the known max number of solidification events
    Kokkos::resize(temperature.LayerTimeTempHistory, grid.DomainSize, MaxSolidificationEvents_Count, 3);
    // Initialize MaxSolidificationEvents to 3 for each layer. LayerTimeTempHistory and NumberOfSolidificationEvents are
    // initialized for each cell on the host and copied to the device
    ViewI_H MaxSolidificationEvents_Host(Kokkos::ViewAllocateWithoutInitializing("MaxSolidificationEvents_Host"),
                                         grid.NumberOfLayers);
    MaxSolidificationEvents_Host(0) = MaxSolidificationEvents_Count;
    MaxSolidificationEvents_Host(1) = MaxSolidificationEvents_Count;
    // Cells solidify 1, 2, or 3 times, depending on their X coordinate
    Kokkos::parallel_for(
        "NumSolidificationEventsInit", grid.nz_layer, KOKKOS_LAMBDA(const int &coord_z) {
            for (int coord_x = 0; coord_x < grid.nx; coord_x++) {
                for (int coord_y = 0; coord_y < grid.ny_local; coord_y++) {
                    int index = grid.get1Dindex(coord_x, coord_y, coord_z);
                    if (coord_x < grid.nx / 2 - 1)
                        temperature.NumberOfSolidificationEvents(index) = 3;
                    else if (coord_x < grid.nx / 2)
                        temperature.NumberOfSolidificationEvents(index) = 2;
                    else
                        temperature.NumberOfSolidificationEvents(index) = 1;
                }
            }
        });
    Kokkos::fence();
    Kokkos::parallel_for(
        "LayerTimeTempHistoryInit", MaxSolidificationEvents_Count, KOKKOS_LAMBDA(const int &n) {
            for (int coord_z = 0; coord_z < grid.nz_layer; coord_z++) {
                for (int coord_x = 0; coord_x < grid.nx; coord_x++) {
                    for (int coord_y = 0; coord_y < grid.ny_local; coord_y++) {
                        int index = grid.get1Dindex(coord_x, coord_y, coord_z);
                        int coord_z_AllLayers = coord_z + grid.z_layer_bottom;
                        if (n < temperature.NumberOfSolidificationEvents(index)) {
                            // melting time step depends on solidification event number
                            temperature.LayerTimeTempHistory(index, n, 0) =
                                coord_z_AllLayers + coord_y + grid.y_offset + (grid.DomainSize * n);
                            // liquidus time stemp depends on solidification event number
                            temperature.LayerTimeTempHistory(index, n, 1) =
                                coord_z_AllLayers + coord_y + grid.y_offset + 1 + (grid.DomainSize * n);
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
    int EstimatedNuclei_ThisRankThisLayer = inputs.nucleation.NMax * pow(grid.deltax, 3) * grid.DomainSize;
    Nucleation<memory_space> nucleation(
        EstimatedNuclei_ThisRankThisLayer, grid.deltax, inputs.nucleation,
        100); // NucleiGrainID should start at -101 - supply optional input arg to constructor
    // Ensure nucleation inputs in nucleation struct were correctly initialized
    EXPECT_DOUBLE_EQ(inputs.nucleation.NMax, nucleation._inputs.NMax);
    EXPECT_DOUBLE_EQ(inputs.nucleation.dTN, nucleation._inputs.dTN);
    EXPECT_DOUBLE_EQ(inputs.nucleation.dTsigma, nucleation._inputs.dTsigma);

    // Fill in nucleation data structures, and assign nucleation undercooling values to potential nucleation events
    // Potential nucleation grains are only associated with liquid cells in layer 1 - they will be initialized for each
    // successive layer when layer 1 in complete
    nucleation.placeNuclei(temperature, inputs.RNGSeed, 1, grid, id);

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
        int coord_z = grid.getCoordZ(index);
        int coord_y = grid.getCoordY(index);
        int coord_z_AllLayers = coord_z + grid.z_layer_bottom;
        // Expected nucleation time with remelting can be one of 3 possibilities, depending on the associated
        // solidification event
        int Expected_NucleationTimeNoRM = coord_z_AllLayers + coord_y + grid.y_offset + 2;
        int AssociatedSEvent = nucleation.NucleationTimes_Host(n) / grid.DomainSize;
        int Expected_NucleationTimeRM = Expected_NucleationTimeNoRM + AssociatedSEvent * grid.DomainSize;
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

    // default inputs struct
    Inputs inputs;
    // manually set grid
    Grid grid;
    // Create views - 125 cells, 75 of which are part of the active layer of the domain (Z = 2-5)
    grid.nx = 5;
    grid.ny_local = 5;
    grid.nz = 5;
    grid.z_layer_bottom = 2;
    grid.nz_layer = 3;
    grid.deltax = 1.0;
    grid.DomainSize = grid.nx * grid.ny_local * grid.nz_layer;
    grid.DomainSize_AllLayers = grid.nx * grid.ny_local * grid.nz;
    grid.BottomOfCurrentLayer = grid.getBottomOfCurrentLayer();
    grid.TopOfCurrentLayer = grid.getTopOfCurrentLayer();
    grid.LayerRange = std::make_pair(grid.BottomOfCurrentLayer, grid.TopOfCurrentLayer);

    // All cells have GrainID of 0, CellType of Liquid - with the exception of the locations where the nucleation events
    // are unable to occur
    CellData<memory_space> cellData(grid.DomainSize_AllLayers, inputs.substrate);
    Kokkos::deep_copy(cellData.CellType_AllLayers, Liquid);
    auto CellType = cellData.getCellTypeSubview(grid);
    auto CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
    auto GrainID = cellData.getGrainIDSubview(grid);
    auto GrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);

    // Create test nucleation data - 10 possible events
    int PossibleNuclei = 10;
    Nucleation<memory_space> nucleation(PossibleNuclei, grid.deltax, inputs.nucleation);
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

    // Interface struct
    Interface<memory_space> interface(grid.DomainSize);
    // Take enough time steps such that every nucleation event has a chance to occur
    for (int cycle = 0; cycle < 10; cycle++) {
        nucleation.nucleate_grain(cycle, grid, cellData, interface);
    }

    // Copy CellType, SteeringVector, numSteer, GrainID back to host to check nucleation results
    CellType = cellData.getCellTypeSubview(grid);
    CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
    auto SteeringVector_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.SteeringVector);
    auto numSteer_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), interface.numSteer);
    GrainID = cellData.getGrainIDSubview(grid);
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
