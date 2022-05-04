#include <Kokkos_Core.hpp>

#include "CAtypes.hpp"
#include "CAupdate.hpp"

#include <gtest/gtest.h>

#include "mpi.h"

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
void testNucleation() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    int SuccessfulNucEvents_ThisRank = 0; // nucleation event counter
    // Counters for nucleation events (successful or not) - host and device
    int NucleationCounter = 0;

    // Create views - each rank has 125 cells, 75 of which are part of the active region of the domain
    int MyXSlices = 5;
    int MyYSlices = 5;
    int nz = 5;
    int ZBound_Low = 2;
    int nzActive = 3;
    int LocalActiveDomainSize = MyXSlices * MyYSlices * nzActive;
    int LocalDomainSize = MyXSlices * MyYSlices * nz;

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
        NucleationTimes_Host(n) = 2 * id + n;
        NucleiLocation_Host(n) = ZBound_Low * MyXSlices * MyYSlices + n;
        // Give these nucleation events grain IDs based on their order, starting with -1 and counting down
        NucleiGrainID_Host(n) = -(n + 1);
    }
    // Include the case where 2 potential nucleation events (3 and 4) happen on the same time step - both successful
    NucleationTimes_Host(3) = NucleationTimes_Host(4);

    // Include the case where a potential nucleation event (2) is unsuccessful
    int UnsuccessfulLocA = ZBound_Low * MyXSlices * MyYSlices + 2;
    CellType_Host(UnsuccessfulLocA) = Active;
    GrainID_Host(UnsuccessfulLocA) = 1;

    // Include the case where 2 potential nucleation events (5 and 6) happen on the same time step - the first
    // successful, the second unsuccessful
    NucleationTimes_Host(5) = NucleationTimes_Host(6);
    int UnsuccessfulLocC = ZBound_Low * MyXSlices * MyYSlices + 6;
    CellType_Host(UnsuccessfulLocC) = Active;
    GrainID_Host(UnsuccessfulLocC) = 2;

    // Include the case where 2 potential nucleation events (8 and 9) happen on the same time step - the first
    // unsuccessful, the second successful
    NucleationTimes_Host(8) = NucleationTimes_Host(9);
    int UnsuccessfulLocB = ZBound_Low * MyXSlices * MyYSlices + 8;
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
    ViewI_H numSteer_Host("SteeringVectorSize_Host", 1);
    ViewI numSteer = Kokkos::create_mirror_view_and_copy(memory_space(), numSteer_Host);

    // Take enough time steps such that every nucleation event has a chance to occur
    for (int cycle = 0; cycle <= (2 * id + 11); cycle++) {
        Nucleation(cycle, SuccessfulNucEvents_ThisRank, NucleationCounter, PossibleNuclei_ThisRankThisLayer,
                   NucleationTimes_Host, NucleiLocation, NucleiGrainID, CellType, GrainID, ZBound_Low, MyXSlices,
                   MyYSlices, SteeringVector, numSteer);
    }

    // Copy CellType, SteeringVector, numSteer, GrainID back to host to check nucleation results
    CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
    ViewI_H SteeringVector_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), SteeringVector);
    numSteer_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), numSteer);
    GrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);

    // Check that all 10 possible nucleation events were attempted
    EXPECT_EQ(NucleationCounter, 10);
    // Check that 7 of the 10 nucleation events were successful
    EXPECT_EQ(SuccessfulNucEvents_ThisRank, 7);
    EXPECT_EQ(numSteer_Host(0), 7);

    // Ensure that the 3 events that should not have occurred, did not occur
    // These cells should be untouched - active type, and same GrainID they was initialized with
    EXPECT_EQ(CellType_Host(ZBound_Low * MyXSlices * MyYSlices + 2), Active);
    EXPECT_EQ(GrainID_Host(ZBound_Low * MyXSlices * MyYSlices + 2), 1);
    EXPECT_EQ(CellType_Host(ZBound_Low * MyXSlices * MyYSlices + 6), Active);
    EXPECT_EQ(GrainID_Host(ZBound_Low * MyXSlices * MyYSlices + 6), 2);
    EXPECT_EQ(CellType_Host(ZBound_Low * MyXSlices * MyYSlices + 8), Active);
    EXPECT_EQ(GrainID_Host(ZBound_Low * MyXSlices * MyYSlices + 8), 3);

    // Check that the successful nucleation events occurred as expected
    std::vector<int> SuccessfulNuc_GrainIDs{-1, -2, -4, -5, -6, -8, -10};
    for (int n = 0; n < 7; n++) {
        // Cell's location within the active layer
        int CellLocation_ThisLayer = SteeringVector_Host(n);
        int RankZ = CellLocation_ThisLayer / (MyXSlices * MyYSlices);
        int Rem = CellLocation_ThisLayer % (MyXSlices * MyYSlices);
        int RankX = Rem / MyYSlices;
        int RankY = Rem % MyYSlices;
        int GlobalZ = RankZ + ZBound_Low;
        // Cell location within the overall domain
        int CellLocation_AllLayers = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
        // Check that this cell type is marked as "FutureActive", and the GrainID matches the expected value
        EXPECT_EQ(CellType_Host(CellLocation_AllLayers), FutureActive);
        EXPECT_EQ(GrainID_Host(CellLocation_AllLayers), SuccessfulNuc_GrainIDs[n]);
    }
}

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, cell_update_tests) { testNucleation(); }

} // end namespace Test
