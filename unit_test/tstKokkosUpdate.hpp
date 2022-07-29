#include <Kokkos_Core.hpp>

#include "CAinitialize.hpp"
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
    // Initialize steering vector size to 0
    ViewI numSteer("SteeringVector", 1);

    // Take enough time steps such that every nucleation event has a chance to occur
    for (int cycle = 0; cycle <= (2 * id + 11); cycle++) {
        Nucleation(cycle, SuccessfulNucEvents_ThisRank, NucleationCounter, PossibleNuclei_ThisRankThisLayer,
                   NucleationTimes_Host, NucleiLocation, NucleiGrainID, CellType, GrainID, ZBound_Low, MyXSlices,
                   MyYSlices, SteeringVector, numSteer);
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

void testFillSteeringVector_Remelt() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // Create views - each rank has 125 cells, 75 of which are part of the active region of the domain
    int MyXSlices = 5;
    int MyYSlices = 5;
    int MyYOffset = 5 * id;
    int nz = 5;
    int ZBound_Low = 2;
    int nzActive = 3;
    int LocalActiveDomainSize = MyXSlices * MyYSlices * nzActive;
    int LocalDomainSize = MyXSlices * MyYSlices * nz;

    // MPI rank locations relative to the global grid
    int DecompositionStrategy = 1; // 1D decomposition in the Y direction
    bool AtNorthBoundary = false;
    bool AtSouthBoundary = false;
    bool AtEastBoundary = false;
    bool AtWestBoundary = false;

    // Buffers for ghost node data (fixed size)
    int BufSizeX = MyXSlices;
    int BufSizeY = 0;
    int BufSizeZ = nzActive;

    // Send/recv buffers for ghost node data - initialize with values of 1.0
    Buffer2D BufferSouthSend(Kokkos::ViewAllocateWithoutInitializing("BufferSouthSend"), BufSizeX * BufSizeZ, 5);
    Buffer2D BufferNorthSend(Kokkos::ViewAllocateWithoutInitializing("BufferNorthSend"), BufSizeX * BufSizeZ, 5);
    Buffer2D BufferEastSend(Kokkos::ViewAllocateWithoutInitializing("BufferEastSend"), BufSizeY * BufSizeZ, 5);
    Buffer2D BufferWestSend(Kokkos::ViewAllocateWithoutInitializing("BufferWestSend"), BufSizeY * BufSizeZ, 5);
    Buffer2D BufferNorthEastSend(Kokkos::ViewAllocateWithoutInitializing("BufferNorthEastSend"), BufSizeZ, 5);
    Buffer2D BufferNorthWestSend(Kokkos::ViewAllocateWithoutInitializing("BufferNorthWestSend"), BufSizeZ, 5);
    Buffer2D BufferSouthEastSend(Kokkos::ViewAllocateWithoutInitializing("BufferSouthEastSend"), BufSizeZ, 5);
    Buffer2D BufferSouthWestSend(Kokkos::ViewAllocateWithoutInitializing("BufferSouthWestSend"), BufSizeZ, 5);
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
        int GlobalZ = i / (MyXSlices * MyYSlices);
        int Rem = i % (MyXSlices * MyYSlices);
        int RankY = Rem % MyYSlices;
        // Let cells be assigned GrainIDs based on the rank ID
        // Cells on rank 0 have grain ID 0, rank 1 have grain ID 1, etc
        GrainID_Host(i) = id;
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
        UndercoolingChange_Host(i) = 0.2 * id;
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
        FillSteeringVector_Remelt(cycle, LocalActiveDomainSize, MyXSlices, MyYSlices, NeighborX, NeighborY, NeighborZ,
                                  CritTimeStep, UndercoolingCurrent, UndercoolingChange, CellType, GrainID, ZBound_Low,
                                  nzActive, SteeringVector, numSteer, numSteer_Host, MeltTimeStep, BufSizeX, BufSizeY,
                                  AtNorthBoundary, AtSouthBoundary, AtEastBoundary, AtWestBoundary, BufferWestSend,
                                  BufferEastSend, BufferNorthSend, BufferSouthSend, BufferNorthEastSend,
                                  BufferNorthWestSend, BufferSouthEastSend, BufferSouthWestSend, DecompositionStrategy);
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
        int GlobalZ = i / (MyXSlices * MyYSlices);
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
                EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(i), (numcycles - CritTimeStep_Host(i)) * 0.2 * id);
                if ((id == 0) || (GlobalZ == 4))
                    EXPECT_EQ(CellType_Host(i), Liquid);
                else {
                    EXPECT_EQ(CellType_Host(i), FutureActive);
                    FutureActiveCells++;
                }
            }
        }
    }
    std::cout << "Future active cells rank " << id << " : " << FutureActiveCells << std::endl;
    // Check the steering vector values on the host
    EXPECT_EQ(FutureActiveCells, numSteer_Host(0));
    for (int i = 0; i < FutureActiveCells; i++) {
        // This cell should correspond to a cell at GlobalZ = 3 (RankZ = 1), and some X and Y
        int LowerBoundCellLocation = MyXSlices * MyYSlices - 1;
        int UpperBoundCellLocation = 2 * MyXSlices * MyYSlices;
        EXPECT_GT(SteeringVector_Host(i), LowerBoundCellLocation);
        EXPECT_LT(SteeringVector_Host(i), UpperBoundCellLocation);
    }

    // Check that the buffer values were either appropriately set to zeros (if the cell underwent melting) or remained
    // at 1.0
    for (int i = 0; i < BufSizeX * BufSizeZ; i++) {
        int RankZ = i / BufSizeX;
        int GlobalZ = RankZ + ZBound_Low;
        int RankX = i % BufSizeX;
        int NorthCellCoordinate = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + (MyYSlices - 1);
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
        int SouthCellCoordinate = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + (MyYSlices - 1);
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

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, cell_update_tests) {
    testNucleation();
    testFillSteeringVector_Remelt();
}

} // end namespace Test
