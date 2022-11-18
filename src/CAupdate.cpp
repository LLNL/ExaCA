// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAupdate.hpp"
#include "CAprint.hpp"
#include "mpi.h"

#include <cmath>

//*****************************************************************************/
void Nucleation(int cycle, int &SuccessfulNucEvents_ThisRank, int &NucleationCounter, int PossibleNuclei_ThisRank,
                ViewI_H NucleationTimes_H, ViewI NucleiLocations, ViewI NucleiGrainID, ViewI CellType, ViewI GrainID,
                int ZBound_Low, int nx, int MyYSlices, ViewI SteeringVector, ViewI numSteer_G) {

    // Is there nucleation left in this layer to check?
    if (NucleationCounter < PossibleNuclei_ThisRank) {
        // Is there at least one potential nucleation event on this rank, at this time step?
        if (cycle == NucleationTimes_H(NucleationCounter)) {
            bool NucleationCheck = true;
            int FirstEvent = NucleationCounter; // first potential nucleation event to check
            // Are there any other nucleation events this time step to check?
            while (NucleationCheck) {
                NucleationCounter++;
                // If the previous nucleation event was the last one for this layer of the simulation, exit loop
                if (NucleationCounter == PossibleNuclei_ThisRank)
                    break;
                // If the next nucleation event corresponds to a future time step, finish check
                if (cycle != NucleationTimes_H(NucleationCounter))
                    NucleationCheck = false;
            }
            int LastEvent = NucleationCounter;
            // parallel_reduce checks each potential nucleation event this time step (FirstEvent, up to but not
            // including LastEvent)
            int NucleationThisDT = 0; // return number of successful event from parallel_reduce
            // Launch kokkos kernel - check if the corresponding CA cell location is liquid
            Kokkos::parallel_reduce(
                "NucleiUpdateLoop", Kokkos::RangePolicy<>(FirstEvent, LastEvent),
                KOKKOS_LAMBDA(const int NucleationCounter_Device, int &update) {
                    int NucleationEventLocation_GlobalGrid = NucleiLocations(NucleationCounter_Device);
                    int update_val =
                        FutureActive; // added to steering vector to become a new active cell as part of cellcapture
                    int old_val = Liquid;
                    int OldCellTypeValue = Kokkos::atomic_compare_exchange(
                        &CellType(NucleationEventLocation_GlobalGrid), old_val, update_val);
                    if (OldCellTypeValue == Liquid) {
                        // Successful nucleation event - atomic update of cell type, proceeded if the atomic
                        // exchange is successful (cell was liquid) Add future active cell location to steering
                        // vector and change cell type, assign new Grain ID
                        GrainID(NucleationEventLocation_GlobalGrid) = NucleiGrainID(NucleationCounter_Device);
                        int GlobalZ = NucleationEventLocation_GlobalGrid / (nx * MyYSlices);
                        int Rem = NucleationEventLocation_GlobalGrid % (nx * MyYSlices);
                        int RankX = Rem / MyYSlices;
                        int RankY = Rem % MyYSlices;
                        int RankZ = GlobalZ - ZBound_Low;
                        int NucleationEventLocation_LocalGrid = RankZ * nx * MyYSlices + RankX * MyYSlices + RankY;
                        SteeringVector(Kokkos::atomic_fetch_add(&numSteer_G(0), 1)) = NucleationEventLocation_LocalGrid;
                        // This undercooled liquid cell is now a nuclei (no nuclei are in the ghost nodes - halo
                        // exchange routine GhostNodes1D or GhostNodes2D is used to fill these)
                        update++;
                    }
                },
                NucleationThisDT);
            // Update the number of successful nuclei counter with the number of successful nucleation events from this
            // time step (NucleationThisDT)
            SuccessfulNucEvents_ThisRank += NucleationThisDT;
        }
    }
}

//*****************************************************************************/
// Determine which cells are associated with the "steering vector" of cells that are either active, or becoming active
// this time step
void FillSteeringVector_NoRemelt(int cycle, int LocalActiveDomainSize, int nx, int MyYSlices, ViewI CritTimeStep,
                                 ViewF UndercoolingCurrent, ViewF UndercoolingChange, ViewI CellType, int ZBound_Low,
                                 int layernumber, ViewI LayerID, ViewI SteeringVector, ViewI numSteer,
                                 ViewI_H numSteer_Host) {

    // Cells associated with this layer that are not solid type but have passed the liquidus (crit time step) have their
    // undercooling values updated Cells that meet the aforementioned criteria and are active type should be added to
    // the steering vector
    Kokkos::parallel_for(
        "FillSV", LocalActiveDomainSize, KOKKOS_LAMBDA(const int &D3D1ConvPosition) {
            // Cells of interest for the CA
            int RankZ = D3D1ConvPosition / (nx * MyYSlices);
            int Rem = D3D1ConvPosition % (nx * MyYSlices);
            int RankX = Rem / MyYSlices;
            int RankY = Rem % MyYSlices;
            int GlobalZ = RankZ + ZBound_Low;
            int GlobalD3D1ConvPosition = GlobalZ * nx * MyYSlices + RankX * MyYSlices + RankY;
            int cellType = CellType(GlobalD3D1ConvPosition);

            int layerCheck = (LayerID(GlobalD3D1ConvPosition) <= layernumber);
            int isNotSolid = (cellType != Solid);
            int pastCritTime = (cycle > CritTimeStep(GlobalD3D1ConvPosition));

            int cell_Liquid = (cellType == Liquid);
            int cell_Active = (cellType == Active);

            if (layerCheck && isNotSolid && pastCritTime) {
                UndercoolingCurrent(GlobalD3D1ConvPosition) +=
                    UndercoolingChange(GlobalD3D1ConvPosition) * (cell_Liquid + cell_Active);
                if (cell_Active) {
                    SteeringVector(Kokkos::atomic_fetch_add(&numSteer(0), 1)) = D3D1ConvPosition;
                }
            }
        });
    Kokkos::deep_copy(numSteer_Host, numSteer);
}

//*****************************************************************************/
// Determine which cells are associated with the "steering vector" of cells that are either active, or becoming active
// this time step - version with remelting
void FillSteeringVector_Remelt(int cycle, int LocalActiveDomainSize, int nx, int MyYSlices, NList NeighborX,
                               NList NeighborY, NList NeighborZ, ViewI CritTimeStep, ViewF UndercoolingCurrent,
                               ViewF UndercoolingChange, ViewI CellType, ViewI GrainID, int ZBound_Low, int nzActive,
                               ViewI SteeringVector, ViewI numSteer, ViewI_H numSteer_Host, ViewI MeltTimeStep,
                               int BufSizeX, bool AtNorthBoundary, bool AtSouthBoundary, Buffer2D BufferNorthSend,
                               Buffer2D BufferSouthSend) {

    Kokkos::parallel_for(
        "FillSV_RM", LocalActiveDomainSize, KOKKOS_LAMBDA(const int &D3D1ConvPosition) {
            // Coordinate of this cell on the "global" (all cells in the Z direction) grid
            int RankZ = D3D1ConvPosition / (nx * MyYSlices);
            int Rem = D3D1ConvPosition % (nx * MyYSlices);
            int RankX = Rem / MyYSlices;
            int RankY = Rem % MyYSlices;
            int GlobalZ = RankZ + ZBound_Low;
            int GlobalD3D1ConvPosition = GlobalZ * nx * MyYSlices + RankX * MyYSlices + RankY;

            int cellType = CellType(GlobalD3D1ConvPosition);
            bool isNotSolid = ((cellType != TempSolid) && (cellType != Solid));
            bool atMeltTime = (cycle == MeltTimeStep(GlobalD3D1ConvPosition));
            bool pastCritTime = (cycle > CritTimeStep(GlobalD3D1ConvPosition));
            if ((atMeltTime) && ((cellType == TempSolid) || (cellType == Active))) {
                // This cell should be a liquid cell
                CellType(GlobalD3D1ConvPosition) = Liquid;
                // Reset current undercooling to zero
                UndercoolingCurrent(GlobalD3D1ConvPosition) = 0.0;
                // Remove solid cell data from the buffer
                loadghostnodes(0, 0, 0, 0, 0, BufSizeX, MyYSlices, RankX, RankY, RankZ, AtNorthBoundary,
                               AtSouthBoundary, BufferSouthSend, BufferNorthSend);
            }
            else if ((isNotSolid) && (pastCritTime)) {
                // Update cell undercooling
                UndercoolingCurrent(GlobalD3D1ConvPosition) += UndercoolingChange(GlobalD3D1ConvPosition);
                if (cellType == Active) {
                    // Add active cells below liquidus to steering vector
                    SteeringVector(Kokkos::atomic_fetch_add(&numSteer(0), 1)) = D3D1ConvPosition;
                }
                else if ((cellType == Liquid) && (GrainID(GlobalD3D1ConvPosition) != 0)) {
                    // If this cell borders at least one solid/tempsolid cell and is part of a grain, it should become
                    // active
                    for (int l = 0; l < 26; l++) {
                        // "l" correpsponds to the specific neighboring cell
                        // Local coordinates of adjacent cell center
                        int MyNeighborX = RankX + NeighborX[l];
                        int MyNeighborY = RankY + NeighborY[l];
                        int MyNeighborZ = RankZ + NeighborZ[l];
                        if ((MyNeighborX >= 0) && (MyNeighborX < nx) && (MyNeighborY >= 0) &&
                            (MyNeighborY < MyYSlices) && (MyNeighborZ < nzActive) && (MyNeighborZ >= 0)) {
                            int GlobalNeighborD3D1ConvPosition =
                                (MyNeighborZ + ZBound_Low) * nx * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
                            if ((CellType(GlobalNeighborD3D1ConvPosition) == TempSolid) ||
                                (CellType(GlobalNeighborD3D1ConvPosition) == Solid) || (RankZ == 0)) {
                                // Cell activation to be performed as part of steering vector
                                l = 26;
                                SteeringVector(Kokkos::atomic_fetch_add(&numSteer(0), 1)) = D3D1ConvPosition;
                                CellType(GlobalD3D1ConvPosition) =
                                    FutureActive; // this cell cannot be captured - is being activated
                            }
                        }
                    }
                }
            }
        });
    Kokkos::fence();

    // Copy size of steering vector (containing positions of undercooled liquid/active cells) to the host
    Kokkos::deep_copy(numSteer_Host, numSteer);
}

//*****************************************************************************/
// Jump to the next time step with work to be done, if nothing left to do in the near future
// Without remelting, the cells of interest are undercooled liquid cells, and the view checked for future work is
// CritTimeStep With remelting, the cells of interest are active cells, and the view checked for future work is
// MeltTimeStep Print intermediate output during this jump if PrintIdleMovieFrames = true
void JumpTimeStep(int &cycle, unsigned long int RemainingCellsOfInterest, unsigned long int LocalIncompleteCells,
                  ViewI FutureWorkView, int LocalActiveDomainSize, int MyYSlices, int ZBound_Low, bool RemeltingYN,
                  ViewI CellType, ViewI LayerID, int id, int layernumber, int np, int nx, int ny, int nz, int MyYOffset,
                  ViewI GrainID, ViewI CritTimeStep, ViewF GrainUnitVector, ViewF UndercoolingChange,
                  ViewF UndercoolingCurrent, std::string OutputFile, int NGrainOrientations, std::string PathToOutput,
                  int &IntermediateFileCounter, int nzActive, double deltax, float XMin, float YMin, float ZMin,
                  int NumberOfLayers, int &XSwitch, std::string TemperatureDataType, bool PrintIdleMovieFrames,
                  int MovieFrameInc, bool PrintBinary, int FinishTimeStep = 0) {

    MPI_Bcast(&RemainingCellsOfInterest, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    if (RemainingCellsOfInterest == 0) {
        // If this rank still has cells that will later undergo transformation (LocalIncompleteCells > 0), check when
        // the next superheated cells go below the liquidus (no remelting) or check when the next solid cells go above
        // the liquidus (remelting) Otherwise, assign the largest possible time step as the next time work needs to be
        // done on the rank
        unsigned long int NextWorkTimeStep;
        if (LocalIncompleteCells > 0) {
            Kokkos::parallel_reduce(
                "CheckNextTSForWork", LocalActiveDomainSize,
                KOKKOS_LAMBDA(const int &D3D1ConvPosition, unsigned long int &tempv) {
                    int RankZ = D3D1ConvPosition / (nx * MyYSlices);
                    int Rem = D3D1ConvPosition % (nx * MyYSlices);
                    int RankX = Rem / MyYSlices;
                    int RankY = Rem % MyYSlices;
                    int GlobalZ = RankZ + ZBound_Low;
                    int GlobalD3D1ConvPosition = GlobalZ * nx * MyYSlices + RankX * MyYSlices + RankY;
                    unsigned long int NextWorkTimeStep_ThisCell =
                        (unsigned long int)(FutureWorkView(GlobalD3D1ConvPosition));
                    // remelting/no remelting criteria for a cell to be associated with future work
                    if (((!(RemeltingYN)) && (CellType(GlobalD3D1ConvPosition) == Liquid) &&
                         (LayerID(GlobalD3D1ConvPosition) == layernumber)) ||
                        ((RemeltingYN) && (CellType(GlobalD3D1ConvPosition) == TempSolid))) {
                        if (NextWorkTimeStep_ThisCell < tempv)
                            tempv = NextWorkTimeStep_ThisCell;
                    }
                },
                Kokkos::Min<unsigned long int>(NextWorkTimeStep));
        }
        else
            NextWorkTimeStep = INT_MAX;

        unsigned long int GlobalNextWorkTimeStep;
        MPI_Allreduce(&NextWorkTimeStep, &GlobalNextWorkTimeStep, 1, MPI_UNSIGNED_LONG, MPI_MIN, MPI_COMM_WORLD);
        if ((GlobalNextWorkTimeStep - cycle) > 5000) {
            if (PrintIdleMovieFrames) {
                // Print any movie frames that occur during the skipped time steps
                for (unsigned long int cycle_jump = cycle + 1; cycle_jump < GlobalNextWorkTimeStep; cycle_jump++) {
                    if (cycle_jump % MovieFrameInc == 0) {
                        // Print current state of ExaCA simulation (up to and including the current layer's data)
                        // Host mirrors of CellType and GrainID are not maintained - pass device views and perform
                        // copy inside of subroutine
                        PrintExaCAData(id, layernumber, np, nx, ny, nz, MyYSlices, MyYOffset, GrainID, CritTimeStep,
                                       GrainUnitVector, LayerID, CellType, UndercoolingChange, UndercoolingCurrent,
                                       OutputFile, NGrainOrientations, PathToOutput, 0, false, false, false, true,
                                       false, IntermediateFileCounter, ZBound_Low, nzActive, deltax, XMin, YMin, ZMin,
                                       NumberOfLayers, PrintBinary);
                        IntermediateFileCounter++;
                    }
                }
            }
            // Jump to next time step when solidification starts again
            cycle = GlobalNextWorkTimeStep - 1;
            if (id == 0)
                std::cout << "Jumping to cycle " << cycle + 1 << std::endl;
        }
        // If all cells have cooled below solidus, jump to next layer (no remelting/using input temp data only)
        if ((!(RemeltingYN)) && (TemperatureDataType == "R") && (cycle >= FinishTimeStep))
            XSwitch = 1;
    }
}

//*****************************************************************************/
// Prints intermediate code output to stdout, checks to see if solidification is complete
void IntermediateOutputAndCheck(int id, int np, int &cycle, int MyYSlices, int MyYOffset, int LocalDomainSize,
                                int LocalActiveDomainSize, int nx, int ny, int nz, int nzActive, double deltax,
                                float XMin, float YMin, float ZMin, int SuccessfulNucEvents_ThisRank, int &XSwitch,
                                ViewI CellType, ViewI CritTimeStep, ViewI GrainID, std::string TemperatureDataType,
                                int *FinishTimeStep, int layernumber, int, int ZBound_Low, int NGrainOrientations,
                                ViewI LayerID, ViewF GrainUnitVector, ViewF UndercoolingChange,
                                ViewF UndercoolingCurrent, std::string PathToOutput, std::string OutputFile,
                                bool PrintIdleMovieFrames, int MovieFrameInc, int &IntermediateFileCounter,
                                int NumberOfLayers, bool PrintBinary) {

    unsigned long int LocalSuperheatedCells;
    unsigned long int LocalUndercooledCells;
    unsigned long int LocalActiveCells;
    unsigned long int LocalSolidCells;
    unsigned long int LocalRemainingLiquidCells;
    Kokkos::parallel_reduce(
        LocalDomainSize,
        KOKKOS_LAMBDA(const int &D3D1ConvPosition, unsigned long int &sum_superheated,
                      unsigned long int &sum_undercooled, unsigned long int &sum_active, unsigned long int &sum_solid,
                      unsigned long int &sum_remaining_liquid) {
            if (LayerID(D3D1ConvPosition) == layernumber) {
                if (CellType(D3D1ConvPosition) == Liquid) {
                    if (CritTimeStep(D3D1ConvPosition) > cycle)
                        sum_superheated += 1;
                    else
                        sum_undercooled += 1;
                }
                else if (CellType(D3D1ConvPosition) == Active)
                    sum_active += 1;
                else if (CellType(D3D1ConvPosition) == Solid)
                    sum_solid += 1;
            }
            else {
                if (CellType(D3D1ConvPosition) == Liquid)
                    sum_remaining_liquid += 1;
            }
        },
        LocalSuperheatedCells, LocalUndercooledCells, LocalActiveCells, LocalSolidCells, LocalRemainingLiquidCells);

    unsigned long int Global_SuccessfulNucEvents_ThisRank = 0;
    unsigned long int GlobalSuperheatedCells, GlobalUndercooledCells, GlobalActiveCells, GlobalSolidCells,
        GlobalRemainingLiquidCells;
    MPI_Reduce(&LocalSuperheatedCells, &GlobalSuperheatedCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalUndercooledCells, &GlobalUndercooledCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalActiveCells, &GlobalActiveCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalSolidCells, &GlobalSolidCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalRemainingLiquidCells, &GlobalRemainingLiquidCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&SuccessfulNucEvents_ThisRank, &Global_SuccessfulNucEvents_ThisRank, 1, MPI_INT, MPI_SUM, 0,
               MPI_COMM_WORLD);

    if (id == 0) {
        std::cout << "cycle = " << cycle << " Superheated liquid cells = " << GlobalSuperheatedCells
                  << " Undercooled liquid cells = " << GlobalUndercooledCells
                  << " Number of nucleation events this layer " << Global_SuccessfulNucEvents_ThisRank
                  << " Remaining liquid cells in future layers of this simulation = " << GlobalRemainingLiquidCells
                  << std::endl;
        if (GlobalSuperheatedCells + GlobalUndercooledCells == 0)
            XSwitch = 1;
    }
    MPI_Bcast(&XSwitch, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // If an appropraite problem type/solidification is not finished, jump to the next time step with work to be done,
    // if nothing left to do in the near future
    if ((XSwitch == 0) && ((TemperatureDataType == "R") || (TemperatureDataType == "S")))
        JumpTimeStep(cycle, GlobalUndercooledCells, LocalSuperheatedCells, CritTimeStep, LocalActiveDomainSize,
                     MyYSlices, ZBound_Low, false, CellType, LayerID, id, layernumber, np, nx, ny, nz, MyYOffset,
                     GrainID, CritTimeStep, GrainUnitVector, UndercoolingChange, UndercoolingCurrent, OutputFile,
                     NGrainOrientations, PathToOutput, IntermediateFileCounter, nzActive, deltax, XMin, YMin, ZMin,
                     NumberOfLayers, XSwitch, TemperatureDataType, PrintIdleMovieFrames, MovieFrameInc, PrintBinary,
                     FinishTimeStep[layernumber]);
}

//*****************************************************************************/
// Prints intermediate code output to stdout (intermediate output collected and printed is different than without
// remelting) and checks to see if solidification is complete in the case where cells can solidify multiple times
void IntermediateOutputAndCheck_Remelt(
    int id, int np, int &cycle, int MyYSlices, int MyYOffset, int LocalActiveDomainSize, int nx, int ny, int nz,
    int nzActive, double deltax, float XMin, float YMin, float ZMin, int SuccessfulNucEvents_ThisRank, int &XSwitch,
    ViewI CellType, ViewI CritTimeStep, ViewI GrainID, std::string TemperatureDataType, int layernumber, int,
    int ZBound_Low, int NGrainOrientations, ViewI LayerID, ViewF GrainUnitVector, ViewF UndercoolingChange,
    ViewF UndercoolingCurrent, std::string PathToOutput, std::string OutputFile, bool PrintIdleMovieFrames,
    int MovieFrameInc, int &IntermediateFileCounter, int NumberOfLayers, ViewI MeltTimeStep, bool PrintBinary) {

    unsigned long int LocalSuperheatedCells;
    unsigned long int LocalUndercooledCells;
    unsigned long int LocalActiveCells;
    unsigned long int LocalTempSolidCells;
    unsigned long int LocalFinishedSolidCells;
    Kokkos::parallel_reduce(
        LocalActiveDomainSize,
        KOKKOS_LAMBDA(const int &D3D1ConvPosition, unsigned long int &sum_superheated,
                      unsigned long int &sum_undercooled, unsigned long int &sum_active,
                      unsigned long int &sum_temp_solid, unsigned long int &sum_finished_solid) {
            int GlobalD3D1ConvPosition = D3D1ConvPosition + ZBound_Low * nx * MyYSlices;
            if (CellType(GlobalD3D1ConvPosition) == Liquid) {
                if (CritTimeStep(GlobalD3D1ConvPosition) > cycle)
                    sum_superheated += 1;
                else
                    sum_undercooled += 1;
            }
            else if (CellType(GlobalD3D1ConvPosition) == Active)
                sum_active += 1;
            else if (CellType(GlobalD3D1ConvPosition) == TempSolid)
                sum_temp_solid += 1;
            else if (CellType(GlobalD3D1ConvPosition) == Solid)
                sum_finished_solid += 1;
        },
        LocalSuperheatedCells, LocalUndercooledCells, LocalActiveCells, LocalTempSolidCells, LocalFinishedSolidCells);

    unsigned long int Global_SuccessfulNucEvents_ThisRank = 0;
    unsigned long int GlobalSuperheatedCells, GlobalUndercooledCells, GlobalActiveCells, GlobalTempSolidCells,
        GlobalFinishedSolidCells;
    MPI_Reduce(&LocalSuperheatedCells, &GlobalSuperheatedCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalUndercooledCells, &GlobalUndercooledCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalActiveCells, &GlobalActiveCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalTempSolidCells, &GlobalTempSolidCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalFinishedSolidCells, &GlobalFinishedSolidCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&SuccessfulNucEvents_ThisRank, &Global_SuccessfulNucEvents_ThisRank, 1, MPI_INT, MPI_SUM, 0,
               MPI_COMM_WORLD);

    if (id == 0) {
        std::cout << "cycle = " << cycle << " in layer " << layernumber
                  << ": Superheated liquid cells = " << GlobalSuperheatedCells
                  << " Undercooled liquid cells = " << GlobalUndercooledCells
                  << " Number of nucleation events this layer " << Global_SuccessfulNucEvents_ThisRank
                  << " cells to undergo at least one more solidification event = " << GlobalTempSolidCells
                  << " cells that are finished with solidification = " << GlobalFinishedSolidCells << std::endl;
        if (GlobalSuperheatedCells + GlobalUndercooledCells + GlobalActiveCells + GlobalTempSolidCells == 0)
            XSwitch = 1;
    }
    MPI_Bcast(&XSwitch, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if ((XSwitch == 0) && ((TemperatureDataType == "R") || (TemperatureDataType == "S")))
        JumpTimeStep(cycle, GlobalActiveCells, LocalTempSolidCells, MeltTimeStep, LocalActiveDomainSize, MyYSlices,
                     ZBound_Low, true, CellType, LayerID, id, layernumber, np, nx, ny, nz, MyYOffset, GrainID,
                     CritTimeStep, GrainUnitVector, UndercoolingChange, UndercoolingCurrent, OutputFile,
                     NGrainOrientations, PathToOutput, IntermediateFileCounter, nzActive, deltax, XMin, YMin, ZMin,
                     NumberOfLayers, XSwitch, TemperatureDataType, PrintIdleMovieFrames, MovieFrameInc, PrintBinary);
}
