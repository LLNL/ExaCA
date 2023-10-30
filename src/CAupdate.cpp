// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAupdate.hpp"

#include "mpi.h"

#include <cmath>

// Using for compatibility with device math functions.
using std::max;
using std::min;

//*****************************************************************************/
// Jump to the next time step with work to be done, if nothing left to do in the near future
// The cells of interest are active cells, and the view checked for future work is
// MeltTimeStep Print intermediate output during this jump if PrintIdleMovieFrames = true
void JumpTimeStep(int &cycle, unsigned long int RemainingCellsOfInterest, unsigned long int LocalTempSolidCells,
                  Temperature<device_memory_space> &temperature, Grid &grid, CellData<device_memory_space> &cellData,
                  int id, int layernumber, int np, ViewF GrainUnitVector, Print print, int NGrainOrientations) {

    auto CellType = cellData.getCellTypeSubview(grid);
    MPI_Bcast(&RemainingCellsOfInterest, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    if (RemainingCellsOfInterest == 0) {
        // If this rank still has cells that will later undergo transformation (LocalIncompleteCells > 0), check when
        // the next solid cells go above the liquidus (remelting) Otherwise, assign the largest possible time step as
        // the next time work needs to be done on the rank
        unsigned long int NextMeltTimeStep;
        if (LocalTempSolidCells > 0) {
            Kokkos::parallel_reduce(
                "CheckNextTSForWork", grid.DomainSize,
                KOKKOS_LAMBDA(const int &index, unsigned long int &tempv) {
                    // criteria for a cell to be associated with future work
                    if (CellType(index) == TempSolid) {
                        int SolidificationCounter_ThisCell = temperature.SolidificationEventCounter(index);
                        unsigned long int NextMeltTimeStep_ThisCell = static_cast<unsigned long int>(
                            temperature.LayerTimeTempHistory(index, SolidificationCounter_ThisCell, 0));
                        if (NextMeltTimeStep_ThisCell < tempv)
                            tempv = NextMeltTimeStep_ThisCell;
                    }
                },
                Kokkos::Min<unsigned long int>(NextMeltTimeStep));
        }
        else
            NextMeltTimeStep = INT_MAX;

        unsigned long int GlobalNextMeltTimeStep;
        MPI_Allreduce(&NextMeltTimeStep, &GlobalNextMeltTimeStep, 1, MPI_UNSIGNED_LONG, MPI_MIN, MPI_COMM_WORLD);
        if ((GlobalNextMeltTimeStep - cycle) > 5000) {
            // Print current grain misorientations (up to and including the current layer's data) for any of the time
            // steps between now and when melting/solidification occurs again, if the print option for idle frame
            // printing was toggled
            print.printIdleIntermediateGrainMisorientation(id, np, cycle, grid, cellData.GrainID_AllLayers,
                                                           cellData.CellType_AllLayers, GrainUnitVector,
                                                           NGrainOrientations, layernumber, GlobalNextMeltTimeStep);
            // Jump to next time step when solidification starts again
            cycle = GlobalNextMeltTimeStep - 1;
            if (id == 0)
                std::cout << "Jumping to cycle " << cycle + 1 << std::endl;
        }
    }
}

//*****************************************************************************/
// Prints intermediate code output to stdout (intermediate output collected and printed is different than without
// remelting) and checks to see if solidification is complete in the case where cells can solidify multiple times
void IntermediateOutputAndCheck(int id, int np, int &cycle, Grid &grid, int SuccessfulNucEvents_ThisRank, int &XSwitch,
                                CellData<device_memory_space> &cellData, Temperature<device_memory_space> &temperature,
                                std::string SimulationType, int layernumber, int NGrainOrientations,
                                ViewF GrainUnitVector, Print print) {

    auto CellType = cellData.getCellTypeSubview(grid);
    auto GrainID = cellData.getGrainIDSubview(grid);
    auto LayerID = cellData.getLayerIDSubview(grid);
    unsigned long int LocalSuperheatedCells;
    unsigned long int LocalUndercooledCells;
    unsigned long int LocalActiveCells;
    unsigned long int LocalTempSolidCells;
    unsigned long int LocalFinishedSolidCells;
    Kokkos::parallel_reduce(
        "IntermediateOutput", grid.DomainSize,
        KOKKOS_LAMBDA(const int &index, unsigned long int &sum_superheated, unsigned long int &sum_undercooled,
                      unsigned long int &sum_active, unsigned long int &sum_temp_solid,
                      unsigned long int &sum_finished_solid) {
            if (CellType(index) == Liquid) {
                int CritTimeStep = temperature.getCritTimeStep(index);
                if (CritTimeStep > cycle)
                    sum_superheated += 1;
                else
                    sum_undercooled += 1;
            }
            else if (CellType(index) == Active)
                sum_active += 1;
            else if (CellType(index) == TempSolid)
                sum_temp_solid += 1;
            else if (CellType(index) == Solid)
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
    // Cells of interest are those currently undergoing a melting-solidification cycle
    unsigned long int RemainingCellsOfInterest = GlobalActiveCells + GlobalSuperheatedCells + GlobalUndercooledCells;
    if ((XSwitch == 0) && ((SimulationType == "R") || (SimulationType == "S")))
        JumpTimeStep(cycle, RemainingCellsOfInterest, LocalTempSolidCells, temperature, grid, cellData, id, layernumber,
                     np, GrainUnitVector, print, NGrainOrientations);
}

//*****************************************************************************/
// Prints intermediate code output to stdout and checks to see the single grain simulation end condition (the grain has
// reached a domain edge) has been satisfied
void IntermediateOutputAndCheck(int id, int cycle, Grid &grid, int &XSwitch, ViewI CellType_AllLayers) {

    unsigned long int LocalLiquidCells, LocalActiveCells, LocalSolidCells;
    ViewB2D EdgesReached(Kokkos::ViewAllocateWithoutInitializing("EdgesReached"), 3, 2); // init to false
    Kokkos::deep_copy(EdgesReached, false);

    Kokkos::parallel_reduce(
        "IntermediateOutput", grid.DomainSize,
        KOKKOS_LAMBDA(const int &index, unsigned long int &sum_liquid, unsigned long int &sum_active,
                      unsigned long int &sum_solid) {
            if (CellType_AllLayers(index) == Liquid)
                sum_liquid += 1;
            else if (CellType_AllLayers(index) == Active) {
                sum_active += 1;
                // Did this cell reach a domain edge?
                int coord_x = grid.getCoordX(index);
                int coord_y = grid.getCoordY(index);
                int coord_z = grid.getCoordZ(index);
                int coord_y_global = coord_y + grid.y_offset;
                if (coord_x == 0)
                    EdgesReached(0, 0) = true;
                if (coord_x == grid.nx - 1)
                    EdgesReached(0, 1) = true;
                if (coord_y_global == 0)
                    EdgesReached(1, 0) = true;
                if (coord_y_global == grid.ny - 1)
                    EdgesReached(1, 1) = true;
                if (coord_z == 0)
                    EdgesReached(2, 0) = true;
                if (coord_z == grid.nz - 1)
                    EdgesReached(2, 1) = true;
            }
            else if (CellType_AllLayers(index) == Solid)
                sum_solid += 1;
        },
        LocalLiquidCells, LocalActiveCells, LocalSolidCells);

    unsigned long int GlobalLiquidCells, GlobalActiveCells, GlobalSolidCells;
    MPI_Reduce(&LocalLiquidCells, &GlobalLiquidCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalActiveCells, &GlobalActiveCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalSolidCells, &GlobalSolidCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "cycle = " << cycle << " : Liquid cells = " << GlobalLiquidCells
                  << " Active cells = " << GlobalActiveCells << " Solid cells = " << GlobalSolidCells << std::endl;

    // Each rank checks to see if a global domain boundary was reached
    ViewB2D_H EdgesReached_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), EdgesReached);
    int XSwitchLocal = 0;
    std::vector<std::string> EdgeDims = {"X", "Y", "Z"};
    std::vector<std::string> EdgeNames = {"Lower", "Upper"};
    for (int edgedim = 0; edgedim < 3; edgedim++) {
        for (int edgename = 0; edgename < 2; edgename++) {
            if (EdgesReached_Host(edgedim, edgename)) {
                std::cout << EdgeNames[edgename] << " edge of domain in the " << EdgeDims[edgedim]
                          << " direction was reached on rank " << id << " and cycle " << cycle
                          << "; simulation is complete" << std::endl;
                XSwitchLocal = 1;
            }
        }
    }
    // Simulation ends if a global domain boundary was reached on any rank
    MPI_Allreduce(&XSwitchLocal, &XSwitch, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
}
