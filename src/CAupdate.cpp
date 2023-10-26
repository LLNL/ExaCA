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
// Determine which cells are associated with the "steering vector" of cells that are either active, or becoming active
// this time step
void FillSteeringVector_NoRemelt(int cycle, Grid &grid, Temperature<device_memory_space> &temperature,
                                 CellData<device_memory_space> &cellData, ViewI SteeringVector, ViewI numSteer,
                                 ViewI_H numSteer_Host) {

    // Cells associated with this layer that are not solid type but have passed the liquidus (crit time step) have their
    // undercooling values updated Cells that meet the aforementioned criteria and are active type should be added to
    // the steering vector
    auto CellType = cellData.getCellTypeSubview(grid);
    Kokkos::parallel_for(
        "FillSV", grid.DomainSize, KOKKOS_LAMBDA(const int &index) {
            int cellType = CellType(index);
            int isNotSolid = (cellType != Solid);
            int CritTimeStep = temperature.LayerTimeTempHistory(index, 0, 1);
            int pastCritTime = (cycle > CritTimeStep);
            int cell_Active = ((cellType == Active) || (cellType == FutureActive));
            if (isNotSolid && pastCritTime) {
                temperature.update_undercooling(index);
                if (cell_Active) {
                    SteeringVector(Kokkos::atomic_fetch_add(&numSteer(0), 1)) = index;
                }
            }
        });
    Kokkos::deep_copy(numSteer_Host, numSteer);
}

//*****************************************************************************/
// Determine which cells are associated with the "steering vector" of cells that are either active, or becoming active
// this time step - version with remelting
void FillSteeringVector_Remelt(int cycle, Grid &grid, Temperature<device_memory_space> &temperature,
                               CellData<device_memory_space> &cellData, ViewI SteeringVector, ViewI numSteer,
                               ViewI_H numSteer_Host) {

    auto CellType = cellData.getCellTypeSubview(grid);
    auto GrainID = cellData.getGrainIDSubview(grid);
    Kokkos::parallel_for(
        "FillSV_RM", grid.DomainSize, KOKKOS_LAMBDA(const int &index) {
            int cellType = CellType(index);
            // Only iterate over cells that are not Solid type
            if (cellType != Solid) {
                int MeltTimeStep = temperature.getMeltTimeStep(cycle, index);
                int CritTimeStep = temperature.getCritTimeStep(index);
                bool atMeltTime = (cycle == MeltTimeStep);
                bool atCritTime = (cycle == CritTimeStep);
                bool pastCritTime = (cycle > CritTimeStep);
                if (atMeltTime) {
                    // Cell melts, undercooling is reset to 0 from the previous value, if any
                    CellType(index) = Liquid;
                    temperature.reset_undercooling(index);
                    if (cellType != TempSolid) {
                        // This cell either hasn't started or hasn't finished the previous solidification event, but has
                        // undergone melting - increment the solidification counter to move on to the next
                        // melt-solidification event
                        temperature.update_solidification_counter(index);
                    }
                    // Any adjacent active cells should also be remelted, as these cells are more likely heating up than
                    // cooling down These are converted to the temporary FutureLiquid state, to be later iterated over
                    // and loaded into the steering vector as necessary Get the x, y, z coordinates of the cell on this
                    // MPI rank
                    int coord_x = grid.getCoordX(index);
                    int coord_y = grid.getCoordY(index);
                    int coord_z = grid.getCoordZ(index);
                    for (int l = 0; l < 26; l++) {
                        // "l" correpsponds to the specific neighboring cell
                        // Local coordinates of adjacent cell center
                        int neighbor_coord_x = coord_x + grid.NeighborX[l];
                        int neighbor_coord_y = coord_y + grid.NeighborY[l];
                        int neighbor_coord_z = coord_z + grid.NeighborZ[l];
                        if ((neighbor_coord_x >= 0) && (neighbor_coord_x < grid.nx) && (neighbor_coord_y >= 0) &&
                            (neighbor_coord_y < grid.ny_local) && (neighbor_coord_z < grid.nz_layer) &&
                            (neighbor_coord_z >= 0)) {
                            int neighbor_index = grid.get1Dindex(neighbor_coord_x, neighbor_coord_y, neighbor_coord_z);
                            if (CellType(neighbor_index) == Active) {
                                // Mark adjacent active cells to this as cells that should be converted into liquid, as
                                // they are more likely heating than cooling
                                CellType(neighbor_index) = FutureLiquid;
                                SteeringVector(Kokkos::atomic_fetch_add(&numSteer(0), 1)) = neighbor_index;
                            }
                        }
                    }
                }
                else if ((cellType != TempSolid) && (pastCritTime)) {
                    // Update cell undercooling
                    temperature.update_undercooling(index);
                    if (cellType == Active) {
                        // Add active cells below liquidus to steering vector
                        SteeringVector(Kokkos::atomic_fetch_add(&numSteer(0), 1)) = index;
                    }
                }
                else if ((atCritTime) && (cellType == Liquid) && (GrainID(index) != 0)) {
                    // Get the x, y, z coordinates of the cell on this MPI rank
                    int coord_x = grid.getCoordX(index);
                    int coord_y = grid.getCoordY(index);
                    int coord_z = grid.getCoordZ(index);
                    // If this cell has cooled to the liquidus temperature, borders at least one solid/tempsolid cell,
                    // and is part of a grain, it should become active. This only needs to be checked on the time step
                    // where the cell reaches the liquidus, not every time step beyond this
                    for (int l = 0; l < 26; l++) {
                        // "l" correpsponds to the specific neighboring cell
                        // Local coordinates of adjacent cell center
                        int neighbor_coord_x = coord_x + grid.NeighborX[l];
                        int neighbor_coord_y = coord_y + grid.NeighborY[l];
                        int neighbor_coord_z = coord_z + grid.NeighborZ[l];
                        if ((neighbor_coord_x >= 0) && (neighbor_coord_x < grid.nx) && (neighbor_coord_y >= 0) &&
                            (neighbor_coord_y < grid.ny_local) && (neighbor_coord_z < grid.nz_layer) &&
                            (neighbor_coord_z >= 0)) {
                            int neighbor_index = grid.get1Dindex(neighbor_coord_x, neighbor_coord_y, neighbor_coord_z);
                            if ((CellType(neighbor_index) == TempSolid) || (CellType(neighbor_index) == Solid) ||
                                (coord_z == 0)) {
                                // Cell activation to be performed as part of steering vector
                                l = 26;
                                SteeringVector(Kokkos::atomic_fetch_add(&numSteer(0), 1)) = index;
                                CellType(index) = FutureActive; // this cell cannot be captured - is being activated
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

// Decentered octahedron algorithm for the capture of new interface cells by grains
void CellCapture(int, int np, int, Grid &grid, InterfacialResponseFunction irf, ViewF GrainUnitVector,
                 ViewF CritDiagonalLength, ViewF DiagonalLength, CellData<device_memory_space> &cellData,
                 Temperature<device_memory_space> &temperature, ViewF DOCenter, int NGrainOrientations,
                 Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, ViewI SendSizeNorth, ViewI SendSizeSouth,
                 ViewI SteeringVector, ViewI numSteer, ViewI_H numSteer_Host, int &BufSize) {

    // Loop over list of active and soon-to-be active cells, potentially performing cell capture events and updating
    // cell types
    auto CellType = cellData.getCellTypeSubview(grid);
    auto GrainID = cellData.getGrainIDSubview(grid);
    Kokkos::parallel_for(
        "CellCapture", numSteer_Host(0), KOKKOS_LAMBDA(const int &num) {
            numSteer(0) = 0;
            int index = SteeringVector(num);
            // Get the x, y, z coordinates of the cell on this MPI rank
            int coord_x = grid.getCoordX(index);
            int coord_y = grid.getCoordY(index);
            int coord_z = grid.getCoordZ(index);
            // Cells of interest for the CA - active cells and future active/liquid cells
            if (CellType(index) == Active) {
                // Update local diagonal length of active cell
                double LocU = temperature.UndercoolingCurrent(index);
                LocU = min(210.0, LocU);
                double V = irf.compute(LocU);
                DiagonalLength(index) += min(0.045, V); // Max amount the diagonal can grow per time step
                // Cycle through all neigboring cells on this processor to see if they have been captured
                // Cells in ghost nodes cannot capture cells on other processors
                bool DeactivateCell = true; // switch that becomes false if the cell has at least 1 liquid type neighbor
                // Which neighbors should be iterated over?
                for (int l = 0; l < 26; l++) {
                    // Local coordinates of adjacent cell center
                    int neighbor_coord_x = coord_x + grid.NeighborX[l];
                    int neighbor_coord_y = coord_y + grid.NeighborY[l];
                    int neighbor_coord_z = coord_z + grid.NeighborZ[l];
                    // Check if neighbor is in bounds
                    if ((neighbor_coord_x >= 0) && (neighbor_coord_x < grid.nx) && (neighbor_coord_y >= 0) &&
                        (neighbor_coord_y < grid.ny_local) && (neighbor_coord_z < grid.nz_layer) &&
                        (neighbor_coord_z >= 0)) {
                        int neighbor_index = grid.get1Dindex(neighbor_coord_x, neighbor_coord_y, neighbor_coord_z);
                        if (CellType(neighbor_index) == Liquid)
                            DeactivateCell = false;
                        // Capture of cell located at "NeighborD3D1ConvPosition" if this condition is satisfied
                        if ((DiagonalLength(index) >= CritDiagonalLength(26 * index + l)) &&
                            (CellType(neighbor_index) == Liquid)) {
                            // Use of atomic_compare_exchange
                            // (https://github.com/kokkos/kokkos/wiki/Kokkos%3A%3Aatomic_compare_exchange) old_val =
                            // atomic_compare_exchange(ptr_to_value,comparison_value, new_value); Atomicly sets the
                            // value at the address given by ptr_to_value to new_value if the current value at
                            // ptr_to_value is equal to comparison_value Returns the previously stored value at the
                            // address independent on whether the exchange has happened. If this cell's is a liquid
                            // cell, change it to "TemporaryUpdate" type and return a value of "liquid" If this cell has
                            // already been changed to "TemporaryUpdate" type, return a value of "0"
                            int update_val = TemporaryUpdate;
                            int old_val = Liquid;
                            int OldCellTypeValue =
                                Kokkos::atomic_compare_exchange(&CellType(neighbor_index), old_val, update_val);
                            // Only proceed if CellType was previously liquid (this current thread changed the value to
                            // TemporaryUpdate)
                            if (OldCellTypeValue == Liquid) {
                                int h = GrainID(index);
                                int MyOrientation = getGrainOrientation(h, NGrainOrientations);

                                // The new cell is captured by this cell's growing octahedron (Grain "h")
                                GrainID(neighbor_index) = h;

                                // (cxold, cyold, czold) are the coordiantes of this decentered octahedron
                                float cxold = DOCenter(3 * index);
                                float cyold = DOCenter(3 * index + 1);
                                float czold = DOCenter(3 * index + 2);

                                // (xp,yp,zp) are the global coordinates of the new cell's center
                                // Note that the Y coordinate is relative to the domain origin to keep the coordinate
                                // system continuous across ranks
                                float xp = coord_x + grid.NeighborX[l] + 0.5;
                                float yp = coord_y + grid.y_offset + grid.NeighborY[l] + 0.5;
                                float zp = coord_z + grid.NeighborZ[l] + 0.5;

                                // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the image
                                // of the center of the new cell
                                float x0 = xp - cxold;
                                float y0 = yp - cyold;
                                float z0 = zp - czold;

                                // mag0 is the magnitude of (x0,y0,z0)
                                float mag0 = sqrtf(x0 * x0 + y0 * y0 + z0 * z0);

                                // Calculate unit vectors for the octahedron that intersect the new cell center
                                float Angle1 = (GrainUnitVector(9 * MyOrientation) * x0 +
                                                GrainUnitVector(9 * MyOrientation + 1) * y0 +
                                                GrainUnitVector(9 * MyOrientation + 2) * z0) /
                                               mag0;
                                float Angle2 = (GrainUnitVector(9 * MyOrientation + 3) * x0 +
                                                GrainUnitVector(9 * MyOrientation + 4) * y0 +
                                                GrainUnitVector(9 * MyOrientation + 5) * z0) /
                                               mag0;
                                float Angle3 = (GrainUnitVector(9 * MyOrientation + 6) * x0 +
                                                GrainUnitVector(9 * MyOrientation + 7) * y0 +
                                                GrainUnitVector(9 * MyOrientation + 8) * z0) /
                                               mag0;
                                float Diag1X = GrainUnitVector(9 * MyOrientation) * (2 * (Angle1 < 0) - 1);
                                float Diag1Y = GrainUnitVector(9 * MyOrientation + 1) * (2 * (Angle1 < 0) - 1);
                                float Diag1Z = GrainUnitVector(9 * MyOrientation + 2) * (2 * (Angle1 < 0) - 1);

                                float Diag2X = GrainUnitVector(9 * MyOrientation + 3) * (2 * (Angle2 < 0) - 1);
                                float Diag2Y = GrainUnitVector(9 * MyOrientation + 4) * (2 * (Angle2 < 0) - 1);
                                float Diag2Z = GrainUnitVector(9 * MyOrientation + 5) * (2 * (Angle2 < 0) - 1);

                                float Diag3X = GrainUnitVector(9 * MyOrientation + 6) * (2 * (Angle3 < 0) - 1);
                                float Diag3Y = GrainUnitVector(9 * MyOrientation + 7) * (2 * (Angle3 < 0) - 1);
                                float Diag3Z = GrainUnitVector(9 * MyOrientation + 8) * (2 * (Angle3 < 0) - 1);

                                float U1[3], U2[3];
                                U1[0] = Diag2X - Diag1X;
                                U1[1] = Diag2Y - Diag1Y;
                                U1[2] = Diag2Z - Diag1Z;
                                U2[0] = Diag3X - Diag1X;
                                U2[1] = Diag3Y - Diag1Y;
                                U2[2] = Diag3Z - Diag1Z;
                                float UU[3];
                                UU[0] = U1[1] * U2[2] - U1[2] * U2[1];
                                UU[1] = U1[2] * U2[0] - U1[0] * U2[2];
                                UU[2] = U1[0] * U2[1] - U1[1] * U2[0];
                                float NDem = sqrtf(UU[0] * UU[0] + UU[1] * UU[1] + UU[2] * UU[2]);
                                float Norm[3];
                                Norm[0] = UU[0] / NDem;
                                Norm[1] = UU[1] / NDem;
                                Norm[2] = UU[2] / NDem;
                                // normal to capturing plane
                                double norm[3], TriangleX[3], TriangleY[3], TriangleZ[3], ParaT;
                                norm[0] = Norm[0];
                                norm[1] = Norm[1];
                                norm[2] = Norm[2];
                                ParaT = (norm[0] * x0 + norm[1] * y0 + norm[2] * z0) /
                                        (norm[0] * Diag1X + norm[1] * Diag1Y + norm[2] * Diag1Z);

                                TriangleX[0] = cxold + ParaT * Diag1X;
                                TriangleY[0] = cyold + ParaT * Diag1Y;
                                TriangleZ[0] = czold + ParaT * Diag1Z;

                                TriangleX[1] = cxold + ParaT * Diag2X;
                                TriangleY[1] = cyold + ParaT * Diag2Y;
                                TriangleZ[1] = czold + ParaT * Diag2Z;

                                TriangleX[2] = cxold + ParaT * Diag3X;
                                TriangleY[2] = cyold + ParaT * Diag3Y;
                                TriangleZ[2] = czold + ParaT * Diag3Z;

                                // Determine which of the 3 corners of the capturing face is closest to the captured
                                // cell center
                                float DistToCorner[3];
                                DistToCorner[0] = sqrtf(((TriangleX[0] - xp) * (TriangleX[0] - xp)) +
                                                        ((TriangleY[0] - yp) * (TriangleY[0] - yp)) +
                                                        ((TriangleZ[0] - zp) * (TriangleZ[0] - zp)));
                                DistToCorner[1] = sqrtf(((TriangleX[1] - xp) * (TriangleX[1] - xp)) +
                                                        ((TriangleY[1] - yp) * (TriangleY[1] - yp)) +
                                                        ((TriangleZ[1] - zp) * (TriangleZ[1] - zp)));
                                DistToCorner[2] = sqrtf(((TriangleX[2] - xp) * (TriangleX[2] - xp)) +
                                                        ((TriangleY[2] - yp) * (TriangleY[2] - yp)) +
                                                        ((TriangleZ[2] - zp) * (TriangleZ[2] - zp)));

                                int x, y, z;
                                x = (DistToCorner[0] < DistToCorner[1]);
                                y = (DistToCorner[1] < DistToCorner[2]);
                                z = (DistToCorner[2] < DistToCorner[0]);

                                int idx = 2 * (z - y) * z + (y - x) * y;
                                float mindisttocorner = DistToCorner[idx];
                                float xc = TriangleX[idx];
                                float yc = TriangleY[idx];
                                float zc = TriangleZ[idx];

                                float x1 = TriangleX[(idx + 1) % 3];
                                float y1 = TriangleY[(idx + 1) % 3];
                                float z1 = TriangleZ[(idx + 1) % 3];
                                float x2 = TriangleX[(idx + 2) % 3];
                                float y2 = TriangleY[(idx + 2) % 3];
                                float z2 = TriangleZ[(idx + 2) % 3];

                                float D1 =
                                    sqrtf(((xp - x2) * (xp - x2)) + ((yp - y2) * (yp - y2)) + ((zp - z2) * (zp - z2)));
                                float D2 =
                                    sqrtf(((xc - x2) * (xc - x2)) + ((yc - y2) * (yc - y2)) + ((zc - z2) * (zc - z2)));
                                float D3 =
                                    sqrtf(((xp - x1) * (xp - x1)) + ((yp - y1) * (yp - y1)) + ((zp - z1) * (zp - z1)));
                                float D4 =
                                    sqrtf(((xc - x1) * (xc - x1)) + ((yc - y1) * (yc - y1)) + ((zc - z1) * (zc - z1)));

                                float I1 = 0;
                                float I2 = D2;
                                float J1 = 0;
                                float J2 = D4;
                                // If minimum distance to corner = 0, the octahedron corner captured the new cell center
                                if (mindisttocorner != 0) {
                                    I1 = D1 * ((xp - x2) * (xc - x2) + (yp - y2) * (yc - y2) + (zp - z2) * (zc - z2)) /
                                         (D1 * D2);
                                    I2 = D2 - I1;
                                    J1 = D3 * ((xp - x1) * (xc - x1) + (yp - y1) * (yc - y1) + (zp - z1) * (zc - z1)) /
                                         (D3 * D4);
                                    J2 = D4 - J1;
                                }
                                float L12 = 0.5 * (fmin(I1, sqrtf(3.0)) + fmin(I2, sqrtf(3.0)));
                                float L13 = 0.5 * (fmin(J1, sqrtf(3.0)) + fmin(J2, sqrtf(3.0)));
                                float NewODiagL = sqrtf(2.0) * fmax(L12, L13); // half diagonal length of new octahedron

                                DiagonalLength(neighbor_index) = NewODiagL;
                                // Calculate coordinates of new decentered octahedron center
                                float CaptDiag[3];
                                CaptDiag[0] = xc - cxold;
                                CaptDiag[1] = yc - cyold;
                                CaptDiag[2] = zc - czold;

                                float CaptDiagMagnitude = sqrt(CaptDiag[0] * CaptDiag[0] + CaptDiag[1] * CaptDiag[1] +
                                                               CaptDiag[2] * CaptDiag[2]);
                                float CaptDiagUV[3];
                                CaptDiagUV[0] = CaptDiag[0] / CaptDiagMagnitude;
                                CaptDiagUV[1] = CaptDiag[1] / CaptDiagMagnitude;
                                CaptDiagUV[2] = CaptDiag[2] / CaptDiagMagnitude;
                                // (cx, cy, cz) are the coordiantes of the new active cell's decentered octahedron
                                float cx = xc - NewODiagL * CaptDiagUV[0];
                                float cy = yc - NewODiagL * CaptDiagUV[1];
                                float cz = zc - NewODiagL * CaptDiagUV[2];

                                DOCenter(3 * neighbor_index) = cx;
                                DOCenter(3 * neighbor_index + 1) = cy;
                                DOCenter(3 * neighbor_index + 2) = cz;

                                // Get new critical diagonal length values for the newly activated cell (at array
                                // position "neighbor_index")
                                calcCritDiagonalLength(neighbor_index, xp, yp, zp, cx, cy, cz, grid.NeighborX,
                                                       grid.NeighborY, grid.NeighborZ, MyOrientation, GrainUnitVector,
                                                       CritDiagonalLength);

                                if (np > 1) {
                                    // TODO: Test loading ghost nodes in a separate kernel, potentially adopting this
                                    // change if the slowdown is minor
                                    int GhostGID = h;
                                    float GhostDOCX = cx;
                                    float GhostDOCY = cy;
                                    float GhostDOCZ = cz;
                                    float GhostDL = NewODiagL;
                                    // Collect data for the ghost nodes, if necessary
                                    // Data loaded into the ghost nodes is for the cell that was just captured
                                    bool DataFitsInBuffer = loadghostnodes(
                                        GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth,
                                        SendSizeSouth, grid.ny_local, neighbor_coord_x, neighbor_coord_y,
                                        neighbor_coord_z, grid.AtNorthBoundary, grid.AtSouthBoundary, BufferSouthSend,
                                        BufferNorthSend, NGrainOrientations, BufSize);
                                    if (!(DataFitsInBuffer)) {
                                        // This cell's data did not fit in the buffer with current size BufSize - mark
                                        // with temporary type
                                        CellType(neighbor_index) = ActiveFailedBufferLoad;
                                    }
                                    else {
                                        // Cell activation is now finished - cell type can be changed from
                                        // TemporaryUpdate to Active
                                        CellType(neighbor_index) = Active;
                                    }
                                } // End if statement for serial/parallel code
                                else {
                                    // Only update the new cell's type once Critical Diagonal Length, Triangle Index,
                                    // and Diagonal Length values have been assigned to it Avoids the race condition in
                                    // which the new cell is activated, and another thread acts on the new active cell
                                    // before the cell's new critical diagonal length/triangle index/diagonal length
                                    // values are assigned
                                    CellType(neighbor_index) = Active;
                                }
                            } // End if statement within locked capture loop
                        }     // End if statement for outer capture loop
                    }         // End if statement over neighbors on the active grid
                }             // End loop over all neighbors of this active cell
                if (DeactivateCell) {
                    // This active cell has no more neighboring cells to be captured
                    // Update the counter for the number of times this cell went from liquid to active to solid
                    bool SolidificationCompleteYN = temperature.update_check_solidification_counter(index);
                    // Did the cell solidify for the last time in the layer?
                    // If so, this cell is solid - ignore until next layer (if needed)
                    // If not, this cell is tempsolid, will become liquid again
                    if (SolidificationCompleteYN)
                        CellType(index) = Solid;
                    else
                        CellType(index) = TempSolid;
                }
            }
            else if (CellType(index) == FutureActive) {
                // Successful nucleation event - this cell is becoming a new active cell
                CellType(index) = TemporaryUpdate; // avoid operating on the new active cell before its
                                                   // associated octahedron data is initialized
                int MyGrainID = GrainID(index);    // GrainID was assigned as part of Nucleation

                // Initialize new octahedron
                createNewOctahedron(index, DiagonalLength, DOCenter, coord_x, coord_y, grid.y_offset, coord_z);
                // The orientation for the new grain will depend on its Grain ID (nucleated grains have negative GrainID
                // values)
                int MyOrientation = getGrainOrientation(MyGrainID, NGrainOrientations);
                // Octahedron center is at (cx, cy, cz) - note that the Y coordinate is relative to the domain origin to
                // keep the coordinate system continuous across ranks
                float cx = coord_x + 0.5;
                float cy = coord_y + grid.y_offset + 0.5;
                float cz = coord_z + 0.5;
                // Calculate critical values at which this active cell leads to the activation of a neighboring liquid
                // cell. Octahedron center and cell center overlap for octahedra created as part of a new grain
                calcCritDiagonalLength(index, cx, cy, cz, cx, cy, cz, grid.NeighborX, grid.NeighborY, grid.NeighborZ,
                                       MyOrientation, GrainUnitVector, CritDiagonalLength);
                if (np > 1) {
                    // TODO: Test loading ghost nodes in a separate kernel, potentially adopting this change if the
                    // slowdown is minor
                    int GhostGID = MyGrainID;
                    float GhostDOCX = cx;
                    float GhostDOCY = cy;
                    float GhostDOCZ = cz;
                    float GhostDL = 0.01;
                    // Collect data for the ghost nodes, if necessary
                    bool DataFitsInBuffer = loadghostnodes(
                        GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth, SendSizeSouth, grid.ny_local,
                        coord_x, coord_y, coord_z, grid.AtNorthBoundary, grid.AtSouthBoundary, BufferSouthSend,
                        BufferNorthSend, NGrainOrientations, BufSize);
                    if (!(DataFitsInBuffer)) {
                        // This cell's data did not fit in the buffer with current size BufSize - mark with temporary
                        // type
                        CellType(index) = ActiveFailedBufferLoad;
                    }
                    else {
                        // Cell activation is now finished - cell type can be changed from TemporaryUpdate to Active
                        CellType(index) = Active;
                    }
                } // End if statement for serial/parallel code
                else {
                    // Cell activation is now finished - cell type can be changed from TemporaryUpdate to Active
                    CellType(index) = Active;
                } // End if statement for serial/parallel code
            }
            else if (CellType(index) == FutureLiquid) {
                // This type was assigned to a cell that was recently transformed from active to liquid, due to its
                // bordering of a cell above the liquidus. This information may need to be sent to other MPI ranks
                // Dummy values for first 4 arguments (Grain ID and octahedron center coordinates), 0 for diagonal
                // length
                bool DataFitsInBuffer =
                    loadghostnodes(-1, -1.0, -1.0, -1.0, 0.0, SendSizeNorth, SendSizeSouth, grid.ny_local, coord_x,
                                   coord_y, coord_z, grid.AtNorthBoundary, grid.AtSouthBoundary, BufferSouthSend,
                                   BufferNorthSend, NGrainOrientations, BufSize);
                if (!(DataFitsInBuffer)) {
                    // This cell's data did not fit in the buffer with current size BufSize - mark with temporary
                    // type
                    CellType(index) = LiquidFailedBufferLoad;
                }
                else {
                    // Cell activation is now finished - cell type can be changed from FutureLiquid to Active
                    CellType(index) = Liquid;
                }
            }
        });
    Kokkos::fence();
}

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
