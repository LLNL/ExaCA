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

// For the case where all cells solidify once, determine which cells are associated with the "steering vector" of
// cells that are either active, or becoming active this time step
void fill_steering_vector_no_remelt(const int cycle, const Grid &grid, CellData<device_memory_space> &cellData,
                                    Temperature<device_memory_space> &temperature,
                                    Interface<device_memory_space> &interface) {

    // Cells associated with this layer that are not solid type but have passed the liquidus (crit time step) have
    // their undercooling values updated Cells that meet the aforementioned criteria and are active type should be
    // added to the steering vector
    auto CellType = cellData.getCellTypeSubview(grid);
    Kokkos::parallel_for(
        "FillSV", grid.domain_size, KOKKOS_LAMBDA(const int &index) {
            int cellType = CellType(index);
            int isNotSolid = (cellType != Solid);
            int CritTimeStep = temperature.LayerTimeTempHistory(index, 0, 1);
            int pastCritTime = (cycle > CritTimeStep);
            int cell_Active = ((cellType == Active) || (cellType == FutureActive));
            if (isNotSolid && pastCritTime) {
                temperature.update_undercooling(index);
                if (cell_Active) {
                    interface.steering_vector(Kokkos::atomic_fetch_add(&interface.num_steer(0), 1)) = index;
                }
            }
        });
    Kokkos::deep_copy(interface.num_steer_host, interface.num_steer);
}

// For the case where cells may melt and solidify multiple times, determine which cells are associated with the
// "steering vector" of cells that are either active, or becoming active this time step - version with remelting
void fill_steering_vector_remelt(const int cycle, const Grid &grid, CellData<device_memory_space> &cellData,
                                 Temperature<device_memory_space> &temperature,
                                 Interface<device_memory_space> &interface) {

    auto CellType = cellData.getCellTypeSubview(grid);
    auto GrainID = cellData.getGrainIDSubview(grid);
    Kokkos::parallel_for(
        "FillSV_RM", grid.domain_size, KOKKOS_LAMBDA(const int &index) {
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
                        // This cell either hasn't started or hasn't finished the previous solidification event, but
                        // has undergone melting - increment the solidification counter to move on to the next
                        // melt-solidification event
                        temperature.update_solidification_counter(index);
                    }
                    // Any adjacent active cells should also be remelted, as these cells are more likely heating up
                    // than cooling down These are converted to the temporary FutureLiquid state, to be later
                    // iterated over and loaded into the steering vector as necessary Get the x, y, z coordinates of
                    // the cell on this MPI rank
                    int coord_x = grid.get_coord_X(index);
                    int coord_y = grid.get_coord_Y(index);
                    int coord_z = grid.get_coord_Z(index);
                    for (int l = 0; l < 26; l++) {
                        // "l" correpsponds to the specific neighboring cell
                        // Local coordinates of adjacent cell center
                        int neighbor_coord_x = coord_x + interface.neighbor_x[l];
                        int neighbor_coord_y = coord_y + interface.neighbor_y[l];
                        int neighbor_coord_z = coord_z + interface.neighbor_z[l];
                        if ((neighbor_coord_x >= 0) && (neighbor_coord_x < grid.nx) && (neighbor_coord_y >= 0) &&
                            (neighbor_coord_y < grid.ny_local) && (neighbor_coord_z < grid.nz_layer) &&
                            (neighbor_coord_z >= 0)) {
                            int neighbor_index =
                                grid.get_1D_index(neighbor_coord_x, neighbor_coord_y, neighbor_coord_z);
                            if (CellType(neighbor_index) == Active) {
                                // Mark adjacent active cells to this as cells that should be converted into liquid,
                                // as they are more likely heating than cooling
                                CellType(neighbor_index) = FutureLiquid;
                                interface.steering_vector(Kokkos::atomic_fetch_add(&interface.num_steer(0), 1)) =
                                    neighbor_index;
                            }
                        }
                    }
                }
                else if ((cellType != TempSolid) && (pastCritTime)) {
                    // Update cell undercooling
                    temperature.update_undercooling(index);
                    if (cellType == Active) {
                        // Add active cells below liquidus to steering vector
                        interface.steering_vector(Kokkos::atomic_fetch_add(&interface.num_steer(0), 1)) = index;
                    }
                }
                else if ((atCritTime) && (cellType == Liquid) && (GrainID(index) != 0)) {
                    // Get the x, y, z coordinates of the cell on this MPI rank
                    int coord_x = grid.get_coord_X(index);
                    int coord_y = grid.get_coord_Y(index);
                    int coord_z = grid.get_coord_Z(index);
                    // If this cell has cooled to the liquidus temperature, borders at least one solid/tempsolid
                    // cell, and is part of a grain, it should become active. This only needs to be checked on the
                    // time step where the cell reaches the liquidus, not every time step beyond this
                    for (int l = 0; l < 26; l++) {
                        // "l" correpsponds to the specific neighboring cell
                        // Local coordinates of adjacent cell center
                        int neighbor_coord_x = coord_x + interface.neighbor_x[l];
                        int neighbor_coord_y = coord_y + interface.neighbor_y[l];
                        int neighbor_coord_z = coord_z + interface.neighbor_z[l];
                        if ((neighbor_coord_x >= 0) && (neighbor_coord_x < grid.nx) && (neighbor_coord_y >= 0) &&
                            (neighbor_coord_y < grid.ny_local) && (neighbor_coord_z < grid.nz_layer) &&
                            (neighbor_coord_z >= 0)) {
                            int neighbor_index =
                                grid.get_1D_index(neighbor_coord_x, neighbor_coord_y, neighbor_coord_z);
                            if ((CellType(neighbor_index) == TempSolid) || (CellType(neighbor_index) == Solid) ||
                                (coord_z == 0)) {
                                // Cell activation to be performed as part of steering vector
                                l = 26;
                                interface.steering_vector(Kokkos::atomic_fetch_add(&interface.num_steer(0), 1)) = index;
                                CellType(index) = FutureActive; // this cell cannot be captured - is being activated
                            }
                        }
                    }
                }
            }
        });
    Kokkos::fence();

    // Copy size of steering vector (containing positions of undercooled liquid/active cells) to the host
    Kokkos::deep_copy(interface.num_steer_host, interface.num_steer);
}

// Decentered octahedron algorithm for the capture of new interface cells by grains
void cell_capture(const int, const int np, const Grid &grid, const InterfacialResponseFunction &irf,
                  CellData<device_memory_space> &cellData, Temperature<device_memory_space> &temperature,
                  Interface<device_memory_space> &interface, const ViewF GrainUnitVector,
                  const int NGrainOrientations) {

    // Get subviews for this layer
    auto CellType = cellData.getCellTypeSubview(grid);
    auto GrainID = cellData.getGrainIDSubview(grid);
    // Loop over list of active and soon-to-be active cells, potentially performing cell capture events and updating
    // cell types
    Kokkos::parallel_for(
        "CellCapture", interface.num_steer_host(0), KOKKOS_LAMBDA(const int &num) {
            interface.num_steer(0) = 0;
            int index = interface.steering_vector(num);
            // Get the x, y, z coordinates of the cell on this MPI rank
            int coord_x = grid.get_coord_X(index);
            int coord_y = grid.get_coord_Y(index);
            int coord_z = grid.get_coord_Z(index);

            // Cells of interest for the CA - active cells and future active/liquid cells
            if (CellType(index) == Active) {
                // Update local diagonal length of active cell
                double LocU = temperature.UndercoolingCurrent(index);
                LocU = min(210.0, LocU);
                double V = irf.compute(LocU);
                interface.diagonal_length(index) += min(0.045, V); // Max amount the diagonal can grow per time step
                // Cycle through all neigboring cells on this processor to see if they have been captured
                // Cells in ghost nodes cannot capture cells on other processors
                bool DeactivateCell = true; // switch that becomes false if the cell has at least 1 liquid type neighbor
                // Which neighbors should be iterated over?
                for (int l = 0; l < 26; l++) {
                    // Local coordinates of adjacent cell center
                    int neighbor_coord_x = coord_x + interface.neighbor_x[l];
                    int neighbor_coord_y = coord_y + interface.neighbor_y[l];
                    int neighbor_coord_z = coord_z + interface.neighbor_z[l];
                    // Check if neighbor is in bounds
                    if ((neighbor_coord_x >= 0) && (neighbor_coord_x < grid.nx) && (neighbor_coord_y >= 0) &&
                        (neighbor_coord_y < grid.ny_local) && (neighbor_coord_z < grid.nz_layer) &&
                        (neighbor_coord_z >= 0)) {
                        int neighbor_index = grid.get_1D_index(neighbor_coord_x, neighbor_coord_y, neighbor_coord_z);
                        if (CellType(neighbor_index) == Liquid)
                            DeactivateCell = false;
                        // Capture of cell located at "NeighborD3D1ConvPosition" if this condition is satisfied
                        if ((interface.diagonal_length(index) >= interface.crit_diagonal_length(26 * index + l)) &&
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
                                float cxold = interface.octahedron_center(3 * index);
                                float cyold = interface.octahedron_center(3 * index + 1);
                                float czold = interface.octahedron_center(3 * index + 2);

                                // (xp,yp,zp) are the global coordinates of the new cell's center
                                // Note that the Y coordinate is relative to the domain origin to keep the coordinate
                                // system continuous across ranks
                                float xp = coord_x + interface.neighbor_x[l] + 0.5;
                                float yp = coord_y + grid.y_offset + interface.neighbor_y[l] + 0.5;
                                float zp = coord_z + interface.neighbor_z[l] + 0.5;

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
                                // If minimum distance to corner = 0, the octahedron corner captured the new cell
                                // center
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

                                interface.diagonal_length(neighbor_index) = NewODiagL;
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

                                interface.octahedron_center(3 * neighbor_index) = cx;
                                interface.octahedron_center(3 * neighbor_index + 1) = cy;
                                interface.octahedron_center(3 * neighbor_index + 2) = cz;

                                // Get new critical diagonal length values for the newly activated cell (at array
                                // position "neighbor_index")
                                interface.calc_crit_diagonal_length(neighbor_index, xp, yp, zp, cx, cy, cz,
                                                                    MyOrientation, GrainUnitVector);

                                if (np > 1) {
                                    // TODO: Test loading ghost nodes in a separate kernel, potentially adopting
                                    // this change if the slowdown is minor
                                    int ghost_grain_id = h;
                                    float ghost_octahedron_center_x = cx;
                                    float ghost_octahedron_center_y = cy;
                                    float ghost_octahedron_center_z = cz;
                                    float ghost_diagonal_length = NewODiagL;
                                    // Collect data for the ghost nodes, if necessary
                                    // Data loaded into the ghost nodes is for the cell that was just captured
                                    bool data_fits_in_buffer = interface.load_ghost_nodes(
                                        ghost_grain_id, ghost_octahedron_center_x, ghost_octahedron_center_y,
                                        ghost_octahedron_center_z, ghost_diagonal_length, grid.ny_local,
                                        neighbor_coord_x, neighbor_coord_y, neighbor_coord_z, grid.at_north_boundary,
                                        grid.at_south_boundary, NGrainOrientations);
                                    if (!(data_fits_in_buffer)) {
                                        // This cell's data did not fit in the buffer with current size buf_size -
                                        // mark with temporary type
                                        CellType(neighbor_index) = ActiveFailedBufferLoad;
                                    }
                                    else {
                                        // Cell activation is now finished - cell type can be changed from
                                        // TemporaryUpdate to Active
                                        CellType(neighbor_index) = Active;
                                    }
                                } // End if statement for serial/parallel code
                                else {
                                    // Only update the new cell's type once Critical Diagonal Length, Triangle
                                    // Index, and Diagonal Length values have been assigned to it Avoids the race
                                    // condition in which the new cell is activated, and another thread acts on the
                                    // new active cell before the cell's new critical diagonal length/triangle
                                    // index/diagonal length values are assigned
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
                interface.create_new_octahedron(index, coord_x, coord_y, grid.y_offset, coord_z);
                // The orientation for the new grain will depend on its Grain ID (nucleated grains have negative
                // GrainID values)
                int MyOrientation = getGrainOrientation(MyGrainID, NGrainOrientations);
                // Octahedron center is at (cx, cy, cz) - note that the Y coordinate is relative to the domain
                // origin to keep the coordinate system continuous across ranks
                float cx = coord_x + 0.5;
                float cy = coord_y + grid.y_offset + 0.5;
                float cz = coord_z + 0.5;
                // Calculate critical values at which this active cell leads to the activation of a neighboring
                // liquid cell. Octahedron center and cell center overlap for octahedra created as part of a new
                // grain
                interface.calc_crit_diagonal_length(index, cx, cy, cz, cx, cy, cz, MyOrientation, GrainUnitVector);
                if (np > 1) {
                    // TODO: Test loading ghost nodes in a separate kernel, potentially adopting this change if the
                    // slowdown is minor
                    int ghost_grain_id = MyGrainID;
                    float ghost_octahedron_center_x = cx;
                    float ghost_octahedron_center_y = cy;
                    float ghost_octahedron_center_z = cz;
                    float ghost_diagonal_length = 0.01;
                    // Collect data for the ghost nodes, if necessary
                    bool data_fits_in_buffer = interface.load_ghost_nodes(
                        ghost_grain_id, ghost_octahedron_center_x, ghost_octahedron_center_y, ghost_octahedron_center_z,
                        ghost_diagonal_length, grid.ny_local, coord_x, coord_y, coord_z, grid.at_north_boundary,
                        grid.at_south_boundary, NGrainOrientations);
                    if (!(data_fits_in_buffer)) {
                        // This cell's data did not fit in the buffer with current size buf_size - mark with
                        // temporary type
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
                bool data_fits_in_buffer =
                    interface.load_ghost_nodes(-1, -1.0, -1.0, -1.0, 0.0, grid.ny_local, coord_x, coord_y, coord_z,
                                               grid.at_north_boundary, grid.at_south_boundary, NGrainOrientations);
                if (!(data_fits_in_buffer)) {
                    // This cell's data did not fit in the buffer with current size buf_size - mark with temporary
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

// Check buffers for overflow and resize/refill as necessary
void check_buffers(const int id, const Grid &grid, CellData<device_memory_space> &cellData,
                   Interface<device_memory_space> &interface, const int n_grain_orientations) {
    // Count the number of cells' in halo regions where the data did not fit into the send buffers
    // Reduce across all ranks, as the same buf_size should be maintained across all ranks
    // If any rank overflowed its buffer size, resize all buffers to the new size plus 10% padding
    int old_buf_size = interface.buf_size;
    int buf_size_local = interface.resize_buffers();
    if (old_buf_size != buf_size_local) {
        if (id == 0)
            std::cout << "Resized number of cells stored in send/recv buffers from " << old_buf_size << " to "
                      << buf_size_local << std::endl;
        refill_buffers(grid, cellData, interface, n_grain_orientations);
    }
}

// Refill the buffers as necessary starting from the old count size, using the data from cells marked with type
// ActiveFailedBufferLoad
void refill_buffers(const Grid &grid, CellData<device_memory_space> &cellData,
                    Interface<device_memory_space> &interface, const int n_grain_orientations) {

    auto CellType = cellData.getCellTypeSubview(grid);
    auto GrainID = cellData.getGrainIDSubview(grid);
    Kokkos::parallel_for(
        "FillSendBuffersOverflow", grid.nx, KOKKOS_LAMBDA(const int &coord_x) {
            for (int coord_z = 0; coord_z < grid.nz_layer; coord_z++) {
                int index_south_buffer = grid.get_1D_index(coord_x, 1, coord_z);
                int index_north_buffer = grid.get_1D_index(coord_x, grid.ny_local - 2, coord_z);
                if (CellType(index_south_buffer) == ActiveFailedBufferLoad) {
                    int ghost_grain_id = GrainID(index_south_buffer);
                    float ghost_octahedron_center_x = interface.octahedron_center(3 * index_south_buffer);
                    float ghost_octahedron_center_y = interface.octahedron_center(3 * index_south_buffer + 1);
                    float ghost_octahedron_center_z = interface.octahedron_center(3 * index_south_buffer + 2);
                    float ghost_diagonal_length = interface.diagonal_length(index_south_buffer);
                    // Collect data for the ghost nodes, if necessary
                    // Data loaded into the ghost nodes is for the cell that was just captured
                    bool data_fits_in_buffer = interface.load_ghost_nodes(
                        ghost_grain_id, ghost_octahedron_center_x, ghost_octahedron_center_y, ghost_octahedron_center_z,
                        ghost_diagonal_length, grid.ny_local, coord_x, 1, coord_z, grid.at_north_boundary,
                        grid.at_south_boundary, n_grain_orientations);
                    CellType(index_south_buffer) = Active;
                    // If data doesn't fit in the buffer after the resize, warn that buffer data may have been lost
                    if (!(data_fits_in_buffer))
                        printf("Error: Send/recv buffer resize failed to include all necessary data, predicted "
                               "results at MPI processor boundaries may be inaccurate\n");
                }
                else if (CellType(index_south_buffer) == LiquidFailedBufferLoad) {
                    // Dummy values for first 4 arguments (Grain ID and octahedron center coordinates), 0 for
                    // diagonal length
                    bool data_fits_in_buffer = interface.load_ghost_nodes(-1, -1.0, -1.0, -1.0, 0.0, grid.ny_local,
                                                                          coord_x, 1, coord_z, grid.at_north_boundary,
                                                                          grid.at_south_boundary, n_grain_orientations);
                    CellType(index_south_buffer) = Liquid;
                    // If data doesn't fit in the buffer after the resize, warn that buffer data may have been lost
                    if (!(data_fits_in_buffer))
                        printf("Error: Send/recv buffer resize failed to include all necessary data, predicted "
                               "results at MPI processor boundaries may be inaccurate\n");
                }
                if (CellType(index_north_buffer) == ActiveFailedBufferLoad) {
                    int ghost_grain_id = GrainID(index_north_buffer);
                    float ghost_octahedron_center_x = interface.octahedron_center(3 * index_north_buffer);
                    float ghost_octahedron_center_y = interface.octahedron_center(3 * index_north_buffer + 1);
                    float ghost_octahedron_center_z = interface.octahedron_center(3 * index_north_buffer + 2);
                    float ghost_diagonal_length = interface.diagonal_length(index_north_buffer);
                    // Collect data for the ghost nodes, if necessary
                    // Data loaded into the ghost nodes is for the cell that was just captured
                    bool data_fits_in_buffer = interface.load_ghost_nodes(
                        ghost_grain_id, ghost_octahedron_center_x, ghost_octahedron_center_y, ghost_octahedron_center_z,
                        ghost_diagonal_length, grid.ny_local, coord_x, grid.ny_local - 2, coord_z,
                        grid.at_north_boundary, grid.at_south_boundary, n_grain_orientations);
                    CellType(index_north_buffer) = Active;
                    // If data doesn't fit in the buffer after the resize, warn that buffer data may have been lost
                    if (!(data_fits_in_buffer))
                        printf("Error: Send/recv buffer resize failed to include all necessary data, predicted "
                               "results at MPI processor boundaries may be inaccurate\n");
                }
                else if (CellType(index_north_buffer) == LiquidFailedBufferLoad) {
                    // Dummy values for first 4 arguments (Grain ID and octahedron center coordinates), 0 for
                    // diagonal length
                    bool data_fits_in_buffer = interface.load_ghost_nodes(
                        -1, -1.0, -1.0, -1.0, 0.0, grid.ny_local, coord_x, grid.ny_local - 2, coord_z,
                        grid.at_north_boundary, grid.at_south_boundary, n_grain_orientations);
                    CellType(index_north_buffer) = Liquid;
                    // If data doesn't fit in the buffer after the resize, warn that buffer data may have been lost
                    if (!(data_fits_in_buffer))
                        printf("Error: Send/recv buffer resize failed to include all necessary data, predicted "
                               "results at MPI processor boundaries may be inaccurate\n");
                }
            }
        });
    Kokkos::fence();
}

// 1D domain decomposition: update ghost nodes with new cell data from Nucleation and CellCapture routines
void halo_update(const int, const int, const Grid &grid, CellData<device_memory_space> &cellData,
                 Interface<device_memory_space> &interface, const int n_grain_orientations,
                 const ViewF GrainUnitVector) {

    std::vector<MPI_Request> SendRequests(2, MPI_REQUEST_NULL);
    std::vector<MPI_Request> RecvRequests(2, MPI_REQUEST_NULL);

    // Send data to each other rank (MPI_Isend)
    MPI_Isend(interface.buffer_south_send.data(), interface.buf_components * interface.buf_size, MPI_FLOAT,
              grid.neighbor_rank_south, 0, MPI_COMM_WORLD, &SendRequests[0]);
    MPI_Isend(interface.buffer_north_send.data(), interface.buf_components * interface.buf_size, MPI_FLOAT,
              grid.neighbor_rank_north, 0, MPI_COMM_WORLD, &SendRequests[1]);

    // Receive buffers for all neighbors (MPI_Irecv)
    MPI_Irecv(interface.buffer_south_recv.data(), interface.buf_components * interface.buf_size, MPI_FLOAT,
              grid.neighbor_rank_south, 0, MPI_COMM_WORLD, &RecvRequests[0]);
    MPI_Irecv(interface.buffer_north_recv.data(), interface.buf_components * interface.buf_size, MPI_FLOAT,
              grid.neighbor_rank_north, 0, MPI_COMM_WORLD, &RecvRequests[1]);

    // unpack in any order
    bool unpack_complete = false;
    auto CellType = cellData.getCellTypeSubview(grid);
    auto GrainID = cellData.getGrainIDSubview(grid);
    while (!unpack_complete) {
        // Get the next buffer to unpack from rank "unpack_index"
        int unpack_index = MPI_UNDEFINED;
        MPI_Waitany(2, RecvRequests.data(), &unpack_index, MPI_STATUS_IGNORE);
        // If there are no more buffers to unpack, leave the while loop
        if (MPI_UNDEFINED == unpack_index) {
            unpack_complete = true;
        }
        // Otherwise unpack the next buffer.
        else {
            Kokkos::parallel_for(
                "BufferUnpack", interface.buf_size, KOKKOS_LAMBDA(const int &buf_position) {
                    int coord_x, coord_y, coord_z, index, new_grain_id;
                    float new_octahedron_center_x, new_octahedron_center_y, new_octahedron_center_z,
                        new_diagonal_length;
                    bool Place = false;
                    // Which rank was the data received from? Is there valid data at this position in the buffer
                    // (i.e., not set to -1.0)?
                    if ((unpack_index == 0) && (interface.buffer_south_recv(buf_position, 0) != -1.0) &&
                        (grid.neighbor_rank_south != MPI_PROC_NULL)) {
                        // Data receieved from South
                        coord_x = static_cast<int>(interface.buffer_south_recv(buf_position, 0));
                        coord_y = 0;
                        coord_z = static_cast<int>(interface.buffer_south_recv(buf_position, 1));
                        index = grid.get_1D_index(coord_x, coord_y, coord_z);
                        // Two possibilities: buffer data with non-zero diagonal length was loaded, and a liquid
                        // cell may have to be updated to active - or zero diagonal length data was loaded, and an
                        // active cell may have to be updated to liquid
                        if (CellType(index) == Liquid) {
                            Place = true;
                            int MyGrainOrientation = static_cast<int>(interface.buffer_south_recv(buf_position, 2));
                            int MyGrainNumber = static_cast<int>(interface.buffer_south_recv(buf_position, 3));
                            new_grain_id = getGrainID(n_grain_orientations, MyGrainOrientation, MyGrainNumber);
                            new_octahedron_center_x = interface.buffer_south_recv(buf_position, 4);
                            new_octahedron_center_y = interface.buffer_south_recv(buf_position, 5);
                            new_octahedron_center_z = interface.buffer_south_recv(buf_position, 6);
                            new_diagonal_length = interface.buffer_south_recv(buf_position, 7);
                        }
                        else if ((CellType(index) == Active) && (interface.buffer_south_recv(buf_position, 7) == 0.0)) {
                            CellType(index) = Liquid;
                        }
                    }
                    else if ((unpack_index == 1) && (interface.buffer_north_recv(buf_position, 0) != -1.0) &&
                             (grid.neighbor_rank_north != MPI_PROC_NULL)) {
                        // Data received from North
                        coord_x = static_cast<int>(interface.buffer_north_recv(buf_position, 0));
                        coord_y = grid.ny_local - 1;
                        coord_z = static_cast<int>(interface.buffer_north_recv(buf_position, 1));
                        index = grid.get_1D_index(coord_x, coord_y, coord_z);
                        // Two possibilities: buffer data with non-zero diagonal length was loaded, and a liquid
                        // cell may have to be updated to active - or zero diagonal length data was loaded, and an
                        // active cell may have to be updated to liquid
                        if (CellType(index) == Liquid) {
                            Place = true;
                            int MyGrainOrientation = static_cast<int>(interface.buffer_north_recv(buf_position, 2));
                            int MyGrainNumber = static_cast<int>(interface.buffer_north_recv(buf_position, 3));
                            new_grain_id = getGrainID(n_grain_orientations, MyGrainOrientation, MyGrainNumber);
                            new_octahedron_center_x = interface.buffer_north_recv(buf_position, 4);
                            new_octahedron_center_y = interface.buffer_north_recv(buf_position, 5);
                            new_octahedron_center_z = interface.buffer_north_recv(buf_position, 6);
                            new_diagonal_length = interface.buffer_north_recv(buf_position, 7);
                        }
                        else if ((CellType(index) == Active) && (interface.buffer_north_recv(buf_position, 7) == 0.0)) {
                            CellType(index) = Liquid;
                        }
                    }
                    if (Place) {
                        // Update this ghost node cell's information with data from other rank
                        GrainID(index) = new_grain_id;
                        interface.octahedron_center(3 * index) = new_octahedron_center_x;
                        interface.octahedron_center(3 * index + 1) = new_octahedron_center_y;
                        interface.octahedron_center(3 * index + 2) = new_octahedron_center_z;
                        int MyOrientation = getGrainOrientation(GrainID(index), n_grain_orientations);
                        interface.diagonal_length(index) = static_cast<float>(new_diagonal_length);
                        // Cell center - note that the Y coordinate is relative to the domain origin to keep the
                        // coordinate system continuous across ranks
                        double xp = coord_x + 0.5;
                        double yp = coord_y + grid.y_offset + 0.5;
                        double zp = coord_z + 0.5;
                        // Calculate critical values at which this active cell leads to the activation of a
                        // neighboring liquid cell
                        interface.calc_crit_diagonal_length(index, xp, yp, zp, new_octahedron_center_x,
                                                            new_octahedron_center_y, new_octahedron_center_z,
                                                            MyOrientation, GrainUnitVector);
                        CellType(index) = Active;
                    }
                });
        }
    }

    // Reset send buffer data to -1 (used as placeholder) and reset the number of cells stored in the buffers to 0
    interface.reset_buffers();
    // Wait on send requests
    MPI_Waitall(2, SendRequests.data(), MPI_STATUSES_IGNORE);
    Kokkos::fence();
}
//*****************************************************************************/
// Jump to the next time step with work to be done, if nothing left to do in the near future
// The cells of interest are active cells, and the view checked for future work is
// MeltTimeStep Print intermediate output during this jump if PrintIdleMovieFrames = true
void JumpTimeStep(int &cycle, unsigned long int RemainingCellsOfInterest, unsigned long int LocalTempSolidCells,
                  Temperature<device_memory_space> &temperature, const Grid &grid,
                  CellData<device_memory_space> &cellData, const int id, const int layernumber, const int np,
                  const ViewF GrainUnitVector, Print print, const int NGrainOrientations) {

    auto CellType = cellData.getCellTypeSubview(grid);
    MPI_Bcast(&RemainingCellsOfInterest, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    if (RemainingCellsOfInterest == 0) {
        // If this rank still has cells that will later undergo transformation (LocalIncompleteCells > 0), check when
        // the next solid cells go above the liquidus (remelting) Otherwise, assign the largest possible time step as
        // the next time work needs to be done on the rank
        unsigned long int NextMeltTimeStep;
        if (LocalTempSolidCells > 0) {
            Kokkos::parallel_reduce(
                "CheckNextTSForWork", grid.domain_size,
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
void IntermediateOutputAndCheck(const int id, const int np, int &cycle, const Grid &grid,
                                int SuccessfulNucEvents_ThisRank, int &XSwitch, CellData<device_memory_space> &cellData,
                                Temperature<device_memory_space> &temperature, std::string SimulationType,
                                const int layernumber, const int NGrainOrientations, const ViewF GrainUnitVector,
                                Print print) {

    auto CellType = cellData.getCellTypeSubview(grid);
    auto GrainID = cellData.getGrainIDSubview(grid);
    auto LayerID = cellData.getLayerIDSubview(grid);
    unsigned long int LocalSuperheatedCells;
    unsigned long int LocalUndercooledCells;
    unsigned long int LocalActiveCells;
    unsigned long int LocalTempSolidCells;
    unsigned long int LocalFinishedSolidCells;
    Kokkos::parallel_reduce(
        "IntermediateOutput", grid.domain_size,
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
void IntermediateOutputAndCheck(const int id, int cycle, const Grid &grid, int &XSwitch, ViewI CellType_AllLayers) {

    unsigned long int LocalLiquidCells, LocalActiveCells, LocalSolidCells;
    ViewB2D EdgesReached(Kokkos::ViewAllocateWithoutInitializing("EdgesReached"), 3, 2); // init to false
    Kokkos::deep_copy(EdgesReached, false);

    Kokkos::parallel_reduce(
        "IntermediateOutput", grid.domain_size,
        KOKKOS_LAMBDA(const int &index, unsigned long int &sum_liquid, unsigned long int &sum_active,
                      unsigned long int &sum_solid) {
            if (CellType_AllLayers(index) == Liquid)
                sum_liquid += 1;
            else if (CellType_AllLayers(index) == Active) {
                sum_active += 1;
                // Did this cell reach a domain edge?
                int coord_x = grid.get_coord_X(index);
                int coord_y = grid.get_coord_Y(index);
                int coord_z = grid.get_coord_Z(index);
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
