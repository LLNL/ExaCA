// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_UPDATE_HPP
#define EXACA_UPDATE_HPP

#include "CAcelldata.hpp"
#include "CAgrid.hpp"
#include "CAinterface.hpp"
#include "CAinterfacialresponse.hpp"
#include "CAorientation.hpp"
#include "CAprint.hpp"
#include "CAtemperature.hpp"

#include <Kokkos_Core.hpp>

#include <string>

// For the case where all cells solidify once, determine which cells are associated with the "steering vector" of
// cells that are either active, or becoming active this time step
template <typename MemorySpace>
void fillSteeringVector_NoRemelt(const int cycle, const Grid &grid, CellData<MemorySpace> &celldata,
                                 Temperature<MemorySpace> &temperature, Interface<MemorySpace> &interface) {

    // Cells associated with this layer that are not solid type but have passed the liquidus (crit time step) have
    // their undercooling values updated Cells that meet the aforementioned criteria and are active type should be
    // added to the steering vector
    Kokkos::parallel_for(
        "FillSV", grid.domain_size, KOKKOS_LAMBDA(const int &index) {
            int cell_type = celldata.cell_type(index);
            bool is_not_solid = (cell_type != Solid);
            int crit_time_step = temperature.layer_time_temp_history(index, 0, 1);
            bool past_crit_time = (cycle > crit_time_step);
            bool cell_active = ((cell_type == Active) || (cell_type == FutureActive));
            if (is_not_solid && past_crit_time) {
                temperature.updateUndercooling(index);
                if (cell_active) {
                    interface.steering_vector(Kokkos::atomic_fetch_add(&interface.num_steer(0), 1)) = index;
                }
            }
        });
    Kokkos::deep_copy(interface.num_steer_host, interface.num_steer);
}

// For the case where cells may melt and solidify multiple times, determine which cells are associated with the
// "steering vector" of cells that are either active, or becoming active this time step - version with remelting
template <typename MemorySpace>
void fillSteeringVector_Remelt(const int cycle, const Grid &grid, CellData<MemorySpace> &celldata,
                               Temperature<MemorySpace> &temperature, Interface<MemorySpace> &interface) {

    auto grain_id = celldata.getGrainIDSubview(grid);
    Kokkos::parallel_for(
        "FillSV_RM", grid.domain_size, KOKKOS_LAMBDA(const int &index) {
            int celltype = celldata.cell_type(index);
            // Only iterate over cells that are not Solid type
            if (celltype != Solid) {
                int melt_time_step = temperature.getMeltTimeStep(cycle, index);
                int crit_time_step = temperature.getCritTimeStep(index);
                bool at_melt_time = (cycle == melt_time_step);
                bool at_crit_time = (cycle == crit_time_step);
                bool past_crit_time = (cycle > crit_time_step);
                if (at_melt_time) {
                    // Cell melts, undercooling is reset to 0 from the previous value, if any
                    celldata.cell_type(index) = Liquid;
                    temperature.resetUndercooling(index);
                    if (celltype != TempSolid) {
                        // This cell either hasn't started or hasn't finished the previous solidification event, but
                        // has undergone melting - increment the solidification counter to move on to the next
                        // melt-solidification event
                        temperature.updateSolidificationCounter(index);
                    }
                    // Any adjacent active cells should also be remelted, as these cells are more likely heating up
                    // than cooling down These are converted to the temporary FutureLiquid state, to be later
                    // iterated over and loaded into the steering vector as necessary Get the x, y, z coordinates of
                    // the cell on this MPI rank
                    int coord_x = grid.getCoordX(index);
                    int coord_y = grid.getCoordY(index);
                    int coord_z = grid.getCoordZ(index);
                    for (int l = 0; l < 26; l++) {
                        // "l" correpsponds to the specific neighboring cell
                        // Local coordinates of adjacent cell center
                        int neighbor_coord_x = coord_x + interface.neighbor_x[l];
                        int neighbor_coord_y = coord_y + interface.neighbor_y[l];
                        int neighbor_coord_z = coord_z + interface.neighbor_z[l];
                        if ((neighbor_coord_x >= 0) && (neighbor_coord_x < grid.nx) && (neighbor_coord_y >= 0) &&
                            (neighbor_coord_y < grid.ny_local) && (neighbor_coord_z < grid.nz_layer) &&
                            (neighbor_coord_z >= 0)) {
                            int neighbor_index = grid.get1DIndex(neighbor_coord_x, neighbor_coord_y, neighbor_coord_z);
                            if (celldata.cell_type(neighbor_index) == Active) {
                                // Mark adjacent active cells to this as cells that should be converted into liquid,
                                // as they are more likely heating than cooling
                                celldata.cell_type(neighbor_index) = FutureLiquid;
                                interface.steering_vector(Kokkos::atomic_fetch_add(&interface.num_steer(0), 1)) =
                                    neighbor_index;
                            }
                        }
                    }
                }
                else if ((celltype != TempSolid) && (past_crit_time)) {
                    // Update cell undercooling
                    temperature.updateUndercooling(index);
                    if (celltype == Active) {
                        // Add active cells below liquidus to steering vector
                        interface.steering_vector(Kokkos::atomic_fetch_add(&interface.num_steer(0), 1)) = index;
                    }
                }
                else if ((at_crit_time) && (celltype == Liquid) && (grain_id(index) != 0)) {
                    // Get the x, y, z coordinates of the cell on this MPI rank
                    int coord_x = grid.getCoordX(index);
                    int coord_y = grid.getCoordY(index);
                    int coord_z = grid.getCoordZ(index);
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
                            int neighbor_index = grid.get1DIndex(neighbor_coord_x, neighbor_coord_y, neighbor_coord_z);
                            if ((celldata.cell_type(neighbor_index) == TempSolid) ||
                                (celldata.cell_type(neighbor_index) == Solid) || (coord_z == 0)) {
                                // Cell activation to be performed as part of steering vector
                                l = 26;
                                interface.steering_vector(Kokkos::atomic_fetch_add(&interface.num_steer(0), 1)) = index;
                                celldata.cell_type(index) =
                                    FutureActive; // this cell cannot be captured - is being activated
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
template <typename MemorySpace>
void cellCapture(const int, const int np, const Grid &grid, const InterfacialResponseFunction &irf,
                 CellData<MemorySpace> &celldata, Temperature<MemorySpace> &temperature,
                 Interface<MemorySpace> &interface, Orientation<MemorySpace> &orientation) {

    // Get grain_id subview for this layer
    auto grain_id = celldata.getGrainIDSubview(grid);
    // Loop over list of active and soon-to-be active cells, potentially performing cell capture events and updating
    // cell types
    Kokkos::parallel_for(
        "CellCapture", interface.num_steer_host(0), KOKKOS_LAMBDA(const int &num) {
            interface.num_steer(0) = 0;
            int index = interface.steering_vector(num);
            // Get the x, y, z coordinates of the cell on this MPI rank
            int coord_x = grid.getCoordX(index);
            int coord_y = grid.getCoordY(index);
            int coord_z = grid.getCoordZ(index);

            // Cells of interest for the CA - active cells and future active/liquid cells
            if (celldata.cell_type(index) == Active) {
                // Get undercooling of active cell
                float loc_u = temperature.undercooling_current(index);
                // Update diagonal length of octahedron based on local undercooling and interfacial response function
                interface.diagonal_length(index) += irf.compute(loc_u);
                // Cycle through all neigboring cells on this processor to see if they have been captured
                // Cells in ghost nodes cannot capture cells on other processors
                bool deactivate_cell =
                    true; // switch that becomes false if the cell has at least 1 liquid type neighbor
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
                        int neighbor_index = grid.get1DIndex(neighbor_coord_x, neighbor_coord_y, neighbor_coord_z);
                        if (celldata.cell_type(neighbor_index) == Liquid)
                            deactivate_cell = false;
                        // Capture of cell located at "neighbor_index" if this condition is satisfied
                        if ((interface.diagonal_length(index) >= interface.crit_diagonal_length(26 * index + l)) &&
                            (celldata.cell_type(neighbor_index) == Liquid)) {
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
                            int old_cell_type_value = Kokkos::atomic_compare_exchange(
                                &celldata.cell_type(neighbor_index), old_val, update_val);
                            // Only proceed if cell_type was previously liquid (this current thread changed the value to
                            // TemporaryUpdate)
                            if (old_cell_type_value == Liquid) {
                                int h = grain_id(index);
                                int my_orientation = getGrainOrientation(h, orientation.n_grain_orientations);

                                // The new cell is captured by this cell's growing octahedron (Grain "h")
                                grain_id(neighbor_index) = h;

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
                                float mag0 = Kokkos::sqrtf(x0 * x0 + y0 * y0 + z0 * z0);

                                // Calculate unit vectors for the octahedron that intersect the new cell center
                                // TODO: Potential inline function in orientation struct that returns Diag1, Diag2, or
                                // Diag3
                                float angle_1 = (orientation.grain_unit_vector(9 * my_orientation) * x0 +
                                                 orientation.grain_unit_vector(9 * my_orientation + 1) * y0 +
                                                 orientation.grain_unit_vector(9 * my_orientation + 2) * z0) /
                                                mag0;
                                float angle_2 = (orientation.grain_unit_vector(9 * my_orientation + 3) * x0 +
                                                 orientation.grain_unit_vector(9 * my_orientation + 4) * y0 +
                                                 orientation.grain_unit_vector(9 * my_orientation + 5) * z0) /
                                                mag0;
                                float angle_3 = (orientation.grain_unit_vector(9 * my_orientation + 6) * x0 +
                                                 orientation.grain_unit_vector(9 * my_orientation + 7) * y0 +
                                                 orientation.grain_unit_vector(9 * my_orientation + 8) * z0) /
                                                mag0;
                                float diag_1x =
                                    orientation.grain_unit_vector(9 * my_orientation) * (2 * (angle_1 < 0) - 1);
                                float diag_1y =
                                    orientation.grain_unit_vector(9 * my_orientation + 1) * (2 * (angle_1 < 0) - 1);
                                float diag_1z =
                                    orientation.grain_unit_vector(9 * my_orientation + 2) * (2 * (angle_1 < 0) - 1);

                                float diag_2x =
                                    orientation.grain_unit_vector(9 * my_orientation + 3) * (2 * (angle_2 < 0) - 1);
                                float diag_2y =
                                    orientation.grain_unit_vector(9 * my_orientation + 4) * (2 * (angle_2 < 0) - 1);
                                float diag_2z =
                                    orientation.grain_unit_vector(9 * my_orientation + 5) * (2 * (angle_2 < 0) - 1);

                                float diag_3x =
                                    orientation.grain_unit_vector(9 * my_orientation + 6) * (2 * (angle_3 < 0) - 1);
                                float diag_3y =
                                    orientation.grain_unit_vector(9 * my_orientation + 7) * (2 * (angle_3 < 0) - 1);
                                float diag_3z =
                                    orientation.grain_unit_vector(9 * my_orientation + 8) * (2 * (angle_3 < 0) - 1);

                                float u1[3], u2[3];
                                u1[0] = diag_2x - diag_1x;
                                u1[1] = diag_2y - diag_1y;
                                u1[2] = diag_2z - diag_1z;
                                u2[0] = diag_3x - diag_1x;
                                u2[1] = diag_3y - diag_1y;
                                u2[2] = diag_3z - diag_1z;
                                float uu[3];
                                uu[0] = u1[1] * u2[2] - u1[2] * u2[1];
                                uu[1] = u1[2] * u2[0] - u1[0] * u2[2];
                                uu[2] = u1[0] * u2[1] - u1[1] * u2[0];
                                float n_dem = Kokkos::hypot(uu[0], uu[1], uu[2]);
                                float norm[3];
                                norm[0] = uu[0] / n_dem;
                                norm[1] = uu[1] / n_dem;
                                norm[2] = uu[2] / n_dem;
                                // normal to capturing plane
                                float norm_plane[3], triangle_x[3], triangle_y[3], triangle_z[3], para_t;
                                norm_plane[0] = norm[0];
                                norm_plane[1] = norm[1];
                                norm_plane[2] = norm[2];
                                para_t = (norm_plane[0] * x0 + norm_plane[1] * y0 + norm_plane[2] * z0) /
                                         (norm_plane[0] * diag_1x + norm_plane[1] * diag_1y + norm_plane[2] * diag_1z);

                                triangle_x[0] = cxold + para_t * diag_1x;
                                triangle_y[0] = cyold + para_t * diag_1y;
                                triangle_z[0] = czold + para_t * diag_1z;

                                triangle_x[1] = cxold + para_t * diag_2x;
                                triangle_y[1] = cyold + para_t * diag_2y;
                                triangle_z[1] = czold + para_t * diag_2z;

                                triangle_x[2] = cxold + para_t * diag_3x;
                                triangle_y[2] = cyold + para_t * diag_3y;
                                triangle_z[2] = czold + para_t * diag_3z;

                                // Determine which of the 3 corners of the capturing face is closest to the captured
                                // cell center
                                float dist_to_corner[3];
                                dist_to_corner[0] =
                                    Kokkos::hypot(triangle_x[0] - xp, triangle_y[0] - yp, triangle_z[0] - zp);
                                dist_to_corner[1] =
                                    Kokkos::hypot(triangle_x[1] - xp, triangle_y[1] - yp, triangle_z[1] - zp);
                                dist_to_corner[2] =
                                    Kokkos::hypot(triangle_x[2] - xp, triangle_y[2] - yp, triangle_z[2] - zp);

                                int x, y, z;
                                x = (dist_to_corner[0] < dist_to_corner[1]);
                                y = (dist_to_corner[1] < dist_to_corner[2]);
                                z = (dist_to_corner[2] < dist_to_corner[0]);

                                int idx = 2 * (z - y) * z + (y - x) * y;
                                float mindist_to_corner = dist_to_corner[idx];
                                float xc = triangle_x[idx];
                                float yc = triangle_y[idx];
                                float zc = triangle_z[idx];

                                float x1 = triangle_x[(idx + 1) % 3];
                                float y1 = triangle_y[(idx + 1) % 3];
                                float z1 = triangle_z[(idx + 1) % 3];
                                float x2 = triangle_x[(idx + 2) % 3];
                                float y2 = triangle_y[(idx + 2) % 3];
                                float z2 = triangle_z[(idx + 2) % 3];

                                float d1 = Kokkos::hypot(xp - x2, yp - y2, zp - z2);
                                float d2 = Kokkos::hypot(xc - x2, yc - y2, zc - z2);
                                float d3 = Kokkos::hypot(xp - x1, yp - y1, zp - z1);
                                float d4 = Kokkos::hypot(xc - x1, yc - y1, zc - z1);

                                float i_1 = 0;
                                float i_2 = d2;
                                float j_1 = 0;
                                float j_2 = d4;
                                // If minimum distance to corner = 0, the octahedron corner captured the new cell
                                // center
                                if (mindist_to_corner != 0) {
                                    i_1 = d1 * ((xp - x2) * (xc - x2) + (yp - y2) * (yc - y2) + (zp - z2) * (zc - z2)) /
                                          (d1 * d2);
                                    i_2 = d2 - i_1;
                                    j_1 = d3 * ((xp - x1) * (xc - x1) + (yp - y1) * (yc - y1) + (zp - z1) * (zc - z1)) /
                                          (d3 * d4);
                                    j_2 = d4 - j_1;
                                }
                                float l_12 = 0.5 * (Kokkos::fmin(i_1, Kokkos::sqrt(3.0f)) +
                                                    Kokkos::fmin(i_2, Kokkos::sqrt(3.0f)));
                                float l_13 = 0.5 * (Kokkos::fmin(j_1, Kokkos::sqrt(3.0f)) +
                                                    Kokkos::fmin(j_2, Kokkos::sqrt(3.0f)));
                                // half diagonal length of new octahedron
                                float new_o_diag_l = Kokkos::sqrt(2.0f) * Kokkos::fmax(l_12, l_13);

                                interface.diagonal_length(neighbor_index) = new_o_diag_l;
                                // Calculate coordinates of new decentered octahedron center
                                float capt_diag[3];
                                capt_diag[0] = xc - cxold;
                                capt_diag[1] = yc - cyold;
                                capt_diag[2] = zc - czold;

                                float capt_diag_magnitude =
                                    Kokkos::sqrt(capt_diag[0] * capt_diag[0] + capt_diag[1] * capt_diag[1] +
                                                 capt_diag[2] * capt_diag[2]);
                                float capt_diag_uv[3];
                                capt_diag_uv[0] = capt_diag[0] / capt_diag_magnitude;
                                capt_diag_uv[1] = capt_diag[1] / capt_diag_magnitude;
                                capt_diag_uv[2] = capt_diag[2] / capt_diag_magnitude;
                                // (cx, cy, cz) are the coordiantes of the new active cell's decentered octahedron
                                float cx = xc - new_o_diag_l * capt_diag_uv[0];
                                float cy = yc - new_o_diag_l * capt_diag_uv[1];
                                float cz = zc - new_o_diag_l * capt_diag_uv[2];

                                interface.octahedron_center(3 * neighbor_index) = cx;
                                interface.octahedron_center(3 * neighbor_index + 1) = cy;
                                interface.octahedron_center(3 * neighbor_index + 2) = cz;

                                // Get new critical diagonal length values for the newly activated cell (at array
                                // position "neighbor_index")
                                interface.calcCritDiagonalLength(neighbor_index, xp, yp, zp, cx, cy, cz, my_orientation,
                                                                 orientation.grain_unit_vector);

                                if (np > 1) {
                                    // TODO: Test loading ghost nodes in a separate kernel, potentially adopting
                                    // this change if the slowdown is minor
                                    int ghost_grain_id = h;
                                    float ghost_octahedron_center_x = cx;
                                    float ghost_octahedron_center_y = cy;
                                    float ghost_octahedron_center_z = cz;
                                    float ghost_diagonal_length = new_o_diag_l;
                                    // Collect data for the ghost nodes, if necessary
                                    // Data loaded into the ghost nodes is for the cell that was just captured
                                    bool data_fits_in_buffer = interface.loadGhostNodes(
                                        ghost_grain_id, ghost_octahedron_center_x, ghost_octahedron_center_y,
                                        ghost_octahedron_center_z, ghost_diagonal_length, grid.ny_local,
                                        neighbor_coord_x, neighbor_coord_y, neighbor_coord_z, grid.at_north_boundary,
                                        grid.at_south_boundary, orientation.n_grain_orientations);
                                    if (!(data_fits_in_buffer)) {
                                        // This cell's data did not fit in the buffer with current size buf_size -
                                        // mark with temporary type
                                        celldata.cell_type(neighbor_index) = ActiveFailedBufferLoad;
                                    }
                                    else {
                                        // Cell activation is now finished - cell type can be changed from
                                        // TemporaryUpdate to Active
                                        celldata.cell_type(neighbor_index) = Active;
                                    }
                                } // End if statement for serial/parallel code
                                else {
                                    // Only update the new cell's type once Critical Diagonal Length, Triangle
                                    // Index, and Diagonal Length values have been assigned to it Avoids the race
                                    // condition in which the new cell is activated, and another thread acts on the
                                    // new active cell before the cell's new critical diagonal length/triangle
                                    // index/diagonal length values are assigned
                                    celldata.cell_type(neighbor_index) = Active;
                                }
                            } // End if statement within locked capture loop
                        }     // End if statement for outer capture loop
                    }         // End if statement over neighbors on the active grid
                }             // End loop over all neighbors of this active cell
                if (deactivate_cell) {
                    // This active cell has no more neighboring cells to be captured
                    // Update the counter for the number of times this cell went from liquid to active to solid
                    bool solidification_complete_y_n = temperature.updateCheckSolidificationCounter(index);
                    // Did the cell solidify for the last time in the layer?
                    // If so, this cell is solid - ignore until next layer (if needed)
                    // If not, this cell is tempsolid, will become liquid again
                    if (solidification_complete_y_n)
                        celldata.cell_type(index) = Solid;
                    else
                        celldata.cell_type(index) = TempSolid;
                }
            }
            else if (celldata.cell_type(index) == FutureActive) {
                // Successful nucleation event - this cell is becoming a new active cell
                celldata.cell_type(index) = TemporaryUpdate; // avoid operating on the new active cell before its
                                                             // associated octahedron data is initialized
                int my_grain_id = grain_id(index);           // grain_id was assigned as part of Nucleation

                // Initialize new octahedron
                interface.createNewOctahedron(index, coord_x, coord_y, grid.y_offset, coord_z);
                // The orientation for the new grain will depend on its Grain ID (nucleated grains have negative
                // grain_id values)
                int my_orientation = getGrainOrientation(my_grain_id, orientation.n_grain_orientations);
                // Octahedron center is at (cx, cy, cz) - note that the Y coordinate is relative to the domain
                // origin to keep the coordinate system continuous across ranks
                float cx = coord_x + 0.5;
                float cy = coord_y + grid.y_offset + 0.5;
                float cz = coord_z + 0.5;
                // Calculate critical values at which this active cell leads to the activation of a neighboring
                // liquid cell. Octahedron center and cell center overlap for octahedra created as part of a new
                // grain
                interface.calcCritDiagonalLength(index, cx, cy, cz, cx, cy, cz, my_orientation,
                                                 orientation.grain_unit_vector);
                if (np > 1) {
                    // TODO: Test loading ghost nodes in a separate kernel, potentially adopting this change if the
                    // slowdown is minor
                    int ghost_grain_id = my_grain_id;
                    float ghost_octahedron_center_x = cx;
                    float ghost_octahedron_center_y = cy;
                    float ghost_octahedron_center_z = cz;
                    float ghost_diagonal_length = 0.01;
                    // Collect data for the ghost nodes, if necessary
                    bool data_fits_in_buffer = interface.loadGhostNodes(
                        ghost_grain_id, ghost_octahedron_center_x, ghost_octahedron_center_y, ghost_octahedron_center_z,
                        ghost_diagonal_length, grid.ny_local, coord_x, coord_y, coord_z, grid.at_north_boundary,
                        grid.at_south_boundary, orientation.n_grain_orientations);
                    if (!(data_fits_in_buffer)) {
                        // This cell's data did not fit in the buffer with current size buf_size - mark with
                        // temporary type
                        celldata.cell_type(index) = ActiveFailedBufferLoad;
                    }
                    else {
                        // Cell activation is now finished - cell type can be changed from TemporaryUpdate to Active
                        celldata.cell_type(index) = Active;
                    }
                } // End if statement for serial/parallel code
                else {
                    // Cell activation is now finished - cell type can be changed from TemporaryUpdate to Active
                    celldata.cell_type(index) = Active;
                } // End if statement for serial/parallel code
            }
            else if (celldata.cell_type(index) == FutureLiquid) {
                // This type was assigned to a cell that was recently transformed from active to liquid, due to its
                // bordering of a cell above the liquidus. This information may need to be sent to other MPI ranks
                // Dummy values for first 4 arguments (Grain ID and octahedron center coordinates), 0 for diagonal
                // length
                bool data_fits_in_buffer = interface.loadGhostNodes(
                    -1, -1.0, -1.0, -1.0, 0.0, grid.ny_local, coord_x, coord_y, coord_z, grid.at_north_boundary,
                    grid.at_south_boundary, orientation.n_grain_orientations);
                if (!(data_fits_in_buffer)) {
                    // This cell's data did not fit in the buffer with current size buf_size - mark with temporary
                    // type
                    celldata.cell_type(index) = LiquidFailedBufferLoad;
                }
                else {
                    // Cell activation is now finished - cell type can be changed from FutureLiquid to Active
                    celldata.cell_type(index) = Liquid;
                }
            }
        });
    Kokkos::fence();
}

// Check buffers for overflow and resize/refill as necessary
template <typename MemorySpace>
void checkBuffers(const int id, const int cycle, const Grid &grid, CellData<MemorySpace> &celldata,
                  Interface<MemorySpace> &interface, const int n_grain_orientations) {
    // Count the number of cells' in halo regions where the data did not fit into the send buffers
    // Reduce across all ranks, as the same buf_size should be maintained across all ranks
    // If any rank overflowed its buffer size, resize all buffers to the new size plus 10% padding
    bool resize_performed = interface.resizeBuffers(id, cycle);
    if (resize_performed)
        refillBuffers(grid, celldata, interface, n_grain_orientations);
}

// Refill the buffers as necessary starting from the old count size, using the data from cells marked with type
// ActiveFailedBufferLoad
template <typename MemorySpace>
void refillBuffers(const Grid &grid, CellData<MemorySpace> &celldata, Interface<MemorySpace> &interface,
                   const int n_grain_orientations) {

    auto grain_id = celldata.getGrainIDSubview(grid);
    Kokkos::parallel_for(
        "FillSendBuffersOverflow", grid.nx, KOKKOS_LAMBDA(const int &coord_x) {
            for (int coord_z = 0; coord_z < grid.nz_layer; coord_z++) {
                int index_south_buffer = grid.get1DIndex(coord_x, 1, coord_z);
                int index_north_buffer = grid.get1DIndex(coord_x, grid.ny_local - 2, coord_z);
                if (celldata.cell_type(index_south_buffer) == ActiveFailedBufferLoad) {
                    int ghost_grain_id = grain_id(index_south_buffer);
                    float ghost_octahedron_center_x = interface.octahedron_center(3 * index_south_buffer);
                    float ghost_octahedron_center_y = interface.octahedron_center(3 * index_south_buffer + 1);
                    float ghost_octahedron_center_z = interface.octahedron_center(3 * index_south_buffer + 2);
                    float ghost_diagonal_length = interface.diagonal_length(index_south_buffer);
                    // Collect data for the ghost nodes, if necessary
                    // Data loaded into the ghost nodes is for the cell that was just captured
                    bool data_fits_in_buffer = interface.loadGhostNodes(
                        ghost_grain_id, ghost_octahedron_center_x, ghost_octahedron_center_y, ghost_octahedron_center_z,
                        ghost_diagonal_length, grid.ny_local, coord_x, 1, coord_z, grid.at_north_boundary,
                        grid.at_south_boundary, n_grain_orientations);
                    celldata.cell_type(index_south_buffer) = Active;
                    // If data doesn't fit in the buffer after the resize, warn that buffer data may have been lost
                    if (!(data_fits_in_buffer))
                        printf("Error: Send/recv buffer resize failed to include all necessary data, predicted "
                               "results at MPI processor boundaries may be inaccurate\n");
                }
                else if (celldata.cell_type(index_south_buffer) == LiquidFailedBufferLoad) {
                    // Dummy values for first 4 arguments (Grain ID and octahedron center coordinates), 0 for
                    // diagonal length
                    bool data_fits_in_buffer =
                        interface.loadGhostNodes(-1, -1.0, -1.0, -1.0, 0.0, grid.ny_local, coord_x, 1, coord_z,
                                                 grid.at_north_boundary, grid.at_south_boundary, n_grain_orientations);
                    celldata.cell_type(index_south_buffer) = Liquid;
                    // If data doesn't fit in the buffer after the resize, warn that buffer data may have been lost
                    if (!(data_fits_in_buffer))
                        printf("Error: Send/recv buffer resize failed to include all necessary data, predicted "
                               "results at MPI processor boundaries may be inaccurate\n");
                }
                if (celldata.cell_type(index_north_buffer) == ActiveFailedBufferLoad) {
                    int ghost_grain_id = grain_id(index_north_buffer);
                    float ghost_octahedron_center_x = interface.octahedron_center(3 * index_north_buffer);
                    float ghost_octahedron_center_y = interface.octahedron_center(3 * index_north_buffer + 1);
                    float ghost_octahedron_center_z = interface.octahedron_center(3 * index_north_buffer + 2);
                    float ghost_diagonal_length = interface.diagonal_length(index_north_buffer);
                    // Collect data for the ghost nodes, if necessary
                    // Data loaded into the ghost nodes is for the cell that was just captured
                    bool data_fits_in_buffer = interface.loadGhostNodes(
                        ghost_grain_id, ghost_octahedron_center_x, ghost_octahedron_center_y, ghost_octahedron_center_z,
                        ghost_diagonal_length, grid.ny_local, coord_x, grid.ny_local - 2, coord_z,
                        grid.at_north_boundary, grid.at_south_boundary, n_grain_orientations);
                    celldata.cell_type(index_north_buffer) = Active;
                    // If data doesn't fit in the buffer after the resize, warn that buffer data may have been lost
                    if (!(data_fits_in_buffer))
                        printf("Error: Send/recv buffer resize failed to include all necessary data, predicted "
                               "results at MPI processor boundaries may be inaccurate\n");
                }
                else if (celldata.cell_type(index_north_buffer) == LiquidFailedBufferLoad) {
                    // Dummy values for first 4 arguments (Grain ID and octahedron center coordinates), 0 for
                    // diagonal length
                    bool data_fits_in_buffer = interface.loadGhostNodes(
                        -1, -1.0, -1.0, -1.0, 0.0, grid.ny_local, coord_x, grid.ny_local - 2, coord_z,
                        grid.at_north_boundary, grid.at_south_boundary, n_grain_orientations);
                    celldata.cell_type(index_north_buffer) = Liquid;
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
template <typename MemorySpace>
void haloUpdate(const int, const int, const Grid &grid, CellData<MemorySpace> &celldata,
                Interface<MemorySpace> &interface, Orientation<MemorySpace> &orientation) {

    std::vector<MPI_Request> send_requests(2, MPI_REQUEST_NULL);
    std::vector<MPI_Request> recv_requests(2, MPI_REQUEST_NULL);

    // Send data to each other rank (MPI_Isend)
    MPI_Isend(interface.buffer_south_send.data(), interface.buf_components * interface.buf_size, MPI_FLOAT,
              grid.neighbor_rank_south, 0, MPI_COMM_WORLD, &send_requests[0]);
    MPI_Isend(interface.buffer_north_send.data(), interface.buf_components * interface.buf_size, MPI_FLOAT,
              grid.neighbor_rank_north, 0, MPI_COMM_WORLD, &send_requests[1]);

    // Receive buffers for all neighbors (MPI_Irecv)
    MPI_Irecv(interface.buffer_south_recv.data(), interface.buf_components * interface.buf_size, MPI_FLOAT,
              grid.neighbor_rank_south, 0, MPI_COMM_WORLD, &recv_requests[0]);
    MPI_Irecv(interface.buffer_north_recv.data(), interface.buf_components * interface.buf_size, MPI_FLOAT,
              grid.neighbor_rank_north, 0, MPI_COMM_WORLD, &recv_requests[1]);

    // unpack in any order
    bool unpack_complete = false;
    auto grain_id = celldata.getGrainIDSubview(grid);
    while (!unpack_complete) {
        // Get the next buffer to unpack from rank "unpack_index"
        int unpack_index = MPI_UNDEFINED;
        MPI_Waitany(2, recv_requests.data(), &unpack_index, MPI_STATUS_IGNORE);
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
                    bool place = false;
                    // Which rank was the data received from? Is there valid data at this position in the buffer
                    // (i.e., not set to -1.0)?
                    if ((unpack_index == 0) && (interface.buffer_south_recv(buf_position, 0) != -1.0) &&
                        (grid.neighbor_rank_south != MPI_PROC_NULL)) {
                        // Data receieved from South
                        coord_x = static_cast<int>(interface.buffer_south_recv(buf_position, 0));
                        coord_y = 0;
                        coord_z = static_cast<int>(interface.buffer_south_recv(buf_position, 1));
                        index = grid.get1DIndex(coord_x, coord_y, coord_z);
                        // Two possibilities: buffer data with non-zero diagonal length was loaded, and a liquid
                        // cell may have to be updated to active - or zero diagonal length data was loaded, and an
                        // active cell may have to be updated to liquid
                        if (celldata.cell_type(index) == Liquid) {
                            place = true;
                            int my_grain_orientation = static_cast<int>(interface.buffer_south_recv(buf_position, 2));
                            int my_grain_number = static_cast<int>(interface.buffer_south_recv(buf_position, 3));
                            new_grain_id =
                                getGrainID(my_grain_orientation, my_grain_number, orientation.n_grain_orientations);
                            new_octahedron_center_x = interface.buffer_south_recv(buf_position, 4);
                            new_octahedron_center_y = interface.buffer_south_recv(buf_position, 5);
                            new_octahedron_center_z = interface.buffer_south_recv(buf_position, 6);
                            new_diagonal_length = interface.buffer_south_recv(buf_position, 7);
                        }
                        else if ((celldata.cell_type(index) == Active) &&
                                 (interface.buffer_south_recv(buf_position, 7) == 0.0)) {
                            celldata.cell_type(index) = Liquid;
                        }
                    }
                    else if ((unpack_index == 1) && (interface.buffer_north_recv(buf_position, 0) != -1.0) &&
                             (grid.neighbor_rank_north != MPI_PROC_NULL)) {
                        // Data received from North
                        coord_x = static_cast<int>(interface.buffer_north_recv(buf_position, 0));
                        coord_y = grid.ny_local - 1;
                        coord_z = static_cast<int>(interface.buffer_north_recv(buf_position, 1));
                        index = grid.get1DIndex(coord_x, coord_y, coord_z);
                        // Two possibilities: buffer data with non-zero diagonal length was loaded, and a liquid
                        // cell may have to be updated to active - or zero diagonal length data was loaded, and an
                        // active cell may have to be updated to liquid
                        if (celldata.cell_type(index) == Liquid) {
                            place = true;
                            int my_grain_orientation = static_cast<int>(interface.buffer_north_recv(buf_position, 2));
                            int my_grain_number = static_cast<int>(interface.buffer_north_recv(buf_position, 3));
                            new_grain_id =
                                getGrainID(my_grain_orientation, my_grain_number, orientation.n_grain_orientations);
                            new_octahedron_center_x = interface.buffer_north_recv(buf_position, 4);
                            new_octahedron_center_y = interface.buffer_north_recv(buf_position, 5);
                            new_octahedron_center_z = interface.buffer_north_recv(buf_position, 6);
                            new_diagonal_length = interface.buffer_north_recv(buf_position, 7);
                        }
                        else if ((celldata.cell_type(index) == Active) &&
                                 (interface.buffer_north_recv(buf_position, 7) == 0.0)) {
                            celldata.cell_type(index) = Liquid;
                        }
                    }
                    if (place) {
                        // Update this ghost node cell's information with data from other rank
                        grain_id(index) = new_grain_id;
                        interface.octahedron_center(3 * index) = new_octahedron_center_x;
                        interface.octahedron_center(3 * index + 1) = new_octahedron_center_y;
                        interface.octahedron_center(3 * index + 2) = new_octahedron_center_z;
                        int my_orientation = getGrainOrientation(grain_id(index), orientation.n_grain_orientations);
                        interface.diagonal_length(index) = static_cast<float>(new_diagonal_length);
                        // Cell center - note that the Y coordinate is relative to the domain origin to keep the
                        // coordinate system continuous across ranks
                        float xp = coord_x + 0.5;
                        float yp = coord_y + grid.y_offset + 0.5;
                        float zp = coord_z + 0.5;
                        // Calculate critical values at which this active cell leads to the activation of a
                        // neighboring liquid cell
                        interface.calcCritDiagonalLength(index, xp, yp, zp, new_octahedron_center_x,
                                                         new_octahedron_center_y, new_octahedron_center_z,
                                                         my_orientation, orientation.grain_unit_vector);
                        celldata.cell_type(index) = Active;
                    }
                });
        }
    }

    // Reset send buffer data to -1 (used as placeholder) and reset the number of cells stored in the buffers to 0
    interface.resetBuffers();
    // Wait on send requests
    MPI_Waitall(2, send_requests.data(), MPI_STATUSES_IGNORE);
    Kokkos::fence();
}
//*****************************************************************************/
// Jump to the next time step with work to be done, if nothing left to do in the near future
// The cells of interest are active cells, and the view checked for future work is
// melt_time_step Print intermediate output during this jump if PrintIdleMovieFrames = true
template <typename MemorySpace>
void jumpTimeStep(int &cycle, int remaining_cells_of_interest, const int local_temp_solid_cells,
                  Temperature<MemorySpace> &temperature, const Grid &grid, CellData<MemorySpace> &celldata,
                  const int id, const int layernumber, const int np, Orientation<MemorySpace> &orientation, Print print,
                  const double deltat, Interface<MemorySpace> &interface) {

    MPI_Bcast(&remaining_cells_of_interest, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (remaining_cells_of_interest == 0) {
        // If this rank still has cells that will later undergo transformation (LocalIncompleteCells > 0), check when
        // the next solid cells go above the liquidus (remelting) Otherwise, assign the largest possible time step as
        // the next time work needs to be done on the rank
        int next_melt_time_step;
        if (local_temp_solid_cells > 0) {
            Kokkos::parallel_reduce(
                "CheckNextTSForWork", grid.domain_size,
                KOKKOS_LAMBDA(const int &index, int &tempv) {
                    // criteria for a cell to be associated with future work (checking this layer's cells only)
                    if (celldata.cell_type(index) == TempSolid) {
                        int solidification_counter_this_cell = temperature.solidification_event_counter(index);
                        int next_melt_time_step_this_cell = static_cast<int>(
                            temperature.layer_time_temp_history(index, solidification_counter_this_cell, 0));
                        if (next_melt_time_step_this_cell < tempv)
                            tempv = next_melt_time_step_this_cell;
                    }
                },
                Kokkos::Min<int>(next_melt_time_step));
        }
        else
            next_melt_time_step = INT_MAX;

        int global_next_melt_time_step;
        MPI_Allreduce(&next_melt_time_step, &global_next_melt_time_step, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if ((global_next_melt_time_step - cycle) > 5000) {
            // Print current grain misorientations (up to and including the current layer's data) for any of the time
            // steps between now and when melting/solidification occurs again, if the print option for idle frame
            // printing was toggled

            print.printIdleIntralayer(id, np, layernumber, deltat, cycle, grid, celldata, temperature, interface,
                                      orientation, global_next_melt_time_step);
            // Jump to next time step when solidification starts again
            cycle = global_next_melt_time_step - 1;
            if (id == 0)
                std::cout << "Jumping to cycle " << cycle + 1 << std::endl;
        }
    }
}

//*****************************************************************************/
// Prints intermediate code output to stdout (intermediate output collected and printed is different than without
// remelting) and checks to see if solidification is complete in the case where cells can solidify multiple times
template <typename MemorySpace>
void intermediateOutputAndCheck(const int id, const int np, int &cycle, const Grid &grid,
                                int successful_nuc_events_this_rank, int &x_switch, CellData<MemorySpace> &celldata,
                                Temperature<MemorySpace> &temperature, std::string simulation_type,
                                const int layernumber, Orientation<MemorySpace> &orientation, Print print,
                                const double deltat, Interface<MemorySpace> &interface) {

    auto grain_id = celldata.getGrainIDSubview(grid);
    int local_superheated_cells, local_undercooled_cells, local_active_cells, local_temp_solid_cells,
        local_finished_solid_cells;
    Kokkos::parallel_reduce(
        "IntermediateOutput", grid.domain_size,
        KOKKOS_LAMBDA(const int &index, int &sum_superheated, int &sum_undercooled, int &sum_active,
                      int &sum_temp_solid, int &sum_finished_solid) {
            int cell_type_this_cell = celldata.cell_type(index);
            if (cell_type_this_cell == Liquid) {
                int crit_time_step = temperature.getCritTimeStep(index);
                if (crit_time_step > cycle)
                    sum_superheated += 1;
                else
                    sum_undercooled += 1;
            }
            else if (cell_type_this_cell == Active)
                sum_active += 1;
            else if (cell_type_this_cell == TempSolid)
                sum_temp_solid += 1;
            else if (cell_type_this_cell == Solid)
                sum_finished_solid += 1;
        },
        local_superheated_cells, local_undercooled_cells, local_active_cells, local_temp_solid_cells,
        local_finished_solid_cells);

    int global_successful_nuc_events_this_rank = 0;
    int global_superheated_cells, global_undercooled_cells, global_active_cells, global_temp_solid_cells,
        global_finished_solid_cells;
    MPI_Reduce(&local_superheated_cells, &global_superheated_cells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_undercooled_cells, &global_undercooled_cells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_active_cells, &global_active_cells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_temp_solid_cells, &global_temp_solid_cells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_finished_solid_cells, &global_finished_solid_cells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&successful_nuc_events_this_rank, &global_successful_nuc_events_this_rank, 1, MPI_INT, MPI_SUM, 0,
               MPI_COMM_WORLD);

    if (id == 0) {
        std::cout << "cycle = " << cycle << " in layer " << layernumber
                  << ": Superheated liquid cells = " << global_superheated_cells
                  << " Undercooled liquid cells = " << global_undercooled_cells
                  << " Number of nucleation events this layer " << global_successful_nuc_events_this_rank
                  << " cells to undergo at least one more solidification event = " << global_temp_solid_cells
                  << " cells that are finished with solidification = " << global_finished_solid_cells << std::endl;
        if (global_superheated_cells + global_undercooled_cells + global_active_cells + global_temp_solid_cells == 0)
            x_switch = 1;
    }
    MPI_Bcast(&x_switch, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Cells of interest are those currently undergoing a melting-solidification cycle
    int remaining_cells_of_interest = global_active_cells + global_superheated_cells + global_undercooled_cells;
    if ((x_switch == 0) && ((simulation_type == "R") || (simulation_type == "S")))
        jumpTimeStep(cycle, remaining_cells_of_interest, local_temp_solid_cells, temperature, grid, celldata, id,
                     layernumber, np, orientation, print, deltat, interface);
}

//*****************************************************************************/
// Prints intermediate code output to stdout and checks to see the single grain simulation end condition (the grain has
// reached a domain edge) has been satisfied
template <typename ViewTypeInt>
void intermediateOutputAndCheck(const int id, int cycle, const Grid &grid, int &x_switch, ViewTypeInt cell_type) {

    int local_liquid_cells, local_active_cells, local_solid_cells;
    using memory_space = typename ViewTypeInt::memory_space;
    Kokkos::View<bool **, memory_space> edges_reached(Kokkos::ViewAllocateWithoutInitializing("edges_reached"), 3,
                                                      2); // init to false
    Kokkos::deep_copy(edges_reached, false);

    Kokkos::parallel_reduce(
        "IntermediateOutput", grid.domain_size,
        KOKKOS_LAMBDA(const int &index, int &sum_liquid, int &sum_active, int &sum_solid) {
            if (cell_type(index) == Liquid)
                sum_liquid += 1;
            else if (cell_type(index) == Active) {
                sum_active += 1;
                // Did this cell reach a domain edge?
                int coord_x = grid.getCoordX(index);
                int coord_y = grid.getCoordY(index);
                int coord_z = grid.getCoordZ(index);
                int coord_y_global = coord_y + grid.y_offset;
                if (coord_x == 0)
                    edges_reached(0, 0) = true;
                if (coord_x == grid.nx - 1)
                    edges_reached(0, 1) = true;
                if (coord_y_global == 0)
                    edges_reached(1, 0) = true;
                if (coord_y_global == grid.ny - 1)
                    edges_reached(1, 1) = true;
                if (coord_z == 0)
                    edges_reached(2, 0) = true;
                if (coord_z == grid.nz - 1)
                    edges_reached(2, 1) = true;
            }
            else if (cell_type(index) == Solid)
                sum_solid += 1;
        },
        local_liquid_cells, local_active_cells, local_solid_cells);

    int global_liquid_cells, global_active_cells, global_solid_cells;
    MPI_Reduce(&local_liquid_cells, &global_liquid_cells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_active_cells, &global_active_cells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_solid_cells, &global_solid_cells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "cycle = " << cycle << " : Liquid cells = " << global_liquid_cells
                  << " Active cells = " << global_active_cells << " Solid cells = " << global_solid_cells << std::endl;

    // Each rank checks to see if a global domain boundary was reached
    auto edges_reached_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), edges_reached);
    int x_switch_local = 0;
    std::vector<std::string> edge_dims = {"X", "Y", "Z"};
    std::vector<std::string> edge_names = {"Lower", "Upper"};
    for (int edgedim = 0; edgedim < 3; edgedim++) {
        for (int edgename = 0; edgename < 2; edgename++) {
            if (edges_reached_host(edgedim, edgename)) {
                std::cout << edge_names[edgename] << " edge of domain in the " << edge_dims[edgedim]
                          << " direction was reached on rank " << id << " and cycle " << cycle
                          << "; simulation is complete" << std::endl;
                x_switch_local = 1;
            }
        }
    }
    // Simulation ends if a global domain boundary was reached on any rank
    MPI_Allreduce(&x_switch_local, &x_switch, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
}

#endif
