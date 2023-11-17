// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_INTERFACE_HPP
#define EXACA_INTERFACE_HPP

#include "CAcelldata.hpp"
#include "CAconfig.hpp"
#include "CAfunctions.hpp"
#include "CAgrid.hpp"
#include "CAinputs.hpp"
#include "CAparsefiles.hpp"
#include "CAtemperature.hpp"
#include "CAtypes.hpp"
#include "mpi.h"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

// Using for compatibility with device math functions.
using std::max;
using std::min;

// Data representing the active cells at the solid-liquid interrface, including MPI buffers
template <typename MemorySpace>
struct Interface {

    using memory_space = MemorySpace;
    using view_type_buffer = Kokkos::View<float **, memory_space>;
    using view_type_float = Kokkos::View<float *, memory_space>;
    using view_type_int = Kokkos::View<int *, memory_space>;
    using view_type_int_host = typename view_type_int::HostMirror;
    using neighbor_list_type = Kokkos::Array<int, 26>;

    // Using the default exec space for this memory space.
    using execution_space = typename memory_space::execution_space;

    // Size of send/recv buffers before and after the potential resize
    int OldBufSize, BufSize, BufComponents;
    view_type_float DiagonalLength, DOCenter, CritDiagonalLength;
    view_type_buffer BufferSouthSend, BufferNorthSend, BufferSouthRecv, BufferNorthRecv;
    view_type_int SendSizeSouth, SendSizeNorth, SteeringVector, numSteer;
    view_type_int_host SendSizeSouth_Host, SendSizeNorth_Host, numSteer_Host;

    // Neighbor lists
    neighbor_list_type NeighborX, NeighborY, NeighborZ;

    // Parallel dispatch tags.
    struct CellCaptureTag {};

    // Store the system details.
    int np;
    Grid grid;
    InterfacialResponseFunction irf;
    Temperature<memory_space> temperature;
    CellData<memory_space> cellData;
    ViewF GrainUnitVector;
    int NGrainOrientations;

    // Constructor for views and view bounds for current layer
    // Use default initialization to 0 for numSteer_host and numSteer and buffer counts
    Interface(int np_, Grid &grid_, InterfacialResponseFunction &irf_, CellData<memory_space> &cell_,
              Temperature<memory_space> &temp_, ViewF &unit_vector_, int n_orientations_,
              int BufSizeInitialEstimate = 25, int BufComponents_temp = 8)
        : np(np_)
        , grid(grid_)
        , irf(irf_)
        , temperature(temp_)
        , cellData(cell_)
        , GrainUnitVector(unit_vector_)
        , NGrainOrientations(n_orientations_) {
        DiagonalLength = view_type_float(Kokkos::ViewAllocateWithoutInitializing("DiagonalLength"), grid.DomainSize);
        DOCenter = view_type_float(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * grid.DomainSize);
        CritDiagonalLength =
            view_type_float(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength"), 26 * grid.DomainSize);
        BufferSouthSend = view_type_buffer(Kokkos::ViewAllocateWithoutInitializing("BufferSouthSend"),
                                           BufSizeInitialEstimate, BufComponents_temp);
        BufferNorthSend = view_type_buffer(Kokkos::ViewAllocateWithoutInitializing("BufferNorthSend"),
                                           BufSizeInitialEstimate, BufComponents_temp);
        BufferSouthRecv = view_type_buffer(Kokkos::ViewAllocateWithoutInitializing("BufferSouthRecv"),
                                           BufSizeInitialEstimate, BufComponents_temp);
        BufferNorthRecv = view_type_buffer(Kokkos::ViewAllocateWithoutInitializing("BufferNorthRecv"),
                                           BufSizeInitialEstimate, BufComponents_temp);
        SendSizeSouth = view_type_int("SendSizeSouth", 1);
        SendSizeNorth = view_type_int("SendSizeNorth", 1);
        SteeringVector = view_type_int(Kokkos::ViewAllocateWithoutInitializing("SteeringVector"), grid.DomainSize);
        numSteer = view_type_int("SteeringVectorSize", 1);
        SendSizeSouth_Host = view_type_int_host("SendSizeSouth_Host", 1);
        SendSizeNorth_Host = view_type_int_host("SendSizeNorth_Host", 1);
        numSteer_Host = view_type_int_host("SteeringVectorSize_Host", 1);

        // Set initial buffer size to the estimate
        BufSize = BufSizeInitialEstimate;
        // Set number of components in the buffer
        BufComponents = BufComponents_temp;
        // Send/recv buffers for ghost node data should be initialized with -1s in the first index as placeholders for
        // empty positions in the buffer, and with send size counts of 0
        reset_buffers();
        // Initialize neighbor lists for iterating over active cells
        neighbor_list_init();
    }

    // Set first index in send buffers to -1 (placeholder) for all cells in the buffer, and reset the counts of number
    // of cells contained in buffers to 0s
    void reset_buffers() {

        auto BufferNorthSend_local = BufferNorthSend;
        auto BufferSouthSend_local = BufferSouthSend;
        auto SendSizeNorth_local = SendSizeNorth;
        auto SendSizeSouth_local = SendSizeSouth;
        Kokkos::parallel_for(
            "BufferReset", BufSize, KOKKOS_LAMBDA(const int &i) {
                BufferNorthSend_local(i, 0) = -1.0;
                BufferSouthSend_local(i, 0) = -1.0;
            });
        Kokkos::parallel_for(
            "HaloCountReset", 1, KOKKOS_LAMBDA(const int) {
                SendSizeNorth_local(0) = 0;
                SendSizeSouth_local(0) = 0;
            });
    }

    // Intialize neighbor list structures (NeighborX, NeighborY, NeighborZ)
    void neighbor_list_init() {

        // Assignment of neighbors around a cell "X" is as follows (in order of closest to furthest from cell "X")
        // Neighbors 0 through 8 are in the -Y direction
        // Neighbors 9 through 16 are in the XY plane with cell X
        // Neighbors 17 through 25 are in the +Y direction
        NeighborX = {0, 1, -1, 0, 0, -1, 1, -1, 1, 0, 0, 1, -1, 1, -1, 1, -1, 0, 1, -1, 0, 0, 1, -1, 1, -1};
        NeighborY = {-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        NeighborZ = {0, 0, 0, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 0, 0, 0, 0, 0, 1, -1, 1, 1, -1, -1};
    }

    // For the case where all cells solidify once, determine which cells are associated with the "steering vector" of
    // cells that are either active, or becoming active this time step
    void fill_steering_vector_no_remelt(int cycle) {

        // Cells associated with this layer that are not solid type but have passed the liquidus (crit time step) have
        // their undercooling values updated Cells that meet the aforementioned criteria and are active type should be
        // added to the steering vector
        auto CellType = cellData.getCellTypeSubview(grid);
        auto SteeringVector_local = SteeringVector;
        auto numSteer_local = numSteer;
        auto numSteer_Host_local = numSteer_Host;
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
                        SteeringVector_local(Kokkos::atomic_fetch_add(&numSteer_local(0), 1)) = index;
                    }
                }
            });
        Kokkos::deep_copy(numSteer_Host_local, numSteer_local);
    }

    // For the case where cells may melt and solidify multiple times, determine which cells are associated with the
    // "steering vector" of cells that are either active, or becoming active this time step - version with remelting
    void fill_steering_vector_remelt(int cycle) {

        auto CellType = cellData.getCellTypeSubview(grid);
        auto GrainID = cellData.getGrainIDSubview(grid);
        auto SteeringVector_local = SteeringVector;
        auto numSteer_local = numSteer;
        auto numSteer_Host_local = numSteer_Host;
        auto NeighborX_local = NeighborX;
        auto NeighborY_local = NeighborY;
        auto NeighborZ_local = NeighborZ;
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
                            // This cell either hasn't started or hasn't finished the previous solidification event, but
                            // has undergone melting - increment the solidification counter to move on to the next
                            // melt-solidification event
                            temperature.update_solidification_counter(index);
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
                            int neighbor_coord_x = coord_x + NeighborX_local[l];
                            int neighbor_coord_y = coord_y + NeighborY_local[l];
                            int neighbor_coord_z = coord_z + NeighborZ_local[l];
                            if ((neighbor_coord_x >= 0) && (neighbor_coord_x < grid.nx) && (neighbor_coord_y >= 0) &&
                                (neighbor_coord_y < grid.ny_local) && (neighbor_coord_z < grid.nz_layer) &&
                                (neighbor_coord_z >= 0)) {
                                int neighbor_index =
                                    grid.get1Dindex(neighbor_coord_x, neighbor_coord_y, neighbor_coord_z);
                                if (CellType(neighbor_index) == Active) {
                                    // Mark adjacent active cells to this as cells that should be converted into liquid,
                                    // as they are more likely heating than cooling
                                    CellType(neighbor_index) = FutureLiquid;
                                    SteeringVector_local(Kokkos::atomic_fetch_add(&numSteer_local(0), 1)) =
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
                            SteeringVector_local(Kokkos::atomic_fetch_add(&numSteer_local(0), 1)) = index;
                        }
                    }
                    else if ((atCritTime) && (cellType == Liquid) && (GrainID(index) != 0)) {
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
                            int neighbor_coord_x = coord_x + NeighborX_local[l];
                            int neighbor_coord_y = coord_y + NeighborY_local[l];
                            int neighbor_coord_z = coord_z + NeighborZ_local[l];
                            if ((neighbor_coord_x >= 0) && (neighbor_coord_x < grid.nx) && (neighbor_coord_y >= 0) &&
                                (neighbor_coord_y < grid.ny_local) && (neighbor_coord_z < grid.nz_layer) &&
                                (neighbor_coord_z >= 0)) {
                                int neighbor_index =
                                    grid.get1Dindex(neighbor_coord_x, neighbor_coord_y, neighbor_coord_z);
                                if ((CellType(neighbor_index) == TempSolid) || (CellType(neighbor_index) == Solid) ||
                                    (coord_z == 0)) {
                                    // Cell activation to be performed as part of steering vector
                                    l = 26;
                                    SteeringVector_local(Kokkos::atomic_fetch_add(&numSteer_local(0), 1)) = index;
                                    CellType(index) = FutureActive; // this cell cannot be captured - is being activated
                                }
                            }
                        }
                    }
                }
            });
        Kokkos::fence();

        // Copy size of steering vector (containing positions of undercooled liquid/active cells) to the host
        Kokkos::deep_copy(numSteer_Host_local, numSteer_local);
    }

    // Decentered octahedron algorithm for the capture of new interface cells by grains
    void cell_capture() {

        // Loop over list of active and soon-to-be active cells, potentially performing cell capture events and updating
        // cell types
        auto policy = Kokkos::RangePolicy<CellCaptureTag>(0, numSteer_Host(0));
        Kokkos::parallel_for("CellCapture", policy, *this);
        Kokkos::fence();
    }

    void operator()(CellCaptureTag, const int num) const {
        numSteer(0) = 0;
        int index = SteeringVector(num);
        // Get the x, y, z coordinates of the cell on this MPI rank
        int coord_x = grid.getCoordX(index);
        int coord_y = grid.getCoordY(index);
        int coord_z = grid.getCoordZ(index);

        // Get subviews for this layer
        auto CellType = cellData.getCellTypeSubview(grid);
        auto GrainID = cellData.getGrainIDSubview(grid);

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
                int neighbor_coord_x = coord_x + NeighborX[l];
                int neighbor_coord_y = coord_y + NeighborY[l];
                int neighbor_coord_z = coord_z + NeighborZ[l];
                // Check if neighbor is in bounds
                if ((neighbor_coord_x >= 0) && (neighbor_coord_x < grid.nx) && (neighbor_coord_y >= 0) &&
                    (neighbor_coord_y < grid.ny) && (neighbor_coord_z < grid.nz_layer) && (neighbor_coord_z >= 0)) {
                    int neighbor_index = grid.get1Dindex(neighbor_coord_x, neighbor_coord_y, neighbor_coord_z);
                    if (CellType(neighbor_index) == Liquid)
                        DeactivateCell = false;
                    // Capture of cell located at "NeighborD3D1ConvPosition" if this condition is satisfied
                    if ((DiagonalLength(index) >= CritDiagonalLength(26 * index + l)) &&
                        (CellType(neighbor_index) == Liquid)) {
                        // Use of atomic_compare_exchange
                        // (https://github.com/kokkos/kokkos/wiki/Kokkos%3A%3Aatomic_compare_exchange) old_val =
                        // atomic_compare_exchange(ptr_to_value,comparison_value, new_value); Atomicly sets the value at
                        // the address given by ptr_to_value to new_value if the current value at ptr_to_value is equal
                        // to comparison_value Returns the previously stored value at the address independent on whether
                        // the exchange has happened. If this cell's is a liquid cell, change it to "TemporaryUpdate"
                        // type and return a value of "liquid" If this cell has already been changed to
                        // "TemporaryUpdate" type, return a value of "0"
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
                            // Note that the Y coordinate is relative to the domain origin to keep the coordinate system
                            // continuous across ranks
                            float xp = coord_x + NeighborX[l] + 0.5;
                            float yp = coord_y + grid.y_offset + NeighborY[l] + 0.5;
                            float zp = coord_z + NeighborZ[l] + 0.5;

                            // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the image of
                            // the center of the new cell
                            float x0 = xp - cxold;
                            float y0 = yp - cyold;
                            float z0 = zp - czold;

                            // mag0 is the magnitude of (x0,y0,z0)
                            float mag0 = sqrtf(x0 * x0 + y0 * y0 + z0 * z0);

                            // Calculate unit vectors for the octahedron that intersect the new cell center
                            float Angle1 =
                                (GrainUnitVector(9 * MyOrientation) * x0 + GrainUnitVector(9 * MyOrientation + 1) * y0 +
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

                            DiagonalLength(neighbor_index) = NewODiagL;
                            // Calculate coordinates of new decentered octahedron center
                            float CaptDiag[3];
                            CaptDiag[0] = xc - cxold;
                            CaptDiag[1] = yc - cyold;
                            CaptDiag[2] = zc - czold;

                            float CaptDiagMagnitude =
                                sqrt(CaptDiag[0] * CaptDiag[0] + CaptDiag[1] * CaptDiag[1] + CaptDiag[2] * CaptDiag[2]);
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
                            calc_crit_diagonal_length(neighbor_index, xp, yp, zp, cx, cy, cz, NeighborX, NeighborY,
                                                      NeighborZ, MyOrientation, GrainUnitVector, CritDiagonalLength);

                            if (np > 1) {
                                // TODO: Test loading ghost nodes in a separate kernel, potentially adopting
                                // this change if the slowdown is minor
                                int GhostGID = h;
                                float GhostDOCX = cx;
                                float GhostDOCY = cy;
                                float GhostDOCZ = cz;
                                float GhostDL = NewODiagL;
                                // Collect data for the ghost nodes, if necessary
                                // Data loaded into the ghost nodes is for the cell that was just captured
                                bool DataFitsInBuffer =
                                    load_ghost_nodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth,
                                                     SendSizeSouth, grid.ny, neighbor_coord_x, neighbor_coord_y,
                                                     neighbor_coord_z, grid.AtNorthBoundary, grid.AtSouthBoundary,
                                                     BufferSouthSend, BufferNorthSend, NGrainOrientations, BufSize);
                                if (!(DataFitsInBuffer)) {
                                    // This cell's data did not fit in the buffer with current size BufSize -
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
            create_new_octahedron(index, DiagonalLength, DOCenter, coord_x, coord_y, grid.y_offset, coord_z);
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
            calc_crit_diagonal_length(index, cx, cy, cz, cx, cy, cz, NeighborX, NeighborY, NeighborZ, MyOrientation,
                                      GrainUnitVector, CritDiagonalLength);
            if (np > 1) {
                // TODO: Test loading ghost nodes in a separate kernel, potentially adopting this change if the
                // slowdown is minor
                int GhostGID = MyGrainID;
                float GhostDOCX = cx;
                float GhostDOCY = cy;
                float GhostDOCZ = cz;
                float GhostDL = 0.01;
                // Collect data for the ghost nodes, if necessary
                bool DataFitsInBuffer =
                    load_ghost_nodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth, SendSizeSouth,
                                     grid.ny, coord_x, coord_y, coord_z, grid.AtNorthBoundary, grid.AtSouthBoundary,
                                     BufferSouthSend, BufferNorthSend, NGrainOrientations, BufSize);
                if (!(DataFitsInBuffer)) {
                    // This cell's data did not fit in the buffer with current size BufSize - mark with
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
            bool DataFitsInBuffer =
                load_ghost_nodes(-1, -1.0, -1.0, -1.0, 0.0, SendSizeNorth, SendSizeSouth, grid.ny, coord_x, coord_y,
                                 coord_z, grid.AtNorthBoundary, grid.AtSouthBoundary, BufferSouthSend, BufferNorthSend,
                                 NGrainOrientations, BufSize);
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
    }

    void check_buffers(int id, Grid &grid, CellData<memory_space> &cellData, int NGrainOrientations) {
        // Count the number of cells' in halo regions where the data did not fit into the send buffers
        // Reduce across all ranks, as the same BufSize should be maintained across all ranks
        // If any rank overflowed its buffer size, resize all buffers to the new size plus 10% padding
        OldBufSize = BufSize;
        BufSize = resize_buffers();
        if (OldBufSize != BufSize) {
            if (id == 0)
                std::cout << "Resized number of cells stored in send/recv buffers from " << OldBufSize << " to "
                          << BufSize << std::endl;
            refill_buffers(grid, cellData, NGrainOrientations);
        }
    }

    // Count the number of cells' in halo regions where the data did not fit into the send buffers, and resize the
    // buffers if necessary, returning the new buffer size
    int resize_buffers(int NumCellsBufferPadding = 25) {

        int NewBufSize;
        SendSizeNorth_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), SendSizeNorth);
        SendSizeSouth_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), SendSizeSouth);
        int max_count_local = max(SendSizeNorth_Host(0), SendSizeSouth_Host(0));
        int max_count_global;
        MPI_Allreduce(&max_count_local, &max_count_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (max_count_global > OldBufSize) {
            // Increase buffer size to fit all data
            // Add numcells_buffer_padding (defaults to 25) cells as additional padding
            NewBufSize = max_count_global + NumCellsBufferPadding;
            Kokkos::resize(BufferNorthSend, NewBufSize, BufComponents);
            Kokkos::resize(BufferSouthSend, NewBufSize, BufComponents);
            Kokkos::resize(BufferNorthRecv, NewBufSize, BufComponents);
            Kokkos::resize(BufferSouthRecv, NewBufSize, BufComponents);
            // Reset count variables on device to the old buffer size
            auto SendSizeNorth_local = SendSizeNorth;
            auto SendSizeSouth_local = SendSizeSouth;
            auto OldBufSize_local = OldBufSize;
            auto BufferNorthSend_local = BufferNorthSend;
            auto BufferSouthSend_local = BufferSouthSend;
            auto BufComponents_local = BufComponents;
            Kokkos::parallel_for(
                "ResetCounts", 1, KOKKOS_LAMBDA(const int &) {
                    SendSizeNorth_local(0) = OldBufSize_local;
                    SendSizeSouth_local(0) = OldBufSize_local;
                });
            // Set -1 values for the new (currently empty) positions in the resized buffer
            Kokkos::parallel_for(
                "InitNewBufCapacity", Kokkos::RangePolicy<>(OldBufSize_local, NewBufSize),
                KOKKOS_LAMBDA(const int &BufPosition) {
                    for (int BufComp = 0; BufComp < BufComponents_local; BufComp++) {
                        BufferNorthSend_local(BufPosition, BufComp) = -1.0;
                        BufferSouthSend_local(BufPosition, BufComp) = -1.0;
                    }
                });
        }
        else
            NewBufSize = OldBufSize;
        return NewBufSize;
    }

    // Refill the buffers as necessary starting from the old count size, using the data from cells marked with type
    // ActiveFailedBufferLoad
    void refill_buffers(Grid &grid, CellData<memory_space> &cellData, int NGrainOrientations) {

        auto CellType = cellData.getCellTypeSubview(grid);
        auto GrainID = cellData.getGrainIDSubview(grid);
        auto DOCenter_local = DOCenter;
        auto DiagonalLength_local = DiagonalLength;
        auto SendSizeNorth_local = SendSizeNorth;
        auto SendSizeSouth_local = SendSizeSouth;
        auto BufferSouthSend_local = BufferSouthSend;
        auto BufferNorthSend_local = BufferNorthSend;
        auto BufSize_local = BufSize;
        Kokkos::parallel_for(
            "FillSendBuffersOverflow", grid.nx, KOKKOS_LAMBDA(const int &coord_x) {
                for (int coord_z = 0; coord_z < grid.nz_layer; coord_z++) {
                    int CellCoordinateSouth = grid.get1Dindex(coord_x, 1, coord_z);
                    int CellCoordinateNorth = grid.get1Dindex(coord_x, grid.ny_local - 2, coord_z);
                    if (CellType(CellCoordinateSouth) == ActiveFailedBufferLoad) {
                        int GhostGID = GrainID(CellCoordinateSouth);
                        float GhostDOCX = DOCenter_local(3 * CellCoordinateSouth);
                        float GhostDOCY = DOCenter_local(3 * CellCoordinateSouth + 1);
                        float GhostDOCZ = DOCenter_local(3 * CellCoordinateSouth + 2);
                        float GhostDL = DiagonalLength_local(CellCoordinateSouth);
                        // Collect data for the ghost nodes, if necessary
                        // Data loaded into the ghost nodes is for the cell that was just captured
                        bool DataFitsInBuffer =
                            load_ghost_nodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth_local,
                                             SendSizeSouth_local, grid.ny_local, coord_x, 1, coord_z,
                                             grid.AtNorthBoundary, grid.AtSouthBoundary, BufferSouthSend_local,
                                             BufferNorthSend_local, NGrainOrientations, BufSize_local);
                        CellType(CellCoordinateSouth) = Active;
                        // If data doesn't fit in the buffer after the resize, warn that buffer data may have been lost
                        if (!(DataFitsInBuffer))
                            printf("Error: Send/recv buffer resize failed to include all necessary data, predicted "
                                   "results at MPI processor boundaries may be inaccurate\n");
                    }
                    else if (CellType(CellCoordinateSouth) == LiquidFailedBufferLoad) {
                        // Dummy values for first 4 arguments (Grain ID and octahedron center coordinates), 0 for
                        // diagonal length
                        bool DataFitsInBuffer = load_ghost_nodes(
                            -1, -1.0, -1.0, -1.0, 0.0, SendSizeNorth_local, SendSizeSouth_local, grid.ny_local, coord_x,
                            1, coord_z, grid.AtNorthBoundary, grid.AtSouthBoundary, BufferSouthSend_local,
                            BufferNorthSend_local, NGrainOrientations, BufSize_local);
                        CellType(CellCoordinateSouth) = Liquid;
                        // If data doesn't fit in the buffer after the resize, warn that buffer data may have been lost
                        if (!(DataFitsInBuffer))
                            printf("Error: Send/recv buffer resize failed to include all necessary data, predicted "
                                   "results at MPI processor boundaries may be inaccurate\n");
                    }
                    if (CellType(CellCoordinateNorth) == ActiveFailedBufferLoad) {
                        int GhostGID = GrainID(CellCoordinateNorth);
                        float GhostDOCX = DOCenter(3 * CellCoordinateNorth);
                        float GhostDOCY = DOCenter(3 * CellCoordinateNorth + 1);
                        float GhostDOCZ = DOCenter(3 * CellCoordinateNorth + 2);
                        float GhostDL = DiagonalLength(CellCoordinateNorth);
                        // Collect data for the ghost nodes, if necessary
                        // Data loaded into the ghost nodes is for the cell that was just captured
                        bool DataFitsInBuffer =
                            load_ghost_nodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth_local,
                                             SendSizeSouth_local, grid.ny_local, coord_x, grid.ny_local - 2, coord_z,
                                             grid.AtNorthBoundary, grid.AtSouthBoundary, BufferSouthSend_local,
                                             BufferNorthSend_local, NGrainOrientations, BufSize_local);
                        CellType(CellCoordinateNorth) = Active;
                        // If data doesn't fit in the buffer after the resize, warn that buffer data may have been lost
                        if (!(DataFitsInBuffer))
                            printf("Error: Send/recv buffer resize failed to include all necessary data, predicted "
                                   "results at MPI processor boundaries may be inaccurate\n");
                    }
                    else if (CellType(CellCoordinateNorth) == LiquidFailedBufferLoad) {
                        // Dummy values for first 4 arguments (Grain ID and octahedron center coordinates), 0 for
                        // diagonal length
                        bool DataFitsInBuffer = load_ghost_nodes(
                            -1, -1.0, -1.0, -1.0, 0.0, SendSizeNorth_local, SendSizeSouth_local, grid.ny_local, coord_x,
                            grid.ny_local - 2, coord_z, grid.AtNorthBoundary, grid.AtSouthBoundary,
                            BufferSouthSend_local, BufferNorthSend_local, NGrainOrientations, BufSize_local);
                        CellType(CellCoordinateNorth) = Liquid;
                        // If data doesn't fit in the buffer after the resize, warn that buffer data may have been lost
                        if (!(DataFitsInBuffer))
                            printf("Error: Send/recv buffer resize failed to include all necessary data, predicted "
                                   "results at MPI processor boundaries may be inaccurate\n");
                    }
                }
            });
        Kokkos::fence();
    }

    // 1D domain decomposition: update ghost nodes with new cell data from Nucleation and CellCapture routines
    void halo_update(int, int, Grid &grid, CellData<device_memory_space> &cellData, int NGrainOrientations,
                     ViewF GrainUnitVector) {

        std::vector<MPI_Request> SendRequests(2, MPI_REQUEST_NULL);
        std::vector<MPI_Request> RecvRequests(2, MPI_REQUEST_NULL);

        // Send data to each other rank (MPI_Isend)
        MPI_Isend(BufferSouthSend.data(), BufComponents * BufSize, MPI_FLOAT, grid.NeighborRank_South, 0,
                  MPI_COMM_WORLD, &SendRequests[0]);
        MPI_Isend(BufferNorthSend.data(), BufComponents * BufSize, MPI_FLOAT, grid.NeighborRank_North, 0,
                  MPI_COMM_WORLD, &SendRequests[1]);

        // Receive buffers for all neighbors (MPI_Irecv)
        MPI_Irecv(BufferSouthRecv.data(), BufComponents * BufSize, MPI_FLOAT, grid.NeighborRank_South, 0,
                  MPI_COMM_WORLD, &RecvRequests[0]);
        MPI_Irecv(BufferNorthRecv.data(), BufComponents * BufSize, MPI_FLOAT, grid.NeighborRank_North, 0,
                  MPI_COMM_WORLD, &RecvRequests[1]);

        // unpack in any order
        bool unpack_complete = false;
        auto CellType = cellData.getCellTypeSubview(grid);
        auto GrainID = cellData.getGrainIDSubview(grid);
        auto DiagonalLength_local = DiagonalLength;
        auto CritDiagonalLength_local = CritDiagonalLength;
        auto DOCenter_local = DOCenter;
        auto BufferSouthRecv_local = BufferSouthRecv;
        auto BufferNorthRecv_local = BufferNorthRecv;
        auto BufSize_local = BufSize;
        auto NeighborX_local = NeighborX;
        auto NeighborY_local = NeighborY;
        auto NeighborZ_local = NeighborZ;
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
                    "BufferUnpack", BufSize_local, KOKKOS_LAMBDA(const int &BufPosition) {
                        int coord_x, coord_y, coord_z, index, NewGrainID;
                        float DOCenterX, DOCenterY, DOCenterZ, NewDiagonalLength;
                        bool Place = false;
                        // Which rank was the data received from? Is there valid data at this position in the buffer
                        // (i.e., not set to -1.0)?
                        if ((unpack_index == 0) && (BufferSouthRecv_local(BufPosition, 0) != -1.0) &&
                            (grid.NeighborRank_South != MPI_PROC_NULL)) {
                            // Data receieved from South
                            coord_x = static_cast<int>(BufferSouthRecv_local(BufPosition, 0));
                            coord_y = 0;
                            coord_z = static_cast<int>(BufferSouthRecv_local(BufPosition, 1));
                            index = grid.get1Dindex(coord_x, coord_y, coord_z);
                            // Two possibilities: buffer data with non-zero diagonal length was loaded, and a liquid
                            // cell may have to be updated to active - or zero diagonal length data was loaded, and an
                            // active cell may have to be updated to liquid
                            if (CellType(index) == Liquid) {
                                Place = true;
                                int MyGrainOrientation = static_cast<int>(BufferSouthRecv_local(BufPosition, 2));
                                int MyGrainNumber = static_cast<int>(BufferSouthRecv_local(BufPosition, 3));
                                NewGrainID = getGrainID(NGrainOrientations, MyGrainOrientation, MyGrainNumber);
                                DOCenterX = BufferSouthRecv_local(BufPosition, 4);
                                DOCenterY = BufferSouthRecv_local(BufPosition, 5);
                                DOCenterZ = BufferSouthRecv_local(BufPosition, 6);
                                NewDiagonalLength = BufferSouthRecv_local(BufPosition, 7);
                            }
                            else if ((CellType(index) == Active) && (BufferSouthRecv_local(BufPosition, 7) == 0.0)) {
                                CellType(index) = Liquid;
                            }
                        }
                        else if ((unpack_index == 1) && (BufferNorthRecv_local(BufPosition, 0) != -1.0) &&
                                 (grid.NeighborRank_North != MPI_PROC_NULL)) {
                            // Data received from North
                            coord_x = static_cast<int>(BufferNorthRecv_local(BufPosition, 0));
                            coord_y = grid.ny_local - 1;
                            coord_z = static_cast<int>(BufferNorthRecv_local(BufPosition, 1));
                            index = grid.get1Dindex(coord_x, coord_y, coord_z);
                            // Two possibilities: buffer data with non-zero diagonal length was loaded, and a liquid
                            // cell may have to be updated to active - or zero diagonal length data was loaded, and an
                            // active cell may have to be updated to liquid
                            if (CellType(index) == Liquid) {
                                Place = true;
                                int MyGrainOrientation = static_cast<int>(BufferNorthRecv_local(BufPosition, 2));
                                int MyGrainNumber = static_cast<int>(BufferNorthRecv_local(BufPosition, 3));
                                NewGrainID = getGrainID(NGrainOrientations, MyGrainOrientation, MyGrainNumber);
                                DOCenterX = BufferNorthRecv_local(BufPosition, 4);
                                DOCenterY = BufferNorthRecv_local(BufPosition, 5);
                                DOCenterZ = BufferNorthRecv_local(BufPosition, 6);
                                NewDiagonalLength = BufferNorthRecv_local(BufPosition, 7);
                            }
                            else if ((CellType(index) == Active) && (BufferNorthRecv_local(BufPosition, 7) == 0.0)) {
                                CellType(index) = Liquid;
                            }
                        }
                        if (Place) {
                            // Update this ghost node cell's information with data from other rank
                            GrainID(index) = NewGrainID;
                            DOCenter_local(3 * index) = DOCenterX;
                            DOCenter_local(3 * index + 1) = DOCenterY;
                            DOCenter_local(3 * index + 2) = DOCenterZ;
                            int MyOrientation = getGrainOrientation(GrainID(index), NGrainOrientations);
                            DiagonalLength_local(index) = static_cast<float>(NewDiagonalLength);
                            // Cell center - note that the Y coordinate is relative to the domain origin to keep the
                            // coordinate system continuous across ranks
                            double xp = coord_x + 0.5;
                            double yp = coord_y + grid.y_offset + 0.5;
                            double zp = coord_z + 0.5;
                            // Calculate critical values at which this active cell leads to the activation of a
                            // neighboring liquid cell
                            calc_crit_diagonal_length(index, xp, yp, zp, DOCenterX, DOCenterY, DOCenterZ,
                                                      NeighborX_local, NeighborY_local, NeighborZ_local, MyOrientation,
                                                      GrainUnitVector, CritDiagonalLength_local);
                            CellType(index) = Active;
                        }
                    });
            }
        }

        // Reset send buffer data to -1 (used as placeholder) and reset the number of cells stored in the buffers to 0
        reset_buffers();
        // Wait on send requests
        MPI_Waitall(2, SendRequests.data(), MPI_STATUSES_IGNORE);
        Kokkos::fence();
    }

    // Resize and reinitialize structs governing the active cells before the next layer of a multilayer problem
    void init_next_layer(int DomainSize) {

        // Realloc steering vector as LocalActiveDomainSize may have changed (old values aren't needed)
        Kokkos::realloc(SteeringVector, DomainSize);

        // Realloc active cell data structure and halo regions on device (old values not needed)
        Kokkos::realloc(DiagonalLength, DomainSize);
        Kokkos::realloc(DOCenter, 3 * DomainSize);
        Kokkos::realloc(CritDiagonalLength, 26 * DomainSize);

        // Reset active cell data structures on device
        Kokkos::deep_copy(DiagonalLength, 0);
        Kokkos::deep_copy(DOCenter, 0);
        Kokkos::deep_copy(CritDiagonalLength, 0);
    }

    // Assign octahedron a small initial size, and a center location
    // Note that the Y coordinate is relative to the domain origin to keep the coordinate system continuous across ranks
    KOKKOS_INLINE_FUNCTION
    void create_new_octahedron(const int index, view_type_float DiagonalLength_local, view_type_float DOCenter_local,
                               const int coord_x, const int coord_y, const int y_offset, const int coord_z) const {
        DiagonalLength_local(index) = 0.01;
        DOCenter_local(3 * index) = coord_x + 0.5;
        DOCenter_local(3 * index + 1) = coord_y + y_offset + 0.5;
        DOCenter_local(3 * index + 2) = coord_z + 0.5;
    }

    // For the newly active cell located at 1D array position D3D1ConvPosition (3D center coordinate of xp, yp, zp),
    // update CritDiagonalLength values for cell capture of neighboring cells. The octahedron has a center located at
    // (cx, cy, cz) Note that yp and cy are relative to the domain origin to keep the coordinate system continuous
    // across ranks
    template <typename ViewType>
    KOKKOS_INLINE_FUNCTION void
    calc_crit_diagonal_length(int index, float xp, float yp, float zp, float cx, float cy, float cz,
                              neighbor_list_type NeighborX_local, neighbor_list_type NeighborY_local,
                              neighbor_list_type NeighborZ_local, int MyOrientation, ViewType GrainUnitVector,
                              view_type_float CritDiagonalLength_local) const {
        // Calculate critical octahedron diagonal length to activate nearest neighbor.
        // First, calculate the unique planes (4) associated with all octahedron faces (8)
        // Then just look at distance between face and the point of interest (cell center of
        // neighbor). The critical diagonal length will be the maximum of these (since all other
        // planes will have passed over the point by then
        // ... meaning it must be in the octahedron)
        double Fx[4], Fy[4], Fz[4];

        Fx[0] = GrainUnitVector(9 * MyOrientation) + GrainUnitVector(9 * MyOrientation + 3) +
                GrainUnitVector(9 * MyOrientation + 6);
        Fx[1] = GrainUnitVector(9 * MyOrientation) - GrainUnitVector(9 * MyOrientation + 3) +
                GrainUnitVector(9 * MyOrientation + 6);
        Fx[2] = GrainUnitVector(9 * MyOrientation) + GrainUnitVector(9 * MyOrientation + 3) -
                GrainUnitVector(9 * MyOrientation + 6);
        Fx[3] = GrainUnitVector(9 * MyOrientation) - GrainUnitVector(9 * MyOrientation + 3) -
                GrainUnitVector(9 * MyOrientation + 6);

        Fy[0] = GrainUnitVector(9 * MyOrientation + 1) + GrainUnitVector(9 * MyOrientation + 4) +
                GrainUnitVector(9 * MyOrientation + 7);
        Fy[1] = GrainUnitVector(9 * MyOrientation + 1) - GrainUnitVector(9 * MyOrientation + 4) +
                GrainUnitVector(9 * MyOrientation + 7);
        Fy[2] = GrainUnitVector(9 * MyOrientation + 1) + GrainUnitVector(9 * MyOrientation + 4) -
                GrainUnitVector(9 * MyOrientation + 7);
        Fy[3] = GrainUnitVector(9 * MyOrientation + 1) - GrainUnitVector(9 * MyOrientation + 4) -
                GrainUnitVector(9 * MyOrientation + 7);

        Fz[0] = GrainUnitVector(9 * MyOrientation + 2) + GrainUnitVector(9 * MyOrientation + 5) +
                GrainUnitVector(9 * MyOrientation + 8);
        Fz[1] = GrainUnitVector(9 * MyOrientation + 2) - GrainUnitVector(9 * MyOrientation + 5) +
                GrainUnitVector(9 * MyOrientation + 8);
        Fz[2] = GrainUnitVector(9 * MyOrientation + 2) + GrainUnitVector(9 * MyOrientation + 5) -
                GrainUnitVector(9 * MyOrientation + 8);
        Fz[3] = GrainUnitVector(9 * MyOrientation + 2) - GrainUnitVector(9 * MyOrientation + 5) -
                GrainUnitVector(9 * MyOrientation + 8);

        for (int n = 0; n < 26; n++) {
            float x0 = xp + NeighborX_local[n] - cx;
            float y0 = yp + NeighborY_local[n] - cy;
            float z0 = zp + NeighborZ_local[n] - cz;
            float D0 = x0 * Fx[0] + y0 * Fy[0] + z0 * Fz[0];
            float D1 = x0 * Fx[1] + y0 * Fy[1] + z0 * Fz[1];
            float D2 = x0 * Fx[2] + y0 * Fy[2] + z0 * Fz[2];
            float D3 = x0 * Fx[3] + y0 * Fy[3] + z0 * Fz[3];
            float Dfabs = fmax(fmax(fabs(D0), fabs(D1)), fmax(fabs(D2), fabs(D3)));
            CritDiagonalLength_local(26 * index + n) = Dfabs;
        }
    }

    // Load data (GrainID, DOCenter, DiagonalLength) into ghost nodes if the given RankY is associated with a 1D halo
    // region Uses check to ensure that the buffer position does not reach the buffer size - if it does, return false
    // (otherwise return true) but keep incrementing the send size counters for use resizing the buffers in the future
    KOKKOS_INLINE_FUNCTION
    bool load_ghost_nodes(const int GhostGID, const float GhostDOCX, const float GhostDOCY, const float GhostDOCZ,
                          const float GhostDL, view_type_int SendSizeNorth_local, view_type_int SendSizeSouth_local,
                          const int ny_local, const int coord_x, const int coord_y, const int coord_z,
                          const bool AtNorthBoundary, const bool AtSouthBoundary,
                          view_type_buffer BufferSouthSend_local, view_type_buffer BufferNorthSend_local,
                          int NGrainOrientations, int BufSize_local) const {
        bool DataFitsInBuffer = true;
        if ((coord_y == 1) && (!(AtSouthBoundary))) {
            int GNPositionSouth = Kokkos::atomic_fetch_add(&SendSizeSouth_local(0), 1);
            if (GNPositionSouth >= BufSize_local)
                DataFitsInBuffer = false;
            else {
                BufferSouthSend_local(GNPositionSouth, 0) = static_cast<float>(coord_x);
                BufferSouthSend_local(GNPositionSouth, 1) = static_cast<float>(coord_z);
                BufferSouthSend_local(GNPositionSouth, 2) =
                    static_cast<float>(getGrainOrientation(GhostGID, NGrainOrientations, false));
                BufferSouthSend_local(GNPositionSouth, 3) =
                    static_cast<float>(getGrainNumber(GhostGID, NGrainOrientations));
                BufferSouthSend_local(GNPositionSouth, 4) = GhostDOCX;
                BufferSouthSend_local(GNPositionSouth, 5) = GhostDOCY;
                BufferSouthSend_local(GNPositionSouth, 6) = GhostDOCZ;
                BufferSouthSend_local(GNPositionSouth, 7) = GhostDL;
            }
        }
        else if ((coord_y == ny_local - 2) && (!(AtNorthBoundary))) {
            int GNPositionNorth = Kokkos::atomic_fetch_add(&SendSizeNorth_local(0), 1);
            if (GNPositionNorth >= BufSize_local)
                DataFitsInBuffer = false;
            else {
                BufferNorthSend_local(GNPositionNorth, 0) = static_cast<float>(coord_x);
                BufferNorthSend_local(GNPositionNorth, 1) = static_cast<float>(coord_z);
                BufferNorthSend_local(GNPositionNorth, 2) =
                    static_cast<float>(getGrainOrientation(GhostGID, NGrainOrientations, false));
                BufferNorthSend_local(GNPositionNorth, 3) =
                    static_cast<float>(getGrainNumber(GhostGID, NGrainOrientations));
                BufferNorthSend_local(GNPositionNorth, 4) = GhostDOCX;
                BufferNorthSend_local(GNPositionNorth, 5) = GhostDOCY;
                BufferNorthSend_local(GNPositionNorth, 6) = GhostDOCZ;
                BufferNorthSend_local(GNPositionNorth, 7) = GhostDL;
            }
        }
        return DataFitsInBuffer;
    }
};

#endif
