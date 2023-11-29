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
    struct RefillBuffersTag {};

    // Constructor for views and view bounds for current layer
    // Use default initialization to 0 for numSteer_host and numSteer and buffer counts
    Interface(int DomainSize, int BufSizeInitialEstimate = 25, int BufComponents_temp = 8)
        : DiagonalLength(view_type_float(Kokkos::ViewAllocateWithoutInitializing("DiagonalLength"), DomainSize))
        , DOCenter(view_type_float(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * DomainSize))
        , CritDiagonalLength(
              view_type_float(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength"), 26 * DomainSize))
        , BufferSouthSend(view_type_buffer(Kokkos::ViewAllocateWithoutInitializing("BufferSouthSend"),
                                           BufSizeInitialEstimate, BufComponents_temp))
        , BufferNorthSend(view_type_buffer(Kokkos::ViewAllocateWithoutInitializing("BufferNorthSend"),
                                           BufSizeInitialEstimate, BufComponents_temp))
        , BufferSouthRecv(view_type_buffer(Kokkos::ViewAllocateWithoutInitializing("BufferSouthRecv"),
                                           BufSizeInitialEstimate, BufComponents_temp))
        , BufferNorthRecv(view_type_buffer(Kokkos::ViewAllocateWithoutInitializing("BufferNorthRecv"),
                                           BufSizeInitialEstimate, BufComponents_temp))
        , SendSizeSouth(view_type_int("SendSizeSouth", 1))
        , SendSizeNorth(view_type_int("SendSizeNorth", 1))
        , SteeringVector(view_type_int(Kokkos::ViewAllocateWithoutInitializing("SteeringVector"), DomainSize))
        , numSteer(view_type_int("SteeringVectorSize", 1))
        , SendSizeSouth_Host(view_type_int_host("SendSizeSouth_Host", 1))
        , SendSizeNorth_Host(view_type_int_host("SendSizeNorth_Host", 1))
        , numSteer_Host(view_type_int_host("SteeringVectorSize_Host", 1)) {

        // Set initial buffer size to the estimate
        BufSize = BufSizeInitialEstimate;
        OldBufSize = BufSize;
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
    void create_new_octahedron(const int index, const int coord_x, const int coord_y, const int y_offset,
                               const int coord_z) const {
        DiagonalLength(index) = 0.01;
        DOCenter(3 * index) = coord_x + 0.5;
        DOCenter(3 * index + 1) = coord_y + y_offset + 0.5;
        DOCenter(3 * index + 2) = coord_z + 0.5;
    }

    // For the newly active cell located at 1D array position D3D1ConvPosition (3D center coordinate of xp, yp, zp),
    // update CritDiagonalLength values for cell capture of neighboring cells. The octahedron has a center located at
    // (cx, cy, cz) Note that yp and cy are relative to the domain origin to keep the coordinate system continuous
    // across ranks
    template <typename ViewType>
    KOKKOS_INLINE_FUNCTION void calc_crit_diagonal_length(int index, float xp, float yp, float zp, float cx, float cy,
                                                          float cz, int MyOrientation, ViewType GrainUnitVector) const {
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
            float x0 = xp + NeighborX[n] - cx;
            float y0 = yp + NeighborY[n] - cy;
            float z0 = zp + NeighborZ[n] - cz;
            float D0 = x0 * Fx[0] + y0 * Fy[0] + z0 * Fz[0];
            float D1 = x0 * Fx[1] + y0 * Fy[1] + z0 * Fz[1];
            float D2 = x0 * Fx[2] + y0 * Fy[2] + z0 * Fz[2];
            float D3 = x0 * Fx[3] + y0 * Fy[3] + z0 * Fz[3];
            float Dfabs = fmax(fmax(fabs(D0), fabs(D1)), fmax(fabs(D2), fabs(D3)));
            CritDiagonalLength(26 * index + n) = Dfabs;
        }
    }

    // Load data (GrainID, DOCenter, DiagonalLength) into ghost nodes if the given RankY is associated with a 1D halo
    // region Uses check to ensure that the buffer position does not reach the buffer size - if it does, return false
    // (otherwise return true) but keep incrementing the send size counters for use resizing the buffers in the future
    KOKKOS_INLINE_FUNCTION
    bool load_ghost_nodes(const int GhostGID, const float GhostDOCX, const float GhostDOCY, const float GhostDOCZ,
                          const float GhostDL, const int ny_local, const int coord_x, const int coord_y,
                          const int coord_z, const bool AtNorthBoundary, const bool AtSouthBoundary,
                          int NGrainOrientations) const {
        bool DataFitsInBuffer = true;
        if ((coord_y == 1) && (!(AtSouthBoundary))) {
            int GNPositionSouth = Kokkos::atomic_fetch_add(&SendSizeSouth(0), 1);
            if (GNPositionSouth >= BufSize)
                DataFitsInBuffer = false;
            else {
                BufferSouthSend(GNPositionSouth, 0) = static_cast<float>(coord_x);
                BufferSouthSend(GNPositionSouth, 1) = static_cast<float>(coord_z);
                BufferSouthSend(GNPositionSouth, 2) =
                    static_cast<float>(getGrainOrientation(GhostGID, NGrainOrientations, false));
                BufferSouthSend(GNPositionSouth, 3) = static_cast<float>(getGrainNumber(GhostGID, NGrainOrientations));
                BufferSouthSend(GNPositionSouth, 4) = GhostDOCX;
                BufferSouthSend(GNPositionSouth, 5) = GhostDOCY;
                BufferSouthSend(GNPositionSouth, 6) = GhostDOCZ;
                BufferSouthSend(GNPositionSouth, 7) = GhostDL;
            }
        }
        else if ((coord_y == ny_local - 2) && (!(AtNorthBoundary))) {
            int GNPositionNorth = Kokkos::atomic_fetch_add(&SendSizeNorth(0), 1);
            if (GNPositionNorth >= BufSize)
                DataFitsInBuffer = false;
            else {
                BufferNorthSend(GNPositionNorth, 0) = static_cast<float>(coord_x);
                BufferNorthSend(GNPositionNorth, 1) = static_cast<float>(coord_z);
                BufferNorthSend(GNPositionNorth, 2) =
                    static_cast<float>(getGrainOrientation(GhostGID, NGrainOrientations, false));
                BufferNorthSend(GNPositionNorth, 3) = static_cast<float>(getGrainNumber(GhostGID, NGrainOrientations));
                BufferNorthSend(GNPositionNorth, 4) = GhostDOCX;
                BufferNorthSend(GNPositionNorth, 5) = GhostDOCY;
                BufferNorthSend(GNPositionNorth, 6) = GhostDOCZ;
                BufferNorthSend(GNPositionNorth, 7) = GhostDL;
            }
        }
        return DataFitsInBuffer;
    }
};

#endif
