// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAghostnodes.hpp"
#include "CAfunctions.hpp"
#include "CAupdate.hpp"
#include "mpi.h"

#include <algorithm>
#include <cmath>
#include <vector>

// Set first index in send buffers to -1 (placeholder) for all cells in the buffer, and reset the counts of number of
// cells contained in buffers to 0s
void ResetSendBuffers(int BufSize, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, ViewI SendSizeNorth,
                      ViewI SendSizeSouth) {

    Kokkos::parallel_for(
        "BufferReset", BufSize, KOKKOS_LAMBDA(const int &i) {
            BufferNorthSend(i, 0) = -1.0;
            BufferSouthSend(i, 0) = -1.0;
        });
    Kokkos::parallel_for(
        "HaloCountReset", 1, KOKKOS_LAMBDA(const int) {
            SendSizeNorth(0) = 0;
            SendSizeSouth(0) = 0;
        });
}

// Count the number of cells' in halo regions where the data did not fit into the send buffers, and resize the buffers
// if necessary, returning the new buffer size
int ResizeBuffers(Buffer2D &BufferNorthSend, Buffer2D &BufferSouthSend, Buffer2D &BufferNorthRecv,
                  Buffer2D &BufferSouthRecv, ViewI SendSizeNorth, ViewI SendSizeSouth, ViewI_H SendSizeNorth_Host,
                  ViewI_H SendSizeSouth_Host, int OldBufSize, int NumCellsBufferPadding) {

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
        Kokkos::resize(BufferNorthSend, NewBufSize, 8);
        Kokkos::resize(BufferSouthSend, NewBufSize, 8);
        Kokkos::resize(BufferNorthRecv, NewBufSize, 8);
        Kokkos::resize(BufferSouthRecv, NewBufSize, 8);
        // Reset count variables on device to the old buffer size
        Kokkos::parallel_for(
            "ResetCounts", 1, KOKKOS_LAMBDA(const int &) {
                SendSizeNorth(0) = OldBufSize;
                SendSizeSouth(0) = OldBufSize;
            });
    }
    else
        NewBufSize = OldBufSize;
    return NewBufSize;
}

// Reset the buffer sizes to a set value (defaulting to 25, which was the initial size) preserving the existing values
void ResetBufferCapacity(Buffer2D &BufferNorthSend, Buffer2D &BufferSouthSend, Buffer2D &BufferNorthRecv,
                         Buffer2D &BufferSouthRecv, int NewBufSize) {
    Kokkos::resize(BufferNorthSend, NewBufSize, 8);
    Kokkos::resize(BufferSouthSend, NewBufSize, 8);
    Kokkos::resize(BufferNorthRecv, NewBufSize, 8);
    Kokkos::resize(BufferSouthRecv, NewBufSize, 8);
}

// Refill the buffers as necessary starting from the old count size, using the data from cells marked with type
// ActiveFailedBufferLoad
void RefillBuffers(int nx, int nzActive, int MyYSlices, int ZBound_Low, ViewI CellType, Buffer2D BufferNorthSend,
                   Buffer2D BufferSouthSend, ViewI SendSizeNorth, ViewI SendSizeSouth, bool AtNorthBoundary,
                   bool AtSouthBoundary, ViewI GrainID, ViewF DOCenter, ViewF DiagonalLength, int NGrainOrientations,
                   int BufSize) {

    Kokkos::parallel_for(
        "FillSendBuffersOverflow", nx, KOKKOS_LAMBDA(const int &i) {
            for (int k = 0; k < nzActive; k++) {
                int GlobalCellCoordinateSouth = (k + ZBound_Low) * nx * MyYSlices + i * MyYSlices + 1;
                int GlobalCellCoordinateNorth = (k + ZBound_Low) * nx * MyYSlices + i * MyYSlices + MyYSlices - 2;
                if (CellType(GlobalCellCoordinateSouth) == ActiveFailedBufferLoad) {
                    int ActiveLayerCoordinateSouth = k * nx * MyYSlices + i * MyYSlices + 1;
                    int GhostGID = GrainID(GlobalCellCoordinateSouth);
                    float GhostDOCX = DOCenter(3 * ActiveLayerCoordinateSouth);
                    float GhostDOCY = DOCenter(3 * ActiveLayerCoordinateSouth + 1);
                    float GhostDOCZ = DOCenter(3 * ActiveLayerCoordinateSouth + 2);
                    float GhostDL = DiagonalLength(ActiveLayerCoordinateSouth);
                    // Collect data for the ghost nodes, if necessary
                    // Data loaded into the ghost nodes is for the cell that was just captured
                    bool DataFitsInBuffer =
                        loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth, SendSizeSouth,
                                       MyYSlices, i, 1, k, AtNorthBoundary, AtSouthBoundary, BufferSouthSend,
                                       BufferNorthSend, NGrainOrientations, BufSize);
                    CellType(GlobalCellCoordinateSouth) = Active;
                    // If data doesn't fit in the buffer after the resize, warn that buffer data may have been lost
                    if (!(DataFitsInBuffer))
                        printf("Warning: Send/recv buffer resize failed to include all necessary data, predicted "
                               "results at MPI processor boundaries may be inaccurate\n");
                }
                else if (CellType(GlobalCellCoordinateSouth) == LiquidFailedBufferLoad) {
                    bool DataFitsInBuffer =
                        loadghostnodes_liquid(SendSizeNorth, SendSizeSouth, MyYSlices, i, 1, k, AtNorthBoundary,
                                              AtSouthBoundary, BufferSouthSend, BufferNorthSend, BufSize);
                    CellType(GlobalCellCoordinateSouth) = Liquid;
                    // If data doesn't fit in the buffer after the resize, warn that buffer data may have been lost
                    if (!(DataFitsInBuffer))
                        printf("Warning: Send/recv buffer resize failed to include all necessary data, predicted "
                               "results at MPI processor boundaries may be inaccurate\n");
                }
                if (CellType(GlobalCellCoordinateNorth) == ActiveFailedBufferLoad) {
                    int ActiveLayerCoordinateNorth = k * nx * MyYSlices + i * MyYSlices + MyYSlices - 2;
                    int GhostGID = GrainID(GlobalCellCoordinateNorth);
                    float GhostDOCX = DOCenter(3 * ActiveLayerCoordinateNorth);
                    float GhostDOCY = DOCenter(3 * ActiveLayerCoordinateNorth + 1);
                    float GhostDOCZ = DOCenter(3 * ActiveLayerCoordinateNorth + 2);
                    float GhostDL = DiagonalLength(ActiveLayerCoordinateNorth);
                    // Collect data for the ghost nodes, if necessary
                    // Data loaded into the ghost nodes is for the cell that was just captured
                    bool DataFitsInBuffer =
                        loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth, SendSizeSouth,
                                       MyYSlices, i, MyYSlices - 2, k, AtNorthBoundary, AtSouthBoundary,
                                       BufferSouthSend, BufferNorthSend, NGrainOrientations, BufSize);
                    CellType(GlobalCellCoordinateNorth) = Active;
                    // If data doesn't fit in the buffer after the resize, warn that buffer data may have been lost
                    if (!(DataFitsInBuffer))
                        printf("Warning: Send/recv buffer resize failed to include all necessary data, predicted "
                               "results at MPI processor boundaries may be inaccurate\n");
                }
                else if (CellType(GlobalCellCoordinateNorth) == LiquidFailedBufferLoad) {
                    bool DataFitsInBuffer = loadghostnodes_liquid(SendSizeNorth, SendSizeSouth, MyYSlices, i,
                                                                  MyYSlices - 2, k, AtNorthBoundary, AtSouthBoundary,
                                                                  BufferSouthSend, BufferNorthSend, BufSize);
                    CellType(GlobalCellCoordinateNorth) = Liquid;
                    // If data doesn't fit in the buffer after the resize, warn that buffer data may have been lost
                    if (!(DataFitsInBuffer))
                        printf("Warning: Send/recv buffer resize failed to include all necessary data, predicted "
                               "results at MPI processor boundaries may be inaccurate\n");
                }
            }
        });
    Kokkos::fence();
}

//*****************************************************************************/
// 1D domain decomposition: update ghost nodes with new cell data from Nucleation and CellCapture routines
void GhostNodes1D(int, int, int NeighborRank_North, int NeighborRank_South, int nx, int MyYSlices, int MyYOffset,
                  NList NeighborX, NList NeighborY, NList NeighborZ, ViewI CellType, ViewF DOCenter, ViewI GrainID,
                  ViewF GrainUnitVector, ViewF DiagonalLength, ViewF CritDiagonalLength, int NGrainOrientations,
                  Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, Buffer2D BufferNorthRecv,
                  Buffer2D BufferSouthRecv, int BufSize, int ZBound_Low, ViewI SendSizeNorth, ViewI SendSizeSouth) {

    std::vector<MPI_Request> SendRequests(2, MPI_REQUEST_NULL);
    std::vector<MPI_Request> RecvRequests(2, MPI_REQUEST_NULL);

    // Send data to each other rank (MPI_Isend)
    MPI_Isend(BufferSouthSend.data(), 8 * BufSize, MPI_FLOAT, NeighborRank_South, 0, MPI_COMM_WORLD, &SendRequests[0]);
    MPI_Isend(BufferNorthSend.data(), 8 * BufSize, MPI_FLOAT, NeighborRank_North, 0, MPI_COMM_WORLD, &SendRequests[1]);

    // Receive buffers for all neighbors (MPI_Irecv)
    MPI_Irecv(BufferSouthRecv.data(), 8 * BufSize, MPI_FLOAT, NeighborRank_South, 0, MPI_COMM_WORLD, &RecvRequests[0]);
    MPI_Irecv(BufferNorthRecv.data(), 8 * BufSize, MPI_FLOAT, NeighborRank_North, 0, MPI_COMM_WORLD, &RecvRequests[1]);

    // unpack in any order
    bool unpack_complete = false;
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
                "BufferUnpack", BufSize, KOKKOS_LAMBDA(const int &BufPosition) {
                    int RankX, RankY, RankZ, NewGrainID;
                    long int CellLocation;
                    float DOCenterX, DOCenterY, DOCenterZ, NewDiagonalLength;
                    bool Place = false;
                    // Which rank was the data received from? Is there valid data at this position in the buffer (i.e.,
                    // not set to -1.0)?
                    if ((unpack_index == 0) && (BufferSouthRecv(BufPosition, 0) != -1.0) &&
                        (NeighborRank_South != MPI_PROC_NULL)) {
                        // Data receieved from South
                        RankX = static_cast<int>(BufferSouthRecv(BufPosition, 0));
                        RankY = 0;
                        RankZ = static_cast<int>(BufferSouthRecv(BufPosition, 1));
                        CellLocation = RankZ * nx * MyYSlices + MyYSlices * RankX + RankY;
                        int GlobalCellLocation = CellLocation + ZBound_Low * nx * MyYSlices;
                        // Two possibilities: buffer data with non-zero diagonal length was loaded, and a liquid cell
                        // may have to be updated to active - or zero diagonal length data was loaded, and an active
                        // cell may have to be updated to liquid
                        if (CellType(GlobalCellLocation) == Liquid) {
                            Place = true;
                            int MyGrainOrientation = static_cast<int>(BufferSouthRecv(BufPosition, 2));
                            int MyGrainNumber = static_cast<int>(BufferSouthRecv(BufPosition, 3));
                            NewGrainID = getGrainID(NGrainOrientations, MyGrainOrientation, MyGrainNumber);
                            DOCenterX = BufferSouthRecv(BufPosition, 4);
                            DOCenterY = BufferSouthRecv(BufPosition, 5);
                            DOCenterZ = BufferSouthRecv(BufPosition, 6);
                            NewDiagonalLength = BufferSouthRecv(BufPosition, 7);
                        }
                        else if ((CellType(GlobalCellLocation) == Active) && (BufferSouthRecv(BufPosition, 7) == 0.0)) {
                            CellType(GlobalCellLocation) = Liquid;
                        }
                    }
                    else if ((unpack_index == 1) && (BufferNorthRecv(BufPosition, 0) != -1.0) &&
                             (NeighborRank_North != MPI_PROC_NULL)) {
                        // Data received from North
                        RankX = static_cast<int>(BufferNorthRecv(BufPosition, 0));
                        RankY = MyYSlices - 1;
                        RankZ = static_cast<int>(BufferNorthRecv(BufPosition, 1));
                        CellLocation = RankZ * nx * MyYSlices + MyYSlices * RankX + RankY;
                        int GlobalCellLocation = CellLocation + ZBound_Low * nx * MyYSlices;
                        // Two possibilities: buffer data with non-zero diagonal length was loaded, and a liquid cell
                        // may have to be updated to active - or zero diagonal length data was loaded, and an active
                        // cell may have to be updated to liquid
                        if (CellType(GlobalCellLocation) == Liquid) {
                            Place = true;
                            int MyGrainOrientation = static_cast<int>(BufferNorthRecv(BufPosition, 2));
                            int MyGrainNumber = static_cast<int>(BufferNorthRecv(BufPosition, 3));
                            NewGrainID = getGrainID(NGrainOrientations, MyGrainOrientation, MyGrainNumber);
                            DOCenterX = BufferNorthRecv(BufPosition, 4);
                            DOCenterY = BufferNorthRecv(BufPosition, 5);
                            DOCenterZ = BufferNorthRecv(BufPosition, 6);
                            NewDiagonalLength = BufferNorthRecv(BufPosition, 7);
                        }
                        else if ((CellType(GlobalCellLocation) == Active) && (BufferNorthRecv(BufPosition, 7) == 0.0)) {
                            CellType(GlobalCellLocation) = Liquid;
                        }
                    }
                    if (Place) {
                        int GlobalZ = RankZ + ZBound_Low;
                        int GlobalCellLocation = GlobalZ * nx * MyYSlices + RankX * MyYSlices + RankY;
                        // Update this ghost node cell's information with data from other rank
                        GrainID(GlobalCellLocation) = NewGrainID;
                        DOCenter((long int)(3) * CellLocation) = DOCenterX;
                        DOCenter((long int)(3) * CellLocation + (long int)(1)) = DOCenterY;
                        DOCenter((long int)(3) * CellLocation + (long int)(2)) = DOCenterZ;
                        int MyOrientation = getGrainOrientation(GrainID(GlobalCellLocation), NGrainOrientations);
                        DiagonalLength(CellLocation) = static_cast<float>(NewDiagonalLength);
                        // Global coordinates of cell center
                        double xp = RankX + 0.5;
                        double yp = RankY + MyYOffset + 0.5;
                        double zp = GlobalZ + 0.5;
                        // Calculate critical values at which this active cell leads to the activation of a neighboring
                        // liquid cell
                        calcCritDiagonalLength(CellLocation, xp, yp, zp, DOCenterX, DOCenterY, DOCenterZ, NeighborX,
                                               NeighborY, NeighborZ, MyOrientation, GrainUnitVector,
                                               CritDiagonalLength);
                        CellType(GlobalCellLocation) = Active;
                    }
                });
        }
    }

    // Reset send buffer data to -1 (used as placeholder) and reset the number of cells stored in the buffers to 0
    ResetSendBuffers(BufSize, BufferNorthSend, BufferSouthSend, SendSizeNorth, SendSizeSouth);
    // Wait on send requests
    MPI_Waitall(2, SendRequests.data(), MPI_STATUSES_IGNORE);
    Kokkos::fence();
}
