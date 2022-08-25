// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAghostnodes.hpp"
#include "CAfunctions.hpp"
#include "CAupdate.hpp"
#include "mpi.h"

#include <cmath>
#include <vector>

//*****************************************************************************/
// 2D domain decomposition: update ghost nodes with new cell data from Nucleation and CellCapture routines
void GhostNodes2D(int, int, int NeighborRank_North, int NeighborRank_South, int NeighborRank_East,
                  int NeighborRank_West, int NeighborRank_NorthEast, int NeighborRank_NorthWest,
                  int NeighborRank_SouthEast, int NeighborRank_SouthWest, int MyXSlices, int MyYSlices, int MyXOffset,
                  int MyYOffset, NList NeighborX, NList NeighborY, NList NeighborZ, ViewI CellType, ViewF DOCenter,
                  ViewI GrainID, ViewF GrainUnitVector, ViewF DiagonalLength, ViewF CritDiagonalLength,
                  int NGrainOrientations, Buffer2D BufferWestSend, Buffer2D BufferEastSend, Buffer2D BufferNorthSend,
                  Buffer2D BufferSouthSend, Buffer2D BufferNorthEastSend, Buffer2D BufferNorthWestSend,
                  Buffer2D BufferSouthEastSend, Buffer2D BufferSouthWestSend, Buffer2D BufferWestRecv,
                  Buffer2D BufferEastRecv, Buffer2D BufferNorthRecv, Buffer2D BufferSouthRecv,
                  Buffer2D BufferNorthEastRecv, Buffer2D BufferNorthWestRecv, Buffer2D BufferSouthEastRecv,
                  Buffer2D BufferSouthWestRecv, int BufSizeX, int BufSizeY, int BufSizeZ, int ZBound_Low) {

    Kokkos::fence();
    MPI_Barrier(MPI_COMM_WORLD);

    // Allocate requests.
    std::vector<MPI_Request> SendRequests(8, MPI_REQUEST_NULL);
    std::vector<MPI_Request> RecvRequests(8, MPI_REQUEST_NULL);

    // Send data to each other rank (MPI_Isend)
    MPI_Isend(BufferSouthSend.data(), 5 * BufSizeX * BufSizeZ, MPI_DOUBLE, NeighborRank_South, 0, MPI_COMM_WORLD,
              &SendRequests[0]);
    MPI_Isend(BufferNorthSend.data(), 5 * BufSizeX * BufSizeZ, MPI_DOUBLE, NeighborRank_North, 0, MPI_COMM_WORLD,
              &SendRequests[1]);
    MPI_Isend(BufferEastSend.data(), 5 * BufSizeY * BufSizeZ, MPI_DOUBLE, NeighborRank_East, 0, MPI_COMM_WORLD,
              &SendRequests[2]);
    MPI_Isend(BufferWestSend.data(), 5 * BufSizeY * BufSizeZ, MPI_DOUBLE, NeighborRank_West, 0, MPI_COMM_WORLD,
              &SendRequests[3]);
    MPI_Isend(BufferNorthWestSend.data(), 5 * BufSizeZ, MPI_DOUBLE, NeighborRank_NorthWest, 0, MPI_COMM_WORLD,
              &SendRequests[4]);
    MPI_Isend(BufferNorthEastSend.data(), 5 * BufSizeZ, MPI_DOUBLE, NeighborRank_NorthEast, 0, MPI_COMM_WORLD,
              &SendRequests[5]);
    MPI_Isend(BufferSouthWestSend.data(), 5 * BufSizeZ, MPI_DOUBLE, NeighborRank_SouthWest, 0, MPI_COMM_WORLD,
              &SendRequests[6]);
    MPI_Isend(BufferSouthEastSend.data(), 5 * BufSizeZ, MPI_DOUBLE, NeighborRank_SouthEast, 0, MPI_COMM_WORLD,
              &SendRequests[7]);

    // Receive buffers for all neighbors (MPI_Irecv)
    MPI_Irecv(BufferSouthRecv.data(), 5 * BufSizeX * BufSizeZ, MPI_DOUBLE, NeighborRank_South, 0, MPI_COMM_WORLD,
              &RecvRequests[0]);
    MPI_Irecv(BufferNorthRecv.data(), 5 * BufSizeX * BufSizeZ, MPI_DOUBLE, NeighborRank_North, 0, MPI_COMM_WORLD,
              &RecvRequests[1]);
    MPI_Irecv(BufferEastRecv.data(), 5 * BufSizeY * BufSizeZ, MPI_DOUBLE, NeighborRank_East, 0, MPI_COMM_WORLD,
              &RecvRequests[2]);
    MPI_Irecv(BufferWestRecv.data(), 5 * BufSizeY * BufSizeZ, MPI_DOUBLE, NeighborRank_West, 0, MPI_COMM_WORLD,
              &RecvRequests[3]);
    MPI_Irecv(BufferNorthWestRecv.data(), 5 * BufSizeZ, MPI_DOUBLE, NeighborRank_NorthWest, 0, MPI_COMM_WORLD,
              &RecvRequests[4]);
    MPI_Irecv(BufferNorthEastRecv.data(), 5 * BufSizeZ, MPI_DOUBLE, NeighborRank_NorthEast, 0, MPI_COMM_WORLD,
              &RecvRequests[5]);
    MPI_Irecv(BufferSouthWestRecv.data(), 5 * BufSizeZ, MPI_DOUBLE, NeighborRank_SouthWest, 0, MPI_COMM_WORLD,
              &RecvRequests[6]);
    MPI_Irecv(BufferSouthEastRecv.data(), 5 * BufSizeZ, MPI_DOUBLE, NeighborRank_SouthEast, 0, MPI_COMM_WORLD,
              &RecvRequests[7]);

    // unpack in any order
    bool unpack_complete = false;
    while (!unpack_complete) {
        // Get the next buffer to unpack from rank "unpack_index"
        int unpack_index = MPI_UNDEFINED;
        MPI_Waitany(8, RecvRequests.data(), &unpack_index, MPI_STATUS_IGNORE);

        // If there are no more buffers to unpack, leave the while loop
        if (MPI_UNDEFINED == unpack_index) {
            unpack_complete = true;
        }
        // Otherwise unpack the next buffer.
        else {
            int RecvBufSize;
            if (unpack_index <= 1) {
                RecvBufSize = BufSizeX * BufSizeZ;
            }
            else if (unpack_index >= 4) {
                RecvBufSize = BufSizeZ;
            }
            else {
                RecvBufSize = BufSizeY * BufSizeZ;
            }
            Kokkos::parallel_for(
                "BufferUnpack", RecvBufSize, KOKKOS_LAMBDA(const int &BufPosition) {
                    int RankX, RankY, RankZ, NewGrainID;
                    long int CellLocation;
                    bool Place = false;
                    double DOCenterX, DOCenterY, DOCenterZ, NewDiagonalLength;
                    // Which rank was the data received from?
                    if ((unpack_index == 0) && (NeighborRank_South != MPI_PROC_NULL)) {
                        // Data receieved from South
                        // Adjust X Position by +1 since X = 0 is not included in this buffer
                        RankX = BufPosition % BufSizeX + 1;
                        RankY = 0;
                        RankZ = BufPosition / BufSizeX;
                        CellLocation = RankZ * MyXSlices * MyYSlices + MyYSlices * RankX + RankY;
                        int GlobalCellLocation = CellLocation + ZBound_Low * MyXSlices * MyYSlices;
                        if ((BufferSouthRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferSouthRecv(BufPosition, 0));
                            DOCenterX = BufferSouthRecv(BufPosition, 1);
                            DOCenterY = BufferSouthRecv(BufPosition, 2);
                            DOCenterZ = BufferSouthRecv(BufPosition, 3);
                            NewDiagonalLength = BufferSouthRecv(BufPosition, 4);
                        }
                    }
                    else if ((unpack_index == 1) && (NeighborRank_North != MPI_PROC_NULL)) {
                        // Data received from North
                        // Adjust X Position by +1 since X = 0 is not included in this buffer
                        RankX = BufPosition % BufSizeX + 1;
                        RankY = MyYSlices - 1;
                        RankZ = BufPosition / BufSizeX;
                        CellLocation = RankZ * MyXSlices * MyYSlices + MyYSlices * RankX + RankY;
                        int GlobalCellLocation = CellLocation + ZBound_Low * MyXSlices * MyYSlices;
                        if ((BufferNorthRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferNorthRecv(BufPosition, 0));
                            DOCenterX = BufferNorthRecv(BufPosition, 1);
                            DOCenterY = BufferNorthRecv(BufPosition, 2);
                            DOCenterZ = BufferNorthRecv(BufPosition, 3);
                            NewDiagonalLength = BufferNorthRecv(BufPosition, 4);
                        }
                    }
                    else if ((unpack_index == 2) && (NeighborRank_East != MPI_PROC_NULL)) {
                        // Data received from East
                        // Adjust Y Position by +1 since Y = 0 is not included in this buffer
                        RankX = MyXSlices - 1;
                        RankY = BufPosition % BufSizeY + 1;
                        RankZ = BufPosition / BufSizeY;
                        CellLocation = RankZ * MyXSlices * MyYSlices + MyYSlices * RankX + RankY;
                        int GlobalCellLocation = CellLocation + ZBound_Low * MyXSlices * MyYSlices;
                        if ((BufferEastRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferEastRecv(BufPosition, 0));
                            DOCenterX = BufferEastRecv(BufPosition, 1);
                            DOCenterY = BufferEastRecv(BufPosition, 2);
                            DOCenterZ = BufferEastRecv(BufPosition, 3);
                            NewDiagonalLength = BufferEastRecv(BufPosition, 4);
                        }
                    }
                    else if ((unpack_index == 3) && (NeighborRank_West != MPI_PROC_NULL)) {
                        // Data received from West
                        // Adjust Y Position by +1 since Y = 0 is not included in this buffer
                        RankX = 0;
                        RankY = BufPosition % BufSizeY + 1;
                        RankZ = BufPosition / BufSizeY;
                        CellLocation = RankZ * MyXSlices * MyYSlices + MyYSlices * RankX + RankY;
                        int GlobalCellLocation = CellLocation + ZBound_Low * MyXSlices * MyYSlices;
                        if ((BufferWestRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferWestRecv(BufPosition, 0));
                            DOCenterX = BufferWestRecv(BufPosition, 1);
                            DOCenterY = BufferWestRecv(BufPosition, 2);
                            DOCenterZ = BufferWestRecv(BufPosition, 3);
                            NewDiagonalLength = BufferWestRecv(BufPosition, 4);
                        }
                    }
                    else if ((unpack_index == 4) && (NeighborRank_NorthWest != MPI_PROC_NULL)) {
                        // Data received from NorthWest
                        RankX = 0;
                        RankY = MyYSlices - 1;
                        RankZ = BufPosition;
                        CellLocation = RankZ * MyXSlices * MyYSlices + MyYSlices * RankX + RankY;
                        int GlobalCellLocation = CellLocation + ZBound_Low * MyXSlices * MyYSlices;
                        if ((BufferNorthWestRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferNorthWestRecv(BufPosition, 0));
                            DOCenterX = BufferNorthWestRecv(BufPosition, 1);
                            DOCenterY = BufferNorthWestRecv(BufPosition, 2);
                            DOCenterZ = BufferNorthWestRecv(BufPosition, 3);
                            NewDiagonalLength = BufferNorthWestRecv(BufPosition, 4);
                        }
                    }
                    else if ((unpack_index == 5) && (NeighborRank_NorthEast != MPI_PROC_NULL)) {
                        // Data received from NorthEast
                        RankX = MyXSlices - 1;
                        RankY = MyYSlices - 1;
                        RankZ = BufPosition;
                        CellLocation = RankZ * MyXSlices * MyYSlices + MyYSlices * RankX + RankY;
                        int GlobalCellLocation = CellLocation + ZBound_Low * MyXSlices * MyYSlices;
                        if ((BufferNorthEastRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferNorthEastRecv(BufPosition, 0));
                            DOCenterX = BufferNorthEastRecv(BufPosition, 1);
                            DOCenterY = BufferNorthEastRecv(BufPosition, 2);
                            DOCenterZ = BufferNorthEastRecv(BufPosition, 3);
                            NewDiagonalLength = BufferNorthEastRecv(BufPosition, 4);
                        }
                    }
                    else if ((unpack_index == 6) && (NeighborRank_SouthWest != MPI_PROC_NULL)) {
                        // Data received from SouthWest
                        RankX = 0;
                        RankY = 0;
                        RankZ = BufPosition;
                        CellLocation = RankZ * MyXSlices * MyYSlices + MyYSlices * RankX + RankY;
                        int GlobalCellLocation = CellLocation + ZBound_Low * MyXSlices * MyYSlices;
                        if ((BufferSouthWestRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferSouthWestRecv(BufPosition, 0));
                            DOCenterX = BufferSouthWestRecv(BufPosition, 1);
                            DOCenterY = BufferSouthWestRecv(BufPosition, 2);
                            DOCenterZ = BufferSouthWestRecv(BufPosition, 3);
                            NewDiagonalLength = BufferSouthWestRecv(BufPosition, 4);
                        }
                    }
                    else if ((unpack_index == 7) && (NeighborRank_SouthEast != MPI_PROC_NULL)) {
                        // Data received from SouthEast
                        RankX = MyXSlices - 1;
                        RankY = 0;
                        RankZ = BufPosition;
                        CellLocation = RankZ * MyXSlices * MyYSlices + MyYSlices * RankX + RankY;
                        int GlobalCellLocation = CellLocation + ZBound_Low * MyXSlices * MyYSlices;
                        if ((BufferSouthEastRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferSouthEastRecv(BufPosition, 0));
                            DOCenterX = BufferSouthEastRecv(BufPosition, 1);
                            DOCenterY = BufferSouthEastRecv(BufPosition, 2);
                            DOCenterZ = BufferSouthEastRecv(BufPosition, 3);
                            NewDiagonalLength = BufferSouthEastRecv(BufPosition, 4);
                        }
                    }

                    if (Place) {
                        // Update this ghost node cell's information with data from other rank
                        int GlobalZ = RankZ + ZBound_Low;
                        int GlobalCellLocation = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                        GrainID(GlobalCellLocation) = NewGrainID;
                        DOCenter((long int)(3) * CellLocation) = static_cast<float>(DOCenterX);
                        DOCenter((long int)(3) * CellLocation + (long int)(1)) = static_cast<float>(DOCenterY);
                        DOCenter((long int)(3) * CellLocation + (long int)(2)) = static_cast<float>(DOCenterZ);
                        int MyOrientation = getGrainOrientation(GrainID(GlobalCellLocation), NGrainOrientations);
                        DiagonalLength(CellLocation) = static_cast<float>(NewDiagonalLength);
                        // Global coordinates of cell center
                        double xp = RankX + MyXOffset + 0.5;
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
    // Wait on send requests
    MPI_Waitall(8, SendRequests.data(), MPI_STATUSES_IGNORE);
    Kokkos::fence();
}

//*****************************************************************************/
// 1D domain decomposition: update ghost nodes with new cell data from Nucleation and CellCapture routines
void GhostNodes1D(int, int, int NeighborRank_North, int NeighborRank_South, int MyXSlices, int MyYSlices, int MyXOffset,
                  int MyYOffset, NList NeighborX, NList NeighborY, NList NeighborZ, ViewI CellType, ViewF DOCenter,
                  ViewI GrainID, ViewF GrainUnitVector, ViewF DiagonalLength, ViewF CritDiagonalLength,
                  int NGrainOrientations, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, Buffer2D BufferNorthRecv,
                  Buffer2D BufferSouthRecv, int BufSizeX, int, int BufSizeZ, int ZBound_Low) {

    std::vector<MPI_Request> SendRequests(2, MPI_REQUEST_NULL);
    std::vector<MPI_Request> RecvRequests(2, MPI_REQUEST_NULL);

    // Send data to each other rank (MPI_Isend)
    MPI_Isend(BufferSouthSend.data(), 5 * BufSizeX * BufSizeZ, MPI_DOUBLE, NeighborRank_South, 0, MPI_COMM_WORLD,
              &SendRequests[0]);
    MPI_Isend(BufferNorthSend.data(), 5 * BufSizeX * BufSizeZ, MPI_DOUBLE, NeighborRank_North, 0, MPI_COMM_WORLD,
              &SendRequests[1]);

    // Receive buffers for all neighbors (MPI_Irecv)
    MPI_Irecv(BufferSouthRecv.data(), 5 * BufSizeX * BufSizeZ, MPI_DOUBLE, NeighborRank_South, 0, MPI_COMM_WORLD,
              &RecvRequests[0]);
    MPI_Irecv(BufferNorthRecv.data(), 5 * BufSizeX * BufSizeZ, MPI_DOUBLE, NeighborRank_North, 0, MPI_COMM_WORLD,
              &RecvRequests[1]);

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
            int RecvBufSize = BufSizeX * BufSizeZ;
            Kokkos::parallel_for(
                "BufferUnpack", RecvBufSize, KOKKOS_LAMBDA(const int &BufPosition) {
                    int RankX, RankY, RankZ, NewGrainID;
                    long int CellLocation;
                    double DOCenterX, DOCenterY, DOCenterZ, NewDiagonalLength;
                    bool Place = false;
                    RankZ = BufPosition / BufSizeX;
                    RankX = BufPosition % BufSizeX;
                    // Which rank was the data received from?
                    if ((unpack_index == 0) && (NeighborRank_South != MPI_PROC_NULL)) {
                        // Data receieved from South
                        RankY = 0;
                        CellLocation = RankZ * MyXSlices * MyYSlices + MyYSlices * RankX + RankY;
                        int GlobalCellLocation = CellLocation + ZBound_Low * MyXSlices * MyYSlices;
                        if ((BufferSouthRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferSouthRecv(BufPosition, 0));
                            DOCenterX = BufferSouthRecv(BufPosition, 1);
                            DOCenterY = BufferSouthRecv(BufPosition, 2);
                            DOCenterZ = BufferSouthRecv(BufPosition, 3);
                            NewDiagonalLength = BufferSouthRecv(BufPosition, 4);
                        }
                    }
                    else if ((unpack_index == 1) && (NeighborRank_North != MPI_PROC_NULL)) {
                        // Data received from North
                        RankY = MyYSlices - 1;
                        CellLocation = RankZ * MyXSlices * MyYSlices + MyYSlices * RankX + RankY;
                        int GlobalCellLocation = CellLocation + ZBound_Low * MyXSlices * MyYSlices;
                        if ((BufferNorthRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferNorthRecv(BufPosition, 0));
                            DOCenterX = BufferNorthRecv(BufPosition, 1);
                            DOCenterY = BufferNorthRecv(BufPosition, 2);
                            DOCenterZ = BufferNorthRecv(BufPosition, 3);
                            NewDiagonalLength = BufferNorthRecv(BufPosition, 4);
                        }
                    }
                    if (Place) {
                        int GlobalZ = RankZ + ZBound_Low;
                        int GlobalCellLocation = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;

                        // Update this ghost node cell's information with data from other rank
                        GrainID(GlobalCellLocation) = NewGrainID;
                        DOCenter((long int)(3) * CellLocation) = static_cast<float>(DOCenterX);
                        DOCenter((long int)(3) * CellLocation + (long int)(1)) = static_cast<float>(DOCenterY);
                        DOCenter((long int)(3) * CellLocation + (long int)(2)) = static_cast<float>(DOCenterZ);
                        int MyOrientation = getGrainOrientation(GrainID(GlobalCellLocation), NGrainOrientations);
                        DiagonalLength(CellLocation) = static_cast<float>(NewDiagonalLength);
                        // Global coordinates of cell center
                        double xp = RankX + MyXOffset + 0.5;
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

    // Wait on send requests
    MPI_Waitall(2, SendRequests.data(), MPI_STATUSES_IGNORE);
    Kokkos::fence();
}
