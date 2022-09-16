// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
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
// 1D domain decomposition: update ghost nodes with new cell data from Nucleation and CellCapture routines
void GhostNodes1D(int, int, int NeighborRank_North, int NeighborRank_South, int nx, int MyYSlices, int MyYOffset,
                  NList NeighborX, NList NeighborY, NList NeighborZ, ViewI CellType, ViewF DOCenter, ViewI GrainID,
                  ViewF GrainUnitVector, ViewF DiagonalLength, ViewF CritDiagonalLength, int NGrainOrientations,
                  Halo &halo, int BufSizeX, int BufSizeZ, int ZBound_Low) {

    std::vector<MPI_Request> SendRequests(2, MPI_REQUEST_NULL);
    std::vector<MPI_Request> RecvRequests(2, MPI_REQUEST_NULL);

    // Send data to each other rank (MPI_Isend)
    MPI_Isend(halo.BufferSouthSend.data(), BufSizeX * BufSizeZ * sizeof(CommData), MPI_BYTE, NeighborRank_South, 0,
              MPI_COMM_WORLD, &SendRequests[0]);
    MPI_Isend(halo.BufferNorthSend.data(), BufSizeX * BufSizeZ * sizeof(CommData), MPI_BYTE, NeighborRank_North, 0,
              MPI_COMM_WORLD, &SendRequests[1]);

    // Receive buffers for all neighbors (MPI_Irecv)
    MPI_Irecv(halo.BufferSouthRecv.data(), BufSizeX * BufSizeZ * sizeof(CommData), MPI_BYTE, NeighborRank_South, 0,
              MPI_COMM_WORLD, &RecvRequests[0]);
    MPI_Irecv(halo.BufferNorthRecv.data(), BufSizeX * BufSizeZ * sizeof(CommData), MPI_BYTE, NeighborRank_North, 0,
              MPI_COMM_WORLD, &RecvRequests[1]);

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
                        CellLocation = RankZ * nx * MyYSlices + MyYSlices * RankX + RankY;
                        int GlobalCellLocation = CellLocation + ZBound_Low * nx * MyYSlices;
                        if ((halo.BufferSouthRecv(BufPosition).DiagonalLength > 0) &&
                            (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = halo.BufferSouthRecv(BufPosition).GrainID;
                            DOCenterX = halo.BufferSouthRecv(BufPosition).DOCenterX;
                            DOCenterY = halo.BufferSouthRecv(BufPosition).DOCenterY;
                            DOCenterZ = halo.BufferSouthRecv(BufPosition).DOCenterZ;
                            NewDiagonalLength = halo.BufferSouthRecv(BufPosition).DiagonalLength;
                        }
                    }
                    else if ((unpack_index == 1) && (NeighborRank_North != MPI_PROC_NULL)) {
                        // Data received from North
                        RankY = MyYSlices - 1;
                        CellLocation = RankZ * nx * MyYSlices + MyYSlices * RankX + RankY;
                        int GlobalCellLocation = CellLocation + ZBound_Low * nx * MyYSlices;
                        if ((halo.BufferNorthRecv(BufPosition).DiagonalLength > 0) &&
                            (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = halo.BufferNorthRecv(BufPosition).GrainID;
                            DOCenterX = halo.BufferNorthRecv(BufPosition).DOCenterX;
                            DOCenterY = halo.BufferNorthRecv(BufPosition).DOCenterY;
                            DOCenterZ = halo.BufferNorthRecv(BufPosition).DOCenterZ;
                            NewDiagonalLength = halo.BufferNorthRecv(BufPosition).DiagonalLength;
                        }
                    }
                    if (Place) {
                        int GlobalZ = RankZ + ZBound_Low;
                        int GlobalCellLocation = GlobalZ * nx * MyYSlices + RankX * MyYSlices + RankY;

                        // Update this ghost node cell's information with data from other rank
                        GrainID(GlobalCellLocation) = NewGrainID;
                        DOCenter((long int)(3) * CellLocation) = static_cast<float>(DOCenterX);
                        DOCenter((long int)(3) * CellLocation + (long int)(1)) = static_cast<float>(DOCenterY);
                        DOCenter((long int)(3) * CellLocation + (long int)(2)) = static_cast<float>(DOCenterZ);
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

    // Wait on send requests
    MPI_Waitall(2, SendRequests.data(), MPI_STATUSES_IGNORE);
    Kokkos::fence();
}
