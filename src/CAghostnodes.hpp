// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_GHOST_HPP
#define EXACA_GHOST_HPP

#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

// Load data (GrainID, DOCenter, DiagonalLength) into ghost nodes if the given RankY is associated with a 1D halo region
KOKKOS_INLINE_FUNCTION void loadghostnodes(const double GhostGID, const double GhostDOCX, const double GhostDOCY,
                                           const double GhostDOCZ, const double GhostDL, const int BufSizeX,
                                           const int MyYSlices, const int RankX, const int RankY, const int RankZ,
                                           const bool AtNorthBoundary, const bool AtSouthBoundary,
                                           Buffer2D BufferSouthSend, Buffer2D BufferNorthSend) {

    if ((RankY == 1) && (!(AtSouthBoundary))) {
        int GNPosition = RankZ * BufSizeX + RankX;
        BufferSouthSend(GNPosition, 0) = GhostGID;
        BufferSouthSend(GNPosition, 1) = GhostDOCX;
        BufferSouthSend(GNPosition, 2) = GhostDOCY;
        BufferSouthSend(GNPosition, 3) = GhostDOCZ;
        BufferSouthSend(GNPosition, 4) = GhostDL;
    }
    else if ((RankY == MyYSlices - 2) && (!(AtNorthBoundary))) {
        int GNPosition = RankZ * BufSizeX + RankX;
        BufferNorthSend(GNPosition, 0) = GhostGID;
        BufferNorthSend(GNPosition, 1) = GhostDOCX;
        BufferNorthSend(GNPosition, 2) = GhostDOCY;
        BufferNorthSend(GNPosition, 3) = GhostDOCZ;
        BufferNorthSend(GNPosition, 4) = GhostDL;
    }
}

// Load data (GrainID, DOCenter, DiagonalLength) into ghost nodes if the given RankX and RankY is associated with a 2D
// halo region
KOKKOS_INLINE_FUNCTION void
loadghostnodes(const double GhostGID, const double GhostDOCX, const double GhostDOCY, const double GhostDOCZ,
               const double GhostDL, const int BufSizeX, const int BufSizeY, const int MyXSlices, const int MyYSlices,
               const int RankX, const int RankY, const int RankZ, const bool AtNorthBoundary,
               const bool AtSouthBoundary, const bool AtWestBoundary, const bool AtEastBoundary,
               Buffer2D BufferSouthSend, Buffer2D BufferNorthSend, Buffer2D BufferWestSend, Buffer2D BufferEastSend,
               Buffer2D BufferNorthEastSend, Buffer2D BufferSouthEastSend, Buffer2D BufferSouthWestSend,
               Buffer2D BufferNorthWestSend) {

    if ((RankY == 1) && (RankX > 1) && (RankX < MyXSlices - 2) && (!(AtSouthBoundary))) {
        int GNPosition = RankZ * BufSizeX + RankX - 1;
        BufferSouthSend(GNPosition, 0) = GhostGID;
        BufferSouthSend(GNPosition, 1) = GhostDOCX;
        BufferSouthSend(GNPosition, 2) = GhostDOCY;
        BufferSouthSend(GNPosition, 3) = GhostDOCZ;
        BufferSouthSend(GNPosition, 4) = GhostDL;
    }
    if ((RankY == MyYSlices - 2) && (RankX > 1) && (RankX < MyXSlices - 2) && (!(AtNorthBoundary))) {
        int GNPosition = RankZ * BufSizeX + RankX - 1;
        BufferNorthSend(GNPosition, 0) = GhostGID;
        BufferNorthSend(GNPosition, 1) = GhostDOCX;
        BufferNorthSend(GNPosition, 2) = GhostDOCY;
        BufferNorthSend(GNPosition, 3) = GhostDOCZ;
        BufferNorthSend(GNPosition, 4) = GhostDL;
    }
    if ((RankX == 1) && (RankY > 1) && (RankY < MyYSlices - 2) && (!(AtWestBoundary))) {
        int GNPosition = RankZ * BufSizeY + RankY - 1;
        BufferWestSend(GNPosition, 0) = GhostGID;
        BufferWestSend(GNPosition, 1) = GhostDOCX;
        BufferWestSend(GNPosition, 2) = GhostDOCY;
        BufferWestSend(GNPosition, 3) = GhostDOCZ;
        BufferWestSend(GNPosition, 4) = GhostDL;
    }
    if ((RankX == MyXSlices - 2) && (RankY > 1) && (RankY < MyYSlices - 2) && (!(AtEastBoundary))) {
        int GNPosition = RankZ * BufSizeY + RankY - 1;
        BufferEastSend(GNPosition, 0) = GhostGID;
        BufferEastSend(GNPosition, 1) = GhostDOCX;
        BufferEastSend(GNPosition, 2) = GhostDOCY;
        BufferEastSend(GNPosition, 3) = GhostDOCZ;
        BufferEastSend(GNPosition, 4) = GhostDL;
    }
    if ((RankY == 1) && (RankX == MyXSlices - 2) && (!(AtSouthBoundary)) && (!(AtEastBoundary))) {
        int GNPosition = RankZ;
        BufferSouthEastSend(GNPosition, 0) = GhostGID;
        BufferSouthEastSend(GNPosition, 1) = GhostDOCX;
        BufferSouthEastSend(GNPosition, 2) = GhostDOCY;
        BufferSouthEastSend(GNPosition, 3) = GhostDOCZ;
        BufferSouthEastSend(GNPosition, 4) = GhostDL;
    }
    if ((RankX == 1) && (RankY == 1) && (!(AtSouthBoundary)) && (!(AtWestBoundary))) {
        int GNPosition = RankZ;
        BufferSouthWestSend(GNPosition, 0) = GhostGID;
        BufferSouthWestSend(GNPosition, 1) = GhostDOCX;
        BufferSouthWestSend(GNPosition, 2) = GhostDOCY;
        BufferSouthWestSend(GNPosition, 3) = GhostDOCZ;
        BufferSouthWestSend(GNPosition, 4) = GhostDL;
    }
    if ((RankX == MyXSlices - 2) && (RankY == MyYSlices - 2) && (!(AtNorthBoundary)) && (!(AtEastBoundary))) {
        int GNPosition = RankZ;
        BufferNorthEastSend(GNPosition, 0) = GhostGID;
        BufferNorthEastSend(GNPosition, 1) = GhostDOCX;
        BufferNorthEastSend(GNPosition, 2) = GhostDOCY;
        BufferNorthEastSend(GNPosition, 3) = GhostDOCZ;
        BufferNorthEastSend(GNPosition, 4) = GhostDL;
    }
    if ((RankX == MyXSlices - 2) && (RankY == 1) && (!(AtNorthBoundary)) && (!(AtWestBoundary))) {
        int GNPosition = RankZ;
        BufferNorthWestSend(GNPosition, 0) = GhostGID;
        BufferNorthWestSend(GNPosition, 1) = GhostDOCX;
        BufferNorthWestSend(GNPosition, 2) = GhostDOCY;
        BufferNorthWestSend(GNPosition, 3) = GhostDOCZ;
        BufferNorthWestSend(GNPosition, 4) = GhostDL;
    }
}
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
                  Buffer2D BufferSouthWestRecv, int BufSizeX, int BufSizeY, int BufSizeZ, int ZBound_Low);
void GhostNodes1D(int, int, int NeighborRank_North, int NeighborRank_South, int MyXSlices, int MyYSlices, int MyXOffset,
                  int MyYOffset, NList NeighborX, NList NeighborY, NList NeighborZ, ViewI CellType, ViewF DOCenter,
                  ViewI GrainID, ViewF GrainUnitVector, ViewF DiagonalLength, ViewF CritDiagonalLength,
                  int NGrainOrientations, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, Buffer2D BufferNorthRecv,
                  Buffer2D BufferSouthRecv, int BufSizeX, int, int BufSizeZ, int ZBound_Low);

#endif
