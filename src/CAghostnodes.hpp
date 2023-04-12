// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_GHOST_HPP
#define EXACA_GHOST_HPP

#include "CAfunctions.hpp"
#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

// Load data (GrainID, DOCenter, DiagonalLength) into ghost nodes if the given RankY is associated with a 1D halo region
KOKKOS_INLINE_FUNCTION void loadghostnodes(const int GhostGID, const float GhostDOCX, const float GhostDOCY,
                                           const float GhostDOCZ, const float GhostDL, const int BufSizeX,
                                           const int MyYSlices, const int RankX, const int RankY, const int RankZ,
                                           const bool AtNorthBoundary, const bool AtSouthBoundary,
                                           Buffer2D BufferSouthSend, Buffer2D BufferNorthSend, int NGrainOrientations) {

    if ((RankY == 1) && (!(AtSouthBoundary))) {
        int GNPosition = RankZ * BufSizeX + RankX;
        int MyGrainOrientation = getGrainOrientation(GhostGID, NGrainOrientations, false);
        BufferSouthSend(GNPosition, 0) = static_cast<float>(MyGrainOrientation);
        BufferSouthSend(GNPosition, 1) = getGrainNumber(GhostGID, NGrainOrientations);
        BufferSouthSend(GNPosition, 2) = GhostDOCX;
        BufferSouthSend(GNPosition, 3) = GhostDOCY;
        BufferSouthSend(GNPosition, 4) = GhostDOCZ;
        BufferSouthSend(GNPosition, 5) = GhostDL;
    }
    else if ((RankY == MyYSlices - 2) && (!(AtNorthBoundary))) {
        int GNPosition = RankZ * BufSizeX + RankX;
        int MyGrainOrientation = getGrainOrientation(GhostGID, NGrainOrientations, false);
        BufferNorthSend(GNPosition, 0) = static_cast<float>(MyGrainOrientation);
        BufferNorthSend(GNPosition, 1) = getGrainNumber(GhostGID, NGrainOrientations);
        BufferNorthSend(GNPosition, 2) = GhostDOCX;
        BufferNorthSend(GNPosition, 3) = GhostDOCY;
        BufferNorthSend(GNPosition, 4) = GhostDOCZ;
        BufferNorthSend(GNPosition, 5) = GhostDL;
    }
}
void GhostNodes1D(int, int, int NeighborRank_North, int NeighborRank_South, int nx, int MyYSlices, int MyYOffset,
                  NList NeighborX, NList NeighborY, NList NeighborZ, ViewI CellType, ViewF DOCenter, ViewI GrainID,
                  ViewF GrainUnitVector, ViewF DiagonalLength, ViewF CritDiagonalLength, int NGrainOrientations,
                  Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, Buffer2D BufferNorthRecv,
                  Buffer2D BufferSouthRecv, int BufSizeX, int BufSizeZ, int ZBound_Low);

#endif
