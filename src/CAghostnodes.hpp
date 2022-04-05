// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_GHOST_HPP
#define EXACA_GHOST_HPP

#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

void GhostNodesInit(int, int, int DecompositionStrategy, int NeighborRank_North, int NeighborRank_South,
                    int NeighborRank_East, int NeighborRank_West, int NeighborRank_NorthEast,
                    int NeighborRank_NorthWest, int NeighborRank_SouthEast, int NeighborRank_SouthWest, int MyXSlices,
                    int MyYSlices, int MyXOffset, int MyYOffset, int ZBound_Low, int nzActive, int nz,
                    int LocalActiveDomainSize, int NGrainOrientations, ViewI NeighborX, ViewI NeighborY,
                    ViewI NeighborZ, ViewF GrainUnitVector, ViewI GrainID, ViewI CellType, ViewF DOCenter,
                    ViewF DiagonalLength, ViewF CritDiagonalLength);
void GhostNodes2D(int, int, int NeighborRank_North, int NeighborRank_South, int NeighborRank_East,
                  int NeighborRank_West, int NeighborRank_NorthEast, int NeighborRank_NorthWest,
                  int NeighborRank_SouthEast, int NeighborRank_SouthWest, int MyXSlices, int MyYSlices, int MyXOffset,
                  int MyYOffset, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewI CellType, ViewF DOCenter,
                  ViewI GrainID, ViewF GrainUnitVector, ViewF DiagonalLength, ViewF CritDiagonalLength,
                  int NGrainOrientations, Buffer2D BufferWestSend, Buffer2D BufferEastSend, Buffer2D BufferNorthSend,
                  Buffer2D BufferSouthSend, Buffer2D BufferNorthEastSend, Buffer2D BufferNorthWestSend,
                  Buffer2D BufferSouthEastSend, Buffer2D BufferSouthWestSend, Buffer2D BufferWestRecv,
                  Buffer2D BufferEastRecv, Buffer2D BufferNorthRecv, Buffer2D BufferSouthRecv,
                  Buffer2D BufferNorthEastRecv, Buffer2D BufferNorthWestRecv, Buffer2D BufferSouthEastRecv,
                  Buffer2D BufferSouthWestRecv, int BufSizeX, int BufSizeY, int BufSizeZ, int ZBound_Low);
void GhostNodes1D(int, int, int NeighborRank_North, int NeighborRank_South, int MyXSlices, int MyYSlices, int MyXOffset,
                  int MyYOffset, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewI CellType, ViewF DOCenter,
                  ViewI GrainID, ViewF GrainUnitVector, ViewF DiagonalLength, ViewF CritDiagonalLength,
                  int NGrainOrientations, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, Buffer2D BufferNorthRecv,
                  Buffer2D BufferSouthRecv, int BufSizeX, int, int BufSizeZ, int ZBound_Low);

#endif
