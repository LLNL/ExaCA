// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_GHOST_HPP
#define EXACA_GHOST_HPP

#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

void GhostExchange(int, int, int NeighborRank_North, int NeighborRank_South, int MyXSlices,
                    int MyYSlices, int ZBound_Low, int nzActive, ViewI GrainID, ViewI CellType,
                   ViewF DOCenter, ViewI CaptureTimeStep, ViewF3D BufferSouthSend, ViewF3D BufferNorthSend, ViewF3D BufferSouthRecv, ViewF3D BufferNorthRecv);
//void GhostNodes2D(int, int, int NeighborRank_North, int NeighborRank_South, int NeighborRank_East,
//                  int NeighborRank_West, int NeighborRank_NorthEast, int NeighborRank_NorthWest,
//                  int NeighborRank_SouthEast, int NeighborRank_SouthWest, int MyXSlices, int MyYSlices, int MyXOffset,
//                  int MyYOffset, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewI CellType, ViewF DOCenter,
//                  ViewI GrainID, ViewF GrainUnitVector, ViewI GrainOrientation, ViewF DiagonalLength,
//                  ViewF CritDiagonalLength, int NGrainOrientations, Buffer2D BufferWestSend, Buffer2D BufferEastSend,
//                  Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, Buffer2D BufferNorthEastSend,
//                  Buffer2D BufferNorthWestSend, Buffer2D BufferSouthEastSend, Buffer2D BufferSouthWestSend,
//                  Buffer2D BufferWestRecv, Buffer2D BufferEastRecv, Buffer2D BufferNorthRecv, Buffer2D BufferSouthRecv,
//                  Buffer2D BufferNorthEastRecv, Buffer2D BufferNorthWestRecv, Buffer2D BufferSouthEastRecv,
//                  Buffer2D BufferSouthWestRecv, int BufSizeX, int BufSizeY, int nzActive, int ZBound_Low);
//void GhostNodes1D(int, int, int NeighborRank_North, int NeighborRank_South, int MyXSlices, int MyYSlices, int MyXOffset,
//                  int MyYOffset, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewI CellType, ViewF DOCenter,
//                  ViewI GrainID, ViewF GrainUnitVector, ViewI GrainOrientation, ViewF DiagonalLength,
//                  ViewF CritDiagonalLength, int NGrainOrientations, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend,
//                  Buffer2D BufferNorthRecv, Buffer2D BufferSouthRecv, int BufSizeX, int, int nzActive, int ZBound_Low, int cycle, ViewF UndercoolingCurrent, ViewF UndercoolingChange, double AConst, double BConst, double CConst, double DConst, ViewI CaptureTimeStep);

#endif
