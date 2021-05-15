// Copyright 2021 Lawrence Livermore National Security, LLC and other
// ExaCA Project Developers. See the top-level LICENSE file for details.
// SPDX-License-Identifier: MIT

#ifndef EXACA_GHOST_HPP
#define EXACA_GHOST_HPP

#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

void GhostNodesInit(int id, int np, int DecompositionStrategy, int MyLeft, int MyRight, int MyIn, int MyOut,
                    int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, int MyXSlices, int MyYSlices,
                    int MyXOffset, int MyYOffset, int ZBound_Low, int nzActive, int LocalActiveDomainSize,
                    int NGrainOrientations, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewF GrainUnitVector,
                    ViewI GrainOrientation, ViewI GrainID, ViewI CellType, ViewF DOCenter, ViewF DiagonalLength,
                    ViewF CritDiagonalLength);
void GhostNodes1D(int cycle, int id, int MyLeft, int MyRight, int MyXSlices, int MyYSlices, int MyXOffset,
                  int MyYOffset, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewI CellType, ViewF DOCenter,
                  ViewI GrainID, ViewF GrainUnitVector, ViewI GrainOrientation, ViewF DiagonalLength,
                  ViewF CritDiagonalLength, int NGrainOrientations, Buffer2D BufferA, Buffer2D BufferB,
                  Buffer2D BufferAR, Buffer2D BufferBR, int BufSizeX, int BufSizeY, int BufSizeZ, int ZBound_Low);
void GhostNodes2D(int cycle, int id, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn,
                  int MyLeftOut, int MyRightOut, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset,
                  ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewI CellType, ViewF DOCenter, ViewI GrainID,
                  ViewF GrainUnitVector, ViewI GrainOrientation, ViewF DiagonalLength, ViewF CritDiagonalLength,
                  int NGrainOrientations, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferC, Buffer2D BufferD,
                  Buffer2D BufferE, Buffer2D BufferF, Buffer2D BufferG, Buffer2D BufferH, Buffer2D BufferAR,
                  Buffer2D BufferBR, Buffer2D BufferCR, Buffer2D BufferDR, Buffer2D BufferER, Buffer2D BufferFR,
                  Buffer2D BufferGR, Buffer2D BufferHR, int BufSizeX, int BufSizeY, int BufSizeZ, int ZBound_Low);

#endif
