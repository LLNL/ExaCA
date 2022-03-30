// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_UPDATE_HPP
#define EXACA_UPDATE_HPP

#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

#include <string>

KOKKOS_INLINE_FUNCTION void loadghostnodes_1D(const int GhostGID, const float GhostDOCX, const float GhostDOCY,
                                              const float GhostDOCZ, const float GhostDL, const int BufSizeX,
                                              const int MyYSlices, const int RankX, const int RankY, const int RankZ,
                                              const bool AtNorthBoundary, const bool AtSouthBoundary,
                                              Buffer2D BufferSouthSend, Buffer2D BufferNorthSend);
KOKKOS_INLINE_FUNCTION void
loadghostnodes_2D(const int GhostGID, const float GhostDOCX, const float GhostDOCY, const float GhostDOCZ,
                  const float GhostDL, const int BufSizeX, const int BufSizeY, const int MyXSlices, const int MyYSlices,
                  const int RankX, const int RankY, const int RankZ, const bool AtNorthBoundary,
                  const bool AtSouthBoundary, const bool AtWestBoundary, const bool AtEastBoundary,
                  Buffer2D BufferSouthSend, Buffer2D BufferNorthSend, Buffer2D BufferWestSend, Buffer2D BufferEastSend,
                  Buffer2D BufferNorthEastSend, Buffer2D BufferSouthEastSend, Buffer2D BufferSouthWestSend,
                  Buffer2D BufferNorthWestSend);
void Nucleation(int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int cycle, int &nn, ViewI CellType,
                ViewI NucleiLocations, ViewI NucleationTimes, ViewI GrainID, ViewI GrainOrientation, ViewF DOCenter,
                ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewF GrainUnitVector, ViewF CritDiagonalLength,
                ViewF DiagonalLength, int NGrainOrientations, int PossibleNuclei_ThisRank, int ZBound_Low,
                int layernumber, ViewI LayerID);
void CellCapture(int np, int cycle, int DecompositionStrategy, int LocalActiveDomainSize, int LocalDomainSize,
                 int MyXSlices, int MyYSlices, double AConst, double BConst, double CConst, double DConst,
                 int MyXOffset, int MyYOffset, ViewI2D ItList, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ,
                 ViewI CritTimeStep, ViewF UndercoolingCurrent, ViewF UndercoolingChange, ViewF GrainUnitVector,
                 ViewF CritDiagonalLength, ViewF DiagonalLength, ViewI GrainOrientation, ViewI CellType, ViewF DOCenter,
                 ViewI GrainID, int NGrainOrientations, Buffer2D BufferWestSend, Buffer2D BufferEastSend,
                 Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, Buffer2D BufferNorthEastSend,
                 Buffer2D BufferNorthWestSend, Buffer2D BufferSouthEastSend, Buffer2D BufferSouthWestSend, int BufSizeX,
                 int BufSizeY, int ZBound_Low, int nzActive, int nz, int layernumber, ViewI LayerID,
                 ViewI SteeringVector, ViewI numSteer_G, ViewI_H numSteer_H, bool AtNorthBoundary, bool AtSouthBoundary,
                 bool AtEastBoundary, bool AtWestBoundary);
void IntermediateOutputAndCheck(int id, int np, int &cycle, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset,
                                int LocalDomainSize, int LocalActiveDomainSize, int nx, int ny, int nz, int nzActive,
                                double deltax, float XMin, float YMin, float ZMin, int DecompositionStrategy,
                                int ProcessorsInXDirection, int ProcessorsInYDirection, int nn, int &XSwitch,
                                ViewI CellType, ViewI_H CellType_H, ViewI CritTimeStep, ViewI_H CritTimeStep_H,
                                ViewI GrainID, ViewI_H GrainID_H, std::string TemperatureDataType, int *FinishTimeStep,
                                int layernumber, int, int ZBound_Low, int NGrainOrientations, bool *Melted,
                                ViewI LayerID, ViewI_H LayerID_H, ViewI_H GrainOrientation_H, ViewF_H GrainUnitVector_H,
                                ViewF_H UndercoolingChange_H, ViewF_H UndercoolingCurrent_H, std::string PathToOutput,
                                std::string OutputFile, bool PrintIdleMovieFrames, int MovieFrameInc,
                                int &IntermediateFileCounter);

#endif
