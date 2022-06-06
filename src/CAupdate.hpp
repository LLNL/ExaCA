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
void Nucleation(int cycle, int &SuccessfulNucEvents_ThisRank, int &NucleationCounter, int PossibleNuclei_ThisRank,
                ViewI_H NucleationTimes_H, ViewI NucleiLocations, ViewI NucleiGrainID, ViewI CellType, ViewI GrainID,
                int ZBound_Low, int MyXSlices, int MyYSlices, ViewI SteeringVector, ViewI numSteer_G);
void CellCapture(int np, int cycle, int DecompositionStrategy, int LocalActiveDomainSize, int, int MyXSlices,
                 int MyYSlices, double AConst, double BConst, double CConst, double DConst, int MyXOffset,
                 int MyYOffset, NList NeighborX, NList NeighborY, NList NeighborZ, ViewI CritTimeStep,
                 ViewF UndercoolingCurrent, ViewF UndercoolingChange, ViewF GrainUnitVector, ViewF CritDiagonalLength,
                 ViewF DiagonalLength, ViewI CellType, ViewF DOCenter, ViewI GrainID, int NGrainOrientations,
                 Buffer2D BufferWestSend, Buffer2D BufferEastSend, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend,
                 Buffer2D BufferNorthEastSend, Buffer2D BufferNorthWestSend, Buffer2D BufferSouthEastSend,
                 Buffer2D BufferSouthWestSend, int BufSizeX, int BufSizeY, int ZBound_Low, int nzActive, int,
                 int layernumber, ViewI LayerID, ViewI SteeringVector, ViewI numSteer_G, ViewI_H numSteer_H,
                 bool AtNorthBoundary, bool AtSouthBoundary, bool AtEastBoundary, bool AtWestBoundary);
void IntermediateOutputAndCheck(int id, int np, int &cycle, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset,
                                int LocalDomainSize, int LocalActiveDomainSize, int nx, int ny, int nz, int nzActive,
                                double deltax, float XMin, float YMin, float ZMin, int DecompositionStrategy,
                                int ProcessorsInXDirection, int ProcessorsInYDirection,
                                int SuccessfulNucEvents_ThisRank, int &XSwitch, ViewI CellType, ViewI CritTimeStep,
                                ViewI GrainID, std::string TemperatureDataType, int *FinishTimeStep, int layernumber,
                                int, int ZBound_Low, int NGrainOrientations, bool *Melted, ViewI LayerID,
                                ViewF GrainUnitVector, ViewF UndercoolingChange, ViewF UndercoolingCurrent,
                                std::string PathToOutput, std::string OutputFile, bool PrintIdleMovieFrames,
                                int MovieFrameInc, int &IntermediateFileCounter, int NumberOfLayers);

#endif
