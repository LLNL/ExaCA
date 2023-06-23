// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_UPDATE_HPP
#define EXACA_UPDATE_HPP

#include "CAcelldata.hpp"
#include "CAinterfacialresponse.hpp"
#include "CAprint.hpp"
#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

#include <string>

void FillSteeringVector_NoRemelt(int cycle, int LocalActiveDomainSize, int nx, int MyYSlices, ViewI CritTimeStep,
                                 ViewF UndercoolingCurrent, ViewF UndercoolingChange, CellData &cellData,
                                 int ZBound_Low, int layernumber, ViewI SteeringVector, ViewI numSteer_G,
                                 ViewI_H numSteer_H);
void FillSteeringVector_Remelt(int cycle, int LocalActiveDomainSize, int nx, int MyYSlices, NList NeighborX,
                               NList NeighborY, NList NeighborZ, ViewI CritTimeStep, ViewF UndercoolingCurrent,
                               ViewF UndercoolingChange, CellData &cellData, int ZBound_Low, int nzActive,
                               ViewI SteeringVector, ViewI numSteer, ViewI_H numSteer_Host, ViewI MeltTimeStep,
                               ViewI SolidificationEventCounter, ViewI NumberOfSolidificationEvents,
                               ViewF3D LayerTimeTempHistory);
void CellCapture(int id, int np, int cycle, int LocalActiveDomainSize, int LocalDomainSize, int nx, int MyYSlices,
                 InterfacialResponseFunction irf, int MyYOffset, NList NeighborX, NList NeighborY, NList NeighborZ,
                 ViewI CritTimeStep, ViewF UndercoolingCurrent, ViewF UndercoolingChange, ViewF GrainUnitVector,
                 ViewF CritDiagonalLength, ViewF DiagonalLength, CellData &cellData, ViewF DOCenter,
                 int NGrainOrientations, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, ViewI SendSizeNorth,
                 ViewI SendSizeSouth, int ZBound_Low, int nzActive, int nz, ViewI SteeringVector, ViewI numSteer_G,
                 ViewI_H numSteer_H, bool AtNorthBoundary, bool AtSouthBoundary, ViewI SolidificationEventCounter,
                 ViewF3D LayerTimeTempHistory, ViewI NumberOfSolidificationEvents, int &BufSize);
void JumpTimeStep(int &cycle, unsigned long int RemainingCellsOfInterest, unsigned long int LocalTempSolidCells,
                  ViewI MeltTimeStep, int LocalActiveDomainSize, int MyYSlices, int ZBound_Low, CellData &cellData,
                  int id, int layernumber, int np, int nx, int ny, ViewF GrainUnitVector, Print print,
                  int NGrainOrientations, int nzActive, double deltax, double XMin, double YMin, double ZMin);
void IntermediateOutputAndCheck(int id, int np, int &cycle, int MyYSlices, int LocalActiveDomainSize, int nx, int ny,
                                int nzActive, double deltax, double XMin, double YMin, double ZMin,
                                int SuccessfulNucEvents_ThisRank, int &XSwitch, CellData &cellData, ViewI CritTimeStep,
                                std::string TemperatureDataType, int layernumber, int, int ZBound_Low,
                                int NGrainOrientations, ViewF GrainUnitVector, Print print, ViewI MeltTimeStep);

#endif
