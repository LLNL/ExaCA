// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_UPDATE_HPP
#define EXACA_UPDATE_HPP

#include "CAcelldata.hpp"
#include "CAfunctions.hpp"
#include "CAghostnodes.hpp"
#include "CAgrid.hpp"
#include "CAinterfacialresponse.hpp"
#include "CAprint.hpp"
#include "CAtemperature.hpp"
#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

#include <string>

void FillSteeringVector_NoRemelt(int cycle, Grid &grid, Temperature<device_memory_space> &temperature,
                                 CellData<device_memory_space> &cellData, ViewI SteeringVector, ViewI numSteer,
                                 ViewI_H numSteer_Host);
void FillSteeringVector_Remelt(int cycle, Grid &grid, Temperature<device_memory_space> &temperature,
                               CellData<device_memory_space> &cellData, ViewI SteeringVector, ViewI numSteer,
                               ViewI_H numSteer_Host);
void CellCapture(int, int np, int, Grid &grid, InterfacialResponseFunction irf, ViewF GrainUnitVector,
                 ViewF CritDiagonalLength, ViewF DiagonalLength, CellData<device_memory_space> &cellData,
                 Temperature<device_memory_space> &temperature, ViewF DOCenter, int NGrainOrientations,
                 Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, ViewI SendSizeNorth, ViewI SendSizeSouth,
                 ViewI SteeringVector, ViewI numSteer, ViewI_H numSteer_Host, int &BufSize);
void JumpTimeStep(int &cycle, unsigned long int RemainingCellsOfInterest, unsigned long int LocalTempSolidCells,
                  Temperature<device_memory_space> &temperature, Grid &grid, CellData<device_memory_space> &cellData,
                  int id, int layernumber, int np, ViewF GrainUnitVector, Print print, int NGrainOrientations);
void IntermediateOutputAndCheck(int id, int np, int &cycle, Grid &grid, int SuccessfulNucEvents_ThisRank, int &XSwitch,
                                CellData<device_memory_space> &cellData, Temperature<device_memory_space> &temperature,
                                std::string SimulationType, int layernumber, int NGrainOrientations,
                                ViewF GrainUnitVector, Print print);
void IntermediateOutputAndCheck(int id, int cycle, Grid &grid, int &XSwitch, ViewI CellType_AllLayers);

#endif
