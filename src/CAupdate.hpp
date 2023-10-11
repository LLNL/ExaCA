// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_UPDATE_HPP
#define EXACA_UPDATE_HPP

#include "CAcelldata.hpp"
#include "CAfunctions.hpp"
#include "CAghostnodes.hpp"
#include "CAinterfacialresponse.hpp"
#include "CAprint.hpp"
#include "CAtemperature.hpp"
#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

#include <string>

void FillSteeringVector_NoRemelt(int cycle, int DomainSize, Temperature<device_memory_space> &temperature,
                                 CellData<device_memory_space> &cellData, ViewI SteeringVector, ViewI numSteer,
                                 ViewI_H numSteer_Host);
void FillSteeringVector_Remelt(int cycle, int DomainSize, int nx, int ny_local, NList NeighborX, NList NeighborY,
                               NList NeighborZ, Temperature<device_memory_space> &temperature,
                               CellData<device_memory_space> &cellData, int nz_layer, ViewI SteeringVector,
                               ViewI numSteer, ViewI_H numSteer_Host);
void CellCapture(int, int np, int, int nx, int ny_local, InterfacialResponseFunction irf, int y_offset, NList NeighborX,
                 NList NeighborY, NList NeighborZ, ViewF GrainUnitVector, ViewF CritDiagonalLength,
                 ViewF DiagonalLength, CellData<device_memory_space> &cellData,
                 Temperature<device_memory_space> &temperature, ViewF DOCenter, int NGrainOrientations,
                 Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, ViewI SendSizeNorth, ViewI SendSizeSouth,
                 int nz_layer, ViewI SteeringVector, ViewI numSteer, ViewI_H numSteer_Host, bool AtNorthBoundary,
                 bool AtSouthBoundary, int &BufSize);
void JumpTimeStep(int &cycle, unsigned long int RemainingCellsOfInterest, unsigned long int LocalTempSolidCells,
                  Temperature<device_memory_space> &temperature, int DomainSize, int ny_local, int z_layer_bottom,
                  CellData<device_memory_space> &cellData, int id, int layernumber, int np, int nx, int ny,
                  ViewF GrainUnitVector, Print print, int NGrainOrientations, int nz_layer, int nz, double deltax,
                  double XMin, double YMin, double ZMin);
void IntermediateOutputAndCheck(int id, int np, int &cycle, int ny_local, int DomainSize, int nx, int ny, int nz,
                                int nz_layer, int z_layer_bottom, double deltax, double XMin, double YMin, double ZMin,
                                int SuccessfulNucEvents_ThisRank, int &XSwitch, CellData<device_memory_space> &cellData,
                                Temperature<device_memory_space> &temperature, std::string SimulationType,
                                int layernumber, int NGrainOrientations, ViewF GrainUnitVector, Print print);
void IntermediateOutputAndCheck(int id, int cycle, int ny_local, int y_offset, int DomainSize, int nx, int ny, int nz,
                                int &XSwitch, ViewI CellType_AllLayers);

#endif
