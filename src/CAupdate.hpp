// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_UPDATE_HPP
#define EXACA_UPDATE_HPP

#include "CAcelldata.hpp"
#include "CAfunctions.hpp"
#include "CAgrid.hpp"
#include "CAinterface.hpp"
#include "CAinterfacialresponse.hpp"
#include "CAprint.hpp"
#include "CAtemperature.hpp"
#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

#include <string>

void fill_steering_vector_no_remelt(const int cycle, const Grid &grid, CellData<device_memory_space> &cellData,
                                    Temperature<device_memory_space> &temperature,
                                    Interface<device_memory_space> &interface);
void fill_steering_vector_remelt(const int cycle, const Grid &grid, CellData<device_memory_space> &cellData,
                                 Temperature<device_memory_space> &temperature,
                                 Interface<device_memory_space> &interface);
void cell_capture(const int, const int np, const Grid &grid, const InterfacialResponseFunction &irf,
                  CellData<device_memory_space> &cellData, Temperature<device_memory_space> &temperature,
                  Interface<device_memory_space> &interface, const ViewF GrainUnitVector, const int NGrainOrientations);
void check_buffers(const int id, const Grid &grid, CellData<device_memory_space> &cellData,
                   Interface<device_memory_space> &interface, const int NGrainOrientations);
void refill_buffers(const Grid &grid, CellData<device_memory_space> &cellData,
                    Interface<device_memory_space> &interface, const int NGrainOrientations);
void halo_update(const int, const int, const Grid &grid, CellData<device_memory_space> &cellData,
                 Interface<device_memory_space> &interface, const int NGrainOrientations, const ViewF GrainUnitVector);
void JumpTimeStep(int &cycle, unsigned long int RemainingCellsOfInterest, unsigned long int LocalTempSolidCells,
                  Temperature<device_memory_space> &temperature, const Grid &grid,
                  CellData<device_memory_space> &cellData, const int id, const int layernumber, const int np,
                  const ViewF GrainUnitVector, Print print, const int NGrainOrientations);
void IntermediateOutputAndCheck(const int id, const int np, int &cycle, const Grid &grid,
                                int SuccessfulNucEvents_ThisRank, int &XSwitch, CellData<device_memory_space> &cellData,
                                Temperature<device_memory_space> &temperature, std::string SimulationType,
                                const int layernumber, const int NGrainOrientations, const ViewF GrainUnitVector,
                                Print print);
void IntermediateOutputAndCheck(const int id, int cycle, const Grid &grid, int &XSwitch, ViewI CellType_AllLayers);

#endif
