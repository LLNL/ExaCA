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

void fill_steering_vector_no_remelt(int cycle, Grid &grid, CellData<device_memory_space> &cellData,
                                    Temperature<device_memory_space> &temperature,
                                    Interface<device_memory_space> &interface);
void fill_steering_vector_remelt(int cycle, Grid &grid, CellData<device_memory_space> &cellData,
                                 Temperature<device_memory_space> &temperature,
                                 Interface<device_memory_space> &interface);
void cell_capture(int, int np, Grid &grid, InterfacialResponseFunction &irf, CellData<device_memory_space> &cellData,
                  Temperature<device_memory_space> &temperature, Interface<device_memory_space> &interface,
                  ViewF GrainUnitVector, int NGrainOrientations);
void check_buffers(int id, Grid &grid, CellData<device_memory_space> &cellData,
                   Interface<device_memory_space> &interface, int NGrainOrientations);
void refill_buffers(Grid &grid, CellData<device_memory_space> &cellData, Interface<device_memory_space> &interface,
                    int NGrainOrientations);
void halo_update(int, int, Grid &grid, CellData<device_memory_space> &cellData,
                 Interface<device_memory_space> &interface, int NGrainOrientations, ViewF GrainUnitVector);
void JumpTimeStep(int &cycle, unsigned long int RemainingCellsOfInterest, unsigned long int LocalTempSolidCells,
                  Temperature<device_memory_space> &temperature, Grid &grid, CellData<device_memory_space> &cellData,
                  int id, int layernumber, int np, ViewF GrainUnitVector, Print print, int NGrainOrientations);
void IntermediateOutputAndCheck(int id, int np, int &cycle, Grid &grid, int SuccessfulNucEvents_ThisRank, int &XSwitch,
                                CellData<device_memory_space> &cellData, Temperature<device_memory_space> &temperature,
                                std::string SimulationType, int layernumber, int NGrainOrientations,
                                ViewF GrainUnitVector, Print print);
void IntermediateOutputAndCheck(int id, int cycle, Grid &grid, int &XSwitch, ViewI CellType_AllLayers);

#endif
