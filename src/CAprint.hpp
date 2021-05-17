// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_PRINT_HPP
#define EXACA_PRINT_HPP

#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

#include <string>
#include <vector>

void CollectGrainData(int pid, int np, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection,
                      int ProcessorsInYDirection, ViewI_H GrainID, ViewI_H GrainOrientation, ViewF_H GrainUnitVector,
                      std::string OutputFile, int DecompositionStrategy, int NGrainOrientations, bool *Melted,
                      std::string PathToOutput, bool FilesToPrint[4], double deltax);
void PrintTempValues(int pid, int np, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection,
                     int ProcessorsInYDirection, ViewI_H CritTimeStep, int DecompositionStrategy,
                     std::string PathToOutput);
void PrintCT(int pid, int np, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection,
             int ProcessorsInYDirection, ViewI_H CellType, std::string OutputFile, int DecompositionStrategy);
void PrintOrientations(std::string FName, int nx, int ny, int nz,
                       std::vector<std::vector<std::vector<bool>>> Melted_WholeDomain,
                       std::vector<std::vector<std::vector<int>>> GrainID_WholeDomain, int NGrainOrientations);
void PrintGrainIDs(std::string FName, int nx, int ny, int nz,
                   std::vector<std::vector<std::vector<bool>>> Melted_WholeDomain,
                   std::vector<std::vector<std::vector<int>>> GrainID_WholeDomain);
void PrintGrainIDsForExaConstit(std::string FName, int nx, int ny, int nz,
                                std::vector<std::vector<std::vector<int>>> GrainID_WholeDomain, double deltax);
void PrintParaview(std::string FName, int nx, int ny, int nz,
                   std::vector<std::vector<std::vector<bool>>> Melted_WholeDomain,
                   std::vector<std::vector<std::vector<int>>> GrainID_WholeDomain, ViewI_H GrainOrientation,
                   ViewF_H GrainUnitVector, int NGrainOrientations);
void PrintGrainAreas(std::string FName, double deltax, int nx, int ny, int nz,
                     std::vector<std::vector<std::vector<int>>> GrainID_WholeDomain);
void PrintWeightedGrainAreas(std::string FName, double deltax, int nx, int ny, int nz,
                             std::vector<std::vector<std::vector<int>>> GrainID_WholeDomain);

#endif
