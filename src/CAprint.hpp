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

void PrintGrainIDsForExaConstit(std::string FName, int nx, int ny, int nz,
                                std::vector<std::vector<std::vector<int>>> GrainID_WholeDomain, double deltax);
void CollectIntField(std::vector<std::vector<std::vector<int>>> &IntVar_WholeDomain, ViewI_H IntVar, int nx, int ny,
                     int nz, int MyXSlices, int MyYSlices, int np, int *RecvXOffset, int *RecvYOffset, int *RecvXSlices,
                     int *RecvYSlices, int *RBufSize);
void CollectFloatField(std::vector<std::vector<std::vector<float>>> &FloatVar_WholeDomain, ViewF_H FloatVar, int nx,
                       int ny, int nz, int MyXSlices, int MyYSlices, int np, int *RecvXOffset, int *RecvYOffset,
                       int *RecvXSlices, int *RecvYSlices, int *RBufSize);
void CollectBoolField(std::vector<std::vector<std::vector<int>>> &IntVar_WholeDomain, bool *BoolVar, int nx, int ny,
                      int nz, int MyXSlices, int MyYSlices, int np, int *RecvXOffset, int *RecvYOffset,
                      int *RecvXSlices, int *RecvYSlices, int *RBufSize);
void SendIntField(ViewI_H VarToSend, int nz, int MyXSlices, int MyYSlices, int SendBufSize);
void SendFloatField(ViewF_H VarToSend, int nz, int MyXSlices, int MyYSlices, int SendBufSize);
void SendBoolField(bool *VarToSend, int nz, int MyXSlices, int MyYSlices, int SendBufSize);
void PrintExaCAData(int id, int np, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection,
                    int ProcessorsInYDirection, ViewI_H GrainID, ViewI_H GrainOrientation, ViewI_H CritTimeStep,
                    ViewF_H GrainUnitVector, ViewI_H LayerID, ViewI_H CellType, ViewF_H UndercoolingChange,
                    ViewF_H UndercoolingCurrent, std::string BaseFileName, int DecompositionStrategy,
                    int NGrainOrientations, bool *Melted, std::string PathToOutput, int PrintDebug,
                    bool PrintMisorientation, bool PrintFullOutput);
void PrintExaCALog(int id, int np, std::string InputFile, std::string SimulationType, int DecompositionStrategy,
                   int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, double AConst, double BConst,
                   double CConst, double DConst, double FreezingRange, double deltax, double NMax, double dTN,
                   double dTsigma, std::string tempfile, int TempFilesInSeries, double HT_deltax, bool RemeltingYN,
                   double deltat, int NumberOfLayers, int LayerHeight, std::string SubstrateFileName,
                   double SubstrateGrainSpacing, bool SubstrateFile, double G, double R, int nx, int ny, int nz,
                   double FractSurfaceSitesActive, std::string PathToOutput, int NSpotsX, int NSpotsY, int SpotOffset,
                   int SpotRadius, std::string BaseFileName, double InitTime, double RunTime, double OutTime);
void PrintParaviewGeneric(int nx, int ny, int nz, std::vector<std::vector<std::vector<int>>> GrainID_WholeDomain,
                          std::vector<std::vector<std::vector<int>>> LayerID_WholeDomain,
                          std::vector<std::vector<std::vector<int>>> CritTimeStep_WholeDomain,
                          std::vector<std::vector<std::vector<int>>> CellType_WholeDomain,
                          std::vector<std::vector<std::vector<float>>> UndercoolingChange_WholeDomain,
                          std::vector<std::vector<std::vector<float>>> UndercoolingCurrent_WholeDomain,
                          std::vector<std::vector<std::vector<int>>> Melted_WholeDomain, std::string PathToOutput,
                          std::string BaseFileName, int PrintDebug, bool PrintFullOutput);
void PrintParaview(std::string BaseFileName, std::string PathToOutput, int nx, int ny, int nz,
                   std::vector<std::vector<std::vector<int>>> Melted_WholeDomain,
                   std::vector<std::vector<std::vector<int>>> GrainID_WholeDomain, ViewI_H GrainOrientation,
                   ViewF_H GrainUnitVector, int NGrainOrientations);

#endif
