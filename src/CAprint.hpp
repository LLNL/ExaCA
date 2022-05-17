// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_PRINT_HPP
#define EXACA_PRINT_HPP

#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

#include <string>

void PrintGrainIDsForExaConstit(std::string FName, int nx, int ny, int nz, ViewI3D_H GrainID_WholeDomain,
                                double deltax);
void CollectIntField(ViewI3D_H IntVar_WholeDomain, ViewI_H IntVar, int nx, int ny, int nz, int MyXSlices, int MyYSlices,
                     int np, ViewI_H RecvXOffset, ViewI_H RecvYOffset, ViewI_H RecvXSlices, ViewI_H RecvYSlices,
                     ViewI_H RBufSize);
void CollectFloatField(ViewF3D_H FloatVar_WholeDomain, ViewF_H FloatVar, int nx, int ny, int nz, int MyXSlices,
                       int MyYSlices, int np, ViewI_H RecvXOffset, ViewI_H RecvYOffset, ViewI_H RecvXSlices,
                       ViewI_H RecvYSlices, ViewI_H RBufSize);
void CollectBoolField(ViewI3D_H IntVar_WholeDomain, bool *BoolVar, int nx, int ny, int nz, int MyXSlices, int MyYSlices,
                      int np, ViewI_H RecvXOffset, ViewI_H RecvYOffset, ViewI_H RecvXSlices, ViewI_H RecvYSlices,
                      ViewI_H RBufSize);
void SendIntField(ViewI_H VarToSend, int nz, int MyXSlices, int MyYSlices, int SendBufSize, int SendBufStartX,
                  int SendBufEndX, int SendBufStartY, int SendBufEndY);
void SendFloatField(ViewF_H VarToSend, int nz, int MyXSlices, int MyYSlices, int SendBufSize, int SendBufStartX,
                    int SendBufEndX, int SendBufStartY, int SendBufEndY);
void SendBoolField(bool *VarToSend, int nz, int MyXSlices, int MyYSlices, int SendBufSize, int SendBufStartX,
                   int SendBufEndX, int SendBufStartY, int SendBufEndY);
void PrintExaCAData(int id, int layernumber, int np, int nx, int ny, int nz, int MyXSlices, int MyYSlices,
                    int MyXOffset, int MyYOffset, int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI GrainID,
                    ViewI CritTimeStep, ViewF GrainUnitVector, ViewI LayerID, ViewI CellType, ViewF UndercoolingChange,
                    ViewF UndercoolingCurrent, std::string BaseFileName, int DecompositionStrategy,
                    int NGrainOrientations, bool *Melted, std::string PathToOutput, int PrintDebug,
                    bool PrintMisorientation, bool PrintFinalUndercooling, bool PrintFullOutput, bool PrintTimeSeries,
                    bool PrintDefaultRVE, int IntermediateFileCounter, int ZBound_Low, int nzActive, double deltax,
                    float XMin, float YMin, float ZMin, int NumberOfLayers);
void PrintExaCALog(int id, int np, std::string InputFile, std::string SimulationType, int DecompositionStrategy,
                   int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, double AConst, double BConst,
                   double CConst, double DConst, double FreezingRange, double deltax, double NMax, double dTN,
                   double dTsigma, std::string tempfile, int TempFilesInSeries, double HT_deltax, bool RemeltingYN,
                   double deltat, int NumberOfLayers, int LayerHeight, std::string SubstrateFileName,
                   double SubstrateGrainSpacing, bool SubstrateFile, double G, double R, int nx, int ny, int nz,
                   double FractSurfaceSitesActive, std::string PathToOutput, int NSpotsX, int NSpotsY, int SpotOffset,
                   int SpotRadius, std::string BaseFileName, double InitTime, double RunTime, double OutTime, int cycle,
                   double InitMaxTime, double InitMinTime, double NuclMaxTime, double NuclMinTime,
                   double CaptureMaxTime, double CaptureMinTime, double GhostMaxTime, double GhostMinTime,
                   double OutMaxTime, double OutMinTime);
void PrintCAFields(int nx, int ny, int nz, ViewI3D_H GrainID_WholeDomain, ViewI3D_H LayerID_WholeDomain,
                   ViewI3D_H CritTimeStep_WholeDomain, ViewI3D_H CellType_WholeDomain,
                   ViewF3D_H UndercoolingChange_WholeDomain, ViewF3D_H UndercoolingCurrent_WholeDomain,
                   ViewI3D_H Melted_WholeDomain, std::string PathToOutput, std::string BaseFileName, int PrintDebug,
                   bool PrintFullOutput, double deltax, float XMin, float YMin, float ZMin);
void PrintGrainMisorientations(std::string BaseFileName, std::string PathToOutput, int nx, int ny, int nz,
                               ViewI3D_H Melted_WholeDomain, ViewI3D_H GrainID_WholeDomain, ViewF_H GrainUnitVector,
                               int NGrainOrientations, double deltax, float XMin, float YMin, float ZMin);
void PrintFinalUndercooling(std::string BaseFileName, std::string PathToOutput, int nx, int ny, int nz,
                            ViewI3D_H Melted_WholeDomain, ViewF3D_H UndercoolingCurrent_WholeDomain, double deltax,
                            float XMin, float YMin, float ZMin);
void PrintExaConstitDefaultRVE(std::string BaseFileName, std::string PathToOutput, int nx, int ny, int nz,
                               ViewI3D_H LayerID_WholeDomain, ViewI3D_H GrainID_WholeDomain, double deltax,
                               int NumberOfLayers);
void PrintIntermediateExaCAState(int IntermediateFileCounter, int layernumber, std::string BaseFileName,
                                 std::string PathToOutput, int ZBound_Low, int nzActive, int nx, int ny,
                                 ViewI3D_H GrainID_WholeDomain, ViewI3D_H CellType_WholeDomain, ViewF_H GrainUnitVector,
                                 int NGrainOrientations, double deltax, float XMin, float YMin, float ZMin);

#endif
