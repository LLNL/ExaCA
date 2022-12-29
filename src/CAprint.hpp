// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_PRINT_HPP
#define EXACA_PRINT_HPP

#include "CAinterfacialresponse.hpp"
#include "CAparsefiles.hpp"
#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

void WriteHeader(std::ofstream &ParaviewOutputStream, std::string FName, bool PrintBinary, int nx, int ny, int nz,
                 double deltax, double XMin, double YMin, double ZMin);
void CollectIntField(ViewI3D_H IntVar_WholeDomain, ViewI_H IntVar, int nx, int ny, int nz, int MyYSlices, int np,
                     ViewI_H RecvYOffset, ViewI_H RecvYSlices, ViewI_H RBufSize);
void CollectFloatField(ViewF3D_H FloatVar_WholeDomain, ViewF_H FloatVar, int nx, int ny, int nz, int MyYSlices, int np,
                       ViewI_H RecvYOffset, ViewI_H RecvYSlices, ViewI_H RBufSize);
void SendIntField(ViewI_H VarToSend, int nz, int nx, int MyYSlices, int SendBufSize, int SendBufStartY,
                  int SendBufEndY);
void SendFloatField(ViewF_H VarToSend, int nz, int nx, int MyYSlices, int SendBufSize, int SendBufStartY,
                    int SendBufEndY);
void PrintExaCAData(int id, int layernumber, int np, int nx, int ny, int nz, int MyYSlices, int MyYOffset,
                    ViewI GrainID, ViewI CritTimeStep, ViewF GrainUnitVector, ViewI LayerID, ViewI CellType,
                    ViewF UndercoolingChange, ViewF UndercoolingCurrent, std::string BaseFileName,
                    int NGrainOrientations, std::string PathToOutput, int PrintDebug, bool PrintMisorientation,
                    bool PrintFinalUndercooling, bool PrintFullOutput, bool PrintTimeSeries, bool PrintDefaultRVE,
                    int IntermediateFileCounter, int ZBound_Low, int nzActive, double deltax, double XMin, double YMin,
                    double ZMin, int NumberOfLayers, bool PrintBinary, int RVESize = 0);
void PrintExaCALog(int id, int np, std::string InputFile, std::string SimulationType, int MyYSlices, int MyYOffset,
                   InterfacialResponseFunction irf, double deltax, double NMax, double dTN, double dTsigma,
                   std::vector<std::string> temp_paths, int TempFilesInSeries, double HT_deltax, bool RemeltingYN,
                   double deltat, int NumberOfLayers, int LayerHeight, std::string SubstrateFileName,
                   double SubstrateGrainSpacing, bool SubstrateFile, double G, double R, int nx, int ny, int nz,
                   double FractSurfaceSitesActive, std::string PathToOutput, int NSpotsX, int NSpotsY, int SpotOffset,
                   int SpotRadius, std::string BaseFileName, double InitTime, double RunTime, double OutTime, int cycle,
                   double InitMaxTime, double InitMinTime, double NuclMaxTime, double NuclMinTime,
                   double CreateSVMinTime, double CreateSVMaxTime, double CaptureMaxTime, double CaptureMinTime,
                   double GhostMaxTime, double GhostMinTime, double OutMaxTime, double OutMinTime, double XMin,
                   double XMax, double YMin, double YMax, double ZMin, double ZMax, int FirstPowderGrainID,
                   double PowderDensity);
void PrintCAFields(int nx, int ny, int nz, ViewI3D_H GrainID_WholeDomain, ViewI3D_H LayerID_WholeDomain,
                   ViewI3D_H CritTimeStep_WholeDomain, ViewI3D_H CellType_WholeDomain,
                   ViewF3D_H UndercoolingChange_WholeDomain, ViewF3D_H UndercoolingCurrent_WholeDomain,
                   std::string PathToOutput, std::string BaseFileName, int PrintDebug, bool PrintFullOutput,
                   double deltax, double XMin, double YMin, double ZMin, bool PrintBinary);
void PrintGrainMisorientations(std::string BaseFileName, std::string PathToOutput, int nx, int ny, int nz,
                               ViewI3D_H LayerID_WholeDomain, ViewI3D_H GrainID_WholeDomain, ViewF_H GrainUnitVector,
                               int NGrainOrientations, double deltax, double XMin, double YMin, double ZMin,
                               bool PrintBinary);
void PrintFinalUndercooling(std::string BaseFileName, std::string PathToOutput, int nx, int ny, int nz,
                            ViewF3D_H UndercoolingCurrent_WholeDomain, ViewI3D_H LayerID_WholeDomain, double deltax,
                            double XMin, double YMin, double ZMin, bool PrintBinary);
void PrintExaConstitDefaultRVE(std::string BaseFileName, std::string PathToOutput, int nx, int ny, int nz,
                               ViewI3D_H LayerID_WholeDomain, ViewI3D_H GrainID_WholeDomain, double deltax,
                               int NumberOfLayers, int RVESize);
void PrintIntermediateExaCAState(int IntermediateFileCounter, int layernumber, std::string BaseFileName,
                                 std::string PathToOutput, int ZBound_Low, int nzActive, int nx, int ny,
                                 ViewI3D_H GrainID_WholeDomain, ViewI3D_H CellType_WholeDomain, ViewF_H GrainUnitVector,
                                 int NGrainOrientations, double deltax, double XMin, double YMin, double ZMin,
                                 bool PrintBinary);
std::string version();
std::string gitCommitHash();

// Write data of type PrintType as ascii or binary, with option to convert between big and small endian binary
template <typename PrintType>
void WriteData(std::ofstream &outstream, PrintType PrintValue, bool PrintBinary, bool SwapEndianYN = false) {
    if (PrintBinary) {
        if (SwapEndianYN)
            SwapEndian(PrintValue);
        int varSize = sizeof(PrintType);
        outstream.write((char *)&PrintValue, varSize);
    }
    else
        outstream << PrintValue << " ";
}

#endif
