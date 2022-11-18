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

#include "mpi.h"

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
                    int IntermediateFileCounter, int ZBound_Low, int nzActive, double deltax, float XMin, float YMin,
                    float ZMin, int NumberOfLayers, bool PrintBinary, int RVESize = 0);
void PrintCAFields(int nx, int ny, int nz, ViewI3D_H GrainID_WholeDomain, ViewI3D_H LayerID_WholeDomain,
                   ViewI3D_H CritTimeStep_WholeDomain, ViewI3D_H CellType_WholeDomain,
                   ViewF3D_H UndercoolingChange_WholeDomain, ViewF3D_H UndercoolingCurrent_WholeDomain,
                   std::string PathToOutput, std::string BaseFileName, int PrintDebug, bool PrintFullOutput,
                   double deltax, float XMin, float YMin, float ZMin, bool PrintBinary);
void PrintGrainMisorientations(std::string BaseFileName, std::string PathToOutput, int nx, int ny, int nz,
                               ViewI3D_H LayerID_WholeDomain, ViewI3D_H GrainID_WholeDomain, ViewF_H GrainUnitVector,
                               int NGrainOrientations, double deltax, float XMin, float YMin, float ZMin,
                               bool PrintBinary);
void PrintFinalUndercooling(std::string BaseFileName, std::string PathToOutput, int nx, int ny, int nz,
                            ViewF3D_H UndercoolingCurrent_WholeDomain, ViewI3D_H LayerID_WholeDomain, double deltax,
                            float XMin, float YMin, float ZMin, bool PrintBinary);
void PrintExaConstitDefaultRVE(std::string BaseFileName, std::string PathToOutput, int nx, int ny, int nz,
                               ViewI3D_H LayerID_WholeDomain, ViewI3D_H GrainID_WholeDomain, double deltax,
                               int NumberOfLayers, int RVESize);
void PrintIntermediateExaCAState(int IntermediateFileCounter, int layernumber, std::string BaseFileName,
                                 std::string PathToOutput, int ZBound_Low, int nzActive, int nx, int ny,
                                 ViewI3D_H GrainID_WholeDomain, ViewI3D_H CellType_WholeDomain, ViewF_H GrainUnitVector,
                                 int NGrainOrientations, double deltax, float XMin, float YMin, float ZMin,
                                 bool PrintBinary);

//*****************************************************************************/
// Print a log file for this ExaCA run, containing information about the run parameters used
// from the input file as well as the decomposition scheme
template <typename IRFtype>
void PrintExaCALog(int id, int np, std::string InputFile, std::string SimulationType, int MyYSlices, int MyYOffset,
                   IRFtype irf, double deltax, double NMax, double dTN, double dTsigma,
                   std::vector<std::string> temp_paths, int TempFilesInSeries, double HT_deltax, bool RemeltingYN,
                   double deltat, int NumberOfLayers, int LayerHeight, std::string SubstrateFileName,
                   double SubstrateGrainSpacing, bool SubstrateFile, double G, double R, int nx, int ny, int nz,
                   double FractSurfaceSitesActive, std::string PathToOutput, int NSpotsX, int NSpotsY, int SpotOffset,
                   int SpotRadius, std::string BaseFileName, double InitTime, double RunTime, double OutTime, int cycle,
                   double InitMaxTime, double InitMinTime, double NuclMaxTime, double NuclMinTime,
                   double CreateSVMinTime, double CreateSVMaxTime, double CaptureMaxTime, double CaptureMinTime,
                   double GhostMaxTime, double GhostMinTime, double OutMaxTime, double OutMinTime, double XMin,
                   double XMax, double YMin, double YMax, double ZMin, double ZMax) {

    int *YSlices = new int[np];
    int *YOffset = new int[np];
    MPI_Gather(&MyYSlices, 1, MPI_INT, YSlices, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&MyYOffset, 1, MPI_INT, YOffset, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (id == 0) {

        std::string FName = PathToOutput + BaseFileName + ".log";
        std::cout << "Printing ExaCA log file" << std::endl;
        std::ofstream ExaCALog;
        ExaCALog.open(FName);
        ExaCALog << "log file for a simulation run with input file " << InputFile << "  run on " << np
                 << " MPI ranks, output written at cycle " << cycle << std::endl;
        ExaCALog << "This simulation took " << InitTime + RunTime + OutTime
                 << " seconds to run, with the init/run/output breakdown as " << InitTime << "/" << RunTime << "/"
                 << OutTime << std::endl;
        ExaCALog << "This simulation was type: " << SimulationType << std::endl;
        ExaCALog << "Domain size in x: " << nx << std::endl;
        ExaCALog << "Domain size in y: " << ny << std::endl;
        ExaCALog << "Domain size in z: " << nz << std::endl;
        ExaCALog << "Cell size: " << deltax << " microns" << std::endl;
        ExaCALog << "Time step: " << deltat << " microseconds" << std::endl;
        ExaCALog << "Lower bound of domain in x: " << XMin << std::endl;
        ExaCALog << "Lower bound of domain in y: " << YMin << std::endl;
        ExaCALog << "Lower bound of domain in z: " << ZMin << std::endl;
        ExaCALog << "Upper bound of domain in x: " << XMax << std::endl;
        ExaCALog << "Upper bound of domain in y: " << YMax << std::endl;
        ExaCALog << "Upper bound of domain in z: " << ZMax << std::endl;
        ExaCALog << "Nucleation density was " << NMax << " m^-3 , mean nucleation undercooling was " << dTN
                 << " K, and standard deviation of nucleation undercooling was " << dTsigma << " K" << std::endl;
        ExaCALog << irf->print() << std::endl;
        if (SimulationType == "C") {
            ExaCALog << "The thermal gradient was " << G << " K/m, and the cooling rate " << R << " K/s" << std::endl;
            ExaCALog << "The fraction of surface sites active was " << FractSurfaceSitesActive << std::endl;
        }
        else if (SimulationType == "S") {
            ExaCALog << "A total of " << NSpotsX << " in X and " << NSpotsY << " in Y were considered" << std::endl;
            ExaCALog << "The spots were offset by " << SpotOffset << " microns, and had radii of " << SpotRadius
                     << " microns" << std::endl;
            ExaCALog << "This pattern was repeated for " << NumberOfLayers << " layers" << std::endl;
        }
        else {
            ExaCALog << NumberOfLayers << " layers were simulated, with an offset of " << LayerHeight << " cells"
                     << std::endl;
            if (RemeltingYN)
                ExaCALog << "Remelting was included" << std::endl;
            else
                ExaCALog << "Remelting was not included" << std::endl;
            if (SubstrateFile)
                ExaCALog << "The substrate file was " << SubstrateFileName << std::endl;
            else
                ExaCALog << "The mean substrate grain size was " << SubstrateGrainSpacing << " microns" << std::endl;
            if (SimulationType == "R") {
                ExaCALog << "The " << TempFilesInSeries << " temperature file(s) repeated in the " << NumberOfLayers
                         << " layer simulation were: " << std::endl;
                for (int i = 0; i < TempFilesInSeries; i++) {
                    std::cout << temp_paths[i] << std::endl;
                }
                ExaCALog << "The temperature data resolution was " << HT_deltax << " microns" << std::endl;
            }
        }
        ExaCALog << "The decomposition scheme used was a 1D decomposition along the Y direction" << std::endl;
        for (int i = 0; i < np; i++) {
            ExaCALog << "Rank " << i << " contained " << YSlices[i] << " cells in y; subdomain was offset by "
                     << YOffset[i] << " in y" << std::endl;
        }
        ExaCALog << "Max/min rank time initializing data  = " << InitMaxTime << " / " << InitMinTime << " s"
                 << std::endl;
        ExaCALog << "Max/min rank time in CA nucleation   = " << NuclMaxTime << " / " << NuclMinTime << " s"
                 << std::endl;
        ExaCALog << "Max/min rank time in CA steering vector creation = " << CreateSVMaxTime << " / " << CreateSVMinTime
                 << " s" << std::endl;
        ExaCALog << "Max/min rank time in CA cell capture = " << CaptureMaxTime << " / " << CaptureMinTime << " s"
                 << std::endl;
        ExaCALog << "Max/min rank time in CA ghosting     = " << GhostMaxTime << " / " << GhostMinTime << " s"
                 << std::endl;
        ExaCALog << "Max/min rank time exporting data     = " << OutMaxTime << " / " << OutMinTime << " s\n"
                 << std::endl;
        ExaCALog.close();
        // Also print this log information to the console
        std::cout << "===================================================================================" << std::endl;
        std::cout << "Having run with = " << np << " processors" << std::endl;
        std::cout << "Output written at cycle = " << cycle << std::endl;
        std::cout << "Total time = " << InitTime + RunTime + OutTime << std::endl;
        std::cout << "Time spent initializing data = " << InitTime << " s" << std::endl;
        std::cout << "Time spent performing CA calculations = " << RunTime << " s" << std::endl;
        std::cout << "Time spent collecting and printing output data = " << OutTime << " s\n" << std::endl;

        std::cout << "Max/min rank time initializing data  = " << InitMaxTime << " / " << InitMinTime << " s"
                  << std::endl;
        std::cout << "Max/min rank time in CA nucleation   = " << NuclMaxTime << " / " << NuclMinTime << " s"
                  << std::endl;
        std::cout << "Max/min rank time in CA steering vector creation = " << CreateSVMaxTime << " / "
                  << CreateSVMinTime << " s" << std::endl;
        std::cout << "Max/min rank time in CA cell capture = " << CaptureMaxTime << " / " << CaptureMinTime << " s"
                  << std::endl;
        std::cout << "Max/min rank time in CA ghosting     = " << GhostMaxTime << " / " << GhostMinTime << " s"
                  << std::endl;
        std::cout << "Max/min rank time exporting data     = " << OutMaxTime << " / " << OutMinTime << " s\n"
                  << std::endl;

        std::cout << "===================================================================================" << std::endl;
    }
}

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
