// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "CAtypes.hpp"
#include <Kokkos_Core.hpp>

//*****************************************************************************/

// These are duplicated from CAinitialize.cpp - will not need if the analysis executable is
// linked properly to the ExaCA-Kokkos executable
std::string parseCoordinatePair(std::string line, int val);
int FindTopOrBottom(int ***LayerID, int XLow, int XHigh, int YLow, int YHigh, int nz, int L, std::string HighLow);

// These are used in reading/parsing ExaCA microstructure data
void ParseLogFile(std::string LogFile, int &nx, int &ny, int &nz, double &deltax, int &NumberOfLayers);
void ReadField(std::ifstream &InputDataStream, int nx, int ny, int nz, ViewI3D_H FieldOfInterest);
void ParseFilenames(std::string BaseFileName, std::string &AnalysisFile, std::string &LogFile, std::string &MicrostructureFile, std::string &RotationFilename, std::string &EulerFilename);
void InitializeData(std::string MicrostructureFile, int nx, int ny, int nz, ViewI3D_H GrainID, ViewI3D_H LayerID, ViewI3D_H Melted);
void ParseAnalysisFile(std::string AnalysisFilename, std::string RotationFilename, int &NumberOfOrientations, bool *AnalysisTypes, std::vector<int> &XMin_RVE,
                       std::vector<int> &XMax_RVE, std::vector<int> &YMin_RVE, std::vector<int> &YMax_RVE,
                       std::vector<int> &ZMin_RVE, std::vector<int> &ZMax_RVE, int &NumberOfRVEs,
                       std::vector<int> &ANGCrossSectionPlane, std::vector<int> &ANGCrossSectionLocation,
                       int &NumberOfANGCrossSections, int &XMin, int &XMax, int &YMin, int &YMax, int &ZMin, int &ZMax,
                       int nx, int ny, int nz, ViewI3D_H LayerID, ViewI3D_H Melted, int NumberOfLayers);
void ParseGrainOrientationFiles(std::string RotationFilename, std::string EulerFilename, int NumberOfOrientations,
                                ViewF3D_H GrainUnitVector, ViewF2D_H GrainEulerAngles);

void PrintExaConstitRVEData(int NumberOfRVEs, std::string BaseFileName, int nx, int ny, int nz, double deltax,
                            ViewI3D_H GrainID, std::vector<int> XLow_RVE, std::vector<int> XHigh_RVE,
                            std::vector<int> YLow_RVE, std::vector<int> YHigh_RVE, std::vector<int> ZLow_RVE,
                            std::vector<int> ZHigh_RVE);
void PrintInversePoleFigureCrossSections(int NumberOfCrossSections, std::string BaseFileName,
                                         std::vector<int> CrossSectionPlane, std::vector<int> CrossSectionLocation,
                                         int nx, int ny, int nz, int NumberOfOrientations, ViewI3D_H GrainID,
                                         ViewF2D_H GrainEulerAngles);
void PrintMisorientationData(bool *AnalysisTypes, std::string BaseFileName, int XMin, int XMax, int YMin, int YMax,
                             int ZMin, int ZMax, ViewI3D_H Melted, ViewF3D_H GrainUnitVector, ViewI3D_H GrainID,
                             int NumberOfOrientations);
void PrintSizeData(bool *AnalysisTypes, std::string BaseFileName, int XMin, int XMax, int YMin, int YMax, int ZMin,
                   int ZMax, int nx, int ny, int nz, ViewI3D_H Melted, ViewI3D_H GrainID_WholeDomain, double deltax);
void PrintGrainAreaData(bool *AnalysisTypes, std::string BaseFileName, double deltax, int XMin, int XMax, int YMin,
                        int YMax, int ZMin, int ZMax, ViewI3D_H GrainID);
void PrintPoleFigureData(bool *AnalysisTypes, std::string BaseFileName, int NumberOfOrientations,
                         ViewF2D_H GrainEulerAngles, int XMin, int XMax, int YMin, int YMax, int ZMin, int ZMax,
                         ViewI3D_H GrainID, ViewI3D_H Melted);
