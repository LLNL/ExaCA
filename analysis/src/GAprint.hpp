// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef GA_PRINT_HPP
#define GA_PRINT_HPP

#include "CAtypes.hpp"
#include <Kokkos_Core.hpp>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

void PrintExaConstitRVEData(int NumberOfRVEs, std::string BaseFileName, int nx, int ny, int nz, double deltax,
                            ViewI3D_H GrainID, std::vector<int> XLow_RVE, std::vector<int> XHigh_RVE,
                            std::vector<int> YLow_RVE, std::vector<int> YHigh_RVE, std::vector<int> ZLow_RVE,
                            std::vector<int> ZHigh_RVE);
void PrintCrossSectionOrientationData(int NumberOfCrossSections, std::string BaseFileName,
                                      std::vector<int> CrossSectionPlane, std::vector<int> CrossSectionLocation, int nx,
                                      int ny, int nz, int NumberOfOrientations, ViewI3D_H GrainID,
                                      std::vector<bool> PrintSectionPF, std::vector<bool> PrintSectionIPF,
                                      bool NewOrientationFormatYN, double deltax, ViewF_H GrainUnitVectors,
                                      ViewF_H GrainEulerAngles);
void PrintMisorientationData(bool *AnalysisTypes, std::string BaseFileName, int XMin, int XMax, int YMin, int YMax,
                             int ZMin, int ZMax, ViewI3D_H LayerID, ViewF_H GrainUnitVector, ViewI3D_H GrainID,
                             int NumberOfOrientations);
void PrintSizeData(bool *AnalysisTypes, std::string BaseFileName, int XMin, int XMax, int YMin, int YMax, int ZMin,
                   int ZMax, int nx, int ny, int nz, ViewI3D_H LayerID, ViewI3D_H GrainID_WholeDomain, double deltax);
void PrintGrainAreaData(bool *AnalysisTypes, std::string BaseFileName, double deltax, int XMin, int XMax, int YMin,
                        int YMax, int ZMin, int ZMax, ViewI3D_H GrainID);
void PrintPoleFigureData(bool *AnalysisTypes, std::string BaseFileName, int NumberOfOrientations, int XMin, int XMax,
                         int YMin, int YMax, int ZMin, int ZMax, ViewI3D_H GrainID, ViewI3D_H LayerID,
                         bool NewOrientationFormatYN, ViewF_H GrainEulerAngles);
void WritePoleFigureDataToFile(std::string Filename, int NumberOfOrientations, ViewF_H GrainEulerAngles,
                               ViewI_H GOHistogram);
void WritePoleFigureDataToFile_OldFormat(std::string Filename, int NumberOfOrientations, ViewI_H GOHistogram);

#endif
