// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef GA_PRINT_HPP
#define GA_PRINT_HPP

#include <Kokkos_Core.hpp>

#include "ExaCA.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

void writeExaConstitRVE(int NumberOfRVEs, std::string BaseFileName, int nx, int ny, int nz, double deltax,
                        ViewI3D_H GrainID, std::vector<int> XLow_RVE, std::vector<int> XHigh_RVE,
                        std::vector<int> YLow_RVE, std::vector<int> YHigh_RVE, std::vector<int> ZLow_RVE,
                        std::vector<int> ZHigh_RVE);
void AnalyzeCrossSection_Unimodal(std::ofstream &QoIs, std::string BaseFileName, std::string ThisCrossSectionPlane,
                                  double deltax, int NumberOfGrains, int CrossSectionSize, std::vector<int> GrainAreas,
                                  float MinGrainSize_microns);
void AnalyzeCrossSection_Bimodal(std::ofstream &QoIs, std::string BaseFileName, std::string ThisCrossSectionPlane,
                                 double deltax, int NumberOfGrains, int CrossSectionSize,
                                 std::vector<int> UniqueGrainIDs, std::vector<int> GrainAreas, int NumberOfOrientations,
                                 ViewF_H GrainUnitVector, ViewF_H GrainRGBValues, float MinGrainSize,
                                 float SmallLargeCutoff_microns);
void printCrossSectionData(int NumberOfCrossSections, std::string BaseFileName,
                           std::vector<std::string> CrossSectionPlane, std::vector<int> CrossSectionLocation, int nx,
                           int ny, int nz, int NumberOfOrientations, ViewI3D_H GrainID,
                           std::vector<bool> PrintSectionPF, std::vector<bool> PrintSectionIPF,
                           std::vector<bool> BimodalAnalysis, double deltax, ViewF_H GrainUnitVector,
                           ViewF_H GrainEulerAngles, ViewF_H GrainRGBValues, std::vector<std::string> CSLabels);
// TODO: Combine in future release, with printAnalysisHeader calling separate length, area, or volume based print
// functions
void printAnalysisHeader(std::ofstream &QoIs, const int XLow, const int XHigh, const int YLow, const int YHigh,
                         const int ZLow, const int ZHigh, std::vector<double> XYZBounds);
void printAnalysisHeader(std::ofstream &QoIs, const int XLow_cells, const int XHigh_cells, const int YLow_cells,
                         const int YHigh_cells, const int ZLow_cells, const int ZHigh_cells, const double XLow_microns,
                         const double XHigh_microns, const double YLow_microns, const double YHigh_microns,
                         const double ZLow_microns, const double ZHigh_microns, const std::string regionName,
                         const std::string regionOrientation);
//****
// Functions for printing data for a region of various type
void printGrainTypeFractions(std::ofstream &QoIs, const int XLow, const int XHigh, const int YLow, const int YHigh,
                             const int ZLow, const int ZHigh, ViewI3D_H GrainID, ViewI3D_H LayerID,
                             int RepresentativeRegionSize);
void printMeanMisorientations(std::ofstream &QoIs, int NumberOfGrains, std::vector<float> GrainMisorientationXVector,
                              std::vector<float> GrainMisorientationYVector,
                              std::vector<float> GrainMisorientationZVector, std::vector<float> GrainSizeVector,
                              double RepresentativeRegionSize_Microns);
void printMisorientationDataOld(int XMin, int XMax, int YMin, int YMax, int ZMin, int ZMax, ViewI3D_H LayerID,
                                ViewF_H GrainUnitVector, ViewI3D_H GrainID, int NumberOfOrientations);
void printMeanSize(std::ofstream &QoIs, int NumberOfGrains, double RepresentativeRegionSize_Microns,
                   std::string RegionType);
void printMeanExtent(std::ofstream &QoIs, std::vector<float> GrainExtent, std::string Direction, int NumberOfGrains);
//****
// Functions for printing data for a volume
void printMeanBuildTransAspectRatio(std::ofstream &QoIs, std::vector<float> GrainExtentX,
                                    std::vector<float> GrainExtentY, std::vector<float> GrainExtentZ,
                                    std::vector<float> GrainSizeVector, double RepresentativeRegionSize_Microns,
                                    int NumberOfGrains);
void printSizeOld(std::string BaseFileName, int NumberOfGrains, std::vector<float> GrainExtentX,
                  std::vector<float> GrainExtentY, const int XMin, const int XMax, const int YMin, const int YMax,
                  const int ZMax, double deltax, ViewI3D_H GrainID);
//****
// Functions for printing a series of data for multiple XY cross-sections
void writeAreaSeries(bool PrintWeightedAreas, bool PrintUnweightedAreas, std::string BaseFileName, double deltax,
                     int XMin, int XMax, int YMin, int YMax, int ZMin, int ZMax, ViewI3D_H GrainID,
                     double ZMin_Coordinate);
//****
void writePerGrainStats(std::string OutputFileName, std::string RegionType, std::vector<int> UniqueGrainIDVector,
                        std::vector<float> GrainMisorientationXVector, std::vector<float> GrainMisorientationYVector,
                        std::vector<float> GrainMisorientationZVector, std::vector<float> GrainSizeVector,
                        std::vector<float> GrainExtentX, std::vector<float> GrainExtentY,
                        std::vector<float> GrainExtentZ, std::vector<float> BuildTransAspectRatio, bool *AnalysisTypes,
                        int NumberOfGrains, bool PrintIPFRGB, std::vector<float> GrainRed,
                        std::vector<float> GrainGreen, std::vector<float> GrainBlue);
void writeIPFColoredCrossSection(std::string BaseFileName, std::string ThisCrossSectionPlane, std::string Plane,
                                 int Index1Low, int Index1High, int Index2Low, int Index2High,
                                 int CrossSectionOutOfPlaneLocation, ViewI3D_H GrainID, ViewF_H GrainEulerAngles,
                                 double deltax, int NumberOfOrientations);
void writePoleFigure(std::string BaseFileName, std::string RegionLabel, int NumberOfOrientations,
                     ViewF_H GrainEulerAngles, ViewI_H GOHistogram);

#endif
