// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef GA_PRINT_HPP
#define GA_PRINT_HPP

#include <Kokkos_Core.hpp>

#include "ExaCA.hpp"
#include "GArepresentativeregion.hpp"

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

void printAnalysisHeader_Area(std::ofstream &QoIs, const int XLow_cells, const int XHigh_cells, const int YLow_cells,
                              const int YHigh_cells, const int ZLow_cells, const int ZHigh_cells,
                              const double XLow_microns, const double XHigh_microns, const double YLow_microns,
                              const double YHigh_microns, const double ZLow_microns, const double ZHigh_microns,
                              const std::string regionName, const std::string regionOrientation);
//****
// Functions for printing data for a region of various type
void printGrainTypeFractions(std::ofstream &QoIs, const int XLow, const int XHigh, const int YLow, const int YHigh,
                             const int ZLow, const int ZHigh, ViewI3D_H GrainID, ViewS3D_H LayerID,
                             int RepresentativeRegionSize);
void printMeanMisorientations(std::ofstream &QoIs, int NumberOfGrains, std::vector<float> GrainMisorientationXVector,
                              std::vector<float> GrainMisorientationYVector,
                              std::vector<float> GrainMisorientationZVector, std::vector<float> GrainSizeVector,
                              double RepresentativeRegionSize_Microns);
void printMeanSize(std::ofstream &QoIs, int NumberOfGrains, double RepresentativeRegionSize_Microns,
                   std::string RegionType, std::string Units);
void printMeanExtent(std::ofstream &QoIs, std::vector<float> GrainExtent, std::string Direction, int NumberOfGrains);
//****
// Functions for printing data for a volume
void printMeanBuildTransAspectRatio(std::ofstream &QoIs, std::vector<float> GrainExtentX,
                                    std::vector<float> GrainExtentY, std::vector<float> GrainExtentZ,
                                    std::vector<float> GrainSizeVector, double RepresentativeRegionSize_Microns,
                                    int NumberOfGrains);
//****
// Functions for printing a series of data for multiple XY cross-sections
void writeAreaSeries(bool PrintWeightedAreas, bool PrintUnweightedAreas, std::string BaseFileName, double deltax,
                     int XMin, int XMax, int YMin, int YMax, int ZMin, int ZMax, ViewI3D_H GrainID,
                     double ZMin_Coordinate);
//****

#endif
