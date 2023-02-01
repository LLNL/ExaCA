// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef GA_UTIL_HPP
#define GA_UTIL_HPP

#include "ExaCA.hpp"

#include <Kokkos_Core.hpp>

#ifdef ExaCA_ENABLE_JSON
#include <nlohmann/json.hpp>
#endif

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

std::string parseCoordinatePair(std::string line, int val);
int FindTopOrBottom(int ***LayerID, int XLow, int XHigh, int YLow, int YHigh, int nz, int L, std::string HighLow);

// These are used in reading/parsing ExaCA microstructure data
int checkLogFormat(std::string LogFile);
#ifdef ExaCA_ENABLE_JSON
void ParseLogFile(std::string LogFile, int &nx, int &ny, int &nz, double &deltax, int &NumberOfLayers,
                  std::vector<double> &XYZBounds, std::string &RotationFilename, std::string &EulerAnglesFilename,
                  std::string &RGBFilename, bool OrientationFilesInInput);
#endif
void ParseLogFile_Old(std::string LogFile, int &nx, int &ny, int &nz, double &deltax, int &NumberOfLayers,
                      std::vector<double> &XYZBounds, std::string &RotationFilename, std::string &EulerAnglesFilename,
                      std::string &RGBFilename, bool OrientationFilesInInput);
void ParseLogFile_OldNoColon(std::string LogFile, int &nx, int &ny, int &nz, double &deltax, int &NumberOfLayers,
                             bool UseXYZBounds, std::vector<double> &XYZBounds);
void ReadASCIIField(std::ifstream &InputDataStream, int nx, int ny, int nz, ViewI3D_H FieldOfInterest);
void ReadBinaryField(std::ifstream &InputDataStream, int nx, int ny, int nz, ViewI3D_H FieldOfInterest,
                     std::string FieldName);
void ParseFilenames(std::string AnalysisFile, std::string &LogFile, std::string &MicrostructureFile,
                    std::string &RotationFilename, std::string &OutputFileName, std::string &EulerAnglesFilename,
                    std::string &RGBFilename, bool &OrientationFilesInInput);
void InitializeData(std::string MicrostructureFile, int nx, int ny, int nz, ViewI3D_H GrainID, ViewI3D_H LayerID);
void ParseAnalysisFile(std::string AnalysisFile, std::string RotationFilename, int &NumberOfOrientations,
                       bool *AnalysisTypes, std::vector<int> &XLow_RVE, std::vector<int> &XHigh_RVE,
                       std::vector<int> &YLow_RVE, std::vector<int> &YHigh_RVE, std::vector<int> &ZLow_RVE,
                       std::vector<int> &ZHigh_RVE, int &NumberOfRVEs, std::vector<std::string> &CrossSectionPlane,
                       std::vector<int> &CrossSectionLocation, int &NumberOfCrossSections, int &XMin, int &XMax,
                       int &YMin, int &YMax, int &ZMin, int &ZMax, int nx, int ny, int nz, ViewI3D_H LayerID,
                       int NumberOfLayers, std::vector<bool> &PrintSectionPF, std::vector<bool> &PrintSectionIPF,
                       std::vector<bool> &BimodalAnalysis, std::vector<std::string> &CSLabels);
ViewI_H createOrientationHistogram(int NumberOfOrientations, ViewI3D_H GrainID, ViewI3D_H LayerID, int XMin, int XMax,
                                   int YMin, int YMax, int ZMin, int ZMax);
ViewI_H createOrientationHistogram(int NumberOfOrientations, std::vector<int> GrainIDVector,
                                   int RepresentativeRegionSize_Cells);
std::vector<int> FindUniqueGrains(const std::vector<int> GrainIDVector);

void CheckInputFiles(std::string &LogFile, std::string MicrostructureFile, std::string &RotationFilename,
                     std::string &RGBFilename, std::string &EulerAnglesFilename);
std::vector<int> FindUniqueGrains(std::vector<int> GrainIDVector);
template <typename ReturnType, typename FirstType, typename SecondType>
ReturnType DivideCast(FirstType Int1, SecondType Int2) {
    return static_cast<ReturnType>(Int1) / static_cast<ReturnType>(Int2);
}

#endif
