// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef GA_UTIL_HPP
#define GA_UTIL_HPP

#include "ExaCA.hpp"

#include <Kokkos_Core.hpp>

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
void ParseLogFile(std::string LogFile, int &nx, int &ny, int &nz, double &deltax, int &NumberOfLayers,
                  bool UseXYZBounds, std::vector<double> &XYZBounds);
void ReadASCIIField(std::ifstream &InputDataStream, int nx, int ny, int nz, ViewI3D_H FieldOfInterest);
void ReadBinaryField(std::ifstream &InputDataStream, int nx, int ny, int nz, ViewI3D_H FieldOfInterest,
                     std::string FieldName);
void ParseFilenames(std::string AnalysisFile, std::string &LogFile, std::string &MicrostructureFile,
                    std::string &RotationFilename, std::string &OutputFileName, std::string &EulerAnglesFilename,
                    bool &NewOrientationFormatYN, std::string &RGBFilename);
void InitializeData(std::string MicrostructureFile, int nx, int ny, int nz, ViewI3D_H GrainID, ViewI3D_H LayerID);
void ParseAnalysisFile(std::string AnalysisFile, std::string RotationFilename, int &NumberOfOrientations,
                       bool *AnalysisTypes, std::vector<int> &XLow_RVE, std::vector<int> &XHigh_RVE,
                       std::vector<int> &YLow_RVE, std::vector<int> &YHigh_RVE, std::vector<int> &ZLow_RVE,
                       std::vector<int> &ZHigh_RVE, int &NumberOfRVEs, std::vector<std::string> &CrossSectionPlane,
                       std::vector<int> &CrossSectionLocation, int &NumberOfCrossSections, int &XMin, int &XMax,
                       int &YMin, int &YMax, int &ZMin, int &ZMax, int nx, int ny, int nz, ViewI3D_H LayerID,
                       int NumberOfLayers, std::vector<bool> &PrintSectionPF, std::vector<bool> &PrintSectionIPF,
                       std::vector<bool> &BimodalAnalysis, bool NewOrientationFormatYN,
                       std::vector<std::string> &CSLabels);
std::vector<int> FindUniqueGrains(const std::vector<int> GrainIDVector);

void CheckInputFiles(std::string BaseFileName, std::string &LogFile, std::string &MicrostructureFile,
                     std::string &RotationFilename, std::string &RGBFilename, std::string &EulerAnglesFilename);
std::vector<int> FindUniqueGrains(std::vector<int> GrainIDVector);
template <typename ReturnType, typename FirstType, typename SecondType>
ReturnType DivideCast(FirstType Int1, SecondType Int2) {
    return static_cast<ReturnType>(Int1) / static_cast<ReturnType>(Int2);
}

#endif
