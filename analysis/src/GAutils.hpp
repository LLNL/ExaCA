// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef GA_UTIL_HPP
#define GA_UTIL_HPP

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

std::string parseCoordinatePair(std::string line, int val);
int FindTopOrBottom(int ***LayerID, int XLow, int XHigh, int YLow, int YHigh, int nz, int L, std::string HighLow);

// These are used in reading/parsing ExaCA microstructure data
void ParseLogFile(std::string LogFile, int &nx, int &ny, int &nz, double &deltax, int &NumberOfLayers);
void ReadField(std::ifstream &InputDataStream, int nx, int ny, int nz, ViewI3D_H FieldOfInterest);
void ParseFilenames(std::string AnalysisFile, std::string &LogFile, std::string &MicrostructureFile,
                    std::string &RotationFilename, std::string &OutputFileName, std::string &EulerAnglesFilename,
                    bool &NewOrientationFormatYN);
void InitializeData(std::string MicrostructureFile, int nx, int ny, int nz, ViewI3D_H GrainID, ViewI3D_H LayerID,
                    ViewI3D_H Melted);
void ParseAnalysisFile(std::string AnalysisFile, std::string RotationFilename, int &NumberOfOrientations,
                       bool *AnalysisTypes, std::vector<int> &XLow_RVE, std::vector<int> &XHigh_RVE,
                       std::vector<int> &YLow_RVE, std::vector<int> &YHigh_RVE, std::vector<int> &ZLow_RVE,
                       std::vector<int> &ZHigh_RVE, int &NumberOfRVEs, std::vector<int> &CrossSectionPlane,
                       std::vector<int> &CrossSectionLocation, int &NumberOfCrossSections, int &XMin, int &XMax,
                       int &YMin, int &YMax, int &ZMin, int &ZMax, int nx, int ny, int nz, ViewI3D_H LayerID, ViewI3D_H,
                       int NumberOfLayers, std::vector<bool> &PrintSectionPF, std::vector<bool> &PrintSectionIPF,
                       bool NewOrientationFormatYN);

#endif
