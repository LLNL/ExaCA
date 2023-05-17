// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
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

int FindTopOrBottom(int ***LayerID, int XLow, int XHigh, int YLow, int YHigh, int nz, int L, std::string HighLow);

// These are used in reading/parsing ExaCA microstructure data
void ParseLogFile(std::string LogFile, int &nx, int &ny, int &nz, double &deltax, int &NumberOfLayers,
                  std::vector<double> &XYZBounds, std::string &RotationFilename, std::string &EulerAnglesFilename,
                  std::string &RGBFilename, bool OrientationFilesInInput);
void ReadASCIIField(std::ifstream &InputDataStream, int nx, int ny, int nz, ViewI3D_H FieldOfInterest);
void ReadBinaryField(std::ifstream &InputDataStream, int nx, int ny, int nz, ViewI3D_H FieldOfInterest,
                     std::string FieldName);
void InitializeData(std::string MicrostructureFile, int nx, int ny, int nz, ViewI3D_H GrainID, ViewI3D_H LayerID);
void CheckInputFiles(std::string &LogFile, std::string MicrostructureFile, std::string &RotationFilename,
                     std::string &RGBFilename, std::string &EulerAnglesFilename);
double convertToMicrons(double deltax, std::string RegionType);
double convertToCells(double deltax, std::string RegionType);
std::vector<float> getGrainMisorientation(std::string Direction, ViewF_H GrainUnitVector,
                                          std::vector<int> UniqueGrainIDVector, int NumberOfOrientations,
                                          int NumberOfGrains);
void dual_print(std::string temp, std::ostream &stream1, std::ostream &stream2);
template <typename ReturnType, typename FirstType, typename SecondType>
ReturnType DivideCast(FirstType Int1, SecondType Int2) {
    return static_cast<ReturnType>(Int1) / static_cast<ReturnType>(Int2);
}

#endif
