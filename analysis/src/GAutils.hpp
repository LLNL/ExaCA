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
void InitializeData(std::string MicrostructureFile, int nx, int ny, int nz, ViewI3D_H &GrainID, ViewS3D_H &LayerID);
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

// Reads portion of a paraview file and places data in the appropriate data structure
// ASCII data at each Z value is separated by a newline
template <typename read_view_type_3d_host>
read_view_type_3d_host ReadASCIIField(std::ifstream &InputDataStream, int nx, int ny, int nz, std::string label) {
    read_view_type_3d_host FieldOfInterest(Kokkos::ViewAllocateWithoutInitializing(label), nz, nx, ny);
    using value_type = typename read_view_type_3d_host::value_type;
    for (int k = 0; k < nz; k++) {
        // Get line from file
        std::string line;
        getline(InputDataStream, line);
        // Parse string at spaces
        std::istringstream ss(line);
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                FieldOfInterest(k, i, j) = ParseASCIIData<value_type>(ss);
            }
        }
    }
    return FieldOfInterest;
}

// Reads binary string of type read_datatype from a paraview file, converts field to the appropriate type to match
// read_view_type_3d_host (i.e, value_type), and place data in the appropriate data structure Each field consists of a
// single binary string (no newlines) Store converted values in view - LayerID data is a short int, GrainID data is an
// int In some older vtk files, LayerID may have been stored as an int and should be converted
template <typename read_view_type_3d_host, typename read_datatype>
read_view_type_3d_host ReadBinaryField(std::ifstream &InputDataStream, int nx, int ny, int nz, std::string label) {
    read_view_type_3d_host FieldOfInterest(Kokkos::ViewAllocateWithoutInitializing(label), nz, nx, ny);
    using value_type = typename read_view_type_3d_host::value_type;
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                read_datatype parsed_value = ReadBinaryData<read_datatype>(InputDataStream, true);
                FieldOfInterest(k, i, j) = static_cast<value_type>(parsed_value);
            }
        }
    }
    return FieldOfInterest;
}

void ReadIgnoreASCIIField(std::ifstream &InputDataStream, int nx, int ny, int nz);

// Reads and discards binary string of type read_datatype from a paraview file
template <typename read_datatype>
void ReadIgnoreBinaryField(std::ifstream &InputDataStream, int nx, int ny, int nz) {
    unsigned char temp[sizeof(read_datatype)];
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                InputDataStream.read(reinterpret_cast<char *>(temp), sizeof(read_datatype));
            }
        }
    }
}

#endif
