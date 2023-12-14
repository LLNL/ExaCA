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

// These are used in reading/parsing ExaCA microstructure data
void ParseLogFile(std::string LogFile, int &nx, int &ny, int &nz, double &deltax, int &NumberOfLayers,
                  std::vector<double> &XYZBounds, std::string &RotationFilename, std::string &EulerAnglesFilename,
                  std::string &RGBFilename, bool OrientationFilesInInput);
void CheckInputFiles(std::string &LogFile, std::string MicrostructureFile, std::string &RotationFilename,
                     std::string &RGBFilename, std::string &EulerAnglesFilename);
double convertToMicrons(double deltax, std::string RegionType);
double convertToCells(double deltax, std::string RegionType);
std::vector<float> getGrainMisorientation(std::string Direction,
                                          Kokkos::View<float *, Kokkos::HostSpace> GrainUnitVector,
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

// Search a specified region of LayerID in x and y for either the smallest Z that doesn't contain any layer "L",
// or the largest Z that doesn't contain any layer "L"
template <typename ViewTypeShort3dHost>
int FindTopOrBottom(ViewTypeShort3dHost LayerID, int XLow, int XHigh, int YLow, int YHigh, int nz, int L,
                    std::string HighLow) {

    int TopBottomZ = -1;
    bool SearchingForTop = true;
    int i = XLow;
    int j = YLow;
    int k = nz - 2;
    if (HighLow == "Low")
        k = 3;
    while (SearchingForTop) {
        if (LayerID(k, i, j) != L) {
            j++;
            if (j > YHigh) {
                if (i > XHigh) {
                    // All X and Y coordinates at this Z have been checked, and this Z coordinate is not part of the
                    // last depositied layer
                    TopBottomZ = k;
                    SearchingForTop = false;
                }
                else {
                    i++;
                    j = YLow;
                }
            }
        }
        else {
            // Search next Z coordinate
            if (HighLow == "High")
                k--;
            else if (HighLow == "Low")
                k++;
            i = XLow;
            j = YLow;
            if ((k == 0) || (k == nz - 1))
                SearchingForTop = false;
        }
    }
    if ((HighLow == "High") && (k == 0))
        TopBottomZ = nz - 2;
    else if ((HighLow == "Low") && (k == nz - 1))
        TopBottomZ = 1;
    return TopBottomZ;
}

template <typename ViewTypeInt3dHost, typename ViewTypeShort3dHost>
void InitializeData(std::string MicrostructureFile, int nx, int ny, int nz, ViewTypeInt3dHost &GrainID,
                    ViewTypeShort3dHost &LayerID) {

    std::ifstream InputDataStream;
    InputDataStream.open(MicrostructureFile);
    if (!(InputDataStream))
        throw std::runtime_error("Error: Cannot find ExaCA microstructure file");
    std::string line;
    bool Binaryvtk = false;
    // 8 header lines at beginning of Paraview file from header
    // 3rd line tells whether data is in ASCII or binary format
    for (int HeaderLines = 0; HeaderLines < 8; HeaderLines++) {
        getline(InputDataStream, line);
        if (HeaderLines == 2)
            if (line.find("BINARY") != std::string::npos)
                Binaryvtk = true;
    }
    for (int field = 0; field < 3; field++) {
        // This line says which variable appears next in the file, along with its type
        // A blank line ends the file read
        std::string read_datatype_string, read_fieldname;
        getline(InputDataStream, line);
        if (line.empty())
            break;
        else {
            // Thus far, every ExaCA field printed to VTK files has been one of these three types - ensure that this
            // field is one of them
            std::vector<std::string> possible_vtk_datatypes = {"short", "int", "float"};
            int num_possible_vtk_types = possible_vtk_datatypes.size();
            for (auto n = 0; n < num_possible_vtk_types; n++) {
                if (line.find(possible_vtk_datatypes[n]) != std::string::npos)
                    read_datatype_string = possible_vtk_datatypes[n];
            }
            if (read_datatype_string.empty())
                throw std::runtime_error("Error: Unknown data type for a field in the vtk file");

            // GrainID and LayerID are the only fields currently used by the analysis, other fields will be
            // skipped/ignored
            std::vector<std::string> possible_fieldnames = {"GrainID", "LayerID"};
            int num_possible_fieldnames = possible_fieldnames.size();
            for (auto n = 0; n < num_possible_fieldnames; n++) {
                if (line.find(possible_fieldnames[n]) != std::string::npos)
                    read_fieldname = possible_fieldnames[n];
            }
            if (read_fieldname.empty()) {
                std::cout << "Reading and ignoring field described by" << line << std::endl;
                if (Binaryvtk) {
                    if (read_datatype_string == "short")
                        ReadIgnoreBinaryField<short>(InputDataStream, nx, ny, nz);
                    else if (read_datatype_string == "int")
                        ReadIgnoreBinaryField<int>(InputDataStream, nx, ny, nz);
                    else if (read_datatype_string == "float")
                        ReadIgnoreBinaryField<float>(InputDataStream, nx, ny, nz);
                }
                else
                    ReadIgnoreASCIIField(InputDataStream, nx, ny, nz);
            }
            else {
                // Valid fieldname to be read
                // 1 more unused line
                getline(InputDataStream, line);
                // Place appropriate data
                if (read_fieldname == "GrainID") {
                    // GrainID data should be type int
                    if (read_datatype_string != "int")
                        throw std::runtime_error("Error: Field GrainID should be data of type int");
                    if (Binaryvtk)
                        GrainID = ReadBinaryField<Kokkos::View<int ***, Kokkos::HostSpace>, int>(InputDataStream, nx,
                                                                                                 ny, nz, "GrainID");
                    else
                        GrainID = ReadASCIIField<Kokkos::View<int ***, Kokkos::HostSpace>>(InputDataStream, nx, ny, nz,
                                                                                           "GrainID");
                }
                else if (read_fieldname == "LayerID") {
                    // LayerID may be int or short, but is stored as type short
                    if ((read_datatype_string != "int") && (read_datatype_string != "short"))
                        throw std::runtime_error("Error: Field LayerID should be data of type int or short");
                    if (!Binaryvtk)
                        LayerID = ReadASCIIField<Kokkos::View<short ***, Kokkos::HostSpace>>(InputDataStream, nx, ny,
                                                                                             nz, "LayerID");
                    else {
                        if (read_datatype_string == "int")
                            LayerID = ReadBinaryField<Kokkos::View<short ***, Kokkos::HostSpace>, int>(
                                InputDataStream, nx, ny, nz, "LayerID");
                        else
                            LayerID = ReadBinaryField<Kokkos::View<short ***, Kokkos::HostSpace>, short>(
                                InputDataStream, nx, ny, nz, "LayerID");
                    }
                }
                std::cout << "Data field " << read_fieldname << " read" << std::endl;
            }
        }
    }
    InputDataStream.close();
}

#endif
