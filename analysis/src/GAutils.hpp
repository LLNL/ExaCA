// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef GA_UTIL_HPP
#define GA_UTIL_HPP

#include "CAparsefiles.hpp"

#include <Kokkos_Core.hpp>

#include <nlohmann/json.hpp>

#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// These are used in reading/parsing ExaCA microstructure data
void parseLogFile(std::string logfile, int &nx, int &ny, int &nz, double &deltax, int &number_of_layers,
                  std::vector<double> &xyz_bounds, std::string &rotation_filename, bool orientation_files_in_input);
double convertToMicrons(double deltax, std::string region_type);
double convertToCells(double deltax, std::string region_type);
void dualPrint(std::string temp, std::ostream &stream1, std::ostream &stream2);
template <typename ReturnType, typename FirstType, typename SecondType>
ReturnType divideCast(FirstType int1, SecondType int2) {
    return static_cast<ReturnType>(int1) / static_cast<ReturnType>(int2);
}

// Search a specified region of layer_id in x and y for either the smallest Z that doesn't contain any layer "L",
// or the largest Z that doesn't contain any layer "L"
template <typename ViewTypeShort3dHost>
int findTopOrBottom(ViewTypeShort3dHost layer_id, int x_low, int x_high, int y_low, int y_high, int nz, int l,
                    std::string high_low) {

    int top_bottom_z = -1;
    bool searching_for_top = true;
    int i = x_low;
    int j = y_low;
    int k = nz - 2;
    if (high_low == "Low")
        k = 3;
    while (searching_for_top) {
        if (layer_id(k, i, j) != l) {
            j++;
            if (j > y_high) {
                if (i > x_high) {
                    // All X and Y coordinates at this Z have been checked, and this Z coordinate is not part of the
                    // last depositied layer
                    top_bottom_z = k;
                    searching_for_top = false;
                }
                else {
                    i++;
                    j = y_low;
                }
            }
        }
        else {
            // Search next Z coordinate
            if (high_low == "High")
                k--;
            else if (high_low == "Low")
                k++;
            i = x_low;
            j = y_low;
            if ((k == 0) || (k == nz - 1))
                searching_for_top = false;
        }
    }
    if ((high_low == "High") && (k == 0))
        top_bottom_z = nz - 2;
    else if ((high_low == "Low") && (k == nz - 1))
        top_bottom_z = 1;
    return top_bottom_z;
}

template <typename ViewTypeInt3dHost, typename ViewTypeShort3dHost>
void initializeData(std::string microstructure_file, int nx, int ny, int nz, ViewTypeInt3dHost &grain_id,
                    ViewTypeShort3dHost &layer_id) {

    std::ifstream input_data_stream;
    input_data_stream.open(microstructure_file);
    if (!(input_data_stream))
        throw std::runtime_error("Error: Cannot find ExaCA microstructure file");
    std::string line;
    bool binary_vtk = false;
    // 8 header lines at beginning of Paraview file from header
    // 3rd line tells whether data is in ASCII or binary format
    for (int header_lines = 0; header_lines < 8; header_lines++) {
        getline(input_data_stream, line);
        if (header_lines == 2)
            if (line.find("BINARY") != std::string::npos)
                binary_vtk = true;
    }
    for (int field = 0; field < 3; field++) {
        // This line says which variable appears next in the file, along with its type
        // A blank line ends the file read
        std::string read_datatype_string, read_fieldname;
        getline(input_data_stream, line);
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

            // grain_id and layer_id are the only fields currently used by the analysis, other fields will be
            // skipped/ignored
            std::vector<std::string> possible_fieldnames = {"GrainID", "LayerID"};
            int num_possible_fieldnames = possible_fieldnames.size();
            for (auto n = 0; n < num_possible_fieldnames; n++) {
                if (line.find(possible_fieldnames[n]) != std::string::npos)
                    read_fieldname = possible_fieldnames[n];
            }
            if (read_fieldname.empty()) {
                std::cout << "Reading and ignoring field described by" << line << std::endl;
                if (binary_vtk) {
                    if (read_datatype_string == "short")
                        readIgnoreBinaryField<short>(input_data_stream, nx, ny, nz);
                    else if (read_datatype_string == "int")
                        readIgnoreBinaryField<int>(input_data_stream, nx, ny, nz);
                    else if (read_datatype_string == "float")
                        readIgnoreBinaryField<float>(input_data_stream, nx, ny, nz);
                }
                else
                    readIgnoreASCIIField(input_data_stream, nx, ny, nz);
            }
            else {
                // Valid fieldname to be read
                // 1 more unused line
                getline(input_data_stream, line);
                // Place appropriate data
                if (read_fieldname == "GrainID") {
                    // grain_id data should be type int
                    if (read_datatype_string != "int")
                        throw std::runtime_error("Error: Field GrainID should be data of type int");
                    if (binary_vtk)
                        grain_id = readBinaryField<Kokkos::View<int ***, Kokkos::HostSpace>, int>(input_data_stream, nx,
                                                                                                  ny, nz, "GrainID");
                    else
                        grain_id = readASCIIField<Kokkos::View<int ***, Kokkos::HostSpace>>(input_data_stream, nx, ny,
                                                                                            nz, "GrainID");
                }
                else if (read_fieldname == "LayerID") {
                    // layer_id may be int or short, but is stored as type short
                    if ((read_datatype_string != "int") && (read_datatype_string != "short"))
                        throw std::runtime_error("Error: Field LayerID should be data of type int or short");
                    if (!binary_vtk)
                        layer_id = readASCIIField<Kokkos::View<short ***, Kokkos::HostSpace>>(input_data_stream, nx, ny,
                                                                                              nz, "LayerID");
                    else {
                        if (read_datatype_string == "int")
                            layer_id = readBinaryField<Kokkos::View<short ***, Kokkos::HostSpace>, int>(
                                input_data_stream, nx, ny, nz, "LayerID");
                        else
                            layer_id = readBinaryField<Kokkos::View<short ***, Kokkos::HostSpace>, short>(
                                input_data_stream, nx, ny, nz, "LayerID");
                    }
                }
                std::cout << "Data field " << read_fieldname << " read" << std::endl;
            }
        }
    }
    input_data_stream.close();
}

#endif
