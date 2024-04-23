// Copyright 2021-2024 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_PARSE_HPP
#define EXACA_PARSE_HPP

#include "CAconfig.hpp"

#include <Kokkos_Core.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

std::string removeWhitespace(std::string line, int pos = -1);
bool getInputBool(std::string val_input);
int getInputInt(std::string val_input);
float getInputFloat(std::string val_input, int factor = 0);
double getInputDouble(std::string val_input, int factor = 0);
void splitString(std::string line, std::vector<std::string> &parsed_line, std::size_t expected_num_values,
                 char separator = ',');
std::size_t checkForHeaderValues(std::string header_line);
bool checkFileExists(const std::string path, const int id, const bool error = true);
std::string checkFileInstalled(const std::string name, const int id);
void checkFileNotEmpty(std::string testfilename);
bool checkTemperatureFileFormat(std::string tempfile_thislayer);
std::size_t checkForHeaderValues(std::string header_line);

// Swaps bits for a variable of type SwapType
template <typename SwapType>
void swapEndian(SwapType &var) {
    // Cast var into a char array (bit values)
    char *var_array = reinterpret_cast<char *>(&var);
    // Size of char array
    int var_size = sizeof(var);
    // Swap the "ith" bit with the bit "i" from the end of the array
    for (long i = 0; i < static_cast<long>(var_size / 2); i++)
        std::swap(var_array[var_size - 1 - i], var_array[i]);
}
// Reads binary data of the type ReadType, optionally swapping the endian format
template <typename ReadType>
ReadType readBinaryData(std::ifstream &instream, bool swap_endian_yn = false) {
    unsigned char temp[sizeof(ReadType)];
    instream.read(reinterpret_cast<char *>(temp), sizeof(ReadType));
    if (swap_endian_yn)
        swapEndian(temp);
    ReadType read_value = reinterpret_cast<ReadType &>(temp);
    return read_value;
}
// Parse space-separated ASCII data loaded into the string stream
template <typename ReadType>
ReadType parseASCIIData(std::istringstream &ss) {
    ReadType read_value;
    ss >> read_value;
    return read_value;
}

// Reads portion of a paraview file and places data in the appropriate data structure
// ASCII data at each Z value is separated by a newline
template <typename read_view_type_3d_host>
read_view_type_3d_host readASCIIField(std::ifstream &input_data_stream, int nx, int ny, int nz, std::string label) {
    read_view_type_3d_host field_of_interest(Kokkos::ViewAllocateWithoutInitializing(label), nz, nx, ny);
    using value_type = typename read_view_type_3d_host::value_type;
    for (int k = 0; k < nz; k++) {
        // Get line from file
        std::string line;
        getline(input_data_stream, line);
        // Parse string at spaces
        std::istringstream ss(line);
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                field_of_interest(k, i, j) = parseASCIIData<value_type>(ss);
            }
        }
    }
    return field_of_interest;
}

// Reads binary string of type read_datatype from a paraview file, converts field to the appropriate type to match
// read_view_type_3d_host (i.e, value_type), and place data in the appropriate data structure Each field consists of a
// single binary string (no newlines) Store converted values in view - LayerID data is a short int, GrainID data is an
// int In some older vtk files, LayerID may have been stored as an int and should be converted
template <typename read_view_type_3d_host, typename read_datatype>
read_view_type_3d_host readBinaryField(std::ifstream &input_data_stream, int nx, int ny, int nz, std::string label) {
    read_view_type_3d_host field_of_interest(Kokkos::ViewAllocateWithoutInitializing(label), nz, nx, ny);
    using value_type = typename read_view_type_3d_host::value_type;
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                read_datatype parsed_value = readBinaryData<read_datatype>(input_data_stream, true);
                field_of_interest(k, i, j) = static_cast<value_type>(parsed_value);
            }
        }
    }
    return field_of_interest;
}

void readIgnoreASCIIField(std::ifstream &input_data_stream, int nx, int ny, int nz);

// Reads and discards binary string of type read_datatype from a paraview file
template <typename read_datatype>
void readIgnoreBinaryField(std::ifstream &input_data_stream, int nx, int ny, int nz) {
    unsigned char temp[sizeof(read_datatype)];
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                input_data_stream.read(reinterpret_cast<char *>(temp), sizeof(read_datatype));
            }
        }
    }
}

#endif
