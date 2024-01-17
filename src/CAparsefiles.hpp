// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_PARSE_HPP
#define EXACA_PARSE_HPP

#include "CAconfig.hpp"

#include <Kokkos_Core.hpp>

#include <nlohmann/json.hpp>

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
void SwapEndian(SwapType &var) {
    // Cast var into a char array (bit values)
    char *varArray = reinterpret_cast<char *>(&var);
    // Size of char array
    int varSize = sizeof(var);
    // Swap the "ith" bit with the bit "i" from the end of the array
    for (long i = 0; i < static_cast<long>(varSize / 2); i++)
        std::swap(varArray[varSize - 1 - i], varArray[i]);
}
// Reads binary data of the type ReadType, optionally swapping the endian format
template <typename ReadType>
ReadType ReadBinaryData(std::ifstream &instream, bool SwapEndianYN = false) {
    unsigned char temp[sizeof(ReadType)];
    instream.read(reinterpret_cast<char *>(temp), sizeof(ReadType));
    if (SwapEndianYN)
        SwapEndian(temp);
    ReadType readValue = reinterpret_cast<ReadType &>(temp);
    return readValue;
}
// Parse space-separated ASCII data loaded into the string stream
template <typename ReadType>
ReadType ParseASCIIData(std::istringstream &ss) {
    ReadType readValue;
    ss >> readValue;
    return readValue;
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
