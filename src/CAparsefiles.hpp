// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_PARSE_HPP
#define EXACA_PARSE_HPP

#include "CAconfig.hpp"
#include "CAfunctions.hpp"
#include "CAtypes.hpp"

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
std::vector<bool> getPrintFieldValues(nlohmann::json inputdata, std::string Fieldtype,
                                      std::vector<std::string> Fieldnames_key);
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
std::array<double, 6> parseTemperatureCoordinateMinMax(std::string tempfile_thislayer, bool BinaryInputData);
// Read and parse the temperature file (double precision values in a comma-separated, ASCII format with a header line -
// or a binary string of double precision values), storing the x, y, z, tm, tl, cr values in the RawData vector. Each
// rank only contains the points corresponding to cells within the associated Y bounds. NumberOfTemperatureDataPoints is
// incremented on each rank as data is added to RawData
template <typename view_type_double_host>
void parseTemperatureData(std::string tempfile_thislayer, double YMin, double deltax, int LowerYBound, int UpperYBound,
                          int &NumberOfTemperatureDataPoints, bool BinaryInputData,
                          view_type_double_host &RawTemperatureData) {

    std::ifstream TemperatureFilestream;
    TemperatureFilestream.open(tempfile_thislayer);
    if (BinaryInputData) {
        while (!TemperatureFilestream.eof()) {
            double XTemperaturePoint = ReadBinaryData<double>(TemperatureFilestream);
            double YTemperaturePoint = ReadBinaryData<double>(TemperatureFilestream);
            // If no data was extracted from the stream, the end of the file was reached
            if (!(TemperatureFilestream))
                break;
            // Check the y value from ParsedLine, to check if this point is stored on this rank
            // Check the CA grid positions of the data point to see which rank(s) should store it
            int YInt = round((YTemperaturePoint - YMin) / deltax);
            if ((YInt >= LowerYBound) && (YInt <= UpperYBound)) {
                // This data point is inside the bounds of interest for this MPI rank
                // Store the x and y values in RawData
                RawTemperatureData(NumberOfTemperatureDataPoints) = XTemperaturePoint;
                NumberOfTemperatureDataPoints++;
                RawTemperatureData(NumberOfTemperatureDataPoints) = YTemperaturePoint;
                NumberOfTemperatureDataPoints++;
                // Parse the remaining 4 components (z, tm, tl, cr) from the line and store in RawData
                for (int component = 2; component < 6; component++) {
                    RawTemperatureData(NumberOfTemperatureDataPoints) = ReadBinaryData<double>(TemperatureFilestream);
                    NumberOfTemperatureDataPoints++;
                }
                int RawTemperatureData_extent = RawTemperatureData.extent(0);
                // Adjust size of RawData if it is near full
                if (NumberOfTemperatureDataPoints >= RawTemperatureData_extent - 6) {
                    Kokkos::resize(RawTemperatureData, RawTemperatureData_extent + 1000000);
                }
            }
            else {
                // This data point is inside the bounds of interest for this MPI rank
                // ignore the z, tm, tl, cr values associated with it
                unsigned char temp[4 * sizeof(double)];
                TemperatureFilestream.read(reinterpret_cast<char *>(temp), 4 * sizeof(double));
            }
        }
    }
    else {
        // Get number of columns in this temperature file
        std::string HeaderLine;
        getline(TemperatureFilestream, HeaderLine);
        int vals_per_line = checkForHeaderValues(HeaderLine);
        while (!TemperatureFilestream.eof()) {
            std::vector<std::string> ParsedLine(6); // Each line has an x, y, z, tm, tl, cr
            std::string ReadLine;
            if (!getline(TemperatureFilestream, ReadLine))
                break;
            // Only parse the first 6 columns of the temperature data
            splitString(ReadLine, ParsedLine, vals_per_line);
            // Check the y value from ParsedLine, to check if this point is stored on this rank
            double YTemperaturePoint = getInputDouble(ParsedLine[1]);
            // Check the CA grid positions of the data point to see which rank(s) should store it
            int YInt = round((YTemperaturePoint - YMin) / deltax);
            if ((YInt >= LowerYBound) && (YInt <= UpperYBound)) {
                // This data point is inside the bounds of interest for this MPI rank: Store the x, z, tm, tl, and cr
                // vals inside of RawData, incrementing with each value added
                for (int component = 0; component < 6; component++) {
                    RawTemperatureData(NumberOfTemperatureDataPoints) = getInputDouble(ParsedLine[component]);
                    NumberOfTemperatureDataPoints++;
                }
                // Adjust size of RawData if it is near full
                int RawTemperatureData_extent = RawTemperatureData.extent(0);
                // Adjust size of RawData if it is near full
                if (NumberOfTemperatureDataPoints >= RawTemperatureData_extent - 6) {
                    Kokkos::resize(RawTemperatureData, RawTemperatureData_extent + 1000000);
                }
            }
        }
    }
}

#endif
