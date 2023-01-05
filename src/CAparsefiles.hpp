// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_PARSE_HPP
#define EXACA_PARSE_HPP

#include "CAtypes.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

void skipLines(std::ifstream &stream, std::string seperator);
std::string removeWhitespace(std::string line, int pos = -1);
bool parseInputFromList(std::string line, std::vector<std::string> InputKeys, std::vector<std::string> &ParsedInputs,
                        int NumInputs = 1);
std::string parseInput(std::ifstream &stream, std::string key);
bool getInputBool(std::string val_input);
int getInputInt(std::string val_input);
float getInputFloat(std::string val_input, int factor = 0);
double getInputDouble(std::string val_input, int factor = 0);
void splitString(std::string line, std::vector<std::string> &parsed_line, int expected_num_values,
                 char separator = ',');
void checkForHeaderValues(std::string header_line);
void getTemperatureDataPoint(std::string s, std::vector<double> &XYZTemperaturePoint);
void parseTInstuctionsFile(int id, const std::string TFieldInstructions, int &TempFilesInSeries, int &NumberOfLayers,
                           int &LayerHeight, double deltax, double &HT_deltax, std::vector<std::string> &temp_paths,
                           bool &RemeltingYN, bool &LayerwiseTempRead);
bool checkFileExists(const std::string path, const int id, const bool error = true);
std::string checkFileInstalled(const std::string name, const int id);
void checkFileNotEmpty(std::string testfilename);
void parseMaterialFile(std::string MaterialFile, double &AConst, double &BConst, double &CConst, double &DConst,
                       double &FreezingRange);
std::string parseCoordinatePair(std::string line, int val);
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
void parseTemperatureData(std::string tempfile_thislayer, double YMin, double deltax, int LowerYBound, int UpperYBound,
                          std::vector<double> &RawData, unsigned int &NumberOfTemperatureDataPoints,
                          bool BinaryInputData);

#endif
