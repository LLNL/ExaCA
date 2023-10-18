// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAparsefiles.hpp"

#include "mpi.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <regex>

// Functions that are used to simplify the parsing of input files, either by ExaCA or related utilities

//*****************************************************************************/
// Remove whitespace from "line", optional argument to take only portion of the line after position "pos"
std::string removeWhitespace(std::string line, int pos) {

    std::string val = line.substr(pos + 1, std::string::npos);
    std::regex r("\\s+");
    val = std::regex_replace(val, r, "");
    return val;
}

// Check if a string is Y (true) or N (false)
bool getInputBool(std::string val_input) {
    std::string val = removeWhitespace(val_input);
    if (val == "N") {
        return false;
    }
    else if (val == "Y") {
        return true;
    }
    else {
        std::string error = "Input \"" + val + "\" must be \"Y\" or \"N\".";
        throw std::runtime_error(error);
    }
}

// Convert string "val_input" to base 10 integer
int getInputInt(std::string val_input) {
    int IntFromString = stoi(val_input, nullptr, 10);
    return IntFromString;
}

// Convert string "val_input" to float value multipled by 10^(factor)
float getInputFloat(std::string val_input, int factor) {
    float FloatFromString = atof(val_input.c_str()) * pow(10, factor);
    return FloatFromString;
}

// Convert string "val_input" to double value multipled by 10^(factor)
double getInputDouble(std::string val_input, int factor) {
    double DoubleFromString = std::stod(val_input.c_str()) * pow(10, factor);
    return DoubleFromString;
}

// Given a string ("line"), parse at "separator" (commas used by default)
// Modifies "parsed_line" to hold the separated values
// expected_num_values may be larger than parsed_line_size, if only a portion of the line is being parsed
void splitString(std::string line, std::vector<std::string> &parsed_line, std::size_t expected_num_values,
                 char separator) {
    // Make sure the right number of values are present on the line - one more than the number of separators
    std::size_t actual_num_values = std::count(line.begin(), line.end(), separator) + 1;
    if (expected_num_values != actual_num_values) {
        std::string error = "Error: Expected " + std::to_string(expected_num_values) +
                            " values while reading file; but " + std::to_string(actual_num_values) + " were found";
        throw std::runtime_error(error);
    }
    // Separate the line into its components, now that the number of values has been checked
    std::size_t parsed_line_size = parsed_line.size();
    for (std::size_t n = 0; n < parsed_line_size - 1; n++) {
        std::size_t pos = line.find(separator);
        parsed_line[n] = line.substr(0, pos);
        line = line.substr(pos + 1, std::string::npos);
    }
    parsed_line[parsed_line_size - 1] = line;
}

// Check to make sure that the 6 expected column names appear in the correct order in the header for this temperature
// file Return the number of columns present - ignore any columns after the 6 of interest
std::size_t checkForHeaderValues(std::string header_line) {

    // Header values from file - number of commas plus one is the size of the header
    std::size_t header_size = std::count(header_line.begin(), header_line.end(), ',') + 1;
    std::vector<std::string> header_values(header_size, "");
    splitString(header_line, header_values, header_size);

    std::vector<std::vector<std::string>> expected_values = {{"x"}, {"y"}, {"z"}, {"tm"}, {"tl", "ts"}, {"r", "cr"}};
    std::size_t num_expected_values = expected_values.size();
    if (num_expected_values > header_size)
        throw std::runtime_error("Error: Fewer values than expected found in temperature file header");

    // Case insensitive comparison
    for (std::size_t n = 0; n < num_expected_values; n++) {
        auto val = removeWhitespace(header_values[n]);
        std::transform(val.begin(), val.end(), val.begin(), ::tolower);
        // Check each header column label against the expected value(s) - throw error if no match
        std::size_t options_size = expected_values[n].size();
        for (std::size_t e = 0; e < options_size; e++) {
            auto ev = expected_values[n][e];
            if (val == ev)
                break;
            else if (e == options_size - 1)
                throw std::runtime_error(ev + " not found in temperature file header");
        }
    }
    return header_size;
}

bool checkFileExists(const std::string path, const int id, const bool error) {
    std::ifstream stream;
    stream.open(path);
    if (!(stream.is_open())) {
        stream.close();
        if (error)
            throw std::runtime_error("Could not locate/open \"" + path + "\"");
        else
            return false;
    }
    stream.close();
    if (id == 0)
        std::cout << "Opened \"" << path << "\"" << std::endl;
    return true;
}

std::string checkFileInstalled(const std::string name, const int id) {
    // Path to file. Prefer installed location; if not installed use source location.
    std::string path = ExaCA_DATA_INSTALL;
    std::string file = path + "/" + name;
    bool files_installed = checkFileExists(file, id, false);
    if (!files_installed) {
        // If full file path, just use it.
        if (name.substr(0, 1) == "/") {
            file = name;
        }
        // If a relative path, it has to be with respect to the source path.
        else {
            path = ExaCA_DATA_SOURCE;
            file = path + "/" + name;
        }
        checkFileExists(file, id);
    }
    return file;
}

// Make sure file contains data
void checkFileNotEmpty(std::string testfilename) {
    std::ifstream testfilestream;
    testfilestream.open(testfilename);
    std::string testline;
    std::getline(testfilestream, testline);
    if (testline.empty())
        throw std::runtime_error("First line of file " + testfilename + " appears empty");
    testfilestream.close();
}

// Check the field names from the given input (Fieldtype = PrintFieldsInit or PrintFieldsFinal) against the possible
// fieldnames listed in Fieldnames_key. Fill the vector PrintFields_given with true or false values depending on whether
// the corresponding field name from Fieldnames_key appeared in the input or not
std::vector<bool> getPrintFieldValues(nlohmann::json inputdata, std::string Fieldtype,
                                      std::vector<std::string> Fieldnames_key) {
    int NumFields_key = Fieldnames_key.size();
    int NumFields_given = inputdata["Printing"][Fieldtype].size();
    std::vector<bool> PrintFields_given(NumFields_key, false);
    // Check each given field against each possible input field name
    for (int field_given = 0; field_given < NumFields_given; field_given++) {
        for (int field_key = 0; field_key < NumFields_key; field_key++) {
            if (inputdata["Printing"][Fieldtype][field_given] == Fieldnames_key[field_key])
                PrintFields_given[field_key] = true;
        }
    }
    return PrintFields_given;
}

// Read x, y, z coordinates in tempfile_thislayer (temperature file in either an ASCII or binary format) and return the
// min and max values
std::array<double, 6> parseTemperatureCoordinateMinMax(std::string tempfile_thislayer, bool BinaryInputData) {

    std::array<double, 6> XYZMinMax;
    std::ifstream TemperatureFilestream;
    TemperatureFilestream.open(tempfile_thislayer);
    std::size_t vals_per_line;

    // Binary temperature data should contain only the six columns of interest
    // Comma-separated double type values may contain additional columns after the 6 used by ExaCA
    if (BinaryInputData)
        vals_per_line = 6;
    else {
        // Read the header line data
        // Make sure the first line contains all required column names: x, y, z, tm, tl, cr
        std::string HeaderLine;
        getline(TemperatureFilestream, HeaderLine);
        vals_per_line = checkForHeaderValues(HeaderLine);
    }

    // Units are assumed to be in meters, meters, seconds, seconds, and K/second
    int XYZPointCount_Estimate = 1000000;
    std::vector<double> XCoordinates(XYZPointCount_Estimate), YCoordinates(XYZPointCount_Estimate),
        ZCoordinates(XYZPointCount_Estimate);
    long unsigned int XYZPointCounter = 0;
    if (BinaryInputData) {
        while (!TemperatureFilestream.eof()) {
            // Get x from the binary string, or, if no data is left, exit the file read
            double XValue = ReadBinaryData<double>(TemperatureFilestream);
            if (!(TemperatureFilestream))
                break;
            // Store the x value that was read, and parse the y and z values
            XCoordinates[XYZPointCounter] = XValue;
            YCoordinates[XYZPointCounter] = ReadBinaryData<double>(TemperatureFilestream);
            ZCoordinates[XYZPointCounter] = ReadBinaryData<double>(TemperatureFilestream);
            // Ignore the tm, tl, cr values associated with this x, y, z
            unsigned char temp[3 * sizeof(double)];
            TemperatureFilestream.read(reinterpret_cast<char *>(temp), 3 * sizeof(double));
            XYZPointCounter++;
            if (XYZPointCounter == XCoordinates.size()) {
                XCoordinates.resize(XYZPointCounter + XYZPointCount_Estimate);
                YCoordinates.resize(XYZPointCounter + XYZPointCount_Estimate);
                ZCoordinates.resize(XYZPointCounter + XYZPointCount_Estimate);
            }
        }
    }
    else {
        while (!TemperatureFilestream.eof()) {
            std::vector<std::string> ParsedLine(3); // Get x, y, z - ignore tm, tl, cr
            std::string ReadLine;
            if (!getline(TemperatureFilestream, ReadLine))
                break;
            splitString(ReadLine, ParsedLine, vals_per_line);
            // Only get x, y, and z values from ParsedLine
            XCoordinates[XYZPointCounter] = getInputDouble(ParsedLine[0]);
            YCoordinates[XYZPointCounter] = getInputDouble(ParsedLine[1]);
            ZCoordinates[XYZPointCounter] = getInputDouble(ParsedLine[2]);
            XYZPointCounter++;
            if (XYZPointCounter == XCoordinates.size()) {
                XCoordinates.resize(XYZPointCounter + XYZPointCount_Estimate);
                YCoordinates.resize(XYZPointCounter + XYZPointCount_Estimate);
                ZCoordinates.resize(XYZPointCounter + XYZPointCount_Estimate);
            }
        }
    }

    XCoordinates.resize(XYZPointCounter);
    YCoordinates.resize(XYZPointCounter);
    ZCoordinates.resize(XYZPointCounter);
    TemperatureFilestream.close();

    // Min/max x, y, and z coordinates from this layer's data
    XYZMinMax[0] = *min_element(XCoordinates.begin(), XCoordinates.end());
    XYZMinMax[1] = *max_element(XCoordinates.begin(), XCoordinates.end());
    XYZMinMax[2] = *min_element(YCoordinates.begin(), YCoordinates.end());
    XYZMinMax[3] = *max_element(YCoordinates.begin(), YCoordinates.end());
    XYZMinMax[4] = *min_element(ZCoordinates.begin(), ZCoordinates.end());
    XYZMinMax[5] = *max_element(ZCoordinates.begin(), ZCoordinates.end());
    return XYZMinMax;
}
