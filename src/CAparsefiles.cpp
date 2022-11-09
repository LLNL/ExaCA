// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAparsefiles.hpp"

#include "CAconfig.hpp"
#include "CAfunctions.hpp"

#include "mpi.h"

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

// Skip initial lines in input files.
void skipLines(std::ifstream &stream, std::string seperator) {
    std::string line;
    while (getline(stream, line)) {
        // remove any whitespace from line
        std::string linechars = removeWhitespace(line);
        if (linechars == seperator)
            break;
    }
}

// From a line read from the input file, check against NumInputs possible matches to see if there is a match with a key
// from InputKeys Store the appropriate portion of the line read from the file in ParsedInputs Return whether a match
// was found for the input or not
bool parseInputFromList(std::string line, std::vector<std::string> InputKeys, std::vector<std::string> &ParsedInputs,
                        int NumInputs) {

    // Only attempt to parse the line if it isn't empty - ignore empty lines
    if (line.empty())
        return false;

    // First, check that there is a colon separating the parameter name from the value
    std::size_t colon = line.find(":");
    if (colon == std::string::npos) {
        // No colon on this line - throw error
        std::string error = "Unable to parse line " + line + " ; no colon separating variable name from value";
        throw std::runtime_error(error);
    }
    std::string BeforeColon = line.substr(0, colon);
    // Check line against the "NumInputs" number of possible matches from the vector "InputKeys"
    // Break when first match is found
    bool FoundInput = false;
    for (int inputnumber = 0; inputnumber < NumInputs; inputnumber++) {
        if (BeforeColon.find(InputKeys[inputnumber]) != std::string::npos) {
            // This is required input number "inputnumber"
            FoundInput = true;
            // Take the potion of the line following the colon and store in the appropriate position in "ParsedInputs"
            std::string AfterColon = line.substr(colon + 1, std::string::npos);
            ParsedInputs[inputnumber] = removeWhitespace(AfterColon);
            break;
        }
    }
    return FoundInput;
}

// Get a line from a file, verify the required input was included with the correct format, and parse the input
std::string parseInput(std::ifstream &stream, std::string key) {

    std::string line;
    std::getline(stream, line);
    std::vector<std::string> ParsedInputs(1);
    std::vector<std::string> InputKeys = {key};
    bool FoundInput = parseInputFromList(line, InputKeys, ParsedInputs);
    if (!(FoundInput)) {
        // Keyword not found
        std::string error = "Required input not present: \"" + key + "\" not found in the input file";
        throw std::runtime_error(error);
    }
    return ParsedInputs[0];
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
void splitString(std::string line, std::vector<std::string> &parsed_line, int expected_num_values, char separator) {
    // Make sure the right number of values are present on the line - one more than the number of separators
    int actual_num_values = std::count(line.begin(), line.end(), separator) + 1;
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

// Check to make sure that all expected column names appear in the header for this temperature file
void checkForHeaderValues(std::string header_line) {

    // Header values from file
    std::size_t header_size = 6;
    std::vector<std::string> header_values(header_size, "");
    splitString(header_line, header_values, 6);

    std::vector<std::vector<std::string>> expected_values = {{"x"}, {"y"}, {"z"}, {"tm"}, {"tl", "ts"}, {"r", "cr"}};

    // Case insensitive comparison
    for (std::size_t n = 0; n < header_size; n++) {
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
}

// Function for obtaining the path to temperature data using the old input file format - TODO: remove in future release
void getTemperatureFilePaths_Old(const std::string path, const std::string name, std::vector<std::string> &all_paths) {
    std::size_t num_files = all_paths.size();
    for (std::size_t i = 1; i < num_files + 1; i++) {
        std::string curr_path = path + "/";
        if (num_files > 1)
            curr_path += std::to_string(i) + name;
        else
            curr_path += name;
        all_paths[i - 1] = curr_path;
    }
}
// Function for parsing temperature input data using the file format - TODO: remove in future release
void parseTemperatureInput_Old(std::vector<std::string> DeprecatedInputs, std::vector<std::string> DeprecatedInputsRead,
                               int NumDeprecatedInputs, std::vector<bool> DeprecatedInputs_RequiredYN,
                               int &TempFilesInSeries, int &NumberOfLayers, int &LayerHeight, double deltax,
                               double &HT_deltax, std::vector<std::string> &temp_paths) {

    // Check that the required deprecated temperature input information exists
    for (int i = 0; i < NumDeprecatedInputs; i++) {
        if ((DeprecatedInputsRead[i].empty()) && (DeprecatedInputs_RequiredYN[i])) {
            std::string error = "Error: Missing " + DeprecatedInputs[i] + " line in input file";
            throw std::runtime_error(error);
        }
    }

    // Parse temperature information
    std::string temppath, tempfile;
    // Path to temperature files may be given, otherwise assumed to be in Temperatures folder
    if (DeprecatedInputsRead[0].empty())
        temppath = "examples/Temperatures";
    else
        temppath = DeprecatedInputsRead[0];
    tempfile = DeprecatedInputsRead[1];
    TempFilesInSeries = getInputInt(DeprecatedInputsRead[2]);
    NumberOfLayers = getInputInt(DeprecatedInputsRead[3]);
    LayerHeight = getInputInt(DeprecatedInputsRead[4]);
    // Input temperature data spacing may be given, otherwise it is assumed HT_deltax = CA cell size
    if (DeprecatedInputsRead[5].empty())
        HT_deltax = deltax;
    else
        HT_deltax = getInputDouble(DeprecatedInputsRead[5], -6);

    // Parse deprecated temperature input information into "temp_paths" vector
    temp_paths.resize(TempFilesInSeries, "");
    getTemperatureFilePaths_Old(temppath, tempfile, temp_paths);
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

void parseTInstuctionsFile(int id, const std::string TFieldInstructions, int &TempFilesInSeries, int &NumberOfLayers,
                           int &LayerHeight, double deltax, double &HT_deltax, std::vector<std::string> &temp_paths,
                           bool &RemeltingYN, bool &LayerwiseTempRead) {

    std::ifstream TemperatureData;
    // Check that file exists and contains data
    checkFileExists(TFieldInstructions, id);
    checkFileNotEmpty(TFieldInstructions);
    // Open and read temperature instuctions file
    TemperatureData.open(TFieldInstructions);
    std::string val;
    // These required inputs should be present in the temperature file
    std::vector<std::string> RequiredTemperatureInputs = {
        "Number of layers",
        "Offset between layers",
    };
    // These may or may not be in the temperature file
    std::vector<std::string> OptionalTemperatureInputs = {
        "Discard temperature data and reread temperature files after each layer",
        "Heat transport data mesh size",
    };
    int NumRequiredTemperatureInputs = RequiredTemperatureInputs.size();
    int NumOptionalTemperatureInputs = OptionalTemperatureInputs.size();
    std::vector<std::string> RequiredTemperatureInputsRead(NumRequiredTemperatureInputs);
    std::vector<std::string> OptionalTemperatureInputsRead(NumOptionalTemperatureInputs);
    // Read first portion of the file to get Number of layers, Offset between layers, Heat transport data mesh size (if
    // given)
    bool ReadingArgs = true;
    while (ReadingArgs) {
        std::getline(TemperatureData, val);
        // Line with ***** indicates separation first and second portions of file
        std::string FirstChar = val.substr(0, 1);
        if (FirstChar.compare("*") == 0) {
            ReadingArgs = false;
        }
        else {
            bool FoundArg = parseInputFromList(val, RequiredTemperatureInputs, RequiredTemperatureInputsRead,
                                               NumRequiredTemperatureInputs);
            if (!(FoundArg)) {
                FoundArg = parseInputFromList(val, OptionalTemperatureInputs, OptionalTemperatureInputsRead,
                                              NumOptionalTemperatureInputs);
                if (!(FoundArg))
                    std::cout << "Ignoring unknown line " << val << " in temperature instructions file "
                              << TFieldInstructions << std::endl;
            }
        }
        if (!(TemperatureData.is_open()))
            throw std::runtime_error("Error: Required separator not found in temperature instructions file");
    }
    NumberOfLayers = getInputInt(RequiredTemperatureInputsRead[0]);
    LayerHeight = getInputInt(RequiredTemperatureInputsRead[1]);
    // If this input was not given, default to reading/storing all temperature data during initialization
    if (OptionalTemperatureInputsRead[0].empty())
        LayerwiseTempRead = false;
    else
        LayerwiseTempRead = getInputBool(OptionalTemperatureInputsRead[0]);
    // If this input was not given, default to using CA cell size as the assumed temperature data resolution
    if (OptionalTemperatureInputsRead[1].empty())
        HT_deltax = deltax;
    else
        HT_deltax = getInputDouble(OptionalTemperatureInputsRead[1], -6);

    if ((!(RemeltingYN)) && (LayerwiseTempRead)) {
        if (id == 0)
            std::cout << "Warning: ability to read temperature files one at a time during initialization of "
                         "new layers is only supported for simulations with remelting"
                      << std::endl;
        LayerwiseTempRead = false;
    }
    // Read second part of file to get the paths/names of the temperature data files used
    TempFilesInSeries = 0;
    while (TemperatureData.is_open()) {
        std::string temppath;
        std::getline(TemperatureData, temppath);
        // Ignore any blank lines at the end of the file
        if (!(temppath.empty())) {
            temp_paths.push_back(temppath);
            TempFilesInSeries++;
        }
        else
            break;
    }
    TemperatureData.close();
    if (TempFilesInSeries == 0)
        throw std::runtime_error("Error: No temperature files listed in the temperature instructions file");
}

// Parse a line that looks like [x,y], returning either val = 0 (x) or val = 1 (y)
std::string parseCoordinatePair(std::string line, int val) {
    size_t linestart = line.find("[");
    size_t linebreak = line.find(",");
    size_t lineend = line.find("]");
    std::string SplitStr;
    if (val == 0)
        SplitStr = line.substr(linestart + 1, linebreak - linestart - 1);
    else if (val == 1)
        SplitStr = line.substr(linebreak + 1, lineend - linebreak - 1);
    std::regex r("\\s+");
    SplitStr = std::regex_replace(SplitStr, r, "");
    return SplitStr;
}
