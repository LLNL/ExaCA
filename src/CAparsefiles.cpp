// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

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
std::string removeWhitespace(std::string line, int pos = -1) {

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
                        int NumInputs = 1) {

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
float getInputFloat(std::string val_input, int factor = 0) {
    float FloatFromString = atof(val_input.c_str()) * pow(10, factor);
    return FloatFromString;
}

// Convert string "val_input" to double value multipled by 10^(factor)
double getInputDouble(std::string val_input, int factor = 0) {
    double DoubleFromString = std::stod(val_input.c_str()) * pow(10, factor);
    return DoubleFromString;
}

// Given a line "s", parse at the commas and return the parsed values as strings in "ParsedLine"
// If AllColumns = true, return all 6 values; otherwise, only parse the first 3 commas/values
void splitString(std::string line, std::vector<std::string> &parsed_line, std::string separator = ",") {
    std::size_t line_size = parsed_line.size();
    for (std::size_t n = 0; n < line_size - 1; n++) {
        std::size_t pos = line.find(separator);
        parsed_line[n] = line.substr(0, pos);
        line = line.substr(pos + 1, std::string::npos);
    }
    parsed_line[line_size - 1] = line;
}

// Check to make sure that all expected column names appear in the header for this temperature file
void checkForHeaderValues(std::string header_line) {

    // Header values from file
    std::size_t header_size = 6;
    std::vector<std::string> header_values(header_size, "");
    splitString(header_line, header_values);

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

// From comma separated data on this line, obtain the x, y, and z coordinates
// if AllColumns = true, also obtain the melting, liquidus, and cooling rate values
void getTemperatureDataPoint(std::string s, std::vector<double> &XYZTemperaturePoint) {

    // temperature values from file as strings
    std::size_t point_size = XYZTemperaturePoint.size();
    std::vector<std::string> TemperatureValues_Read(point_size, "");
    splitString(s, TemperatureValues_Read);

    // convert to double
    for (std::size_t n = 0; n < point_size; n++)
        XYZTemperaturePoint[n] = stod(TemperatureValues_Read[n]);
}

void getTemperatureFilePaths(const std::string path, const std::string name, std::vector<std::string> &all_paths) {
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

bool checkFileExists(const std::string path, const std::string type, const int id, const bool error = true) {
    std::ifstream stream;
    stream.open(path);
    if (!(stream.is_open())) {
        stream.close();
        if (error)
            throw std::runtime_error("Could not locate/open " + type + " file");
        else
            return false;
    }
    stream.close();
    if (id == 0)
        std::cout << type + " file " << path << " opened" << std::endl;
    return true;
}

std::string checkFileInstalled(const std::string name, const std::string type, const int id) {
    // Path to file. Prefer installed location; if not installed use source location.
    std::string path = ExaCA_DATA_INSTALL;
    // Note type must match directory.
    std::string file = path + "/" + type + "/" + name;
    bool files_installed = checkFileExists(file, type, id, false);
    if (!files_installed) {
        path = ExaCA_DATA_SOURCE;
        file = path + "/" + type + "/" + name;
        checkFileExists(file, type, id);
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

void parseMaterialFile(std::string MaterialFile, double &AConst, double &BConst, double &CConst, double &DConst,
                       double &FreezingRange) {
    std::ifstream MaterialData;
    MaterialData.open(MaterialFile);
    skipLines(MaterialData, "*****");
    std::string val;
    // Interfacial response function A, B, C, D, and the solidification range for the alloy
    // The order of these is important: "Alloy freezing range" should be before "A", as a search for "A" in either
    // string will return true
    std::vector<std::string> MaterialInputs = {
        "Alloy freezing range", // Required input 0
        "A",                    // Required input 1
        "B",                    // Required input 2
        "C",                    // Required input 3
        "D",                    // Required input 4
    };
    int NumMaterialInputs = MaterialInputs.size();
    std::vector<std::string> MaterialInputsRead(NumMaterialInputs);
    while (std::getline(MaterialData, val)) {
        // Check if this is one of the expected inputs - otherwise throw an error
        bool FoundInput = parseInputFromList(val, MaterialInputs, MaterialInputsRead, NumMaterialInputs);
        if (!(FoundInput)) {
            std::string error = "Error: Unexpected line " + val + " present in material file " + MaterialFile +
                                " : file should only contain A, B, C, D, and Alloy freezing range";
            throw std::runtime_error(error);
        }
    }

    FreezingRange = getInputDouble(MaterialInputsRead[0]);
    AConst = getInputDouble(MaterialInputsRead[1]);
    BConst = getInputDouble(MaterialInputsRead[2]);
    CConst = getInputDouble(MaterialInputsRead[3]);
    DConst = getInputDouble(MaterialInputsRead[4]);

    MaterialData.close();
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
