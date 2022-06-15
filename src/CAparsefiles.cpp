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

void parseTInstuctionsFile(int id, const std::string TFieldInstructions, int &TempFilesInSeries, int &NumberOfLayers,
                           int &LayerHeight, double deltax, double &HT_deltax, std::vector<std::string> &temp_paths) {

    std::ifstream TemperatureData;
    // Check that file exists and contains data
    checkFileExists(TFieldInstructions, "Temperature instruction", id);
    checkFileNotEmpty(TFieldInstructions);
    // Open and read temperature instuctions file
    TemperatureData.open(TFieldInstructions);
    std::string val;
    // Three required inputs should be present in the temperature file
    std::vector<std::string> TemperatureInputs = {
        "Number of layers",
        "Offset between layers",
        "Heat transport data mesh size",
    };
    int NumTemperatureInputs = TemperatureInputs.size();
    std::vector<std::string> TemperatureInputsRead(NumTemperatureInputs);
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
            bool FoundArg = parseInputFromList(val, TemperatureInputs, TemperatureInputsRead, NumTemperatureInputs);
            if (!(FoundArg))
                std::cout << "Ignoring unknown line " << val << " in temperature instructions file "
                          << TFieldInstructions << std::endl;
        }
        if (!(TemperatureData.is_open()))
            throw std::runtime_error("Error: Required separator not found in temperature instructions file");
    }
    NumberOfLayers = getInputInt(TemperatureInputsRead[0]);
    LayerHeight = getInputInt(TemperatureInputsRead[1]);
    if (TemperatureInputsRead[2].empty())
        HT_deltax = deltax;
    else
        HT_deltax = getInputDouble(TemperatureInputsRead[2], -6);
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

// Initialize grain orientations and unit vectors
void OrientationInit(int, int &NGrainOrientations, ViewF &GrainOrientationData, std::string GrainOrientationFile, int ValsPerLine) {

    // Read file of grain orientations
    std::ifstream O;
    O.open(GrainOrientationFile);

    // Line 1 is the number of orientation values to read (if not specified already)
    std::string ValueRead;
    getline(O, ValueRead);
    NGrainOrientations = getInputInt(ValueRead);
    
    // Temporary host view for storing grain orientations read from file
    ViewF_H GrainOrientationData_Host(Kokkos::ViewAllocateWithoutInitializing("GrainOrientationData_H"), ValsPerLine * NGrainOrientations);
    // Populate data structure for grain orientation data
    for (int i = 0; i < NGrainOrientations; i++) {
        std::string s;
        if (!getline(O, s))
            break;
        std::istringstream ss(s);
        int Comp = 0;
        while (ss) { // This is the 3 grain orientation angles or 9 rotation matrix components
            std::string s;
            if (!getline(ss, s, ','))
                break;
            float ReadGO = atof(s.c_str());
            GrainOrientationData_Host(ValsPerLine * i + Comp) = ReadGO;
            Comp++;
        }
    }
    O.close();

    // Resize device view and orientation data to device
    Kokkos::realloc(GrainOrientationData, ValsPerLine * NGrainOrientations);
    GrainOrientationData = Kokkos::create_mirror_view_and_copy(device_memory_space(), GrainOrientationData_Host);
}
