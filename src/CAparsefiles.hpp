// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_PARSE_HPP
#define EXACA_PARSE_HPP

#include "CAtypes.hpp"

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
void splitString(std::string line, std::vector<std::string> &parsed_line, std::string separator = ",");
void checkForHeaderValues(std::string header_line);
void getTemperatureDataPoint(std::string s, std::vector<double> &XYZTemperaturePoint);
void parseTemperatureInput_Old(std::vector<std::string> DeprecatedInputs, std::vector<std::string> DeprecatedInputsRead,
                               int NumDeprecatedInputs, std::vector<bool> DeprecatedInputs_RequiredYN,
                               int &TempFilesInSeries, int &NumberOfLayers, int &LayerHeight, double deltax,
                               double &HT_deltax, std::vector<std::string> &temp_paths);
void parseTInstuctionsFile(int id, const std::string TFieldInstructions, int &TempFilesInSeries, int &NumberOfLayers,
                           int &LayerHeight, double deltax, double &HT_deltax, std::vector<std::string> &temp_paths);
bool checkFileExists(const std::string path, const std::string type, const int id, const bool error = true);
std::string checkFileInstalled(const std::string name, const std::string type, const int id);
void checkFileNotEmpty(std::string testfilename);
void parseMaterialFile(std::string MaterialFile, double &AConst, double &BConst, double &CConst, double &DConst,
                       double &FreezingRange);
std::string parseCoordinatePair(std::string line, int val);
void parseCommaSeparatedArgs(std::string AnalysisFile, std::string ReadLine, std::vector<std::string> &ParsedLine);

#endif
