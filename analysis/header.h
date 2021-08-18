// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <stdexcept>
#include <string>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

//*****************************************************************************/

// These are duplicated from CAinitialize.cpp - will not need if the analysis executable is
// linked properly to the ExaCA-Kokkos executable
void skipLines(std::ifstream &stream);
std::string getKey(std::ifstream &stream, std::string &line, std::size_t &colon);
std::string removeWhitespace(std::string line, std::size_t colon);
std::string parseInput(std::ifstream &stream, std::string key);
std::string parseInputMultiple(std::ifstream &stream, std::string key1, std::string key2, int &WhichKey);
bool parseInputBool(std::ifstream &stream, std::string key);

// These are used in reading/parsing ExaCA microstructure data
void ParseLogFile(std::string BaseFileName, int &nx, int &ny, int &nz, double &deltax, int &NumberOfLayers);
void InitializeData(std::string InputFile, int &nx, int &ny, int &nz, int*** GrainID, int*** LayerID, int*** Melted);
void ParseAnalysisFile(std::string OutputAnalysisFile, std::string &RotationFilename, std::string &EulerFilename, int &NumberOfOrientations, bool* AnalysisTypes, std::vector <int> &CenterX_RVE, std::vector <int> &CenterY_RVE, std::vector <int> &CenterZ_RVE, std::vector <int> &Size_RVE, std::vector <bool> &ExaConstitPrint, std::vector <bool> &pyEBSDPrint, int &NumberOfRVEs, std::vector <int> &ANGCrossSectionPlane, std::vector <int> &ANGCrossSectionLocation, int &NumberOfANGCrossSections, int &XMin, int &XMax, int &YMin, int &YMax, int nx, int ny, int nz, int*** LayerID, int*** Melted, int NumberOfLayers);
void ParseGrainOrientationFiles(std::string RotationFilename, std::string EulerFilename, int NumberOfOrientations, float*** GrainUnitVector, float** GrainEulerAngles);


void PrintMisorientationData(bool* AnalysisTypes, std::string BaseFileName, int XMin, int XMax, int YMin, int YMax, int nz, int*** Melted, float*** GrainUnitVector, int*** GrainID, int NumberOfOrientations, int &NumberOfMeltedCells);
void PrintRVEData(int NumberOfRVEs, std::string BaseFileName, int nx, int ny, int nz, double deltax, int*** GrainID, std::vector <int> CenterX_RVE, std::vector <int> CenterY_RVE, std::vector <int> CenterZ_RVE, std::vector <int> Size_RVE, std::vector <bool> ExaConstitPrint, std::vector <bool> pyEBSDPrint);
void PrintVolumeData(bool* AnalysisTypes, int XMin, int XMax, int YMin, int YMax, int nz, int NumberOfMeltedCells, int*** Melted, int*** GrainID_WholeDomain);
void PrintGrainAreaData(bool* AnalysisTypes, std::string BaseFileName, double deltax, int XMin, int XMax, int YMin, int YMax, int nz, int*** GrainID);

