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
std::string parseCoordinatePair(std::string line, int val);
int FindTopOrBottom(int*** LayerID, int XLow, int XHigh, int YLow, int YHigh, int nz, int L, std::string HighLow);

// These are used in reading/parsing ExaCA microstructure data
void ParseLogFile(std::string BaseFileName, int &nx, int &ny, int &nz, double &deltax, int &NumberOfLayers);
void ReadField(std::ifstream &InputDataStream, int nx, int ny, int nz, int*** FieldOfInterest);
void InitializeData(std::string InputFile, int nx, int ny, int nz, int*** GrainID, int*** LayerID, int*** Melted);
void ParseAnalysisFile(std::string OutputAnalysisFile, std::string &RotationFilename, std::string &EulerFilename, int &NumberOfOrientations, bool* AnalysisTypes, std::vector <int> &XMin_RVE, std::vector <int> &XMax_RVE, std::vector <int> &YMin_RVE, std::vector <int> &YMax_RVE, std::vector <int> &ZMin_RVE, std::vector <int> &ZMax_RVE, int &NumberOfRVEs, std::vector <int> &ANGCrossSectionPlane, std::vector <int> &ANGCrossSectionLocation, int &NumberOfANGCrossSections, int &XMin, int &XMax, int &YMin, int &YMax, int &ZMin, int &ZMax, int nx, int ny, int nz, int*** LayerID, int*** Melted, int NumberOfLayers);
void ParseGrainOrientationFiles(std::string RotationFilename, std::string EulerFilename, int NumberOfOrientations, float*** GrainUnitVector, float** GrainEulerAngles);

void PrintExaConstitRVEData(int NumberOfRVEs, std::string BaseFileName, int nx, int ny, int nz, double deltax, int*** GrainID, std::vector <int> XLow_RVE, std::vector <int> XHigh_RVE, std::vector <int> YLow_RVE, std::vector <int> YHigh_RVE, std::vector <int> ZLow_RVE, std::vector <int> ZHigh_RVE);
void PrintInversePoleFigureCrossSections(int NumberOfCrossSections, std::string BaseFileName, std::vector <int> CrossSectionPlane, std::vector <int> CrossSectionLocation, int nx, int ny, int nz, int NumberOfOrientations, int*** GrainID, float** GrainEulerAngles);
void PrintMisorientationData(bool* AnalysisTypes, std::string BaseFileName, int XMin, int XMax, int YMin, int YMax, int ZMin, int ZMax, int*** Melted, float*** GrainUnitVector, int*** GrainID, int NumberOfOrientations);
void PrintSizeData(bool* AnalysisTypes, std::string BaseFileName, int XMin, int XMax, int YMin, int YMax, int ZMin, int ZMax, int nx, int ny, int nz, int*** Melted, int*** GrainID_WholeDomain, double deltax);
void PrintGrainAreaData(bool* AnalysisTypes, std::string BaseFileName, double deltax, int XMin, int XMax, int YMin, int YMax, int ZMin, int ZMax, int*** GrainID);
void PrintPoleFigureData(bool* AnalysisTypes, std::string BaseFileName, int NumberOfOrientations, float** GrainEulerAngles, int XMin, int XMax, int YMin, int YMax, int ZMin, int ZMax, int*** GrainID, int*** Melted);
