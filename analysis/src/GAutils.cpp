// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "GAutils.hpp"

#include <CAparsefiles.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>

// Search a specified region of LayerID in x and y for either the smallest Z that doesn't contain any layer "L",
// or the largest Z that doesn't contain any layer "L"
int FindTopOrBottom(ViewI3D_H LayerID, int XLow, int XHigh, int YLow, int YHigh, int nz, int L, std::string HighLow) {

    int TopBottomZ = -1;
    bool SearchingForTop = true;
    int i = XLow;
    int j = YLow;
    int k = nz - 2;
    if (HighLow == "Low")
        k = 3;
    while (SearchingForTop) {
        if (LayerID(k, i, j) != L) {
            j++;
            if (j > YHigh) {
                if (i > XHigh) {
                    // All X and Y coordinates at this Z have been checked, and this Z coordinate is not part of the
                    // last depositied layer
                    TopBottomZ = k;
                    SearchingForTop = false;
                }
                else {
                    i++;
                    j = YLow;
                }
            }
        }
        else {
            // Search next Z coordinate
            if (HighLow == "High")
                k--;
            else if (HighLow == "Low")
                k++;
            i = XLow;
            j = YLow;
            if ((k == 0) || (k == nz - 1))
                SearchingForTop = false;
        }
    }
    if ((HighLow == "High") && (k == 0))
        TopBottomZ = nz - 2;
    else if ((HighLow == "Low") && (k == nz - 1))
        TopBottomZ = 1;
    return TopBottomZ;
}

void ParseLogFile(std::string LogFile, int &nx, int &ny, int &nz, double &deltax, int &NumberOfLayers) {

    std::ifstream InputDataStream;
    InputDataStream.open(LogFile);
    if (!(InputDataStream))
        throw std::runtime_error("Error: Cannot find ExaCA log file");
    std::string line;
    // 2 unused lines at beginning of log file
    for (int DummyLines = 0; DummyLines < 2; DummyLines++) {
        getline(InputDataStream, line);
    }
    std::string SimulationType = parseInput(InputDataStream, "simulation was type");

    // Line containing x, y, and z dimensions of domain
    std::string nxString = parseInput(InputDataStream, "Domain size in x");
    std::string nyString = parseInput(InputDataStream, "Domain size in y");
    std::string nzString = parseInput(InputDataStream, "Domain size in z");
    std::string deltaxString = parseInput(InputDataStream, "Cell size");
    // Convert strings to integers
    nx = std::stoi(nxString);
    ny = std::stoi(nyString);
    nz = std::stoi(nzString);
    deltax = atof(deltaxString.c_str());
    std::cout << "Dimensions of data are: nx = " << nx << ", ny = " << ny << ", nz = " << nz << std::endl;
    std::cout << "Cell size is " << deltax << " microns" << std::endl;
    if (SimulationType == "C")
        NumberOfLayers = 1;
    else {
        // Keep reading to determine if this is a multilayer simulation of not
        bool FindingNumLayers = true;
        while (FindingNumLayers) {
            getline(InputDataStream, line);
            std::size_t NLayersFound = line.find("layers");
            if (NLayersFound != std::string::npos) {
                std::string NumLayersString = line.substr(0, NLayersFound);
                NumberOfLayers = stoi(NumLayersString, nullptr, 10);
                FindingNumLayers = false;
            }
        }
    }
    std::cout << "Number of layers in microstructure data: " << NumberOfLayers << std::endl;
    InputDataStream.close();
}

// Reads portion of a paraview file and places data in the appropriate data structure
void ReadField(std::ifstream &InputDataStream, int nx, int ny, int nz, ViewI3D_H FieldOfInterest) {
    std::string line;
    for (int k = 0; k < nz; k++) {
        getline(InputDataStream, line);
        // Used to split string around spaces
        std::istringstream ss(line);
        std::string GIDVal; // for storing each separated string
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int val;
                ss >> val;
                FieldOfInterest(k, i, j) = val;
            }
        }
    }
}

// Read the analysis file to determine the file names/paths for the microstructure and the orientations
void ParseFilenames(std::string AnalysisFile, std::string &LogFile, std::string &MicrostructureFile,
                    std::string &RotationFilename, std::string &OutputFileName, std::string &EulerAnglesFilename,
                    bool &NewOrientationFormatYN) {

    // The analysis file should be in the analysis subdirectory of ExaCA
    std::cout << "Looking for analysis file: " << AnalysisFile << std::endl;
    std::ifstream Analysis;
    Analysis.open(AnalysisFile);
    if (!(Analysis))
        throw std::runtime_error("Error: Cannot find ExaCA analysis file");
    skipLines(Analysis, "*****");

    std::vector<std::string> NamesOfFiles = {"name of log (.log) file", "name of microstructure (.vtk) file",
                                             "file of grain orientations (rotation matrix form)",
                                             "name of data files of output resulting from this analysis",
                                             "file of grain orientations (Bunge Euler angle form ZXZ)"};
    std::vector<std::string> FilesRead(5);

    for (int i = 0; i < 5; i++) {
        std::string line;
        std::getline(Analysis, line);
        if (line == "**********")
            break;
        bool LineFound = parseInputFromList(line, NamesOfFiles, FilesRead, 5);
        if (!(LineFound))
            std::cout << "Warning: the line " << line
                      << " did not match any analysis file options expected by ExaCA and will be ignored" << std::endl;
    }
    // Check that all required arguments were present
    for (int i = 0; i < 4; i++) {
        if (FilesRead[i].empty()) {
            std::string error = "Error: Required input " + NamesOfFiles[i] + " not found in analysis file";
            throw std::runtime_error(error);
        }
    }
    LogFile = FilesRead[0];
    MicrostructureFile = FilesRead[1];
    RotationFilename = FilesRead[2];
    OutputFileName = FilesRead[3];

    // If the file of bunge convention euler angles was given, initialize this and set NewOrientationFormatYN to true;
    // otherwise, set NewOrientationFormatYN to false
    if (FilesRead[4].empty()) {
        NewOrientationFormatYN = false;
    }
    else {
        NewOrientationFormatYN = true;
        EulerAnglesFilename = FilesRead[4];
    }
    Analysis.close();
    // Path to file of grain orientations based on install/source location
    RotationFilename = checkFileInstalled(RotationFilename, "Substrate", 0);
    // Check that files are not empty
    checkFileNotEmpty(LogFile);
    checkFileNotEmpty(MicrostructureFile);
    checkFileNotEmpty(RotationFilename);
    // Same checks for EulerAnglesFilename, if given
    if (NewOrientationFormatYN) {
        EulerAnglesFilename = checkFileInstalled(EulerAnglesFilename, "Substrate", 0);
        checkFileNotEmpty(EulerAnglesFilename);
    }
}

void InitializeData(std::string MicrostructureFile, int nx, int ny, int nz, ViewI3D_H GrainID, ViewI3D_H LayerID,
                    ViewI3D_H Melted) {

    std::ifstream InputDataStream;
    InputDataStream.open(MicrostructureFile);
    if (!(InputDataStream))
        throw std::runtime_error("Error: Cannot find ExaCA microstructure file");
    std::string line;
    // 8 unused lines at beginning of Paraview file from header
    for (int DummyLines = 0; DummyLines < 8; DummyLines++) {
        getline(InputDataStream, line);
    }
    for (int field = 0; field < 3; field++) {
        // This line says which variable appears next in the file
        getline(InputDataStream, line);
        size_t FoundGrainID, FoundLayerID, FoundMelted;
        FoundGrainID = line.find("GrainID");
        FoundLayerID = line.find("LayerID");
        FoundMelted = line.find("Melted");
        // 1 more unused line
        getline(InputDataStream, line);
        // Place appropriate data
        if (FoundGrainID != std::string::npos)
            ReadField(InputDataStream, nx, ny, nz, GrainID);
        else if (FoundLayerID != std::string::npos)
            ReadField(InputDataStream, nx, ny, nz, LayerID);
        else if (FoundMelted != std::string::npos)
            ReadField(InputDataStream, nx, ny, nz, Melted);
        else
            throw std::runtime_error(
                "Error: unexpected data field in ExaCA microstructure file (not GrainID, LayerID, nor Melted");
        std::cout << "Data field " << field + 1 << " of 3 read" << std::endl;
    }
    InputDataStream.close();
}

// Read the analysis file to determine which analysis operations will be performed on the given ExaCA data
void ParseAnalysisFile(std::string AnalysisFile, std::string RotationFilename, int &NumberOfOrientations,
                       bool *AnalysisTypes, std::vector<int> &XLow_RVE, std::vector<int> &XHigh_RVE,
                       std::vector<int> &YLow_RVE, std::vector<int> &YHigh_RVE, std::vector<int> &ZLow_RVE,
                       std::vector<int> &ZHigh_RVE, int &NumberOfRVEs, std::vector<int> &CrossSectionPlane,
                       std::vector<int> &CrossSectionLocation, int &NumberOfCrossSections, int &XMin, int &XMax,
                       int &YMin, int &YMax, int &ZMin, int &ZMax, int nx, int ny, int nz, ViewI3D_H LayerID, ViewI3D_H,
                       int NumberOfLayers, std::vector<bool> &PrintSectionPF, std::vector<bool> &PrintSectionIPF,
                       bool NewOrientationFormatYN) {

    int FullDomainCenterX = nx / 2;
    int FullDomainCenterY = ny / 2;
    int FullDomainCenterZ = nz / 2;
    std::ifstream Analysis;
    Analysis.open(AnalysisFile);
    // Skip lines that were already read as part of ParseFilenames
    skipLines(Analysis, "**********");

    // Get the number of orientations from the file of rotation matrices
    std::ifstream GrainO;
    GrainO.open(RotationFilename);
    if (!(GrainO))
        throw std::runtime_error("Error: Cannot find ExaCA rotations file");
    std::string NumberOfOrientationsString;
    getline(GrainO, NumberOfOrientationsString);
    std::cout << "Number of orientations " << NumberOfOrientationsString << std::endl;
    NumberOfOrientations = stoi(NumberOfOrientationsString, nullptr, 10);
    GrainO.close();
    skipLines(Analysis, "*****");

    // Read some number of lines of RVE data
    bool ReadingRVEs = true;
    while (ReadingRVEs) {
        std::string RVEline;
        getline(Analysis, RVEline);
        if (RVEline == "**********")
            break;
        std::string RVE[4];
        // Should be 4 inputs on this line - break it into its components
        size_t compstart = 0, compend = 0;
        for (int comp = 0; comp < 4; comp++) {
            if (comp == 0)
                compstart = RVEline.find("[");
            else
                compstart = compend; // from last parsed string
            if (comp == 3)
                compend = RVEline.find("]");
            else
                compend = RVEline.find(",", compstart + 1);
            RVE[comp] = RVEline.substr(compstart + 1, compend - compstart - 1);
            // remove any whitespace from parsed string
            std::regex r("\\s+");
            RVE[comp] = std::regex_replace(RVE[comp], r, "");
        }
        // RVE[0] = Size, RVE[1] = Center X, RVE[2] = CenterY, RVE[3] = CenterZ
        // Default option for Center X or Center Y takes RVE X/Y center to be at the center of the data
        // Default option for Center Z (if single layer dataset) takes the RVE Z center to be at the center
        // Default option for Center Z (if multilayer dataset) takes Z center such that we avoid the last layer if
        // possible, but are as close to the top as possible RVE Size
        int Size_RVE = stoi(RVE[0], nullptr, 10);
        int HalfSize_RVE = Size_RVE / 2;
        // Center X
        int CenterX;
        if (RVE[1] == "D")
            CenterX = FullDomainCenterX;
        else
            CenterX = stoi(RVE[1], nullptr, 10);
        // Center Y
        int CenterY;
        if (RVE[2] == "D")
            CenterY = FullDomainCenterY;
        else
            CenterY = stoi(RVE[2], nullptr, 10);
        std::cout << "Center of RVE in X was selected to be at " << CenterX << std::endl;
        std::cout << "Center of RVE in Y was selected to be at " << CenterY << std::endl;
        // Obtain XLow/YLow/XHigh/YHigh for this RVE based on CenterX/CenterY/Size
        XLow_RVE.push_back(CenterX - HalfSize_RVE);
        YLow_RVE.push_back(CenterY - HalfSize_RVE);
        XHigh_RVE.push_back(Size_RVE + XLow_RVE[NumberOfRVEs] - 1);
        YHigh_RVE.push_back(Size_RVE + YLow_RVE[NumberOfRVEs] - 1);
        if (XLow_RVE[NumberOfRVEs] < 0) {
            XLow_RVE[NumberOfRVEs] = 0;
            std::cout << "WARNING: RVE number " << NumberOfRVEs
                      << " lower x bound adjusted to 0 to avoid out of bounds error" << std::endl;
        }
        if (YLow_RVE[NumberOfRVEs] < 0) {
            YLow_RVE[NumberOfRVEs] = 0;
            std::cout << "WARNING: RVE number " << NumberOfRVEs
                      << " lower y bound adjusted to 0 to avoid out of bounds error" << std::endl;
        }
        if (XHigh_RVE[NumberOfRVEs] > nx - 1) {
            XHigh_RVE[NumberOfRVEs] = nx - 1;
            std::cout << "WARNING: RVE number " << NumberOfRVEs << " upper x bound adjusted to " << nx - 1
                      << " to avoid out of bounds error" << std::endl;
        }
        if (YHigh_RVE[NumberOfRVEs] > ny - 1) {
            YHigh_RVE[NumberOfRVEs] = ny - 1;
            std::cout << "WARNING: RVE number " << NumberOfRVEs << " upper y bound adjusted to " << ny - 1
                      << " to avoid out of bounds error" << std::endl;
        }
        // Center Z
        if (RVE[3] != "D") {
            int CenterZ = stoi(RVE[3], nullptr, 10);
            ZLow_RVE.push_back(CenterZ - HalfSize_RVE);
            ZHigh_RVE.push_back(Size_RVE + ZLow_RVE[NumberOfRVEs] - 1);
        }
        else {
            int RVETop =
                FindTopOrBottom(LayerID, XLow_RVE[NumberOfRVEs], XHigh_RVE[NumberOfRVEs], YLow_RVE[NumberOfRVEs],
                                YHigh_RVE[NumberOfRVEs], nz, NumberOfLayers - 1, "High");
            int RVEBottom = RVETop - Size_RVE + 1;
            if (RVEBottom < 0) {
                // We can't avoid getting at least part of the top layer in the RVE, so default to taking the RVE from
                // the center of the domain in z
                ZLow_RVE.push_back(FullDomainCenterZ - HalfSize_RVE);
                ZHigh_RVE.push_back(Size_RVE + ZLow_RVE[NumberOfRVEs] - 1);
            }
            else {
                ZLow_RVE.push_back(RVEBottom);
                ZHigh_RVE.push_back(RVETop);
            }
        }
        std::cout << "RVE number " << NumberOfRVEs << " has X bounds " << XLow_RVE[NumberOfRVEs] << "/"
                  << XHigh_RVE[NumberOfRVEs] << "; Y bounds " << YLow_RVE[NumberOfRVEs] << "/"
                  << YHigh_RVE[NumberOfRVEs] << "; Z bounds " << ZLow_RVE[NumberOfRVEs] << "/"
                  << ZHigh_RVE[NumberOfRVEs] << std::endl;
        NumberOfRVEs++;
    }
    skipLines(Analysis, "*****");
    // Read some number of lines of cross-section data
    bool ReadingCrossSections = true;
    while (ReadingCrossSections) {
        std::string CSline;
        getline(Analysis, CSline);
        if (CSline == "**********")
            break;

        // Check for "old" format ([section type, coordinate out of plane]
        // vs "new" format (section type, coordinate out of plane, print pole figure data Y/N, print IPF-colored
        // cross-section data Y/N)
        std::vector<std::string> CS(4);
        if (NewOrientationFormatYN) {
            // Should be 4 inputs on this line - break it into its components
            splitString(CSline, CS, 4);
            if (CS[2] == "Y")
                PrintSectionPF.push_back(true);
            else
                PrintSectionPF.push_back(false);
            if (CS[3] == "Y")
                PrintSectionIPF.push_back(true);
            else
                PrintSectionIPF.push_back(false);
        }
        else {
            // Should be 2 inputs on this line - break it into its components
            CS[0] = parseCoordinatePair(CSline, 0);
            CS[1] = parseCoordinatePair(CSline, 1);
            PrintSectionPF.push_back(false);
            PrintSectionIPF.push_back(true);
        }
        if (CS[0] == "XZ") {
            CrossSectionPlane.push_back(0);
            if (CS[1] == "END")
                CrossSectionLocation.push_back(ny - 2);
            else if (CS[1] == "D")
                CrossSectionLocation.push_back(FullDomainCenterY);
            else {
                CrossSectionLocation.push_back(stoi(CS[1], nullptr, 10));
                if (CrossSectionLocation[NumberOfCrossSections] >= ny) {
                    CrossSectionLocation[NumberOfCrossSections] = ny - 1;
                    std::cout << "Adjusting cross-section number " << NumberOfCrossSections + 1
                              << " y value, which was larger than the domain size in y" << std::endl;
                }
            }
        }
        else if (CS[0] == "YZ") {
            CrossSectionPlane.push_back(1);
            if (CS[1] == "END")
                CrossSectionLocation.push_back(nx - 2);
            else if (CS[1] == "D")
                CrossSectionLocation.push_back(FullDomainCenterX);
            else {
                CrossSectionLocation.push_back(stoi(CS[1], nullptr, 10));
                if (CrossSectionLocation[NumberOfCrossSections] >= nx) {
                    CrossSectionLocation[NumberOfCrossSections] = nx - 1;
                    std::cout << "Adjusting .ang (pyEBSD) cross-section number " << NumberOfCrossSections + 1
                              << " x value, which was larger than the domain size in x" << std::endl;
                }
            }
        }
        else {
            CrossSectionPlane.push_back(2);
            if (CS[1] == "END")
                CrossSectionLocation.push_back(nz - 2);
            else if (CS[1] == "D")
                CrossSectionLocation.push_back(FullDomainCenterZ);
            else {
                CrossSectionLocation.push_back(stoi(CS[1], nullptr, 10));
                if (CrossSectionLocation[NumberOfCrossSections] >= nz) {
                    CrossSectionLocation[NumberOfCrossSections] = nz - 1;
                    std::cout << "Adjusting .ang (pyEBSD) cross-section number " << NumberOfCrossSections + 1
                              << " z value, which was larger than the domain size in z" << std::endl;
                }
            }
        }
        std::cout << "Cross-section number " << NumberOfCrossSections << " has parsed string data "
                  << CrossSectionPlane[NumberOfCrossSections] << " " << CrossSectionLocation[NumberOfCrossSections]
                  << std::endl;
        NumberOfCrossSections++;
    }
    skipLines(Analysis, "*****");

    // Should be 2 inputs on this line (x lower and upper bounds)- break it into its components
    std::string XYZCoordLine;
    getline(Analysis, XYZCoordLine);
    std::string XMinString = parseCoordinatePair(XYZCoordLine, 0);
    std::string XMaxString = parseCoordinatePair(XYZCoordLine, 1);
    // Should be 2 inputs on this line (y lower and upper bounds)- break it into its components
    getline(Analysis, XYZCoordLine);
    std::string YMinString = parseCoordinatePair(XYZCoordLine, 0);
    std::string YMaxString = parseCoordinatePair(XYZCoordLine, 1);
    // Should be 2 inputs on this line (z lower and upper bounds)- break it into its components
    getline(Analysis, XYZCoordLine);
    std::string ZMinString = parseCoordinatePair(XYZCoordLine, 0);
    std::string ZMaxString = parseCoordinatePair(XYZCoordLine, 1);

    if (XMinString == "D")
        XMin = 50;
    else
        XMin = stoi(XMinString, nullptr, 10);

    if (XMaxString == "D")
        XMax = nx - 50;
    else
        XMax = stoi(XMaxString, nullptr, 10);

    if (YMinString == "D")
        YMin = 50;
    else
        YMin = stoi(YMinString, nullptr, 10);

    if (YMaxString == "D")
        YMax = ny - 50;
    else
        YMax = stoi(YMaxString, nullptr, 10);
    if (ZMinString == "D") {
        ZMin = FindTopOrBottom(LayerID, XMin, XMax, YMin, YMax, nz, 0, "Low");
    }
    else
        ZMin = stoi(ZMinString, nullptr, 10);
    if (ZMaxString == "D") {
        ZMax = FindTopOrBottom(LayerID, XMin, XMax, YMin, YMax, nz, NumberOfLayers - 1, "High");
    }
    else
        ZMax = stoi(ZMaxString, nullptr, 10);

    std::cout << "XYZ region for grain statistics is X = " << XMin << " through " << XMax << "; Y = " << YMin
              << " through " << YMax << "; Z = " << ZMin << " through " << ZMax << std::endl;
    // Read Y/N options for analyzing microstructure data
    std::vector<std::string> AnalysisOptions = {
        "misorientation",           "grain volume",      "mean aspect ratio",        "mean grain area",
        "mean weighted grain area", "mean grain height", "grain width distribution", "volume for pole figure analysis"};
    int NumAnalysisOptions = AnalysisOptions.size();
    std::vector<std::string> AnalysisOptionsYN(NumAnalysisOptions);

    // Read files to get Y/N values for possible analysis options from the list "AnalysisOptions", storing results in
    // AnalysisOptionsYN
    for (int i = 0; i < NumAnalysisOptions; i++) {
        std::string line;
        std::getline(Analysis, line);
        bool AnalysisOptionFound = parseInputFromList(line, AnalysisOptions, AnalysisOptionsYN, NumAnalysisOptions);
        if (!(AnalysisOptionFound))
            std::cout << "Warning: unknown analysis option " << line << std::endl;
    }
    // Check which analysis options were found in the file and store their results as bools
    for (int i = 0; i < NumAnalysisOptions; i++) {
        if (!(AnalysisOptionsYN[i].empty()))
            AnalysisTypes[i] = getInputBool(AnalysisOptionsYN[i]);
        else
            AnalysisTypes[i] = false;
    }
    Analysis.close();
}
