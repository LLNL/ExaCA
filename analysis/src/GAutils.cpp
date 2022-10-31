// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
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

void ParseLogFile(std::string LogFile, int &nx, int &ny, int &nz, double &deltax, int &NumberOfLayers,
                  bool UseXYZBounds, std::vector<double> &XYZBounds) {

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

    if (UseXYZBounds) {
        // Time step not used by analysis script
        getline(InputDataStream, line);

        // Get domain bounds
        std::vector<std::string> XYZBoundsString(6);

        // Check that the bounds are present in the file
        // If this line does not contain the lower bound of the domain in x, throw an error that the input data set is
        // not compatible with this executable
        std::string testline;
        std::string expectedline = "Lower bound of domain in x:";
        std::getline(InputDataStream, testline);
        if (testline.find(expectedline) == std::string::npos)
            throw std::runtime_error("Error: input data set not compatible with amb analysis executable: domain bounds "
                                     "are not listed in the log file");
        else {
            // Parse this line as normal
            std::size_t colon = testline.find(":");
            XYZBoundsString[0] = testline.substr(colon + 1, std::string::npos);
            XYZBoundsString[0] = removeWhitespace(XYZBoundsString[0]);
        }
        // Parse remaining bounds
        XYZBoundsString[1] = parseInput(InputDataStream, "Lower bound of domain in y");
        XYZBoundsString[2] = parseInput(InputDataStream, "Lower bound of domain in z");
        XYZBoundsString[3] = parseInput(InputDataStream, "Upper bound of domain in x");
        XYZBoundsString[4] = parseInput(InputDataStream, "Upper bound of domain in y");
        XYZBoundsString[5] = parseInput(InputDataStream, "Upper bound of domain in z");

        // If log file contained the X, Y, Z bounds of the domain and they're used in the analysis, convert them
        for (int n = 0; n < 6; n++)
            XYZBounds[n] = getInputDouble(XYZBoundsString[n]);
        std::cout << "X bounds of data are " << XYZBounds[0] << ", " << XYZBounds[3] << "; Y bounds of data are "
                  << XYZBounds[1] << ", " << XYZBounds[4] << "; Z bounds of data are " << XYZBounds[2] << ", "
                  << XYZBounds[5] << std::endl;
    }

    // Get number of layers
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

// Reads and ignored field data
void ReadIgnoreField(std::ifstream &InputDataStream, int nx, int ny, int nz) {
    std::string line;
    for (int i = 0; i < nx * ny * nz; i++)
        getline(InputDataStream, line);
}

// Read the analysis file to determine the file names/paths for the microstructure and the orientations
void ParseFilenames(std::string AnalysisFile, std::string &LogFile, std::string &MicrostructureFile,
                    std::string &RotationFilename, std::string &OutputFileName, std::string &EulerAnglesFilename,
                    bool &NewOrientationFormatYN, std::string &RGBFilename) {

    // The analysis file should be in the analysis subdirectory of ExaCA
    std::cout << "Looking for analysis file: " << AnalysisFile << std::endl;
    std::ifstream Analysis;
    Analysis.open(AnalysisFile);
    if (!(Analysis))
        throw std::runtime_error("Error: Cannot find ExaCA analysis file");
    skipLines(Analysis, "*****");

    std::vector<std::string> NamesOfFiles = {
        "name of log (.log) file",
        "name of microstructure (.vtk) file",
        "file of grain orientations (rotation matrix form)",
        "name of data files of output resulting from this analysis",
        "file of grain orientations (Bunge Euler angle form ZXZ)",
        "file corresponding to grain orientation IPF-Z colors as fractional RGB values"};
    std::vector<std::string> FilesRead(6);

    for (int i = 0; i < 6; i++) {
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
        std::cout << "Warning: the name of the file of Euler angles was not provided in the analysis input file, but "
                     "will be required in a future release"
                  << std::endl;
    }
    else {
        NewOrientationFormatYN = true;
        EulerAnglesFilename = FilesRead[4];
    }
    Analysis.close();
    // Path to file of grain orientations based on install/source location
    RotationFilename = checkFileInstalled(RotationFilename, 0);
    // Check that files are not empty
    checkFileNotEmpty(LogFile);
    checkFileNotEmpty(MicrostructureFile);
    checkFileNotEmpty(RotationFilename);
    // Same checks for EulerAnglesFilename, if given
    if (NewOrientationFormatYN) {
        EulerAnglesFilename = checkFileInstalled(EulerAnglesFilename, 0);
        checkFileNotEmpty(EulerAnglesFilename);
    }
    // RGB file - default value or value from file
    if (FilesRead[5].empty()) {
        std::cout << "Warning: the name of the file of RGB colors for each orientation corresponding to the IPF-Z "
                     "color was not required, but will be required in a future release. The default file "
                     "GrainOrientationRGB_IPF-Z.csv will be used to map colors to orientations"
                  << std::endl;
        RGBFilename = "GrainOrientationRGB_IPF-Z.csv";
    }
    else
        RGBFilename = FilesRead[5];
    // Check that file exists and was successfully installed
    RGBFilename = checkFileInstalled(RGBFilename, 0);
}

// Ensure that the appropriate files exist - orientation files should be installed, other files should start with the
// BaseFileName value given from the command line
void CheckInputFiles_AMB(std::string BaseFileName, std::string &LogFile, std::string &MicrostructureFile,
                         std::string &RotationFilename, std::string &EulerAnglesFilename, std::string &RGBFilename) {

    // Names of files
    LogFile = BaseFileName + ".log";
    MicrostructureFile = BaseFileName + ".vtk";
    RotationFilename = "GrainOrientationVectors.csv";
    RotationFilename = checkFileInstalled(RotationFilename, 0);
    EulerAnglesFilename = "GrainOrientationEulerAnglesBungeZXZ.csv";
    EulerAnglesFilename = checkFileInstalled(EulerAnglesFilename, 0);
    RGBFilename = "GrainOrientationRGB_IPF-Z.csv";
    RGBFilename = checkFileInstalled(RGBFilename, 0);

    // Check that files are not empty
    checkFileNotEmpty(LogFile);
    checkFileNotEmpty(MicrostructureFile);
    checkFileNotEmpty(RotationFilename);
    checkFileNotEmpty(EulerAnglesFilename);
    checkFileNotEmpty(RGBFilename);
}
void InitializeData(std::string MicrostructureFile, int nx, int ny, int nz, ViewI3D_H GrainID, ViewI3D_H LayerID) {

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
        // A blank line ends the file read
        getline(InputDataStream, line);
        if (line.empty())
            break;
        size_t FoundGrainID = line.find("GrainID");
        size_t FoundLayerID = line.find("LayerID");
        size_t FoundMelted = line.find("Melted");
        // 1 more unused line
        getline(InputDataStream, line);
        // Place appropriate data
        if (FoundGrainID != std::string::npos)
            ReadField(InputDataStream, nx, ny, nz, GrainID);
        else if (FoundLayerID != std::string::npos)
            ReadField(InputDataStream, nx, ny, nz, LayerID);
        else if (FoundMelted != std::string::npos) {
            std::cout << "Note: Melted data is no longer used and will be ignored" << std::endl;
            ReadIgnoreField(InputDataStream, nx, ny, nz);
        }
        else
            throw std::runtime_error(
                "Error: unexpected data field in ExaCA microstructure file (not GrainID or LayerID)");
        std::string Fieldname;
        if (FoundGrainID != std::string::npos)
            Fieldname = "GrainID";
        else if (FoundLayerID != std::string::npos)
            Fieldname = "LayerID";
        else if (FoundMelted != std::string::npos)
            Fieldname = "Melted";
        std::cout << "Data field " << Fieldname << " read" << std::endl;
    }
    InputDataStream.close();
}

// Read the analysis file to determine which analysis operations will be performed on the given ExaCA data
void ParseAnalysisFile(std::string AnalysisFile, std::string RotationFilename, int &NumberOfOrientations,
                       bool *AnalysisTypes, std::vector<int> &XLow_RVE, std::vector<int> &XHigh_RVE,
                       std::vector<int> &YLow_RVE, std::vector<int> &YHigh_RVE, std::vector<int> &ZLow_RVE,
                       std::vector<int> &ZHigh_RVE, int &NumberOfRVEs, std::vector<std::string> &CrossSectionPlane,
                       std::vector<int> &CrossSectionLocation, int &NumberOfCrossSections, int &XMin, int &XMax,
                       int &YMin, int &YMax, int &ZMin, int &ZMax, int nx, int ny, int nz, ViewI3D_H LayerID,
                       int NumberOfLayers, std::vector<bool> &PrintSectionPF, std::vector<bool> &PrintSectionIPF,
                       std::vector<bool> &BimodalAnalysis, bool NewOrientationFormatYN,
                       std::vector<std::string> &CSLabels) {

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
        // For the "new" format, there are either 4 or 5 comma-separated components on this line, depending on whether
        // the option of analyzing small and large area grains separately is given. If it is not given, a single
        // distribution of grain areas in the cross-section will be analyzed
        int NumComponents;
        if (NewOrientationFormatYN)
            NumComponents = std::count(CSline.begin(), CSline.end(), ',') + 1;
        else
            NumComponents = 2;
        std::vector<std::string> CS(NumComponents);
        if (NewOrientationFormatYN) {
            // Should be 4 or 5 inputs on this line - break it into its components
            splitString(CSline, CS, NumComponents);
            // Remove whitespace from parsed cross-section data
            for (int n = 0; n < NumComponents; n++) {
                CS[n] = removeWhitespace(CS[n]);
            }
            if (CS[2] == "Y")
                PrintSectionPF.push_back(true);
            else
                PrintSectionPF.push_back(false);
            if (CS[3] == "Y")
                PrintSectionIPF.push_back(true);
            else
                PrintSectionIPF.push_back(false);
            if (NumComponents == 4)
                BimodalAnalysis.push_back(false);
            else {
                if (CS[4] == "Y")
                    BimodalAnalysis.push_back(true);
                else
                    BimodalAnalysis.push_back(false);
            }
        }
        else {
            // Should be 2 inputs on this line - break it into its components
            CS[0] = parseCoordinatePair(CSline, 0);
            CS[1] = parseCoordinatePair(CSline, 1);
            PrintSectionPF.push_back(false);
            PrintSectionIPF.push_back(true);
            BimodalAnalysis.push_back(false);
            std::cout << "Warning: old format (2 comma-separated values per line) for Section 2 of analysis input file "
                         "detected, the new format will be required in a future release"
                      << std::endl;
        }
        if (CS[0] == "XZ") {
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
        std::string ThisCrossSectionPlane = CS[0] + std::to_string(CrossSectionLocation[NumberOfCrossSections]);
        CrossSectionPlane.push_back(ThisCrossSectionPlane);
        std::cout << "Cross-section number " << NumberOfCrossSections << " has parsed string data "
                  << CrossSectionPlane[NumberOfCrossSections] << " " << CrossSectionLocation[NumberOfCrossSections]
                  << std::endl;
        std::string CrossSectionLabel = CS[0] + " cross-section located at out-of-plane coordinate " +
                                        std::to_string(CrossSectionLocation[NumberOfCrossSections]);
        CSLabels.push_back(CrossSectionLabel);
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
    if (!(AnalysisOptionsYN[7].empty()))
        std::cout << "Note: ability to print pole figure data for the specified volume will be deprecated in a future "
                     "release, option will be available for cross-sections only"
                  << std::endl;
    Analysis.close();
}

// Given an input vector of integer Grain ID values, return an output vector consisting of the unique Grain ID values,
// sorted from lowest to highest
std::vector<int> FindUniqueGrains(const std::vector<int> GrainIDVector) {
    std::vector<int> UniqueGrainIDVector = GrainIDVector;
    std::sort(UniqueGrainIDVector.begin(), UniqueGrainIDVector.end());
    std::vector<int>::iterator it;
    it = std::unique(UniqueGrainIDVector.begin(), UniqueGrainIDVector.end());
    UniqueGrainIDVector.resize(std::distance(UniqueGrainIDVector.begin(), it));
    return UniqueGrainIDVector;
}
