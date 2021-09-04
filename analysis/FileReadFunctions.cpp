// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

// These are duplicated from CAinitialize.cpp - will not need if the analysis executable is
// linked properly to the ExaCA-Kokkos executable
#include "header.h"
#include <stdexcept>
#include <string>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <regex>

// Skip initial lines in input files.
void skipLines(std::ifstream &stream) {
    std::string line;
    while (getline(stream, line)) {
        if (line == "*****")
            break;
    }
}

// Find the colon on a line read from a file, and return the value before the colon
// Throw an error if no colon is found, get a new line and throw a warning if the value before
// the colon is a deprecated input
// Also store the location of the colon within the line
std::string getKey(std::ifstream &stream, std::string &line, std::size_t &colon) {

    bool DeprecatedInputCheck = true;
    std::string actual_key;
    while (DeprecatedInputCheck) {
        std::getline(stream, line);
        colon = line.find(":");
        actual_key = line.substr(0, colon);
        // Check for colon seperator
        if (colon == std::string::npos) {
            std::string error = "Input \"" + actual_key + "\" must be separated from value by \":\".";
            throw std::runtime_error(error);
        }
        std::vector<std::string> deprecated_inputs = {"Burst buffer", "Source of input length unit"};
        for (auto di : deprecated_inputs)
            if (actual_key.find(di) != std::string::npos) {
                std::cout << "WARNING - this input has been deprecated and has no effect, \"" << actual_key << "\""
                          << std::endl;
                // Ignore this line and get another line
            }
            else {
                // Input was not a deprecated line - continue
                DeprecatedInputCheck = false;
            }
    }
    return actual_key;
}

// Remove whitespace from "line", taking only the portion of the line that comes after the colon
std::string removeWhitespace(std::string line, std::size_t colon) {

    std::string val = line.substr(colon + 1, std::string::npos);
    std::regex r("\\s+");
    val = std::regex_replace(val, r, "");
    return val;
}

// Verify the required input was included with the correct format, and parse the input (with only 1 possibility for the
// possible input)
std::string parseInput(std::ifstream &stream, std::string key) {

    std::size_t colon;
    std::string line, val;
    std::string actual_key = getKey(stream, line, colon);

    // Check for keyword
    if (actual_key.find(key) == std::string::npos) {
        // Keyword not found
        std::string error = "Required input not present: " + key + " not found in the input file";
        throw std::runtime_error(error);
    }
    // Keyword was found
    val = removeWhitespace(line, colon);
    return val;
}

// Verify the required boolean input was included with the correct format.
bool parseInputBool(std::ifstream &stream, std::string key) {
    std::string val = parseInput(stream, key);
    if (val == "N") {
        return false;
    }
    else if (val == "Y") {
        return true;
    }
    else {
        std::string error = "Input \"" + key + "\" must be \"Y\" or \"N\".";
        throw std::runtime_error(error);
    }
}

// Parse a line that looks like [x,y], returning either val = 0 (x) or val = 1 (y)
std::string parseCoordinatePair(std::string line, int val) {
    size_t linestart = line.find("[");
    size_t linebreak = line.find(",");
    size_t lineend = line.find("]");
    std::string SplitStr;
    if (val == 0) SplitStr = line.substr(linestart + 1, linebreak - linestart - 1);
    else if (val == 1) SplitStr = line.substr(linebreak + 1, lineend - linebreak - 1);
    std::regex r("\\s+");
    SplitStr = std::regex_replace(SplitStr, r, "");
    return SplitStr;
}

// Search a specified region of LayerID in x and y for either the smallest Z that doesn't contain any layer "L",
// or the largest Z that doesn't contain any layer "L"
int FindTopOrBottom(int*** LayerID, int XLow, int XHigh, int YLow, int YHigh, int nz, int L, std::string HighLow) {
    
    int TopBottomZ;
    bool SearchingForTop = true;
    int i = XLow;
    int j = YLow;
    int k;
    if (HighLow == "High") k = nz-2;
    else if (HighLow == "Low") k = 3;
    while (SearchingForTop) {
        if (LayerID[k][i][j] != L) {
            j++;
            if (j > YHigh) {
                if (i > XHigh) {
                    // All X and Y coordinates at this Z have been checked, and this Z coordinate is not part of the last depositied layer
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
            if (HighLow == "High") k--;
            else if (HighLow == "Low") k++;
            i = XLow;
            j = YLow;
            if ((k == 0)||(k == nz-1)) SearchingForTop = false;
        }
    }
    if ((HighLow == "High")&&(k == 0)) TopBottomZ = nz-2;
    else if ((HighLow == "Low")&&(k == nz-1)) TopBottomZ = 1;
    return TopBottomZ;
}

void ParseLogFile(std::string BaseFileName, int &nx, int &ny, int &nz, double &deltax, int &NumberOfLayers) {
    
    std::string LogFile = BaseFileName + ".log";
    std::cout << "Looking for ExaCA log file with name " << LogFile << std::endl;
    std::ifstream InputDataStream;
    InputDataStream.open(LogFile);
    if (!(InputDataStream)) throw std::runtime_error("Error: Cannot find ExaCA log file");
    std::string line;
    // 2 unused lines at beginning of log file
    for (int DummyLines=0; DummyLines<2; DummyLines++) {
        getline(InputDataStream, line);
    }
    std::string SimulationType = parseInput(InputDataStream,"simulation was type");

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
    if (SimulationType == "C") NumberOfLayers = 1;
    else {
        // Keep reading to determine if this is a multilayer simulation of not
        bool FindingNumLayers = true;
        while (FindingNumLayers) {
            getline(InputDataStream, line);
            std::size_t NLayersFound = line.find("layers");
            if (NLayersFound != std::string::npos) {
                std::size_t ForFound = line.find("for");
                std::string NumLayersString = line.substr(ForFound+3, NLayersFound-(ForFound+3));
                NumberOfLayers = stoi(NumLayersString, nullptr, 10);
                FindingNumLayers = false;
            }
        }
    }
    std::cout << "Number of layers in microstructure data: " << NumberOfLayers << std::endl;
    InputDataStream.close();
    
}

// Reads portion of a paraview file and places data in the appropriate data structure
void ReadField(std::ifstream &InputDataStream, int nx, int ny, int nz, int*** FieldOfInterest) {
    std::string line;
    for (int k=0; k<nz; k++) {
        getline(InputDataStream, line);
        // Used to split string around spaces
        std::istringstream ss(line);
        std::string GIDVal; // for storing each separated string
        for (int j=0; j<ny; j++) {
            for (int i=0; i<nx; i++) {
                int val;
                ss >> val;
                FieldOfInterest[k][i][j] = val;
            }
        }
    }
    
}

void InitializeData(std::string InputFile, int nx, int ny, int nz, int*** GrainID, int*** LayerID, int*** Melted) {
    
    std::ifstream InputDataStream;
    InputDataStream.open(InputFile);
    if (!(InputDataStream)) throw std::runtime_error("Error: Cannot find ExaCA microstructure file");
    std::string line;
    // 8 unused lines at beginning of Paraview file from header
    for (int DummyLines=0; DummyLines<8; DummyLines++) {
        getline(InputDataStream, line);
    }
    for (int field=0; field<3; field++) {
    // This line says which variable appears next in the file
        getline(InputDataStream, line);
        size_t FoundGrainID, FoundLayerID, FoundMelted;
        FoundGrainID = line.find("GrainID");
        FoundLayerID = line.find("LayerID");
        FoundMelted = line.find("Melted");
        // 1 more unused line
        getline(InputDataStream, line);
        // Place appropriate data
        if (FoundGrainID != std::string::npos) ReadField(InputDataStream, nx, ny, nz, GrainID);
        else if (FoundLayerID != std::string::npos) ReadField(InputDataStream, nx, ny, nz, LayerID);
        else if (FoundMelted != std::string::npos) ReadField(InputDataStream, nx, ny, nz, Melted);
        else throw std::runtime_error("Error: unexpected data field in ExaCA microstructure file (not GrainID, LayerID, nor Melted");
        std::cout << "Data field " << field+1 << " of 3 read" << std::endl;
    }
    InputDataStream.close();

}

// Read the analysis file to determine which analysis operations will be performed on the given ExaCA data
void ParseAnalysisFile(std::string OutputAnalysisFile, std::string &RotationFilename, std::string &EulerFilename, int &NumberOfOrientations, bool* AnalysisTypes, std::vector <int> &XLow_RVE, std::vector <int> &XHigh_RVE, std::vector <int> &YLow_RVE, std::vector <int> &YHigh_RVE, std::vector <int> &ZLow_RVE, std::vector <int> &ZHigh_RVE, int &NumberOfRVEs, std::vector <int> &CrossSectionPlane, std::vector <int> &CrossSectionLocation, int &NumberOfCrossSections, int &XMin, int &XMax, int &YMin, int &YMax, int &ZMin, int &ZMax, int nx, int ny, int nz, int*** LayerID, int*** Melted, int NumberOfLayers) {
    
    int FullDomainCenterX = nx/2;
    int FullDomainCenterY = ny/2;
    int FullDomainCenterZ = nz/2;
    std::ifstream Analysis;
    Analysis.open(OutputAnalysisFile);
    if (!(Analysis)) throw std::runtime_error("Error: Cannot find ExaCA analysis file");
    skipLines(Analysis);
    
    // Read names of files containing the grain orientations, in rotation matrix and Euler angle forms
    RotationFilename = parseInput(Analysis, "rotation matrix");
    std::ifstream GrainO;
    
    // Get the number of orientations from the file of rotation matrices
    GrainO.open(RotationFilename);
    std::string NumberOfOrientationsString;
    getline(GrainO,NumberOfOrientationsString);
    std::cout << "Number of orientations " << NumberOfOrientationsString << std::endl;
    NumberOfOrientations = stoi(NumberOfOrientationsString, nullptr, 10);
    GrainO.close();
    EulerFilename = parseInput(Analysis, "Euler angle");
    skipLines(Analysis);
    
    // Read some number of lines of RVE data
    bool ReadingRVEs = true;
    while (ReadingRVEs) {
        std::string RVEline;
        getline(Analysis, RVEline);
        if (RVEline == "**********") break;
        std::string RVE[4];
        // Should be 4 inputs on this line - break it into its components
        size_t compstart, compend;
        for (int comp=0; comp<4; comp++) {
            if (comp == 0) compstart = RVEline.find("[");
            else compstart = compend; // from last parsed string
            if (comp == 3) compend = RVEline.find("]");
            else compend = RVEline.find(",", compstart + 1);
            RVE[comp] = RVEline.substr(compstart + 1, compend - compstart - 1);
            // remove any whitespace from parsed string
            std::regex r("\\s+");
            RVE[comp] = std::regex_replace(RVE[comp], r, "");
        }
        // RVE[0] = Size, RVE[1] = Center X, RVE[2] = CenterY, RVE[3] = CenterZ
        // Default option for Center X or Center Y takes RVE X/Y center to be at the center of the data
        // Default option for Center Z (if single layer dataset) takes the RVE Z center to be at the center
        // Default option for Center Z (if multilayer dataset) takes Z center such that we avoid the last layer if possible, but are as close to the top as possible
        // RVE Size
        int Size_RVE = stoi(RVE[0], nullptr, 10);
        int HalfSize_RVE = Size_RVE/2;
        // Center X
        int CenterX;
        if (RVE[1] == "D") CenterX = FullDomainCenterX;
        else CenterX = stoi(RVE[1], nullptr, 10);
        // Center Y
        int CenterY;
        if (RVE[2] == "D") CenterY = FullDomainCenterY;
        else CenterY = stoi(RVE[2], nullptr, 10);
        std::cout << "Center of RVE in X was selected to be at " << CenterX << std::endl;
        std::cout << "Center of RVE in Y was selected to be at " << CenterY << std::endl;
        // Obtain XLow/YLow/XHigh/YHigh for this RVE based on CenterX/CenterY/Size
        XLow_RVE.push_back(CenterX - HalfSize_RVE);
        YLow_RVE.push_back(CenterY - HalfSize_RVE);
        XHigh_RVE.push_back(Size_RVE + XLow_RVE[NumberOfRVEs] - 1);
        YHigh_RVE.push_back(Size_RVE + YLow_RVE[NumberOfRVEs] - 1);
        if (XLow_RVE[NumberOfRVEs] < 0) {
            XLow_RVE[NumberOfRVEs] = 0;
            std::cout << "WARNING: RVE number " << NumberOfRVEs << " lower x bound adjusted to 0 to avoid out of bounds error" << std::endl;
        }
        if (YLow_RVE[NumberOfRVEs] < 0) {
            YLow_RVE[NumberOfRVEs] = 0;
            std::cout << "WARNING: RVE number " << NumberOfRVEs << " lower y bound adjusted to 0 to avoid out of bounds error" << std::endl;
        }
        if (XHigh_RVE[NumberOfRVEs] > nx-1) {
            XHigh_RVE[NumberOfRVEs] = nx-1;
            std::cout << "WARNING: RVE number " << NumberOfRVEs << " upper x bound adjusted to " << nx-1 << " to avoid out of bounds error" << std::endl;
        }
        if (YHigh_RVE[NumberOfRVEs] > ny-1) {
            YHigh_RVE[NumberOfRVEs] = ny-1;
            std::cout << "WARNING: RVE number " << NumberOfRVEs << " upper y bound adjusted to " << ny-1 << " to avoid out of bounds error" << std::endl;
        }
        // Center Z
        if (RVE[3] != "D") {
            int CenterZ = stoi(RVE[2], nullptr, 10);
            ZLow_RVE.push_back(CenterZ - HalfSize_RVE);
            ZHigh_RVE.push_back(Size_RVE + ZLow_RVE[NumberOfRVEs] - 1);
        }
        else {
            int RVETop = FindTopOrBottom(LayerID, XLow_RVE[NumberOfRVEs], XHigh_RVE[NumberOfRVEs], YLow_RVE[NumberOfRVEs], YHigh_RVE[NumberOfRVEs], nz, NumberOfLayers-1, "High");
            int RVEBottom = RVETop - Size_RVE + 1;
            if (RVEBottom < 0) {
                // We can't avoid getting at least part of the top layer in the RVE, so default to taking the RVE from the center of the domain in z
                ZLow_RVE.push_back(FullDomainCenterZ - HalfSize_RVE);
                ZHigh_RVE.push_back(Size_RVE + ZLow_RVE[NumberOfRVEs] - 1);
            }
            else {
                ZLow_RVE.push_back(RVEBottom);
                ZHigh_RVE.push_back(RVETop);
            }
        }
        std::cout << "RVE number " << NumberOfRVEs << " has X bounds " << XLow_RVE[NumberOfRVEs] << "/" << XHigh_RVE[NumberOfRVEs] << "; Y bounds " << YLow_RVE[NumberOfRVEs] << "/" << YHigh_RVE[NumberOfRVEs] << "; Z bounds " << ZLow_RVE[NumberOfRVEs] << "/" << ZHigh_RVE[NumberOfRVEs] << std::endl;
        NumberOfRVEs++;
    }
    skipLines(Analysis);
    // Read some number of lines of cross-section data
    bool ReadingCrossSections = true;
    while (ReadingCrossSections) {
        std::string CSline;
        getline(Analysis, CSline);
        if (CSline == "**********") break;
        
        std::string CS[2];
        // Should be 2 inputs on this line - break it into its components
        CS[0] = parseCoordinatePair(CSline, 0);
        CS[1] = parseCoordinatePair(CSline, 1);
        if (CS[0] == "XZ") {
            CrossSectionPlane.push_back(0);
            if (CS[1] == "END") CrossSectionLocation.push_back(ny-2);
            else if (CS[1] == "D") CrossSectionLocation.push_back(FullDomainCenterY);
            else {
                CrossSectionLocation.push_back(stoi(CS[1], nullptr, 10));
                if (CrossSectionLocation[NumberOfCrossSections] >= ny) {
                    CrossSectionLocation[NumberOfCrossSections] = ny-1;
                    std::cout << "Adjusting cross-section number " << NumberOfCrossSections+1 << " y value, which was larger than the domain size in y" << std::endl;
                }
            }
        }
        else if (CS[0] == "YZ") {
            CrossSectionPlane.push_back(1);
            if (CS[1] == "END") CrossSectionLocation.push_back(nx-2);
            else if (CS[1] == "D") CrossSectionLocation.push_back(FullDomainCenterX);
            else {
                CrossSectionLocation.push_back(stoi(CS[1], nullptr, 10));
                if (CrossSectionLocation[NumberOfCrossSections] >= nx) {
                    CrossSectionLocation[NumberOfCrossSections] = nx-1;
                    std::cout << "Adjusting .ang (pyEBSD) cross-section number " << NumberOfCrossSections+1 << " x value, which was larger than the domain size in x" << std::endl;
                }
            }
        }
        else {
            CrossSectionPlane.push_back(2);
            if (CS[1] == "END") CrossSectionLocation.push_back(nz-2);
            else if (CS[1] == "D") CrossSectionLocation.push_back(FullDomainCenterZ);
            else {
                CrossSectionLocation.push_back(stoi(CS[1], nullptr, 10));
                if (CrossSectionLocation[NumberOfCrossSections] >= nz) {
                    CrossSectionLocation[NumberOfCrossSections] = nz-1;
                    std::cout << "Adjusting .ang (pyEBSD) cross-section number " << NumberOfCrossSections+1 << " z value, which was larger than the domain size in z" << std::endl;
                }
            }
        }
        std::cout << ".ang (pyEBSD) Cross-section number " << NumberOfCrossSections << " has parsed string data " << CrossSectionPlane[NumberOfCrossSections] << " " << CrossSectionLocation[NumberOfCrossSections] << std::endl;
        NumberOfCrossSections++;
    }
    skipLines(Analysis);

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

    if (XMinString == "D") XMin = 50;
    else XMin = stoi(XMinString, nullptr, 10);
    
    if (XMaxString == "D") XMax = nx-50;
    else XMax = stoi(XMaxString, nullptr, 10);
    
    if (YMinString == "D") YMin = 50;
    else YMin = stoi(YMinString, nullptr, 10);
    
    if (YMaxString == "D") YMax = ny-50;
    else YMax = stoi(YMaxString, nullptr, 10);
    if (ZMinString == "D") {
        ZMin = FindTopOrBottom(LayerID, XMin, XMax, YMin, YMax, nz, 0, "Low");
    }
    else ZMin = stoi(ZMinString, nullptr, 10);
    if (ZMaxString == "D") {
        ZMax = FindTopOrBottom(LayerID, XMin, XMax, YMin, YMax, nz, NumberOfLayers-1, "High");
    }
    else ZMax = stoi(ZMaxString, nullptr, 10);

    std::cout << "XYZ region for grain statistics is X = " << XMin << " through " << XMax << "; Y = " << YMin << " through " << YMax << "; Z = " << ZMin << " through " << ZMax << std::endl;
    // Read Y/N options for analyzing microstructure data
    AnalysisTypes[0] = parseInputBool(Analysis, "misorientation");
    AnalysisTypes[1] = parseInputBool(Analysis, "grain volume");
    AnalysisTypes[2] = parseInputBool(Analysis, "mean aspect ratio");
    AnalysisTypes[3] = parseInputBool(Analysis, "mean grain area");
    AnalysisTypes[4] = parseInputBool(Analysis, "mean weighted grain area");
    AnalysisTypes[5] = parseInputBool(Analysis, "mean grain height");
    AnalysisTypes[6] = parseInputBool(Analysis, "grain width distribution");
    AnalysisTypes[7] = parseInputBool(Analysis, "volume for pole figure analysis");
    Analysis.close();
}

void ParseGrainOrientationFiles(std::string RotationFilename, std::string EulerFilename, int NumberOfOrientations, float*** GrainUnitVector, float**) {
    
    // Read file of rotation matrix form grain orientations
    std::ifstream O;
    O.open(RotationFilename);

    // Line 1 is the number of orientation values to read (already known)
    std::string ValueRead;
    getline(O, ValueRead);

    // Populate data structure for grain unit vectors
    for (int i = 0; i < NumberOfOrientations; i++) {
        std::string s;
        if (!getline(O, s))
            break;
        std::istringstream ss(s);
        int Comp = 0;
        int UVNumber = 0;
        while (ss) { // This is the 3 grain orientation angles
            std::string s;
            if (!getline(ss, s, ','))
                break;
            float ReadGO = atof(s.c_str());
            // X,Y,Z of a single unit vector
            GrainUnitVector[i][UVNumber][Comp] = ReadGO;
            Comp++;
            if (Comp > 2) {
                Comp = 0;
                UVNumber++;
            }
        }
    }
    O.close();
    
//    Read file of euler form grain orientations - to be used in the future with pyEBSD data analysis
//    std::ifstream OE;
//    OE.open(EulerFilename);
//    for (int i = 0; i < NumberOfOrientations; i++) {
//        std::string s;
//        if (!getline(O, s))
//            break;
//        std::istringstream ss(s);
//        int Comp = 0;
//        while (ss) { // This is the 3 grain orientation angles
//            std::string s;
//            if (!getline(ss, s, ','))
//                break;
//            float ReadGO = atof(s.c_str());
//            // X,Y,Z of a single unit vector
//            GrainEulerAngles[i][Comp] = ReadGO;
//            Comp++;
//            if (Comp > 2) Comp = 0;
//        }
//    }
//    OE.close();
}