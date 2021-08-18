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

// Verify the required input was included with the correct format, and parse the input
// 2 possible keywords to check against - WhichKey denotes which of the two was found
std::string parseInputMultiple(std::ifstream &stream, std::string key1, std::string key2, int &WhichKey) {

    std::size_t colon;
    std::string line, val;
    std::string actual_key = getKey(stream, line, colon);

    // Check for keywords
    if (actual_key.find(key1) == std::string::npos) {
        // Keyword 1 not found - check for Keyword 2
        if (actual_key.find(key2) == std::string::npos) {
            // Neither possible input found
            std::string error =
                "Required input not present: Neither " + key1 + " nor " + key2 + " was found in the input file";
            throw std::runtime_error(error);
        }
        else {
            // The second keyword was found
            WhichKey = 2;
        }
    }
    else {
        // First keyword was found
        WhichKey = 1;
    }
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
        for (int DummyLines=0; DummyLines<3; DummyLines++) {
            getline(InputDataStream, line);
        }
        getline(InputDataStream, line);
        std::size_t EndString = line.find("layers");
        std::string NumLayersString = line.substr(0, EndString);
        NumberOfLayers = stoi(NumLayersString, nullptr, 10);
    }
    InputDataStream.close();
    
}

void InitializeData(std::string InputFile, int &nx, int &ny, int &nz, int*** GrainID, int*** LayerID, int*** Melted) {
    
    std::ifstream InputDataStream;
    InputDataStream.open(InputFile);
    if (!(InputDataStream)) throw std::runtime_error("Error: Cannot find ExaCA microstructure file");
    std::string line;
    // 10 unused lines at beginning of Paraview file from header
    for (int DummyLines=0; DummyLines<10; DummyLines++) {
        getline(InputDataStream, line);
    }
    // Read Grain ID data (each line is at constant k, varied i/j)
    for (int k=0; k<nz; k++) {
        int i = 0;
        int j = 0;
        getline(InputDataStream, line);
        // Used to split string around spaces
        std::istringstream ss(line);
        std::string GIDVal; // for storing each separated string
        for (int j=0; j<ny; j++) {
            for (int i=0; i<nx; i++) {
                int val;
                ss >> val;
                GrainID[k][i][j] = val;
            }
        }
    }
    std::cout << "Grain ID data read" << std::endl;
    // 2 more unused lines separating Grain ID and Layer ID data
    for (int DummyLines=0; DummyLines<2; DummyLines++) {
        getline(InputDataStream, line);
    }
    // Read Layer ID data (each line is at constant k, varied i/j)
    for (int k=0; k<nz; k++) {
        int i = 0;
        int j = 0;
        getline(InputDataStream, line);
        // Used to split string around spaces
        std::istringstream ss(line);
        std::string GIDVal; // for storing each separated string
        for (int j=0; j<ny; j++) {
            for (int i=0; i<nx; i++) {
                int val;
                ss >> val;
                LayerID[k][i][j] = val;
            }
        }
    }
    std::cout << "Layer ID data read" << std::endl;
    // 2 more unused lines separating Layer ID and Melted data
    for (int DummyLines=0; DummyLines<2; DummyLines++) {
        getline(InputDataStream, line);
    }
    // Read Melted data (each line is at constant k, varied i/j)
    for (int k=0; k<nz; k++) {
        int i = 0;
        int j = 0;
        getline(InputDataStream, line);
        // Used to split string around spaces
        std::istringstream ss(line);
        std::string GIDVal; // for storing each separated string
        for (int j=0; j<ny; j++) {
            for (int i=0; i<nx; i++) {
                int val;
                ss >> val;
                Melted[k][i][j] = val;
            }
        }
    }
    std::cout << "Melting data read" << std::endl;
    InputDataStream.close();
    
}

// Read the analysis file to determine which analysis operations will be performed on the given ExaCA data
void ParseAnalysisFile(std::string OutputAnalysisFile, std::string &RotationFilename, std::string &EulerFilename, int &NumberOfOrientations, bool* AnalysisTypes, std::vector <int> &CenterX_RVE, std::vector <int> &CenterY_RVE, std::vector <int> &CenterZ_RVE, std::vector <int> &Size_RVE, std::vector <bool> &ExaConstitPrint, std::vector <bool> &pyEBSDPrint, int &NumberOfRVEs, std::vector <int> &ANGCrossSectionPlane, std::vector <int> &ANGCrossSectionLocation, int &NumberOfANGCrossSections, int &XMin, int &XMax, int &YMin, int &YMax, int nx, int ny, int nz, int*** LayerID, int*** Melted, int NumberOfLayers) {
    
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
        std::string RVE[6];
        // Should be 6 inputs on this line - break it into its components
        size_t compstart, compend;
        for (int comp=0; comp<6; comp++) {
            if (comp == 0) compstart = RVEline.find("[");
            else compstart = compend; // from last parsed string
            if (comp == 5) compend = RVEline.find("]");
            else compend = RVEline.find(",", compstart + 1);
            RVE[comp] = RVEline.substr(compstart + 1, compend - compstart - 1);
            // remove any whitespace from parsed string
            std::regex r("\\s+");
            RVE[comp] = std::regex_replace(RVE[comp], r, "");
        }
        // RVE[0] = Size, RVE[1] = Center X, RVE[2] = CenterY, RVE[3] = CenterZ
        // RVE[4] = print .ang file Y/N, RVE[5] = print ExaConstit file Y/N
        // Default option for Center X or Center Y takes RVE X/Y center to be at the center of the data
        // Default option for Center Z (if single layer dataset) takes the RVE Z center to be at the center
        // Default option for Center Z (if multilayer dataset) takes Z center such that we avoid the last layer if possible, but are as close to the top as possible
        // RVE Size
        Size_RVE[NumberOfRVEs] = stoi(RVE[0], nullptr, 10);
        // Center X and Center Y
        if (RVE[1] == "D") CenterX_RVE[NumberOfRVEs] = nx/2;
        else CenterX_RVE[NumberOfRVEs] = stoi(RVE[1], nullptr, 10);
        if (RVE[2] == "D") CenterY_RVE[NumberOfRVEs] = ny/2;
        else CenterY_RVE[NumberOfRVEs] = stoi(RVE[2], nullptr, 10);
        // Center Z
        if (RVE[3] != "D") CenterZ_RVE[NumberOfRVEs] = stoi(RVE[2], nullptr, 10);
        else if ((RVE[3] != "D")&&(NumberOfLayers == 1)) CenterZ_RVE[NumberOfRVEs] = nz/2;
        else {
            int nztop;
            bool SearchingForTop = true;
            int i = 0;
            int j = 0;
            int k = nz-1;
            while (SearchingForTop) {
                if (LayerID[k][i][j] != NumberOfLayers-1) {
                    j++;
                    if (j == ny) {
                        if (i == nx-1) {
                            // All X and Y coordinates at this Z have been checked, and this Z coordinate is not part of the last depositied layer
                            nztop = k;
                            SearchingForTop = false;
                        }
                        else {
                            i++;
                            j = 0;
                        }
                    }
                }
                else {
                    // Search next Z coordinate down - this one was part of the last deposited layer
                    k--;
                    i = 0;
                    j = 0;
                    if (k < 0) SearchingForTop = false;
                }
            }
            int nzbottom = nztop - Size_RVE[NumberOfRVEs];
            if (nzbottom < 0) {
                // We can't avoid getting at least part of the top layer in the RVE, so default to taking the RVE from the center of the domain in z
                CenterZ_RVE[NumberOfRVEs] = nz/2;
            }
            else {
                CenterZ_RVE[NumberOfRVEs] = (nztop + nzbottom)/2;
            }
        }
        if (RVE[4] == "Y") ExaConstitPrint[NumberOfRVEs] = true;
        else ExaConstitPrint[NumberOfRVEs] = 0;
        if (RVE[5] == "Y") pyEBSDPrint[NumberOfRVEs] = true;
        else pyEBSDPrint[NumberOfRVEs] = false;
        std::cout << "RVE number " << NumberOfRVEs << " has size " << Size_RVE[NumberOfRVEs] << " and center location " << CenterX_RVE[NumberOfRVEs] << "," << CenterY_RVE[NumberOfRVEs] << " " << CenterZ_RVE[NumberOfRVEs] << std::endl;
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
        size_t linestart = CSline.find("[");
        size_t linebreak = CSline.find(",");
        size_t lineend = CSline.find("]");
        CS[0] = CSline.substr(linestart + 1, linebreak - linestart - 1);
        CS[1] = CSline.substr(linebreak + 1, lineend - linebreak - 1);
        // remove any whitespace from parsed strings
        std::regex r("\\s+");
        CS[0] = std::regex_replace(CS[0], r, "");
        CS[1] = std::regex_replace(CS[1], r, "");
        
        if (CS[0] == "XZ") {
            ANGCrossSectionPlane[NumberOfANGCrossSections] = 0;
            if (CS[1] == "END") ANGCrossSectionLocation[NumberOfANGCrossSections] = ny-1;
            else if (CS[1] == "D") ANGCrossSectionLocation[NumberOfANGCrossSections] = ny/2;
            else {
                ANGCrossSectionLocation[NumberOfANGCrossSections] = stoi(CS[1], nullptr, 10);
                if (ANGCrossSectionLocation[NumberOfANGCrossSections] >= ny) {
                    ANGCrossSectionLocation[NumberOfANGCrossSections] = ny-1;
                    std::cout << "Adjusting .ang (pyEBSD) cross-section number " << NumberOfANGCrossSections+1 << " y value, which was larger than the domain size in y" << std::endl;
                }
            }
        }
        else if (CS[0] == "YZ") {
            ANGCrossSectionPlane[NumberOfANGCrossSections] = 1;
            if (CS[1] == "END") ANGCrossSectionLocation[NumberOfANGCrossSections] = nx-1;
            else if (CS[1] == "D") ANGCrossSectionLocation[NumberOfANGCrossSections] = nx/2;
            else {
                ANGCrossSectionLocation[NumberOfANGCrossSections] = stoi(CS[1], nullptr, 10);
                if (ANGCrossSectionLocation[NumberOfANGCrossSections] >= nx) {
                    ANGCrossSectionLocation[NumberOfANGCrossSections] = nx-1;
                    std::cout << "Adjusting .ang (pyEBSD) cross-section number " << NumberOfANGCrossSections+1 << " x value, which was larger than the domain size in x" << std::endl;
                }
            }
        }
        else {
            ANGCrossSectionPlane[NumberOfANGCrossSections] = 2;
            if (CS[1] == "END") ANGCrossSectionLocation[NumberOfANGCrossSections] = nz-1;
            else if (CS[1] == "D") ANGCrossSectionLocation[NumberOfANGCrossSections] = nz/2;
            else {
                ANGCrossSectionLocation[NumberOfANGCrossSections] = stoi(CS[1], nullptr, 10);
                if (ANGCrossSectionLocation[NumberOfANGCrossSections] >= nz) {
                    ANGCrossSectionLocation[NumberOfANGCrossSections] = nz-1;
                    std::cout << "Adjusting .ang (pyEBSD) cross-section number " << NumberOfANGCrossSections+1 << " z value, which was larger than the domain size in z" << std::endl;
                }
            }
        }
        std::cout << ".ang (pyEBSD) Cross-section number " << NumberOfANGCrossSections+1 << " has parsed string data " << ANGCrossSectionPlane[NumberOfANGCrossSections] << " " << ANGCrossSectionLocation[NumberOfANGCrossSections] << std::endl;
        NumberOfANGCrossSections++;
    }
    skipLines(Analysis);
    std::string XYSizeString = parseInput(Analysis, "XY Size");
    int XYSize;
    if (XYSizeString != "D") XYSize = stoi(XYSizeString, nullptr, 10);
    else {
        if (NumberOfLayers == 1) {
            // XY size just involves cropping out 20 cells from each of the +/- X and Y boundaries
            if ((nx < 60)||(ny < 60)) {
                // This simulation volume is too small for any meaningful grain statistics
                // Will not calculate any of the requested Y/N options
                std::cout << "WARNING: This microstructure is too small to obtain any meaningful grain statistics - ignoring any toggled print options" << std::endl;
                for (int i=0; i<8; i++) {
                    AnalysisTypes[i] = false;
                }
                return;
            }
            else {
                // Selection is as close to the center of the XY cross-section as possible
                if (nx < ny) {
                    XMin = 20;
                    XMax = nx - 20;
                    YMin = ny/2 - ((nx-40)/2);
                    YMax = ny/2 + ((nx-40)/2);
                }
                else {
                    YMin = 20;
                    YMax = ny - 20;
                    XMin = nx/2 - ((ny-40)/2);
                    XMax = nx/2 + ((ny-40)/2);
                }
            }
        }
        // Determine the XY size based on the boundary between substrate/simulated microstructure at the top surface
        // Rough estimation of domain "bottom" - all solidification events at a Z height that were associated
        // with a layer other than the first one
        XMin = 0;
        XMax = nx-1;
        YMin = 0;
        YMax = ny-1;
        for (int i=0; i<nx/2; i++) {
            for (int j=0; j<ny; j++) {
                if (!(Melted[nz-1][i][j])) {
                    if (i > XMin) XMin = i;
                }
            }
        }
        for (int i=nx/2; i<nx; i++) {
            for (int j=0; j<ny; j++) {
                if (!(Melted[nz-1][i][j])) {
                    if (i < XMax) XMax = i;
                }
            }
        }
        for (int i=0; i<nx; i++) {
            for (int j=0; j<ny/2; j++) {
                if (!(Melted[nz-1][i][j])) {
                    if (j > YMin) YMin = j;
                }
            }
        }
        for (int i=0; i<nx; i++) {
            for (int j=ny/2; j<ny; j++) {
                if (!(Melted[nz-1][i][j])) {
                    if (j < YMax) YMax = j;
                }
            }
        }
        std::cout << "XY cross-section for grain statistics is X = " << XMin << " through " << XMax << "; Y = " << YMin << " through " << YMax << std::endl;
    }
    // Read Y/N options for analyzing microstructure data
    AnalysisTypes[0] = parseInputBool(Analysis, "misorientation");
    AnalysisTypes[1] = parseInputBool(Analysis, "grain volume");
    AnalysisTypes[2] = parseInputBool(Analysis, "mean grain area");
    AnalysisTypes[3] = parseInputBool(Analysis, "mean weighted grain area");
    AnalysisTypes[4] = parseInputBool(Analysis, "mean grain height");
    AnalysisTypes[5] = parseInputBool(Analysis, "grain width distribution");
    AnalysisTypes[6] = parseInputBool(Analysis, "grain width in X and Y");
    AnalysisTypes[7] = parseInputBool(Analysis, "mean volume-weighted aspect ratio");

}

void ParseGrainOrientationFiles(std::string RotationFilename, std::string EulerFilename, int NumberOfOrientations, float*** GrainUnitVector, float** GrainEulerAngles) {
    
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
    
    // Read file of euler form grain orientations
    // std::ifstream OE;
    // OE.open(EulerFilename);

    // Populate data structure for euler angles
    // Work in progress - not currently used
    // OE.close();
}
