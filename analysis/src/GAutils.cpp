// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "GAutils.hpp"
#include "CAfunctions.hpp"
#include "CAparsefiles.hpp"

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

// Check the format of the log file:
// 0: text file without colon-separated format
// 1: text file with colon separator on all lines
// 2: json format
int checkLogFormat(std::string LogFile) {
    int LogFormat;
    std::ifstream InputDataStream;
    InputDataStream.open(LogFile);
    if (!(InputDataStream))
        throw std::runtime_error("Error: Cannot find ExaCA log file");
    bool JsonInputFormat = checkInputFileFormat(LogFile, 0);
    if (JsonInputFormat)
        LogFormat = 2;
    else {
        std::string testline;
        // If the version is listed on the first line, this is the newer of the non-json log file formats
        getline(InputDataStream, testline);
        std::size_t found_version = testline.find("ExaCA version");
        if (found_version == std::string::npos)
            LogFormat = 0;
        else
            LogFormat = 1;
    }
    return LogFormat;
}

// Parse log file using json format
#ifdef ExaCA_ENABLE_JSON
void ParseLogFile(std::string LogFile, int &nx, int &ny, int &nz, double &deltax, int &NumberOfLayers,
                  std::vector<double> &XYZBounds, std::string &RotationFilename, std::string &EulerAnglesFilename,
                  std::string &RGBFilename, bool OrientationFilesInInput) {

    std::ifstream InputDataStream;
    InputDataStream.open(LogFile);
    nlohmann::json logdata = nlohmann::json::parse(InputDataStream);

    // X, Y, Z bounds of domain (in microns)
    XYZBounds[0] = logdata["Domain"]["XBounds"][0];
    XYZBounds[1] = logdata["Domain"]["YBounds"][0];
    XYZBounds[2] = logdata["Domain"]["ZBounds"][0];
    XYZBounds[3] = logdata["Domain"]["XBounds"][1];
    XYZBounds[4] = logdata["Domain"]["YBounds"][1];
    XYZBounds[5] = logdata["Domain"]["ZBounds"][1];
    // Cell size (in microns)
    deltax = logdata["Domain"]["CellSize"];
    // Number of cells per domain direction
    nx = logdata["Domain"]["Nx"];
    ny = logdata["Domain"]["Ny"];
    nz = logdata["Domain"]["Nz"];
    std::string SimulationType = logdata["SimulationType"];
    if (SimulationType == "C")
        NumberOfLayers = 1;
    else
        NumberOfLayers = logdata["Domain"]["NumberOfLayers"];
    if (OrientationFilesInInput)
        std::cout << "Note: orientation filename specified in log file will be used, overriding value from analysis "
                     "input file"
                  << std::endl;
    RotationFilename = logdata["GrainOrientationFile"];
    if (RotationFilename.find("GrainOrientationVectors.csv") != std::string::npos) {
        // Default files for euler angles and RGB mapping
        EulerAnglesFilename = "GrainOrientationEulerAnglesBungeZXZ.csv";
        RGBFilename = "GrainOrientationRGB_IPF-Z.csv";
    }
    else {
        // Custom files for euler angles and RGB mapping based on rotation filename
        std::size_t baseorientation_startpos = RotationFilename.find_last_of("/");
        std::size_t endpos = RotationFilename.find_last_of(".");
        std::string baseorientationname =
            RotationFilename.substr(baseorientation_startpos + 1, endpos - baseorientation_startpos - 1);
        std::size_t customname_startpos = baseorientationname.find_last_of("_");
        std::string customname = baseorientationname.substr(customname_startpos + 1, endpos - customname_startpos - 1);
        EulerAnglesFilename = "GrainOrientationEulerAnglesBungeZXZ_" + customname + ".csv";
        RGBFilename = "GrainOrientationRGB_IPF-Z_" + customname + ".csv";
    }
    // Full path to install location was already given for the rotations file, need to check install location for
    // other files
    EulerAnglesFilename = checkFileInstalled(EulerAnglesFilename, 0);
    RGBFilename = checkFileInstalled(RGBFilename, 0);
    checkFileNotEmpty(RotationFilename);
    checkFileNotEmpty(EulerAnglesFilename);
    checkFileNotEmpty(RGBFilename);
    InputDataStream.close();
}
#endif

void ParseLogFile_Old(std::string LogFile, int &nx, int &ny, int &nz, double &deltax, int &NumberOfLayers,
                      std::vector<double> &XYZBounds, std::string &RotationFilename, std::string &EulerAnglesFilename,
                      std::string &RGBFilename, bool OrientationFilesInInput) {

    std::vector<std::string> LogInputs = {
        "Lower bound of domain in x", // Input 0
        "Lower bound of domain in y", // Input 1
        "Lower bound of domain in z", // Input 2
        "Upper bound of domain in x", // Input 3
        "Upper bound of domain in y", // Input 4
        "Upper bound of domain in z", // Input 5
        "Cell size",                  // Input 6
        "Domain size in x",           // Input 7
        "Domain size in y",           // Input 8
        "Domain size in z",           // Input 9
        "Simulation was type",        // Input 10
        "Number of layers simulated", // Input 11
        "Grain orientation file"      // Input 12
    };
    int NumLogInputs = LogInputs.size();
    std::vector<std::string> LogInputsRead(NumLogInputs);

    std::ifstream InputDataStream;
    InputDataStream.open(LogFile);
    std::string line;
    while (std::getline(InputDataStream, line)) {
        // Line with *** denotes end of relevant portion of log file
        if (line.find("***") == std::string::npos) {
            // Check against inputs of interest
            bool FoundYN = parseInputFromList(line, LogInputs, LogInputsRead, NumLogInputs);
            if (FoundYN)
                std::cout << line << std::endl;
        }
        else
            break;
    }
    for (int n = 0; n < 6; n++) {
        XYZBounds[n] = getInputDouble(LogInputsRead[n]);
    }
    deltax = getInputDouble(LogInputsRead[6]);
    nx = getInputInt(LogInputsRead[7]);
    ny = getInputInt(LogInputsRead[8]);
    nz = getInputInt(LogInputsRead[9]);
    std::string SimulationType = LogInputsRead[10];
    if (SimulationType == "C")
        NumberOfLayers = 1;
    else
        NumberOfLayers = getInputInt(LogInputsRead[11]);
    // Was orientation file name specified in the log file? Otherwise get from the analysis file
    if (LogInputsRead[12].empty()) {
        if (OrientationFilesInInput)
            std::cout << "Warning: specification of orientation file names will be required in the log file in a "
                         "future release, will no longer be accepted in the analysis input file"
                      << std::endl;
        else
            throw std::runtime_error(
                "Error: Orientation file names were not given in the log file nor the analysis input file");
    }
    else {
        OrientationFilesInInput = true;
        RotationFilename = LogInputsRead[12];
        if (OrientationFilesInInput)
            std::cout << "Note: using orientation file names from log file, not those specified in analysis input file"
                      << std::endl;
        if (RotationFilename.find("GrainOrientationVectors.csv") != std::string::npos) {
            // Default files for euler angles and RGB mapping
            EulerAnglesFilename = "GrainOrientationEulerAnglesBungeZXZ.csv";
            RGBFilename = "GrainOrientationRGB_IPF-Z.csv";
        }
        else {
            // Custom files for euler angles and RGB mapping based on rotation filename
            std::size_t baseorientation_startpos = RotationFilename.find_last_of("/");
            std::size_t endpos = RotationFilename.find_last_of(".");
            std::string baseorientationname =
                RotationFilename.substr(baseorientation_startpos + 1, endpos - baseorientation_startpos - 1);
            std::size_t customname_startpos = baseorientationname.find_last_of("_");
            std::string customname =
                baseorientationname.substr(customname_startpos + 1, endpos - customname_startpos - 1);
            EulerAnglesFilename = "GrainOrientationEulerAnglesBungeZXZ_" + customname + ".csv";
            RGBFilename = "GrainOrientationRGB_IPF-Z_" + customname + ".csv";
        }
        // Full path to install location was already given for the rotations file, need to check install location for
        // other files
        EulerAnglesFilename = checkFileInstalled(EulerAnglesFilename, 0);
        RGBFilename = checkFileInstalled(RGBFilename, 0);
        checkFileNotEmpty(RotationFilename);
        checkFileNotEmpty(EulerAnglesFilename);
        checkFileNotEmpty(RGBFilename);
    }
}

void ParseLogFile_OldNoColon(std::string LogFile, int &nx, int &ny, int &nz, double &deltax, int &NumberOfLayers,
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
// ASCII data at each Z value is separated by a newline
void ReadASCIIField(std::ifstream &InputDataStream, int nx, int ny, int nz, ViewI3D_H FieldOfInterest) {
    for (int k = 0; k < nz; k++) {
        // Get line from file
        std::string line;
        getline(InputDataStream, line);
        // Parse string at spaces
        std::istringstream ss(line);
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                FieldOfInterest(k, i, j) = ParseASCIIData<int>(ss);
            }
        }
    }
}

// Reads binary string from a paraview file, converts to int field, and places data in the appropriate data structure
// Each field consists of a single binary string (no newlines)
void ReadBinaryField(std::ifstream &InputDataStream, int nx, int ny, int nz, ViewI3D_H FieldOfInterest,
                     std::string FieldName) {
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                // Store converted values in view - LayerID data is a short int, GrainID data is an int
                if (FieldName == "GrainID")
                    FieldOfInterest(k, i, j) = ReadBinaryData<int>(InputDataStream, true);
                else if (FieldName == "LayerID")
                    FieldOfInterest(k, i, j) = ReadBinaryData<short>(InputDataStream, true);
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
                    std::string &RGBFilename, bool &OrientationFilesInInput) {

    // The analysis file should be in the analysis subdirectory of ExaCA
    std::cout << "Looking for analysis file: " << AnalysisFile << std::endl;
    std::ifstream Analysis;
    Analysis.open(AnalysisFile);
    if (!(Analysis))
        throw std::runtime_error("Error: Cannot find ExaCA analysis file");
    skipLines(Analysis, "*****");

    std::vector<std::string> NamesOfFilesRequired = {"name of microstructure (.vtk) file",
                                                     "name of data files of output resulting from this analysis"};
    std::vector<std::string> NamesOfFilesOptional = {
        "Path to/name of log file", "file of grain orientations (rotation matrix form)",
        "file of grain orientations (Bunge Euler angle form ZXZ)",
        "file corresponding to grain orientation IPF-Z colors as fractional RGB values"};
    std::vector<std::string> NamesOfFilesDeprecated = {"Path to/name of log (.log) file"};
    // If orientation file names are given here, store these values but override if given in the log file with a warning
    int NumRequiredInputs = NamesOfFilesRequired.size();
    int NumOptionalInputs = NamesOfFilesOptional.size();
    int NumDeprecatedInputs = NamesOfFilesDeprecated.size();
    std::vector<std::string> RequiredFilesRead(NumRequiredInputs);
    std::vector<std::string> OptionalFilesRead(NumOptionalInputs);
    std::vector<std::string> DeprecatedFilesRead(NumDeprecatedInputs);

    for (int i = 0; i < 6; i++) {
        std::string line;
        std::getline(Analysis, line);
        if (line == "**********")
            break;
        bool LineFound = parseInputFromList(line, NamesOfFilesRequired, RequiredFilesRead, NumRequiredInputs);
        if (!(LineFound)) {
            LineFound = parseInputFromList(line, NamesOfFilesOptional, OptionalFilesRead, NumOptionalInputs);
            if (!(LineFound)) {
                LineFound = parseInputFromList(line, NamesOfFilesDeprecated, DeprecatedFilesRead, NumDeprecatedInputs);
                if (!(LineFound))
                    std::cout << "Warning: the line " << line
                              << " did not match any analysis file options expected by ExaCA and will be ignored"
                              << std::endl;
            }
        }
    }
    Analysis.close();

    // Check that all required arguments were present
    for (int i = 0; i < NumRequiredInputs; i++) {
        if (RequiredFilesRead[i].empty()) {
            std::string error = "Error: Required input " + NamesOfFilesRequired[i] + " not found in analysis file";
            throw std::runtime_error(error);
        }
    }
    // Required input files
    MicrostructureFile = RequiredFilesRead[0];
    OutputFileName = RequiredFilesRead[1];
    // Log file is required - either with the .log extension (deprecated) or with no extension
    if (!(DeprecatedFilesRead[0].empty())) {
        LogFile = DeprecatedFilesRead[0];
        std::cout << "Note: Compatibility with .log files will be removed in a future release, beyond which .json "
                     "format for the log data will be required";
    }
    else if (!(OptionalFilesRead[0].empty())) {
#ifdef ExaCA_ENABLE_JSON
        // If json is enabled, first check if log is in json format
        LogFile = OptionalFilesRead[0] + ".json";
        bool Json_found = checkFileExists(LogFile, 0, false);
        // If json is enabled and the .json file does not exist, try the .log file
        if (!(Json_found))
            LogFile = OptionalFilesRead[0] + ".log";
#else
        LogFile = OptionalFilesRead[0] + ".log";
#endif
    }
    else
        throw std::runtime_error("Error: Required input Path to/name of log file not found in analysis file");

    // Check that required input files are not empty
    checkFileNotEmpty(LogFile);
    checkFileNotEmpty(MicrostructureFile);

    // Check for optional inputs - if not included, get these from the log file
    if ((OptionalFilesRead[1].empty()) && (OptionalFilesRead[2].empty()) && (OptionalFilesRead[3].empty()))
        OrientationFilesInInput = false;
    else {
        OrientationFilesInInput = true;
        RotationFilename = OptionalFilesRead[1];
        if (OptionalFilesRead[2].empty()) {
            // File not given, use default file
            EulerAnglesFilename = "GrainOrientationEulerAnglesBungeZXZ.csv";
            std::cout << "Defaulting to use of GrainOrientationEulerAnglesBungeZXZ.csv for euler angles" << std::endl;
        }
        else {
            // File is given, allow printing of pole figure data using updated format
            EulerAnglesFilename = OptionalFilesRead[2];
        }
        if (OptionalFilesRead[3].empty()) {
            // File not given, use default file
            RGBFilename = "GrainOrientationRGB_IPF-Z.csv";
            std::cout << "The default file GrainOrientationRGB_IPF-Z.csv will be used to map colors to orientations"
                      << std::endl;
        }
        else {
            // File is given, use this value
            RGBFilename = OptionalFilesRead[3];
        }
        RotationFilename = checkFileInstalled(RotationFilename, 0);
        checkFileNotEmpty(RotationFilename);
        EulerAnglesFilename = checkFileInstalled(EulerAnglesFilename, 0);
        checkFileNotEmpty(EulerAnglesFilename);
        RGBFilename = checkFileInstalled(RGBFilename, 0);
        checkFileNotEmpty(RGBFilename);
    }
}

// Ensure that the appropriate files exist - orientation files should be installed, other files should start with the
// BaseFileName value given from the command line
void CheckInputFiles(std::string &LogFile, std::string MicrostructureFile, std::string &RotationFilename,
                     std::string &EulerAnglesFilename, std::string &RGBFilename) {

    // Names of files
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
    bool Binaryvtk = false;
    // 8 header lines at beginning of Paraview file from header
    // 3rd line tells whether data is in ASCII or binary format
    for (int HeaderLines = 0; HeaderLines < 8; HeaderLines++) {
        getline(InputDataStream, line);
        if (HeaderLines == 2)
            if (line.find("BINARY") != std::string::npos)
                Binaryvtk = true;
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
        if (FoundGrainID != std::string::npos) {
            if (Binaryvtk)
                ReadBinaryField(InputDataStream, nx, ny, nz, GrainID, "GrainID");
            else
                ReadASCIIField(InputDataStream, nx, ny, nz, GrainID);
        }
        else if (FoundLayerID != std::string::npos) {
            if (Binaryvtk)
                ReadBinaryField(InputDataStream, nx, ny, nz, LayerID, "LayerID");
            else
                ReadASCIIField(InputDataStream, nx, ny, nz, LayerID);
        }
        else if (FoundMelted != std::string::npos) {
            std::cout << "Note: Melted data is no longer used and will be ignored" << std::endl;
            // Note: version of ExaCA that printed output with melted data is older than the version that prints binary
            // vtk output, melted field will always consist of ASCII characters if present in the output
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
                       std::vector<bool> &BimodalAnalysis, std::vector<std::string> &CSLabels) {

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
        // cross-section data Y/N). This is no longer linked to whether the euler angles file was provided in the
        // analysis file, since the default euler angles file is now used if not otherwise given
        bool NewAnalysisInputFormatYN;
        std::size_t found1 = CSline.find("]");
        std::size_t found2 = CSline.find("[");
        // New format has no brackets, old format has 2 comma separated values with brackets
        if ((found1 == std::string::npos) && (found2 == std::string::npos))
            NewAnalysisInputFormatYN = true;
        else
            NewAnalysisInputFormatYN = false;
        // For the "new" format, there are either 4 or 5 comma-separated components on this line, depending on whether
        // the option of analyzing small and large area grains separately is given. If it is not given, a single
        // distribution of grain areas in the cross-section will be analyzed
        int NumComponents;
        if (NewAnalysisInputFormatYN)
            NumComponents = std::count(CSline.begin(), CSline.end(), ',') + 1;
        else
            NumComponents = 2;
        std::vector<std::string> CS(NumComponents);
        if (NewAnalysisInputFormatYN) {
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

// Get the scaling factor to convert a number of cells into a length, area, or volume
double convertToMicrons(double deltax, std::string RegionType) {
    double CellToSizeScaling;
    if (RegionType == "volume")
        CellToSizeScaling = pow(deltax, 3) * pow(10, 18);
    else if (RegionType == "area")
        CellToSizeScaling = pow(deltax, 2) * pow(10, 12);
    else if (RegionType == "length")
        CellToSizeScaling = deltax * pow(10, 6);
    else
        throw std::runtime_error("Error: unknown region type");
    return CellToSizeScaling;
}

// Get the scaling factor to convert a length, area, or volume into a number of cells
double convertToCells(double deltax, std::string RegionType) {
    double SizeToCellScaling;
    if (RegionType == "volume")
        SizeToCellScaling = 1.0 / (pow(deltax, 3) * pow(10, 18));
    else if (RegionType == "area")
        SizeToCellScaling = 1.0 / (pow(deltax, 2) * pow(10, 12));
    else if (RegionType == "length")
        SizeToCellScaling = 1.0 / (deltax * pow(10, 6));
    else
        throw std::runtime_error("Error: unknown region type");
    return SizeToCellScaling;
}

// Helper functions used to analyze specific regions of microstructure data
// Subroutines starting with "get" return the data specified
// Subroutines starting with "calc" calculate the quantity specified but do not return it
// Store a vector of integer Grain ID values consisting of the values bounded by [XLow, XHigh], [YLow, YHigh], [ZLow,
// ZHigh] from the 3D view "GrainID"
std::vector<int> getRepresentativeRegionGrainIDs(ViewI3D_H GrainID, const int XLow, const int XHigh, const int YLow,
                                                 const int YHigh, const int ZLow, const int ZHigh,
                                                 const int RepresentativeRegionSize_Cells) {
    std::vector<int> GrainIDVector(RepresentativeRegionSize_Cells);
    int count = 0;
    for (int k = ZLow; k <= ZHigh; k++) {
        for (int i = XLow; i <= XHigh; i++) {
            for (int j = YLow; j <= YHigh; j++) {
                GrainIDVector[count] = GrainID(k, i, j);
                count++;
            }
        }
    }
    return GrainIDVector;
}

// Given an input vector of integer Grain ID values, return an output vector consisting of the unique Grain ID values,
// sorted from lowest to highest. Store the number of grains
std::vector<int> getUniqueGrains(const std::vector<int> GrainIDVector, int &NumberOfGrains) {
    std::vector<int> UniqueGrainIDVector = GrainIDVector;
    std::sort(UniqueGrainIDVector.begin(), UniqueGrainIDVector.end());
    std::vector<int>::iterator it;
    it = std::unique(UniqueGrainIDVector.begin(), UniqueGrainIDVector.end());
    UniqueGrainIDVector.resize(std::distance(UniqueGrainIDVector.begin(), it));
    NumberOfGrains = UniqueGrainIDVector.size();
    return UniqueGrainIDVector;
}

// Create a histogram of orientations for texture determination, using the GrainID values in the volume bounded by
// [XMin,XMax], [YMin,YMax], [ZMin,ZMax] and excluding and cells that did not undergo melting (GrainID = -1)
ViewI_H getOrientationHistogram(int NumberOfOrientations, ViewI3D_H GrainID, ViewI3D_H LayerID, int XMin, int XMax,
                                int YMin, int YMax, int ZMin, int ZMax) {

    // Init histogram values to zero
    ViewI_H GOHistogram("GOHistogram", NumberOfOrientations);
    for (int k = ZMin; k <= ZMax; k++) {
        for (int j = YMin; j <= YMax; j++) {
            for (int i = XMin; i <= XMax; i++) {
                if (LayerID(k, i, j) != -1) {
                    int GOVal = getGrainOrientation(GrainID(k, i, j), NumberOfOrientations);
                    GOHistogram(GOVal)++;
                }
            }
        }
    }
    return GOHistogram;
}

// Create a histogram of orientations for texture determination, using the GrainID values in the vector
ViewI_H getOrientationHistogram(int NumberOfOrientations, std::vector<int> GrainIDVector,
                                int RepresentativeRegionSize_Cells) {

    // Init histogram values to zero
    ViewI_H GOHistogram("GOHistogram", NumberOfOrientations);
    for (int n = 0; n < RepresentativeRegionSize_Cells; n++) {
        int GOVal = getGrainOrientation(GrainIDVector[n], NumberOfOrientations);
        GOHistogram(GOVal)++;
    }
    return GOHistogram;
}

// Given an input vector of integer Grain ID values "GrainIDVector", and an input vector of the
// unique Grain ID values "UniqueGrainIDVector" (of size "NumberOfGrains"), return a third vector "GrainSizeVector"
// listing the size of each of the "NumberOfGrains" grains (scaled by the cell size deltax and depending on whether the
// region is 1D, 2D, or 3D)
std::vector<float> getGrainSizes(const std::vector<int> GrainIDVector, const std::vector<int> UniqueGrainIDVector,
                                 const int NumberOfGrains, double deltax, std::string RegionType) {

    std::vector<float> GrainSizeVector(NumberOfGrains);
    double conv = convertToMicrons(deltax, RegionType);
    for (int n = 0; n < NumberOfGrains; n++) {
        int GrainSizeCells = std::count(GrainIDVector.begin(), GrainIDVector.end(), UniqueGrainIDVector[n]);
        // convert to either microns, square microns, or cubic microns
        GrainSizeVector[n] = conv * GrainSizeCells;
    }
    return GrainSizeVector;
}

// Given the 3D grain structure "GrainID", determine the extent in the direction specified of each of the
// "NumberOfGrains" unique grains from the volume bounded by [XLow, XHigh], [YLow, YHigh], [ZLow, ZHigh]. Extent is
// calculated in microns
void calcGrainExtent(std::vector<float> &GrainExtent, ViewI3D_H GrainID, const std::vector<int> UniqueGrainIDVector,
                     std::vector<float> GrainSizeVector, const int NumberOfGrains, const int XLow, const int XHigh,
                     const int YLow, const int YHigh, const int ZLow, const int ZHigh, std::string Direction,
                     double deltax, std::string RegionType) {

    for (int n = 0; n < NumberOfGrains; n++) {
        int ThisGrainID = UniqueGrainIDVector[n];
        int ThisGrainSize = std::round(GrainSizeVector[n] * convertToCells(deltax, RegionType));
        std::vector<int> GrainCoordinate(ThisGrainSize);
        int count = 0;
        for (int k = ZLow; k <= ZHigh; k++) {
            for (int i = XLow; i <= XHigh; i++) {
                for (int j = YLow; j <= YHigh; j++) {
                    if (GrainID(k, i, j) == ThisGrainID) {
                        if (Direction == "X")
                            GrainCoordinate[count] = i;
                        else if (Direction == "Y")
                            GrainCoordinate[count] = j;
                        else if (Direction == "Z")
                            GrainCoordinate[count] = k;
                        count++;
                    }
                }
            }
        }
        int MinCoord = *std::min_element(GrainCoordinate.begin(), GrainCoordinate.end());
        int MaxCoord = *std::max_element(GrainCoordinate.begin(), GrainCoordinate.end());
        GrainExtent[n] = (MaxCoord - MinCoord + 1) * convertToMicrons(deltax, "length");
    }
}

// Calculate misorientation relative to the specified cardinal direction for each grain and store in the vector
std::vector<float> getGrainMisorientation(std::string Direction, ViewF_H GrainUnitVector,
                                          std::vector<int> UniqueGrainIDVector, int NumberOfOrientations,
                                          int NumberOfGrains) {

    std::vector<float> GrainMisorientationVector(NumberOfGrains);
    int direction;
    if (Direction == "X")
        direction = 0;
    else if (Direction == "Y")
        direction = 1;
    else if (Direction == "Z")
        direction = 2;
    else
        throw std::runtime_error("Error: invalid direction specified in calcGrainMisorientation: should be X, Y, or Z");
    ViewF_H GrainMisorientation = MisorientationCalc(NumberOfOrientations, GrainUnitVector, direction);
    for (int n = 0; n < NumberOfGrains; n++) {
        int MyOrientation = getGrainOrientation(UniqueGrainIDVector[n], NumberOfOrientations);
        float MyMisorientation = GrainMisorientation(MyOrientation);
        GrainMisorientationVector[n] = MyMisorientation;
    }
    return GrainMisorientationVector;
}

// From the grain extents in x, y, and z, calcualte the aspect ratio for each grain in the build to the average of the
// transverse directions
void calcBuildTransAspectRatio(std::vector<float> &BuildTransAspectRatio, std::vector<float> GrainExtentX,
                               std::vector<float> GrainExtentY, std::vector<float> GrainExtentZ, int NumberOfGrains) {

    for (int n = 0; n < NumberOfGrains; n++) {
        float AR_XZ = GrainExtentZ[n] / GrainExtentX[n];
        float AR_YZ = GrainExtentZ[n] / GrainExtentY[n];
        BuildTransAspectRatio[n] = 0.5 * (AR_XZ + AR_YZ);
    }
}

// Determine the R, G, or B value for this grain orientation for the IPF-Z inverse pole figure colormap (Color = 0 for
// R, 1 for G, 2 for B)
std::vector<float> getIPFZColor(int Color, std::vector<int> UniqueGrainIDVector, int NumberOfOrientations,
                                ViewF_H GrainRGBValues, int NumberOfGrains) {
    std::vector<float> IPFZColor(NumberOfGrains);
    for (int n = 0; n < NumberOfGrains; n++) {
        int MyOrientation = getGrainOrientation(UniqueGrainIDVector[n], NumberOfOrientations);
        IPFZColor[n] = GrainRGBValues(3 * MyOrientation + Color);
    }
    return IPFZColor;
}
