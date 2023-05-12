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

// Parse log file using json format
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
        std::string pathtoorientations = RotationFilename.substr(0, baseorientation_startpos + 1);
        std::string baseorientationname =
            RotationFilename.substr(baseorientation_startpos + 1, endpos - baseorientation_startpos - 1);
        std::size_t customname_startpos = baseorientationname.find_last_of("_");
        std::string customname = baseorientationname.substr(customname_startpos + 1, endpos - customname_startpos - 1);
        EulerAnglesFilename = pathtoorientations + "GrainOrientationEulerAnglesBungeZXZ_" + customname + ".csv";
        RGBFilename = pathtoorientations + "GrainOrientationRGB_IPF-Z_" + customname + ".csv";
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
// listing the size of each of the "NumberOfGrains" grains (as either microns, square microns, or cubic microns
// depending on the dimensionality of the region type)
std::vector<float> getGrainSizes(const std::vector<int> GrainIDVector, const std::vector<int> UniqueGrainIDVector,
                                 const int NumberOfGrains, double deltax, std::string RegionType) {

    std::vector<float> GrainSizeVector_Microns(NumberOfGrains);
    double conv = convertToMicrons(deltax, RegionType);
    for (int n = 0; n < NumberOfGrains; n++) {
        int GrainSizeCells = std::count(GrainIDVector.begin(), GrainIDVector.end(), UniqueGrainIDVector[n]);
        // convert to either microns, square microns, or cubic microns
        GrainSizeVector_Microns[n] = conv * GrainSizeCells;
    }
    return GrainSizeVector_Microns;
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
