// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
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
int FindTopOrBottom(ViewS3D_H LayerID, int XLow, int XHigh, int YLow, int YHigh, int nz, int L, std::string HighLow) {

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

// Reads and ignored ASCII field data
void ReadIgnoreASCIIField(std::ifstream &InputDataStream, int nx, int ny, int nz) {
    std::string line;
    for (int i = 0; i < nx * ny * nz; i++)
        getline(InputDataStream, line);
}

void InitializeData(std::string MicrostructureFile, int nx, int ny, int nz, ViewI3D_H &GrainID, ViewS3D_H &LayerID) {

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
        // This line says which variable appears next in the file, along with its type
        // A blank line ends the file read
        std::string read_datatype_string, read_fieldname;
        getline(InputDataStream, line);
        if (line.empty())
            break;
        else {
            // Thus far, every ExaCA field printed to VTK files has been one of these three types - ensure that this
            // field is one of them
            std::vector<std::string> possible_vtk_datatypes = {"short", "int", "float"};
            int num_possible_vtk_types = possible_vtk_datatypes.size();
            for (auto n = 0; n < num_possible_vtk_types; n++) {
                if (line.find(possible_vtk_datatypes[n]) != std::string::npos)
                    read_datatype_string = possible_vtk_datatypes[n];
            }
            if (read_datatype_string.empty())
                throw std::runtime_error("Error: Unknown data type for a field in the vtk file");

            // GrainID and LayerID are the only fields currently used by the analysis, other fields will be
            // skipped/ignored
            std::vector<std::string> possible_fieldnames = {"GrainID", "LayerID"};
            int num_possible_fieldnames = possible_fieldnames.size();
            for (auto n = 0; n < num_possible_fieldnames; n++) {
                if (line.find(possible_fieldnames[n]) != std::string::npos)
                    read_fieldname = possible_fieldnames[n];
            }
            if (read_fieldname.empty()) {
                std::cout << "Reading and ignoring field described by" << line << std::endl;
                if (Binaryvtk) {
                    if (read_datatype_string == "short")
                        ReadIgnoreBinaryField<short>(InputDataStream, nx, ny, nz);
                    else if (read_datatype_string == "int")
                        ReadIgnoreBinaryField<int>(InputDataStream, nx, ny, nz);
                    else if (read_datatype_string == "float")
                        ReadIgnoreBinaryField<float>(InputDataStream, nx, ny, nz);
                }
                else
                    ReadIgnoreASCIIField(InputDataStream, nx, ny, nz);
            }
            else {
                // Valid fieldname to be read
                // 1 more unused line
                getline(InputDataStream, line);
                // Place appropriate data
                if (read_fieldname == "GrainID") {
                    // GrainID data should be type int
                    if (read_datatype_string != "int")
                        throw std::runtime_error("Error: Field GrainID should be data of type int");
                    if (Binaryvtk)
                        GrainID = ReadBinaryField<ViewI3D_H, int>(InputDataStream, nx, ny, nz, "GrainID");
                    else
                        GrainID = ReadASCIIField<ViewI3D_H>(InputDataStream, nx, ny, nz, "GrainID");
                }
                else if (read_fieldname == "LayerID") {
                    // LayerID may be int or short, but is stored as type short
                    if ((read_datatype_string != "int") && (read_datatype_string != "short"))
                        throw std::runtime_error("Error: Field LayerID should be data of type int or short");
                    if (!Binaryvtk)
                        LayerID = ReadASCIIField<ViewS3D_H>(InputDataStream, nx, ny, nz, "LayerID");
                    else {
                        if (read_datatype_string == "int")
                            LayerID = ReadBinaryField<ViewS3D_H, int>(InputDataStream, nx, ny, nz, "LayerID");
                        else
                            LayerID = ReadBinaryField<ViewS3D_H, short>(InputDataStream, nx, ny, nz, "LayerID");
                    }
                }
                std::cout << "Data field " << read_fieldname << " read" << std::endl;
            }
        }
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

// Print the string "output" to both output file streams
void dual_print(std::string output, std::ostream &stream1, std::ostream &stream2) {
    stream1 << output;
    stream2 << output;
}
