// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "GAutils.hpp"

#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>

// Parse log file using json format
void ParseLogFile(std::string LogFile, int &nx, int &ny, int &nz, double &deltax, int &NumberOfLayers,
                  std::vector<double> &XYZBounds, std::string &grain_unit_vector_file, bool OrientationFilesInInput) {

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
    grain_unit_vector_file = logdata["GrainOrientationFile"];
    InputDataStream.close();
}

// Reads and ignored ASCII field data
void ReadIgnoreASCIIField(std::ifstream &InputDataStream, int nx, int ny, int nz) {
    std::string line;
    for (int i = 0; i < nx * ny * nz; i++)
        getline(InputDataStream, line);
}

// Get the scaling factor to convert a number of cells into a length, area, or volume
double convertToMicrons(double deltax, std::string RegionType) {
    double CellToSizeScaling;
    if (RegionType == "volume")
        CellToSizeScaling = Kokkos::pow(deltax, 3) * Kokkos::pow(10, 18);
    else if (RegionType == "area")
        CellToSizeScaling = Kokkos::pow(deltax, 2) * Kokkos::pow(10, 12);
    else if (RegionType == "length")
        CellToSizeScaling = deltax * Kokkos::pow(10, 6);
    else
        throw std::runtime_error("Error: unknown region type");
    return CellToSizeScaling;
}

// Get the scaling factor to convert a length, area, or volume into a number of cells
double convertToCells(double deltax, std::string RegionType) {
    double SizeToCellScaling;
    if (RegionType == "volume")
        SizeToCellScaling = 1.0 / (pow(deltax, 3) * Kokkos::pow(10, 18));
    else if (RegionType == "area")
        SizeToCellScaling = 1.0 / (pow(deltax, 2) * Kokkos::pow(10, 12));
    else if (RegionType == "length")
        SizeToCellScaling = 1.0 / (deltax * Kokkos::pow(10, 6));
    else
        throw std::runtime_error("Error: unknown region type");
    return SizeToCellScaling;
}

// Print the string "output" to both output file streams
void dual_print(std::string output, std::ostream &stream1, std::ostream &stream2) {
    stream1 << output;
    stream2 << output;
}
