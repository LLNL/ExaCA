// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAinitialize.hpp"

#include "CAconfig.hpp"
#include "CAfunctions.hpp"
#include "CAparsefiles.hpp"

#include "mpi.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <regex>

// Initializes input parameters, mesh, temperature field, and grain structures for CA simulations

//*****************************************************************************/
// Read ExaCA input file.
void InputReadFromFile(int id, std::string InputFile, std::string &SimulationType, int &DecompositionStrategy,
                       double &AConst, double &BConst, double &CConst, double &DConst, double &FreezingRange,
                       double &deltax, double &NMax, double &dTN, double &dTsigma, std::string &OutputFile,
                       std::string &GrainOrientationFile, std::string &temppath, std::string &tempfile,
                       int &TempFilesInSeries, std::vector<std::string> &temp_paths, bool &ExtraWalls,
                       double &HT_deltax, bool &RemeltingYN, double &deltat, int &NumberOfLayers, int &LayerHeight,
                       std::string &SubstrateFileName, float &SubstrateGrainSpacing, bool &UseSubstrateFile, double &G,
                       double &R, int &nx, int &ny, int &nz, double &FractSurfaceSitesActive, std::string &PathToOutput,
                       int &PrintDebug, bool &PrintMisorientation, bool &PrintFullOutput, int &NSpotsX, int &NSpotsY,
                       int &SpotOffset, int &SpotRadius, bool &PrintTimeSeries, int &TimeSeriesInc,
                       bool &PrintIdleTimeSeriesFrames) {

    // For now, assuming no remelting
    RemeltingYN = false;

    // Required inputs that should be present in the input file, regardless of problem type
    std::vector<std::string> RequiredInputs_General = {
        "Decomposition strategy",                        // Required input 0
        "Material",                                      // Required input 1
        "Cell size",                                     // Required input 2
        "Heterogeneous nucleation density",              // Required input 3
        "Mean nucleation undercooling",                  // Required input 4
        "Standard deviation of nucleation undercooling", // Required input 5
        "Path to output",                                // Required input 6
        "Output file base name",                         // Required input 7
        "File of grain orientations",                    // Required input 8
        "grain misorientation values",                   // Required input 9
        "all ExaCA data"                                 // Required input 10
    };

    // Optional inputs that may be present in the input file, regardless of problem type
    std::vector<std::string> OptionalInputs_General = {
        "Debug check (reduced)",              // Optional input 0
        "Debug check (extensive)",            // Optional input 1
        "Print intermediate output frames",   // Optional input 2
        "separate frames",                    // Optional input 3
        "output even if system is unchanged", // Optional input 4
    };

    // Values used temporarily to store information from the file
    int NRatio = 0;
    float SpotOffsetFloat = 0.0;
    float SpotRadiusFloat = 0.0;
    double TimeSeriesFrameInc_time = 0.0;
    std::string MaterialName, GrainOrientationFile_Read;

    // Read the input file provided on the command line - determine which problem type is being simulated
    checkFileExists(InputFile, "Input", id);
    std::ifstream InputData;
    InputData.open(InputFile);
    skipLines(InputData, "*****"); // Skip lines until past the header lines above the asterisks
    // First line after the header must be the problem type: either R, C, or S
    // Additional required/optional inputs depending on problem type
    SimulationType = parseInput(InputData, "Problem type");
    std::vector<std::string> RequiredInputs_ProblemSpecific, OptionalInputs_ProblemSpecific;
    if (SimulationType == "C") {
        RequiredInputs_ProblemSpecific.resize(7);
        RequiredInputs_ProblemSpecific[0] = "Thermal gradient";
        RequiredInputs_ProblemSpecific[1] = "Cooling rate";
        RequiredInputs_ProblemSpecific[2] = "Velocity";
        RequiredInputs_ProblemSpecific[3] = "Domain size in x";
        RequiredInputs_ProblemSpecific[4] = "Domain size in y";
        RequiredInputs_ProblemSpecific[5] = "Domain size in z";
        RequiredInputs_ProblemSpecific[6] = "surface sites active";
    }
    else if (SimulationType == "S") {
        RequiredInputs_ProblemSpecific.resize(9);
        RequiredInputs_ProblemSpecific[0] = "Thermal gradient";
        RequiredInputs_ProblemSpecific[1] = "Cooling rate";
        RequiredInputs_ProblemSpecific[2] = "Velocity";
        RequiredInputs_ProblemSpecific[3] = "spots in x";
        RequiredInputs_ProblemSpecific[4] = "spots in y";
        RequiredInputs_ProblemSpecific[5] = "Offset between spot centers";
        RequiredInputs_ProblemSpecific[6] = "Radii of spots";
        RequiredInputs_ProblemSpecific[7] = "Number of layers";
        RequiredInputs_ProblemSpecific[8] = "Offset between layers";
        OptionalInputs_ProblemSpecific.resize(2);
        OptionalInputs_ProblemSpecific[0] = "Substrate grain spacing";
        OptionalInputs_ProblemSpecific[1] = "Substrate filename";
    }
    else if (SimulationType == "R") {
        RequiredInputs_ProblemSpecific.resize(6);
        RequiredInputs_ProblemSpecific[0] = "Time step";
        RequiredInputs_ProblemSpecific[1] = "Temperature filename";
        RequiredInputs_ProblemSpecific[2] = "Number of temperature files";
        RequiredInputs_ProblemSpecific[3] = "Number of layers";
        RequiredInputs_ProblemSpecific[4] = "Offset between layers";
        RequiredInputs_ProblemSpecific[5] = "Extra set of wall cells";
        OptionalInputs_ProblemSpecific.resize(4);
        OptionalInputs_ProblemSpecific[0] = "Substrate grain spacing";
        OptionalInputs_ProblemSpecific[1] = "Substrate filename";
        OptionalInputs_ProblemSpecific[2] = "Heat transport data mesh size";
        OptionalInputs_ProblemSpecific[3] = "Path to temperature file(s)";
    }
    else {
        std::string error = "Error: problem type must be C, S, or R: the value given was " + SimulationType;
        throw std::runtime_error(error);
    }
    int NumRequiredInputs_General = RequiredInputs_General.size();
    int NumOptionalInputs_General = OptionalInputs_General.size();
    int NumRequiredInputs_ProblemSpecific = RequiredInputs_ProblemSpecific.size();
    int NumOptionalInputs_ProblemSpecific = OptionalInputs_ProblemSpecific.size();
    std::vector<std::string> RequiredInputsRead_General(NumRequiredInputs_General),
        OptionalInputsRead_General(NumOptionalInputs_General),
        RequiredInputsRead_ProblemSpecific(NumRequiredInputs_ProblemSpecific),
        OptionalInputsRead_ProblemSpecific(NumOptionalInputs_ProblemSpecific);

    // Read the rest of the input file to initialize required and optional input parameters
    std::string line;
    while (std::getline(InputData, line)) {
        // Ignore lines with *, as these are not inputs
        if (line.find("***") == std::string::npos) {
            // Check if this is a general required input
            bool RequiredYN =
                parseInputFromList(line, RequiredInputs_General, RequiredInputsRead_General, NumRequiredInputs_General);
            if (!(RequiredYN)) {
                // If not a general required input, check if this is a problem type specific required input
                RequiredYN = parseInputFromList(line, RequiredInputs_ProblemSpecific,
                                                RequiredInputsRead_ProblemSpecific, NumRequiredInputs_ProblemSpecific);
            }
            if (!(RequiredYN)) {
                // Check if this is an general optional input
                bool OptionalYN = parseInputFromList(line, OptionalInputs_General, OptionalInputsRead_General,
                                                     NumOptionalInputs_General);
                if (!(OptionalYN)) {
                    // If not a general optional input, check if this is a problem type specific optional input
                    OptionalYN =
                        parseInputFromList(line, OptionalInputs_ProblemSpecific, OptionalInputsRead_ProblemSpecific,
                                           NumOptionalInputs_ProblemSpecific);
                }
                if (!(OptionalYN) && (id == 0)) {
                    std::cout << "WARNING: input " << line
                              << " did not match any optional nor required inputs known to ExaCA and will be ignored"
                              << std::endl;
                }
            }
        }
    }
    InputData.close();

    // Ensure that all required inputs were given - general and problem type specific
    for (int i = 0; i < NumRequiredInputs_General; i++) {
        if (RequiredInputsRead_General[i].empty()) {
            std::string error =
                "Error: Required input " + RequiredInputs_General[i] + " was not present in the input file";
            throw std::runtime_error(error);
        }
    }
    for (int i = 0; i < NumRequiredInputs_ProblemSpecific; i++) {
        if (RequiredInputsRead_ProblemSpecific[i].empty()) {
            std::string error =
                "Error: Required input " + RequiredInputs_ProblemSpecific[i] + " was not present in the input file";
            throw std::runtime_error(error);
        }
    }

    // Convert information read from the file into values usable by ExaCA
    // Required inputs for all problems
    DecompositionStrategy = getInputInt(RequiredInputsRead_General[0]);
    MaterialName = RequiredInputsRead_General[1];
    deltax = getInputDouble(RequiredInputsRead_General[2], -6);
    NMax = getInputDouble(RequiredInputsRead_General[3], 12);
    dTN = getInputDouble(RequiredInputsRead_General[4]);
    dTsigma = getInputDouble(RequiredInputsRead_General[5]);
    PathToOutput = RequiredInputsRead_General[6];
    OutputFile = RequiredInputsRead_General[7];
    GrainOrientationFile_Read = RequiredInputsRead_General[8];
    PrintMisorientation = getInputBool(RequiredInputsRead_General[9]);
    PrintFullOutput = getInputBool(RequiredInputsRead_General[10]);
    // Problem type-specific inputs
    if (SimulationType == "R") {
        deltat = getInputDouble(RequiredInputsRead_ProblemSpecific[0], -6);
        tempfile = RequiredInputsRead_ProblemSpecific[1];
        TempFilesInSeries = getInputInt(RequiredInputsRead_ProblemSpecific[2]);
        NumberOfLayers = getInputInt(RequiredInputsRead_ProblemSpecific[3]);
        LayerHeight = getInputInt(RequiredInputsRead_ProblemSpecific[4]);
        ExtraWalls = getInputBool(RequiredInputsRead_ProblemSpecific[5]);
        if (id == 0) {
            std::cout << "CA Simulation using temperature data from file(s)" << std::endl;
            std::cout << "The time step is " << deltat << " seconds" << std::endl;
            if (TempFilesInSeries > 1)
                std::cout << "Temperature data files are *" << tempfile << " , and there are " << TempFilesInSeries
                          << " in the series" << std::endl;
            else
                std::cout << "The temperature data file is " << tempfile << std::endl;
            std::cout << "A total of " << NumberOfLayers << " layers of solidification offset by " << LayerHeight
                      << " CA cells will be simulated" << std::endl;
        }
    }
    else if (SimulationType == "C") {
        G = getInputDouble(RequiredInputsRead_ProblemSpecific[0]);
        R = getInputDouble(RequiredInputsRead_ProblemSpecific[1]);
        NRatio = getInputInt(RequiredInputsRead_ProblemSpecific[2]);
        nx = getInputInt(RequiredInputsRead_ProblemSpecific[3]);
        nx = nx + 2; // Domain size in x, wall cells at boundaries
        ny = getInputInt(RequiredInputsRead_ProblemSpecific[4]);
        ny = ny + 2; // Domain size in y, wall cells at boundaries
        nz = getInputInt(RequiredInputsRead_ProblemSpecific[5]);
        nz = nz + 1; // Domain size in z, wall cells at Z = 0, but not Z = nz-1
        FractSurfaceSitesActive = getInputDouble(RequiredInputsRead_ProblemSpecific[6]);
        deltat = deltax / ((double)(NRatio) * (R / G));
        NumberOfLayers = 1;
        LayerHeight = nz;
        if (id == 0) {
            std::cout << "CA Simulation using a unidirectional, fixed thermal gradient of " << G
                      << " K/m and a cooling rate of " << R << " K/s" << std::endl;
            std::cout << "The time step is " << deltat * pow(10, 6) << " microseconds" << std::endl;
            std::cout << "The fraction of CA cells at the bottom surface that are active is " << FractSurfaceSitesActive
                      << std::endl;
        }
    }
    else if (SimulationType == "S") {
        G = getInputFloat(RequiredInputsRead_ProblemSpecific[0]);
        R = getInputFloat(RequiredInputsRead_ProblemSpecific[1]);
        NRatio = getInputInt(RequiredInputsRead_ProblemSpecific[2]);
        NSpotsX = getInputInt(RequiredInputsRead_ProblemSpecific[3]);
        NSpotsY = getInputInt(RequiredInputsRead_ProblemSpecific[4]);
        SpotOffsetFloat = getInputFloat(RequiredInputsRead_ProblemSpecific[5]);
        SpotRadiusFloat = getInputFloat(RequiredInputsRead_ProblemSpecific[6]);
        NumberOfLayers = getInputInt(RequiredInputsRead_ProblemSpecific[7]);
        LayerHeight = getInputInt(RequiredInputsRead_ProblemSpecific[8]);
        deltat = deltax / (NRatio * (R / G));
        SpotOffset = SpotOffsetFloat * 1.0 * pow(10, -6) / deltax;
        SpotRadius = SpotRadiusFloat * 1.0 * pow(10, -6) / deltax;
        // Calculate nx, ny, and nz based on spot array pattern and number of layers
        // Z = 0 consists of wall cells, Z = 1 is active cells just outside of the melt pool footprint
        // Z = 2 is the lowest Z coordinate with the first layers spot melt data
        nz = SpotRadius + 2 + (NumberOfLayers - 1) * LayerHeight;
        nx = 4 + 2 * SpotRadius + SpotOffset * (NSpotsX - 1);
        ny = 4 + 2 * SpotRadius + SpotOffset * (NSpotsY - 1);
        if (id == 0) {
            std::cout << "CA Simulation using a radial, fixed thermal gradient of " << G
                      << " K/m as a series of hemispherical spots, and a cooling rate of " << R << " K/s" << std::endl;
            std::cout << "A total of " << NumberOfLayers << " spots per layer, with layers offset by " << LayerHeight
                      << " CA cells will be simulated" << std::endl;
            std::cout << "The time step is " << deltat * pow(10, 6) << " microseconds" << std::endl;
        }
    }
    // Optional inputs - should files post-initialization be printed for debugging?
    bool PrintDebugA = false;
    bool PrintDebugB = false;
    if (!(OptionalInputsRead_General[0].empty()))
        PrintDebugA = getInputBool(OptionalInputsRead_General[0]);
    if (!(OptionalInputsRead_General[1].empty()))
        PrintDebugB = getInputBool(OptionalInputsRead_General[1]);
    if ((OptionalInputsRead_General[0].empty()) && (OptionalInputsRead_General[1].empty()))
        PrintDebug = 0;
    if (PrintDebugB)
        PrintDebug = 2;
    else if (PrintDebugA)
        PrintDebug = 1;
    else
        PrintDebug = 0;
    // Should intermediate output be printed?
    if (!(OptionalInputsRead_General[2].empty())) {
        PrintTimeSeries = getInputBool(OptionalInputsRead_General[2]);
        if (PrintTimeSeries) {
            if (!(OptionalInputsRead_General[3].empty())) {
                TimeSeriesFrameInc_time = getInputFloat(OptionalInputsRead_General[3], -6);
                TimeSeriesInc = round(TimeSeriesFrameInc_time / deltat);
                if (!(OptionalInputsRead_General[4].empty()))
                    PrintIdleTimeSeriesFrames = getInputBool(OptionalInputsRead_General[4]);
                else
                    PrintIdleTimeSeriesFrames = false;
            }
            else
                throw std::runtime_error(
                    "Error: Cannot print intermediate output data without a specified output increment");
        }
        else {
            PrintTimeSeries = false;
            PrintIdleTimeSeriesFrames = false;
        }
    }
    else {
        PrintTimeSeries = false;
        PrintIdleTimeSeriesFrames = false;
    }
    // For simulations with substrate grain structures, should an input grain spacing or a substrate file be used?
    if ((SimulationType == "S") || (SimulationType == "R")) {
        // Exactly one of the two inputs "sub grain size" and "sub filename" should be present
        if ((!(OptionalInputsRead_ProblemSpecific[0].empty())) && (!(OptionalInputsRead_ProblemSpecific[1].empty())))
            throw std::runtime_error("Error: only one of substrate grain size and substrate structure filename should "
                                     "be provided in the input file");
        else if ((OptionalInputsRead_ProblemSpecific[0].empty()) && (OptionalInputsRead_ProblemSpecific[1].empty()))
            throw std::runtime_error(
                "Error: neither substrate grain size nor substrate structure filename was provided in the input file");
        else {
            if (!(OptionalInputsRead_ProblemSpecific[0].empty())) {
                UseSubstrateFile = false;
                SubstrateGrainSpacing = getInputFloat(OptionalInputsRead_ProblemSpecific[0]);
            }
            else {
                // Check that substrate file exists
                UseSubstrateFile = true;
                SubstrateFileName = OptionalInputsRead_ProblemSpecific[1];
                checkFileExists(SubstrateFileName, "Substrate", id);
            }
        }
    }
    // Input temperature data spacing or path to temperature data may be given
    // If not given, it is assumed HT_deltax = CA cell size and temperature data is located in examples/Temperatures
    if (SimulationType == "R") {
        if (OptionalInputsRead_ProblemSpecific[2].empty())
            HT_deltax = deltax;
        else
            HT_deltax = getInputDouble(OptionalInputsRead_ProblemSpecific[2], -6);

        if (OptionalInputsRead_ProblemSpecific[3].empty())
            temppath = "examples/Temperatures";
        else
            temppath = OptionalInputsRead_ProblemSpecific[3];

        // Check that temperature file(s) exist
        temp_paths.resize(TempFilesInSeries, "");
        getTemperatureFilePaths(temppath, tempfile, temp_paths);
        for (int i = 0; i < TempFilesInSeries; i++) {
            if (id == 0)
                std::cout << "Checking file " << temp_paths[i] << std::endl;
            checkFileExists(temp_paths[i], "Temperature", id);
        }
    }

    // Path to file of materials constants based on install/source location
    std::string MaterialFile = checkFileInstalled(MaterialName, "Materials", id);
    checkFileNotEmpty(MaterialFile);
    // Read material file (specified from main input file) to obtain values for A, B, C, and D for the interfacial
    // reponse function
    parseMaterialFile(MaterialFile, AConst, BConst, CConst, DConst, FreezingRange);

    // Path to file of grain orientations based on install/source location
    GrainOrientationFile = checkFileInstalled(GrainOrientationFile_Read, "Substrate", id);
    checkFileNotEmpty(GrainOrientationFile);

    if (id == 0) {
        std::cout << "Decomposition Strategy is " << DecompositionStrategy << std::endl;
        std::cout << "Material simulated is " << MaterialName
                  << ", interfacial response function constants are A = " << AConst << ", B = " << BConst
                  << ", C = " << CConst << ", and D = " << DConst << std::endl;
        std::cout << "CA cell size is " << deltax * pow(10, 6) << " microns" << std::endl;
        std::cout << "Nucleation density is " << NMax << " per m^3" << std::endl;
        std::cout << "Mean nucleation undercooling is " << dTN << " K, standard deviation of distribution is "
                  << dTsigma << "K" << std::endl;
        if (PrintTimeSeries)
            std::cout << "Intermediate output for movie frames will be printed every " << TimeSeriesInc
                      << " time steps (or every " << TimeSeriesFrameInc_time << " microseconds)" << std::endl;
    }
}

//*****************************************************************************/
// Intialize neighbor list structures (NeighborX, NeighborY, NeighborZ, and ItList)
void NeighborListInit(ViewI_H NeighborX, ViewI_H NeighborY, ViewI_H NeighborZ, ViewI2D_H ItList) {

    // Assignment of neighbors around a cell "X" is as follows (in order of closest to furthest from cell "X")
    // Neighbors 0 through 8 are in the -Y direction
    // Neighbors 9 through 16 are in the XY plane with cell X
    // Neighbors 17 through 25 are in the +Y direction

    NeighborX(0) = 0;
    NeighborY(0) = -1;
    NeighborZ(0) = 0;
    NeighborX(1) = 1;
    NeighborY(1) = -1;
    NeighborZ(1) = 0;
    NeighborX(2) = -1;
    NeighborY(2) = -1;
    NeighborZ(2) = 0;
    NeighborX(3) = 0;
    NeighborY(3) = -1;
    NeighborZ(3) = 1;
    NeighborX(4) = 0;
    NeighborY(4) = -1;
    NeighborZ(4) = -1;
    NeighborX(5) = -1;
    NeighborY(5) = -1;
    NeighborZ(5) = 1;
    NeighborX(6) = 1;
    NeighborY(6) = -1;
    NeighborZ(6) = 1;
    NeighborX(7) = -1;
    NeighborY(7) = -1;
    NeighborZ(7) = -1;
    NeighborX(8) = 1;
    NeighborY(8) = -1;
    NeighborZ(8) = -1;

    NeighborX(9) = 0;
    NeighborY(9) = 0;
    NeighborZ(9) = 1;
    NeighborX(10) = 0;
    NeighborY(10) = 0;
    NeighborZ(10) = -1;
    NeighborX(11) = 1;
    NeighborY(11) = 0;
    NeighborZ(11) = 1;
    NeighborX(12) = -1;
    NeighborY(12) = 0;
    NeighborZ(12) = 1;
    NeighborX(13) = 1;
    NeighborY(13) = 0;
    NeighborZ(13) = -1;
    NeighborX(14) = -1;
    NeighborY(14) = 0;
    NeighborZ(14) = -1;
    NeighborX(15) = 1;
    NeighborY(15) = 0;
    NeighborZ(15) = 0;
    NeighborX(16) = -1;
    NeighborY(16) = 0;
    NeighborZ(16) = 0;

    NeighborX(17) = 0;
    NeighborY(17) = 1;
    NeighborZ(17) = 0;
    NeighborX(18) = 1;
    NeighborY(18) = 1;
    NeighborZ(18) = 0;
    NeighborX(19) = -1;
    NeighborY(19) = 1;
    NeighborZ(19) = 0;
    NeighborX(20) = 0;
    NeighborY(20) = 1;
    NeighborZ(20) = 1;
    NeighborX(21) = 0;
    NeighborY(21) = 1;
    NeighborZ(21) = -1;
    NeighborX(22) = 1;
    NeighborY(22) = 1;
    NeighborZ(22) = 1;
    NeighborX(23) = -1;
    NeighborY(23) = 1;
    NeighborZ(23) = 1;
    NeighborX(24) = 1;
    NeighborY(24) = 1;
    NeighborZ(24) = -1;
    NeighborX(25) = -1;
    NeighborY(25) = 1;
    NeighborZ(25) = -1;

    // If X and Y coordinates are not on edges, Case 0: iteratation over neighbors 0-25 possible
    for (int i = 0; i <= 25; i++) {
        ItList(0, i) = i;
    }
    // If Y coordinate is on lower edge, Case 1: iteration over only neighbors 9-25 possible
    int Pos = 0;
    for (int i = 0; i <= 25; i++) {
        if (NeighborY(i) != -1) {
            ItList(1, Pos) = i;
            Pos++;
        }
    }
    // If Y coordinate is on upper edge, Case 2: iteration over only neighbors 0-16 possible
    Pos = 0;
    for (int i = 0; i <= 25; i++) {
        if (NeighborY(i) != 1) {
            ItList(2, Pos) = i;
            Pos++;
        }
    }
    // If X coordinate is on lower edge, Case 3: iteration over only neighbors
    // 0,1,3,4,6,8,9,10,11,13,15,17,18,20,21,22,24
    Pos = 0;
    for (int i = 0; i <= 25; i++) {
        if (NeighborX(i) != -1) {
            ItList(3, Pos) = i;
            Pos++;
        }
    }
    // If X coordinate is on upper edge, Case 4: iteration over only neighbors
    // 0,2,3,4,5,7,9,10,12,14,16,17,19,20,21,23,25
    Pos = 0;
    for (int i = 0; i <= 25; i++) {
        if (NeighborX(i) != 1) {
            ItList(4, Pos) = i;
            Pos++;
        }
    }
    // If X/Y coordinates are on lower edge, Case 5: iteration over only neighbors 9,10,11,13,15,17,18,20,21,22,24
    Pos = 0;
    for (int i = 0; i <= 25; i++) {
        if ((NeighborX(i) != -1) && (NeighborY(i) != -1)) {
            ItList(5, Pos) = i;
            Pos++;
        }
    }
    // If X coordinate is on upper edge/Y on lower edge, Case 6:
    Pos = 0;
    for (int i = 0; i <= 25; i++) {
        if ((NeighborX(i) != 1) && (NeighborY(i) != -1)) {
            ItList(6, Pos) = i;
            Pos++;
        }
    }
    // If X coordinate is on lower edge/Y on upper edge, Case 7:
    Pos = 0;
    for (int i = 0; i <= 25; i++) {
        if ((NeighborX(i) != -1) && (NeighborY(i) != 1)) {
            ItList(7, Pos) = i;
            Pos++;
        }
    }
    // If X/Y coordinates are on upper edge, Case 8:
    Pos = 0;
    for (int i = 0; i <= 25; i++) {
        if ((NeighborX(i) != 1) && (NeighborY(i) != 1)) {
            ItList(8, Pos) = i;
            Pos++;
        }
    }
}

// Obtain the physical XYZ bounds of the domain, using either domain size from the input file, or reading temperature
// data files and parsing the coordinates
void FindXYZBounds(std::string SimulationType, int id, double &deltax, int &nx, int &ny, int &nz,
                   std::vector<std::string> &temp_paths, float &XMin, float &XMax, float &YMin, float &YMax,
                   float &ZMin, float &ZMax, int &LayerHeight, int NumberOfLayers, int TempFilesInSeries,
                   float *ZMinLayer, float *ZMaxLayer) {

    if (SimulationType == "R") {
        // Two passes through reading temperature data files- the first pass only reads the headers to
        // determine units and X/Y/Z bounds of the simulaton domain. Using the X/Y/Z bounds of the simulation domain,
        // nx, ny, and nz can be calculated and the domain decomposed among MPI processes. The maximum number of
        // remelting events in the simulation can also be calculated. The second pass reads the actual X/Y/Z/liquidus
        // time/cooling rate data and each rank stores the data relevant to itself in "RawData" - this is done in the
        // subroutine "ReadTemperatureData"
        XMin = 1000000.0;
        YMin = 1000000.0;
        ZMin = 1000000.0;
        XMax = -1000000.0;
        YMax = -1000000.0;
        ZMax = -1000000.0;

        // Read the first temperature file, first line to determine if the "new" OpenFOAM output format (with a 1 line
        // header) is used, or whether the "old" OpenFOAM header (which contains information like the X/Y/Z bounds of
        // the simulation domain) is
        std::ifstream FirstTemperatureFile;
        FirstTemperatureFile.open(temp_paths[0]);
        std::string FirstLineFirstFile;
        getline(FirstTemperatureFile, FirstLineFirstFile);
        std::size_t found = FirstLineFirstFile.find("Number of temperature data points");
        if (found != std::string::npos) {
            // Old temperature data format detected - no longer supported by ExaCA
            std::string error = "Error: Old header and temperature file format no longer supported";
            throw std::runtime_error(error);
        }

        // Read all data files to determine the domain bounds, max number of remelting events
        // for simulations with remelting
        int LayersToRead = std::min(NumberOfLayers, TempFilesInSeries); // was given in input file
        for (int LayerReadCount = 1; LayerReadCount <= LayersToRead; LayerReadCount++) {

            std::string tempfile_thislayer = temp_paths[LayerReadCount - 1];
            std::ifstream TemperatureFile;
            TemperatureFile.open(tempfile_thislayer);

            // Read the header line data
            // Make sure the first line contains all required column names: x, y, z, tm, tl, cr
            // Check 2 mixes of lower and uppercase letters/similar possible column names just in case
            std::string HeaderLine;
            getline(TemperatureFile, HeaderLine);
            checkForHeaderValues(HeaderLine);

            double XMin_ThisLayer = 1000000.0;
            double XMax_ThisLayer = -1000000.0;
            double YMin_ThisLayer = 1000000.0;
            double YMax_ThisLayer = -1000000.0;
            double ZMin_ThisLayer = 1000000.0;
            double ZMax_ThisLayer = -100000.0;

            // Units are assumed to be in meters, meters, seconds, seconds, and K/second
            std::vector<double> XCoordinates(1000000), YCoordinates(1000000), ZCoordinates(1000000);
            long unsigned int XYZPointCounter = 0;
            while (!TemperatureFile.eof()) {
                std::string s;
                getline(TemperatureFile, s);
                if (s.empty())
                    break;
                // Only get x, y, and z values
                std::vector<double> XYZTemperaturePoint(3, 0);
                getTemperatureDataPoint(s, XYZTemperaturePoint);
                XCoordinates[XYZPointCounter] = XYZTemperaturePoint[0];
                YCoordinates[XYZPointCounter] = XYZTemperaturePoint[1];
                ZCoordinates[XYZPointCounter] = XYZTemperaturePoint[2];
                XYZPointCounter++;
                if (XYZPointCounter == XCoordinates.size()) {
                    XCoordinates.resize(XYZPointCounter + 1000000);
                    YCoordinates.resize(XYZPointCounter + 1000000);
                    ZCoordinates.resize(XYZPointCounter + 1000000);
                }
            }
            XCoordinates.resize(XYZPointCounter);
            YCoordinates.resize(XYZPointCounter);
            ZCoordinates.resize(XYZPointCounter);
            TemperatureFile.close();
            // Min/max x, y, and z coordinates from this layer's data
            XMin_ThisLayer = *min_element(XCoordinates.begin(), XCoordinates.end());
            YMin_ThisLayer = *min_element(YCoordinates.begin(), YCoordinates.end());
            ZMin_ThisLayer = *min_element(ZCoordinates.begin(), ZCoordinates.end());
            XMax_ThisLayer = *max_element(XCoordinates.begin(), XCoordinates.end());
            YMax_ThisLayer = *max_element(YCoordinates.begin(), YCoordinates.end());
            ZMax_ThisLayer = *max_element(ZCoordinates.begin(), ZCoordinates.end());

            // Based on the input file's layer offset, adjust ZMin/ZMax from the temperature data coordinate
            // system to the multilayer CA coordinate system Check to see in the XYZ bounds for this layer are
            // also limiting for the entire multilayer CA coordinate system
            ZMin_ThisLayer += deltax * LayerHeight * (LayerReadCount - 1);
            ZMax_ThisLayer += deltax * LayerHeight * (LayerReadCount - 1);
            if (XMin_ThisLayer < XMin)
                XMin = XMin_ThisLayer;
            if (XMax_ThisLayer > XMax)
                XMax = XMax_ThisLayer;
            if (YMin_ThisLayer < YMin)
                YMin = YMin_ThisLayer;
            if (YMax_ThisLayer > YMax)
                YMax = YMax_ThisLayer;
            if (ZMin_ThisLayer < ZMin)
                ZMin = ZMin_ThisLayer;
            if (ZMax_ThisLayer > ZMax)
                ZMax = ZMax_ThisLayer;
            ZMinLayer[LayerReadCount - 1] = ZMin_ThisLayer;
            ZMaxLayer[LayerReadCount - 1] = ZMax_ThisLayer;
            if (id == 0)
                std::cout << "Layer = " << LayerReadCount << " Z Bounds are " << ZMin_ThisLayer << " " << ZMax_ThisLayer
                          << std::endl;
        }
        // Extend domain in Z (build) direction if the number of layers are simulated is greater than the number
        // of temperature files read
        if (NumberOfLayers > TempFilesInSeries) {
            for (int LayerReadCount = TempFilesInSeries; LayerReadCount < NumberOfLayers; LayerReadCount++) {
                if (TempFilesInSeries == 1) {
                    // Only one temperature file was read, so the upper Z bound should account for an additional
                    // "NumberOfLayers-1" worth of data Since all layers have the same temperature data, each
                    // layer's "ZMinLayer" is just translated from that of the first layer
                    ZMinLayer[LayerReadCount] = ZMinLayer[LayerReadCount - 1] + deltax * LayerHeight;
                    ZMaxLayer[LayerReadCount] = ZMaxLayer[LayerReadCount - 1] + deltax * LayerHeight;
                    ZMax += deltax * LayerHeight;
                }
                else {
                    // "TempFilesInSeries" temperature files was read, so the upper Z bound should account for
                    // an additional "NumberOfLayers-TempFilesInSeries" worth of data
                    int RepeatedFile = (LayerReadCount) % TempFilesInSeries;
                    int RepeatUnit = LayerReadCount / TempFilesInSeries;
                    ZMinLayer[LayerReadCount] =
                        ZMinLayer[RepeatedFile] + RepeatUnit * TempFilesInSeries * deltax * LayerHeight;
                    ZMaxLayer[LayerReadCount] =
                        ZMaxLayer[RepeatedFile] + RepeatUnit * TempFilesInSeries * deltax * LayerHeight;
                    ZMax += deltax * LayerHeight;
                }
            }
        }

        // Now at the conclusion of "Loop 0", the decomposition can be performed as the domain bounds are known
        // (all header lines from all files have been read)
        // CA nodes in each direction: (+2 for wall cells at the +/- X/Y boundaries
        // +2 for active/wall cells at the -Z boundary, no wall at the top Z boundary
        nx = round((XMax - XMin) / deltax) + 1 + 4;
        ny = round((YMax - YMin) / deltax) + 1 + 4;
        nz = round((ZMax - ZMin) / deltax) + 1 + 2;
    }
    else {
        // Using fixed G/R values to set temperature field - no temperature data to read
        // Setting physical domain bounds for consistency with problems that use input temperature data
        // Let XMin/YMin/ZMin be equal to 0, XMax/YMax/ZMax equal to nx*deltax, ny*deltax, nz*deltax
        XMin = 0.0;
        YMin = 0.0;
        ZMin = 0.0;
        XMax = nx * deltax;
        YMax = ny * deltax;
        ZMax = nx * deltax;
    }
    if (id == 0) {
        std::cout << "Domain size: " << nx << " by " << ny << " by " << nz << std::endl;
        std::cout << "X Limits of domain: " << XMin << " and " << XMax << std::endl;
        std::cout << "Y Limits of domain: " << YMin << " and " << YMax << std::endl;
        std::cout << "Z Limits of domain: " << ZMin << " and " << ZMax << std::endl;
        std::cout << "================================================================" << std::endl;
    }
}

// Decompose the domain into subdomains on each MPI rank: Each subdomain contains "MyXSlices" cells in X, and
// "MyYSlices" in Y. Each subdomain is offset from the full domain origin by "MyXOffset" cells in X, and "MyYOffset"
// cells in Y
void DomainDecomposition(int DecompositionStrategy, int id, int np, int &MyXSlices, int &MyYSlices, int &MyXOffset,
                         int &MyYOffset, int &NeighborRank_North, int &NeighborRank_South, int &NeighborRank_East,
                         int &NeighborRank_West, int &NeighborRank_NorthEast, int &NeighborRank_NorthWest,
                         int &NeighborRank_SouthEast, int &NeighborRank_SouthWest, int &nx, int &ny, int &nz,
                         int &ProcessorsInXDirection, int &ProcessorsInYDirection, long int &LocalDomainSize) {

    InitialDecomposition(DecompositionStrategy, nx, ny, ProcessorsInXDirection, ProcessorsInYDirection, id, np,
                         NeighborRank_North, NeighborRank_South, NeighborRank_East, NeighborRank_West,
                         NeighborRank_NorthEast, NeighborRank_NorthWest, NeighborRank_SouthEast,
                         NeighborRank_SouthWest);

    MyXOffset = XOffsetCalc(id, nx, ProcessorsInXDirection, ProcessorsInYDirection, DecompositionStrategy);
    MyXSlices = XMPSlicesCalc(id, nx, ProcessorsInXDirection, ProcessorsInYDirection, DecompositionStrategy);

    MyYOffset = YOffsetCalc(id, ny, ProcessorsInYDirection, np, DecompositionStrategy);
    MyYSlices = YMPSlicesCalc(id, ny, ProcessorsInYDirection, np, DecompositionStrategy);

    LocalDomainSize = MyXSlices * MyYSlices * nz; // Number of cells on this MPI rank
}

// Read in temperature data from files, stored in "RawData", with the appropriate MPI ranks storing the appropriate data
void ReadTemperatureData(double deltax, double HT_deltax, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset,
                         int nx, int ny, float XMin, float YMin, std::vector<std::string> &temp_paths,
                         int NumberOfLayers, int TempFilesInSeries, unsigned int &NumberOfTemperatureDataPoints,
                         std::vector<double> &RawData, int *FirstValue, int *LastValue) {

    int HTtoCAratio = HT_deltax / deltax; // OpenFOAM/CA cell size ratio
    int LowerXBound, LowerYBound, UpperXBound, UpperYBound;
    if (MyXOffset <= 2)
        LowerXBound = 2;
    else
        LowerXBound = MyXOffset - ((MyXOffset - 2) % HTtoCAratio);
    if (MyYOffset <= 2)
        LowerYBound = 2;
    else
        LowerYBound = MyYOffset - ((MyYOffset - 2) % HTtoCAratio);

    if (MyXOffset + MyXSlices - 1 >= nx - 3)
        UpperXBound = nx - 3;
    else
        UpperXBound = MyXOffset + MyXSlices - 1 + HTtoCAratio - ((MyXOffset + (MyXSlices - 1) - 2) % HTtoCAratio);
    if (MyYOffset + MyYSlices - 1 >= ny - 3)
        UpperYBound = ny - 3;
    else
        UpperYBound = MyYOffset + MyYSlices - 1 + HTtoCAratio - ((MyYOffset + (MyYSlices - 1) - 2) % HTtoCAratio);

    // Store raw data relevant to each rank in the vector structure RawData
    // With row/col 0 being wall cells and row/col 1 being solid cells outside of the melted area, the
    // domain starts at row/col 2 As the wall cells are not part of the physical domain (solid cells at
    // row/col 1 are defined as X = Y = 0, the melted region domain data starts at X = Y = deltax, with
    // data points at X or Y = deltax + N*HT_deltax through X or Y = nx-3 or ny-3

    // The X and Y bounds are the region (for this MPI rank) of the physical domain that needs to be
    // read extends past the actual spatial extent of the local domain for purposes of interpolating
    // from HT_deltax to deltax
    // Two passes through reading temperature data files- this is the second pass, reading the actual X/Y/Z/liquidus
    // time/cooling rate data and each rank stores the data relevant to itself in "RawData". With remelting
    // (SimulationType == "RM"), this is the same except that some X/Y/Z coordinates may be repeated in a file, and
    // a "melting time" value is stored in addition to liquidus time and cooling rate
    // Second pass through the files - ignore header line
    int LayersToRead = std::min(NumberOfLayers, TempFilesInSeries); // was given in input file
    for (int LayerReadCount = 1; LayerReadCount <= LayersToRead; LayerReadCount++) {

        std::string tempfile_thislayer = temp_paths[LayerReadCount - 1];
        std::ifstream TemperatureFile;
        TemperatureFile.open(tempfile_thislayer);
        FirstValue[LayerReadCount - 1] = NumberOfTemperatureDataPoints;

        std::string DummyLine;
        // ignore header line
        getline(TemperatureFile, DummyLine);

        // Read data from the remaining lines - values should be separated by commas
        // Space separated data is no longer accepted by ExaCA
        while (!TemperatureFile.eof()) {
            std::string s;
            getline(TemperatureFile, s);
            if (s.empty())
                break;
            // Get x, y, z, melting, liquidus, and cooling rate values
            std::vector<double> XYZTemperaturePoint(6, 0);
            getTemperatureDataPoint(s, XYZTemperaturePoint);

            // Check the CA grid positions of the data point to see which rank(s) should store it
            int XInt = 0, YInt = 0;
            XInt = round((XYZTemperaturePoint[0] - XMin) / deltax) + 2;
            YInt = round((XYZTemperaturePoint[1] - YMin) / deltax) + 2;
            if ((XInt >= LowerXBound) && (XInt <= UpperXBound) && (YInt >= LowerYBound) && (YInt <= UpperYBound)) {
                // This data point is inside the bounds of interest for this MPI rank - store inside of RawData
                RawData[NumberOfTemperatureDataPoints] = XYZTemperaturePoint[0];
                NumberOfTemperatureDataPoints++;
                RawData[NumberOfTemperatureDataPoints] = XYZTemperaturePoint[1];
                NumberOfTemperatureDataPoints++;
                RawData[NumberOfTemperatureDataPoints] = XYZTemperaturePoint[2];
                NumberOfTemperatureDataPoints++;
                RawData[NumberOfTemperatureDataPoints] = XYZTemperaturePoint[3];
                NumberOfTemperatureDataPoints++;
                RawData[NumberOfTemperatureDataPoints] = XYZTemperaturePoint[4];
                NumberOfTemperatureDataPoints++;
                RawData[NumberOfTemperatureDataPoints] = XYZTemperaturePoint[5];
                NumberOfTemperatureDataPoints++;
                if (NumberOfTemperatureDataPoints >= RawData.size() - 6) {
                    int OldSize = RawData.size();
                    RawData.resize(OldSize + 1000000);
                }
            }
        }
        LastValue[LayerReadCount - 1] = NumberOfTemperatureDataPoints;
    } // End loop over all files read for all layers
    RawData.resize(NumberOfTemperatureDataPoints);
    // Determine start values for each layer's data within "RawData"
    if (NumberOfLayers > TempFilesInSeries) {
        for (int LayerReadCount = TempFilesInSeries; LayerReadCount < NumberOfLayers; LayerReadCount++) {
            if (TempFilesInSeries == 1) {
                // Since all layers have the same temperature data, each layer's "ZMinLayer" is just
                // translated from that of the first layer
                FirstValue[LayerReadCount] = FirstValue[LayerReadCount - 1];
                LastValue[LayerReadCount] = LastValue[LayerReadCount - 1];
            }
            else {
                // All layers have different temperature data but in a repeating pattern
                int RepeatedFile = (LayerReadCount) % TempFilesInSeries;
                FirstValue[LayerReadCount] = FirstValue[RepeatedFile];
                LastValue[LayerReadCount] = LastValue[RepeatedFile];
            }
        }
    }
}

//*****************************************************************************/
// Initialize temperature data for a constrained solidification test problem
void TempInit_DirSolidification(double G, double R, int id, int &MyXSlices, int &MyYSlices, int &MyXOffset,
                                int &MyYOffset, double deltax, double deltat, int nx, int ny, int nz,
                                ViewI_H CritTimeStep, ViewF_H UndercoolingChange, ViewF_H UndercoolingCurrent,
                                bool *Melted, int &nzActive, int &ZBound_Low, int &ZBound_High, ViewI_H LayerID) {

    // Contrained solidification test problem
    // Since Z = 0 cells will be initialized as wall cells, it isn't necessary to iterate over them
    // Z = 1 through Z = nz - 1 are the cells of interest for solidification
    ZBound_Low = 1;
    ZBound_High = nz - 1;
    nzActive = ZBound_High - ZBound_Low + 1;
    if (id == 0)
        std::cout << "Constrained solidification problem's active domain is from Z = " << ZBound_Low << " through "
                  << ZBound_High << std::endl;

    // Initialize temperature field in Z direction with thermal gradient G set in input file
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < MyXSlices; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                int GlobalD3D1ConvPosition = k * MyXSlices * MyYSlices + i * MyYSlices + j;
                UndercoolingCurrent(GlobalD3D1ConvPosition) = 0;
                LayerID(GlobalD3D1ConvPosition) = 0;
                UndercoolingChange(GlobalD3D1ConvPosition) = R * deltat;
                CritTimeStep(GlobalD3D1ConvPosition) = (int)(((k - 1) * G * deltax) / (R * deltat));
                int GlobalX = i + MyXOffset;
                int GlobalY = j + MyYOffset;
                if ((GlobalX > -1) && (GlobalX < nx) && (GlobalY > -1) && (GlobalY < ny) && (k > 0)) {
                    Melted[GlobalD3D1ConvPosition] = true;
                }
                else {
                    Melted[GlobalD3D1ConvPosition] = false;
                }
            }
        }
    }
}

// Initialize temperature data for an array of overlapping spot melts
void TempInit_SpotMelt(double G, double R, std::string, int id, int &MyXSlices, int &MyYSlices, int &MyXOffset,
                       int &MyYOffset, double deltax, double deltat, int nz, ViewI_H CritTimeStep,
                       ViewF_H UndercoolingChange, ViewF_H UndercoolingCurrent, bool *Melted, int LayerHeight,
                       int NumberOfLayers, int &nzActive, int &ZBound_Low, int &ZBound_High, double FreezingRange,
                       ViewI_H LayerID, int NSpotsX, int NSpotsY, int SpotRadius, int SpotOffset) {

    // First layer: bottom is at Z = 1 (Z = 0 is a wall)
    // Z = 2 through Z = SpotRadius + 1 contains temperature data
    ZBound_Low = 1;
    ZBound_High = SpotRadius + 1;
    nzActive = ZBound_High - ZBound_Low + 1;
    if (id == 0)
        std::cout << "Layer 0's active domain is from Z = " << ZBound_Low << " through " << ZBound_High << std::endl;
    // No cells intitially have undercooling, nor have melting/solidification data
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < MyXSlices; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                int GlobalD3D1ConvPosition = k * MyXSlices * MyYSlices + i * MyYSlices + j;
                UndercoolingCurrent(GlobalD3D1ConvPosition) = 0;
                Melted[GlobalD3D1ConvPosition] = false;
                LayerID(GlobalD3D1ConvPosition) = -1;
            }
        }
    }

    int NumberOfSpots = NSpotsX * NSpotsY;
    // Outer edges of spots are initialized at the liquidus temperature
    // Spots cool at constant rate R, spot thermal gradient = G
    // Time between "start" of next spot is the time it takes for the previous spot
    // to have entirely gone below the solidus temperature
    float IsothermVelocity = (R / G) * deltat / deltax;                                  // in cells per time step
    int TimeBetweenSpots = SpotRadius / IsothermVelocity + (FreezingRange / R) / deltat; // in time steps
    if (id == 0)
        std::cout << "Initializing temperature field for " << NumberOfSpots
                  << " spots, each of which takes approximately " << TimeBetweenSpots << " time steps to solidify"
                  << std::endl;
    for (int layernumber = 0; layernumber < NumberOfLayers; layernumber++) {
        for (int n = 0; n < NumberOfSpots; n++) {
            if (id == 0)
                std::cout << "Initializing spot " << n << " on layer " << layernumber << std::endl;
            // Initialize critical time step/cooling rate values for this spot/this layer
            int XSpotPos = 2 + SpotRadius + (n % NSpotsX) * SpotOffset;
            int YSpotPos = 2 + SpotRadius + (n / NSpotsX) * SpotOffset;
            int ZSpotPos = SpotRadius + 1 + LayerHeight * layernumber;
            for (int k = 2; k <= SpotRadius + 1; k++) {
                // Distance of this cell from the spot center
                float DistZ = (float)(ZSpotPos - (k + LayerHeight * layernumber));
                for (int i = 0; i < MyXSlices; i++) {
                    int XGlobal = i + MyXOffset;
                    float DistX = (float)(XSpotPos - XGlobal);
                    for (int j = 0; j < MyYSlices; j++) {
                        int YGlobal = j + MyYOffset;
                        float DistY = (float)(YSpotPos - YGlobal);
                        float TotDist = sqrt(DistX * DistX + DistY * DistY + DistZ * DistZ);
                        if (TotDist <= SpotRadius) {
                            int GlobalD3D1ConvPosition =
                                (k + layernumber * LayerHeight) * MyXSlices * MyYSlices + i * MyYSlices + j;
                            CritTimeStep(GlobalD3D1ConvPosition) =
                                1 + (int)(((float)(SpotRadius)-TotDist) / IsothermVelocity) + TimeBetweenSpots * n;
                            UndercoolingChange(GlobalD3D1ConvPosition) = R * deltat;
                            LayerID(GlobalD3D1ConvPosition) = layernumber;
                            Melted[GlobalD3D1ConvPosition] = true;
                        }
                    }
                }
            }
        }
    }
}

// Initialize temperature data for a problem using the reduced/sparse data format and input temperature data from
// file(s)
void TempInit_Reduced(int id, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset, double &deltax,
                      double HT_deltax, double deltat, int &nx, int &ny, int &nz, ViewI_H CritTimeStep,
                      ViewF_H UndercoolingChange, ViewF_H UndercoolingCurrent, float XMin, float YMin, float ZMin,
                      bool *Melted, float *ZMinLayer, float *ZMaxLayer, int LayerHeight, int NumberOfLayers,
                      int &nzActive, int &ZBound_Low, int &ZBound_High, int *FinishTimeStep, double FreezingRange,
                      ViewI_H LayerID, int *FirstValue, int *LastValue, std::vector<double> RawData) {

    // Initialize temperature views to 0
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < MyXSlices; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                int Coord3D1D = k * MyXSlices * MyYSlices + i * MyYSlices + j;
                CritTimeStep(Coord3D1D) = 0;
                UndercoolingChange(Coord3D1D) = 0.0;
                UndercoolingCurrent(Coord3D1D) = 0.0;
            }
        }
    }

    // Temperature data read
    double HTtoCAratio_unrounded = HT_deltax / deltax;
    double HTtoCAratio_floor = floor(HTtoCAratio_unrounded);
    if (((HTtoCAratio_unrounded - HTtoCAratio_floor) > 0.0005) && (id == 0)) {
        std::string error = "Error: Temperature data point spacing not evenly divisible by CA cell size";
        throw std::runtime_error(error);
    }
    else if (((HTtoCAratio_unrounded - HTtoCAratio_floor) > 0.000001) && (id == 0)) {
        std::cout << "Note: Adjusting cell size from " << deltax << " to " << HT_deltax / HTtoCAratio_floor
                  << " to "
                     "ensure even divisibility of CA cell size into temperature data spacing"
                  << std::endl;
    }
    // Adjust deltax to exact value based on temperature data spacing and ratio between heat transport/CA cell sizes
    deltax = HT_deltax / HTtoCAratio_floor;
    int HTtoCAratio = round(HT_deltax / deltax); // OpenFOAM/CA cell size ratio

    // With row/col 0 being wall cells and row/col 1 being solid cells outside of the melted area, the domain starts
    // at row/col 2 As the wall cells are not part of the physical domain (solid cells at row/col 1 are defined as X
    // = Y = 0, the melted region domain data starts at X = Y = deltax, with data points at X or Y = deltax +
    // N*HT_deltax through X or Y = nx-3 or ny-3

    // The X and Y bounds are the region (for this MPI rank) of the physical domain that needs to be read
    // Extends past the actual spatial extent of the local domain for purposes of interpolating from HT_deltax to
    // deltax
    int LowerXBound, LowerYBound, UpperXBound, UpperYBound;
    if (MyXOffset <= 2)
        LowerXBound = 2;
    else
        LowerXBound = MyXOffset - ((MyXOffset - 2) % HTtoCAratio);
    if (MyYOffset <= 2)
        LowerYBound = 2;
    else
        LowerYBound = MyYOffset - ((MyYOffset - 2) % HTtoCAratio);

    if (MyXOffset + MyXSlices - 1 >= nx - 3)
        UpperXBound = nx - 3;
    else
        UpperXBound = MyXOffset + MyXSlices - 1 + HTtoCAratio - ((MyXOffset + (MyXSlices - 1) - 2) % HTtoCAratio);
    if (MyYOffset + MyYSlices - 1 >= ny - 3)
        UpperYBound = ny - 3;
    else
        UpperYBound = MyYOffset + MyYSlices - 1 + HTtoCAratio - ((MyYOffset + (MyYSlices - 1) - 2) % HTtoCAratio);

    // No sites have melted yet
    for (int i = 0; i < MyXSlices * MyYSlices * nz; i++) {
        Melted[i] = false;
        LayerID(i) = -1;
    }

    // removed unused LayerwiseTSOffset variable

    for (int LayerCounter = 0; LayerCounter < NumberOfLayers; LayerCounter++) {

        double SmallestTime = 1000000000;
        double SmallestTime_Global = 1000000000;
        double LargestTime = 0;
        double LargestTime_Global = 0;

        // How many CA cells in the vertical direction are needed to hold this layer's temperature data?
        int nzTempValuesThisLayer =
            round((ZMaxLayer[LayerCounter] - ZMinLayer[LayerCounter]) / deltax) +
            1; // (note this doesn't include the 2 rows of wall/active cells at the bottom surface)
        if (id == 0)
            std::cout << "Initializing temporary temperature data structures with " << nzTempValuesThisLayer
                      << " cells in z direction" << std::endl;
        if (id == 0)
            std::cout << "Layer " << LayerCounter << " rank " << id << " ZMin this layer is " << ZMinLayer[LayerCounter]
                      << std::endl;
        std::vector<std::vector<std::vector<double>>> CR, CritTL;
        for (int k = 0; k < nzTempValuesThisLayer; k++) {
            std::vector<std::vector<double>> TemperatureXX;
            for (int i = LowerXBound; i <= UpperXBound; i++) {
                std::vector<double> TemperatureX;
                for (int j = LowerYBound; j <= UpperYBound; j++) {
                    TemperatureX.push_back(-1.0);
                }
                TemperatureXX.push_back(TemperatureX);
            }
            CR.push_back(TemperatureXX);
            CritTL.push_back(TemperatureXX);
        }
        // Data was already read into the "RawData" temporary data structure
        // Determine which section of "RawData" is relevant for this layer of the overall domain
        int StartRange = FirstValue[LayerCounter];
        int EndRange = LastValue[LayerCounter];
        int XInt = -1;
        int YInt = -1;
        int ZInt = -1;
        if (id == 0)
            std::cout << "Range for layer " << LayerCounter << " on rank 0 is " << StartRange << " to " << EndRange
                      << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        for (int i = StartRange; i < EndRange; i++) {
            // Pos = 3 contains melting time data - not currently used as CA does not yet include remelting
            int Pos = i % 6;
            if (Pos == 0) {
                XInt = round((RawData[i] - XMin) / deltax) + 2;
            }
            else if (Pos == 1) {
                YInt = round((RawData[i] - YMin) / deltax) + 2;
            }
            else if (Pos == 2) {
                ZInt = round((RawData[i] + deltax * LayerHeight * LayerCounter - ZMinLayer[LayerCounter]) / deltax);
            }
            else if (Pos == 4) {
                // Liquidus time - only keep the last time that this point went below the liquidus
                if (RawData[i] > CritTL[ZInt][XInt - LowerXBound][YInt - LowerYBound]) {
                    CritTL[ZInt][XInt - LowerXBound][YInt - LowerYBound] = RawData[i];
                    if (RawData[i] < SmallestTime)
                        SmallestTime = RawData[i];
                }
                else {
                    // This is not the last time that the cell is going below the liquidus,
                    // do not store cooling rate data - skip to next point on the list
                    i = i + 1;
                }
            }
            else if (Pos == 5) {
                CR[ZInt][XInt - LowerXBound][YInt - LowerYBound] = RawData[i];
                float SolidusTime = CritTL[ZInt][XInt - LowerXBound][YInt - LowerYBound] +
                                    FreezingRange / CR[ZInt][XInt - LowerXBound][YInt - LowerYBound];
                if (SolidusTime > LargestTime)
                    LargestTime = SolidusTime;
            }
        }
        // If reading data from files without a script, time values start at 0 for each layer
        // If reading data with input from a script time values each layer are continuous, are should be
        // renormalized to 0 for each layer
        MPI_Reduce(&LargestTime, &LargestTime_Global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Bcast(&LargestTime_Global, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Reduce(&SmallestTime, &SmallestTime_Global, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Bcast(&SmallestTime_Global, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (id == 0) {
            std::cout << "Smallest time globally for layer " << LayerCounter << " is " << SmallestTime_Global
                      << std::endl;
            std::cout << "Largest time globally for layer " << LayerCounter << " is " << LargestTime_Global
                      << std::endl;
        }
        // renormalize last time step value with the start of the layer as time step 0
        FinishTimeStep[LayerCounter] = round((LargestTime_Global - SmallestTime_Global) / deltat);
        if (id == 0)
            std::cout << " Layer " << LayerCounter << " FINISH TIME STEP IS " << FinishTimeStep[LayerCounter]
                      << std::endl;
        if (id == 0)
            std::cout << "Layer " << LayerCounter << " temperatures read" << std::endl;

        // Data interpolation between heat transport and CA grids, if necessary
        if (HTtoCAratio != 1) {
            for (int k = 0; k < nzTempValuesThisLayer; k++) {
                int LowZ = k - (k % HTtoCAratio);
                int HighZ = LowZ + HTtoCAratio;
                double FHighZ = (double)(k - LowZ) / (double)(HTtoCAratio);
                double FLowZ = 1.0 - FHighZ;
                if (HighZ > nzTempValuesThisLayer - 1)
                    HighZ = LowZ;
                for (int i = 0; i <= UpperXBound - LowerXBound; i++) {
                    int LowX = i - (i % HTtoCAratio);
                    int HighX = LowX + HTtoCAratio;
                    double FHighX = (double)(i - LowX) / (double)(HTtoCAratio);
                    double FLowX = 1.0 - FHighX;
                    if (HighX >= UpperXBound - LowerXBound)
                        HighX = UpperXBound - LowerXBound;

                    for (int j = 0; j <= UpperYBound - LowerYBound; j++) {
                        int LowY = j - (j % HTtoCAratio);
                        int HighY = LowY + HTtoCAratio;
                        double FHighY = (float)(j - LowY) / (float)(HTtoCAratio);
                        double FLowY = 1.0 - FHighY;
                        if (HighY >= UpperYBound - LowerYBound)
                            HighY = UpperYBound - LowerYBound;
                        double Pt1 = CritTL[LowZ][LowX][LowY];
                        double Pt2 = CritTL[LowZ][HighX][LowY];
                        double Pt12 = FLowX * Pt1 + FHighX * Pt2;
                        double Pt3 = CritTL[LowZ][LowX][HighY];
                        double Pt4 = CritTL[LowZ][HighX][HighY];
                        double Pt34 = FLowX * Pt3 + FHighX * Pt4;
                        double Pt1234 = Pt12 * FLowY + Pt34 * FHighY;
                        double Pt5 = CritTL[HighZ][LowX][LowY];
                        double Pt6 = CritTL[HighZ][HighX][LowY];
                        double Pt56 = FLowX * Pt5 + FHighX * Pt6;
                        double Pt7 = CritTL[HighZ][LowX][HighY];
                        double Pt8 = CritTL[HighZ][HighX][HighY];
                        double Pt78 = FLowX * Pt7 + FHighX * Pt8;
                        double Pt5678 = Pt56 * FLowY + Pt78 * FHighY;
                        if ((Pt1 > 0) && (Pt2 > 0) && (Pt3 > 0) && (Pt4 > 0) && (Pt5 > 0) && (Pt6 > 0) && (Pt7 > 0) &&
                            (Pt8 > 0)) {
                            CritTL[k][i][j] = Pt1234 * FLowZ + Pt5678 * FHighZ;
                        }
                        Pt1 = CR[LowZ][LowX][LowY];
                        Pt2 = CR[LowZ][HighX][LowY];
                        Pt12 = FLowX * Pt1 + FHighX * Pt2;
                        Pt3 = CR[LowZ][LowX][HighY];
                        Pt4 = CR[LowZ][HighX][HighY];
                        Pt34 = FLowX * Pt3 + FHighX * Pt4;
                        Pt1234 = Pt12 * FLowY + Pt34 * FHighY;
                        Pt5 = CR[HighZ][LowX][LowY];
                        Pt6 = CR[HighZ][HighX][LowY];
                        Pt56 = FLowX * Pt5 + FHighX * Pt6;
                        Pt7 = CR[HighZ][LowX][HighY];
                        Pt8 = CR[HighZ][HighX][HighY];
                        Pt78 = FLowX * Pt7 + FHighX * Pt8;
                        Pt5678 = Pt56 * FLowY + Pt78 * FHighY;
                        if ((Pt1 > 0) && (Pt2 > 0) && (Pt3 > 0) && (Pt4 > 0) && (Pt5 > 0) && (Pt6 > 0) && (Pt7 > 0) &&
                            (Pt8 > 0)) {
                            CR[k][i][j] = Pt1234 * FLowZ + Pt5678 * FHighZ;
                        }
                    }
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Interpolation done" << std::endl;

        // Convert CritTL, CritTS matrices into CritTimeStep and UndercoolingChange (change in undercooling with
        // time step) "ZMin" is the global Z coordinate that corresponds to cells at Z = 2 (Z = 0 is the domain's
        // bottom wall, Z = 1 are the active cells just outside of the melt pool) "ZMax" is the global Z coordinate
        // that corresponds to cells at Z = nz-2 (Z = nz-1 is the domain's top wall)
        if (LayerCounter == 0) {
            ZBound_Low = 1;
            ZBound_High = nzTempValuesThisLayer + 1;
            nzActive = ZBound_High - ZBound_Low + 1;
            if (id == 0)
                std::cout << "Active domain for layer 0 is Z = 1 through " << ZBound_High << std::endl;
        }
        if (id == 0)
            std::cout << "Layer " << LayerCounter << " data belongs to global z coordinates of "
                      << round((ZMinLayer[LayerCounter] - ZMin) / deltax) + 2 << " through "
                      << round((ZMinLayer[LayerCounter] - ZMin) / deltax) + 2 + nzTempValuesThisLayer - 1 << std::endl;

        for (int k = 0; k < nzTempValuesThisLayer; k++) {
            for (int ii = LowerXBound; ii <= UpperXBound; ii++) {
                for (int jj = LowerYBound; jj <= UpperYBound; jj++) {
                    if ((ii >= MyXOffset) && (ii < MyXOffset + MyXSlices) && (jj >= MyYOffset) &&
                        (jj < MyYOffset + MyYSlices)) {
                        int Adj_i = ii - MyXOffset;
                        int Adj_j = jj - MyYOffset;
                        // Liquidus time normalized to the time at which the layer started solidifying
                        double CTLiq = CritTL[k][ii - LowerXBound][jj - LowerYBound] - SmallestTime_Global;
                        if (CTLiq > 0) {
                            // Where does this layer's temperature data belong on the global (including all layers)
                            // grid? Adjust Z coordinate by ZMin
                            int ZOffset = round((ZMinLayer[LayerCounter] - ZMin) / deltax) + k + 2;
                            int Coord3D1D = ZOffset * MyXSlices * MyYSlices + Adj_i * MyYSlices + Adj_j;
                            Melted[Coord3D1D] = true;
                            CritTimeStep(Coord3D1D) = round(CTLiq / deltat);
                            LayerID(Coord3D1D) = LayerCounter;
                            UndercoolingChange(Coord3D1D) =
                                std::abs(CR[k][ii - LowerXBound][jj - LowerYBound]) * deltat;
                        }
                    }
                }
            }
        }
    } // End read over all temperature files and placement of data

    if (id == 0)
        std::cout << "First layer Z bounds are " << ZBound_Low << " and " << ZBound_High << std::endl;
}
//*****************************************************************************/
// Initialize grain orientations and unit vectors
void OrientationInit(int id, int NGrainOrientations, ViewI_H GrainOrientation, ViewF_H GrainUnitVector,
                     std::string GrainOrientationFile) {

    // Read file of grain orientations
    std::ifstream O;
    O.open(GrainOrientationFile);

    // Line 1 is the number of orientation values to read (if not specified already)
    std::string ValueRead;
    getline(O, ValueRead);

    // Populate data structure for grain unit vectors
    for (int i = 0; i < NGrainOrientations; i++) {
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
            GrainUnitVector(9 * i + 3 * UVNumber + Comp) = ReadGO;

            Comp++;
            if (Comp > 2) {
                Comp = 0;
                UVNumber++;
            }
        }
    }
    O.close();

    // The grain orientations that correspond to each Grain ID must be the same across all ranks
    // Shuffle list of "NGrainOrientation" orientations
    int *GrainOrientation_master = new int[NGrainOrientations];
    if (id == 0) {
        for (int h = 0; h < NGrainOrientations; h++) {
            GrainOrientation_master[h] = h;
        }
    }
    MPI_Bcast(&GrainOrientation_master[0], NGrainOrientations, MPI_INT, 0, MPI_COMM_WORLD);
    for (int h = 0; h < NGrainOrientations; h++) {
        GrainOrientation(h) = GrainOrientation_master[h];
    }
}

// Initializes cell types and epitaxial Grain ID values where substrate grains are active cells on the bottom surface of
// the constrained domain
void SubstrateInit_ConstrainedGrowth(double FractSurfaceSitesActive, int MyXSlices, int MyYSlices, int nx, int ny,
                                     int nz, int MyXOffset, int MyYOffset, int id, int np, ViewI_H CellType,
                                     ViewI_H GrainID) {

    std::mt19937_64 gen(id);
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // Counter for the number of active cells
    int SubstrateActCells_ThisRank = 0;

    // Wall cells at global domain boundaries (other than the top Z boundary)
    // Other cells are either regions that will melt, or part of the substrate
    for (int k = 1; k < nz; k++) {
        for (int i = 0; i < MyXSlices; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                int GlobalX = i + MyXOffset;
                int GlobalY = j + MyYOffset;
                int CAGridLocation = k * MyXSlices * MyYSlices + i * MyYSlices + j;
                if ((GlobalX == -1) || (GlobalX == nx) || (GlobalY == -1) || (GlobalY == ny) || (k == 0)) {
                    CellType(CAGridLocation) = Wall;
                    GrainID(CAGridLocation) = 0;
                }
                else {
                    if ((k == 1) && (i > 0) && (i < MyXSlices - 1) && (j > 0) && (j < MyYSlices - 1)) {
                        // Randomly locate substrate grain seeds - do not place seeds in the ghost nodes, as this would
                        // double count these cells as possible grain sites
                        double R = dis(gen);
                        if (R < FractSurfaceSitesActive) {
                            SubstrateActCells_ThisRank++;
                            CellType(CAGridLocation) = Active;
                        }
                        else {
                            CellType(CAGridLocation) = Liquid;
                        }
                    }
                    else
                        CellType(CAGridLocation) = Liquid;
                }
            }
        }
    }
    // Assign grain IDs to bottom surface grains
    int FirstEpitaxialGrainID = 1;
    if (np > 1) {
        // Grains for epitaxial growth - determine GrainIDs on each MPI rank
        if (id == 0) {
            int SBuf = FirstEpitaxialGrainID + SubstrateActCells_ThisRank;
            MPI_Send(&SBuf, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        }
        else if (id == np - 1) {
            int RBuf;
            MPI_Recv(&RBuf, 1, MPI_INT, np - 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            FirstEpitaxialGrainID = RBuf;
        }
        else {
            int RBuf;
            MPI_Recv(&RBuf, 1, MPI_INT, id - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            FirstEpitaxialGrainID = RBuf;
            int SBuf = RBuf + SubstrateActCells_ThisRank;
            MPI_Send(&SBuf, 1, MPI_INT, id + 1, 0, MPI_COMM_WORLD);
        }
    }
    for (int i = 0; i < MyXSlices * MyYSlices * nz; i++) {
        if (CellType(i) == Active) {
            GrainID(i) = FirstEpitaxialGrainID;
            FirstEpitaxialGrainID++;
        }
    }
}

// Initializes Grain ID values where the substrate comes from a file
void SubstrateInit_FromFile(std::string SubstrateFileName, int nz, int MyXSlices, int MyYSlices, int MyXOffset,
                            int MyYOffset, int id, ViewI_H CritTimeStep, ViewI_H GrainID) {

    // Assign GrainID values to cells that are part of the substrate
    std::ifstream Substrate;
    Substrate.open(SubstrateFileName);
    int Substrate_LowX = MyXOffset;
    int Substrate_HighX = MyXOffset + MyXSlices;
    int Substrate_LowY = MyYOffset;
    int Substrate_HighY = MyYOffset + MyYSlices;
    int nxS, nyS, nzS;
    std::string s;
    getline(Substrate, s);
    std::size_t found = s.find("=");
    std::string str = s.substr(found + 1, s.length() - 1);
    nzS = stoi(str, nullptr, 10);
    getline(Substrate, s);
    found = s.find("=");
    str = s.substr(found + 1, s.length() - 1);
    nyS = stoi(str, nullptr, 10);
    getline(Substrate, s);
    found = s.find("=");
    str = s.substr(found + 1, s.length() - 1);
    nxS = stoi(str, nullptr, 10);
    if ((id == 0) && (nzS < nz)) {
        // Do not allow simulation if there is inssufficient substrate data in the specified file
        std::string error = "Error: only " + std::to_string(nzS) +
                            " layers of substrate data are present in the file " + SubstrateFileName + " ; at least " +
                            std::to_string(nz) +
                            " layers of substrate data are required to simulate the specified solidification problem";
        throw std::runtime_error(error);
    }

    // Assign GrainID values to cells that are part of the substrate
    // Cells that border the melted region are type active, others are type solid
    for (int k = 1; k < nzS; k++) {
        if (k == nz)
            break;
        for (int j = 0; j < nyS; j++) {
            for (int i = 0; i < nxS; i++) {
                std::string GIDVal;
                getline(Substrate, GIDVal);
                if ((i >= Substrate_LowX) && (i < Substrate_HighX) && (j >= Substrate_LowY) && (j < Substrate_HighY)) {
                    int CAGridLocation;
                    CAGridLocation = k * MyXSlices * MyYSlices + (i - MyXOffset) * MyYSlices + (j - MyYOffset);
                    if (CritTimeStep(CAGridLocation) == 0) {
                        // This cell is part of the substrate
                        GrainID(CAGridLocation) = stoi(GIDVal, nullptr, 10);
                    }
                    else {
                        // This cell is part of at least one layer's melt pool footprint
                        GrainID(CAGridLocation) = 0;
                    }
                }
            }
        }
    }
    Substrate.close();
    if (id == 0)
        std::cout << "Substrate file read complete" << std::endl;
}

// Initializes Grain ID values where the baseplate is generated using an input grain spacing and a Voronoi Tessellation,
// while the remainder of the interface is seeded with CA-cell sized substrate grains (emulating bulk nucleation
// alongside the edges of partially melted powder particles)
void SubstrateInit_FromGrainSpacing(float SubstrateGrainSpacing, int nx, int ny, int nz, int nzActive, int MyXSlices,
                                    int MyYSlices, int MyXOffset, int MyYOffset, int LocalActiveDomainSize, int id,
                                    int np, double deltax, ViewI_H GrainID, ViewI_H CritTimeStep) {

    // Seed random number generator such that each rank generates the same baseplate grain center locations
    double SubstrateSeed = 1.0;
    std::mt19937_64 gen(SubstrateSeed);
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // Probability that a given cell will be the center of a baseplate grain
    double BaseplateGrainProb = (deltax * deltax * deltax) /
                                (SubstrateGrainSpacing * SubstrateGrainSpacing * SubstrateGrainSpacing * pow(10, -18));
    ViewI_H NumBaseplateGrains_H(Kokkos::ViewAllocateWithoutInitializing("NBaseplate"), 1);
    NumBaseplateGrains_H(0) = 0;
    ViewI_H BaseplateGrainX_H(Kokkos::ViewAllocateWithoutInitializing("BaseplateGrainX"), nx * ny * nzActive);
    ViewI_H BaseplateGrainY_H(Kokkos::ViewAllocateWithoutInitializing("BaseplateGrainY"), nx * ny * nzActive);
    ViewI_H BaseplateGrainZ_H(Kokkos::ViewAllocateWithoutInitializing("BaseplateGrainZ"), nx * ny * nzActive);

    // For the entire baseplate (all x and y coordinate, but only layer 0 z coordinates), identify baseplate grain
    // centers This will eventually be done on the device when random number generation with Kokkos is more efficient
    if (id == 0)
        std::cout << "Baseplate spanning domain coordinates Z = 1 through " << nzActive << std::endl;
    for (int k = 1; k <= nzActive; k++) {
        for (int i = 1; i < nx - 1; i++) {
            for (int j = 1; j < ny - 1; j++) {
                double R = dis(gen);
                if (R < BaseplateGrainProb) {
                    int OldIndexValue = NumBaseplateGrains_H(0);
                    BaseplateGrainX_H(OldIndexValue) = i;
                    BaseplateGrainY_H(OldIndexValue) = j;
                    BaseplateGrainZ_H(OldIndexValue) = k;
                    NumBaseplateGrains_H(0)++;
                }
            }
        }
    }
    if (id == 0)
        std::cout << "Number of baseplate grains: " << NumBaseplateGrains_H(0) << std::endl;
    using memory_space = Kokkos::DefaultExecutionSpace::memory_space;
    ViewI NumBaseplateGrains_G = Kokkos::create_mirror_view_and_copy(memory_space(), NumBaseplateGrains_H);
    ViewI BaseplateGrainX_G = Kokkos::create_mirror_view_and_copy(memory_space(), BaseplateGrainX_H);
    ViewI BaseplateGrainY_G = Kokkos::create_mirror_view_and_copy(memory_space(), BaseplateGrainY_H);
    ViewI BaseplateGrainZ_G = Kokkos::create_mirror_view_and_copy(memory_space(), BaseplateGrainZ_H);

    // While eventually all of the initialization routines will be performed on the device (minus reading from files),
    // for now we initialize GrainID on the device, copy back to the host for the remainder of the initialization
    // routines, and then copy back to the device Also need a copy of CritTimeStep on the device to do this This is done
    // because the loop was too slow on the host/in serial
    using memory_space = Kokkos::DefaultExecutionSpace::memory_space;
    ViewI CritTimeStep_G = Kokkos::create_mirror_view_and_copy(memory_space(), CritTimeStep);
    ViewI GrainID_G = Kokkos::create_mirror_view_and_copy(memory_space(), GrainID);
    Kokkos::parallel_for(
        "CellCapture",
        Kokkos::MDRangePolicy<Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>(
            {1, 0, 0}, {nzActive + 1, MyXSlices, MyYSlices}),
        KOKKOS_LAMBDA(const int k, const int i, const int j) {
            int GlobalX = i + MyXOffset;
            int GlobalY = j + MyYOffset;
            int CAGridLocation = k * MyXSlices * MyYSlices + i * MyYSlices + j;
            if (CritTimeStep_G(CAGridLocation) == 0) {
                // This cell is part of the substrate - determine which grain center the cell is closest to, in order to
                // assign it a grain ID If closest to grain "n", assign grain ID "n+1" (grain ID = 0 is not used)
                float MinDistanceToThisGrain = (float)(LocalActiveDomainSize);
                int ClosestGrainIndex = -1;
                for (int n = 0; n < NumBaseplateGrains_G(); n++) {
                    float DistanceToThisGrainX = (float)(abs(BaseplateGrainX_G(n) - GlobalX));
                    float DistanceToThisGrainY = (float)(abs(BaseplateGrainY_G(n) - GlobalY));
                    float DistanceToThisGrainZ = (float)(abs(BaseplateGrainZ_G(n) - k));
                    float DistanceToThisGrain = sqrtf(DistanceToThisGrainX * DistanceToThisGrainX +
                                                      DistanceToThisGrainY * DistanceToThisGrainY +
                                                      DistanceToThisGrainZ * DistanceToThisGrainZ);
                    if (DistanceToThisGrain < MinDistanceToThisGrain) {
                        ClosestGrainIndex = n;
                        MinDistanceToThisGrain = DistanceToThisGrain;
                    }
                }
                GrainID_G(CAGridLocation) = ClosestGrainIndex + 1;
            }
            else {
                // This cell is part of layer 0's melt pool footprint
                GrainID_G(CAGridLocation) = 0;
            }
        });

    // Copy Grain ID back to host
    Kokkos::deep_copy(GrainID, GrainID_G);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Baseplate grain structure initialized" << std::endl;
    // Initialize grain seeds above baseplate to emulate bulk nucleation at edge of melted powder particles
    int PowderLayerActCells_ThisRank = 0;
    for (int k = nzActive + 1; k < nz; k++) {
        for (int i = 1; i < MyXSlices - 1; i++) {
            for (int j = 1; j < MyYSlices - 1; j++) {
                int CAGridLocation = k * MyXSlices * MyYSlices + i * MyYSlices + j;
                if (CritTimeStep(CAGridLocation) == 0) {
                    // This cell is part of the powder layer - count how many of these exist on this rank
                    PowderLayerActCells_ThisRank++;
                }
                else {
                    // This cell is part of one of the layer's melt pool footprint
                    GrainID(CAGridLocation) = 0;
                }
            }
        }
    }
    // Assign grain IDs to bulk grain nuclei
    int FirstEpitaxialGrainID = NumBaseplateGrains_H(0) + 1;
    if (np > 1) {
        // Grains for epitaxial growth - determine GrainIDs on each MPI rank
        if (id == 0) {
            int SBuf = FirstEpitaxialGrainID + PowderLayerActCells_ThisRank;
            MPI_Send(&SBuf, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        }
        else if (id == np - 1) {
            int RBuf;
            MPI_Recv(&RBuf, 1, MPI_INT, np - 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            FirstEpitaxialGrainID = RBuf;
        }
        else {
            int RBuf;
            MPI_Recv(&RBuf, 1, MPI_INT, id - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            FirstEpitaxialGrainID = RBuf;
            int SBuf = RBuf + PowderLayerActCells_ThisRank;
            MPI_Send(&SBuf, 1, MPI_INT, id + 1, 0, MPI_COMM_WORLD);
        }
    }
    for (int D3D1ConvPosition = MyXSlices * MyYSlices * (nzActive + 1); D3D1ConvPosition < MyXSlices * MyYSlices * nz;
         D3D1ConvPosition++) {
        if (CritTimeStep(D3D1ConvPosition) == 0) {
            GrainID(D3D1ConvPosition) = FirstEpitaxialGrainID;
            FirstEpitaxialGrainID++;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Initialized baseplate and powder grain structure" << std::endl;
}

//*****************************************************************************/
// Initializes cells at border of solid and liquid as active type, wall cells along global domain bounds
void ActiveCellWallInit(int id, int MyXSlices, int MyYSlices, int nx, int ny, int nz, int MyXOffset, int MyYOffset,
                        ViewI_H CellType, ViewI_H GrainID, ViewI_H CritTimeStep, ViewI2D_H ItList, ViewI_H NeighborX,
                        ViewI_H NeighborY, ViewI_H NeighborZ, ViewF_H UndercoolingChange, bool ExtraWalls) {

    // Wall cells at global domain boundaries
    // Other cells are either regions that will melt, or part of the substrate
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < MyXSlices; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                int GlobalX = i + MyXOffset;
                int GlobalY = j + MyYOffset;
                int CAGridLocation = k * MyXSlices * MyYSlices + i * MyYSlices + j;
                // Walls at X/Y boundaries, bottom Z boundary but not top Z boundary
                if ((GlobalX == -1) || (GlobalX == nx) || (GlobalY == -1) || (GlobalY == ny) || (k == 0)) {
                    CellType(CAGridLocation) = Wall;
                    GrainID(CAGridLocation) = 0;
                }
                if (ExtraWalls) {
                    if ((GlobalX == 0) || (GlobalX == nx - 1) || (GlobalY == 0) || (GlobalY == ny - 1) ||
                        (GlobalX == 1) || (GlobalX == nx - 2) || (GlobalY == 1) || (GlobalY == ny - 2)) {
                        CellType(CAGridLocation) = Wall;
                        GrainID(CAGridLocation) = 0;
                    }
                }
            }
        }
    }
    // Initialize solid cells where no temperature data exists, and liquid cells where temperature data exists
    // This is done prior to initializing active cells, as active cells are initialized based on neighbor cell types
    for (int k = 1; k < nz; k++) {
        for (int j = 0; j < MyYSlices; j++) {
            for (int i = 0; i < MyXSlices; i++) {
                int CAGridLocation = k * MyXSlices * MyYSlices + i * MyYSlices + j;
                if (CritTimeStep(CAGridLocation) == 0) {
                    CellType(CAGridLocation) = Solid;
                }
                else {
                    // This cell has associated melting data and is initialized as liquid
                    CellType(CAGridLocation) = Liquid;
                }
            }
        }
    }
    // Count number of active cells are at the solid-liquid boundary
    int SubstrateActCells_ThisRank = 0;
    for (int k = 1; k < nz; k++) {
        for (int j = 0; j < MyYSlices; j++) {
            for (int i = 0; i < MyXSlices; i++) {
                int CAGridLocation = k * MyXSlices * MyYSlices + i * MyYSlices + j;
                if (CritTimeStep(CAGridLocation) == 0) {
                    // This is a solid or active cell, depending on whether it is located at the interface of the liquid
                    // Check to see if this site is actually at the solid-liquid interface
                    int LCount = 0;
                    // Which neighbors should be iterated over?
                    int ItBounds = FindItBounds(i, j, MyXSlices, MyYSlices);
                    int NListLength;
                    if (ItBounds == 0) {
                        NListLength = 26;
                    }
                    else if (ItBounds > 4) {
                        NListLength = 11;
                    }
                    else {
                        NListLength = 17;
                    }
                    // "ll" corresponds to the specific position on the list of neighboring cells
                    for (int ll = 0; ll < NListLength; ll++) {
                        // "l" correpsponds to the specific neighboring cell
                        int l = ItList(ItBounds, ll);
                        // Local coordinates of adjacent cell center
                        int MyNeighborX = i + NeighborX(l);
                        int MyNeighborY = j + NeighborY(l);
                        int MyNeighborZ = k + NeighborZ(l);
                        int NeighborD3D1ConvPosition =
                            MyNeighborZ * MyXSlices * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
                        if (MyNeighborZ < nz) {
                            if ((CritTimeStep(NeighborD3D1ConvPosition) > 0) &&
                                (CellType(NeighborD3D1ConvPosition) != Wall)) {
                                LCount++;
                                // At the interface - assign this cell a value for undercooling change taken from one of
                                // its neighbors
                                CellType(CAGridLocation) = Active;
                                UndercoolingChange(CAGridLocation) = UndercoolingChange(NeighborD3D1ConvPosition);
                                SubstrateActCells_ThisRank++;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    int TotalSubstrateActCells;
    MPI_Reduce(&SubstrateActCells_ThisRank, &TotalSubstrateActCells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Number of substrate active cells across all ranks: " << TotalSubstrateActCells << std::endl;
}

//*****************************************************************************/
// Initializes cell types where the substrate comes from a file
void GrainInit(int layernumber, int NGrainOrientations, int DecompositionStrategy, int nz, int LocalActiveDomainSize,
               int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int id, int np, int NeighborRank_North,
               int NeighborRank_South, int NeighborRank_East, int NeighborRank_West, int NeighborRank_NorthEast,
               int NeighborRank_NorthWest, int NeighborRank_SouthEast, int NeighborRank_SouthWest, ViewI2D_H ItList,
               ViewI_H NeighborX, ViewI_H NeighborY, ViewI_H NeighborZ, ViewI_H GrainOrientation,
               ViewF_H GrainUnitVector, ViewF_H DiagonalLength, ViewI_H CellType, ViewI_H GrainID,
               ViewF_H CritDiagonalLength, ViewF_H DOCenter, ViewI_H CritTimeStep, double deltax, double NMax,
               int &NextLayer_FirstNucleatedGrainID, int &PossibleNuclei_ThisRank, int ZBound_High, int ZBound_Low) {

    // RNG for heterogenous nuclei locations in the liquid
    std::mt19937_64 gen(id);
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // Probability that a given liquid site will be a potential nucleus location
    double BulkProb = NMax * deltax * deltax * deltax;

    // Count the number of nucleation events that may potentially occur on this rank (not counting ghost nodes, to avoid
    // double counting cells are potential nucleation sites)
    PossibleNuclei_ThisRank = 0;
    for (int k = 1; k < nz; k++) {
        for (int j = 1; j < MyYSlices - 1; j++) {
            for (int i = 1; i < MyXSlices - 1; i++) {
                int CAGridLocation = k * MyXSlices * MyYSlices + i * MyYSlices + j;
                if (CritTimeStep(CAGridLocation) != 0) {
                    double R = dis(gen);
                    if (R < BulkProb) {
                        PossibleNuclei_ThisRank++;
                        CellType(CAGridLocation) = TemporaryInit;
                    }
                }
            }
        }
    }

    int TotalNucleatedGrains;
    MPI_Reduce(&PossibleNuclei_ThisRank, &TotalNucleatedGrains, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Number of potential nucleated grains: " << TotalNucleatedGrains << std::endl;

    int FirstNucleatedGID_Rank0;
    if (layernumber == -1) {
        FirstNucleatedGID_Rank0 = -1;
    }
    else {
        FirstNucleatedGID_Rank0 = NextLayer_FirstNucleatedGrainID;
    }
    int MyFirstNGrainID; // First GrainID for nuclei on this rank, for this layer
    if (np > 1) {
        // Assign GrainIDs for nucleated grains (negative values)
        // Grains for nucleated growth
        if (id == 0) {
            int SBuf = FirstNucleatedGID_Rank0 - PossibleNuclei_ThisRank;
            MPI_Send(&SBuf, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
            MyFirstNGrainID = FirstNucleatedGID_Rank0;
        }
        else if (id == np - 1) {
            int RBuf;
            MPI_Recv(&RBuf, 1, MPI_INT, np - 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MyFirstNGrainID = RBuf;
            NextLayer_FirstNucleatedGrainID = MyFirstNGrainID - PossibleNuclei_ThisRank;
        }
        else {
            int RBuf;
            MPI_Recv(&RBuf, 1, MPI_INT, id - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MyFirstNGrainID = RBuf;
            int SBuf = MyFirstNGrainID - PossibleNuclei_ThisRank;
            MPI_Send(&SBuf, 1, MPI_INT, id + 1, 0, MPI_COMM_WORLD);
        }
        MPI_Bcast(&NextLayer_FirstNucleatedGrainID, 1, MPI_INT, np - 1, MPI_COMM_WORLD);
    }
    else {
        // No communication among ranks
        NextLayer_FirstNucleatedGrainID = FirstNucleatedGID_Rank0 - PossibleNuclei_ThisRank;
        MyFirstNGrainID = FirstNucleatedGID_Rank0;
    }

    // Assign Grain IDs to nucleated grains
    // Set up active cell data structures for appropriate cells
    int SouthNucleiSendCount = 0;
    int NorthNucleiSendCount = 0;
    int EastNucleiSendCount = 0;
    int WestNucleiSendCount = 0;
    int NorthWestNucleiSendCount = 0;
    int NorthEastNucleiSendCount = 0;
    int SouthWestNucleiSendCount = 0;
    int SouthEastNucleiSendCount = 0;

    // Set up active cell octahedra for growth, mark cell data to be communicated across ranks in ghost nodes
    // Assign GrainIDs to nuclei sites
    int NCounter = MyFirstNGrainID;

    // Nonactive cells should start with diagonal lengths of 0
    for (int i = 0; i < LocalActiveDomainSize; i++) {
        DiagonalLength(i) = 0.0;
    }
    for (int GlobalZ = 1; GlobalZ < nz; GlobalZ++) {
        for (int RankX = 0; RankX < MyXSlices; RankX++) {
            for (int RankY = 0; RankY < MyYSlices; RankY++) {
                long int D3D1ConvPositionGlobal = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                if (CellType(D3D1ConvPositionGlobal) == Active) {

                    // If part of the active domain, calculate Critical diagonal lengths
                    if (GlobalZ <= ZBound_High) {
                        int GlobalX = RankX + MyXOffset;
                        int GlobalY = RankY + MyYOffset;
                        int RankZ = GlobalZ - ZBound_Low;
                        int D3D1ConvPosition = RankZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                        int MyGrainID = GrainID(D3D1ConvPositionGlobal);
                        DiagonalLength(D3D1ConvPosition) = 0.01;
                        DOCenter(3 * D3D1ConvPosition) = GlobalX + 0.5;
                        DOCenter(3 * D3D1ConvPosition + 1) = GlobalY + 0.5;
                        DOCenter(3 * D3D1ConvPosition + 2) = GlobalZ + 0.5;

                        // The orientation for the new grain will depend on its Grain ID
                        int MyOrientation = GrainOrientation(((abs(MyGrainID) - 1) % NGrainOrientations));
                        // Calculate critical values at which this active cell leads to the activation of a neighboring
                        // liquid cell (xp,yp,zp) is the new cell's center on the global grid
                        double xp = GlobalX + 0.5;
                        double yp = GlobalY + 0.5;
                        double zp = GlobalZ + 0.5;

                        float cx = DOCenter((long int)(3 * D3D1ConvPosition));
                        float cy = DOCenter((long int)(3 * D3D1ConvPosition + 1));
                        float cz = DOCenter((long int)(3 * D3D1ConvPosition + 2));

                        // Calculate critical diagonal lengths for the new active cell located at (xp,yp,zp) on the
                        // local grid For each neighbor (l=0 to 25), calculate which octahedron face leads to cell
                        // capture Calculate critical octahedron diagonal length to activate each nearest neighbor, as
                        // well as the coordinates of the triangle vertices on the capturing face
                        for (int n = 0; n < 26; n++) {

                            // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the image of
                            // the center of a neighbor cell
                            double x0 = xp + NeighborX(n) - cx;
                            double y0 = yp + NeighborY(n) - cy;
                            double z0 = zp + NeighborZ(n) - cz;
                            // mag0 is the magnitude of (x0,y0,z0)
                            double mag0 = pow(pow(x0, 2.0) + pow(y0, 2.0) + pow(z0, 2.0), 0.5);

                            // Calculate unit vectors for the octahedron that intersect the new cell center
                            double Diag1X, Diag1Y, Diag1Z, Diag2X, Diag2Y, Diag2Z, Diag3X, Diag3Y, Diag3Z;
                            double Angle1 =
                                (GrainUnitVector(9 * MyOrientation) * x0 + GrainUnitVector(9 * MyOrientation + 1) * y0 +
                                 GrainUnitVector(9 * MyOrientation + 2) * z0) /
                                mag0;
                            if (Angle1 < 0) {
                                Diag1X = GrainUnitVector(9 * MyOrientation);
                                Diag1Y = GrainUnitVector(9 * MyOrientation + 1);
                                Diag1Z = GrainUnitVector(9 * MyOrientation + 2);
                            }
                            else {
                                Diag1X = -GrainUnitVector(9 * MyOrientation);
                                Diag1Y = -GrainUnitVector(9 * MyOrientation + 1);
                                Diag1Z = -GrainUnitVector(9 * MyOrientation + 2);
                            }
                            double Angle2 = (GrainUnitVector(9 * MyOrientation + 3) * x0 +
                                             GrainUnitVector(9 * MyOrientation + 4) * y0 +
                                             GrainUnitVector(9 * MyOrientation + 5) * z0) /
                                            mag0;
                            if (Angle2 < 0) {
                                Diag2X = GrainUnitVector(9 * MyOrientation + 3);
                                Diag2Y = GrainUnitVector(9 * MyOrientation + 4);
                                Diag2Z = GrainUnitVector(9 * MyOrientation + 5);
                            }
                            else {
                                Diag2X = -GrainUnitVector(9 * MyOrientation + 3);
                                Diag2Y = -GrainUnitVector(9 * MyOrientation + 4);
                                Diag2Z = -GrainUnitVector(9 * MyOrientation + 5);
                            }

                            double Angle3 = (GrainUnitVector(9 * MyOrientation + 6) * x0 +
                                             GrainUnitVector(9 * MyOrientation + 7) * y0 +
                                             GrainUnitVector(9 * MyOrientation + 8) * z0) /
                                            mag0;
                            if (Angle3 < 0) {
                                Diag3X = GrainUnitVector(9 * MyOrientation + 6);
                                Diag3Y = GrainUnitVector(9 * MyOrientation + 7);
                                Diag3Z = GrainUnitVector(9 * MyOrientation + 8);
                            }
                            else {
                                Diag3X = -GrainUnitVector(9 * MyOrientation + 6);
                                Diag3Y = -GrainUnitVector(9 * MyOrientation + 7);
                                Diag3Z = -GrainUnitVector(9 * MyOrientation + 8);
                            }

                            double U1[3], U2[3], UU[3], Norm[3];
                            U1[0] = Diag2X - Diag1X;
                            U1[1] = Diag2Y - Diag1Y;
                            U1[2] = Diag2Z - Diag1Z;
                            U2[0] = Diag3X - Diag1X;
                            U2[1] = Diag3Y - Diag1Y;
                            U2[2] = Diag3Z - Diag1Z;
                            UU[0] = U1[1] * U2[2] - U1[2] * U2[1];
                            UU[1] = U1[2] * U2[0] - U1[0] * U2[2];
                            UU[2] = U1[0] * U2[1] - U1[1] * U2[0];
                            double NDem = sqrt(UU[0] * UU[0] + UU[1] * UU[1] + UU[2] * UU[2]);
                            Norm[0] = UU[0] / NDem;
                            Norm[1] = UU[1] / NDem;
                            Norm[2] = UU[2] / NDem;
                            // normal to capturing plane
                            double normx = Norm[0];
                            double normy = Norm[1];
                            double normz = Norm[2];
                            double ParaT = (normx * x0 + normy * y0 + normz * z0) /
                                           (normx * Diag1X + normy * Diag1Y + normz * Diag1Z);
                            float CDLVal = pow(
                                pow(ParaT * Diag1X, 2.0) + pow(ParaT * Diag1Y, 2.0) + pow(ParaT * Diag1Z, 2.0), 0.5);
                            CritDiagonalLength((long int)(26) * D3D1ConvPosition + (long int)(n)) = CDLVal;
                        }
                    }
                }
                else if (CellType(D3D1ConvPositionGlobal) == TemporaryInit) {
                    // Mark and count the number of nucleation events to be sent to other ranks
                    GrainID(D3D1ConvPositionGlobal) = NCounter;
                    NCounter--;
                    if (np > 1) {
                        if (DecompositionStrategy == 1) {
                            if (RankY == 1) {
                                SouthNucleiSendCount++;
                            }
                            else if (RankY == MyYSlices - 2) {
                                NorthNucleiSendCount++;
                            }
                        }
                        else {
                            if (RankY == 1) {
                                // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                                if (RankX == MyXSlices - 2) {
                                    SouthEastNucleiSendCount++;
                                    EastNucleiSendCount++;
                                    SouthNucleiSendCount++;
                                }
                                else if (RankX == 1) {
                                    SouthWestNucleiSendCount++;
                                    WestNucleiSendCount++;
                                    SouthNucleiSendCount++;
                                }
                                else if ((RankX > 1) && (RankX < MyXSlices - 2)) {
                                    // This is being sent to MyLeft
                                    SouthNucleiSendCount++;
                                }
                            }
                            else if (RankY == MyYSlices - 2) {
                                // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                                if (RankX == MyXSlices - 2) {
                                    NorthEastNucleiSendCount++;
                                    EastNucleiSendCount++;
                                    NorthNucleiSendCount++;
                                }
                                else if (RankX == 1) {
                                    NorthWestNucleiSendCount++;
                                    WestNucleiSendCount++;
                                    NorthNucleiSendCount++;
                                }
                                else if ((RankX > 1) && (RankX < MyXSlices - 2)) {
                                    NorthNucleiSendCount++;
                                }
                            }
                            else if ((RankX == 1) && (RankY > 1) && (RankY < MyYSlices - 2)) {
                                WestNucleiSendCount++;
                            }
                            else if ((RankX == MyXSlices - 2) && (RankY > 1) && (RankY < MyYSlices - 2)) {
                                EastNucleiSendCount++;
                            }
                        }
                    }
                }
            }
        }
    }

    // Send/recieve number of nuclei in the ghost nodes so that each rank knows the total number of
    // nuclei in it's domain
    int SouthNucleiRecvCount = 0;
    int NorthNucleiRecvCount = 0;
    int EastNucleiRecvCount = 0;
    int WestNucleiRecvCount = 0;
    int NorthWestNucleiRecvCount = 0;
    int NorthEastNucleiRecvCount = 0;
    int SouthWestNucleiRecvCount = 0;
    int SouthEastNucleiRecvCount = 0;

    // South/North exchange number of nuclei in halo regions
    MPI_Sendrecv(&NorthNucleiSendCount, 1, MPI_INT, NeighborRank_North, 0, &SouthNucleiRecvCount, 1, MPI_INT,
                 NeighborRank_South, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&SouthNucleiSendCount, 1, MPI_INT, NeighborRank_South, 0, &NorthNucleiRecvCount, 1, MPI_INT,
                 NeighborRank_North, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (DecompositionStrategy != 1) {
        // East/West exchange number of nuclei in halo regions
        MPI_Sendrecv(&EastNucleiSendCount, 1, MPI_INT, NeighborRank_East, 0, &WestNucleiRecvCount, 1, MPI_INT,
                     NeighborRank_West, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&WestNucleiSendCount, 1, MPI_INT, NeighborRank_West, 0, &EastNucleiRecvCount, 1, MPI_INT,
                     NeighborRank_East, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // NorthWest/SouthEast exchange number of nuclei in halo regions
        MPI_Sendrecv(&NorthWestNucleiSendCount, 1, MPI_INT, NeighborRank_NorthWest, 0, &SouthEastNucleiRecvCount, 1,
                     MPI_INT, NeighborRank_SouthEast, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&SouthEastNucleiSendCount, 1, MPI_INT, NeighborRank_SouthEast, 0, &NorthWestNucleiRecvCount, 1,
                     MPI_INT, NeighborRank_NorthWest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // NorthEast/SouthWest exchange number of nuclei in halo regions
        MPI_Sendrecv(&NorthEastNucleiSendCount, 1, MPI_INT, NeighborRank_NorthEast, 0, &SouthWestNucleiRecvCount, 1,
                     MPI_INT, NeighborRank_SouthWest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&SouthWestNucleiSendCount, 1, MPI_INT, NeighborRank_SouthWest, 0, &NorthEastNucleiRecvCount, 1,
                     MPI_INT, NeighborRank_NorthEast, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if (NeighborRank_South == MPI_PROC_NULL)
        SouthNucleiRecvCount = 0;
    if (NeighborRank_North == MPI_PROC_NULL)
        NorthNucleiRecvCount = 0;
    if (NeighborRank_East == MPI_PROC_NULL)
        EastNucleiRecvCount = 0;
    if (NeighborRank_West == MPI_PROC_NULL)
        WestNucleiRecvCount = 0;
    if (NeighborRank_NorthWest == MPI_PROC_NULL)
        NorthWestNucleiRecvCount = 0;
    if (NeighborRank_NorthEast == MPI_PROC_NULL)
        NorthEastNucleiRecvCount = 0;
    if (NeighborRank_SouthWest == MPI_PROC_NULL)
        SouthWestNucleiRecvCount = 0;
    if (NeighborRank_SouthEast == MPI_PROC_NULL)
        SouthEastNucleiRecvCount = 0;

    PossibleNuclei_ThisRank +=
        (SouthNucleiRecvCount + NorthNucleiRecvCount + EastNucleiRecvCount + WestNucleiRecvCount +
         NorthWestNucleiRecvCount + NorthEastNucleiRecvCount + SouthWestNucleiRecvCount + SouthEastNucleiRecvCount);

    // Remove delay cells not bordering others
    for (int RankZ = 1; RankZ < nz; RankZ++) {
        for (int RankX = 1; RankX < MyXSlices - 1; RankX++) {
            for (int RankY = 1; RankY < MyYSlices - 1; RankY++) {
                int D3D1ConvPosition = RankZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                if (CellType(D3D1ConvPosition) == Liquid) {
                    // Check to see if this cell is at the interface
                    int LCount = 0;
                    // Which neighbors should be iterated over?
                    int ItBounds = FindItBounds(RankX, RankY, MyXSlices, MyYSlices);
                    int NListLength;
                    if (ItBounds == 0) {
                        NListLength = 26;
                    }
                    else if (ItBounds > 4) {
                        NListLength = 11;
                    }
                    else {
                        NListLength = 17;
                    }
                    // "ll" corresponds to the specific position on the list of neighboring cells
                    for (int ll = 0; ll < NListLength; ll++) {
                        // "l" correpsponds to the specific neighboring cell
                        int l = ItList(ItBounds, ll);
                        // Local coordinates of adjacent cell center
                        int MyNeighborX = RankX + NeighborX(l);
                        int MyNeighborY = RankY + NeighborY(l);
                        int MyNeighborZ = RankZ + NeighborZ(l);
                        if (MyNeighborZ < nz) {
                            int NeighborD3D1ConvPosition =
                                MyNeighborZ * MyXSlices * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
                            if ((CellType(NeighborD3D1ConvPosition) != Solid) &&
                                (CellType(NeighborD3D1ConvPosition) != Wall)) {
                                LCount++;
                            }
                        }
                    }
                    if (LCount == 0) {
                        // This cell is returned to solid type
                        CellType(D3D1ConvPosition) = Solid;
                    }
                }
            }
        }
    }
}

//*****************************************************************************/
// After initializing grain structure and filling ghost nodes, the known potential nucleation sites are placed into the
// nucleation data structures Each nucleation event is assigned a time step, beyond which if the associated cell is not
// solid or actve, the event occurs This data is synced across MPI ranks, for nucleation events that occur in the ghost
// nodes
void NucleiInit(int DecompositionStrategy, int MyXSlices, int MyYSlices, int nz, int id, double dTN, double dTsigma,
                int NeighborRank_North, int NeighborRank_South, int NeighborRank_East, int NeighborRank_West,
                int NeighborRank_NorthEast, int NeighborRank_NorthWest, int NeighborRank_SouthEast,
                int NeighborRank_SouthWest, ViewI_H NucleiLocation, ViewI_H NucleationTimes, ViewI_H CellType,
                ViewI_H GrainID, ViewI_H CritTimeStep, ViewF_H UndercoolingChange) {

    // Counts and buffers for sending/recieving nucleation data from ghost nodes
    int SouthNucleiSendCount = 0;
    int NorthNucleiSendCount = 0;
    int EastNucleiSendCount = 0;
    int WestNucleiSendCount = 0;
    int NorthWestNucleiSendCount = 0;
    int NorthEastNucleiSendCount = 0;
    int SouthWestNucleiSendCount = 0;
    int SouthEastNucleiSendCount = 0;
    std::vector<int> TempNucleiDataSouth, TempNucleiDataNorth, TempNucleiDataEast, TempNucleiDataWest,
        TempNucleiDataNorthWest, TempNucleiDataNorthEast, TempNucleiDataSouthWest, TempNucleiDataSouthEast;

    // Gaussian distribution of nucleation undercooling
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(dTN, dTsigma);
    for (int i = 0; i < 120 * id; i++) {
        distribution(generator);
    }

    // Collect data for ghost nodes' nucleation events
    int NEvent = 0;
    for (int i = 0; i < MyXSlices * MyYSlices * nz; i++) {
        if (CellType(i) == TemporaryInit) {
            NucleiLocation(NEvent) = i;
            // Undercooling for this nucleation event
            double LocNucUnd = distribution(generator);
            // Time steps to reach this undercooling after cell goes below the liquidus
            int TimeToNucUnd = CritTimeStep(i) + round(LocNucUnd / UndercoolingChange(i));
            NucleationTimes(NEvent) = std::max(CritTimeStep(i), TimeToNucUnd);

            // Determine if other MPI ranks need information about this potential nucleation event
            // If so, store the location (X,Y,Z), GrainID, and nucleation time step value to be sent
            int RankZ = i / (MyXSlices * MyYSlices);
            int Rem = i % (MyXSlices * MyYSlices);
            int RankX = Rem / MyYSlices;
            int RankY = Rem % MyYSlices;
            if (DecompositionStrategy == 1) {
                if (RankY == 1) {
                    SouthNucleiSendCount++;
                    TempNucleiDataSouth.push_back(RankX);
                    TempNucleiDataSouth.push_back(RankZ);
                    TempNucleiDataSouth.push_back(TimeToNucUnd);
                    TempNucleiDataSouth.push_back(GrainID(i));
                }
                else if (RankY == MyYSlices - 2) {
                    NorthNucleiSendCount++;
                    TempNucleiDataNorth.push_back(RankX);
                    TempNucleiDataNorth.push_back(RankZ);
                    TempNucleiDataNorth.push_back(TimeToNucUnd);
                    TempNucleiDataNorth.push_back(GrainID(i));
                }
            }
            else {
                if (RankY == 1) {
                    // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                    if (RankX == MyXSlices - 2) {

                        SouthEastNucleiSendCount++;
                        TempNucleiDataSouthEast.push_back(RankZ);
                        TempNucleiDataSouthEast.push_back(TimeToNucUnd);
                        TempNucleiDataSouthEast.push_back(GrainID(i));

                        EastNucleiSendCount++;
                        TempNucleiDataEast.push_back(RankY);
                        TempNucleiDataEast.push_back(RankZ);
                        TempNucleiDataEast.push_back(TimeToNucUnd);
                        TempNucleiDataEast.push_back(GrainID(i));

                        SouthNucleiSendCount++;
                        TempNucleiDataSouth.push_back(RankX);
                        TempNucleiDataSouth.push_back(RankZ);
                        TempNucleiDataSouth.push_back(TimeToNucUnd);
                        TempNucleiDataSouth.push_back(GrainID(i));
                    }
                    else if (RankX == 1) {

                        SouthWestNucleiSendCount++;
                        TempNucleiDataSouthWest.push_back(RankZ);
                        TempNucleiDataSouthWest.push_back(TimeToNucUnd);
                        TempNucleiDataSouthWest.push_back(GrainID(i));

                        WestNucleiSendCount++;
                        TempNucleiDataWest.push_back(RankY);
                        TempNucleiDataWest.push_back(RankZ);
                        TempNucleiDataWest.push_back(TimeToNucUnd);
                        TempNucleiDataWest.push_back(GrainID(i));

                        SouthNucleiSendCount++;
                        TempNucleiDataSouth.push_back(RankX);
                        TempNucleiDataSouth.push_back(RankZ);
                        TempNucleiDataSouth.push_back(TimeToNucUnd);
                        TempNucleiDataSouth.push_back(GrainID(i));
                    }
                    else if ((RankX > 1) && (RankX < MyXSlices - 2)) {
                        // This is being sent to MyLeft
                        SouthNucleiSendCount++;
                        TempNucleiDataSouth.push_back(RankX);
                        TempNucleiDataSouth.push_back(RankZ);
                        TempNucleiDataSouth.push_back(TimeToNucUnd);
                        TempNucleiDataSouth.push_back(GrainID(i));
                    }
                }
                else if (RankY == MyYSlices - 2) {
                    // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                    if (RankX == MyXSlices - 2) {
                        NorthEastNucleiSendCount++;
                        TempNucleiDataNorthEast.push_back(RankZ);
                        TempNucleiDataNorthEast.push_back(TimeToNucUnd);
                        TempNucleiDataNorthEast.push_back(GrainID(i));

                        EastNucleiSendCount++;
                        TempNucleiDataEast.push_back(RankY);
                        TempNucleiDataEast.push_back(RankZ);
                        TempNucleiDataEast.push_back(TimeToNucUnd);
                        TempNucleiDataEast.push_back(GrainID(i));

                        NorthNucleiSendCount++;
                        TempNucleiDataNorth.push_back(RankX);
                        TempNucleiDataNorth.push_back(RankZ);
                        TempNucleiDataNorth.push_back(TimeToNucUnd);
                        TempNucleiDataNorth.push_back(GrainID(i));
                    }
                    else if (RankX == 1) {
                        NorthWestNucleiSendCount++;
                        TempNucleiDataNorthWest.push_back(RankZ);
                        TempNucleiDataNorthWest.push_back(TimeToNucUnd);
                        TempNucleiDataNorthWest.push_back(GrainID(i));

                        WestNucleiSendCount++;
                        TempNucleiDataWest.push_back(RankY);
                        TempNucleiDataWest.push_back(RankZ);
                        TempNucleiDataWest.push_back(TimeToNucUnd);
                        TempNucleiDataWest.push_back(GrainID(i));

                        NorthNucleiSendCount++;
                        TempNucleiDataNorth.push_back(RankX);
                        TempNucleiDataNorth.push_back(RankZ);
                        TempNucleiDataNorth.push_back(TimeToNucUnd);
                        TempNucleiDataNorth.push_back(GrainID(i));
                    }
                    else if ((RankX > 1) && (RankX < MyXSlices - 2)) {
                        NorthNucleiSendCount++;
                        TempNucleiDataNorth.push_back(RankX);
                        TempNucleiDataNorth.push_back(RankZ);
                        TempNucleiDataNorth.push_back(TimeToNucUnd);
                        TempNucleiDataNorth.push_back(GrainID(i));
                    }
                }
                else if ((RankX == 1) && (RankY > 1) && (RankY < MyYSlices - 2)) {
                    WestNucleiSendCount++;
                    TempNucleiDataWest.push_back(RankY);
                    TempNucleiDataWest.push_back(RankZ);
                    TempNucleiDataWest.push_back(TimeToNucUnd);
                    TempNucleiDataWest.push_back(GrainID(i));
                }
                else if ((RankX == MyXSlices - 2) && (RankY > 1) && (RankY < MyYSlices - 2)) {
                    EastNucleiSendCount++;
                    TempNucleiDataEast.push_back(RankY);
                    TempNucleiDataEast.push_back(RankZ);
                    TempNucleiDataEast.push_back(TimeToNucUnd);
                    TempNucleiDataEast.push_back(GrainID(i));
                }
            }
            NEvent++;
        }
    }

    // Determine whether or not ghost node information transfer needs to take place
    int SouthNucleiRecvCount = 0;
    int NorthNucleiRecvCount = 0;
    int EastNucleiRecvCount = 0;
    int WestNucleiRecvCount = 0;
    int NorthWestNucleiRecvCount = 0;
    int NorthEastNucleiRecvCount = 0;
    int SouthWestNucleiRecvCount = 0;
    int SouthEastNucleiRecvCount = 0;

    // South/North exchange number of nuclei in halo regions
    MPI_Sendrecv(&NorthNucleiSendCount, 1, MPI_INT, NeighborRank_North, 0, &SouthNucleiRecvCount, 1, MPI_INT,
                 NeighborRank_South, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&SouthNucleiSendCount, 1, MPI_INT, NeighborRank_South, 0, &NorthNucleiRecvCount, 1, MPI_INT,
                 NeighborRank_North, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (DecompositionStrategy != 1) {
        // East/West exchange number of nuclei in halo regions
        MPI_Sendrecv(&EastNucleiSendCount, 1, MPI_INT, NeighborRank_East, 0, &WestNucleiRecvCount, 1, MPI_INT,
                     NeighborRank_West, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&WestNucleiSendCount, 1, MPI_INT, NeighborRank_West, 0, &EastNucleiRecvCount, 1, MPI_INT,
                     NeighborRank_East, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // NorthWest/SouthEast exchange number of nuclei in halo regions
        MPI_Sendrecv(&NorthWestNucleiSendCount, 1, MPI_INT, NeighborRank_NorthWest, 0, &SouthEastNucleiRecvCount, 1,
                     MPI_INT, NeighborRank_SouthEast, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&SouthEastNucleiSendCount, 1, MPI_INT, NeighborRank_SouthEast, 0, &NorthWestNucleiRecvCount, 1,
                     MPI_INT, NeighborRank_NorthWest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // NorthEast/SouthWest exchange number of nuclei in halo regions
        MPI_Sendrecv(&NorthEastNucleiSendCount, 1, MPI_INT, NeighborRank_NorthEast, 0, &SouthWestNucleiRecvCount, 1,
                     MPI_INT, NeighborRank_SouthWest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&SouthWestNucleiSendCount, 1, MPI_INT, NeighborRank_SouthWest, 0, &NorthEastNucleiRecvCount, 1,
                     MPI_INT, NeighborRank_NorthEast, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if (NeighborRank_South == MPI_PROC_NULL)
        SouthNucleiRecvCount = 0;
    if (NeighborRank_North == MPI_PROC_NULL)
        NorthNucleiRecvCount = 0;
    if (NeighborRank_East == MPI_PROC_NULL)
        EastNucleiRecvCount = 0;
    if (NeighborRank_West == MPI_PROC_NULL)
        WestNucleiRecvCount = 0;
    if (NeighborRank_NorthWest == MPI_PROC_NULL)
        NorthWestNucleiRecvCount = 0;
    if (NeighborRank_NorthEast == MPI_PROC_NULL)
        NorthEastNucleiRecvCount = 0;
    if (NeighborRank_SouthWest == MPI_PROC_NULL)
        SouthWestNucleiRecvCount = 0;
    if (NeighborRank_SouthEast == MPI_PROC_NULL)
        SouthEastNucleiRecvCount = 0;

    // Buffers for recieving nuclei data
    ViewI_H NucleiRecvBufferSouth("NucleiRecvBufferSouth", 4 * SouthNucleiRecvCount);
    ViewI_H NucleiRecvBufferNorth("NucleiRecvBufferNorth", 4 * NorthNucleiRecvCount);
    ViewI_H NucleiRecvBufferEast("NucleiRecvBufferEast", 4 * EastNucleiRecvCount);
    ViewI_H NucleiRecvBufferWest("NucleiRecvBufferWest", 4 * WestNucleiRecvCount);
    ViewI_H NucleiRecvBufferNorthWest("NucleiRecvBufferNorthWest", 3 * NorthWestNucleiRecvCount);
    ViewI_H NucleiRecvBufferNorthEast("NucleiRecvBufferNorthEast", 3 * NorthEastNucleiRecvCount);
    ViewI_H NucleiRecvBufferSouthWest("NucleiRecvBufferSouthWest", 3 * SouthWestNucleiRecvCount);
    ViewI_H NucleiRecvBufferSouthEast("NucleiRecvBufferSouthEast", 3 * SouthEastNucleiRecvCount);

    // Collect ghost node data and send to other ranks- south and north
    if (SouthNucleiSendCount > 0) {
        ViewI_H NucleiSendBufferSouth("NucleiSendBufferSouth", 4 * SouthNucleiSendCount);
        for (int i = 0; i < 4 * SouthNucleiSendCount; i++) {
            NucleiSendBufferSouth(i) = TempNucleiDataSouth[i];
        }
        if (NorthNucleiRecvCount == 0) {
            // Sending data to id = id - 1 only
            MPI_Send(NucleiSendBufferSouth.data(), SouthNucleiSendCount * 4, MPI_INT, NeighborRank_South, 0,
                     MPI_COMM_WORLD);
        }
        else {
            // Sending data to id = id - 1 and recieving data from id = id + 1
            MPI_Sendrecv(NucleiSendBufferSouth.data(), SouthNucleiSendCount * 4, MPI_INT, NeighborRank_South, 0,
                         NucleiRecvBufferNorth.data(), NorthNucleiRecvCount * 4, MPI_INT, NeighborRank_North, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    else if (NorthNucleiRecvCount > 0) {
        // Recieving data from id = id + 1 only
        MPI_Recv(NucleiRecvBufferNorth.data(), NorthNucleiRecvCount * 4, MPI_INT, NeighborRank_North, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    }

    if (NorthNucleiSendCount > 0) {
        ViewI_H NucleiSendBufferNorth("NucleiSendBufferNorth", 4 * NorthNucleiSendCount);
        for (int i = 0; i < 4 * NorthNucleiSendCount; i++) {
            NucleiSendBufferNorth(i) = TempNucleiDataNorth[i];
        }
        if (SouthNucleiRecvCount == 0) {
            // Sending data to id = id + 1 only
            MPI_Send(NucleiSendBufferNorth.data(), NorthNucleiSendCount * 4, MPI_INT, NeighborRank_North, 1,
                     MPI_COMM_WORLD);
        }
        else {
            // Sending data to id = id + 1 and recieving data from id = id - 1
            MPI_Sendrecv(NucleiSendBufferNorth.data(), NorthNucleiSendCount * 4, MPI_INT, NeighborRank_North, 1,
                         NucleiRecvBufferSouth.data(), SouthNucleiRecvCount * 4, MPI_INT, NeighborRank_South, 1,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    else if (SouthNucleiRecvCount > 0) {
        // Recieving data from id = id - 1 only
        MPI_Recv(NucleiRecvBufferSouth.data(), SouthNucleiRecvCount * 4, MPI_INT, NeighborRank_South, 1, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    }

    if (DecompositionStrategy != 1) {
        // Collect ghost node data and send to other ranks- East and West
        if (EastNucleiSendCount > 0) {
            ViewI_H NucleiSendBufferEast("NucleiSendBufferEast", 4 * EastNucleiSendCount);
            for (int i = 0; i < 4 * EastNucleiSendCount; i++) {
                NucleiSendBufferEast(i) = TempNucleiDataEast[i];
            }

            if (WestNucleiRecvCount == 0) {
                // Sending data only
                MPI_Send(NucleiSendBufferEast.data(), EastNucleiSendCount * 4, MPI_INT, NeighborRank_East, 0,
                         MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
                MPI_Sendrecv(NucleiSendBufferEast.data(), EastNucleiSendCount * 4, MPI_INT, NeighborRank_East, 0,
                             NucleiRecvBufferWest.data(), WestNucleiRecvCount * 4, MPI_INT, NeighborRank_West, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        else if (WestNucleiRecvCount > 0) {
            // Recieving data only
            MPI_Recv(NucleiRecvBufferWest.data(), WestNucleiRecvCount * 4, MPI_INT, NeighborRank_West, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (WestNucleiSendCount > 0) {
            ViewI_H NucleiSendBufferWest("NucleiSendBufferWest", 4 * WestNucleiSendCount);
            for (int i = 0; i < 4 * WestNucleiSendCount; i++) {
                NucleiSendBufferWest(i) = TempNucleiDataWest[i];
            }
            if (EastNucleiRecvCount == 0) {
                // Sending data only
                MPI_Send(NucleiSendBufferWest.data(), WestNucleiSendCount * 4, MPI_INT, NeighborRank_West, 1,
                         MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
                MPI_Sendrecv(NucleiSendBufferWest.data(), WestNucleiSendCount * 4, MPI_INT, NeighborRank_West, 1,
                             NucleiRecvBufferEast.data(), EastNucleiRecvCount * 4, MPI_INT, NeighborRank_East, 1,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        else if (EastNucleiRecvCount > 0) {
            // Recieving data only
            MPI_Recv(NucleiRecvBufferEast.data(), EastNucleiRecvCount * 4, MPI_INT, NeighborRank_East, 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Collect ghost node data and send to other ranks- NorthWest and SouthEast
        if (NorthWestNucleiSendCount > 0) {
            ViewI_H NucleiSendBufferNorthWest("NucleiSendBufferNorthWest", 3 * NorthWestNucleiSendCount);
            for (int i = 0; i < 3 * NorthWestNucleiSendCount; i++) {
                NucleiSendBufferNorthWest(i) = TempNucleiDataNorthWest[i];
            }
            if (SouthEastNucleiRecvCount == 0) {
                // Sending data only
                MPI_Send(NucleiSendBufferNorthWest.data(), NorthWestNucleiSendCount * 3, MPI_INT,
                         NeighborRank_NorthWest, 0, MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
                MPI_Sendrecv(NucleiSendBufferNorthWest.data(), NorthWestNucleiSendCount * 3, MPI_INT,
                             NeighborRank_NorthWest, 0, NucleiRecvBufferSouthEast.data(), SouthEastNucleiRecvCount * 3,
                             MPI_INT, NeighborRank_SouthEast, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        else if (SouthEastNucleiRecvCount > 0) {
            // Recieving data only
            MPI_Recv(NucleiRecvBufferSouthEast.data(), SouthEastNucleiRecvCount * 3, MPI_INT, NeighborRank_SouthEast, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (SouthEastNucleiSendCount > 0) {
            ViewI_H NucleiSendBufferSouthEast("bufferH", 3 * SouthEastNucleiSendCount);
            for (int i = 0; i < 3 * SouthEastNucleiSendCount; i++) {
                NucleiSendBufferSouthEast(i) = TempNucleiDataSouthEast[i];
            }
            if (NorthWestNucleiRecvCount == 0) {
                // Sending data only
                MPI_Send(NucleiSendBufferSouthEast.data(), SouthEastNucleiSendCount * 3, MPI_INT,
                         NeighborRank_SouthEast, 0, MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
                MPI_Sendrecv(NucleiSendBufferSouthEast.data(), SouthEastNucleiSendCount * 3, MPI_INT,
                             NeighborRank_SouthEast, 0, NucleiRecvBufferNorthWest.data(), NorthWestNucleiRecvCount * 3,
                             MPI_INT, NeighborRank_NorthWest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        else if (NorthWestNucleiRecvCount > 0) {
            // Recieving data only
            MPI_Recv(NucleiRecvBufferNorthWest.data(), NorthWestNucleiRecvCount * 3, MPI_INT, NeighborRank_NorthWest, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Collect ghost node data and send to other ranks- NorthEast and SouthWest
        if (NorthEastNucleiSendCount > 0) {
            ViewI_H NucleiSendBufferNorthEast("NucleiSendBufferNorthEast", 4 * NorthEastNucleiSendCount);
            for (int i = 0; i < 4 * NorthEastNucleiSendCount; i++) {
                NucleiSendBufferNorthEast(i) = TempNucleiDataNorthEast[i];
            }
            if (SouthWestNucleiRecvCount == 0) {
                // Sending data only
                MPI_Send(NucleiSendBufferNorthEast.data(), NorthEastNucleiSendCount * 3, MPI_INT,
                         NeighborRank_NorthEast, 1, MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
                MPI_Sendrecv(NucleiSendBufferNorthEast.data(), NorthEastNucleiSendCount * 3, MPI_INT,
                             NeighborRank_NorthEast, 1, NucleiRecvBufferSouthWest.data(), SouthWestNucleiRecvCount * 3,
                             MPI_INT, NeighborRank_SouthWest, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        else if (SouthWestNucleiRecvCount > 0) {
            // Recieving data only
            MPI_Recv(NucleiRecvBufferSouthWest.data(), SouthWestNucleiRecvCount * 3, MPI_INT, NeighborRank_SouthWest, 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (SouthWestNucleiSendCount > 0) {
            ViewI_H NucleiSendBufferSouthWest("bufferG", 3 * SouthWestNucleiSendCount);
            for (int i = 0; i < 3 * SouthWestNucleiSendCount; i++) {
                NucleiSendBufferSouthWest(i) = TempNucleiDataSouthWest[i];
            }
            if (NorthEastNucleiRecvCount == 0) {
                // Sending data only
                MPI_Send(NucleiSendBufferSouthWest.data(), SouthWestNucleiSendCount * 3, MPI_INT,
                         NeighborRank_SouthWest, 1, MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
                MPI_Sendrecv(NucleiRecvBufferSouthWest.data(), SouthWestNucleiSendCount * 3, MPI_INT,
                             NeighborRank_SouthWest, 1, NucleiRecvBufferNorthEast.data(), NorthEastNucleiRecvCount * 3,
                             MPI_INT, NeighborRank_NorthEast, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        else if (NorthEastNucleiRecvCount > 0) {
            // Recieving data only
            MPI_Recv(NucleiRecvBufferNorthEast.data(), NorthEastNucleiRecvCount * 3, MPI_INT, NeighborRank_NorthEast, 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // Place ghost node data recieved from the south (if needed)
    if (SouthNucleiRecvCount > 0) {
        for (int i = 0; i < SouthNucleiRecvCount; i++) {
            int RankX = NucleiRecvBufferSouth(4 * i);
            int RankY = 0;
            int RankZ = NucleiRecvBufferSouth(4 * i + 1);
            int CellLocation = RankZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
            NucleiLocation(NEvent) = CellLocation;
            NucleationTimes(NEvent) = NucleiRecvBufferSouth(4 * i + 2);
            CellType(CellLocation) = TemporaryInit;
            GrainID(CellLocation) = NucleiRecvBufferSouth(4 * i + 3);
            NEvent++;
        }
    }

    // Place ghost node data recieved from the north (if needed)
    if (NorthNucleiRecvCount > 0) {
        for (int i = 0; i < NorthNucleiRecvCount; i++) {
            int RankX = NucleiRecvBufferNorth(4 * i);
            int RankY = MyYSlices - 1;
            int RankZ = NucleiRecvBufferNorth(4 * i + 1);
            int CellLocation = RankZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
            NucleiLocation(NEvent) = CellLocation;
            NucleationTimes(NEvent) = NucleiRecvBufferNorth(4 * i + 2);
            CellType(CellLocation) = TemporaryInit;
            GrainID(CellLocation) = NucleiRecvBufferNorth(4 * i + 3);
            NEvent++;
        }
    }

    if (DecompositionStrategy != 1) {
        // Place ghost node data recieved from the east (if needed)
        if (EastNucleiRecvCount > 0) {
            for (int i = 0; i < EastNucleiRecvCount; i++) {
                int RankX = MyXSlices - 1;
                int RankY = NucleiRecvBufferEast(4 * i);
                int RankZ = NucleiRecvBufferEast(4 * i + 1);
                int CellLocation = RankZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                NucleiLocation(NEvent) = CellLocation;
                NucleationTimes(NEvent) = NucleiRecvBufferEast(4 * i + 2);
                CellType(CellLocation) = TemporaryInit;
                GrainID(CellLocation) = NucleiRecvBufferEast(4 * i + 3);
                NEvent++;
            }
        }

        // Place ghost node data recieved from the west (if needed)
        if (WestNucleiRecvCount > 0) {
            for (int i = 0; i < WestNucleiRecvCount; i++) {
                int RankX = 0;
                int RankY = NucleiRecvBufferWest(4 * i);
                int RankZ = NucleiRecvBufferWest(4 * i + 1);
                int CellLocation = RankZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                NucleiLocation(NEvent) = CellLocation;
                NucleationTimes(NEvent) = NucleiRecvBufferWest(4 * i + 2);
                CellType(CellLocation) = TemporaryInit;
                GrainID(CellLocation) = NucleiRecvBufferWest(4 * i + 3);
                NEvent++;
            }
        }

        // Place ghost node data recieved from the northwest (if needed)
        if (NorthWestNucleiRecvCount > 0) {
            for (int i = 0; i < NorthWestNucleiRecvCount; i++) {
                int RankX = 0;
                int RankY = MyYSlices - 1;
                int RankZ = NucleiRecvBufferNorthWest(3 * i);
                int CellLocation = RankZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                NucleiLocation(NEvent) = CellLocation;
                NucleationTimes(NEvent) = NucleiRecvBufferNorthWest(3 * i + 1);
                CellType(CellLocation) = TemporaryInit;
                GrainID(CellLocation) = NucleiRecvBufferNorthWest(3 * i + 2);
                NEvent++;
            }
        }

        // Place ghost node data recieved from the northeast (if needed)
        if (NorthEastNucleiRecvCount > 0) {
            for (int i = 0; i < NorthEastNucleiRecvCount; i++) {
                int RankX = MyXSlices - 1;
                int RankY = MyYSlices - 1;
                int RankZ = NucleiRecvBufferNorthEast(3 * i);
                int CellLocation = RankZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                NucleiLocation(NEvent) = CellLocation;
                NucleationTimes(NEvent) = NucleiRecvBufferNorthEast(3 * i + 1);
                CellType(CellLocation) = TemporaryInit;
                GrainID(CellLocation) = NucleiRecvBufferNorthEast(3 * i + 2);
                NEvent++;
            }
        }

        // Place ghost node data recieved from the southwest (if needed)
        if (SouthWestNucleiRecvCount > 0) {
            for (int i = 0; i < SouthWestNucleiRecvCount; i++) {
                int RankX = 0;
                int RankY = 0;
                int RankZ = NucleiRecvBufferSouthWest(3 * i);
                int CellLocation = RankZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                NucleiLocation(NEvent) = CellLocation;
                NucleationTimes(NEvent) = NucleiRecvBufferSouthWest(3 * i + 1);
                CellType(CellLocation) = TemporaryInit;
                GrainID(CellLocation) = NucleiRecvBufferSouthWest(3 * i + 2);
                NEvent++;
            }
        }

        // Place ghost node data recieved from the southeast (if needed)
        if (SouthEastNucleiRecvCount > 0) {
            for (int i = 0; i < SouthEastNucleiRecvCount; i++) {
                int RankX = MyXSlices - 1;
                int RankY = 0;
                int RankZ = NucleiRecvBufferSouthEast(3 * i);
                int CellLocation = RankZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                NucleiLocation(NEvent) = CellLocation;
                NucleationTimes(NEvent) = NucleiRecvBufferSouthEast(3 * i + 1);
                CellType(CellLocation) = TemporaryInit;
                GrainID(CellLocation) = NucleiRecvBufferSouthEast(3 * i + 2);
                NEvent++;
            }
        }
    }

    // Replace the temporary "TemporaryInit" cell type with the "Liquid" cell type - as GrainID is used to mark the
    // potential nucleus orientation, TemporaryInit is unneeded
    for (int i = 0; i < MyXSlices * MyYSlices * nz; i++) {
        if (CellType(i) == TemporaryInit)
            CellType(i) = Liquid;
    }
    std::cout << "(" << id << ": " << NEvent << ") " << std::flush;
}

//*****************************************************************************/
void DomainShiftAndResize(int id, int MyXSlices, int MyYSlices, int &ZShift, int &ZBound_Low, int &ZBound_High,
                          int &nzActive, int LocalDomainSize, int &LocalActiveDomainSize, int &BufSizeZ,
                          int LayerHeight, ViewI CellType, int layernumber, ViewI LayerID) {

    int ZBound_LowOld = ZBound_Low;

    // The top "top" of the active domain is a shift of "LayerHeight" from the previous domain top
    ZBound_High += LayerHeight;

    // The new "bottom" of the active domain is located just below the lowest active cells remaining in the domain
    int NewMin;
    Kokkos::parallel_reduce(
        "MinReduce", LocalDomainSize,
        KOKKOS_LAMBDA(const int &D3D1ConvPosition, int &lmin) {
            if (CellType(D3D1ConvPosition) == Active) {
                // Check Z position of this active cell
                int RankZ = D3D1ConvPosition / (MyXSlices * MyYSlices);
                if (RankZ < lmin)
                    lmin = RankZ;
            }
        },
        Kokkos::Min<int>(NewMin));
    NewMin--;
    MPI_Allreduce(&NewMin, &ZBound_Low, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    // Shift in +Z direction for the bottom of the active region
    ZShift = ZBound_Low - ZBound_LowOld;

    if (id == 0)
        std::cout << "New domain bottom at Z = " << ZBound_Low << ", a shift of " << ZShift << " cells" << std::endl;

    // The new "top" of the active domain is located at the highest location with cells solidifying during the next
    // layer
    int NewMax;
    Kokkos::parallel_reduce(
        "MaxReduce", LocalDomainSize,
        KOKKOS_LAMBDA(const int &D3D1ConvPosition, int &lmax) {
            if (LayerID(D3D1ConvPosition) == layernumber + 1) {
                // Check Z position of this active cell
                int RankZ = D3D1ConvPosition / (MyXSlices * MyYSlices);
                if (RankZ > lmax)
                    lmax = RankZ;
            }
        },
        Kokkos::Max<int>(NewMax));
    MPI_Allreduce(&NewMax, &ZBound_High, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "New domain top at Z = " << ZBound_High << std::endl;
    // Change in active region data structures' sizes
    nzActive = ZBound_High - ZBound_Low + 1;
    LocalActiveDomainSize = MyXSlices * MyYSlices * nzActive;

    // Change in height of buffers
    BufSizeZ = nzActive;

    if (id == 0)
        std::cout << "New active domain height is " << nzActive << std::endl;
}

//*****************************************************************************/
void LayerSetup(int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int LocalActiveDomainSize,
                ViewI GrainOrientation, int NGrainOrientations, ViewF GrainUnitVector, ViewI NeighborX, ViewI NeighborY,
                ViewI NeighborZ, ViewF DiagonalLength, ViewI CellType, ViewI GrainID, ViewF CritDiagonalLength,
                ViewF DOCenter, int DecompositionStrategy, Buffer2D BufferWestSend, Buffer2D BufferEastSend,
                Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, Buffer2D BufferNorthEastSend,
                Buffer2D BufferNorthWestSend, Buffer2D BufferSouthEastSend, Buffer2D BufferSouthWestSend,
                Buffer2D BufferWestRecv, Buffer2D BufferEastRecv, Buffer2D BufferNorthRecv, Buffer2D BufferSouthRecv,
                Buffer2D BufferNorthEastRecv, Buffer2D BufferNorthWestRecv, Buffer2D BufferSouthEastRecv,
                Buffer2D BufferSouthWestRecv, int &ZBound_Low) {

    // Reset active cell data structures
    Kokkos::deep_copy(DiagonalLength, 0);
    Kokkos::deep_copy(DOCenter, 0);
    Kokkos::deep_copy(CritDiagonalLength, 0);
    Kokkos::deep_copy(BufferSouthSend, 0.0);
    Kokkos::deep_copy(BufferSouthRecv, 0.0);
    Kokkos::deep_copy(BufferNorthSend, 0.0);
    Kokkos::deep_copy(BufferNorthRecv, 0.0);
    if (DecompositionStrategy != 1) {
        Kokkos::deep_copy(BufferEastSend, 0.0);
        Kokkos::deep_copy(BufferEastRecv, 0.0);
        Kokkos::deep_copy(BufferWestSend, 0.0);
        Kokkos::deep_copy(BufferWestRecv, 0.0);
        Kokkos::deep_copy(BufferNorthWestSend, 0.0);
        Kokkos::deep_copy(BufferNorthWestRecv, 0.0);
        Kokkos::deep_copy(BufferNorthEastSend, 0.0);
        Kokkos::deep_copy(BufferNorthEastRecv, 0.0);
        Kokkos::deep_copy(BufferSouthWestSend, 0.0);
        Kokkos::deep_copy(BufferSouthWestRecv, 0.0);
        Kokkos::deep_copy(BufferSouthEastSend, 0.0);
        Kokkos::deep_copy(BufferSouthEastRecv, 0.0);
    }

    Kokkos::parallel_for(
        "NewActiveCellInit", LocalActiveDomainSize, KOKKOS_LAMBDA(const int &D3D1ConvPosition) {
            // Initialize active cell data structures for those that are now part of the active domain
            int RankZ = D3D1ConvPosition / (MyXSlices * MyYSlices);
            int Rem = D3D1ConvPosition % (MyXSlices * MyYSlices);
            int RankX = Rem / MyYSlices;
            int RankY = Rem % MyYSlices;
            int GlobalZ = RankZ + ZBound_Low;
            int GlobalD3D1ConvPosition = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
            if (CellType(GlobalD3D1ConvPosition) == Active) {

                int GlobalX = RankX + MyXOffset;
                int GlobalY = RankY + MyYOffset;
                int MyGrainID = GrainID(GlobalD3D1ConvPosition);
                DiagonalLength(D3D1ConvPosition) = 0.01;
                DOCenter(3 * D3D1ConvPosition) = GlobalX + 0.5;
                DOCenter(3 * D3D1ConvPosition + 1) = GlobalY + 0.5;
                DOCenter(3 * D3D1ConvPosition + 2) = GlobalZ + 0.5;

                // The orientation for the new grain will depend on its Grain ID
                int MyOrientation = GrainOrientation(((abs(MyGrainID) - 1) % NGrainOrientations));
                // Calculate critical values at which this active cell leads to the activation of a neighboring liquid
                // cell (xp,yp,zp) is the new cell's center on the global grid
                double xp = GlobalX + 0.5;
                double yp = GlobalY + 0.5;
                double zp = GlobalZ + 0.5;

                float cx = DOCenter((long int)(3 * D3D1ConvPosition));
                float cy = DOCenter((long int)(3 * D3D1ConvPosition + 1));
                float cz = DOCenter((long int)(3 * D3D1ConvPosition + 2));

                // Calculate critical diagonal lengths for the new active cell located at (xp,yp,zp) on the local grid
                // For each neighbor (l=0 to 25), calculate which octahedron face leads to cell capture
                // Calculate critical octahedron diagonal length to activate each nearest neighbor, as well as the
                // coordinates of the triangle vertices on the capturing face
                for (int n = 0; n < 26; n++) {

                    // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the image of the center
                    // of a neighbor cell
                    double x0 = xp + NeighborX(n) - cx;
                    double y0 = yp + NeighborY(n) - cy;
                    double z0 = zp + NeighborZ(n) - cz;
                    // mag0 is the magnitude of (x0,y0,z0)
                    double mag0 = pow(pow(x0, 2.0) + pow(y0, 2.0) + pow(z0, 2.0), 0.5);

                    // Calculate unit vectors for the octahedron that intersect the new cell center
                    double Diag1X, Diag1Y, Diag1Z, Diag2X, Diag2Y, Diag2Z, Diag3X, Diag3Y, Diag3Z;
                    double Angle1 =
                        (GrainUnitVector(9 * MyOrientation) * x0 + GrainUnitVector(9 * MyOrientation + 1) * y0 +
                         GrainUnitVector(9 * MyOrientation + 2) * z0) /
                        mag0;
                    if (Angle1 < 0) {
                        Diag1X = GrainUnitVector(9 * MyOrientation);
                        Diag1Y = GrainUnitVector(9 * MyOrientation + 1);
                        Diag1Z = GrainUnitVector(9 * MyOrientation + 2);
                    }
                    else {
                        Diag1X = -GrainUnitVector(9 * MyOrientation);
                        Diag1Y = -GrainUnitVector(9 * MyOrientation + 1);
                        Diag1Z = -GrainUnitVector(9 * MyOrientation + 2);
                    }
                    double Angle2 =
                        (GrainUnitVector(9 * MyOrientation + 3) * x0 + GrainUnitVector(9 * MyOrientation + 4) * y0 +
                         GrainUnitVector(9 * MyOrientation + 5) * z0) /
                        mag0;
                    if (Angle2 < 0) {
                        Diag2X = GrainUnitVector(9 * MyOrientation + 3);
                        Diag2Y = GrainUnitVector(9 * MyOrientation + 4);
                        Diag2Z = GrainUnitVector(9 * MyOrientation + 5);
                    }
                    else {
                        Diag2X = -GrainUnitVector(9 * MyOrientation + 3);
                        Diag2Y = -GrainUnitVector(9 * MyOrientation + 4);
                        Diag2Z = -GrainUnitVector(9 * MyOrientation + 5);
                    }

                    double Angle3 =
                        (GrainUnitVector(9 * MyOrientation + 6) * x0 + GrainUnitVector(9 * MyOrientation + 7) * y0 +
                         GrainUnitVector(9 * MyOrientation + 8) * z0) /
                        mag0;
                    if (Angle3 < 0) {
                        Diag3X = GrainUnitVector(9 * MyOrientation + 6);
                        Diag3Y = GrainUnitVector(9 * MyOrientation + 7);
                        Diag3Z = GrainUnitVector(9 * MyOrientation + 8);
                    }
                    else {
                        Diag3X = -GrainUnitVector(9 * MyOrientation + 6);
                        Diag3Y = -GrainUnitVector(9 * MyOrientation + 7);
                        Diag3Z = -GrainUnitVector(9 * MyOrientation + 8);
                    }

                    double U1[3], U2[3], UU[3], Norm[3];
                    U1[0] = Diag2X - Diag1X;
                    U1[1] = Diag2Y - Diag1Y;
                    U1[2] = Diag2Z - Diag1Z;
                    U2[0] = Diag3X - Diag1X;
                    U2[1] = Diag3Y - Diag1Y;
                    U2[2] = Diag3Z - Diag1Z;
                    UU[0] = U1[1] * U2[2] - U1[2] * U2[1];
                    UU[1] = U1[2] * U2[0] - U1[0] * U2[2];
                    UU[2] = U1[0] * U2[1] - U1[1] * U2[0];
                    double NDem = sqrt(UU[0] * UU[0] + UU[1] * UU[1] + UU[2] * UU[2]);
                    Norm[0] = UU[0] / NDem;
                    Norm[1] = UU[1] / NDem;
                    Norm[2] = UU[2] / NDem;
                    // normal to capturing plane
                    double normx = Norm[0];
                    double normy = Norm[1];
                    double normz = Norm[2];
                    double ParaT =
                        (normx * x0 + normy * y0 + normz * z0) / (normx * Diag1X + normy * Diag1Y + normz * Diag1Z);
                    float CDLVal =
                        pow(pow(ParaT * Diag1X, 2.0) + pow(ParaT * Diag1Y, 2.0) + pow(ParaT * Diag1Z, 2.0), 0.5);
                    CritDiagonalLength((long int)(26) * D3D1ConvPosition + (long int)(n)) = CDLVal;
                }
            }
        });
}
