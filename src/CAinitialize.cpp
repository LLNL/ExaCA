// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAinitialize.hpp"

#include "CAconfig.hpp"
#include "CAfunctions.hpp"
#include "CAghostnodes.hpp"
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
                       int &TempFilesInSeries, std::vector<std::string> &temp_paths, double &HT_deltax,
                       bool &RemeltingYN, double &deltat, int &NumberOfLayers, int &LayerHeight,
                       std::string &SubstrateFileName, float &SubstrateGrainSpacing, bool &UseSubstrateFile, double &G,
                       double &R, int &nx, int &ny, int &nz, double &FractSurfaceSitesActive, std::string &PathToOutput,
                       int &PrintDebug, bool &PrintMisorientation, bool &PrintFinalUndercoolingVals,
                       bool &PrintFullOutput, int &NSpotsX, int &NSpotsY, int &SpotOffset, int &SpotRadius,
                       bool &PrintTimeSeries, int &TimeSeriesInc, bool &PrintIdleTimeSeriesFrames,
                       bool &PrintDefaultRVE, double &RNGSeed, double &FractPowderSitesActive) {

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
        "Debug check (reduced)",                        // Optional input 0
        "Debug check (extensive)",                      // Optional input 1
        "Print intermediate output frames",             // Optional input 2
        "separate frames",                              // Optional input 3
        "output even if system is unchanged",           // Optional input 4
        "file of final undercooling values",            // Optional input 5
        "Random seed for grains and nuclei generation", // Optional input 6
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
    // An "M" after the problem type ("SM" or "RM") indicates that the problem uses remelting logic
    // Additional required/optional inputs depending on problem type
    SimulationType = parseInput(InputData, "Problem type");
    if (SimulationType == "RM") {
        // Simulation using external temperature data ("R") with remelting
        SimulationType = "R";
        RemeltingYN = true;
    }
    else if (SimulationType == "SM") {
        // Simulation using hemispherical spot melt data ("S") with remelting
        SimulationType = "S";
        RemeltingYN = true;
    }
    else {
        // Simulation does not including remelting logic
        RemeltingYN = false;
    }
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
        OptionalInputs_ProblemSpecific.resize(3);
        OptionalInputs_ProblemSpecific[0] = "Substrate grain spacing";
        OptionalInputs_ProblemSpecific[1] = "Substrate filename";
        OptionalInputs_ProblemSpecific[2] = "Fraction of powder sites active";
    }
    else if (SimulationType == "R") {
        RequiredInputs_ProblemSpecific.resize(5);
        RequiredInputs_ProblemSpecific[0] = "Time step";
        RequiredInputs_ProblemSpecific[1] = "Temperature filename";
        RequiredInputs_ProblemSpecific[2] = "Number of temperature files";
        RequiredInputs_ProblemSpecific[3] = "Number of layers";
        RequiredInputs_ProblemSpecific[4] = "Offset between layers";
        OptionalInputs_ProblemSpecific.resize(7);
        OptionalInputs_ProblemSpecific[0] = "Substrate grain spacing";
        OptionalInputs_ProblemSpecific[1] = "Substrate filename";
        OptionalInputs_ProblemSpecific[2] = "Heat transport data mesh size";
        OptionalInputs_ProblemSpecific[3] = "Path to temperature file(s)";
        OptionalInputs_ProblemSpecific[4] = "Extra set of wall cells";
        OptionalInputs_ProblemSpecific[5] = "default RVE output";
        OptionalInputs_ProblemSpecific[6] = "Fraction of powder sites active";
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
            if (RemeltingYN)
                std::cout << "This simulation includes logic for cells melting and multiple solidification events"
                          << std::endl;
        }
    }
    else if (SimulationType == "C") {
        G = getInputDouble(RequiredInputsRead_ProblemSpecific[0]);
        R = getInputDouble(RequiredInputsRead_ProblemSpecific[1]);
        NRatio = getInputInt(RequiredInputsRead_ProblemSpecific[2]);
        nx = getInputInt(RequiredInputsRead_ProblemSpecific[3]);
        ny = getInputInt(RequiredInputsRead_ProblemSpecific[4]);
        nz = getInputInt(RequiredInputsRead_ProblemSpecific[5]);
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
        nz = SpotRadius + 1 + (NumberOfLayers - 1) * LayerHeight;
        nx = 2 * SpotRadius + 1 + SpotOffset * (NSpotsX - 1);
        ny = 2 * SpotRadius + 1 + SpotOffset * (NSpotsY - 1);
        // Fraction of powder sites active, if given, should be set. Otherwise, the default value is 1.0
        if (OptionalInputsRead_ProblemSpecific[2].empty())
            FractPowderSitesActive = 1.0;
        else
            FractPowderSitesActive = getInputDouble(OptionalInputsRead_ProblemSpecific[2]);
        if (id == 0) {
            std::cout << "CA Simulation using a radial, fixed thermal gradient of " << G
                      << " K/m as a series of hemispherical spots, and a cooling rate of " << R << " K/s" << std::endl;
            std::cout << "A total of " << NumberOfLayers << " spots per layer, with layers offset by " << LayerHeight
                      << " CA cells will be simulated" << std::endl;
            std::cout << "The time step is " << deltat * pow(10, 6) << " microseconds" << std::endl;
            if (RemeltingYN)
                std::cout << "This simulation includes logic for cells melting and multiple solidification events"
                          << std::endl;
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
    // Should the final undercooling of all cells be printed?
    if (OptionalInputsRead_General[5].empty())
        PrintFinalUndercoolingVals = false;
    else
        PrintFinalUndercoolingVals = getInputBool(OptionalInputsRead_General[5]);
    // Should the RNG seed for grain/nuclei initialization use the default, or a custom value?
    if (OptionalInputsRead_General[6].empty())
        RNGSeed = 0.0;
    else
        RNGSeed = getInputDouble(OptionalInputsRead_General[6]);
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

        if ((RemeltingYN) && (HT_deltax != deltax)) {
            throw std::runtime_error("Error: For simulations with external temperature data and remelting logic, CA "
                                     "cell size and input temperature data resolution must be equivalent");
        }
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
        if ((!(OptionalInputsRead_ProblemSpecific[4].empty())) &&
            (id == 0)) // Fixme: Remove this optional input eventually
            std::cout << "Note: optional input ExtraWalls is no longer used, all simulations by default have no walls "
                         "cells along domain boundaries"
                      << std::endl;
        // Should optional RVE data be printed for the standard location (center of domain in X and Y, as close to the
        // top of the domain in Z as possiblewithout including the last layer's microstructure)?
        PrintDefaultRVE = false;
        if (!(OptionalInputsRead_ProblemSpecific[5].empty()))
            PrintDefaultRVE = getInputBool(OptionalInputsRead_ProblemSpecific[5]);

        // Fraction of powder sites active, if given, should be set. Otherwise, the default value is 1.0
        if (OptionalInputsRead_ProblemSpecific[6].empty())
            FractPowderSitesActive = 1.0;
        else
            FractPowderSitesActive = getInputDouble(OptionalInputsRead_ProblemSpecific[6]);
    }
    else {
        // RVE data print option is only for simulation type R
        PrintDefaultRVE = false;
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
// Intialize neighbor list structures (NeighborX, NeighborY, NeighborZ)
void NeighborListInit(ViewI &NeighborX, ViewI &NeighborY, ViewI &NeighborZ) {

    // Temporary host views for initialization
    ViewI_H NeighborX_Host(Kokkos::ViewAllocateWithoutInitializing("NeighborX_Host"), 26);
    ViewI_H NeighborY_Host(Kokkos::ViewAllocateWithoutInitializing("NeighborY_Host"), 26);
    ViewI_H NeighborZ_Host(Kokkos::ViewAllocateWithoutInitializing("NeighborZ_Host"), 26);

    // Assignment of neighbors around a cell "X" is as follows (in order of closest to furthest from cell "X")
    // Neighbors 0 through 8 are in the -Y direction
    // Neighbors 9 through 16 are in the XY plane with cell X
    // Neighbors 17 through 25 are in the +Y direction

    NeighborX_Host(0) = 0;
    NeighborY_Host(0) = -1;
    NeighborZ_Host(0) = 0;
    NeighborX_Host(1) = 1;
    NeighborY_Host(1) = -1;
    NeighborZ_Host(1) = 0;
    NeighborX_Host(2) = -1;
    NeighborY_Host(2) = -1;
    NeighborZ_Host(2) = 0;
    NeighborX_Host(3) = 0;
    NeighborY_Host(3) = -1;
    NeighborZ_Host(3) = 1;
    NeighborX_Host(4) = 0;
    NeighborY_Host(4) = -1;
    NeighborZ_Host(4) = -1;
    NeighborX_Host(5) = -1;
    NeighborY_Host(5) = -1;
    NeighborZ_Host(5) = 1;
    NeighborX_Host(6) = 1;
    NeighborY_Host(6) = -1;
    NeighborZ_Host(6) = 1;
    NeighborX_Host(7) = -1;
    NeighborY_Host(7) = -1;
    NeighborZ_Host(7) = -1;
    NeighborX_Host(8) = 1;
    NeighborY_Host(8) = -1;
    NeighborZ_Host(8) = -1;

    NeighborX_Host(9) = 0;
    NeighborY_Host(9) = 0;
    NeighborZ_Host(9) = 1;
    NeighborX_Host(10) = 0;
    NeighborY_Host(10) = 0;
    NeighborZ_Host(10) = -1;
    NeighborX_Host(11) = 1;
    NeighborY_Host(11) = 0;
    NeighborZ_Host(11) = 1;
    NeighborX_Host(12) = -1;
    NeighborY_Host(12) = 0;
    NeighborZ_Host(12) = 1;
    NeighborX_Host(13) = 1;
    NeighborY_Host(13) = 0;
    NeighborZ_Host(13) = -1;
    NeighborX_Host(14) = -1;
    NeighborY_Host(14) = 0;
    NeighborZ_Host(14) = -1;
    NeighborX_Host(15) = 1;
    NeighborY_Host(15) = 0;
    NeighborZ_Host(15) = 0;
    NeighborX_Host(16) = -1;
    NeighborY_Host(16) = 0;
    NeighborZ_Host(16) = 0;

    NeighborX_Host(17) = 0;
    NeighborY_Host(17) = 1;
    NeighborZ_Host(17) = 0;
    NeighborX_Host(18) = 1;
    NeighborY_Host(18) = 1;
    NeighborZ_Host(18) = 0;
    NeighborX_Host(19) = -1;
    NeighborY_Host(19) = 1;
    NeighborZ_Host(19) = 0;
    NeighborX_Host(20) = 0;
    NeighborY_Host(20) = 1;
    NeighborZ_Host(20) = 1;
    NeighborX_Host(21) = 0;
    NeighborY_Host(21) = 1;
    NeighborZ_Host(21) = -1;
    NeighborX_Host(22) = 1;
    NeighborY_Host(22) = 1;
    NeighborZ_Host(22) = 1;
    NeighborX_Host(23) = -1;
    NeighborY_Host(23) = 1;
    NeighborZ_Host(23) = 1;
    NeighborX_Host(24) = 1;
    NeighborY_Host(24) = 1;
    NeighborZ_Host(24) = -1;
    NeighborX_Host(25) = -1;
    NeighborY_Host(25) = 1;
    NeighborZ_Host(25) = -1;

    // Copy data back to device views
    NeighborX = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), NeighborX_Host);
    NeighborY = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), NeighborY_Host);
    NeighborZ = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), NeighborZ_Host);
}

// Obtain the physical XYZ bounds of the domain, using either domain size from the input file, or reading temperature
// data files and parsing the coordinates
void FindXYZBounds(std::string SimulationType, int id, double &deltax, int &nx, int &ny, int &nz,
                   std::vector<std::string> &temp_paths, float &XMin, float &XMax, float &YMin, float &YMax,
                   float &ZMin, float &ZMax, int &LayerHeight, int NumberOfLayers, int TempFilesInSeries,
                   float *ZMinLayer, float *ZMaxLayer, int SpotRadius) {

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
        // CA cells in each direction span from the lower to the higher bound of the temperature data - without wall
        // cells or padding around the simulation edges
        nx = round((XMax - XMin) / deltax) + 1;
        ny = round((YMax - YMin) / deltax) + 1;
        nz = round((ZMax - ZMin) / deltax) + 1;
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
        ZMax = nz * deltax;
        // If this is a spot melt problem, also set the ZMin/ZMax for each layer
        for (int n = 0; n < NumberOfLayers; n++) {
            ZMinLayer[n] = deltax * (LayerHeight * n);
            ZMaxLayer[n] = deltax * (SpotRadius + LayerHeight * n);
        }
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
                         int &ProcessorsInXDirection, int &ProcessorsInYDirection, long int &LocalDomainSize,
                         bool &AtNorthBoundary, bool &AtSouthBoundary, bool &AtEastBoundary, bool &AtWestBoundary) {

    // Determine which subdomains are at which locations on the grid relative to the others
    InitialDecomposition(DecompositionStrategy, nx, ny, ProcessorsInXDirection, ProcessorsInYDirection, id, np,
                         NeighborRank_North, NeighborRank_South, NeighborRank_East, NeighborRank_West,
                         NeighborRank_NorthEast, NeighborRank_NorthWest, NeighborRank_SouthEast, NeighborRank_SouthWest,
                         AtNorthBoundary, AtSouthBoundary, AtEastBoundary, AtWestBoundary);
    // Determine, for each MPI process id, the local grid size in x and y (and the offsets in x and y relative to the
    // overall simulation domain)
    MyXOffset = XOffsetCalc(id, nx, ProcessorsInXDirection, ProcessorsInYDirection, DecompositionStrategy);
    MyXSlices = XMPSlicesCalc(id, nx, ProcessorsInXDirection, ProcessorsInYDirection, DecompositionStrategy);

    MyYOffset = YOffsetCalc(id, ny, ProcessorsInYDirection, np, DecompositionStrategy);
    MyYSlices = YMPSlicesCalc(id, ny, ProcessorsInYDirection, np, DecompositionStrategy);

    // Add ghost nodes at subdomain overlaps
    AddGhostNodes(DecompositionStrategy, NeighborRank_West, NeighborRank_East, NeighborRank_North, NeighborRank_South,
                  MyXSlices, MyXOffset, MyYSlices, MyYOffset);

    LocalDomainSize = MyXSlices * MyYSlices * nz; // Number of cells on this MPI rank
}

// Read in temperature data from files, stored in "RawData", with the appropriate MPI ranks storing the appropriate data
void ReadTemperatureData(int id, double &deltax, double HT_deltax, int &HTtoCAratio, int MyXSlices, int MyYSlices,
                         int MyXOffset, int MyYOffset, float XMin, float YMin, std::vector<std::string> &temp_paths,
                         int NumberOfLayers, int TempFilesInSeries, unsigned int &NumberOfTemperatureDataPoints,
                         std::vector<double> &RawData, int *FirstValue, int *LastValue, bool RemeltingYN,
                         ViewI &MaxSolidificationEvents, float *ZMinLayer, float *ZMaxLayer, int LayerHeight) {

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
    HTtoCAratio = round(HT_deltax / deltax); // OpenFOAM/CA cell size ratio
    // If HTtoCAratio > 1, an interpolation of input temperature data is needed
    // The X and Y bounds are the region (for this MPI rank) of the physical domain that needs to be
    // read extends past the actual spatial extent of the local domain for purposes of interpolating
    // from HT_deltax to deltax
    int LowerXBound = MyXOffset - (MyXOffset % HTtoCAratio);
    int LowerYBound = MyYOffset - (MyYOffset % HTtoCAratio);
    int UpperXBound, UpperYBound;
    if (HTtoCAratio == 1) {
        UpperXBound = MyXOffset + MyXSlices - 1;
        UpperYBound = MyYOffset + MyYSlices - 1;
    }
    else {
        UpperXBound = MyXOffset + MyXSlices - 1 + HTtoCAratio - (MyXOffset + MyXSlices - 1) % HTtoCAratio;
        UpperYBound = MyYOffset + MyYSlices - 1 + HTtoCAratio - (MyYOffset + MyYSlices - 1) % HTtoCAratio;
    }

    // Maximum number of times a cell in a given layer will undergo solidification
    ViewI_H MaxSolidificationEvents_Host(Kokkos::ViewAllocateWithoutInitializing("MaxSolidificationEvents_Host"),
                                         NumberOfLayers);

    // Store raw data relevant to each rank in the vector structure RawData
    // Two passes through reading temperature data files- this is the second pass, reading the actual X/Y/Z/liquidus
    // time/cooling rate data and each rank stores the data relevant to itself in "RawData". With remelting
    // (SimulationType == "RM"), this is the same except that some X/Y/Z coordinates may be repeated in a file, and
    // a "melting time" value is stored in addition to liquidus time and cooling rate
    // Second pass through the files - ignore header line
    int LayersToRead = std::min(NumberOfLayers, TempFilesInSeries); // was given in input file
    for (int LayerReadCount = 1; LayerReadCount <= LayersToRead; LayerReadCount++) {

        // View to store, on each rank for this file, the max number of times a cell goes above/below the liquidus
        int TempSizeX = 0;
        int TempSizeY = 0;
        int TempSizeZ = 0;
        if (RemeltingYN) {
            TempSizeX = UpperXBound - LowerXBound + 1;
            TempSizeY = UpperYBound - LowerYBound + 1;
            TempSizeZ = round((ZMaxLayer[LayerReadCount - 1] - ZMinLayer[LayerReadCount - 1]) / deltax) + 1;
        }
        // Init to 0
        ViewI3D_H TempMeltCount("TempMeltCount", TempSizeZ, TempSizeX, TempSizeY);

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
            XInt = round((XYZTemperaturePoint[0] - XMin) / deltax);
            YInt = round((XYZTemperaturePoint[1] - YMin) / deltax);
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
                // Increment view for number of times this cell coordinate has appeared in the file
                // Since ZMinLayer is not the raw "smallest Z value" from the temperature file, but rather is the
                // "CA-adjusted" smallest Z value with each layer being offset by deltax * LayerHeight meters, this same
                // adjustment must be made to the raw Z value read from the file here
                if (RemeltingYN) {
                    int ZInt = round((XYZTemperaturePoint[2] + deltax * LayerHeight * (LayerReadCount - 1) -
                                      ZMinLayer[LayerReadCount - 1]) /
                                     deltax);
                    TempMeltCount(ZInt, XInt - LowerXBound, YInt - LowerYBound)++;
                }
                if (NumberOfTemperatureDataPoints >= RawData.size() - 6) {
                    int OldSize = RawData.size();
                    RawData.resize(OldSize + 1000000);
                }
            }
        }
        LastValue[LayerReadCount - 1] = NumberOfTemperatureDataPoints;
        // Determine max number of remelting events for each layer (may be different on each MPI rank)
        // If no remelting, this value is 1 for each layer (all cells solidify once, since they can't remelt)
        if (RemeltingYN) {
            int MaxCount = 0;
            for (int k = 0; k < TempSizeZ; k++) {
                for (int i = 0; i < TempSizeX; i++) {
                    for (int j = 0; j < TempSizeY; j++) {
                        if (TempMeltCount(k, i, j) > MaxCount)
                            MaxCount = TempMeltCount(k, i, j);
                    }
                }
            }
            int MaxCountGlobal;
            MPI_Allreduce(&MaxCount, &MaxCountGlobal, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            MaxSolidificationEvents_Host(LayerReadCount - 1) = MaxCountGlobal;
        }
        else {
            MaxSolidificationEvents_Host(LayerReadCount - 1) = 1;
        }
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
                MaxSolidificationEvents_Host(LayerReadCount) = MaxSolidificationEvents_Host(LayerReadCount - 1);
            }
            else {
                // All layers have different temperature data but in a repeating pattern
                int RepeatedFile = (LayerReadCount) % TempFilesInSeries;
                FirstValue[LayerReadCount] = FirstValue[RepeatedFile];
                LastValue[LayerReadCount] = LastValue[RepeatedFile];
                MaxSolidificationEvents_Host(LayerReadCount) = MaxSolidificationEvents_Host(LayerReadCount - 1);
            }
        }
    }
    // Copy MaxSolidificationEvents back to device
    MaxSolidificationEvents =
        Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), MaxSolidificationEvents_Host);
}

//*****************************************************************************/
// Get the Z coordinate of the lower bound of iteration
// If using remelting, this is called before initializing temperature data views for each layer
// If not using remelting, this is only called before initializing temperature data view for the first layer, since the
// existing DomainShiftAndResize subroutine already exists for other layers. These functions should be unified in a
// future update
int calcZBound_Low(bool RemeltingYN, std::string SimulationType, int LayerHeight, int layernumber, float *ZMinLayer,
                   float ZMin, double deltax) {

    int ZBound_Low = 0; // set to zero if no remelting (was previously done in individual TempInit_[ProblemType]
                        // subroutines - implicitly assumes bottom of layer 0 is the bottom of the overall domain - this
                        // should be fixed in the future for edge cases where this isn't true)
    if (RemeltingYN) {
        if (SimulationType == "S") {
            // lower bound of domain is an integer multiple of the layer spacing, since the temperature field is the
            // same for every layer
            ZBound_Low = LayerHeight * layernumber;
        }
        else if (SimulationType == "R") {
            // lower bound of domain is based on the data read from the file(s)
            ZBound_Low = round((ZMinLayer[layernumber] - ZMin) / deltax);
        }
        else
            throw std::runtime_error("Error: simulations with remelting must be simulation type SM or RM");
    }
    return ZBound_Low;
}
//*****************************************************************************/
// Get the Z coordinate of the upper bound of iteration
// If using remelting, this is called before initializing temperature data views for each layer
// If not using remelting, this is only called before initializing temperature data view for the first layer, since the
// existing DomainShiftAndResize subroutine already exists for other layers. These functions should be unified in a
// future update
int calcZBound_High(std::string SimulationType, int SpotRadius, int LayerHeight, int layernumber, float ZMin,
                    double deltax, int nz, float *ZMaxLayer) {

    int ZBound_High = -1; // assign dummy initial value
    if (SimulationType == "C") {
        // Not a multilayer problem, top of "layer" is the top of the overall simulation domain
        ZBound_High = nz - 1;
    }
    else if (SimulationType == "S") {
        // Top of layer is equal to the spot radius for a problem of hemispherical spot solidification, plus an offset
        // depending on the layer number
        ZBound_High = SpotRadius + LayerHeight * layernumber;
    }
    else if (SimulationType == "R") {
        // Top of layer comes from the layer's file data (implicitly assumes bottom of layer 0 is the bottom of the
        // overall domain - this should be fixed in the future for edge cases where this isn't true)
        ZBound_High = round((ZMaxLayer[layernumber] - ZMin) / deltax);
    }
    if (ZBound_High == -1)
        throw std::runtime_error(
            "Error: ZBound_High went uninitialized, problem type may not have been a valid option");
    return ZBound_High;
}
//*****************************************************************************/
// Calculate the size of the active domain in Z
int calcnzActive(int ZBound_Low, int ZBound_High, int id, int layernumber) {
    int nzActive = ZBound_High - ZBound_Low + 1;
    if (id == 0)
        std::cout << "Layer " << layernumber << "'s active domain is from Z = " << ZBound_Low << " through "
                  << ZBound_High << " (" << nzActive << ") cells" << std::endl;
    return nzActive;
}
//*****************************************************************************/
// Calculate the size of the domain, as a number of cells
int calcLocalActiveDomainSize(int MyXSlices, int MyYSlices, int nzActive) {
    int LocalActiveDomainSize = MyXSlices * MyYSlices * nzActive;
    return LocalActiveDomainSize;
}
//*****************************************************************************/
// Initialize temperature data for a constrained solidification test problem
void TempInit_DirSolidification(double G, double R, int, int &MyXSlices, int &MyYSlices, double deltax, double deltat,
                                int nz, int LocalDomainSize, ViewI &CritTimeStep, ViewF &UndercoolingChange,
                                bool *Melted, ViewI &LayerID, ViewI &MaxSolidificationEvents) {

    // These views are initialized on the host, filled with data, and then copied to the device for layer "layernumber"
    // This view is initialized with zeros
    ViewI_H LayerID_Host("LayerID_H", LocalDomainSize);
    // These views will be filled with non-zero values
    ViewI_H CritTimeStep_Host(Kokkos::ViewAllocateWithoutInitializing("CritTimeStep_H"), LocalDomainSize);
    ViewF_H UndercoolingChange_Host(Kokkos::ViewAllocateWithoutInitializing("UndercoolingChange_H"), LocalDomainSize);
    ViewI_H MaxSolidificationEvents_Host(Kokkos::ViewAllocateWithoutInitializing("MaxSolidificationEvents_H"), 1);

    // Initialize temperature field in Z direction with thermal gradient G set in input file
    // Cells at the bottom surface (Z = 0) are at the liquidus at time step 0 (no wall cells at the bottom boundary)
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < MyXSlices; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                int GlobalD3D1ConvPosition = k * MyXSlices * MyYSlices + i * MyYSlices + j;
                UndercoolingChange_Host(GlobalD3D1ConvPosition) = R * deltat;
                CritTimeStep_Host(GlobalD3D1ConvPosition) = (int)((k * G * deltax) / (R * deltat));
                Melted[GlobalD3D1ConvPosition] = true;
            }
        }
    }
    // No remelting, only one layer: max number of times a cell will become solid is 1
    MaxSolidificationEvents_Host(0) = 1;

    // Copy initialized host data back to device
    CritTimeStep = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), CritTimeStep_Host);
    LayerID = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), LayerID_Host);
    UndercoolingChange = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), UndercoolingChange_Host);
    MaxSolidificationEvents =
        Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), MaxSolidificationEvents_Host);
}

// Initialize temperature data for an array of overlapping spot melts (done during simulation initialization, no
// remelting)
void TempInit_SpotMelt(double G, double R, std::string, int id, int &MyXSlices, int &MyYSlices, int &MyXOffset,
                       int &MyYOffset, double deltax, double deltat, int nz, int LocalDomainSize, ViewI &CritTimeStep,
                       ViewF &UndercoolingChange, bool *Melted, int LayerHeight, int NumberOfLayers,
                       double FreezingRange, ViewI &LayerID, int NSpotsX, int NSpotsY, int SpotRadius, int SpotOffset,
                       ViewI &MaxSolidificationEvents) {

    // This view will be filled with non-zero values on the host, and later copied to the device
    ViewI_H LayerID_Host(Kokkos::ViewAllocateWithoutInitializing("LayerID_H"), LocalDomainSize);
    // These views are initialized with zero values on the host, populated with other data, and later copied to the
    // device
    ViewI_H CritTimeStep_Host("CritTimeStep_H", LocalDomainSize);
    ViewF_H UndercoolingChange_Host("UndercoolingChange_H", LocalDomainSize);
    ViewI_H MaxSolidificationEvents_Host("MaxSolidificationEvents_H", NumberOfLayers);

    // No cells intitially have undercooling, nor have melting/solidification data
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < MyXSlices; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                int GlobalD3D1ConvPosition = k * MyXSlices * MyYSlices + i * MyYSlices + j;
                Melted[GlobalD3D1ConvPosition] = false;
                LayerID_Host(GlobalD3D1ConvPosition) = -1;
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
        // Without remelting, max number of times a cell will become solid is 1, regardless of the layer number
        MaxSolidificationEvents_Host(layernumber) = 1;
        for (int n = 0; n < NumberOfSpots; n++) {
            if (id == 0)
                std::cout << "Initializing spot " << n << " on layer " << layernumber << std::endl;
            // Initialize critical time step/cooling rate values for this spot/this layer
            int XSpotPos = SpotRadius + (n % NSpotsX) * SpotOffset;
            int YSpotPos = SpotRadius + (n / NSpotsX) * SpotOffset;
            int ZSpotPos = SpotRadius + LayerHeight * layernumber;
            for (int k = 0; k <= SpotRadius; k++) {
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
                            CritTimeStep_Host(GlobalD3D1ConvPosition) =
                                1 + (int)(((float)(SpotRadius)-TotDist) / IsothermVelocity) + TimeBetweenSpots * n;
                            UndercoolingChange_Host(GlobalD3D1ConvPosition) = R * deltat;
                            LayerID_Host(GlobalD3D1ConvPosition) = layernumber;
                            Melted[GlobalD3D1ConvPosition] = true;
                        }
                    }
                }
            }
        }
    }

    // Copy initialized host data back to device
    CritTimeStep = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), CritTimeStep_Host);
    LayerID = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), LayerID_Host);
    UndercoolingChange = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), UndercoolingChange_Host);
    MaxSolidificationEvents =
        Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), MaxSolidificationEvents_Host);
}

// For an overlapping spot melt pattern, determine the maximum number of times a cell will melt/solidify as part of a
// layer
int calcMaxSolidificationEventsSpot(int MyXSlices, int MyYSlices, int NumberOfSpots, int NSpotsX, int SpotRadius,
                                    int SpotOffset, int MyXOffset, int MyYOffset) {

    ViewI2D_H MaxSolidificationEvents_Temp("SEvents_Temp", MyXSlices, MyYSlices);
    for (int n = 0; n < NumberOfSpots; n++) {
        int XSpotPos = SpotRadius + (n % NSpotsX) * SpotOffset;
        int YSpotPos = SpotRadius + (n / NSpotsX) * SpotOffset;
        for (int i = 0; i < MyXSlices; i++) {
            int XGlobal = i + MyXOffset;
            float DistX = (float)(XSpotPos - XGlobal);
            for (int j = 0; j < MyYSlices; j++) {
                int YGlobal = j + MyYOffset;
                float DistY = (float)(YSpotPos - YGlobal);
                float TotDist = sqrt(DistX * DistX + DistY * DistY);
                if (TotDist <= SpotRadius) {
                    MaxSolidificationEvents_Temp(i, j)++;
                }
            }
        }
    }
    int TempMax = 0;
    for (int i = 0; i < MyXSlices; i++) {
        for (int j = 0; j < MyYSlices; j++) {
            if (MaxSolidificationEvents_Temp(i, j) > TempMax) {
                TempMax = MaxSolidificationEvents_Temp(i, j);
            }
        }
    }
    // Max solidification events should be the same on each rank (and each layer, since the pattern for all layers is
    // identical)
    int GlobalMaxSEvents;
    MPI_Allreduce(&TempMax, &GlobalMaxSEvents, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    return GlobalMaxSEvents;
}

// Initialize temperature data for an array of overlapping spot melts (done at the start of each layer, with remelting)
void TempInit_SpotMeltRemelting(int layernumber, double G, double R, std::string, int id, int &MyXSlices,
                                int &MyYSlices, int &MyXOffset, int &MyYOffset, double deltax, double deltat,
                                int ZBound_Low, int nz, int LocalActiveDomainSize, int LocalDomainSize,
                                ViewI &CritTimeStep, ViewF &UndercoolingChange, ViewF &UndercoolingCurrent,
                                bool *Melted, int, double FreezingRange, ViewI &LayerID, int NSpotsX, int NSpotsY,
                                int SpotRadius, int SpotOffset, ViewF3D &LayerTimeTempHistory,
                                ViewI &NumberOfSolidificationEvents, ViewI &MeltTimeStep,
                                ViewI &MaxSolidificationEvents, ViewI &SolidificationEventCounter) {

    int NumberOfSpots = NSpotsX * NSpotsY;

    // Temporary host view for the maximum number of times a cell in a given layer will solidify
    ViewI_H MaxSolidificationEvents_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), MaxSolidificationEvents);
    MaxSolidificationEvents_Host(layernumber) = calcMaxSolidificationEventsSpot(
        MyXSlices, MyYSlices, NumberOfSpots, NSpotsX, SpotRadius, SpotOffset, MyXOffset, MyYOffset);

    // These views are initialized to zeros on the host (requires knowing MaxSolidificationEvents first), filled with
    // data, and then copied to the device for layer "layernumber"
    ViewF3D_H LayerTimeTempHistory_Host("TimeTempHistory_H", LocalActiveDomainSize,
                                        MaxSolidificationEvents_Host(layernumber), 3);
    ViewI_H NumberOfSolidificationEvents_Host("NumSEvents_H", LocalActiveDomainSize);

    // Resize device views for active domain size if initializing first layer (don't resize after that, as all layers
    // are the same)
    if (layernumber == 0) {
        Kokkos::resize(LayerTimeTempHistory, LocalActiveDomainSize, MaxSolidificationEvents_Host(0), 3);
        Kokkos::resize(NumberOfSolidificationEvents, LocalActiveDomainSize);
        Kokkos::resize(SolidificationEventCounter, LocalActiveDomainSize);
        Kokkos::resize(MeltTimeStep, LocalDomainSize);
    }

    // Temporary host views for storing initialized temperature data for active region data structures
    // No resize of device views necessary, as multilayer spot melt simulations have the same active domain size for all
    // layers These views are copied to the host, updated for layer "layernumber", and later copied back to the device
    ViewI_H LayerID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), LayerID);
    ViewI_H MeltTimeStep_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), MeltTimeStep);
    ViewI_H CritTimeStep_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CritTimeStep);
    ViewF_H UndercoolingChange_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), UndercoolingChange);

    if (layernumber == 0) {
        // No cells intitially have undercooling, nor have melting/solidification data
        for (int k = 0; k < nz; k++) {
            for (int i = 0; i < MyXSlices; i++) {
                for (int j = 0; j < MyYSlices; j++) {
                    int GlobalD3D1ConvPosition = k * MyXSlices * MyYSlices + i * MyYSlices + j;
                    Melted[GlobalD3D1ConvPosition] = false;
                    LayerID_Host(GlobalD3D1ConvPosition) = -1;
                }
            }
        }
    }

    // Outer edges of spots are initialized at the liquidus temperature
    // Spots cool at constant rate R, spot thermal gradient = G
    // Time between "start" of next spot is the time it takes for the previous spot
    // to have entirely gone below the solidus temperature
    float IsothermVelocity = (R / G) * deltat / deltax;                                  // in cells per time step
    int TimeBetweenSpots = SpotRadius / IsothermVelocity + (FreezingRange / R) / deltat; // in time steps

    if (id == 0)
        std::cout << "Initializing temperature field for " << NumberOfSpots << " spots on layer " << layernumber
                  << ", each of which takes approximately " << TimeBetweenSpots << " time steps to solidify"
                  << std::endl;

    int Count = 0;
    for (int n = 0; n < NumberOfSpots; n++) {
        if (id == 0)
            std::cout << "Initializing spot " << n << " on layer " << layernumber << std::endl;
        // Initialize LayerTimeTempHistory data values for this spot/this layer - relative to the layer bottom
        int XSpotPos = SpotRadius + (n % NSpotsX) * SpotOffset;
        int YSpotPos = SpotRadius + (n / NSpotsX) * SpotOffset;
        for (int k = 0; k <= SpotRadius; k++) {
            // Distance of this cell from the spot center
            float DistZ = (float)(SpotRadius - k);
            for (int i = 0; i < MyXSlices; i++) {
                int XGlobal = i + MyXOffset;
                float DistX = (float)(XSpotPos - XGlobal);
                for (int j = 0; j < MyYSlices; j++) {
                    int YGlobal = j + MyYOffset;
                    float DistY = (float)(YSpotPos - YGlobal);
                    float TotDist = sqrt(DistX * DistX + DistY * DistY + DistZ * DistZ);
                    if (TotDist <= SpotRadius) {
                        int D3D1ConvPosition = k * MyXSlices * MyYSlices + i * MyYSlices + j;
                        // Melt time
                        LayerTimeTempHistory_Host(D3D1ConvPosition, NumberOfSolidificationEvents_Host(D3D1ConvPosition),
                                                  0) = 1 + TimeBetweenSpots * n;
                        // Liquidus time
                        LayerTimeTempHistory_Host(D3D1ConvPosition, NumberOfSolidificationEvents_Host(D3D1ConvPosition),
                                                  1) =
                            1 + (int)(((float)(SpotRadius)-TotDist) / IsothermVelocity) + TimeBetweenSpots * n;
                        // Cooling rate
                        LayerTimeTempHistory_Host(D3D1ConvPosition, NumberOfSolidificationEvents_Host(D3D1ConvPosition),
                                                  2) = R * deltat;
                        NumberOfSolidificationEvents_Host(D3D1ConvPosition)++;
                        Count++;
                    }
                }
            }
        }
    }

    // Initialize data for first melt-solidification event for all cells with data in this layer
    for (int k = 0; k <= SpotRadius; k++) {
        for (int i = 0; i < MyXSlices; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                int D3D1ConvPosition = k * MyXSlices * MyYSlices + i * MyYSlices + j;
                int GlobalD3D1ConvPosition = (k + ZBound_Low) * MyXSlices * MyYSlices + i * MyYSlices + j;
                if (NumberOfSolidificationEvents_Host(D3D1ConvPosition) > 0) {
                    MeltTimeStep_Host(GlobalD3D1ConvPosition) = LayerTimeTempHistory_Host(D3D1ConvPosition, 0, 0);
                    CritTimeStep_Host(GlobalD3D1ConvPosition) = LayerTimeTempHistory_Host(D3D1ConvPosition, 0, 1);
                    UndercoolingChange_Host(GlobalD3D1ConvPosition) = LayerTimeTempHistory_Host(D3D1ConvPosition, 0, 2);
                    LayerID_Host(GlobalD3D1ConvPosition) = layernumber;
                    Melted[GlobalD3D1ConvPosition] = true;
                }
                else {
                    MeltTimeStep_Host(GlobalD3D1ConvPosition) = 0;
                    CritTimeStep_Host(GlobalD3D1ConvPosition) = 0;
                    UndercoolingChange_Host(GlobalD3D1ConvPosition) = 0.0;
                }
            }
        }
    }

    // Initial undercooling of all cells is 0, solidification event counter is 0 at the start of each layer
    Kokkos::deep_copy(UndercoolingCurrent, 0.0);
    Kokkos::deep_copy(SolidificationEventCounter, 0.0);

    // Copy host view data back to device
    MaxSolidificationEvents =
        Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), MaxSolidificationEvents_Host);
    LayerID = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), LayerID_Host);
    MeltTimeStep = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), MeltTimeStep_Host);
    CritTimeStep = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), CritTimeStep_Host);
    UndercoolingChange = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), UndercoolingChange_Host);
    LayerTimeTempHistory =
        Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), LayerTimeTempHistory_Host);
    NumberOfSolidificationEvents =
        Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), NumberOfSolidificationEvents_Host);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Spot melt temperature field with remelting for layer " << layernumber
                  << " initialized; each cell will solidify up to " << MaxSolidificationEvents_Host(layernumber)
                  << " times" << std::endl;
}

// Initialize temperature data for a problem using the reduced/sparse data format and input temperature data from
// file(s)
void TempInit_Reduced(int id, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset, double deltax,
                      int HTtoCAratio, double deltat, int nz, int LocalDomainSize, ViewI &CritTimeStep,
                      ViewF &UndercoolingChange, float XMin, float YMin, float ZMin, bool *Melted, float *ZMinLayer,
                      float *ZMaxLayer, int LayerHeight, int NumberOfLayers, int *FinishTimeStep, double FreezingRange,
                      ViewI &LayerID, int *FirstValue, int *LastValue, std::vector<double> RawData) {

    // These views are initialized to zeros on the host, filled with data, and then copied to the device for layer
    // "layernumber"
    ViewI_H LayerID_Host("LayerID_H", LocalDomainSize);
    ViewI_H CritTimeStep_Host("CritTimeStep_H", LocalDomainSize);
    ViewF_H UndercoolingChange_Host("UndercoolingChange_H", LocalDomainSize);

    // Temperature data read
    // If HTtoCAratio > 1, an interpolation of input temperature data is needed
    // The X and Y bounds are the region (for this MPI rank) of the physical domain that needs to be
    // read extends past the actual spatial extent of the local domain for purposes of interpolating
    // from HT_deltax to deltax
    int LowerXBound = MyXOffset - (MyXOffset % HTtoCAratio);
    int LowerYBound = MyYOffset - (MyYOffset % HTtoCAratio);
    int UpperXBound = MyXOffset + MyXSlices - 1 + HTtoCAratio - ((MyXOffset + MyYSlices - 1) % HTtoCAratio);
    int UpperYBound = MyYOffset + MyYSlices - 1 + HTtoCAratio - ((MyYOffset + MyYSlices - 1) % HTtoCAratio);

    // No sites have melted yet
    for (int i = 0; i < MyXSlices * MyYSlices * nz; i++) {
        Melted[i] = false;
        LayerID_Host(i) = -1;
    }

    // removed unused LayerwiseTSOffset variable

    for (int LayerCounter = 0; LayerCounter < NumberOfLayers; LayerCounter++) {

        double SmallestTime = 1000000000;
        double SmallestTime_Global = 1000000000;
        double LargestTime = 0;
        double LargestTime_Global = 0;

        // How many CA cells in the vertical direction are needed to hold this layer's temperature data?
        int nzTempValuesThisLayer = round((ZMaxLayer[LayerCounter] - ZMinLayer[LayerCounter]) / deltax) + 1;
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
                XInt = round((RawData[i] - XMin) / deltax);
            }
            else if (Pos == 1) {
                YInt = round((RawData[i] - YMin) / deltax);
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
        // Previously normalized "time 0" to be the first time a cell goes below the liquidus; this is removed such that
        // time is consistent with simulations that include remelting, which don't normalize "time 0"
        FinishTimeStep[LayerCounter] = round(LargestTime_Global / deltat);
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
        // time step) "ZMin" is the global Z coordinate that corresponds to cells at Z = 0, "ZMax" is the global Z
        // coordinate that corresponds to cells at Z = nz-1
        if (id == 0)
            std::cout << "Layer " << LayerCounter << " data belongs to global z coordinates of "
                      << round((ZMinLayer[LayerCounter] - ZMin) / deltax) << " through "
                      << round((ZMinLayer[LayerCounter] - ZMin) / deltax) + nzTempValuesThisLayer - 1 << std::endl;

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
                            int ZOffset = round((ZMinLayer[LayerCounter] - ZMin) / deltax) + k;
                            int Coord3D1D = ZOffset * MyXSlices * MyYSlices + Adj_i * MyYSlices + Adj_j;
                            Melted[Coord3D1D] = true;
                            CritTimeStep_Host(Coord3D1D) = round(CTLiq / deltat);
                            LayerID_Host(Coord3D1D) = LayerCounter;
                            UndercoolingChange_Host(Coord3D1D) =
                                std::abs(CR[k][ii - LowerXBound][jj - LowerYBound]) * deltat;
                        }
                    }
                }
            }
        }
    } // End read over all temperature files and placement of data

    // Copy initialized host data back to device
    CritTimeStep = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), CritTimeStep_Host);
    LayerID = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), LayerID_Host);
    UndercoolingChange = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), UndercoolingChange_Host);
}

// Initialize temperature fields for this layer if remelting is considered and data comes from files
void TempInit_Remelt(int layernumber, int id, int MyXSlices, int MyYSlices, int nz, int LocalActiveDomainSize,
                     int LocalDomainSize, int MyXOffset, int MyYOffset, double &deltax, double deltat,
                     double FreezingRange, ViewF3D &LayerTimeTempHistory, ViewI &NumberOfSolidificationEvents,
                     ViewI MaxSolidificationEvents, ViewI &MeltTimeStep, ViewI &CritTimeStep, ViewF &UndercoolingChange,
                     ViewF &UndercoolingCurrent, float XMin, float YMin, bool *Melted, float *ZMinLayer,
                     int LayerHeight, int nzActive, int ZBound_Low, int *FinishTimeStep, ViewI &LayerID,
                     int *FirstValue, int *LastValue, std::vector<double> RawData, ViewI &SolidificationEventCounter) {

    // Resize device views to have sizes compatible with the temporary host views
    // Copy this view back to the host as part of resizing LayerTimeTempHistory
    ViewI_H MaxSolidificationEvents_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), MaxSolidificationEvents);
    Kokkos::resize(LayerTimeTempHistory, LocalActiveDomainSize, MaxSolidificationEvents_Host(0), 3);
    Kokkos::resize(NumberOfSolidificationEvents, LocalActiveDomainSize);
    Kokkos::resize(SolidificationEventCounter, LocalActiveDomainSize);
    if (layernumber == 0) {
        // Only needs to be resized during initialization of the first layer, as LocalDomainSize is constant while
        // LocalActiveDomainSize is not
        Kokkos::resize(MeltTimeStep, LocalDomainSize);
    }

    // These views are copied to the host, updated for layer "layernumber", and later copied back to the device
    ViewI_H LayerID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), LayerID);
    ViewI_H MeltTimeStep_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), MeltTimeStep);
    ViewI_H CritTimeStep_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CritTimeStep);
    ViewF_H UndercoolingChange_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), UndercoolingChange);

    // These views are initialized to zeros on the host, filled with data, and then copied to the device for layer
    // "layernumber"
    ViewF3D_H LayerTimeTempHistory_Host("TimeTempHistory_H", LocalActiveDomainSize, MaxSolidificationEvents_Host(0), 3);
    ViewI_H NumberOfSolidificationEvents_Host("NumSEvents_H", LocalActiveDomainSize);

    if (layernumber == 0) {
        // No sites have melted yet
        for (int i = 0; i < MyXSlices * MyYSlices * nz; i++) {
            Melted[i] = false;
            LayerID_Host(i) = -1;
        }
    }

    // Data was already read into the "RawData" temporary data structure
    // Determine which section of "RawData" is relevant for this layer of the overall domain
    int StartRange = FirstValue[layernumber];
    int EndRange = LastValue[layernumber];
    int XInt = -1;
    int YInt = -1;
    int ZInt = -1;
    int D3D1ConvPosition = 0;
    double LargestTime = 0;
    double LargestTime_Global = 0;
    if (id == 0)
        std::cout << "Range of raw data for layer " << layernumber << " on rank 0 is " << StartRange << " to "
                  << EndRange << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = StartRange; i < EndRange; i++) {

        int Pos = i % 6;
        if (Pos == 0) {
            XInt = round((RawData[i] - XMin) / deltax);
        }
        else if (Pos == 1) {
            YInt = round((RawData[i] - YMin) / deltax);
        }
        else if (Pos == 2) {
            ZInt = round((RawData[i] + deltax * LayerHeight * layernumber - ZMinLayer[layernumber]) / deltax);
        }
        else if (Pos == 3) {
            // Melt time step (smallest possible value is time step 1, if time = 0)
            if ((XInt >= MyXOffset) && (XInt < MyXOffset + MyXSlices) && (YInt >= MyYOffset) &&
                (YInt < MyYOffset + MyYSlices)) {
                D3D1ConvPosition = ZInt * MyXSlices * MyYSlices + (XInt - MyXOffset) * MyYSlices + (YInt - MyYOffset);
                LayerTimeTempHistory_Host(D3D1ConvPosition, NumberOfSolidificationEvents_Host(D3D1ConvPosition), 0) =
                    round(RawData[i] / deltat) + 1;
            }
            else {
                // skip to next data point
                i = i + 2;
            }
        }
        else if (Pos == 4) {
            // Crit (liquidus) time step
            LayerTimeTempHistory_Host(D3D1ConvPosition, NumberOfSolidificationEvents_Host(D3D1ConvPosition), 1) =
                round(RawData[i] / deltat) + 1;
        }
        else if (Pos == 5) {
            // Cooling rate per time step
            LayerTimeTempHistory_Host(D3D1ConvPosition, NumberOfSolidificationEvents_Host(D3D1ConvPosition), 2) =
                std::abs(RawData[i]) * deltat;
            NumberOfSolidificationEvents_Host(D3D1ConvPosition)++;
            double SolidusTime = RawData[i - 1] + FreezingRange / RawData[i];
            if (SolidusTime > LargestTime)
                LargestTime = SolidusTime;
        }
    }
    MPI_Allreduce(&LargestTime, &LargestTime_Global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Largest time globally for layer " << layernumber << " is " << LargestTime_Global << std::endl;
    FinishTimeStep[layernumber] = round((LargestTime_Global) / deltat);
    if (id == 0)
        std::cout << " Layer " << layernumber << " FINISH TIME STEP IS " << FinishTimeStep[layernumber] << std::endl;
    if (id == 0)
        std::cout << "Layer " << layernumber << " temperatures read" << std::endl;

    // Reorder solidification events in LayerTimeTempHistory(location,event number,component) so that they are in order
    // based on the melting time values (component = 0)
    for (int n = 0; n < LocalActiveDomainSize; n++) {
        if (NumberOfSolidificationEvents_Host(n) > 0) {
            for (int i = 0; i < NumberOfSolidificationEvents_Host(n) - 1; i++) {
                for (int j = (i + 1); j < NumberOfSolidificationEvents_Host(n); j++) {
                    if (LayerTimeTempHistory_Host(n, i, 0) > LayerTimeTempHistory_Host(n, j, 0)) {
                        // Swap these two points - melting event "j" happens before event "i"
                        float OldMeltVal = LayerTimeTempHistory_Host(n, i, 0);
                        float OldLiqVal = LayerTimeTempHistory_Host(n, i, 1);
                        float OldCRVal = LayerTimeTempHistory_Host(n, i, 2);
                        LayerTimeTempHistory_Host(n, i, 0) = LayerTimeTempHistory_Host(n, j, 0);
                        LayerTimeTempHistory_Host(n, i, 1) = LayerTimeTempHistory_Host(n, j, 1);
                        LayerTimeTempHistory_Host(n, i, 2) = LayerTimeTempHistory_Host(n, j, 2);
                        LayerTimeTempHistory_Host(n, j, 0) = OldMeltVal;
                        LayerTimeTempHistory_Host(n, j, 1) = OldLiqVal;
                        LayerTimeTempHistory_Host(n, j, 2) = OldCRVal;
                    }
                }
            }
        }
    }
    // If a cell melts twice before reaching the liquidus temperature, this is a double counted solidification
    // event and should be removed
    for (int n = 0; n < LocalActiveDomainSize; n++) {
        if (NumberOfSolidificationEvents_Host(n) > 1) {
            for (int i = 0; i < NumberOfSolidificationEvents_Host(n) - 1; i++) {
                if (LayerTimeTempHistory_Host(n, i + 1, 0) < LayerTimeTempHistory_Host(n, i, 1)) {
                    std::cout << "Cell " << n << " removing anomalous event " << i + 1 << " out of "
                              << NumberOfSolidificationEvents_Host(n) - 1 << std::endl;
                    // Keep whichever event has the larger liquidus time
                    if (LayerTimeTempHistory_Host(n, i + 1, 1) > LayerTimeTempHistory_Host(n, i, 1)) {
                        LayerTimeTempHistory_Host(n, i, 0) = LayerTimeTempHistory_Host(n, i + 1, 0);
                        LayerTimeTempHistory_Host(n, i, 1) = LayerTimeTempHistory_Host(n, i + 1, 1);
                        LayerTimeTempHistory_Host(n, i, 2) = LayerTimeTempHistory_Host(n, i + 1, 2);
                    }
                    LayerTimeTempHistory_Host(n, i + 1, 0) = 0.0;
                    LayerTimeTempHistory_Host(n, i + 1, 1) = 0.0;
                    LayerTimeTempHistory_Host(n, i + 1, 2) = 0.0;
                    // Reshuffle other solidification events over if needed
                    for (int ii = (i + 1); ii < NumberOfSolidificationEvents_Host(n) - 1; ii++) {
                        LayerTimeTempHistory_Host(n, ii, 0) = LayerTimeTempHistory_Host(n, ii + 1, 0);
                        LayerTimeTempHistory_Host(n, ii, 1) = LayerTimeTempHistory_Host(n, ii + 1, 1);
                        LayerTimeTempHistory_Host(n, ii, 2) = LayerTimeTempHistory_Host(n, ii + 1, 2);
                    }
                    NumberOfSolidificationEvents_Host(n)--;
                }
            }
        }
    }
    // First melt-solidification event from LayerTimeTempHistory to happen is initialized
    for (int k = 0; k < nzActive; k++) {
        int GlobalZ = k + ZBound_Low;
        for (int i = 0; i < MyXSlices; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                int D3D1ConvPosition = k * MyXSlices * MyYSlices + i * MyYSlices + j;
                int GlobalD3D1ConvPosition = GlobalZ * MyXSlices * MyYSlices + i * MyYSlices + j;
                if (LayerTimeTempHistory_Host(D3D1ConvPosition, 0, 0) > 0) {
                    // This cell undergoes solidification in layer "layernumber" at least once
                    Melted[GlobalD3D1ConvPosition] = true;
                    LayerID_Host(GlobalD3D1ConvPosition) = layernumber;
                    MeltTimeStep_Host(GlobalD3D1ConvPosition) =
                        (int)(LayerTimeTempHistory_Host(D3D1ConvPosition, 0, 0));
                    CritTimeStep_Host(GlobalD3D1ConvPosition) =
                        (int)(LayerTimeTempHistory_Host(D3D1ConvPosition, 0, 1));
                    UndercoolingChange_Host(GlobalD3D1ConvPosition) = LayerTimeTempHistory_Host(D3D1ConvPosition, 0, 2);
                }
                else {
                    // This cell does not undergo solidification in layer "layernumber"
                    MeltTimeStep_Host(GlobalD3D1ConvPosition) = 0;
                    CritTimeStep_Host(GlobalD3D1ConvPosition) = 0;
                    UndercoolingChange_Host(GlobalD3D1ConvPosition) = 0.0;
                }
            }
        }
    }

    // Initial undercooling of all cells is 0, solidification event counter is 0 at the start of each layer
    Kokkos::deep_copy(UndercoolingCurrent, 0.0);
    Kokkos::deep_copy(SolidificationEventCounter, 0.0);

    // Copy host view data back to device
    MaxSolidificationEvents =
        Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), MaxSolidificationEvents_Host);
    LayerID = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), LayerID_Host);
    MeltTimeStep = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), MeltTimeStep_Host);
    CritTimeStep = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), CritTimeStep_Host);
    UndercoolingChange = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), UndercoolingChange_Host);
    LayerTimeTempHistory =
        Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), LayerTimeTempHistory_Host);
    NumberOfSolidificationEvents =
        Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), NumberOfSolidificationEvents_Host);

    if (id == 0)
        std::cout << "Layer " << layernumber << " temperature field is from Z = " << ZBound_Low << " through "
                  << nzActive + ZBound_Low - 1 << " of the global domain" << std::endl;
}

//*****************************************************************************/
// Initialize grain orientations and unit vectors
void OrientationInit(int, int NGrainOrientations, ViewF &GrainUnitVector, std::string GrainOrientationFile) {

    // Temporary host view for storing grain orientations read from file
    ViewF_H GrainUnitVector_Host(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector_H"), 9 * NGrainOrientations);

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
            GrainUnitVector_Host(9 * i + 3 * UVNumber + Comp) = ReadGO;

            Comp++;
            if (Comp > 2) {
                Comp = 0;
                UVNumber++;
            }
        }
    }
    O.close();

    // Copy orientation data to device
    GrainUnitVector = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), GrainUnitVector_Host);
}

// Initializes cell types and epitaxial Grain ID values where substrate grains are active cells on the bottom surface of
// the constrained domain. Also initialize active cell data structures associated with the substrate grains
void SubstrateInit_ConstrainedGrowth(int id, double FractSurfaceSitesActive, int MyXSlices, int MyYSlices, int nx,
                                     int ny, int MyXOffset, int MyYOffset, ViewI NeighborX, ViewI NeighborY,
                                     ViewI NeighborZ, ViewF GrainUnitVector, int NGrainOrientations, ViewI CellType,
                                     ViewI GrainID, ViewF DiagonalLength, ViewF DOCenter, ViewF CritDiagonalLength,
                                     double RNGSeed, int np, int DecompositionStrategy, Buffer2D BufferWestSend,
                                     Buffer2D BufferEastSend, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend,
                                     Buffer2D BufferNorthEastSend, Buffer2D BufferNorthWestSend,
                                     Buffer2D BufferSouthEastSend, Buffer2D BufferSouthWestSend, int BufSizeX,
                                     int BufSizeY, bool AtNorthBoundary, bool AtSouthBoundary, bool AtEastBoundary,
                                     bool AtWestBoundary) {

    // Calls to Xdist(gen) and Y dist(gen) return random locations for grain seeds
    // Since X = 0 and X = nx-1 are the cell centers of the last cells in X, locations are evenly scattered between X =
    // -0.49999 and X = nx - 0.5, as the cells have a half width of 0.5
    std::mt19937_64 gen(RNGSeed);
    std::uniform_real_distribution<double> Xdist(-0.49999, nx - 0.5);
    std::uniform_real_distribution<double> Ydist(-0.49999, ny - 0.5);

    // Determine number of active cells from the fraction of sites active and the number of sites on the bottom domain
    // surface
    int SubstrateActCells = std::round(FractSurfaceSitesActive * nx * ny);

    // On all ranks, generate active site locations - same list on every rank
    // TODO: Generate random numbers on GPU, instead of using host view and copying over - ensure that locations are in
    // the same order every time based on the RNGSeed
    ViewI_H ActCellX_Host(Kokkos::ViewAllocateWithoutInitializing("ActCellX_Host"), SubstrateActCells);
    ViewI_H ActCellY_Host(Kokkos::ViewAllocateWithoutInitializing("ActCellY_Host"), SubstrateActCells);

    // Randomly locate substrate grain seeds for cells in the interior of this subdomain (at the k = 0 bottom surface)
    for (int n = 0; n < SubstrateActCells; n++) {
        double XLocation = Xdist(gen);
        double YLocation = Ydist(gen);
        // Randomly select integer coordinates between 0 and nx-1 or ny-1
        ActCellX_Host(n) = round(XLocation);
        ActCellY_Host(n) = round(YLocation);
    }

    // Copy views of substrate grain locations back to the device
    using memory_space = Kokkos::DefaultExecutionSpace::memory_space;
    ViewI ActCellX_Device = Kokkos::create_mirror_view_and_copy(memory_space(), ActCellX_Host);
    ViewI ActCellY_Device = Kokkos::create_mirror_view_and_copy(memory_space(), ActCellY_Host);

    // Start with all cells as liquid prior to locating substrate grain seeds
    Kokkos::deep_copy(CellType, Liquid);

    // Determine which grains/active cells belong to which MPI ranks
    Kokkos::parallel_for(
        "ConstrainedGrainInit", SubstrateActCells, KOKKOS_LAMBDA(const int &n) {
            // What are the X and Y coordinates of this active cell relative to the X and Y bounds of this rank?
            if ((ActCellX_Device(n) >= MyXOffset) && (ActCellX_Device(n) < MyXOffset + MyXSlices) &&
                (ActCellY_Device(n) >= MyYOffset) && (ActCellY_Device(n) < MyYOffset + MyYSlices)) {
                // Convert X and Y coordinates to values relative to this MPI rank's grid (Z = 0 for these active cells,
                // at bottom surface) GrainIDs come from the position on the list of substrate active cells to avoid
                // reusing the same value
                int LocalX = ActCellX_Device(n) - MyXOffset;
                int LocalY = ActCellY_Device(n) - MyYOffset;
                int D3D1ConvPosition = LocalX * MyYSlices + LocalY;
                CellType(D3D1ConvPosition) = Active;
                GrainID(D3D1ConvPosition) = n + 1; // assign GrainID > 0 to epitaxial seeds
                // Initialize active cell data structures
                int GlobalX = LocalX + MyXOffset;
                int GlobalY = LocalY + MyYOffset;
                int GlobalZ = 0;
                DiagonalLength(D3D1ConvPosition) = 0.01;
                DOCenter(3 * D3D1ConvPosition) = GlobalX + 0.5;
                DOCenter(3 * D3D1ConvPosition + 1) = GlobalY + 0.5;
                DOCenter(3 * D3D1ConvPosition + 2) = GlobalZ + 0.5;

                // The orientation for the new grain will depend on its Grain ID
                int MyOrientation = getGrainOrientation(n + 1, NGrainOrientations);
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
                for (int l = 0; l < 26; l++) {

                    // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the image of
                    // the center of a neighbor cell
                    double x0 = xp + NeighborX(l) - cx;
                    double y0 = yp + NeighborY(l) - cy;
                    double z0 = zp + NeighborZ(l) - cz;
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
                    CritDiagonalLength((long int)(26) * D3D1ConvPosition + (long int)(l)) = CDLVal;
                } // end loop over 26 diagonals
                // If this new active cell is in the halo region, load the send buffers
                if (np > 1) {

                    float GhostGID = GrainID(D3D1ConvPosition);
                    float GhostDOCX = GlobalX + 0.5;
                    float GhostDOCY = GlobalY + 0.5;
                    float GhostDOCZ = GlobalZ + 0.5;
                    float GhostDL = 0.01;
                    // Collect data for the ghost nodes, if necessary
                    if (DecompositionStrategy == 1)
                        loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, BufSizeX, MyYSlices, LocalX,
                                       LocalY, 0, AtNorthBoundary, AtSouthBoundary, BufferSouthSend, BufferNorthSend);
                    else
                        loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, BufSizeX, BufSizeY,
                                       MyXSlices, MyYSlices, LocalX, LocalY, 0, AtNorthBoundary, AtSouthBoundary,
                                       AtWestBoundary, AtEastBoundary, BufferSouthSend, BufferNorthSend, BufferWestSend,
                                       BufferEastSend, BufferNorthEastSend, BufferSouthEastSend, BufferSouthWestSend,
                                       BufferNorthWestSend);
                } // End if statement for serial/parallel code
            }
        });
    if (id == 0)
        std::cout << "Number of substrate active cells across all ranks: " << SubstrateActCells << std::endl;
}

// Initializes Grain ID values where the substrate comes from a file
void SubstrateInit_FromFile(std::string SubstrateFileName, int nz, int MyXSlices, int MyYSlices, int MyXOffset,
                            int MyYOffset, int id, ViewI &GrainID_Device) {

    // Assign GrainID values to cells that are part of the substrate - read values from file and initialize using
    // temporary host view
    ViewI_H GrainID_Host(Kokkos::ViewAllocateWithoutInitializing("GrainID_Host"), MyXSlices * MyYSlices * nz);
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

    // Assign GrainID values to all cells - cells that will be part of the melt pool footprint may still need their
    // initial GrainID
    for (int k = 0; k < nzS; k++) {
        if (k == nz)
            break;
        for (int j = 0; j < nyS; j++) {
            for (int i = 0; i < nxS; i++) {
                std::string GIDVal;
                getline(Substrate, GIDVal);
                if ((i >= Substrate_LowX) && (i < Substrate_HighX) && (j >= Substrate_LowY) && (j < Substrate_HighY)) {
                    int CAGridLocation;
                    CAGridLocation = k * MyXSlices * MyYSlices + (i - MyXOffset) * MyYSlices + (j - MyYOffset);
                    GrainID_Host(CAGridLocation) = stoi(GIDVal, nullptr, 10);
                }
            }
        }
    }
    Substrate.close();
    // Copy GrainIDs read from file to device
    GrainID_Device = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), GrainID_Host);
    if (id == 0)
        std::cout << "Substrate file read complete" << std::endl;
}

// Initializes Grain ID values where the baseplate is generated using an input grain spacing and a Voronoi Tessellation
void BaseplateInit_FromGrainSpacing(float SubstrateGrainSpacing, int nx, int ny, float *ZMinLayer, float *ZMaxLayer,
                                    int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int id, double deltax,
                                    ViewI GrainID, double RNGSeed, int &NextLayer_FirstEpitaxialGrainID) {

    // Seed random number generator such that each rank generates the same baseplate grain center locations
    // Calls to Xdist(gen), Ydist(gen), Zdist(gen) return random locations for grain seeds
    // Since X = 0 and X = nx-1 are the cell centers of the last cells in X, locations are evenly scattered between X =
    // -0.49999 and X = nx - 0.5, as the cells have a half width of 0.5
    // Baseplate size in Z corresponds to the Z coordinates of layer 0
    int BaseplateSizeZ = round((ZMaxLayer[0] - ZMinLayer[0]) / deltax) + 1;
    std::mt19937_64 gen(RNGSeed);
    std::uniform_real_distribution<double> Xdist(-0.49999, nx - 0.5);
    std::uniform_real_distribution<double> Ydist(-0.49999, ny - 0.5);
    std::uniform_real_distribution<double> Zdist(-0.49999, BaseplateSizeZ - 0.5);

    // Based on the baseplate size (domain of first layer) and the substrate grain spacing, determine the number of
    // baseplate grains
    double BaseplateSize = nx * ny * BaseplateSizeZ * pow(deltax, 3) * pow(10, 18); // in cubic microns
    double SubstrateGrainSize = pow(SubstrateGrainSpacing, 3);                      // in cubic microns
    int NumberOfBaseplateGrains = round(BaseplateSize / SubstrateGrainSize);
    // Need at least 1 baseplate grain
    NumberOfBaseplateGrains = std::max(NumberOfBaseplateGrains, 1);
    ViewI_H NumBaseplateGrains_Host(Kokkos::ViewAllocateWithoutInitializing("NBaseplate_Host"), 1);
    NumBaseplateGrains_Host(0) = NumberOfBaseplateGrains;
    // TODO: Use device RNG to generate baseplate grain locations, instead of host with copy
    // Store grain centers as float coordinates, not integer cell locations, for more precise substrate generation
    ViewF_H BaseplateGrainX_Host(Kokkos::ViewAllocateWithoutInitializing("BaseplateGrainX_Host"),
                                 nx * ny * BaseplateSizeZ);
    ViewF_H BaseplateGrainY_Host(Kokkos::ViewAllocateWithoutInitializing("BaseplateGrainY_Host"),
                                 nx * ny * BaseplateSizeZ);
    ViewF_H BaseplateGrainZ_Host(Kokkos::ViewAllocateWithoutInitializing("BaseplateGrainZ_Host"),
                                 nx * ny * BaseplateSizeZ);

    // For the entire baseplate (all x and y coordinate, but only layer 0 z coordinates), identify baseplate grain
    // centers
    if (id == 0)
        std::cout << "Baseplate spanning domain coordinates Z = 0 through " << BaseplateSizeZ - 1 << std::endl;
    for (int n = 0; n < NumberOfBaseplateGrains; n++) {
        BaseplateGrainX_Host(n) = Xdist(gen);
        BaseplateGrainY_Host(n) = Ydist(gen);
        BaseplateGrainZ_Host(n) = Zdist(gen);
    }
    if (id == 0)
        std::cout << "Number of baseplate grains: " << NumBaseplateGrains_Host(0) << std::endl;
    using memory_space = Kokkos::DefaultExecutionSpace::memory_space;
    ViewI NumBaseplateGrains_Device = Kokkos::create_mirror_view_and_copy(memory_space(), NumBaseplateGrains_Host);
    ViewF BaseplateGrainX_Device = Kokkos::create_mirror_view_and_copy(memory_space(), BaseplateGrainX_Host);
    ViewF BaseplateGrainY_Device = Kokkos::create_mirror_view_and_copy(memory_space(), BaseplateGrainY_Host);
    ViewF BaseplateGrainZ_Device = Kokkos::create_mirror_view_and_copy(memory_space(), BaseplateGrainZ_Host);

    //  Initialize GrainIDs on device for baseplate
    Kokkos::parallel_for(
        "BaseplateGen",
        Kokkos::MDRangePolicy<Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>(
            {0, 0, 0}, {BaseplateSizeZ, MyXSlices, MyYSlices}),
        KOKKOS_LAMBDA(const int k, const int i, const int j) {
            int GlobalX = i + MyXOffset;
            int GlobalY = j + MyYOffset;
            int CAGridLocation = k * MyXSlices * MyYSlices + i * MyYSlices + j;
            // All cells are given GrainID values, even if they're part of the melt pool footprint, as they may need the
            // values later if at the interface This cell is part of the substrate - determine which grain center the
            // cell is closest to, in order to assign it a grain ID If closest to grain "n", assign grain ID "n+1"
            // (grain ID = 0 is not used)
            float MinDistanceToThisGrain = (float)(nx * ny * BaseplateSizeZ);
            int ClosestGrainIndex = -1;
            for (int n = 0; n < NumBaseplateGrains_Device(0); n++) {
                float DistanceToThisGrainX = (float)(std::abs(BaseplateGrainX_Device(n) - GlobalX));
                float DistanceToThisGrainY = (float)(std::abs(BaseplateGrainY_Device(n) - GlobalY));
                float DistanceToThisGrainZ = (float)(std::abs(BaseplateGrainZ_Device(n) - k));
                float DistanceToThisGrain =
                    sqrtf(DistanceToThisGrainX * DistanceToThisGrainX + DistanceToThisGrainY * DistanceToThisGrainY +
                          DistanceToThisGrainZ * DistanceToThisGrainZ);
                if (DistanceToThisGrain < MinDistanceToThisGrain) {
                    ClosestGrainIndex = n;
                    MinDistanceToThisGrain = DistanceToThisGrain;
                }
            }
            GrainID(CAGridLocation) = ClosestGrainIndex + 1;
        });

    NextLayer_FirstEpitaxialGrainID =
        NumberOfBaseplateGrains + 2; // avoid reusing GrainID in next layer's powder grain structure
    if (id == 0)
        std::cout << "Baseplate grain structure initialized" << std::endl;
}

// Each layer's top Z coordinates are seeded with CA-cell sized substrate grains (emulating bulk nucleation alongside
// the edges of partially melted powder particles)
void PowderInit(int layernumber, int nx, int ny, int LayerHeight, float *ZMaxLayer, float ZMin, double deltax,
                int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int id, ViewI GrainID, double RNGSeed,
                int &NextLayer_FirstEpitaxialGrainID, double FractPowderSitesActive) {

    // On all ranks, generate list of powder grain IDs (starting with NextLayer_FirstEpitaxialGrainID, and shuffle them
    // so that their locations aren't sequential and depend on the RNGSeed (different for each layer)
    std::mt19937_64 gen(RNGSeed + layernumber);
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // TODO: This should be performed on the device, rather than the host
    // Number of cells in the powder layer
    int PowderLayerCells = nx * ny * LayerHeight;
    // Counter for the number of cells in the powder layer that are the source of a new grain
    int PowderLayerActiveCells = 0;
    // List of GrainIDs: a non-zero value is assigned to cells that are the source of a new grain
    std::vector<int> PowderGrainIDs(PowderLayerCells, 0);
    for (int n = 0; n < PowderLayerCells; n++) {
        double ProbActive = dis(gen);
        if (ProbActive < FractPowderSitesActive) {
            // This cell is the source of a new grain from the powder layer
            PowderGrainIDs[n] = PowderLayerActiveCells + NextLayer_FirstEpitaxialGrainID;
            PowderLayerActiveCells++;
        }
    }
    std::shuffle(PowderGrainIDs.begin(), PowderGrainIDs.end(), gen);
    // Copy powder layer GrainIDs into a host view, then a device view
    ViewI_H PowderGrainIDs_Host(Kokkos::ViewAllocateWithoutInitializing("PowderGrainIDs_Host"), nx * ny * LayerHeight);
    for (int n = 0; n < PowderLayerCells; n++) {
        PowderGrainIDs_Host(n) = PowderGrainIDs[n];
    }
    using memory_space = Kokkos::DefaultExecutionSpace::memory_space;
    ViewI PowderGrainIDs_Device = Kokkos::create_mirror_view_and_copy(memory_space(), PowderGrainIDs_Host);

    // Associate powder grain IDs with CA cells in the powder layer
    // Use bounds from temperature field for this layer to determine which cells are part of the powder
    int PowderTopZ = round((ZMaxLayer[layernumber] - ZMin) / deltax) + 1;
    int PowderBottomZ = PowderTopZ - LayerHeight;
    if (id == 0)
        std::cout << "Initializing powder layer for Z = " << PowderBottomZ << " through " << PowderTopZ - 1
                  << std::endl;

    int PowderStart = nx * ny * PowderBottomZ;
    int PowderEnd = nx * ny * PowderTopZ;

    Kokkos::parallel_for(
        "PowderGrainInit", Kokkos::RangePolicy<>(PowderStart, PowderEnd), KOKKOS_LAMBDA(const int &n) {
            int GlobalZ = n / (nx * ny);
            int Rem = n % (nx * ny);
            int GlobalX = Rem / ny;
            int GlobalY = Rem % ny;
            // Is this powder coordinate in X and Y in bounds for this rank? Is the grain id of this site unassigned
            // (wasn't captured during solidification of the previous layer)?
            int GlobalD3D1ConvPosition =
                GlobalZ * MyXSlices * MyYSlices + (GlobalX - MyXOffset) * MyYSlices + (GlobalY - MyYOffset);
            if ((GlobalX >= MyXOffset) && (GlobalX < MyXOffset + MyXSlices) && (GlobalY >= MyYOffset) &&
                (GlobalY < MyYOffset + MyYSlices) && (GrainID(GlobalD3D1ConvPosition) == 0))
                GrainID(GlobalD3D1ConvPosition) = PowderGrainIDs_Device(n - PowderStart);
        });
    Kokkos::fence();

    // Update NextLayer_FirstEpitaxialGrainID for next layer
    NextLayer_FirstEpitaxialGrainID += PowderLayerActiveCells;
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Initialized powder grain structure for layer " << layernumber << std::endl;
}

//*****************************************************************************/
// Initializes cells at border of solid and liquid as active type - performed on device
void CellTypeInit(int layernumber, int id, int np, int DecompositionStrategy, int MyXSlices, int MyYSlices,
                  int MyXOffset, int MyYOffset, int ZBound_Low, int nz, int LocalActiveDomainSize, int LocalDomainSize,
                  ViewI CellType, ViewI CritTimeStep, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ,
                  int NGrainOrientations, ViewF GrainUnitVector, ViewF DiagonalLength, ViewI GrainID,
                  ViewF CritDiagonalLength, ViewF DOCenter, ViewI LayerID, Buffer2D BufferWestSend,
                  Buffer2D BufferEastSend, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend,
                  Buffer2D BufferNorthEastSend, Buffer2D BufferNorthWestSend, Buffer2D BufferSouthEastSend,
                  Buffer2D BufferSouthWestSend, int BufSizeX, int BufSizeY, bool AtNorthBoundary, bool AtSouthBoundary,
                  bool AtEastBoundary, bool AtWestBoundary) {

    // Start with all cells as solid for the first layer, with liquid cells where temperature data exists
    if (layernumber == 0) {
        Kokkos::deep_copy(CellType, Solid);
        // This is done prior to initializing active cells, as active cells are initialized based on neighbor cell types
        Kokkos::parallel_for(
            "CellTypeInitSolLiq", LocalDomainSize, KOKKOS_LAMBDA(const int &GlobalD3D1ConvPosition) {
                if (CritTimeStep(GlobalD3D1ConvPosition) != 0)
                    CellType(GlobalD3D1ConvPosition) = Liquid;
            });
        Kokkos::fence();
        int ActCellCount = 0;
        Kokkos::parallel_reduce(
            "CellTypeInitAct", LocalDomainSize,
            KOKKOS_LAMBDA(const int &GlobalD3D1ConvPosition, int &ActCellCount) {
                // Cells of interest for the CA
                int GlobalZ = GlobalD3D1ConvPosition / (MyXSlices * MyYSlices);
                int Rem = GlobalD3D1ConvPosition % (MyXSlices * MyYSlices);
                int RankX = Rem / MyYSlices;
                int RankY = Rem % MyYSlices;
                if (CellType(GlobalD3D1ConvPosition) == Liquid) {
                    // This is a liquid or active cell, depending on whether it is located at the interface of the
                    // solid Check to see if this site is actually at the solid-liquid interface "l" corresponds to
                    // a specific neighbor direction (of 26 possible neighbors)
                    // All liquid cells at the bottom surface (k = 0) are at the interface and should be active
                    // Without solid/wall cell padding below the domain, some of these wouldn't be caught by
                    // the logic in the below loop over all neighbors
                    for (int l = 0; l < 26; l++) {
                        // Global coordinates of adjacent cell center
                        int MyNeighborX = RankX + NeighborX(l);
                        int MyNeighborY = RankY + NeighborY(l);
                        int MyNeighborZ = GlobalZ + NeighborZ(l);
                        int NeighborD3D1ConvPosition =
                            MyNeighborZ * MyXSlices * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
                        if ((MyNeighborX >= 0) && (MyNeighborX < MyXSlices) && (MyNeighborY >= 0) &&
                            (MyNeighborY < MyYSlices) && (MyNeighborZ >= 0) && (MyNeighborZ < nz)) {
                            if ((CellType(NeighborD3D1ConvPosition) == Solid) || (GlobalZ == 0)) {
                                // This cell is at the interface - becomes active type
                                CellType(GlobalD3D1ConvPosition) = Active;
                                ActCellCount++;
                                // exit loop over neighboring cells
                                l = 26;
                            }
                        }
                    }
                }
            },
            ActCellCount);
        Kokkos::fence();
        int TotalSubstrateActCells;
        MPI_Reduce(&ActCellCount, &TotalSubstrateActCells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Number of initial substrate active cells across all ranks and all layers: "
                      << TotalSubstrateActCells << std::endl;
    }

    // Each layer, count number of active cells are at the solid-liquid boundary for this layer's portion of the domain
    Kokkos::parallel_for(
        "CellTypeInitAct", LocalActiveDomainSize, KOKKOS_LAMBDA(const int &D3D1ConvPosition) {
            // Cells of interest for the CA
            int RankZ = D3D1ConvPosition / (MyXSlices * MyYSlices);
            int Rem = D3D1ConvPosition % (MyXSlices * MyYSlices);
            int RankX = Rem / MyYSlices;
            int RankY = Rem % MyYSlices;
            int GlobalZ = RankZ + ZBound_Low;
            int GlobalD3D1ConvPosition = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
            if ((CellType(GlobalD3D1ConvPosition) == Active) && (LayerID(GlobalD3D1ConvPosition) == layernumber)) {
                // This cell was marked as active previously - is there a non-zero GrainID at this site (either from a
                // previous layer's growth or from the initialized powder layer)? If so, initialize active cell data
                // structures. Otherwise, return this cell to the liquid state
                int MyGrainID = GrainID(GlobalD3D1ConvPosition);
                if (MyGrainID == 0)
                    CellType(GlobalD3D1ConvPosition) = Liquid;
                else {
                    int GlobalX = RankX + MyXOffset;
                    int GlobalY = RankY + MyYOffset;
                    int RankZ = GlobalZ - ZBound_Low;
                    int D3D1ConvPosition = RankZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;

                    DiagonalLength(D3D1ConvPosition) = 0.01;
                    DOCenter(3 * D3D1ConvPosition) = GlobalX + 0.5;
                    DOCenter(3 * D3D1ConvPosition + 1) = GlobalY + 0.5;
                    DOCenter(3 * D3D1ConvPosition + 2) = GlobalZ + 0.5;

                    // The orientation for the new grain will depend on its Grain ID
                    int MyOrientation = getGrainOrientation(MyGrainID, NGrainOrientations);
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
                    } // end loop over 26 diagonals
                    // If this new active cell is in the halo region, load the send buffers
                    if (np > 1) {

                        double GhostGID = static_cast<double>(MyGrainID);
                        double GhostDOCX = static_cast<double>(GlobalX + 0.5);
                        double GhostDOCY = static_cast<double>(GlobalY + 0.5);
                        double GhostDOCZ = static_cast<double>(GlobalZ + 0.5);
                        double GhostDL = 0.01;
                        // Collect data for the ghost nodes, if necessary
                        if (DecompositionStrategy == 1)
                            loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, BufSizeX, MyYSlices,
                                           RankX, RankY, RankZ, AtNorthBoundary, AtSouthBoundary, BufferSouthSend,
                                           BufferNorthSend);
                        else
                            loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, BufSizeX, BufSizeY,
                                           MyXSlices, MyYSlices, RankX, RankY, RankZ, AtNorthBoundary, AtSouthBoundary,
                                           AtWestBoundary, AtEastBoundary, BufferSouthSend, BufferNorthSend,
                                           BufferWestSend, BufferEastSend, BufferNorthEastSend, BufferSouthEastSend,
                                           BufferSouthWestSend, BufferNorthWestSend);
                    } // End if statement for serial/parallel code
                }
            }
        });
    Kokkos::fence();
}

//*****************************************************************************/
// Initializes cells for the current layer as either solid (don't resolidify) or tempsolid (will melt and resolidify)
void CellTypeInit_Remelt(int MyXSlices, int MyYSlices, int LocalActiveDomainSize, ViewI CellType, ViewI CritTimeStep,
                         int id, int ZBound_Low) {

    int MeltPoolCellCount;
    Kokkos::parallel_reduce(
        "CellTypeInitSolidRM", LocalActiveDomainSize,
        KOKKOS_LAMBDA(const int &D3D1ConvPosition, int &tmpval) {
            int GlobalD3D1ConvPosition = D3D1ConvPosition + ZBound_Low * MyXSlices * MyYSlices;
            if (CritTimeStep(GlobalD3D1ConvPosition) != 0) {
                CellType(GlobalD3D1ConvPosition) = TempSolid;
                tmpval++;
            }
            else
                CellType(GlobalD3D1ConvPosition) = Solid;
        },
        MeltPoolCellCount);
    int TotalMeltPoolCellCount;
    MPI_Reduce(&MeltPoolCellCount, &TotalMeltPoolCellCount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Number of cells across all ranks to undergo solidification at least once: "
                  << TotalMeltPoolCellCount << std::endl;
}

//*****************************************************************************/
// Initialize nucleation site locations, GrainID values, and time at which nucleation events will potentially occur
// Modified to include multiple possible nucleation events in cells that melt and solidify multiple times
void NucleiInit(int layernumber, double RNGSeed, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nx,
                int ny, int nzActive, int ZBound_Low, int id, double NMax, double dTN, double dTsigma, double deltax,
                ViewI &NucleiLocation, ViewI_H &NucleationTimes_Host, ViewI &NucleiGrainID, ViewI CellType,
                ViewI CritTimeStep, ViewF UndercoolingChange, ViewI LayerID, int &PossibleNuclei_ThisRankThisLayer,
                int &Nuclei_WholeDomain, bool AtNorthBoundary, bool AtSouthBoundary, bool AtEastBoundary,
                bool AtWestBoundary, bool RemeltingYN, int &NucleationCounter, ViewI &MaxSolidificationEvents,
                ViewI NumberOfSolidificationEvents, ViewF3D LayerTimeTempHistory) {

    // TODO: convert this subroutine into kokkos kernels, rather than copying data back to the host, and nucleation data
    // back to the device again. This is currently performed on the device due to heavy usage of standard library
    // algorithm functions

    // Copy temperature data into temporary host views for this subroutine
    ViewI_H LayerID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), LayerID);
    ViewI_H CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
    ViewI_H CritTimeStep_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CritTimeStep);
    ViewF_H UndercoolingChange_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), UndercoolingChange);
    ViewI_H MaxSolidificationEvents_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), MaxSolidificationEvents);
    ViewI_H NumberOfSolidificationEvents_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), NumberOfSolidificationEvents);
    ViewF3D_H LayerTimeTempHistory_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), LayerTimeTempHistory);

    // Three counters tracked here:
    // Nuclei_WholeDomain - tracks all nuclei (regardless of whether an event would be possible based on the layer ID
    // and cell type), used for Grain ID assignment to ensure that no Grain ID get reused - same on all MPI ranks
    // Nuclei_ThisLayer - the same as Nuclei_WholeDomain, except that it is recalculated each layer - same on all MPI
    // ranks PossibleNuclei_ThisRankThisLayer - the subset of Nuclei_ThisLayer that are located within the bounds of a
    // given MPI rank that may possibly occur (nuclei locations are associated with a liquid cell with a layer ID that
    // matches this layer number). Starts at 0 each layer PossibleNuclei_AllRanksThisLayer - the subset of
    // Nuclei_ThisLayer that may possibly occur in the layer (nuclei locations are associated with a liquid cell with a
    // layer ID that matches this layer number). Starts at 0 each layer
    PossibleNuclei_ThisRankThisLayer = 0;
    if (layernumber == 0)
        Nuclei_WholeDomain = 0;

    // Probability that a given liquid site will be a potential nucleus location
    double BulkProb = NMax * deltax * deltax * deltax;

    // Use new RNG seed for each layer
    std::mt19937_64 generator(RNGSeed + layernumber);
    // Uniform distribution for nuclei location assignment
    std::uniform_real_distribution<double> Xdist(-0.49999, nx - 0.5);
    std::uniform_real_distribution<double> Ydist(-0.49999, ny - 0.5);
    std::uniform_real_distribution<double> Zdist(-0.49999, nzActive - 0.5);
    // Gaussian distribution of nucleation undercooling
    std::normal_distribution<double> Gdistribution(dTN, dTsigma);

    // Max number of nucleated grains in this layer
    int Nuclei_ThisLayerSingle = BulkProb * (nx * ny * nzActive); // equivalent to Nuclei_ThisLayer if no remelting
    int Nuclei_ThisLayer = Nuclei_ThisLayerSingle * MaxSolidificationEvents_Host(layernumber);
    // Temporary vectors for storing nucleated grain IDs and undercooling values
    // Nuclei Grain ID are assigned to avoid reusing values from previous layers
    std::vector<int> NucleiGrainID_WholeDomain_V(Nuclei_ThisLayer);
    std::vector<double> NucleiUndercooling_WholeDomain_V(Nuclei_ThisLayer);
    // Views for storing potential nucleated grain coordinates
    ViewI_H NucleiX(Kokkos::ViewAllocateWithoutInitializing("NucleiX"), Nuclei_ThisLayer);
    ViewI_H NucleiY(Kokkos::ViewAllocateWithoutInitializing("NucleiY"), Nuclei_ThisLayer);
    ViewI_H NucleiZ(Kokkos::ViewAllocateWithoutInitializing("NucleiZ"), Nuclei_ThisLayer);
    for (int meltevent = 0; meltevent < MaxSolidificationEvents_Host(layernumber); meltevent++) {
        for (int n = 0; n < Nuclei_ThisLayerSingle; n++) {
            int NEvent = meltevent * Nuclei_ThisLayerSingle + n;
            // Generate possible nuclei locations
            double NucleiX_unrounded = Xdist(generator);
            double NucleiY_unrounded = Ydist(generator);
            double NucleiZ_unrounded = Zdist(generator);
            // Round these coordinates so they're associated with a specific cell on the grid
            NucleiX(NEvent) = std::round(NucleiX_unrounded);
            NucleiY(NEvent) = std::round(NucleiY_unrounded);
            NucleiZ(NEvent) = std::round(NucleiZ_unrounded);
            // Assign each nuclei a Grain ID (negative values used for nucleated grains) and an undercooling
            NucleiGrainID_WholeDomain_V[NEvent] = -(Nuclei_WholeDomain + NEvent + 1); // avoid using grain ID 0
            NucleiUndercooling_WholeDomain_V[NEvent] = Gdistribution(generator);
        }
    }
    // Shuffle these vectors to make sure the same grain IDs and undercooling don't end up in the same spots each layer
    std::shuffle(NucleiGrainID_WholeDomain_V.begin(), NucleiGrainID_WholeDomain_V.end(), generator);
    std::shuffle(NucleiUndercooling_WholeDomain_V.begin(), NucleiUndercooling_WholeDomain_V.end(), generator);

    if ((id == 0) && (Nuclei_ThisLayer > 0))
        std::cout << "Range of Grain IDs from which layer " << layernumber
                  << " nucleation events were selected: " << Nuclei_WholeDomain + 1 << " through "
                  << Nuclei_WholeDomain + Nuclei_ThisLayer << std::endl;
    // Update number of nuclei counter for whole domain based on the number of nuclei in this layer
    Nuclei_WholeDomain += Nuclei_ThisLayer;

    // Loop through nuclei for this layer - each MPI rank storing the nucleation events that are possible (i.e,
    // nucleation event is associated with a CA cell on that MPI rank's subdomain, the cell is liquid type, and the cell
    // is associated with the current layer of the multilayer problem) Don't put nuclei in "ghost" cells - those
    // nucleation events occur on other ranks and the existing halo exchange functionality will handle this
    std::vector<int> NucleiGrainID_MyRank_V(Nuclei_ThisLayer), NucleiLocation_MyRank_V(Nuclei_ThisLayer);
    std::vector<double> NucleationTimes_MyRank_V(Nuclei_ThisLayer);
    for (int meltevent = 0; meltevent < MaxSolidificationEvents_Host(layernumber); meltevent++) {
        for (int n = 0; n < Nuclei_ThisLayerSingle; n++) {
            int NEvent = meltevent * Nuclei_ThisLayerSingle + n;
            if (((NucleiX(NEvent) > MyXOffset) || (AtWestBoundary)) &&
                ((NucleiX(NEvent) < MyXOffset + MyXSlices - 1) || (AtEastBoundary)) &&
                ((NucleiY(NEvent) > MyYOffset) || (AtSouthBoundary)) &&
                ((NucleiY(NEvent) < MyYOffset + MyYSlices - 1) || (AtNorthBoundary))) {
                // Convert 3D location (using global X and Y coordinates) into a 1D location (using local X and Y
                // coordinates) for the possible nucleation event, both as relative to the bottom of this layer as well
                // as relative to the bottom of the overall domain
                int NucleiLocation_ThisLayer = NucleiZ[NEvent] * MyXSlices * MyYSlices +
                                               (NucleiX[NEvent] - MyXOffset) * MyYSlices +
                                               (NucleiY[NEvent] - MyYOffset);
                int NucleiLocation_AllLayers = (NucleiZ(NEvent) + ZBound_Low) * MyXSlices * MyYSlices +
                                               (NucleiX(NEvent) - MyXOffset) * MyYSlices +
                                               (NucleiY(NEvent) - MyYOffset);
                // Criteria for placing a nucleus is different with and without remelting due to problem initialization
                // differences - with remelting, a cell can up to as many nucleation sites as it it goes above/below the
                // liquidus
                if (RemeltingYN) {
                    if (meltevent < NumberOfSolidificationEvents_Host(NucleiLocation_ThisLayer)) {
                        // Nucleation event is possible - cell undergoes solidification at least once, this nucleation
                        // event is associated with one of the time periods during which the associated cell undergoes
                        // solidification
                        NucleiLocation_MyRank_V[PossibleNuclei_ThisRankThisLayer] = NucleiLocation_AllLayers;

                        int CritTimeStep_ThisEvent = LayerTimeTempHistory_Host(NucleiLocation_ThisLayer, meltevent, 1);
                        float UndercoolingChange_ThisEvent =
                            LayerTimeTempHistory_Host(NucleiLocation_ThisLayer, meltevent, 2);
                        int TimeToNucUnd = CritTimeStep_ThisEvent +
                                           round(NucleiUndercooling_WholeDomain_V[n] / UndercoolingChange_ThisEvent);
                        NucleationTimes_MyRank_V[PossibleNuclei_ThisRankThisLayer] =
                            std::max(CritTimeStep_ThisEvent, TimeToNucUnd);
                        // Assign this cell the potential nucleated grain ID
                        NucleiGrainID_MyRank_V[PossibleNuclei_ThisRankThisLayer] = NucleiGrainID_WholeDomain_V[NEvent];
                        // Increment counter on this MPI rank
                        PossibleNuclei_ThisRankThisLayer++;
                    }
                }
                else {
                    if ((CellType_Host(NucleiLocation_AllLayers) == Liquid) &&
                        (LayerID_Host(NucleiLocation_AllLayers) == layernumber)) {
                        // Nucleation event is possible - cell is liquid and associated with this layer of the
                        // multilayer problem - Add nuclei location and undercooling to the list for this rank
                        NucleiLocation_MyRank_V[PossibleNuclei_ThisRankThisLayer] = NucleiLocation_AllLayers;
                        int TimeToNucUnd = CritTimeStep_Host(NucleiLocation_AllLayers) +
                                           round(NucleiUndercooling_WholeDomain_V[NEvent] /
                                                 UndercoolingChange_Host(NucleiLocation_AllLayers));
                        NucleationTimes_MyRank_V[PossibleNuclei_ThisRankThisLayer] =
                            std::max(CritTimeStep_Host(NucleiLocation_AllLayers), TimeToNucUnd);
                        // Assign this cell the potential nucleated grain ID
                        NucleiGrainID_MyRank_V[PossibleNuclei_ThisRankThisLayer] = NucleiGrainID_WholeDomain_V[n];
                        // Increment counter on this MPI rank
                        PossibleNuclei_ThisRankThisLayer++;
                    }
                }
            }
        }
    }

    // How many nucleation events are actually possible (associated with a cell in this layer that will undergo
    // solidification)?
    int PossibleNuclei_AllRanksThisLayer;
    MPI_Reduce(&PossibleNuclei_ThisRankThisLayer, &PossibleNuclei_AllRanksThisLayer, 1, MPI_INT, MPI_SUM, 0,
               MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Number of potential nucleation events in layer " << layernumber << " : "
                  << PossibleNuclei_AllRanksThisLayer << std::endl;

    // Now that the number of nucleation events on each rank is known, resize these vectors
    NucleiLocation_MyRank_V.resize(PossibleNuclei_ThisRankThisLayer);
    NucleiGrainID_MyRank_V.resize(PossibleNuclei_ThisRankThisLayer);
    NucleationTimes_MyRank_V.resize(PossibleNuclei_ThisRankThisLayer);

    // Sort the list of time steps at which nucleation occurs, keeping the time steps paired with the corresponding
    // locations for nucleation events and grain IDs
    std::vector<std::tuple<int, int, int>> NucleationTimeLocID;
    NucleationTimeLocID.reserve(PossibleNuclei_ThisRankThisLayer);
    for (int n = 0; n < PossibleNuclei_ThisRankThisLayer; n++) {
        NucleationTimeLocID.push_back(
            std::make_tuple(NucleationTimes_MyRank_V[n], NucleiLocation_MyRank_V[n], NucleiGrainID_MyRank_V[n]));
    }
    // Sorting from low to high
    std::sort(NucleationTimeLocID.begin(), NucleationTimeLocID.end());

    // With PossibleNuclei_ThisRankThisLayer now known, resize views appropriately
    // Resize nucleation views now that PossibleNuclei_ThisRank is known for all MPI ranks
    Kokkos::resize(NucleationTimes_Host, PossibleNuclei_ThisRankThisLayer);
    Kokkos::resize(NucleiLocation, PossibleNuclei_ThisRankThisLayer);
    Kokkos::resize(NucleiGrainID, PossibleNuclei_ThisRankThisLayer);

    // Create temporary view to store nucleation locations, grain ID data initialized on the host
    // NucleationTimes_H are stored using a host view that is passed to Nucleation subroutine and later used - don't
    // need a temporary host view
    ViewI_H NucleiLocation_Host(Kokkos::ViewAllocateWithoutInitializing("NucleiLocation_Host"),
                                PossibleNuclei_ThisRankThisLayer);
    ViewI_H NucleiGrainID_Host(Kokkos::ViewAllocateWithoutInitializing("NucleiGrainID_Host"),
                               PossibleNuclei_ThisRankThisLayer);
    for (int n = 0; n < PossibleNuclei_ThisRankThisLayer; n++) {
        NucleationTimes_Host(n) = std::get<0>(NucleationTimeLocID[n]);
        NucleiLocation_Host(n) = std::get<1>(NucleationTimeLocID[n]);
        NucleiGrainID_Host(n) = std::get<2>(NucleationTimeLocID[n]);
    }
    // Copy nucleation data to the device
    NucleiLocation = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), NucleiLocation_Host);
    NucleiGrainID = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), NucleiGrainID_Host);

    // Initialize counter for the layer to 0
    NucleationCounter = 0;

    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Nucleation data structures for layer " << layernumber << " initialized" << std::endl;
}

//*****************************************************************************/
void DomainShiftAndResize(int id, int MyXSlices, int MyYSlices, int &ZShift, int &ZBound_Low, int &ZBound_High,
                          int &nzActive, int LocalDomainSize, int &LocalActiveDomainSize, int &BufSizeZ,
                          int LayerHeight, ViewI CellType, int layernumber, ViewI LayerID) {

    int ZBound_LowOld = ZBound_Low;

    // The top "top" of the active domain is a shift of "LayerHeight" from the previous domain top
    ZBound_High += LayerHeight;

    // The new "bottom" of the active domain is located at the Z coordinate of the lowest active cells remaining in the
    // domain
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
void ZeroResetViews(int LocalActiveDomainSize, int BufSizeX, int BufSizeY, int BufSizeZ, ViewF &DiagonalLength,
                    ViewF &CritDiagonalLength, ViewF &DOCenter, int DecompositionStrategy, Buffer2D &BufferWestSend,
                    Buffer2D &BufferEastSend, Buffer2D &BufferNorthSend, Buffer2D &BufferSouthSend,
                    Buffer2D &BufferNorthEastSend, Buffer2D &BufferNorthWestSend, Buffer2D &BufferSouthEastSend,
                    Buffer2D &BufferSouthWestSend, Buffer2D &BufferWestRecv, Buffer2D &BufferEastRecv,
                    Buffer2D &BufferNorthRecv, Buffer2D &BufferSouthRecv, Buffer2D &BufferNorthEastRecv,
                    Buffer2D &BufferNorthWestRecv, Buffer2D &BufferSouthEastRecv, Buffer2D &BufferSouthWestRecv,
                    ViewI &SteeringVector) {

    // Resize active cell data structure and halo regions on device
    // Resize steering vector as LocalActiveDomainSize may have changed
    Kokkos::resize(SteeringVector, LocalActiveDomainSize);

    // Resize active cell data structures - host and device
    Kokkos::resize(DiagonalLength, LocalActiveDomainSize);
    Kokkos::resize(DOCenter, 3 * LocalActiveDomainSize);
    Kokkos::resize(CritDiagonalLength, 26 * LocalActiveDomainSize);
    Kokkos::resize(BufferNorthSend, BufSizeX * BufSizeZ, 5);
    Kokkos::resize(BufferSouthSend, BufSizeX * BufSizeZ, 5);
    Kokkos::resize(BufferEastSend, BufSizeY * BufSizeZ, 5);
    Kokkos::resize(BufferWestSend, BufSizeY * BufSizeZ, 5);
    Kokkos::resize(BufferNorthEastSend, BufSizeZ, 5);
    Kokkos::resize(BufferNorthWestSend, BufSizeZ, 5);
    Kokkos::resize(BufferSouthEastSend, BufSizeZ, 5);
    Kokkos::resize(BufferSouthWestSend, BufSizeZ, 5);

    Kokkos::resize(BufferNorthRecv, BufSizeX * BufSizeZ, 5);
    Kokkos::resize(BufferSouthRecv, BufSizeX * BufSizeZ, 5);
    Kokkos::resize(BufferEastRecv, BufSizeY * BufSizeZ, 5);
    Kokkos::resize(BufferWestRecv, BufSizeY * BufSizeZ, 5);
    Kokkos::resize(BufferNorthEastRecv, BufSizeZ, 5);
    Kokkos::resize(BufferNorthWestRecv, BufSizeZ, 5);
    Kokkos::resize(BufferSouthEastRecv, BufSizeZ, 5);
    Kokkos::resize(BufferSouthWestRecv, BufSizeZ, 5);

    // Reset active cell data structures on device
    Kokkos::deep_copy(DiagonalLength, 0);
    Kokkos::deep_copy(DOCenter, 0);
    Kokkos::deep_copy(CritDiagonalLength, 0);
    // Reset halo region structures on device
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
}
