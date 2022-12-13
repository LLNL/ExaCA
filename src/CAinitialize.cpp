// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAinitialize.hpp"

#include "CAconfig.hpp"
#include "CAfunctions.hpp"
#include "CAghostnodes.hpp"
#include "CAparsefiles.hpp"
#include "CAupdate.hpp"

#include "mpi.h"

#ifdef ExaCA_ENABLE_JSON
#include <nlohmann/json.hpp>
#endif

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <regex>

// Initializes input parameters, mesh, temperature field, and grain structures for CA simulations
//*****************************************************************************/
// Read ExaCA input file.
void InputReadFromFile_Old(int id, std::string InputFile, std::string &SimulationType, double &deltax, double &NMax,
                           double &dTN, double &dTsigma, std::string &OutputFile, std::string &GrainOrientationFile,
                           int &TempFilesInSeries, std::vector<std::string> &temp_paths, double &HT_deltax,
                           bool &RemeltingYN, double &deltat, int &NumberOfLayers, int &LayerHeight,
                           std::string &MaterialFileName, std::string &SubstrateFileName, float &SubstrateGrainSpacing,
                           bool &UseSubstrateFile, double &G, double &R, int &nx, int &ny, int &nz,
                           double &FractSurfaceSitesActive, std::string &PathToOutput, int &PrintDebug,
                           bool &PrintMisorientation, bool &PrintFinalUndercoolingVals, bool &PrintFullOutput,
                           int &NSpotsX, int &NSpotsY, int &SpotOffset, int &SpotRadius, bool &PrintTimeSeries,
                           int &TimeSeriesInc, bool &PrintIdleTimeSeriesFrames, bool &PrintDefaultRVE, double &RNGSeed,
                           bool &BaseplateThroughPowder, double &PowderActiveFraction, int &RVESize,
                           bool &LayerwiseTempRead, bool &PrintBinary) {

    // Required inputs that should be present in the input file, regardless of problem type
    std::vector<std::string> RequiredInputs_General = {
        "Material",                                      // Required input 0
        "Cell size",                                     // Required input 1
        "Heterogeneous nucleation density",              // Required input 2
        "Mean nucleation undercooling",                  // Required input 3
        "Standard deviation of nucleation undercooling", // Required input 4
        "Path to output",                                // Required input 5
        "Output file base name",                         // Required input 6
        "File of grain orientations",                    // Required input 7
        "grain misorientation values",                   // Required input 8
        "all ExaCA data"                                 // Required input 9
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
        "Print vtk data as binary",                     // Optional input 7
    };
    std::vector<std::string> DeprecatedInputs_General = {
        "Decomposition strategy", // Deprecated input 0
    };

    // Values used temporarily to store information from the file
    int NRatio = 0;
    float SpotOffsetFloat = 0.0;
    float SpotRadiusFloat = 0.0;
    double TimeSeriesFrameInc_time = 0.0;
    std::string MaterialName, GrainOrientationFile_Read;

    // Read the input file provided on the command line - determine which problem type is being simulated
    checkFileExists(InputFile, id);
    std::ifstream InputData;
    InputData.open(InputFile);
    skipLines(InputData, "*****"); // Skip lines until past the header lines above the asterisks
    // First line after the header must be the problem type: either R, C, or S
    // An "M" after the problem type ("SM" or "RM") indicates that the problem uses remelting logic
    // Additional required/optional inputs depending on problem type
    SimulationType = parseInput(InputData, "Problem type");
    std::vector<std::string> RequiredInputs_ProblemSpecific, OptionalInputs_ProblemSpecific,
        DeprecatedInputs_ProblemSpecific;
    std::vector<bool> DeprecatedInputs_RequiredYN;
    // Currently no deprecated inputs. Keeping variables in anticipation of more in the future.
    DeprecatedInputs_ProblemSpecific.resize(0);
    // If using deprecated inputs, which ones are required/optional
    DeprecatedInputs_RequiredYN.resize(0);
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
        OptionalInputs_ProblemSpecific.resize(4);
        OptionalInputs_ProblemSpecific[0] = "Substrate grain spacing";
        OptionalInputs_ProblemSpecific[1] = "Substrate filename";
        OptionalInputs_ProblemSpecific[2] = "Extend baseplate through layers";
        OptionalInputs_ProblemSpecific[3] = "Density of powder surface sites active";
    }
    else if (SimulationType == "R") {
        RequiredInputs_ProblemSpecific.resize(2);
        RequiredInputs_ProblemSpecific[0] = "Time step";
        RequiredInputs_ProblemSpecific[1] = "Path to and name of temperature field assembly instructions";
        OptionalInputs_ProblemSpecific.resize(6);
        OptionalInputs_ProblemSpecific[0] = "Substrate grain spacing";
        OptionalInputs_ProblemSpecific[1] = "Substrate filename";
        OptionalInputs_ProblemSpecific[2] = "default RVE output";
        OptionalInputs_ProblemSpecific[3] = "Extend baseplate through layers";
        OptionalInputs_ProblemSpecific[4] = "Density of powder surface sites active";
        OptionalInputs_ProblemSpecific[5] = "Default RVE size, in CA cells";
    }
    else {
        std::string error = "Error: problem type must be C, S, or R: the value given was " + SimulationType;
        throw std::runtime_error(error);
    }
    int NumRequiredInputs_General = RequiredInputs_General.size();
    int NumOptionalInputs_General = OptionalInputs_General.size();
    int NumRequiredInputs_ProblemSpecific = RequiredInputs_ProblemSpecific.size();
    int NumOptionalInputs_ProblemSpecific = OptionalInputs_ProblemSpecific.size();
    int NumDeprecatedInputs_ProblemSpecific = DeprecatedInputs_ProblemSpecific.size();
    int NumDeprecatedInputs_General = DeprecatedInputs_General.size();
    std::vector<std::string> RequiredInputsRead_General(NumRequiredInputs_General),
        OptionalInputsRead_General(NumOptionalInputs_General),
        RequiredInputsRead_ProblemSpecific(NumRequiredInputs_ProblemSpecific),
        OptionalInputsRead_ProblemSpecific(NumOptionalInputs_ProblemSpecific),
        DeprecatedInputsRead_ProblemSpecific(NumDeprecatedInputs_ProblemSpecific),
        DeprecatedInputsRead_General(NumDeprecatedInputs_General);

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
                    if (!(OptionalYN)) {
                        // If not a general/problem type specific optional input, check against deprecated inputs
                        // (general and problem specific)
                        bool DeprecatedYN = parseInputFromList(line, DeprecatedInputs_ProblemSpecific,
                                                               DeprecatedInputsRead_ProblemSpecific,
                                                               NumDeprecatedInputs_ProblemSpecific);
                        if (!(DeprecatedYN))
                            DeprecatedYN =
                                parseInputFromList(line, DeprecatedInputs_General, DeprecatedInputsRead_General,
                                                   NumDeprecatedInputs_General);
                        if (!(DeprecatedYN) && (id == 0)) {
                            std::cout << "WARNING: input " << line
                                      << " did not match any optional, required, nor deprecated inputs known to ExaCA "
                                         "and will be ignored"
                                      << std::endl;
                        }
                    }
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
    MaterialName = RequiredInputsRead_General[0];
    deltax = getInputDouble(RequiredInputsRead_General[1], -6);
    NMax = getInputDouble(RequiredInputsRead_General[2], 12);
    dTN = getInputDouble(RequiredInputsRead_General[3]);
    dTsigma = getInputDouble(RequiredInputsRead_General[4]);
    PathToOutput = RequiredInputsRead_General[5];
    OutputFile = RequiredInputsRead_General[6];
    GrainOrientationFile_Read = RequiredInputsRead_General[7];
    PrintMisorientation = getInputBool(RequiredInputsRead_General[8]);
    PrintFullOutput = getInputBool(RequiredInputsRead_General[9]);
    // Problem type-specific inputs
    if (SimulationType == "R") {
        deltat = getInputDouble(RequiredInputsRead_ProblemSpecific[0], -6);
        // If not using deprecated inputs, open specified file and parse data into "temp_paths" vector, also obtaining
        // TempFilesInSeries, NumberOfLayers, LayerHeight, HT_deltax
        std::string TemperatureFieldInstructions = RequiredInputsRead_ProblemSpecific[1];
        if (!(RequiredInputsRead_ProblemSpecific[1].empty()))
            parseTInstuctionsFile(id, TemperatureFieldInstructions, TempFilesInSeries, NumberOfLayers, LayerHeight,
                                  deltax, HT_deltax, temp_paths, RemeltingYN, LayerwiseTempRead);
        if ((OptionalInputsRead_ProblemSpecific[3].empty()))
            BaseplateThroughPowder = false; // defaults to using baseplate only for layer 0 substrate
        else
            BaseplateThroughPowder = getInputBool(OptionalInputsRead_ProblemSpecific[3]);
        if ((OptionalInputsRead_ProblemSpecific[4].empty()))
            PowderActiveFraction = 1.0; // defaults to a unique grain at each site in the powder layers
        else {
            PowderActiveFraction =
                getInputDouble(OptionalInputsRead_ProblemSpecific[4], 12) *
                pow(deltax, 3); // powder density is given as a density per unit volume, normalized by 10^12 m^-3 -->
                                // convert this into a fraction of sites active on the CA grid (0 to 1)
            if ((PowderActiveFraction < 0.0) || (PowderActiveFraction > 1.0))
                throw std::runtime_error("Error: Density of powder surface sites active must be larger than 0 and less "
                                         "than 1/(CA cell volume)");
            if (BaseplateThroughPowder)
                throw std::runtime_error("Error: if the option to extend the baseplate through through all powder "
                                         "layers is turned on, a powder density cannot be given");
        }
        if (id == 0) {
            std::cout << "CA Simulation using temperature data from file(s)" << std::endl;
            std::cout << "The time step is " << deltat << " seconds" << std::endl;
            std::cout << "The first temperature data file to be read is " << temp_paths[0] << ", and there are "
                      << TempFilesInSeries << " in the series" << std::endl;
            std::cout << "A total of " << NumberOfLayers << " layers of solidification offset by " << LayerHeight
                      << " CA cells will be simulated" << std::endl;
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
        if ((OptionalInputsRead_ProblemSpecific[2].empty()))
            BaseplateThroughPowder = false; // defaults to using baseplate only for layer 0 substrate
        else
            BaseplateThroughPowder = getInputBool(OptionalInputsRead_ProblemSpecific[2]);
        if ((OptionalInputsRead_ProblemSpecific[3].empty()))
            PowderActiveFraction = 1.0; // defaults to a unique grain at each site in the powder layers
        else {
            PowderActiveFraction =
                getInputDouble(OptionalInputsRead_ProblemSpecific[3], 12) *
                pow(deltax, 3); // powder density is given as a density per unit volume, normalized by 10^12 m^-3 -->
                                // convert this into a density of sites active on the CA grid (0 to 1)
            if ((PowderActiveFraction < 0.0) || (PowderActiveFraction > 1.0))
                throw std::runtime_error("Error: Density of powder surface sites active must be larger than 0 and less "
                                         "than 1/(CA cell volume)");
            if (BaseplateThroughPowder)
                throw std::runtime_error("Error: if the option to extend the baseplate through the powder layers it "
                                         "toggled, a powder density cannot be given");
        }
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
    // Should the vtk output files be printed as binary data? (By default, they are printed as ASCII data)
    if (OptionalInputsRead_General[7].empty())
        PrintBinary = false;
    else
        PrintBinary = getInputBool(OptionalInputsRead_General[7]);
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
                checkFileExists(SubstrateFileName, id);
            }
        }
    }
    if (SimulationType == "R") {

        // Check that cell size/temperature data resolution are equivalent for simulations with remelting
        if ((RemeltingYN) && (HT_deltax != deltax)) {
            throw std::runtime_error("Error: For simulations with external temperature data and remelting logic, CA "
                                     "cell size and input temperature data resolution must be equivalent");
        }
        // Check that temperature file(s) exist
        for (int i = 0; i < TempFilesInSeries; i++) {
            if (id == 0)
                std::cout << "Checking file " << temp_paths[i] << std::endl;
            checkFileExists(temp_paths[i], id);
        }

        // Should optional RVE data be printed for the standard location (center of domain in X and Y, as close to the
        // top of the domain in Z as possiblewithout including the last layer's microstructure)?
        PrintDefaultRVE = false;
        if (!(OptionalInputsRead_ProblemSpecific[2].empty()))
            PrintDefaultRVE = getInputBool(OptionalInputsRead_ProblemSpecific[2]);
        // If printing RVE data, specify RVE size (in cells per side). Otherwise, RVE is 0.5 mm per side
        if (!(OptionalInputsRead_ProblemSpecific[5].empty()))
            RVESize = getInputInt(OptionalInputsRead_ProblemSpecific[5]);
        else
            RVESize = 0.0005 / deltax;
        if ((id == 0) && (PrintDefaultRVE))
            std::cout << "RVE data from a default location in the simulation will be printed, consisting of " << RVESize
                      << " cells per side" << std::endl;
    }
    else {
        // RVE data print option is only for simulation type R
        PrintDefaultRVE = false;
    }

    // Path to file of materials constants based on install/source location
    MaterialFileName = checkFileInstalled(MaterialName, id);
    checkFileNotEmpty(MaterialFileName);

    // Path to file of grain orientations based on install/source location
    GrainOrientationFile = checkFileInstalled(GrainOrientationFile_Read, id);
    checkFileNotEmpty(GrainOrientationFile);

    // No printing of debug data if remelting is used - not currently supported
    if ((RemeltingYN) && (PrintDebug > 0)) {
        std::cout << "Printing of views following initialization with use of remelting is not currently supported"
                  << std::endl;
        PrintDebug = 0;
    }
    if (id == 0) {
        std::cout << "Decomposition Strategy is 1D, with the domain partitioned in the Y direction" << std::endl;
        std::cout << "Material simulated is " << MaterialName << std::endl;
        std::cout << "CA cell size is " << deltax * pow(10, 6) << " microns" << std::endl;
        std::cout << "Nucleation density is " << NMax << " per m^3" << std::endl;
        std::cout << "Mean nucleation undercooling is " << dTN << " K, standard deviation of distribution is "
                  << dTsigma << "K" << std::endl;
        if (PrintTimeSeries)
            std::cout << "Intermediate output for movie frames will be printed every " << TimeSeriesInc
                      << " time steps (or every " << TimeSeriesFrameInc_time << " microseconds)" << std::endl;
        if (RemeltingYN)
            std::cout << "This simulation includes logic for cells melting and multiple solidification events"
                      << std::endl;
    }
}

// Read ExaCA input file - new (json) file format
#ifdef ExaCA_ENABLE_JSON
void InputReadFromFile(int id, std::string InputFile, std::string &SimulationType, double &deltax, double &NMax,
                       double &dTN, double &dTsigma, std::string &OutputFile, std::string &GrainOrientationFile,
                       int &TempFilesInSeries, std::vector<std::string> &temp_paths, double &HT_deltax,
                       bool &RemeltingYN, double &deltat, int &NumberOfLayers, int &LayerHeight,
                       std::string &MaterialFileName, std::string &SubstrateFileName, float &SubstrateGrainSpacing,
                       bool &UseSubstrateFile, double &G, double &R, int &nx, int &ny, int &nz,
                       double &FractSurfaceSitesActive, std::string &PathToOutput, int &PrintDebug,
                       bool &PrintMisorientation, bool &PrintFinalUndercoolingVals, bool &PrintFullOutput, int &NSpotsX,
                       int &NSpotsY, int &SpotOffset, int &SpotRadius, bool &PrintTimeSeries, int &TimeSeriesInc,
                       bool &PrintIdleTimeSeriesFrames, bool &PrintDefaultRVE, double &RNGSeed,
                       bool &BaseplateThroughPowder, double &PowderActiveFraction, int &RVESize,
                       bool &LayerwiseTempRead, bool &PrintBinary) {

    std::ifstream InputData(InputFile);
    nlohmann::json inputdata = nlohmann::json::parse(InputData);

    // General inputs
    SimulationType = inputdata["SimulationType"];
    // "C": constrained (directional) solidification
    // "S": array of overlapping hemispherical spots
    // "R": time-temperature history comes from external files
    // Check if simulation type includes remelting ("M" suffix to input problem type)
    if (SimulationType == "RM") {
        SimulationType = "R";
        RemeltingYN = true;
    }
    else if (SimulationType == "SM") {
        SimulationType = "S";
        RemeltingYN = true;
    }
    else {
        // Simulation does not including remelting logic
        RemeltingYN = false;
    }
    // Input files that should be present for all problem types
    std::string MaterialFileName_Read = inputdata["MaterialFileName"];
    std::string GrainOrientationFile_Read = inputdata["GrainOrientationFile"];
    // Path to file of materials constants based on install/source location
    MaterialFileName = checkFileInstalled(MaterialFileName_Read, id);
    checkFileNotEmpty(MaterialFileName);
    // Path to file of grain orientations based on install/source location
    GrainOrientationFile = checkFileInstalled(GrainOrientationFile_Read, id);
    checkFileNotEmpty(GrainOrientationFile);
    // Seed for random number generator (defaults to 0 if not given)
    if (inputdata.contains("RNGSeed"))
        RNGSeed = inputdata["RNGSeed"];
    else
        RNGSeed = 0.0;

    // Domain inputs:
    // Cell size - given in meters, stored in micrometers
    deltax = inputdata["Domain"]["deltax"];
    deltax = deltax * pow(10, -6);
    // Time step - given in seconds, stored in microseconds
    deltat = inputdata["Domain"]["deltat"];
    deltat = deltat * pow(10, -6);
    if (SimulationType == "C") {
        // Domain size, in cells
        nx = inputdata["Domain"]["nx"];
        ny = inputdata["Domain"]["ny"];
        nz = inputdata["Domain"]["nz"];
        NumberOfLayers = 1;
        LayerHeight = nz;
    }
    else {
        // Number of layers, layer height are needed for problem types S and R
        NumberOfLayers = inputdata["Domain"]["numberOfLayers"];
        LayerHeight = inputdata["Domain"]["layerOffset"];
        // Type S needs spot information, which is then used to compute the domain bounds
        if (SimulationType == "S") {
            NSpotsX = inputdata["Domain"]["nSpotsX"];
            NSpotsY = inputdata["Domain"]["nSpotsY"];
            // Radius and offset are given in micrometers, convert to cells
            SpotRadius = inputdata["Domain"]["rSpots"];
            SpotRadius = SpotRadius * pow(10, -6) / deltax;
            SpotOffset = inputdata["Domain"]["spotOffset"];
            SpotOffset = SpotOffset * pow(10, -6) / deltax;
            // Calculate nx, ny, and nz based on spot array pattern and number of layers
            nz = SpotRadius + 1 + (NumberOfLayers - 1) * LayerHeight;
            nx = 2 * SpotRadius + 1 + SpotOffset * (NSpotsX - 1);
            ny = 2 * SpotRadius + 1 + SpotOffset * (NSpotsY - 1);
        }
    }

    // Nucleation inputs:
    // Nucleation density (normalized by 10^12 m^-3), mean nucleation undercooling/st dev undercooling(K)
    NMax = inputdata["Nucleation"]["NMax"];
    NMax = NMax * pow(10, 12);
    dTN = inputdata["Nucleation"]["dTN"];
    dTsigma = inputdata["Nucleation"]["dTsigma"];

    // Temperature inputs:
    if (SimulationType == "R") {
        // Temperature data resolution - default to using CA cell size if the assumed temperature data resolution if not
        // given
        if (inputdata["TemperatureData"].contains("HTdeltax")) {
            HT_deltax = inputdata["TemperatureData"]["HTdeltax"];
            // Value is given in micrometers, convert to meters
            HT_deltax = HT_deltax * pow(10, -6);
        }
        else
            HT_deltax = deltax;
        if ((RemeltingYN) && (HT_deltax != deltax))
            throw std::runtime_error("Error: For simulations with external temperature data and remelting logic, CA "
                                     "cell size and input temperature data resolution must be equivalent");
        // Read all temperature files at once (default), or one at a time?
        if (inputdata["TemperatureData"].contains("layerwiseTempRead")) {
            LayerwiseTempRead = inputdata["TemperatureData"]["layerwiseTempRead"];
            if ((!(RemeltingYN)) && (LayerwiseTempRead)) {
                if (id == 0)
                    std::cout << "Warning: ability to read temperature files one at a time during initialization of "
                                 "new layers is only supported for simulations with remelting"
                              << std::endl;
                LayerwiseTempRead = false;
            }
        }
        else
            LayerwiseTempRead = false;
        // Get the paths/number of/names of the temperature data files used
        TempFilesInSeries = inputdata["TemperatureData"]["TemperatureFiles"].size();
        if (TempFilesInSeries == 0)
            throw std::runtime_error("Error: No temperature files listed in the temperature instructions file");
        else {
            for (int filename = 0; filename < TempFilesInSeries; filename++)
                temp_paths.push_back(inputdata["TemperatureData"]["TemperatureFiles"][filename]);
        }
    }
    else {
        // Temperature data uses fixed thermal gradient (K/m) and cooling rate (K/s)
        G = inputdata["TemperatureData"]["G"];
        R = inputdata["TemperatureData"]["R"];
    }

    // Substrate inputs:
    if (SimulationType == "C") {
        // Fraction of sites at bottom surface active
        FractSurfaceSitesActive = inputdata["Substrate"]["fActive"];
    }
    else {
        // Substrate data - should data come from an initial size or a file?
        if ((inputdata["Substrate"].contains("substrateFilename")) && (inputdata["Substrate"].contains("S0")))
            throw std::runtime_error("Error: only one of substrate grain size and substrate structure filename should "
                                     "be provided in the input file");
        else if (inputdata["Substrate"].contains("substrateFilename")) {
            SubstrateFileName = inputdata["Substrate"]["substrateFilename"];
            UseSubstrateFile = true;
        }
        else if (inputdata["Substrate"].contains("S0")) {
            SubstrateGrainSpacing = inputdata["Substrate"]["S0"];
            UseSubstrateFile = false;
        }
        // Should the baseplate microstructure be extended through the powder layers? Default is false
        if (inputdata["Substrate"].contains("extendSubstrateThroughPower"))
            BaseplateThroughPowder = inputdata["Substrate"]["extendSubstrateThroughPower"];
        else
            BaseplateThroughPowder = false;
        if (inputdata["Substrate"].contains("powderDensity")) {
            // powder density is given as a density per unit volume, normalized by 10^12 m^-3 --> convert this into a
            // density of sites active on the CA grid (0 to 1)
            PowderActiveFraction = inputdata["Substrate"]["powderDensity"];
            PowderActiveFraction = PowderActiveFraction * pow(10, 12) * pow(deltax, 3);
            if ((PowderActiveFraction < 0.0) || (PowderActiveFraction > 1.0))
                throw std::runtime_error("Error: Density of powder surface sites active must be larger than 0 and less "
                                         "than 1/(CA cell volume)");
            if (BaseplateThroughPowder)
                throw std::runtime_error("Error: if the option to extend the baseplate through the powder layers it "
                                         "toggled, a powder density cannot be given");
        }
        else
            PowderActiveFraction = 1.0; // defaults to a unique grain at each site in the powder layers
    }

    // Printing inputs:
    // Path to output data
    PathToOutput = inputdata["Printing"]["pathToOutput"];
    // Name of output data
    OutputFile = inputdata["Printing"]["outputFile"];
    // Should ASCII or binary be used to print vtk data? Defaults to ASCII if not given
    if (inputdata["Printing"].contains("printBinary"))
        PrintBinary = inputdata["Printing"]["printBinary"];
    else
        PrintBinary = false;
    // Should default ExaConstit output be printed after the simulation?
    if (inputdata["Printing"].contains("printExaConstitDefault")) {
        PrintDefaultRVE = inputdata["Printing"]["printExaConstitDefault"];
        // If so, what size RVE should be printed? Defaults to 0.5 by 0.5 by 0.5 mm
        if (PrintDefaultRVE) {
            if (inputdata["Printing"].contains("printExaConstitSize"))
                RVESize = inputdata["Printing"]["printExaConstitSize"];
            else
                RVESize = 0.0005 / deltax;
        }
    }
    else
        PrintDefaultRVE = false;
    // Which fields should be printed at the start and end of the simulation? Certain permutations are not currently
    // supported
    // TODO: redo to customize printing of data, allow any fields to be printed at start or end of simulation, with or
    // without remelting
    // Field names are: GrainID, LayerID, CritTimeStep, UndercoolingCurrent, UndercoolingChange, CellType,
    // GrainMisorientation
    std::vector<std::string> Fieldnames_key = {"GrainID",
                                               "LayerID",
                                               "CritTimeStep",
                                               "UndercoolingCurrent",
                                               "UndercoolingChange",
                                               "CellType",
                                               "GrainMisorientation"};
    std::vector<bool> PrintFieldsInit = getPrintFieldValues(inputdata, "PrintFieldsInit", Fieldnames_key);
    std::vector<bool> PrintFieldsFinal = getPrintFieldValues(inputdata, "PrintFieldsFinal", Fieldnames_key);
    // Printing init data: if any of fields 3, 4, 5 are present, print first six fields. If any of fields 0, 1, 2 are
    // present, print first three fields (for consistency with old input file format)
    if ((PrintFieldsInit[3]) || (PrintFieldsInit[4]) || (PrintFieldsInit[5]))
        PrintDebug = 2;
    else if ((PrintFieldsInit[0]) || (PrintFieldsInit[1]) || (PrintFieldsInit[2]))
        PrintDebug = 1;
    else
        PrintDebug = 0;
    if ((PrintFieldsInit[6]) && (id == 0))
        std::cout
            << "Warning: Option of printing grain misorientation data currently only supported at end of simulation"
            << std::endl;
    // Printing final data: "UndercoolingCurrent" at the end of the simulation is the final undercooling of each cell
    if ((PrintFieldsFinal[0]) || (PrintFieldsFinal[1]))
        PrintFullOutput = true;
    else
        PrintFullOutput = false;
    if (PrintFieldsFinal[3])
        PrintFinalUndercoolingVals = true;
    else
        PrintFinalUndercoolingVals = false;
    if (PrintFieldsFinal[6])
        PrintMisorientation = true;
    else
        PrintMisorientation = false;
    if (((PrintFieldsFinal[2]) || (PrintFieldsFinal[4]) || (PrintFieldsFinal[5])) && (id == 0))
        std::cout << "Warning: Options of printing crit time step, undercooling change, and cell type are only "
                     "supported after simulation initialization"
                  << std::endl;
    // Should intermediate output be printed?
    if (inputdata["Printing"].contains("PrintIntermediateOutput")) {
        // An increment of 0 will set the intermediate file printing to false
        double TimeSeriesFrameInc_time = inputdata["Printing"]["PrintIntermediateOutput"]["Frequency"];
        if (TimeSeriesFrameInc_time == 0) {
            PrintTimeSeries = false;
            PrintIdleTimeSeriesFrames = false;
        }
        else {
            PrintTimeSeries = true;
            // Increment is given in microseconds, convert to seconds
            TimeSeriesFrameInc_time = TimeSeriesFrameInc_time * pow(10, -6);
            TimeSeriesInc = round(TimeSeriesFrameInc_time / deltat);
            // Should the intermediate output be printed even if the simulation was unchanged from the previous output
            // step?
            PrintIdleTimeSeriesFrames = inputdata["Printing"]["PrintIntermediateOutput"]["PrintIdleFrames"];
        }
    }
    else {
        PrintTimeSeries = false;
        PrintIdleTimeSeriesFrames = false;
    }

    // Print information to console about the input file data read
    if (id == 0) {
        std::cout << "Material simulated is " << MaterialFileName << std::endl;
        std::cout << "CA cell size is " << deltax * pow(10, 6) << " microns" << std::endl;
        std::cout << "Nucleation density is " << NMax << " per m^3" << std::endl;
        std::cout << "Mean nucleation undercooling is " << dTN << " K, standard deviation of distribution is "
                  << dTsigma << "K" << std::endl;
        if (PrintTimeSeries)
            std::cout << "Intermediate output for movie frames will be printed every " << TimeSeriesInc
                      << " time steps (or every " << TimeSeriesInc * deltat << " microseconds)" << std::endl;
        if (RemeltingYN)
            std::cout << "This simulation includes logic for cells melting and multiple solidification events"
                      << std::endl;
        if (SimulationType == "C") {
            std::cout << "CA Simulation using a unidirectional, fixed thermal gradient of " << G
                      << " K/m and a cooling rate of " << R << " K/s" << std::endl;
            std::cout << "The time step is " << deltat * pow(10, 6) << " microseconds" << std::endl;
            std::cout << "The fraction of CA cells at the bottom surface that are active is " << FractSurfaceSitesActive
                      << std::endl;
        }
        else if (SimulationType == "S") {
            std::cout << "CA Simulation using a radial, fixed thermal gradient of " << G
                      << " K/m as a series of hemispherical spots, and a cooling rate of " << R << " K/s" << std::endl;
            std::cout << "A total of " << NumberOfLayers << " spots per layer, with layers offset by " << LayerHeight
                      << " CA cells will be simulated" << std::endl;
            std::cout << "The time step is " << deltat * pow(10, 6) << " microseconds" << std::endl;
        }
        else if (SimulationType == "R") {
            std::cout << "CA Simulation using temperature data from file(s)" << std::endl;
            std::cout << "The time step is " << deltat << " seconds" << std::endl;
            std::cout << "The first temperature data file to be read is " << temp_paths[0] << ", and there are "
                      << TempFilesInSeries << " in the series" << std::endl;
            std::cout << "A total of " << NumberOfLayers << " layers of solidification offset by " << LayerHeight
                      << " CA cells will be simulated" << std::endl;
        }
    }
}
#endif

void checkPowderOverflow(int nx, int ny, int LayerHeight, int NumberOfLayers, bool BaseplateThroughPowder,
                         double PowderActiveFraction) {

    // Check to make sure powder grain density is compatible with the number of powder sites
    // If this problem type includes a powder layer of some grain density, ensure that integer overflow won't occur when
    // assigning powder layer GrainIDs
    if (!(BaseplateThroughPowder)) {
        long int NumCellsPowderLayers =
            (long int)(nx) * (long int)(ny) * (long int)(LayerHeight) * (long int)(NumberOfLayers - 1);
        long int NumAssignedCellsPowderLayers =
            std::lround(round(static_cast<double>(NumCellsPowderLayers) * PowderActiveFraction));
        if (NumAssignedCellsPowderLayers > INT_MAX)
            throw std::runtime_error("Error: A smaller value for powder density is required to avoid potential integer "
                                     "overflow when assigning powder layer GrainID");
    }
}

//*****************************************************************************/
// Intialize neighbor list structures (NeighborX, NeighborY, NeighborZ)
void NeighborListInit(NList &NeighborX, NList &NeighborY, NList &NeighborZ) {

    // Assignment of neighbors around a cell "X" is as follows (in order of closest to furthest from cell "X")
    // Neighbors 0 through 8 are in the -Y direction
    // Neighbors 9 through 16 are in the XY plane with cell X
    // Neighbors 17 through 25 are in the +Y direction
    NeighborX = {0, 1, -1, 0, 0, -1, 1, -1, 1, 0, 0, 1, -1, 1, -1, 1, -1, 0, 1, -1, 0, 0, 1, -1, 1, -1};
    NeighborY = {-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    NeighborZ = {0, 0, 0, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 0, 0, 0, 0, 0, 1, -1, 1, 1, -1, -1};
}

// Check if the temperature data is in ASCII or binary format
bool checkTemperatureFileFormat(std::string tempfile_thislayer) {
    bool BinaryInputData;
    std::size_t found = tempfile_thislayer.find(".catemp");
    if (found == std::string::npos)
        BinaryInputData = false;
    else
        BinaryInputData = true;
    return BinaryInputData;
}

// Obtain the physical XYZ bounds of the domain, using either domain size from the input file, or reading temperature
// data files and parsing the coordinates
void FindXYZBounds(std::string SimulationType, int id, double &deltax, int &nx, int &ny, int &nz,
                   std::vector<std::string> &temp_paths, double &XMin, double &XMax, double &YMin, double &YMax,
                   double &ZMin, double &ZMax, int &LayerHeight, int NumberOfLayers, int TempFilesInSeries,
                   double *ZMinLayer, double *ZMaxLayer, int SpotRadius) {

    if (SimulationType == "R") {
        // Two passes through reading temperature data files- the first pass only reads the headers to
        // determine units and X/Y/Z bounds of the simulaton domain. Using the X/Y/Z bounds of the simulation domain,
        // nx, ny, and nz can be calculated and the domain decomposed among MPI processes. The maximum number of
        // remelting events in the simulation can also be calculated. The second pass reads the actual X/Y/Z/liquidus
        // time/cooling rate data and each rank stores the data relevant to itself in "RawData" - this is done in the
        // subroutine "ReadTemperatureData"
        XMin = std::numeric_limits<double>::max();
        YMin = std::numeric_limits<double>::max();
        ZMin = std::numeric_limits<double>::max();
        XMax = std::numeric_limits<double>::lowest();
        YMax = std::numeric_limits<double>::lowest();
        ZMax = std::numeric_limits<double>::lowest();

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
            // Get min and max x coordinates in this file, which can be a binary or ASCII input file
            // binary file type uses extension .catemp, all other file types assumed to be comma-separated ASCII input
            bool BinaryInputData = checkTemperatureFileFormat(tempfile_thislayer);
            // { Xmin, Xmax, Ymin, Ymax, Zmin, Zmax }
            std::array<double, 6> XYZMinMax_ThisLayer =
                parseTemperatureCoordinateMinMax(tempfile_thislayer, BinaryInputData);

            // Based on the input file's layer offset, adjust ZMin/ZMax from the temperature data coordinate
            // system to the multilayer CA coordinate system Check to see in the XYZ bounds for this layer are
            // also limiting for the entire multilayer CA coordinate system
            XYZMinMax_ThisLayer[4] += deltax * LayerHeight * (LayerReadCount - 1);
            XYZMinMax_ThisLayer[5] += deltax * LayerHeight * (LayerReadCount - 1);
            if (XYZMinMax_ThisLayer[0] < XMin)
                XMin = XYZMinMax_ThisLayer[0];
            if (XYZMinMax_ThisLayer[1] > XMax)
                XMax = XYZMinMax_ThisLayer[1];
            if (XYZMinMax_ThisLayer[2] < YMin)
                YMin = XYZMinMax_ThisLayer[2];
            if (XYZMinMax_ThisLayer[3] > YMax)
                YMax = XYZMinMax_ThisLayer[3];
            if (XYZMinMax_ThisLayer[4] < ZMin)
                ZMin = XYZMinMax_ThisLayer[4];
            if (XYZMinMax_ThisLayer[5] > ZMax)
                ZMax = XYZMinMax_ThisLayer[5];
            ZMinLayer[LayerReadCount - 1] = XYZMinMax_ThisLayer[4];
            ZMaxLayer[LayerReadCount - 1] = XYZMinMax_ThisLayer[5];
            if (id == 0)
                std::cout << "Layer = " << LayerReadCount << " Z Bounds are " << XYZMinMax_ThisLayer[4] << " "
                          << XYZMinMax_ThisLayer[5] << std::endl;
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

// Decompose the domain into subdomains on each MPI rank: Calculate MyYSlices and MyYOffset for each rank, where each
// subdomain contains "MyYSlices" in Y, offset from the full domain origin by "MyYOffset" cells in Y
void DomainDecomposition(int id, int np, int &MyYSlices, int &MyYOffset, int &NeighborRank_North,
                         int &NeighborRank_South, int &nx, int &ny, int &nz, long int &LocalDomainSize,
                         bool &AtNorthBoundary, bool &AtSouthBoundary) {

    // Compare total MPI ranks to total Y cells.
    if (np > ny)
        throw std::runtime_error("Error: Cannot run with more MPI ranks than cells in Y (decomposition direction).");

    // Determine which subdomains are at which locations on the grid relative to the others
    InitialDecomposition(id, np, NeighborRank_North, NeighborRank_South, AtNorthBoundary, AtSouthBoundary);
    // Determine, for each MPI process id, the local grid size in x and y (and the offsets in x and y relative to the
    // overall simulation domain)
    MyYOffset = YOffsetCalc(id, ny, np);
    MyYSlices = YMPSlicesCalc(id, ny, np);

    // Add ghost nodes at subdomain overlaps
    AddGhostNodes(NeighborRank_North, NeighborRank_South, MyYSlices, MyYOffset);

    LocalDomainSize = nx * MyYSlices * nz; // Number of cells on this MPI rank
}

// Read in temperature data from files, stored in "RawData", with the appropriate MPI ranks storing the appropriate data
void ReadTemperatureData(int id, double &deltax, double HT_deltax, int &HTtoCAratio, int MyYSlices, int MyYOffset,
                         double YMin, std::vector<std::string> &temp_paths, int NumberOfLayers, int TempFilesInSeries,
                         unsigned int &NumberOfTemperatureDataPoints, std::vector<double> &RawData, int *FirstValue,
                         int *LastValue, bool LayerwiseTempRead, int layernumber) {

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
    // The Y bounds are the region (for this MPI rank) of the physical domain that needs to be
    // read extends past the actual spatial extent of the local domain for purposes of interpolating
    // from HT_deltax to deltax
    int LowerYBound = MyYOffset - (MyYOffset % HTtoCAratio);
    int UpperYBound;
    if (HTtoCAratio == 1)
        UpperYBound = MyYOffset + MyYSlices - 1;
    else
        UpperYBound = MyYOffset + MyYSlices - 1 + HTtoCAratio - (MyYOffset + MyYSlices - 1) % HTtoCAratio;

    // Store raw data relevant to each rank in the vector structure RawData
    // Two passes through reading temperature data files- this is the second pass, reading the actual X/Y/Z/liquidus
    // time/cooling rate data and each rank stores the data relevant to itself in "RawData". With remelting
    // (SimulationType == "RM"), this is the same except that some X/Y/Z coordinates may be repeated in a file, and
    // a "melting time" value is stored in addition to liquidus time and cooling rate
    NumberOfTemperatureDataPoints = 0; // reset to 0 during each call to ReadTemperatureData
    // Second pass through the files - ignore header line
    int FirstLayerToRead, LastLayerToRead;
    if (LayerwiseTempRead) {
        FirstLayerToRead = layernumber;
        LastLayerToRead = layernumber;
    }
    else {
        FirstLayerToRead = 0;
        LastLayerToRead = std::min(NumberOfLayers, TempFilesInSeries) - 1;
    }
    // Which temperature files should be read? Just the one file for layer "layernumber", or all of them?
    for (int LayerReadCount = FirstLayerToRead; LayerReadCount <= LastLayerToRead; LayerReadCount++) {

        std::string tempfile_thislayer;
        if (LayerwiseTempRead) {
            int LayerInSeries = layernumber % TempFilesInSeries;
            tempfile_thislayer = temp_paths[LayerInSeries];
        }
        else
            tempfile_thislayer = temp_paths[LayerReadCount];

        FirstValue[LayerReadCount] = NumberOfTemperatureDataPoints;
        // Read and parse temperature file for either binary or ASCII, storing the appropriate values on each MPI rank
        // within RawData and incrementing NumberOfTemperatureDataPoints appropriately
        bool BinaryInputData = checkTemperatureFileFormat(tempfile_thislayer);
        parseTemperatureData(tempfile_thislayer, YMin, deltax, LowerYBound, UpperYBound, RawData,
                             NumberOfTemperatureDataPoints, BinaryInputData);
        LastValue[LayerReadCount] = NumberOfTemperatureDataPoints;
    } // End loop over all files read for all layers
    RawData.resize(NumberOfTemperatureDataPoints);
    // Determine start values for each layer's data within "RawData", if all layers were read
    if (!(LayerwiseTempRead)) {
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
}

//*****************************************************************************/
// Get the Z coordinate of the lower bound of iteration
int calcZBound_Low(std::string SimulationType, int LayerHeight, int layernumber, double *ZMinLayer, double ZMin,
                   double deltax) {

    int ZBound_Low = -1; // assign dummy initial value
    if (SimulationType == "C") {
        // Not a multilayer problem, top of "layer" is the top of the overall simulation domain
        ZBound_Low = 0;
    }
    else if (SimulationType == "S") {
        // lower bound of domain is an integer multiple of the layer spacing, since the temperature field is the
        // same for every layer
        ZBound_Low = LayerHeight * layernumber;
    }
    else if (SimulationType == "R") {
        // lower bound of domain is based on the data read from the file(s)
        ZBound_Low = round((ZMinLayer[layernumber] - ZMin) / deltax);
    }
    if (ZBound_Low == -1)
        throw std::runtime_error("Error: ZBound_Low went uninitialized, problem type must be C, S, or R");
    return ZBound_Low;
}
//*****************************************************************************/
// Get the Z coordinate of the upper bound of iteration
int calcZBound_High(std::string SimulationType, int SpotRadius, int LayerHeight, int layernumber, double ZMin,
                    double deltax, int nz, double *ZMaxLayer) {

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
        throw std::runtime_error("Error: ZBound_High went uninitialized, problem type must be C, S, or R");
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
int calcLocalActiveDomainSize(int nx, int MyYSlices, int nzActive) {
    int LocalActiveDomainSize = nx * MyYSlices * nzActive;
    return LocalActiveDomainSize;
}
//*****************************************************************************/
// Initialize temperature data for a constrained solidification test problem
void TempInit_DirSolidification(double G, double R, int, int &nx, int &MyYSlices, double deltax, double deltat, int nz,
                                int LocalDomainSize, ViewI &CritTimeStep, ViewF &UndercoolingChange, ViewI &LayerID) {

    // These views are initialized on the host, filled with data, and then copied to the device for layer "layernumber"
    // This view is initialized with zeros
    ViewI_H LayerID_Host("LayerID_H", LocalDomainSize);
    // These views will be filled with non-zero values
    ViewI_H CritTimeStep_Host(Kokkos::ViewAllocateWithoutInitializing("CritTimeStep_H"), LocalDomainSize);
    ViewF_H UndercoolingChange_Host(Kokkos::ViewAllocateWithoutInitializing("UndercoolingChange_H"), LocalDomainSize);

    // Initialize temperature field in Z direction with thermal gradient G set in input file
    // Cells at the bottom surface (Z = 0) are at the liquidus at time step 0 (no wall cells at the bottom boundary)
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                int GlobalD3D1ConvPosition = k * nx * MyYSlices + i * MyYSlices + j;
                UndercoolingChange_Host(GlobalD3D1ConvPosition) = R * deltat;
                CritTimeStep_Host(GlobalD3D1ConvPosition) = (int)((k * G * deltax) / (R * deltat));
            }
        }
    }

    // Copy initialized host data back to device
    CritTimeStep = Kokkos::create_mirror_view_and_copy(device_memory_space(), CritTimeStep_Host);
    LayerID = Kokkos::create_mirror_view_and_copy(device_memory_space(), LayerID_Host);
    UndercoolingChange = Kokkos::create_mirror_view_and_copy(device_memory_space(), UndercoolingChange_Host);
}

// Initialize temperature data for an array of overlapping spot melts (done during simulation initialization, no
// remelting)
void TempInit_SpotNoRemelt(double G, double R, std::string, int id, int &nx, int &MyYSlices, int &MyYOffset,
                           double deltax, double deltat, int, int LocalDomainSize, ViewI &CritTimeStep,
                           ViewF &UndercoolingChange, int LayerHeight, int NumberOfLayers, double FreezingRange,
                           ViewI &LayerID, int NSpotsX, int NSpotsY, int SpotRadius, int SpotOffset) {

    // This view is initialized with -1 for all cells, populated with other data, and later copied to the device
    ViewI_H LayerID_Host(Kokkos::ViewAllocateWithoutInitializing("LayerID_H"), LocalDomainSize);
    Kokkos::deep_copy(LayerID_Host, -1);
    // These views are initialized with zero values on the host, populated with other data, and later copied to the
    // device
    ViewI_H CritTimeStep_Host("CritTimeStep_H", LocalDomainSize);
    ViewF_H UndercoolingChange_Host("UndercoolingChange_H", LocalDomainSize);

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
            int XSpotPos = SpotRadius + (n % NSpotsX) * SpotOffset;
            int YSpotPos = SpotRadius + (n / NSpotsX) * SpotOffset;
            int ZSpotPos = SpotRadius + LayerHeight * layernumber;
            for (int k = 0; k <= SpotRadius; k++) {
                // Distance of this cell from the spot center
                float DistZ = (float)(ZSpotPos - (k + LayerHeight * layernumber));
                for (int i = 0; i < nx; i++) {
                    float DistX = (float)(XSpotPos - i);
                    for (int j = 0; j < MyYSlices; j++) {
                        int YGlobal = j + MyYOffset;
                        float DistY = (float)(YSpotPos - YGlobal);
                        float TotDist = sqrt(DistX * DistX + DistY * DistY + DistZ * DistZ);
                        if (TotDist <= SpotRadius) {
                            int GlobalD3D1ConvPosition =
                                (k + layernumber * LayerHeight) * nx * MyYSlices + i * MyYSlices + j;
                            CritTimeStep_Host(GlobalD3D1ConvPosition) =
                                1 + (int)(((float)(SpotRadius)-TotDist) / IsothermVelocity) + TimeBetweenSpots * n;
                            UndercoolingChange_Host(GlobalD3D1ConvPosition) = R * deltat;
                            LayerID_Host(GlobalD3D1ConvPosition) = layernumber;
                        }
                    }
                }
            }
        }
    }

    // Copy initialized host data back to device
    CritTimeStep = Kokkos::create_mirror_view_and_copy(device_memory_space(), CritTimeStep_Host);
    LayerID = Kokkos::create_mirror_view_and_copy(device_memory_space(), LayerID_Host);
    UndercoolingChange = Kokkos::create_mirror_view_and_copy(device_memory_space(), UndercoolingChange_Host);
}

// For an overlapping spot melt pattern, determine the maximum number of times a cell will melt/solidify as part of a
// layer
int calcMaxSolidificationEventsSpot(int nx, int MyYSlices, int NumberOfSpots, int NSpotsX, int SpotRadius,
                                    int SpotOffset, int MyYOffset) {

    ViewI2D_H MaxSolidificationEvents_Temp("SEvents_Temp", nx, MyYSlices);
    for (int n = 0; n < NumberOfSpots; n++) {
        int XSpotPos = SpotRadius + (n % NSpotsX) * SpotOffset;
        int YSpotPos = SpotRadius + (n / NSpotsX) * SpotOffset;
        for (int i = 0; i < nx; i++) {
            float DistX = (float)(XSpotPos - i);
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
    for (int i = 0; i < nx; i++) {
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
void TempInit_SpotRemelt(int layernumber, double G, double R, std::string, int id, int &nx, int &MyYSlices,
                         int &MyYOffset, double deltax, double deltat, int ZBound_Low, int, int LocalActiveDomainSize,
                         int LocalDomainSize, ViewI &CritTimeStep, ViewF &UndercoolingChange,
                         ViewF &UndercoolingCurrent, int, double FreezingRange, ViewI &LayerID, int NSpotsX,
                         int NSpotsY, int SpotRadius, int SpotOffset, ViewF3D &LayerTimeTempHistory,
                         ViewI &NumberOfSolidificationEvents, ViewI &MeltTimeStep, ViewI &MaxSolidificationEvents,
                         ViewI &SolidificationEventCounter) {

    int NumberOfSpots = NSpotsX * NSpotsY;

    // Temporary host view for the maximum number of times a cell in a given layer will solidify
    ViewI_H MaxSolidificationEvents_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), MaxSolidificationEvents);
    MaxSolidificationEvents_Host(layernumber) =
        calcMaxSolidificationEventsSpot(nx, MyYSlices, NumberOfSpots, NSpotsX, SpotRadius, SpotOffset, MyYOffset);

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

    // First layer - all LayerID values are -1, to later be populated with other values
    if (layernumber == 0)
        Kokkos::deep_copy(LayerID_Host, -1);

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

    for (int n = 0; n < NumberOfSpots; n++) {
        if (id == 0)
            std::cout << "Initializing spot " << n << " on layer " << layernumber << std::endl;
        // Initialize LayerTimeTempHistory data values for this spot/this layer - relative to the layer bottom
        int XSpotPos = SpotRadius + (n % NSpotsX) * SpotOffset;
        int YSpotPos = SpotRadius + (n / NSpotsX) * SpotOffset;
        for (int k = 0; k <= SpotRadius; k++) {
            // Distance of this cell from the spot center
            float DistZ = (float)(SpotRadius - k);
            for (int i = 0; i < nx; i++) {
                float DistX = (float)(XSpotPos - i);
                for (int j = 0; j < MyYSlices; j++) {
                    int YGlobal = j + MyYOffset;
                    float DistY = (float)(YSpotPos - YGlobal);
                    float TotDist = sqrt(DistX * DistX + DistY * DistY + DistZ * DistZ);
                    if (TotDist <= SpotRadius) {
                        int D3D1ConvPosition = k * nx * MyYSlices + i * MyYSlices + j;
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
                    }
                }
            }
        }
    }

    // Initialize data for first melt-solidification event for all cells with data in this layer
    for (int k = 0; k <= SpotRadius; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                int D3D1ConvPosition = k * nx * MyYSlices + i * MyYSlices + j;
                int GlobalD3D1ConvPosition = (k + ZBound_Low) * nx * MyYSlices + i * MyYSlices + j;
                if (NumberOfSolidificationEvents_Host(D3D1ConvPosition) > 0) {
                    MeltTimeStep_Host(GlobalD3D1ConvPosition) = LayerTimeTempHistory_Host(D3D1ConvPosition, 0, 0);
                    CritTimeStep_Host(GlobalD3D1ConvPosition) = LayerTimeTempHistory_Host(D3D1ConvPosition, 0, 1);
                    UndercoolingChange_Host(GlobalD3D1ConvPosition) = LayerTimeTempHistory_Host(D3D1ConvPosition, 0, 2);
                    LayerID_Host(GlobalD3D1ConvPosition) = layernumber;
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
    MaxSolidificationEvents = Kokkos::create_mirror_view_and_copy(device_memory_space(), MaxSolidificationEvents_Host);
    LayerID = Kokkos::create_mirror_view_and_copy(device_memory_space(), LayerID_Host);
    MeltTimeStep = Kokkos::create_mirror_view_and_copy(device_memory_space(), MeltTimeStep_Host);
    CritTimeStep = Kokkos::create_mirror_view_and_copy(device_memory_space(), CritTimeStep_Host);
    UndercoolingChange = Kokkos::create_mirror_view_and_copy(device_memory_space(), UndercoolingChange_Host);
    LayerTimeTempHistory = Kokkos::create_mirror_view_and_copy(device_memory_space(), LayerTimeTempHistory_Host);
    NumberOfSolidificationEvents =
        Kokkos::create_mirror_view_and_copy(device_memory_space(), NumberOfSolidificationEvents_Host);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Spot melt temperature field with remelting for layer " << layernumber
                  << " initialized; each cell will solidify up to " << MaxSolidificationEvents_Host(layernumber)
                  << " times" << std::endl;
}

// Read data from storage, and calculate the normalized x value of the data point
int getTempCoordX(int i, double XMin, double deltax, const std::vector<double> &RawData) {
    int XInt = round((RawData[i] - XMin) / deltax);
    return XInt;
}
// Read data from storage, and calculate the normalized y value of the data point
int getTempCoordY(int i, double YMin, double deltax, const std::vector<double> &RawData) {
    int YInt = round((RawData[i + 1] - YMin) / deltax);
    return YInt;
}
// Read data from storage, and calculate the normalized z value of the data point
int getTempCoordZ(int i, double deltax, const std::vector<double> &RawData, int LayerHeight, int LayerCounter,
                  double *ZMinLayer) {
    int ZInt = round((RawData[i + 2] + deltax * LayerHeight * LayerCounter - ZMinLayer[LayerCounter]) / deltax);
    return ZInt;
}
// Read data from storage, obtain melting time
double getTempCoordTM(int i, const std::vector<double> &RawData) {
    double TMelting = RawData[i + 3];
    return TMelting;
}
// Read data from storage, obtain liquidus time
double getTempCoordTL(int i, const std::vector<double> &RawData) {
    double TLiquidus = RawData[i + 4];
    return TLiquidus;
}
// Read data from storage, obtain cooling rate
double getTempCoordCR(int i, const std::vector<double> &RawData) {
    double CoolingRate = RawData[i + 5];
    return CoolingRate;
}

// Initialize temperature data for a problem using the reduced/sparse data format and input temperature data from
// file(s)
void TempInit_ReadDataNoRemelt(int id, int &nx, int &MyYSlices, int &MyYOffset, double deltax, int HTtoCAratio,
                               double deltat, int, int LocalDomainSize, ViewI &CritTimeStep, ViewF &UndercoolingChange,
                               double XMin, double YMin, double ZMin, double *ZMinLayer, double *ZMaxLayer,
                               int LayerHeight, int NumberOfLayers, int *FinishTimeStep, double FreezingRange,
                               ViewI &LayerID, int *FirstValue, int *LastValue, std::vector<double> RawData, int ny) {

    // These views are initialized to zeros on the host, filled with data, and then copied to the device for layer
    // "layernumber"
    ViewI_H LayerID_Host(Kokkos::ViewAllocateWithoutInitializing("LayerID_H"), LocalDomainSize);
    ViewI_H CritTimeStep_Host("CritTimeStep_H", LocalDomainSize);
    ViewF_H UndercoolingChange_Host("UndercoolingChange_H", LocalDomainSize);

    // Temperature data read
    // If HTtoCAratio > 1, an interpolation of input temperature data is needed
    // The Y bounds are the region (for this MPI rank) of the physical domain that needs to be
    // read extends past the actual spatial extent of the local domain for purposes of interpolating
    // from HT_deltax to deltax
    int LowerYBound = MyYOffset - (MyYOffset % HTtoCAratio);
    int UpperYBound = MyYOffset + MyYSlices - 1 + HTtoCAratio - ((MyYOffset + MyYSlices - 1) % HTtoCAratio);
    // Make sure that upper Y bound doesn't extend beyond simulation domain
    if (UpperYBound >= ny)
        UpperYBound = ny - 1;

    // LayerID = -1 for cells that don't solidify as part of any layer of the multilayer problem
    Kokkos::deep_copy(LayerID_Host, -1);

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
            for (int i = 0; i < nx; i++) {
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
        if (id == 0)
            std::cout << "Range for layer " << LayerCounter << " on rank 0 is " << StartRange << " to " << EndRange
                      << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        for (int i = StartRange; i < EndRange; i += 6) {

            // Get the integer X, Y, Z coordinates associated with this data point, along with TL and CR values
            int XInt = getTempCoordX(i, XMin, deltax, RawData);
            int YInt = getTempCoordY(i, YMin, deltax, RawData);
            int ZInt = getTempCoordZ(i, deltax, RawData, LayerHeight, LayerCounter, ZMinLayer);
            double TLiquidus = getTempCoordTL(i, RawData);
            // Liquidus time/cooling rate - only keep values for the last time that this point went below the liquidus
            if (TLiquidus > CritTL[ZInt][XInt][YInt - LowerYBound]) {
                CritTL[ZInt][XInt][YInt - LowerYBound] = TLiquidus;
                if (TLiquidus < SmallestTime) {
                    // Store smallest read TLiquidus value over all cells
                    SmallestTime = RawData[i];
                }
                double CoolingRate = getTempCoordCR(i, RawData);
                CR[ZInt][XInt][YInt - LowerYBound] = CoolingRate;
                double SolidusTime =
                    CritTL[ZInt][XInt][YInt - LowerYBound] + FreezingRange / CR[ZInt][XInt][YInt - LowerYBound];
                if (SolidusTime > LargestTime) {
                    // Store largest TSolidus value (based on liquidus/cooling rate/freezing range) over all cells
                    LargestTime = SolidusTime;
                }
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
                for (int i = 0; i < nx; i++) {
                    int LowX = i - (i % HTtoCAratio);
                    int HighX = LowX + HTtoCAratio;
                    double FHighX = (double)(i - LowX) / (double)(HTtoCAratio);
                    double FLowX = 1.0 - FHighX;
                    if (HighX >= nx)
                        HighX = LowX;

                    for (int j = 0; j <= UpperYBound - LowerYBound; j++) {
                        int LowY = j - (j % HTtoCAratio);
                        int HighY = LowY + HTtoCAratio;
                        double FHighY = (float)(j - LowY) / (float)(HTtoCAratio);
                        double FLowY = 1.0 - FHighY;
                        if (HighY > UpperYBound - LowerYBound)
                            HighY = LowY;
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
            for (int i = 0; i < nx; i++) {
                for (int jj = LowerYBound; jj <= UpperYBound; jj++) {
                    if ((jj >= MyYOffset) && (jj < MyYOffset + MyYSlices)) {
                        int Adj_j = jj - MyYOffset;
                        // Liquidus time normalized to the time at which the layer started solidifying
                        double CTLiq = CritTL[k][i][jj - LowerYBound] - SmallestTime_Global;
                        if (CTLiq > 0) {
                            // Where does this layer's temperature data belong on the global (including all layers)
                            // grid? Adjust Z coordinate by ZMin
                            int ZOffset = round((ZMinLayer[LayerCounter] - ZMin) / deltax) + k;
                            int Coord3D1D = ZOffset * nx * MyYSlices + i * MyYSlices + Adj_j;
                            CritTimeStep_Host(Coord3D1D) = round(CTLiq / deltat);
                            LayerID_Host(Coord3D1D) = LayerCounter;
                            UndercoolingChange_Host(Coord3D1D) = std::abs(CR[k][i][jj - LowerYBound]) * deltat;
                        }
                    }
                }
            }
        }
    } // End read over all temperature files and placement of data

    // Copy initialized host data back to device
    CritTimeStep = Kokkos::create_mirror_view_and_copy(device_memory_space(), CritTimeStep_Host);
    LayerID = Kokkos::create_mirror_view_and_copy(device_memory_space(), LayerID_Host);
    UndercoolingChange = Kokkos::create_mirror_view_and_copy(device_memory_space(), UndercoolingChange_Host);
}

// Calculate the number of times that a cell in layer "layernumber" undergoes melting/solidification, and store in
// MaxSolidificationEvents_Host
void calcMaxSolidificationEventsR(int id, int layernumber, int TempFilesInSeries, ViewI_H MaxSolidificationEvents_Host,
                                  int StartRange, int EndRange, std::vector<double> RawData, double XMin, double YMin,
                                  double deltax, double *ZMinLayer, int LayerHeight, int nx, int MyYSlices,
                                  int MyYOffset, int LocalActiveDomainSize) {

    if (layernumber > TempFilesInSeries) {
        // Use the value from a previously checked layer, since the time-temperature history is reused
        if (TempFilesInSeries == 1) {
            // All layers have the same temperature data, MaxSolidificationEvents for this layer is the same as the last
            MaxSolidificationEvents_Host(layernumber) = MaxSolidificationEvents_Host(layernumber - 1);
        }
        else {
            // All layers have different temperature data but in a repeating pattern
            int RepeatedFile = layernumber % TempFilesInSeries;
            MaxSolidificationEvents_Host(layernumber) = MaxSolidificationEvents_Host(RepeatedFile);
        }
    }
    else {
        // Need to calculate MaxSolidificationEvents(layernumber) from the values in RawData
        // Init to 0
        ViewI_H TempMeltCount("TempMeltCount", LocalActiveDomainSize);

        for (int i = StartRange; i < EndRange; i += 6) {

            // Get the integer X, Y, Z coordinates associated with this data point
            int XInt = getTempCoordX(i, XMin, deltax, RawData);
            int YInt = getTempCoordY(i, YMin, deltax, RawData);
            int ZInt = getTempCoordZ(i, deltax, RawData, LayerHeight, layernumber, ZMinLayer);
            // Convert to 1D coordinate in the current layer's domain
            int D3D1ConvPosition = ZInt * nx * MyYSlices + XInt * MyYSlices + (YInt - MyYOffset);
            TempMeltCount(D3D1ConvPosition)++;
        }
        int MaxCount = 0;
        for (int i = 0; i < LocalActiveDomainSize; i++) {
            if (TempMeltCount(i) > MaxCount)
                MaxCount = TempMeltCount(i);
        }
        int MaxCountGlobal;
        MPI_Allreduce(&MaxCount, &MaxCountGlobal, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MaxSolidificationEvents_Host(layernumber) = MaxCountGlobal;
    }
    if (id == 0)
        std::cout << "The maximum number of melting/solidification events during layer " << layernumber << " is "
                  << MaxSolidificationEvents_Host(layernumber) << std::endl;
}

// Initialize temperature fields for this layer if remelting is considered and data comes from files
void TempInit_ReadDataRemelt(int layernumber, int id, int nx, int MyYSlices, int, int LocalActiveDomainSize,
                             int LocalDomainSize, int MyYOffset, double &deltax, double deltat, double FreezingRange,
                             ViewF3D &LayerTimeTempHistory, ViewI &NumberOfSolidificationEvents,
                             ViewI &MaxSolidificationEvents, ViewI &MeltTimeStep, ViewI &CritTimeStep,
                             ViewF &UndercoolingChange, ViewF &UndercoolingCurrent, double XMin, double YMin,
                             double *ZMinLayer, int LayerHeight, int nzActive, int ZBound_Low, int *FinishTimeStep,
                             ViewI &LayerID, int *FirstValue, int *LastValue, std::vector<double> RawData,
                             ViewI &SolidificationEventCounter, int TempFilesInSeries) {

    // Data was already read into the "RawData" temporary data structure
    // Determine which section of "RawData" is relevant for this layer of the overall domain
    int StartRange = FirstValue[layernumber];
    int EndRange = LastValue[layernumber];

    // Resize device views to have sizes compatible with the temporary host views
    // Copy MaxSolidificationEvents back to the host as part of resizing LayerTimeTempHistory
    ViewI_H MaxSolidificationEvents_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), MaxSolidificationEvents);

    // Get the maximum number of times a cell in layer "layernumber" will undergo melting/solidification
    // Store in the host view "MaxSolidificationEvents_Host"
    calcMaxSolidificationEventsR(id, layernumber, TempFilesInSeries, MaxSolidificationEvents_Host, StartRange, EndRange,
                                 RawData, XMin, YMin, deltax, ZMinLayer, LayerHeight, nx, MyYSlices, MyYOffset,
                                 LocalActiveDomainSize);
    // With MaxSolidificationEvents_Host(layernumber) known, can resize LayerTimeTempHistory
    Kokkos::resize(LayerTimeTempHistory, LocalActiveDomainSize, MaxSolidificationEvents_Host(layernumber), 3);
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
    ViewF3D_H LayerTimeTempHistory_Host("TimeTempHistory_H", LocalActiveDomainSize,
                                        MaxSolidificationEvents_Host(layernumber), 3);
    ViewI_H NumberOfSolidificationEvents_Host("NumSEvents_H", LocalActiveDomainSize);

    // First layer - all LayerID values are -1, to later be populated with other values
    if (layernumber == 0)
        Kokkos::deep_copy(LayerID_Host, -1);

    double LargestTime = 0;
    double LargestTime_Global = 0;
    if (id == 0)
        std::cout << "Range of raw data for layer " << layernumber << " on rank 0 is " << StartRange << " to "
                  << EndRange << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = StartRange; i < EndRange; i += 6) {

        // Get the integer X, Y, Z coordinates associated with this data point, along with the associated TM, TL, CR
        // values
        int XInt = getTempCoordX(i, XMin, deltax, RawData);
        int YInt = getTempCoordY(i, YMin, deltax, RawData);
        int ZInt = getTempCoordZ(i, deltax, RawData, LayerHeight, layernumber, ZMinLayer);
        double TMelting = getTempCoordTM(i, RawData);
        double TLiquidus = getTempCoordTL(i, RawData);
        double CoolingRate = getTempCoordCR(i, RawData);

        // 1D cell coordinate on this MPI rank's domain
        int D3D1ConvPosition = ZInt * nx * MyYSlices + XInt * MyYSlices + (YInt - MyYOffset);
        // Store TM, TL, CR values for this solidification event in LayerTimeTempHistory
        LayerTimeTempHistory_Host(D3D1ConvPosition, NumberOfSolidificationEvents_Host(D3D1ConvPosition), 0) =
            round(TMelting / deltat) + 1;
        LayerTimeTempHistory_Host(D3D1ConvPosition, NumberOfSolidificationEvents_Host(D3D1ConvPosition), 1) =
            round(TLiquidus / deltat) + 1;
        LayerTimeTempHistory_Host(D3D1ConvPosition, NumberOfSolidificationEvents_Host(D3D1ConvPosition), 2) =
            std::abs(CoolingRate) * deltat;
        // Increment number of solidification events for this cell
        NumberOfSolidificationEvents_Host(D3D1ConvPosition)++;
        // Estimate of the time step where the last possible solidification is expected to occur
        double SolidusTime = TLiquidus + FreezingRange / CoolingRate;
        if (SolidusTime > LargestTime)
            LargestTime = SolidusTime;
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
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                int D3D1ConvPosition = k * nx * MyYSlices + i * MyYSlices + j;
                int GlobalD3D1ConvPosition = GlobalZ * nx * MyYSlices + i * MyYSlices + j;
                if (LayerTimeTempHistory_Host(D3D1ConvPosition, 0, 0) > 0) {
                    // This cell undergoes solidification in layer "layernumber" at least once
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
    MaxSolidificationEvents = Kokkos::create_mirror_view_and_copy(device_memory_space(), MaxSolidificationEvents_Host);
    LayerID = Kokkos::create_mirror_view_and_copy(device_memory_space(), LayerID_Host);
    MeltTimeStep = Kokkos::create_mirror_view_and_copy(device_memory_space(), MeltTimeStep_Host);
    CritTimeStep = Kokkos::create_mirror_view_and_copy(device_memory_space(), CritTimeStep_Host);
    UndercoolingChange = Kokkos::create_mirror_view_and_copy(device_memory_space(), UndercoolingChange_Host);
    LayerTimeTempHistory = Kokkos::create_mirror_view_and_copy(device_memory_space(), LayerTimeTempHistory_Host);
    NumberOfSolidificationEvents =
        Kokkos::create_mirror_view_and_copy(device_memory_space(), NumberOfSolidificationEvents_Host);

    if (id == 0)
        std::cout << "Layer " << layernumber << " temperature field is from Z = " << ZBound_Low << " through "
                  << nzActive + ZBound_Low - 1 << " of the global domain" << std::endl;
}

//*****************************************************************************/
// Initialize grain orientations and unit vectors
void OrientationInit(int, int &NGrainOrientations, ViewF &GrainOrientationData, std::string GrainOrientationFile,
                     int ValsPerLine) {

    // Read file of grain orientations
    std::ifstream O;
    O.open(GrainOrientationFile);

    // Line 1 is the number of orientation values to read (if not specified already)
    std::string ValueRead;
    getline(O, ValueRead);
    NGrainOrientations = getInputInt(ValueRead);

    // Temporary host view for storing grain orientations read from file
    ViewF_H GrainOrientationData_Host(Kokkos::ViewAllocateWithoutInitializing("GrainOrientationData_H"),
                                      ValsPerLine * NGrainOrientations);
    // Populate data structure for grain orientation data
    for (int i = 0; i < NGrainOrientations; i++) {
        std::vector<std::string> ParsedLine(ValsPerLine);
        std::string ReadLine;
        if (!getline(O, ReadLine))
            break;
        splitString(ReadLine, ParsedLine, ValsPerLine);
        // Place the 3 grain orientation angles or 9 rotation matrix components into the orientation data view
        for (int Comp = 0; Comp < ValsPerLine; Comp++) {
            GrainOrientationData_Host(ValsPerLine * i + Comp) = getInputFloat(ParsedLine[Comp]);
        }
    }
    O.close();

    // Resize device view and orientation data to device
    Kokkos::realloc(GrainOrientationData, ValsPerLine * NGrainOrientations);
    GrainOrientationData = Kokkos::create_mirror_view_and_copy(device_memory_space(), GrainOrientationData_Host);
}

// Initializes cell types and epitaxial Grain ID values where substrate grains are active cells on the bottom surface of
// the constrained domain. Also initialize active cell data structures associated with the substrate grains
void SubstrateInit_ConstrainedGrowth(int id, double FractSurfaceSitesActive, int MyYSlices, int nx, int ny,
                                     int MyYOffset, NList NeighborX, NList NeighborY, NList NeighborZ,
                                     ViewF GrainUnitVector, int NGrainOrientations, ViewI CellType, ViewI GrainID,
                                     ViewF DiagonalLength, ViewF DOCenter, ViewF CritDiagonalLength, double RNGSeed,
                                     int np, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, int BufSizeX,
                                     bool AtNorthBoundary, bool AtSouthBoundary) {

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
    ViewI ActCellX_Device = Kokkos::create_mirror_view_and_copy(device_memory_space(), ActCellX_Host);
    ViewI ActCellY_Device = Kokkos::create_mirror_view_and_copy(device_memory_space(), ActCellY_Host);

    // Start with all cells as liquid prior to locating substrate grain seeds
    Kokkos::deep_copy(CellType, Liquid);

    // Determine which grains/active cells belong to which MPI ranks
    Kokkos::parallel_for(
        "ConstrainedGrainInit", SubstrateActCells, KOKKOS_LAMBDA(const int &n) {
            // What are the X and Y coordinates of this active cell relative to the X and Y bounds of this rank?
            if ((ActCellY_Device(n) >= MyYOffset) && (ActCellY_Device(n) < MyYOffset + MyYSlices)) {
                // Convert X and Y coordinates to values relative to this MPI rank's grid (Z = 0 for these active cells,
                // at bottom surface) GrainIDs come from the position on the list of substrate active cells to avoid
                // reusing the same value
                int GlobalX = ActCellX_Device(n);
                int LocalY = ActCellY_Device(n) - MyYOffset;
                int D3D1ConvPosition = GlobalX * MyYSlices + LocalY;
                CellType(D3D1ConvPosition) = Active;
                GrainID(D3D1ConvPosition) = n + 1; // assign GrainID > 0 to epitaxial seeds
                // Initialize active cell data structures
                int GlobalY = LocalY + MyYOffset;
                int GlobalZ = 0;
                // Initialize new octahedron
                createNewOctahedron(D3D1ConvPosition, DiagonalLength, DOCenter, GlobalX, GlobalY, GlobalZ);

                // The orientation for the new grain will depend on its Grain ID
                int MyOrientation = getGrainOrientation(n + 1, NGrainOrientations);
                float cx = GlobalX + 0.5;
                float cy = GlobalY + 0.5;
                float cz = GlobalZ + 0.5;
                // Calculate critical values at which this active cell leads to the activation of a neighboring liquid
                // cell. Octahedron center and cell center overlap for octahedra created as part of a new grain
                calcCritDiagonalLength(D3D1ConvPosition, cx, cy, cz, cx, cy, cz, NeighborX, NeighborY, NeighborZ,
                                       MyOrientation, GrainUnitVector, CritDiagonalLength);
                // If this new active cell is in the halo region, load the send buffers
                if (np > 1) {

                    float GhostGID = GrainID(D3D1ConvPosition);
                    float GhostDOCX = GlobalX + 0.5;
                    float GhostDOCY = GlobalY + 0.5;
                    float GhostDOCZ = GlobalZ + 0.5;
                    float GhostDL = 0.01;
                    // Collect data for the ghost nodes, if necessary
                    loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, BufSizeX, MyYSlices, GlobalX,
                                   LocalY, 0, AtNorthBoundary, AtSouthBoundary, BufferSouthSend, BufferNorthSend);
                } // End if statement for serial/parallel code
            }
        });
    if (id == 0)
        std::cout << "Number of substrate active cells across all ranks: " << SubstrateActCells << std::endl;
}

// Initializes Grain ID values where the substrate comes from a file
void SubstrateInit_FromFile(std::string SubstrateFileName, int nz, int nx, int MyYSlices, int MyYOffset, int id,
                            ViewI &GrainID_Device, int nzActive, bool BaseplateThroughPowder) {

    // Assign GrainID values to cells that are part of the substrate - read values from file and initialize using
    // temporary host view
    ViewI_H GrainID_Host(Kokkos::ViewAllocateWithoutInitializing("GrainID_Host"), nx * MyYSlices * nz);
    std::ifstream Substrate;
    Substrate.open(SubstrateFileName);
    int Substrate_LowY = MyYOffset;
    int Substrate_HighY = MyYOffset + MyYSlices;
    int nxS, nyS, nzS;
    int BaseplateSizeZ; // in CA cells
    if (BaseplateThroughPowder)
        BaseplateSizeZ = nz; // baseplate microstructure used as entire domain's initial condition
    else
        BaseplateSizeZ = nzActive; // baseplate microstructure is layer 0's initial condition
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
    if ((id == 0) && (nzS < BaseplateSizeZ)) {
        // Do not allow simulation if there is inssufficient substrate data in the specified file
        std::string error = "Error: only " + std::to_string(nzS) +
                            " layers of substrate data are present in the file " + SubstrateFileName + " ; at least " +
                            std::to_string(BaseplateSizeZ) +
                            " layers of substrate data are required to simulate the specified solidification problem";
        throw std::runtime_error(error);
    }

    // Assign GrainID values to all cells - cells that will be part of the melt pool footprint may still need their
    // initial GrainID
    for (int k = 0; k < nzS; k++) {
        if (k == BaseplateSizeZ)
            break;
        for (int j = 0; j < nyS; j++) {
            for (int i = 0; i < nxS; i++) {
                std::string GIDVal;
                getline(Substrate, GIDVal);
                if ((j >= Substrate_LowY) && (j < Substrate_HighY)) {
                    int CAGridLocation;
                    CAGridLocation = k * nx * MyYSlices + i * MyYSlices + (j - MyYOffset);
                    GrainID_Host(CAGridLocation) = stoi(GIDVal, nullptr, 10);
                }
            }
        }
    }
    Substrate.close();
    // Copy GrainIDs read from file to device
    GrainID_Device = Kokkos::create_mirror_view_and_copy(device_memory_space(), GrainID_Host);
    if (id == 0)
        std::cout << "Substrate file read complete" << std::endl;
}

// Initializes Grain ID values where the baseplate is generated using an input grain spacing and a Voronoi Tessellation
void BaseplateInit_FromGrainSpacing(float SubstrateGrainSpacing, int nx, int ny, double *ZMinLayer, double *ZMaxLayer,
                                    int MyYSlices, int MyYOffset, int id, double deltax, ViewI GrainID, double RNGSeed,
                                    int &NextLayer_FirstEpitaxialGrainID, int nz, double BaseplateThroughPowder) {

    // Number of cells to assign GrainID
    int BaseplateSizeZ; // in CA cells
    if (BaseplateThroughPowder)
        BaseplateSizeZ = nz; // baseplate microstructure used as entire domain's initial condition
    else
        BaseplateSizeZ = round((ZMaxLayer[0] - ZMinLayer[0]) / deltax) +
                         1; // baseplate microstructure is layer 0's initial condition
    std::mt19937_64 gen(RNGSeed);

    // Based on the baseplate volume (convert to cubic microns to match units) and the substrate grain spacing,
    // determine the number of baseplate grains
    int BaseplateVolume = nx * ny * BaseplateSizeZ;
    double BaseplateVolume_microns = BaseplateVolume * pow(deltax, 3) * pow(10, 18);
    double SubstrateMeanGrainVolume_microns = pow(SubstrateGrainSpacing, 3);
    int NumberOfBaseplateGrains = round(BaseplateVolume_microns / SubstrateMeanGrainVolume_microns);
    // Need at least 1 baseplate grain, cannot have more baseplate grains than cells in the baseplate
    NumberOfBaseplateGrains = std::max(NumberOfBaseplateGrains, 1);
    NumberOfBaseplateGrains = std::min(NumberOfBaseplateGrains, nx * ny * BaseplateSizeZ);
    // TODO: Use device RNG to generate baseplate grain locations, instead of host with copy
    // List of potential grain IDs (starting at 1) - index corresponds to the associated CA cell location
    // Assign positive values for indices 0 through NumberOfBaseplateGrains-1, assign zeros to the remaining indices
    std::vector<int> BaseplateGrainLocations(BaseplateVolume);
    std::vector<int> BaseplateGrainIDs(BaseplateVolume);
    for (int i = 0; i < BaseplateVolume; i++) {
        BaseplateGrainLocations[i] = i;
        if (i < NumberOfBaseplateGrains)
            BaseplateGrainIDs[i] = i + 1;
        else
            BaseplateGrainIDs[i] = 0;
    }
    // Shuffle list of grain IDs
    std::shuffle(BaseplateGrainIDs.begin(), BaseplateGrainIDs.end(), gen);

    // Create views of baseplate grain IDs and locations - copying BaseplateGrainIDs and BaseplateGrainLocations values
    // only for indices where BaseplateGrainIDs(i) != 0 (cells with BaseplateGrainIDs(i) = 0 were not assigned baseplate
    // center locations)
    ViewI_H BaseplateGrainLocations_Host(Kokkos::ViewAllocateWithoutInitializing("BaseplateGrainLocations_Host"),
                                         NumberOfBaseplateGrains);
    ViewI_H BaseplateGrainIDs_Host(Kokkos::ViewAllocateWithoutInitializing("BaseplateGrainIDs_Host"),
                                   NumberOfBaseplateGrains);
    int count = 0;
    for (int i = 0; i < BaseplateVolume; i++) {
        if (BaseplateGrainIDs[i] != 0) {
            BaseplateGrainLocations_Host(count) = BaseplateGrainLocations[i];
            BaseplateGrainIDs_Host(count) = BaseplateGrainIDs[i];
            count++;
        }
    }

    // Copy baseplate views to the device
    ViewI BaseplateGrainIDs_Device = Kokkos::create_mirror_view_and_copy(device_memory_space(), BaseplateGrainIDs_Host);
    ViewI BaseplateGrainLocations_Device =
        Kokkos::create_mirror_view_and_copy(device_memory_space(), BaseplateGrainLocations_Host);
    if (id == 0) {
        std::cout << "Baseplate spanning domain coordinates Z = 0 through " << BaseplateSizeZ - 1 << std::endl;
        std::cout << "Number of baseplate grains: " << NumberOfBaseplateGrains << std::endl;
    }

    // First, assign cells that are associated with grain centers the appropriate non-zero GrainID values (assumes
    // GrainID values were initialized to zeros)
    Kokkos::parallel_for(
        "BaseplateInit", NumberOfBaseplateGrains, KOKKOS_LAMBDA(const int &n) {
            int BaseplateGrainLoc = BaseplateGrainLocations_Device(n);
            // x, y, z associated with baseplate grain "n", at 1D coordinate "BaseplateGrainLoc"
            int z_n = BaseplateGrainLoc / (nx * ny);
            int Rem = BaseplateGrainLoc % (nx * ny);
            int x_n = Rem / ny;
            int y_n = Rem % ny;
            if ((y_n >= MyYOffset) && (y_n < MyYOffset + MyYSlices)) {
                // This grain is associated with a cell on this MPI rank
                int CAGridLocation = z_n * nx * MyYSlices + x_n * MyYSlices + (y_n - MyYOffset);
                GrainID(CAGridLocation) = BaseplateGrainIDs_Device(n);
            }
        });
    Kokkos::fence();
    // For cells that are not associated with grain centers, assign them the GrainID of the nearest grain center
    Kokkos::parallel_for(
        "BaseplateGen",
        Kokkos::MDRangePolicy<Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>(
            {0, 0, 0}, {BaseplateSizeZ, nx, MyYSlices}),
        KOKKOS_LAMBDA(const int k, const int i, const int j) {
            int CAGridLocation = k * nx * MyYSlices + i * MyYSlices + j;
            if (GrainID(CAGridLocation) == 0) {
                // This cell needs to be assigned a GrainID value
                // Check each possible baseplate grain center to find the closest one
                float MinDistanceToThisGrain = nx * ny * BaseplateSizeZ;
                int MinDistanceToThisGrain_GrainID = 0;
                for (int n = 0; n < NumberOfBaseplateGrains; n++) {
                    // Baseplate grain center at x_n, y_n, z_n - how far is the cell at i, j+MyYOffset, k?
                    int z_n = BaseplateGrainLocations_Device(n) / (nx * ny);
                    int Rem = BaseplateGrainLocations_Device(n) % (nx * ny);
                    int x_n = Rem / ny;
                    int y_n = Rem % ny;
                    float DistanceToThisGrainX = i - x_n;
                    float DistanceToThisGrainY = (j + MyYOffset) - y_n;
                    float DistanceToThisGrainZ = k - z_n;
                    float DistanceToThisGrain = sqrtf(DistanceToThisGrainX * DistanceToThisGrainX +
                                                      DistanceToThisGrainY * DistanceToThisGrainY +
                                                      DistanceToThisGrainZ * DistanceToThisGrainZ);
                    if (DistanceToThisGrain < MinDistanceToThisGrain) {
                        // This is the closest grain center to cell at "CAGridLocation" - update values
                        MinDistanceToThisGrain = DistanceToThisGrain;
                        MinDistanceToThisGrain_GrainID = BaseplateGrainIDs_Device(n);
                    }
                }
                // GrainID associated with the closest baseplate grain center
                GrainID(CAGridLocation) = MinDistanceToThisGrain_GrainID;
            }
        });

    NextLayer_FirstEpitaxialGrainID =
        NumberOfBaseplateGrains + 2; // avoid reusing GrainID in next layer's powder grain structure
    if (id == 0)
        std::cout << "Baseplate grain structure initialized" << std::endl;
}

// Each layer's top Z coordinates are seeded with CA-cell sized substrate grains (emulating bulk nucleation alongside
// the edges of partially melted powder particles)
void PowderInit(int layernumber, int nx, int ny, int LayerHeight, double *ZMaxLayer, double ZMin, double deltax,
                int MyYSlices, int MyYOffset, int id, ViewI GrainID, double RNGSeed,
                int &NextLayer_FirstEpitaxialGrainID, double PowderActiveFraction) {

    // On all ranks, generate list of powder grain IDs (starting with NextLayer_FirstEpitaxialGrainID, and shuffle them
    // so that their locations aren't sequential and depend on the RNGSeed (different for each layer)
    std::mt19937_64 gen(RNGSeed + layernumber);
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // TODO: This should be performed on the device, rather than the host
    int PowderLayerCells = nx * ny * LayerHeight;
    int PowderLayerAssignedCells = round(static_cast<double>(PowderLayerCells) * PowderActiveFraction);
    std::vector<int> PowderGrainIDs(PowderLayerCells, 0);
    for (int n = 0; n < PowderLayerAssignedCells; n++) {
        PowderGrainIDs[n] = n + NextLayer_FirstEpitaxialGrainID; // assigned a nonzero GrainID
    }
    std::shuffle(PowderGrainIDs.begin(), PowderGrainIDs.end(), gen);
    // Copy powder layer GrainIDs into a host view, then a device view
    ViewI_H PowderGrainIDs_Host(Kokkos::ViewAllocateWithoutInitializing("PowderGrainIDs_Host"), PowderLayerCells);
    for (int n = 0; n < PowderLayerCells; n++) {
        PowderGrainIDs_Host(n) = PowderGrainIDs[n];
    }
    ViewI PowderGrainIDs_Device = Kokkos::create_mirror_view_and_copy(device_memory_space(), PowderGrainIDs_Host);
    // Associate powder grain IDs with CA cells in the powder layer
    // Use bounds from temperature field for this layer to determine which cells are part of the powder
    int PowderTopZ = round((ZMaxLayer[layernumber] - ZMin) / deltax) + 1;
    int PowderBottomZ = PowderTopZ - LayerHeight;
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Initializing powder layer for Z = " << PowderBottomZ << " through " << PowderTopZ - 1 << " ("
                  << nx * ny * (PowderTopZ - PowderBottomZ) << " cells)" << std::endl;

    int PowderStart = nx * ny * PowderBottomZ;
    int PowderEnd = nx * ny * PowderTopZ;
    if (id == 0)
        std::cout << "Powder layer has " << PowderLayerAssignedCells
                  << " cells assigned new grain ID values, ranging from " << NextLayer_FirstEpitaxialGrainID
                  << " through " << NextLayer_FirstEpitaxialGrainID + PowderLayerAssignedCells - 1 << std::endl;
    Kokkos::parallel_for(
        "PowderGrainInit", Kokkos::RangePolicy<>(PowderStart, PowderEnd), KOKKOS_LAMBDA(const int &n) {
            int GlobalZ = n / (nx * ny);
            int Rem = n % (nx * ny);
            int GlobalX = Rem / ny;
            int GlobalY = Rem % ny;
            // Is this powder coordinate in X and Y in bounds for this rank? Is the grain id of this site unassigned
            // (wasn't captured during solidification of the previous layer)?
            int GlobalD3D1ConvPosition = GlobalZ * nx * MyYSlices + GlobalX * MyYSlices + (GlobalY - MyYOffset);
            if ((GlobalY >= MyYOffset) && (GlobalY < MyYOffset + MyYSlices) && (GrainID(GlobalD3D1ConvPosition) == 0))
                GrainID(GlobalD3D1ConvPosition) = PowderGrainIDs_Device(n - PowderStart);
        });
    Kokkos::fence();

    // Update NextLayer_FirstEpitaxialGrainID for next layer
    NextLayer_FirstEpitaxialGrainID += PowderLayerAssignedCells;
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Initialized powder grain structure for layer " << layernumber << std::endl;
}

//*****************************************************************************/
// Initializes cells at border of solid and liquid as active type - performed on device
void CellTypeInit_NoRemelt(int layernumber, int id, int np, int nx, int MyYSlices, int MyYOffset, int ZBound_Low,
                           int nz, int LocalActiveDomainSize, int LocalDomainSize, ViewI CellType, ViewI CritTimeStep,
                           NList NeighborX, NList NeighborY, NList NeighborZ, int NGrainOrientations,
                           ViewF GrainUnitVector, ViewF DiagonalLength, ViewI GrainID, ViewF CritDiagonalLength,
                           ViewF DOCenter, ViewI LayerID, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend,
                           int BufSizeX, bool AtNorthBoundary, bool AtSouthBoundary) {

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
                int GlobalZ = GlobalD3D1ConvPosition / (nx * MyYSlices);
                int Rem = GlobalD3D1ConvPosition % (nx * MyYSlices);
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
                        int MyNeighborX = RankX + NeighborX[l];
                        int MyNeighborY = RankY + NeighborY[l];
                        int MyNeighborZ = GlobalZ + NeighborZ[l];
                        int NeighborD3D1ConvPosition =
                            MyNeighborZ * nx * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
                        if ((MyNeighborX >= 0) && (MyNeighborX < nx) && (MyNeighborY >= 0) &&
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
            int RankZ = D3D1ConvPosition / (nx * MyYSlices);
            int Rem = D3D1ConvPosition % (nx * MyYSlices);
            int GlobalX = Rem / MyYSlices;
            int RankY = Rem % MyYSlices;
            int GlobalZ = RankZ + ZBound_Low;
            int GlobalD3D1ConvPosition = GlobalZ * nx * MyYSlices + GlobalX * MyYSlices + RankY;
            if ((CellType(GlobalD3D1ConvPosition) == Active) && (LayerID(GlobalD3D1ConvPosition) == layernumber)) {
                // This cell was marked as active previously - initialize active cell data structures
                int GlobalY = RankY + MyYOffset;
                int RankZ = GlobalZ - ZBound_Low;
                int D3D1ConvPosition = RankZ * nx * MyYSlices + GlobalX * MyYSlices + RankY;
                int MyGrainID = GrainID(GlobalD3D1ConvPosition);
                // Initialize new octahedron
                createNewOctahedron(D3D1ConvPosition, DiagonalLength, DOCenter, GlobalX, GlobalY, GlobalZ);

                // The orientation for the new grain will depend on its Grain ID
                int MyOrientation = getGrainOrientation(MyGrainID, NGrainOrientations);
                float cx = GlobalX + 0.5;
                float cy = GlobalY + 0.5;
                float cz = GlobalZ + 0.5;
                // Calculate critical values at which this active cell leads to the activation of a neighboring liquid
                // cell. Octahedron center and cell center overlap for octahedra created as part of a new grain
                calcCritDiagonalLength(D3D1ConvPosition, cx, cy, cz, cx, cy, cz, NeighborX, NeighborY, NeighborZ,
                                       MyOrientation, GrainUnitVector, CritDiagonalLength);
                // If this new active cell is in the halo region, load the send buffers
                if (np > 1) {

                    double GhostGID = static_cast<double>(MyGrainID);
                    double GhostDOCX = static_cast<double>(GlobalX + 0.5);
                    double GhostDOCY = static_cast<double>(GlobalY + 0.5);
                    double GhostDOCZ = static_cast<double>(GlobalZ + 0.5);
                    double GhostDL = 0.01;
                    // Collect data for the ghost nodes, if necessary
                    loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, BufSizeX, MyYSlices, GlobalX,
                                   RankY, RankZ, AtNorthBoundary, AtSouthBoundary, BufferSouthSend, BufferNorthSend);

                } // End if statement for serial/parallel code
            }
        });
    Kokkos::fence();
}

//*****************************************************************************/
// Initializes cells for the current layer as either solid (don't resolidify) or tempsolid (will melt and resolidify)
void CellTypeInit_Remelt(int nx, int MyYSlices, int LocalActiveDomainSize, ViewI CellType, ViewI CritTimeStep, int id,
                         int ZBound_Low) {

    int MeltPoolCellCount;
    Kokkos::parallel_reduce(
        "CellTypeInitSolidRM", LocalActiveDomainSize,
        KOKKOS_LAMBDA(const int &D3D1ConvPosition, int &local_count) {
            int GlobalD3D1ConvPosition = D3D1ConvPosition + ZBound_Low * nx * MyYSlices;
            if (CritTimeStep(GlobalD3D1ConvPosition) != 0) {
                CellType(GlobalD3D1ConvPosition) = TempSolid;
                local_count++;
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
// Determine which nuclei are located on a given MPI rank, and may possibly occur during simulation
// Place the appropriate nuclei data into the structures NucleiLocation, NucleationTimes, and NucleiGrainID
// Case without remelting (each cell can only have 1 nuclei max, each cell solidifies at most, one time)
void placeNucleiData_NoRemelt(int Nuclei_ThisLayerSingle, ViewI_H NucleiX, ViewI_H NucleiY, ViewI_H NucleiZ,
                              int MyYOffset, int nx, int MyYSlices, bool AtNorthBoundary, bool AtSouthBoundary,
                              int ZBound_Low, ViewI_H CellType_Host, ViewI_H LayerID_Host, ViewI_H CritTimeStep_Host,
                              ViewF_H UndercoolingChange_Host, int layernumber,
                              std::vector<int> NucleiGrainID_WholeDomain_V,
                              std::vector<double> NucleiUndercooling_WholeDomain_V,
                              std::vector<int> &NucleiGrainID_MyRank_V, std::vector<int> &NucleiLocation_MyRank_V,
                              std::vector<int> &NucleationTimes_MyRank_V, int &PossibleNuclei_ThisRankThisLayer) {

    for (int NEvent = 0; NEvent < Nuclei_ThisLayerSingle; NEvent++) {
        if (((NucleiY(NEvent) > MyYOffset) || (AtSouthBoundary)) &&
            ((NucleiY(NEvent) < MyYOffset + MyYSlices - 1) || (AtNorthBoundary))) {
            // Convert 3D location (using global X and Y coordinates) into a 1D location (using local X and Y
            // coordinates) for the possible nucleation event, relative to the bottom of the overall domain
            int NucleiLocation_AllLayers = (NucleiZ(NEvent) + ZBound_Low) * nx * MyYSlices +
                                           NucleiX(NEvent) * MyYSlices + (NucleiY(NEvent) - MyYOffset);
            // Nucleus place criteria - cell is initially liquid, associated with the current layer of the problem
            if ((CellType_Host(NucleiLocation_AllLayers) == Liquid) &&
                (LayerID_Host(NucleiLocation_AllLayers) == layernumber)) {
                // Nucleation event is possible - cell is liquid and associated with this layer of the
                // multilayer problem - Add nuclei location and undercooling to the list for this rank
                NucleiLocation_MyRank_V[PossibleNuclei_ThisRankThisLayer] = NucleiLocation_AllLayers;
                int TimeToNucUnd =
                    CritTimeStep_Host(NucleiLocation_AllLayers) +
                    round(NucleiUndercooling_WholeDomain_V[NEvent] / UndercoolingChange_Host(NucleiLocation_AllLayers));
                NucleationTimes_MyRank_V[PossibleNuclei_ThisRankThisLayer] =
                    std::max(CritTimeStep_Host(NucleiLocation_AllLayers), TimeToNucUnd);
                // Assign this cell the potential nucleated grain ID
                NucleiGrainID_MyRank_V[PossibleNuclei_ThisRankThisLayer] = NucleiGrainID_WholeDomain_V[NEvent];
                // Increment counter on this MPI rank
                PossibleNuclei_ThisRankThisLayer++;
            }
        }
    }
}

// Determine which nuclei are located on a given MPI rank, and may possibly occur during simulation
// Place the appropriate nuclei data into the structures NucleiLocation, NucleationTimes, and NucleiGrainID
// Case with remelting (each cell can solidify multiple times, can be the home of multiple nucleation events)
void placeNucleiData_Remelt(int NucleiMultiplier, int Nuclei_ThisLayerSingle, ViewI_H NucleiX, ViewI_H NucleiY,
                            ViewI_H NucleiZ, int MyYOffset, int nx, int MyYSlices, bool AtNorthBoundary,
                            bool AtSouthBoundary, int ZBound_Low, ViewI_H NumberOfSolidificationEvents_Host,
                            ViewF3D_H LayerTimeTempHistory_Host, std::vector<int> NucleiGrainID_WholeDomain_V,
                            std::vector<double> NucleiUndercooling_WholeDomain_V,
                            std::vector<int> &NucleiGrainID_MyRank_V, std::vector<int> &NucleiLocation_MyRank_V,
                            std::vector<int> &NucleationTimes_MyRank_V, int &PossibleNuclei_ThisRankThisLayer) {

    for (int meltevent = 0; meltevent < NucleiMultiplier; meltevent++) {
        for (int n = 0; n < Nuclei_ThisLayerSingle; n++) {
            int NEvent = meltevent * Nuclei_ThisLayerSingle + n;
            if (((NucleiY(NEvent) > MyYOffset) || (AtSouthBoundary)) &&
                ((NucleiY(NEvent) < MyYOffset + MyYSlices - 1) || (AtNorthBoundary))) {
                // Convert 3D location (using global X and Y coordinates) into a 1D location (using local X and Y
                // coordinates) for the possible nucleation event, both as relative to the bottom of this layer as well
                // as relative to the bottom of the overall domain
                int NucleiLocation_ThisLayer =
                    NucleiZ(NEvent) * nx * MyYSlices + NucleiX(NEvent) * MyYSlices + (NucleiY(NEvent) - MyYOffset);
                int NucleiLocation_AllLayers = (NucleiZ(NEvent) + ZBound_Low) * nx * MyYSlices +
                                               NucleiX(NEvent) * MyYSlices + (NucleiY(NEvent) - MyYOffset);
                // Criteria for placing a nucleus - whether or not this nuclei is associated with a solidification event
                if (meltevent < NumberOfSolidificationEvents_Host(NucleiLocation_ThisLayer)) {
                    // Nucleation event is possible - cell undergoes solidification at least once, this nucleation
                    // event is associated with one of the time periods during which the associated cell undergoes
                    // solidification
                    NucleiLocation_MyRank_V[PossibleNuclei_ThisRankThisLayer] = NucleiLocation_AllLayers;

                    int CritTimeStep_ThisEvent = LayerTimeTempHistory_Host(NucleiLocation_ThisLayer, meltevent, 1);
                    float UndercoolingChange_ThisEvent =
                        LayerTimeTempHistory_Host(NucleiLocation_ThisLayer, meltevent, 2);
                    int TimeToNucUnd = CritTimeStep_ThisEvent +
                                       round(NucleiUndercooling_WholeDomain_V[NEvent] / UndercoolingChange_ThisEvent);
                    NucleationTimes_MyRank_V[PossibleNuclei_ThisRankThisLayer] =
                        std::max(CritTimeStep_ThisEvent, TimeToNucUnd);
                    // Assign this cell the potential nucleated grain ID
                    NucleiGrainID_MyRank_V[PossibleNuclei_ThisRankThisLayer] = NucleiGrainID_WholeDomain_V[NEvent];
                    // Increment counter on this MPI rank
                    PossibleNuclei_ThisRankThisLayer++;
                }
            }
        }
    }
}

// Initialize nucleation site locations, GrainID values, and time at which nucleation events will potentially occur
// Modified to include multiple possible nucleation events in cells that melt and solidify multiple times
void NucleiInit(int layernumber, double RNGSeed, int MyYSlices, int MyYOffset, int nx, int ny, int nzActive,
                int ZBound_Low, int id, double NMax, double dTN, double dTsigma, double deltax, ViewI &NucleiLocation,
                ViewI_H &NucleationTimes_Host, ViewI &NucleiGrainID, ViewI CellType, ViewI CritTimeStep,
                ViewF UndercoolingChange, ViewI LayerID, int &PossibleNuclei_ThisRankThisLayer, int &Nuclei_WholeDomain,
                bool AtNorthBoundary, bool AtSouthBoundary, bool RemeltingYN, int &NucleationCounter,
                ViewI &MaxSolidificationEvents, ViewI NumberOfSolidificationEvents, ViewF3D LayerTimeTempHistory) {

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
    // Multiplier for the number of nucleation events per layer, based on the number of solidification events
    int NucleiMultiplier;
    if (RemeltingYN)
        NucleiMultiplier = MaxSolidificationEvents_Host(layernumber);
    else
        NucleiMultiplier = 1;
    int Nuclei_ThisLayer = Nuclei_ThisLayerSingle * NucleiMultiplier;
    // Temporary vectors for storing nucleated grain IDs and undercooling values
    // Nuclei Grain ID are assigned to avoid reusing values from previous layers
    std::vector<int> NucleiGrainID_WholeDomain_V(Nuclei_ThisLayer);
    std::vector<double> NucleiUndercooling_WholeDomain_V(Nuclei_ThisLayer);
    // Views for storing potential nucleated grain coordinates
    ViewI_H NucleiX(Kokkos::ViewAllocateWithoutInitializing("NucleiX"), Nuclei_ThisLayer);
    ViewI_H NucleiY(Kokkos::ViewAllocateWithoutInitializing("NucleiY"), Nuclei_ThisLayer);
    ViewI_H NucleiZ(Kokkos::ViewAllocateWithoutInitializing("NucleiZ"), Nuclei_ThisLayer);
    for (int meltevent = 0; meltevent < NucleiMultiplier; meltevent++) {
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
    std::vector<int> NucleiGrainID_MyRank_V(Nuclei_ThisLayer), NucleiLocation_MyRank_V(Nuclei_ThisLayer),
        NucleationTimes_MyRank_V(Nuclei_ThisLayer);
    if (RemeltingYN)
        placeNucleiData_Remelt(NucleiMultiplier, Nuclei_ThisLayerSingle, NucleiX, NucleiY, NucleiZ, MyYOffset, nx,
                               MyYSlices, AtNorthBoundary, AtSouthBoundary, ZBound_Low,
                               NumberOfSolidificationEvents_Host, LayerTimeTempHistory_Host,
                               NucleiGrainID_WholeDomain_V, NucleiUndercooling_WholeDomain_V, NucleiGrainID_MyRank_V,
                               NucleiLocation_MyRank_V, NucleationTimes_MyRank_V, PossibleNuclei_ThisRankThisLayer);
    else
        placeNucleiData_NoRemelt(Nuclei_ThisLayerSingle, NucleiX, NucleiY, NucleiZ, MyYOffset, nx, MyYSlices,
                                 AtNorthBoundary, AtSouthBoundary, ZBound_Low, CellType_Host, LayerID_Host,
                                 CritTimeStep_Host, UndercoolingChange_Host, layernumber, NucleiGrainID_WholeDomain_V,
                                 NucleiUndercooling_WholeDomain_V, NucleiGrainID_MyRank_V, NucleiLocation_MyRank_V,
                                 NucleationTimes_MyRank_V, PossibleNuclei_ThisRankThisLayer);

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
    NucleiLocation = Kokkos::create_mirror_view_and_copy(device_memory_space(), NucleiLocation_Host);
    NucleiGrainID = Kokkos::create_mirror_view_and_copy(device_memory_space(), NucleiGrainID_Host);

    // Initialize counter for the layer to 0
    NucleationCounter = 0;

    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Nucleation data structures for layer " << layernumber << " initialized" << std::endl;
}

//*****************************************************************************/
void ZeroResetViews(int LocalActiveDomainSize, int BufSizeX, int BufSizeZ, ViewF &DiagonalLength,
                    ViewF &CritDiagonalLength, ViewF &DOCenter, Buffer2D &BufferNorthSend, Buffer2D &BufferSouthSend,
                    Buffer2D &BufferNorthRecv, Buffer2D &BufferSouthRecv, ViewI &SteeringVector) {

    // Realloc steering vector as LocalActiveDomainSize may have changed (old values aren't needed)
    Kokkos::realloc(SteeringVector, LocalActiveDomainSize);

    // Realloc active cell data structure and halo regions on device (old values not needed)
    Kokkos::realloc(DiagonalLength, LocalActiveDomainSize);
    Kokkos::realloc(DOCenter, 3 * LocalActiveDomainSize);
    Kokkos::realloc(CritDiagonalLength, 26 * LocalActiveDomainSize);
    Kokkos::realloc(BufferNorthSend, BufSizeX * BufSizeZ, 5);
    Kokkos::realloc(BufferSouthSend, BufSizeX * BufSizeZ, 5);
    Kokkos::realloc(BufferNorthRecv, BufSizeX * BufSizeZ, 5);
    Kokkos::realloc(BufferSouthRecv, BufSizeX * BufSizeZ, 5);

    // Reset active cell data structures on device
    Kokkos::deep_copy(DiagonalLength, 0);
    Kokkos::deep_copy(DOCenter, 0);
    Kokkos::deep_copy(CritDiagonalLength, 0);
    // Reset halo region structures on device
    Kokkos::deep_copy(BufferSouthSend, 0.0);
    Kokkos::deep_copy(BufferSouthRecv, 0.0);
    Kokkos::deep_copy(BufferNorthSend, 0.0);
    Kokkos::deep_copy(BufferNorthRecv, 0.0);
}
