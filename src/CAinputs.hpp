// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_INPUTS_HPP
#define EXACA_INPUTS_HPP

#include "CAinfo.hpp"
#include "CAinterfacialresponse.hpp"

#include "mpi.h"

#include <nlohmann/json.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Structs to organize data within inputs struct
struct DomainInputs {
    double deltax = 0.0, deltat = 0.0;
    // number of CA cells in each direction only initialized here for problem types C, S, and SingleGrain
    int nx = 0, ny = 0, nz = 0;
    // multilayer problems only
    int NumberOfLayers = 1, LayerHeight = 0;
    // problem type S only
    int NSpotsX = 0, NSpotsY = 0, SpotRadius = 0, SpotOffset = 0;
};

struct NucleationInputs {
    // unused for single grain problem type, zero by default
    double NMax = 0.0;
    double dTN = 0.0;
    double dTsigma = 0.0;
};

struct TemperatureInputs {
    // Used for problem type R (by default, all temperature files read during init)
    bool LayerwiseTempRead = false;
    int TempFilesInSeries = 0;
    std::vector<std::string> temp_paths;
    // Use for problem types other than R (no temperature files to read) - default to no initial undercooling at
    // solidification front
    double G = 0, R = 0;
    double initUndercooling = 0.0;
};

struct SubstrateInputs {
    // problem type C only - one of these two inputs will be given
    double FractSurfaceSitesActive = 0.0;
    std::vector<int> GrainLocationsX, GrainLocationsY, GrainIDs;
    bool CustomGrainLocationsIDs = false;
    bool FillBottomSurface = false;
    // problem type SingleGrain only
    int singleGrainOrientation = 0;
    // problem types S and R only
    bool UseSubstrateFile = false;
    bool BaseplateThroughPowder = false;
    std::string SubstrateFileName = "";
    double SubstrateGrainSpacing = 0.0;
    // defaults to all sites in powder layer initialized with a new grain
    double PowderActiveFraction = 1.0;
    // Top of baseplate assumed at Z = 0 if not otherwise given
    double BaseplateTopZ = 0.0;
};

struct PrintInputs {
    // Base name of CA output
    std::string BaseFileName = "";
    // Path to CA output
    std::string PathToOutput = "";
    // Names of output fields that can be printed to files during or at the end of a simulation
    std::vector<std::string> Fieldnames_key = {"GrainID",
                                               "LayerID",
                                               "GrainMisorientation",
                                               "UndercoolingCurrent",
                                               "MeltTimeStep",
                                               "CritTimeStep",
                                               "UndercoolingChange",
                                               "CellType",
                                               "DiagonalLength",
                                               "SolidificationEventCounter",
                                               "NumberOfSolidificationEvents"};
    // Fields to be printed during a given layer
    bool intralayer = false;
    int intralayer_increment = 1;
    bool intralayer_idle_frames = false;
    bool intralayer_grain_id = false;
    bool intralayer_layer_id = false;
    bool intralayer_grain_misorientation = false;
    bool intralayer_undercooling_current = false;
    bool intralayer_melt_time_step = false;
    bool intralayer_crit_time_step = false;
    bool intralayer_undercooling_change = false;
    bool intralayer_cell_type = false;
    bool intralayer_diagonal_length = false;
    bool intralayer_solidification_event_counter = false;
    bool intralayer_number_of_solidification_events = false;
    // Fields to be printed at end of a layer
    bool interlayer_full = false;
    bool interlayer_current = false;
    bool interlayer_grain_id = false;
    bool interlayer_layer_id = false;
    bool interlayer_grain_misorientation = false;
    bool interlayer_undercooling_current = false;
    bool interlayer_melt_time_step = false;
    bool interlayer_crit_time_step = false;
    bool interlayer_undercooling_change = false;
    bool interlayer_cell_type = false;
    bool interlayer_diagonal_length = false;
    bool interlayer_solidification_event_counter = false;
    bool interlayer_number_of_solidification_events = false;
    // List of layers following which the interlayer fields should be printed (will always include final layer of
    // simulation)
    std::vector<int> print_layer_number;

    // Should binary be used for printing vtk data?
    bool PrintBinary = false;

    // Should the default RVE data for ExaConstit be printed? If so, with what size?
    bool PrintDefaultRVE = false;
    int RVESize = 0;
};

struct Inputs {

    std::string SimulationType = "", MaterialFileName = "", GrainOrientationFile = "";
    unsigned long RNGSeed = 0.0;
    DomainInputs domain;
    NucleationInputs nucleation;
    TemperatureInputs temperature;
    SubstrateInputs substrate;
    PrintInputs print;

    // Creates input struct with uninitialized/default values, used in unit tests
    Inputs(){};

    Inputs(int id, std::string InputFile) {

        // Open and read JSON input file
        std::ifstream InputData(InputFile);
        nlohmann::json inputdata = nlohmann::json::parse(InputData);

        // General inputs
        SimulationType = inputdata["SimulationType"];
        // "C": constrained (directional) solidification
        // "S": array of overlapping hemispherical spots
        // "R": time-temperature history comes from external files
        // Check if simulation type includes remelting ("M" suffix to input problem type) - all simulations now use
        // remelting, so in the absence of this suffix, print warning that the problem will use remelting
        // DirSoldification problems now include remelting logic
        if (SimulationType == "RM") {
            SimulationType = "R";
        }
        else if (SimulationType == "SM") {
            SimulationType = "S";
        }
        else if ((SimulationType == "S") || (SimulationType == "R")) {
            if (id == 0) {
                if (SimulationType == "S")
                    std::cout
                        << "Warning: While the specified problem type did not include remelting, all simulations now "
                           "include remelting"
                        << std::endl;
            }
        }
        if ((SimulationType == "S") && (id == 0))
            std::cout << "Warning: The spot melt array simulation type (Problem type S) is now deprecated and will be "
                         "removed in a future release"
                      << std::endl;

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
        if (inputdata.contains("RandomSeed"))
            RNGSeed = inputdata["RandomSeed"];

        // Domain inputs:
        // Cell size - given in meters, stored in micrometers
        domain.deltax = inputdata["Domain"]["CellSize"];
        domain.deltax = domain.deltax * pow(10, -6);
        // Time step - given in seconds, stored in microseconds
        domain.deltat = inputdata["Domain"]["TimeStep"];
        domain.deltat = domain.deltat * pow(10, -6);
        if ((SimulationType == "C") || (SimulationType == "SingleGrain")) {
            // Domain size, in cells
            domain.nx = inputdata["Domain"]["Nx"];
            domain.ny = inputdata["Domain"]["Ny"];
            domain.nz = inputdata["Domain"]["Nz"];
        }
        else {
            // Number of layers, layer height are needed for problem types S and R
            domain.NumberOfLayers = inputdata["Domain"]["NumberOfLayers"];
            domain.LayerHeight = inputdata["Domain"]["LayerOffset"];
            // Type S needs spot information, which is then used to compute the domain bounds
            if (SimulationType == "S") {
                domain.NSpotsX = inputdata["Domain"]["NSpotsX"];
                domain.NSpotsY = inputdata["Domain"]["NSpotsY"];
                // Radius and offset are given in micrometers, convert to cells
                domain.SpotRadius = inputdata["Domain"]["RSpots"];
                domain.SpotRadius = domain.SpotRadius * pow(10, -6) / domain.deltax;
                domain.SpotOffset = inputdata["Domain"]["SpotOffset"];
                domain.SpotOffset = domain.SpotOffset * pow(10, -6) / domain.deltax;
                // Calculate nx, ny, and nz based on spot array pattern and number of layers
                domain.nz = domain.SpotRadius + 1 + (domain.NumberOfLayers - 1) * domain.LayerHeight;
                domain.nx = 2 * domain.SpotRadius + 1 + domain.SpotOffset * (domain.NSpotsX - 1);
                domain.ny = 2 * domain.SpotRadius + 1 + domain.SpotOffset * (domain.NSpotsY - 1);
            }
        }

        // Nucleation inputs:
        // Nucleation density (normalized by 10^12 m^-3), mean nucleation undercooling/st dev undercooling(K)
        // SingleGrain problem type does not have nucleation, just a grain of a single orientation in the domain center
        if (SimulationType != "SingleGrain") {
            nucleation.NMax = inputdata["Nucleation"]["Density"];
            nucleation.NMax = nucleation.NMax * pow(10, 12);
            nucleation.dTN = inputdata["Nucleation"]["MeanUndercooling"];
            nucleation.dTsigma = inputdata["Nucleation"]["StDev"];
        }

        // Temperature inputs:
        if (SimulationType == "R") {
            if ((inputdata["TemperatureData"].contains("HeatTransferCellSize")) && (id == 0))
                std::cout << "Note: Heat transport data cell size is no longer an input used in ExaCA, temperature "
                             "data must be at the same resolution as the CA cell size"
                          << std::endl;
            // Read all temperature files at once (default), or one at a time?
            if (inputdata["TemperatureData"].contains("LayerwiseTempRead")) {
                temperature.LayerwiseTempRead = inputdata["TemperatureData"]["LayerwiseTempRead"];
            }
            // Get the paths/number of/names of the temperature data files used
            temperature.TempFilesInSeries = inputdata["TemperatureData"]["TemperatureFiles"].size();
            if (temperature.TempFilesInSeries == 0)
                throw std::runtime_error("Error: No temperature files listed in the temperature instructions file");
            else {
                for (int filename = 0; filename < temperature.TempFilesInSeries; filename++)
                    temperature.temp_paths.push_back(inputdata["TemperatureData"]["TemperatureFiles"][filename]);
            }
        }
        else {
            // Temperature data uses fixed thermal gradient (K/m) and cooling rate (K/s)
            temperature.G = inputdata["TemperatureData"]["G"];
            temperature.R = inputdata["TemperatureData"]["R"];
            // Optional initial undercooling for problem type C (required for type SingleGrain). Defaults to 0 for
            // problem type C if not given, should be a positive number
            if (inputdata["TemperatureData"].contains("InitUndercooling"))
                if (inputdata["TemperatureData"]["InitUndercooling"] < 0)
                    throw std::runtime_error("Error: optional temperature data argument InitUndercooling should be "
                                             "greater than or equal to zero");
            if (SimulationType == "SingleGrain")
                temperature.initUndercooling = inputdata["TemperatureData"]["InitUndercooling"];
            else if ((SimulationType == "C") && (inputdata["TemperatureData"].contains("InitUndercooling")))
                temperature.initUndercooling = inputdata["TemperatureData"]["InitUndercooling"];
            if ((temperature.G > 0) && (Kokkos::fabs(temperature.R) < 0.000001)) {
                // Throw error for edge case where the cooling rate is 0, but cells in the domain would be initialized
                // above the liquidus temperature (i.e., cells that would never solidify)
                int location_init_undercooling;
                if (SimulationType == "C")
                    location_init_undercooling = 0;
                else
                    location_init_undercooling = Kokkos::floorf(static_cast<float>(domain.nz) / 2.0);
                int location_liquidus_isotherm =
                    location_init_undercooling +
                    Kokkos::round(temperature.initUndercooling / (temperature.G * domain.deltax));
                if ((temperature.R == 0) && (location_liquidus_isotherm <= domain.nz - 1))
                    throw std::runtime_error("Error: Domain will not fully solidify based on the size in Z, initial "
                                             "undercooling, thermal gradient, and cooling rate in the input file");
            }
        }

        // Substrate inputs:
        if (SimulationType == "C") {
            // Must contain FractionSurfaceSitesActive OR all of GrainLocationsX/GrainLocationsY/GrainIDs, not both
            bool contains_FractionSurfaceSitesActive = inputdata["Substrate"].contains("FractionSurfaceSitesActive");
            bool contains_GrainLocationsX = inputdata["Substrate"].contains("GrainLocationsX");
            bool contains_GrainLocationsY = inputdata["Substrate"].contains("GrainLocationsY");
            bool contains_GrainIDs = inputdata["Substrate"].contains("GrainIDs");
            bool contains_CustomGrainInitInfo =
                contains_GrainLocationsX && contains_GrainLocationsY && contains_GrainIDs;
            if (contains_FractionSurfaceSitesActive) {
                if ((id == 0) && (contains_CustomGrainInitInfo))
                    std::cout << "Warning: FractionSurfaceSitesActive will be ignored as explicit grain locations and "
                                 "IDs are given in the input file"
                              << std::endl;
                else
                    substrate.FractSurfaceSitesActive = inputdata["Substrate"]["FractionSurfaceSitesActive"];
            }
            else if (contains_CustomGrainInitInfo) {
                // Ensure the number of grains is consistent
                std::size_t numGrains = inputdata["Substrate"]["GrainLocationsX"].size();
                if ((inputdata["Substrate"]["GrainLocationsY"].size() != numGrains) ||
                    (inputdata["Substrate"]["GrainIDs"].size() != numGrains))
                    throw std::runtime_error(
                        "Error: GrainLocationsX, GrainLocationsY, and GrainIDs must be the same size");
                else {
                    for (std::size_t graincount = 0; graincount < numGrains; graincount++) {
                        substrate.GrainLocationsX.push_back(inputdata["Substrate"]["GrainLocationsX"][graincount]);
                        substrate.GrainLocationsY.push_back(inputdata["Substrate"]["GrainLocationsY"][graincount]);
                        substrate.GrainIDs.push_back(inputdata["Substrate"]["GrainIDs"][graincount]);
                    }
                    substrate.CustomGrainLocationsIDs = true;
                }
            }
            else
                throw std::runtime_error(
                    "Error: either FractionSurfaceSitesActive or lists of GrainLocationsX/GrainLocationsY/GrainIDs are "
                    "required to initialize this problem type");
            if (inputdata["Substrate"].contains("FillBottomSurface"))
                substrate.FillBottomSurface = inputdata["Substrate"]["FillBottomSurface"];
        }
        else if (SimulationType == "SingleGrain") {
            // Orientation of the single grain at the domain center
            substrate.singleGrainOrientation = inputdata["Substrate"]["GrainOrientation"];
        }
        else {
            // Substrate data - should data come from an initial size or a file?
            if ((inputdata["Substrate"].contains("SubstrateFilename")) && (inputdata["Substrate"].contains("MeanSize")))
                throw std::runtime_error(
                    "Error: only one of substrate grain size and substrate structure filename should "
                    "be provided in the input file");
            else if (inputdata["Substrate"].contains("SubstrateFilename")) {
                substrate.SubstrateFileName = inputdata["Substrate"]["SubstrateFilename"];
                substrate.UseSubstrateFile = true;
            }
            else if (inputdata["Substrate"].contains("MeanSize")) {
                substrate.SubstrateGrainSpacing = inputdata["Substrate"]["MeanSize"];
                substrate.UseSubstrateFile = false;
            }
            // Should the baseplate microstructure be extended through the powder layers? Default is false
            if (inputdata["Substrate"].contains("ExtendSubstrateThroughPower"))
                substrate.BaseplateThroughPowder = inputdata["Substrate"]["ExtendSubstrateThroughPower"];
            // defaults to a unique grain at each site in the powder layers if not given
            if (inputdata["Substrate"].contains("PowderDensity")) {
                // powder density is given as a density per unit volume, normalized by 10^12 m^-3 --> convert this into
                // a density of sites active on the CA grid (0 to 1)
                substrate.PowderActiveFraction = inputdata["Substrate"]["PowderDensity"];
                substrate.PowderActiveFraction = substrate.PowderActiveFraction * pow(10, 12) * pow(domain.deltax, 3);
                if ((substrate.PowderActiveFraction < 0.0) || (substrate.PowderActiveFraction > 1.0))
                    throw std::runtime_error(
                        "Error: Density of powder surface sites active must be larger than 0 and less "
                        "than 1/(CA cell volume)");
            }
            if ((inputdata["Substrate"].contains("PowderFirstLayer")) && (id == 0))
                std::cout << "Warning: PowderFirstLayer input is no longer used, the top of the first layer must be "
                             "specified using BaseplateTopZ (which will otherwise default to Z = 0)"
                          << std::endl;
            // The top of the baseplate is designated using BaseplateTopZ (assumed to be Z = 0 if not given in input
            // file)
            if (inputdata["Substrate"].contains("BaseplateTopZ"))
                substrate.BaseplateTopZ = inputdata["Substrate"]["BaseplateTopZ"];
            else if (SimulationType == "S") {
                // Baseplate top if not otherwise given for spot array problem is at the top of the first layer
                substrate.BaseplateTopZ = domain.deltax * domain.SpotRadius;
            }
            if ((substrate.BaseplateThroughPowder) && (inputdata["Substrate"].contains("PowderDensity")))
                throw std::runtime_error("Error: if the option to extend the baseplate through the powder layers is "
                                         "toggled, a powder layer density cannot be given");
        }

        // Printing inputs
        getPrintDataFromInputFile(inputdata, id, domain.deltat);

        // Print information to console about the input file data read
        if (id == 0) {
            std::cout << "Material simulated is " << MaterialFileName << std::endl;
            std::cout << "CA cell size is " << domain.deltax * pow(10, 6) << " microns" << std::endl;
            std::cout << "Nucleation density is " << nucleation.NMax << " per m^3" << std::endl;
            std::cout << "Mean nucleation undercooling is " << nucleation.dTN
                      << " K, standard deviation of distribution is " << nucleation.dTsigma << "K" << std::endl;
            if (SimulationType == "C") {
                std::cout << "CA Simulation using a unidirectional, fixed thermal gradient of " << temperature.G
                          << " K/m and a cooling rate of " << temperature.R << " K/s" << std::endl;
                std::cout << "The time step is " << domain.deltat * pow(10, 6) << " microseconds" << std::endl;
                if (substrate.CustomGrainLocationsIDs)
                    std::cout << "Input grain locations and ID values will be used at the bottom surface" << std::endl;
                else
                    std::cout << "The fraction of CA cells at the bottom surface that are active is "
                              << substrate.FractSurfaceSitesActive << std::endl;
            }
            else if (SimulationType == "S") {
                std::cout << "CA Simulation using a radial, fixed thermal gradient of " << temperature.G
                          << " K/m as a series of hemispherical spots, and a cooling rate of " << temperature.R
                          << " K/s" << std::endl;
                std::cout << "A total of " << domain.NumberOfLayers << " spots per layer, with layers offset by "
                          << domain.LayerHeight << " CA cells will be simulated" << std::endl;
                std::cout << "The time step is " << domain.deltat * pow(10, 6) << " microseconds" << std::endl;
            }
            else if (SimulationType == "R") {
                std::cout << "CA Simulation using temperature data from file(s)" << std::endl;
                std::cout << "The time step is " << domain.deltat << " seconds" << std::endl;
                std::cout << "The first temperature data file to be read is " << temperature.temp_paths[0]
                          << ", and there are " << temperature.TempFilesInSeries << " in the series" << std::endl;
                std::cout << "A total of " << domain.NumberOfLayers << " layers of solidification offset by "
                          << domain.LayerHeight << " CA cells will be simulated" << std::endl;
            }
        }
    }

    // Using the old print inputs, read the input data file and initialize appropriate variables to non-default values
    // if necessary
    void getPrintDataFromInputFile_Old(nlohmann::json inputdata, int id, double deltat) {
        // Which fields should be printed at the start and end of the simulation?
        std::vector<std::string> InitFieldnames_key = {"GrainID", "LayerID", "MeltTimeStep", "CritTimeStep",
                                                       "UndercoolingChange"};
        std::vector<bool> PrintFieldsInit = getPrintFieldValues_Old(inputdata, "PrintFieldsInit", InitFieldnames_key);
        if (PrintFieldsInit[0])
            print.intralayer_grain_id = true;
        if (PrintFieldsInit[1])
            print.intralayer_layer_id = true;
        if (PrintFieldsInit[2])
            print.intralayer_melt_time_step = true;
        if (PrintFieldsInit[3])
            print.intralayer_crit_time_step = true;
        if (PrintFieldsInit[4])
            print.intralayer_undercooling_change = true;
        int num_print_intralayer_inputs = PrintFieldsInit.size();
        // True if any fields are printed
        for (int n = 0; n < num_print_intralayer_inputs; n++) {
            if (PrintFieldsInit[n])
                print.intralayer = true;
        }

        std::vector<std::string> FinalFieldnames_key = {
            "GrainID",      "LayerID",      "GrainMisorientation", "UndercoolingCurrent",
            "MeltTimeStep", "CritTimeStep", "UndercoolingChange",  "CellType"};
        std::vector<bool> PrintFieldsFinal =
            getPrintFieldValues_Old(inputdata, "PrintFieldsFinal", FinalFieldnames_key);
        if (PrintFieldsFinal[0])
            print.interlayer_grain_id = true;
        if (PrintFieldsFinal[1])
            print.interlayer_layer_id = true;
        if (PrintFieldsFinal[2])
            print.interlayer_grain_misorientation = true;
        if (PrintFieldsFinal[3])
            print.interlayer_undercooling_current = true;
        if (PrintFieldsFinal[4])
            print.interlayer_melt_time_step = true;
        if (PrintFieldsFinal[5])
            print.interlayer_crit_time_step = true;
        if (PrintFieldsFinal[6])
            print.interlayer_undercooling_change = true;
        if (PrintFieldsFinal[7])
            print.interlayer_cell_type = true;
        // Print fields final - only after last layer of simulation. Set increment to the number of layers
        print.print_layer_number.push_back(domain.NumberOfLayers - 1);
        if ((print.interlayer_grain_id) || (print.interlayer_layer_id))
            print.interlayer_full = true;
        for (int n = 3; n < 7; n++) {
            if (PrintFieldsInit[n])
                print.interlayer_current = true;
        }

        // Initial fields should be printed, set intralayer increment to large number so only the first frame will be
        // printed (this is overwritten if PrintIntermediateOutput and a Frequency are toggled)
        if (print.intralayer) {
            print.intralayer = true;
            print.intralayer_increment = INT_MAX;
        }
        // Should intermediate output be printed?
        if (inputdata["Printing"].contains("PrintIntermediateOutput")) {
            // An increment of 0 will set the intermediate file printing to false
            double TimeSeriesFrameInc_time = inputdata["Printing"]["PrintIntermediateOutput"]["Frequency"];
            if (TimeSeriesFrameInc_time != 0) {
                print.intralayer = true;
                // Increment is given in microseconds, convert to seconds
                TimeSeriesFrameInc_time = TimeSeriesFrameInc_time * pow(10, -6);
                // Overwrite INT_MAX value for intralayer_increment if PrintInitFields was toggled
                print.intralayer_increment = Kokkos::round(TimeSeriesFrameInc_time / deltat);
                // Should the intermediate output be printed even if the simulation was unchanged from the previous
                // output step?
                print.intralayer_idle_frames = inputdata["Printing"]["PrintIntermediateOutput"]["PrintIdleFrames"];
                if (id == 0)
                    std::cout << "Intermediate output for movie frames will be printed every "
                              << print.intralayer_increment << " time steps (or every "
                              << print.intralayer_increment * deltat << " microseconds)" << std::endl;
            }
        }
    }

    // Read the input data file and initialize appropriate variables to non-default values if necessary
    void getPrintDataFromInputFile(nlohmann::json inputdata, int id, double deltat) {
        // Path to output data
        print.PathToOutput = inputdata["Printing"]["PathToOutput"];
        // Name of output data
        print.BaseFileName = inputdata["Printing"]["OutputFile"];
        // Should ASCII or binary be used to print vtk data? Defaults to ASCII if not given
        if (inputdata["Printing"].contains("PrintBinary"))
            print.PrintBinary = inputdata["Printing"]["PrintBinary"];
        // Should default ExaConstit output be printed after the simulation? If so, what size RVE?
        // If a size of 0 is given, this is set to false
        if (inputdata["Printing"].contains("PrintExaConstitSize")) {
            print.RVESize = inputdata["Printing"]["PrintExaConstitSize"];
            if (print.RVESize != 0)
                print.PrintDefaultRVE = true;
        }

        if ((inputdata["Printing"].contains("PrintFieldsInit")) ||
            (inputdata["Printing"].contains("PrintIntermediateOutput")) ||
            (inputdata["Printing"].contains("PrintIntermediateOutput"))) {
            if (id == 0)
                std::cout << "Warning: Old print input format detected; compatibility with this will be removed in a "
                             "future release. See examples/README for updated format"
                          << std::endl;
            getPrintDataFromInputFile_Old(inputdata, id, deltat);
        }
        else {
            if (inputdata["Printing"].contains("Intralayer")) {
                // Fields to be printed during a simulation or during a layer of a simulation - if increment
                print.intralayer = true;
                if (inputdata["Printing"]["Intralayer"]["Increment"] == 0) {
                    // Files only printed for the initial state of a given layer
                    print.intralayer_increment = INT_MAX;
                }
                else {
                    // Files to be printed at some interval during each layer
                    print.intralayer_increment = inputdata["Printing"]["Intralayer"]["Increment"];
                    print.intralayer_idle_frames = inputdata["Printing"]["Intralayer"]["PrintIdleFrames"];
                }
                // Which fields should be printed during the layers?
                std::vector<bool> print_fields_intralayer =
                    getPrintFieldValues(inputdata, "Intralayer", print.Fieldnames_key);
                if (print_fields_intralayer[0])
                    print.intralayer_grain_id = true;
                if (print_fields_intralayer[1])
                    print.intralayer_layer_id = true;
                if (print_fields_intralayer[2])
                    print.intralayer_grain_misorientation = true;
                if (print_fields_intralayer[3])
                    print.intralayer_undercooling_current = true;
                if (print_fields_intralayer[4])
                    print.intralayer_melt_time_step = true;
                if (print_fields_intralayer[5])
                    print.intralayer_crit_time_step = true;
                if (print_fields_intralayer[6])
                    print.intralayer_undercooling_change = true;
                if (print_fields_intralayer[7])
                    print.intralayer_cell_type = true;
                if (print_fields_intralayer[8])
                    print.intralayer_diagonal_length = true;
                if (print_fields_intralayer[9])
                    print.intralayer_solidification_event_counter = true;
                if (print_fields_intralayer[10])
                    print.intralayer_number_of_solidification_events = true;
                // True if any fields are printed
                int num_print_intralayer_inputs = print_fields_intralayer.size();
                for (int n = 0; n < num_print_intralayer_inputs; n++) {
                    if (print_fields_intralayer[n])
                        print.intralayer = true;
                }
            }
            // List of layers following which interlayer data should be printed (will always print after last layer by
            // default)
            if (inputdata["Printing"]["Interlayer"].contains("Layers")) {
                if ((id == 0) && (inputdata["Printing"]["Interlayer"].contains("Increment")))
                    std::cout << "Warning: A list of layers to print and a layer increment were both present in the "
                                 "input file print options, the layer increment will be ignored"
                              << std::endl;
                int num_print_layers = inputdata["Printing"]["Interlayer"]["Layers"].size();
                for (int n = 0; n < num_print_layers; n++) {
                    int print_layer_val = inputdata["Printing"]["Interlayer"]["Layers"][n];
                    if (print_layer_val >= domain.NumberOfLayers) {
                        if (id == 0)
                            std::cout << "Note: adjusting layer value of " << print_layer_val << " to "
                                      << domain.NumberOfLayers - 1
                                      << " as the simulation only contains layers 0 through " << domain.NumberOfLayers
                                      << std::endl;
                        print_layer_val = domain.NumberOfLayers - 1;
                    }
                    print.print_layer_number.push_back(print_layer_val);
                }
                // Make sure files print after last layer, even if it wasn't listed
                if (print.print_layer_number[num_print_layers - 1] != domain.NumberOfLayers - 1)
                    print.print_layer_number.push_back(domain.NumberOfLayers - 1);
            }
            else if (inputdata["Printing"]["Interlayer"].contains("Increment")) {
                // Print layer numbers starting at 0 and at interlayer_increment, always including the last layer
                int interlayer_increment = inputdata["Printing"]["Interlayer"]["Increment"];
                for (int n = 0; n < domain.NumberOfLayers - 1; n += interlayer_increment)
                    print.print_layer_number.push_back(n);
                print.print_layer_number.push_back(domain.NumberOfLayers - 1);
            }
            else
                print.print_layer_number.push_back(domain.NumberOfLayers - 1);

            // Which fields should be printed during the layers?
            std::vector<bool> print_fields_interlayer =
                getPrintFieldValues(inputdata, "Interlayer", print.Fieldnames_key);
            if (print_fields_interlayer[0])
                print.interlayer_grain_id = true;
            if (print_fields_interlayer[1])
                print.interlayer_layer_id = true;
            if (print_fields_interlayer[2])
                print.interlayer_grain_misorientation = true;
            if (print_fields_interlayer[3])
                print.interlayer_undercooling_current = true;
            if (print_fields_interlayer[4])
                print.interlayer_melt_time_step = true;
            if (print_fields_interlayer[5])
                print.interlayer_crit_time_step = true;
            if (print_fields_interlayer[6])
                print.interlayer_undercooling_change = true;
            if (print_fields_interlayer[7])
                print.interlayer_cell_type = true;
            if (print_fields_interlayer[8])
                print.interlayer_diagonal_length = true;
            if (print_fields_interlayer[9])
                print.interlayer_solidification_event_counter = true;
            if (print_fields_interlayer[10])
                print.interlayer_number_of_solidification_events = true;
            if ((print.interlayer_grain_id) || (print.interlayer_layer_id) || (print.interlayer_undercooling_current))
                print.interlayer_full = true;
            // First 4 inputs are full domain inputs - check if any of the others were toggled
            int num_interlayer_current_inputs = print_fields_interlayer.size();
            for (int n = 4; n < num_interlayer_current_inputs; n++) {
                if (print_fields_interlayer[n])
                    print.interlayer_current = true;
            }
        }
        if (id == 0)
            std::cout << "Successfully parsed data printing options from input file" << std::endl;
    }

    // Ensure that input powder layer init options are compatible with this domain size, if needed for this problem type
    // TODO: Expand check that inputs are valid for the problem type
    void checkPowderOverflow(int nx, int ny, int LayerHeight, int NumberOfLayers) {
        // Check to make sure powder grain density is compatible with the number of powder sites
        // If this problem type includes a powder layer of some grain density, ensure that integer overflow won't occur
        // when assigning powder layer GrainIDs
        if (!(substrate.BaseplateThroughPowder)) {
            long int NumCellsPowderLayers =
                (long int)(nx) * (long int)(ny) * (long int)(LayerHeight) * (long int)(NumberOfLayers - 1);
            long int NumAssignedCellsPowderLayers =
                std::lround(Kokkos::round(static_cast<double>(NumCellsPowderLayers) * substrate.PowderActiveFraction));
            if (NumAssignedCellsPowderLayers > INT_MAX)
                throw std::runtime_error(
                    "Error: A smaller value for powder density is required to avoid potential integer "
                    "overflow when assigning powder layer GrainID");
        }
    }

    // Print a log file for this ExaCA run in json file format, containing information about the run parameters used
    // from the input file as well as the decomposition scheme
    // Note: Passing external values for inputs like deltax that will later be stored in the grid class, with the grid
    // class passed to this function
    void PrintExaCALog(int id, int np, std::string InputFile, int ny_local, int y_offset,
                       InterfacialResponseFunction irf, double deltax, int NumberOfLayers, int LayerHeight, int nx,
                       int ny, int nz, double InitTime, double RunTime, double OutTime, int cycle, double InitMaxTime,
                       double InitMinTime, double NuclMaxTime, double NuclMinTime, double CreateSVMinTime,
                       double CreateSVMaxTime, double CaptureMaxTime, double CaptureMinTime, double GhostMaxTime,
                       double GhostMinTime, double OutMaxTime, double OutMinTime, double XMin, double XMax, double YMin,
                       double YMax, double ZMin, double ZMax, float VolFractionNucleated) {

        int *ny_local_allranks = new int[np];
        int *y_offset_allranks = new int[np];
        MPI_Gather(&ny_local, 1, MPI_INT, ny_local_allranks, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(&y_offset, 1, MPI_INT, y_offset_allranks, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (id == 0) {
            std::string FName = print.PathToOutput + print.BaseFileName + ".json";
            std::cout << "Printing ExaCA log file" << std::endl;
            std::ofstream ExaCALog;
            ExaCALog.open(FName);
            ExaCALog << "{" << std::endl;
            ExaCALog << "   \"ExaCAVersion\": \"" << version() << "\", " << std::endl;
            ExaCALog << "   \"ExaCACommitHash\": \"" << gitCommitHash() << "\", " << std::endl;
            ExaCALog << "   \"KokkosVersion\": \"" << kokkosVersion() << "\", " << std::endl;
            ExaCALog << "   \"InputFile\": \"" << InputFile << "\", " << std::endl;
            ExaCALog << "   \"TimeStepOfOutput\": " << cycle << "," << std::endl;
            ExaCALog << "   \"SimulationType\": \"" << SimulationType << "\"," << std::endl;
            ExaCALog << "   \"GrainOrientationFile\": \"" << GrainOrientationFile << "\"," << std::endl;
            ExaCALog << "   \"Domain\": {" << std::endl;
            ExaCALog << "      \"Nx\": " << nx << "," << std::endl;
            ExaCALog << "      \"Ny\": " << ny << "," << std::endl;
            ExaCALog << "      \"Nz\": " << nz << "," << std::endl;
            ExaCALog << "      \"CellSize\": " << deltax << "," << std::endl;
            ExaCALog << "      \"TimeStep\": " << domain.deltat << "," << std::endl;
            ExaCALog << "      \"XBounds\": [" << XMin << "," << XMax << "]," << std::endl;
            ExaCALog << "      \"YBounds\": [" << YMin << "," << YMax << "]," << std::endl;
            ExaCALog << "      \"ZBounds\": [" << ZMin << "," << ZMax << "]";
            if (SimulationType != "C") {
                ExaCALog << "," << std::endl;
                ExaCALog << "      \"NumberOfLayers\": " << NumberOfLayers << "," << std::endl;
                ExaCALog << "      \"LayerOffset\": " << LayerHeight;
                if (SimulationType == "S") {
                    ExaCALog << "," << std::endl;
                    ExaCALog << "      \"NSpotsX\": " << domain.NSpotsX << "," << std::endl;
                    ExaCALog << "      \"NSpotsY\": " << domain.NSpotsY << "," << std::endl;
                    ExaCALog << "      \"RSpots\": " << domain.SpotRadius << "," << std::endl;
                    ExaCALog << "      \"SpotOffset\": " << domain.SpotOffset << std::endl;
                }
                else
                    ExaCALog << std::endl;
            }
            else
                ExaCALog << std::endl;
            ExaCALog << "   }," << std::endl;
            ExaCALog << "   \"Nucleation\": {" << std::endl;
            ExaCALog << "      \"Density\": " << nucleation.NMax << "," << std::endl;
            ExaCALog << "      \"MeanUndercooling\": " << nucleation.dTN << "," << std::endl;
            ExaCALog << "      \"StDevUndercooling\": " << nucleation.dTsigma << "," << std::endl;
            ExaCALog << "      \"VolFractionNucleated\": " << VolFractionNucleated << std::endl;
            ExaCALog << "   }," << std::endl;
            ExaCALog << "   \"TemperatureData\": {" << std::endl;
            if (SimulationType == "R") {
                ExaCALog << "       \"TemperatureFiles\": [";
                for (int i = 0; i < temperature.TempFilesInSeries - 1; i++) {
                    ExaCALog << "\"" << temperature.temp_paths[i] << "\", ";
                }
                ExaCALog << "\"" << temperature.temp_paths[temperature.TempFilesInSeries - 1] << "\"]" << std::endl;
            }
            else {
                ExaCALog << "      \"G\": " << temperature.G << "," << std::endl;
                ExaCALog << "      \"R\": " << temperature.R << std::endl;
            }
            ExaCALog << "   }," << std::endl;
            ExaCALog << "   \"Substrate\": {" << std::endl;
            if (SimulationType == "C")
                ExaCALog << "       \"FractionSurfaceSitesActive\": " << substrate.FractSurfaceSitesActive << std::endl;
            else if (SimulationType == "SingleGrain")
                ExaCALog << "       \"GrainOrientation\": " << substrate.singleGrainOrientation << std::endl;
            else {
                if (substrate.UseSubstrateFile)
                    ExaCALog << "       \"SubstrateFilename\": " << substrate.SubstrateFileName << std::endl;
                else
                    ExaCALog << "       \"MeanSize\": " << substrate.SubstrateGrainSpacing << std::endl;
            }
            ExaCALog << "   }," << std::endl;
            ExaCALog << irf.print() << std::endl;
            ExaCALog << "   \"NumberMPIRanks\": " << np << "," << std::endl;
            ExaCALog << "   \"Decomposition\": {" << std::endl;
            ExaCALog << "       \"SubdomainYSize\": [";
            for (int i = 0; i < np - 1; i++)
                ExaCALog << ny_local_allranks[i] << ",";
            ExaCALog << ny_local_allranks[np - 1] << "]," << std::endl;
            ExaCALog << "       \"SubdomainYOffset\": [";
            for (int i = 0; i < np - 1; i++)
                ExaCALog << y_offset_allranks[i] << ",";
            ExaCALog << y_offset_allranks[np - 1] << "]" << std::endl;
            ExaCALog << "   }," << std::endl;
            ExaCALog << "   \"Timing\": {" << std::endl;
            ExaCALog << "       \"Runtime\": " << InitTime + RunTime + OutTime << "," << std::endl;
            ExaCALog << "       \"InitRunOutputBreakdown\": [" << InitTime << "," << RunTime << "," << OutTime << "],"
                     << std::endl;
            ExaCALog << "       \"MaxMinInitTime\": [" << InitMaxTime << "," << InitMinTime << "]," << std::endl;
            ExaCALog << "       \"MaxMinNucleationTime\": [" << NuclMaxTime << "," << NuclMinTime << "]," << std::endl;
            ExaCALog << "       \"MaxMinSteeringVectorCreationTime\": [" << CreateSVMaxTime << "," << CreateSVMinTime
                     << "]," << std::endl;
            ExaCALog << "       \"MaxMinCellCaptureTime\": [" << CaptureMaxTime << "," << CaptureMinTime << "],"
                     << std::endl;
            ExaCALog << "       \"MaxMinGhostExchangeTime\": [" << GhostMaxTime << "," << GhostMinTime << "],"
                     << std::endl;
            ExaCALog << "       \"MaxMinOutputTime\": [" << OutMaxTime << "," << OutMinTime << "]" << std::endl;
            ExaCALog << "   }" << std::endl;
            ExaCALog << "}" << std::endl;
            ExaCALog.close();
        }
    }
};

#endif
