// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_INPUTS_HPP
#define EXACA_INPUTS_HPP

#include "CAinfo.hpp"
#include "CAparsefiles.hpp"
#include "CAtimers.hpp"

#include "mpi.h"

#include <Kokkos_Core.hpp>

#include <nlohmann/json.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Structs to organize data within inputs struct
struct DomainInputs {
    double deltax = 0.0, deltat = 0.0;
    // number of CA cells in each direction only initialized here for problem types C, Spot, and SingleGrain
    int nx = 0, ny = 0, nz = 0;
    // multilayer problems only
    int number_of_layers = 1, layer_height = 0;
    // problem type Spot only
    int spot_radius = 0;
};

struct NucleationInputs {
    // unused for single grain problem type, zero by default
    double n_max = 0.0;
    double dtn = 0.0;
    double dtsigma = 0.0;
};

struct InterfacialResponseInputs {
    float freezing_range;
    float A;
    float B;
    float C;
    float D = 0.0;
    enum IRFtypes {
        cubic = 0,
        quadratic = 1,
        power = 2,
    };
    int function = cubic;
};

struct TemperatureInputs {
    // Used for problem type R (by default, all temperature files read during init)
    bool layerwise_temp_read = false;
    int temp_files_in_series = 0;
    std::vector<std::string> temp_paths;
    // Use for problem types other than R (no temperature files to read) - default to no initial undercooling at
    // solidification front
    double G = 0, R = 0;
    double init_undercooling = 0.0;
};

struct SubstrateInputs {
    // problem type C only
    std::string surface_init_mode = "";
    // Only used for mode (i)
    double fract_surface_sites_active = 0.0;
    // Only used for mode (ii)
    double surface_site_density = 0.0;
    // Only used for mode (iii)
    std::vector<int> grain_locations_x, grain_locations_y, grain_ids;
    bool fill_bottom_surface = false;
    // problem type SingleGrain only
    int single_grain_orientation = 0;
    // problem types Spot and R only
    bool use_substrate_file = false;
    bool baseplate_through_powder = false;
    std::string substrate_filename = "";
    double substrate_grain_spacing = 0.0;
    // defaults to all sites in powder layer initialized with a new grain
    double powder_active_fraction = 1.0;
    // Top of baseplate assumed at Z = 0 if not otherwise given
    double baseplate_top_z = 0.0;
    // Initial size of octahedra during initialization of an active cell
    float init_oct_size = 0.01;
};

struct PrintInputs {
    // Base name of CA output
    std::string base_filename = "";
    // Path to CA output
    std::string path_to_output = "";
    // Names of output fields that can be printed to files during or at the end of a simulation
    std::vector<std::string> fieldnames_key = {"GrainID",
                                               "LayerID",
                                               "GrainMisorientation",
                                               "UndercoolingCurrent",
                                               "UndercoolingSolidificationStart",
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
    bool intralayer_undercooling_solidification_start = false;
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
    bool interlayer_undercooling_solidification_start = false;
    bool interlayer_undercooling_current = false;
    bool interlayer_melt_time_step = false;
    bool interlayer_crit_time_step = false;
    bool interlayer_undercooling_change = false;
    bool interlayer_cell_type = false;
    bool interlayer_diagonal_length = false;
    bool interlayer_solidification_event_counter = false;
    bool interlayer_number_of_solidification_events = false;
    // True if intralayer_undercooling_solidification_start or interlayer_undercooling_solidification_start is true
    bool store_solidification_start = false;
    bool print_front_undercooling = false;
    // List of layers following which the interlayer fields should be printed (will always include final layer of
    // simulation)
    std::vector<int> print_layer_number;

    // Should binary be used for printing vtk data?
    bool print_binary = false;

    // Should the default RVE data for ExaConstit be printed? If so, with what size?
    bool print_default_rve = false;
    int rve_size = 0;
};

// Error if this is not a valid simulation type.
inline void validSimulationType(std::string simulation_type) {
    if (simulation_type != "Directional" && simulation_type != "Spot" && simulation_type != "SingleGrain" &&
        simulation_type != "FromFile" && simulation_type != "FromFinch")
        throw std::runtime_error("Error: unknown problem type \"" + simulation_type + "\".");
}

struct Inputs {

    std::string simulation_type = "", material_filename = "", grain_orientation_file = "";
    unsigned long rng_seed = 0.0;
    DomainInputs domain;
    NucleationInputs nucleation;
    InterfacialResponseInputs irf;
    TemperatureInputs temperature;
    SubstrateInputs substrate;
    PrintInputs print;
    std::string file_name;

    // Creates input struct with uninitialized/default values, used in unit tests
    Inputs(){};

    Inputs(const int id, const std::string input_file)
        : file_name(input_file) {

        // Open and read JSON input file
        std::ifstream input_data_stream(input_file);
        nlohmann::json input_data = nlohmann::json::parse(input_data_stream);

        // General inputs
        simulation_type = input_data["SimulationType"];
        // "Directional": directional solidification
        // "Spot": hemispherical spot with fixed thermal gradient and cooling rate
        // "FromFile": time-temperature history comes from external files
        if (simulation_type == "C") {
            simulation_type = "Directional";
            std::cout << "Warning: Problem type \"C\" is now \"Directional\". Previous name will be removed in a "
                         "future release."
                      << std::endl;
        }
        else if (simulation_type == "RM" || simulation_type == "R") {
            simulation_type = "FromFile";
            std::cout
                << "Warning: Problem type \"R\" is now \"FromFile\". Previous name will be removed in a future release."
                << std::endl;
        }
        else if ((simulation_type == "S") || (simulation_type == "SM")) {
            throw std::runtime_error("Error: The spot melt array simulation type (Problem type S or SM) was removed "
                                     "after version 1.3; simulation of a single hemispherical spot can be performed "
                                     "using problem type Spot. See README for details");
        }

        // Check for valid simulation type.
        if (simulation_type == "FromFunch")
            simulation_type = "FromFinch";
        validSimulationType(simulation_type);

        // Input files that should be present for all problem types
        std::string material_filename_read = input_data["MaterialFileName"];
        std::string grain_orientation_file_read = input_data["GrainOrientationFile"];
        // Path to file of materials constants based on install/source location
        material_filename = checkFileInstalled(material_filename_read, id);
        checkFileNotEmpty(material_filename);
        // Path to file of grain orientations based on install/source location
        grain_orientation_file = checkFileInstalled(grain_orientation_file_read, id);
        checkFileNotEmpty(grain_orientation_file);
        // Seed for random number generator (defaults to 0 if not given)
        if (input_data.contains("RandomSeed"))
            rng_seed = input_data["RandomSeed"];

        // Domain inputs:
        // Cell size - given in meters, stored in micrometers
        domain.deltax = input_data["Domain"]["CellSize"];
        domain.deltax = domain.deltax * pow(10, -6);
        // Time step - given in seconds, stored in microseconds
        domain.deltat = input_data["Domain"]["TimeStep"];
        domain.deltat = domain.deltat * pow(10, -6);
        if ((simulation_type == "Directional") || (simulation_type == "SingleGrain")) {
            // Domain size, in cells
            domain.nx = input_data["Domain"]["Nx"];
            domain.ny = input_data["Domain"]["Ny"];
            domain.nz = input_data["Domain"]["Nz"];
        }
        else if ((simulation_type == "FromFile") || (simulation_type == "FromFinch")) {
            // Number of layers, layer height are needed for problem type and R
            domain.number_of_layers = input_data["Domain"]["NumberOfLayers"];
            domain.layer_height = input_data["Domain"]["LayerOffset"];
        }
        else if (simulation_type == "Spot") {
            std::vector<std::string> unused_spot_inputs = {"NSpotsX", "NSpotsY", "SpotOffset", "RSpots"};
            int num_unused_spot_inputs = unused_spot_inputs.size();
            for (int n = 0; n < num_unused_spot_inputs; n++)
                if ((input_data["Domain"].contains(unused_spot_inputs[n])) && (id == 0))
                    std::cout << "Note: Input " << unused_spot_inputs[n] << " is unused in the Spot problem type"
                              << std::endl;
            // Radius is given in micrometers, convert to cells
            domain.spot_radius = input_data["Domain"]["SpotRadius"];
            // Calculate nx, ny, and nz based on spot radius
            domain.nz = domain.spot_radius + 1;
            domain.nx = 2 * domain.spot_radius + 2;
            domain.ny = 2 * domain.spot_radius + 2;
        }

        // Nucleation inputs:
        // Nucleation density (normalized by 10^12 m^-3), mean nucleation undercooling/st dev undercooling(K)
        // SingleGrain problem type does not have nucleation, just a grain of a single orientation in the domain center
        if (simulation_type != "SingleGrain") {
            nucleation.n_max = input_data["Nucleation"]["Density"];
            nucleation.n_max = nucleation.n_max * pow(10, 12);
            nucleation.dtn = input_data["Nucleation"]["MeanUndercooling"];
            nucleation.dtsigma = input_data["Nucleation"]["StDev"];
        }

        parseIRF(id);

        // Temperature inputs:
        if (simulation_type == "FromFile") {
            if ((input_data["TemperatureData"].contains("HeatTransferCellSize")) && (id == 0))
                std::cout << "Note: Heat transport data cell size is no longer an input used in ExaCA, temperature "
                             "data must be at the same resolution as the CA cell size"
                          << std::endl;
            // Read all temperature files at once (default), or one at a time?
            if (input_data["TemperatureData"].contains("LayerwiseTempRead")) {
                temperature.layerwise_temp_read = input_data["TemperatureData"]["LayerwiseTempRead"];
            }
            // Get the paths/number of/names of the temperature data files used
            temperature.temp_files_in_series = input_data["TemperatureData"]["TemperatureFiles"].size();
            if (temperature.temp_files_in_series == 0)
                throw std::runtime_error("Error: No temperature files listed in the temperature instructions file");
            else {
                for (int filename = 0; filename < temperature.temp_files_in_series; filename++)
                    temperature.temp_paths.push_back(input_data["TemperatureData"]["TemperatureFiles"][filename]);
            }
        }
        else if (simulation_type != "FromFinch") {
            // Temperature data uses fixed thermal gradient (K/m) and cooling rate (K/s)
            temperature.G = input_data["TemperatureData"]["G"];
            temperature.R = input_data["TemperatureData"]["R"];
            // Optional initial undercooling for problem type C (required for type SingleGrain). Defaults to 0 for
            // problem type C if not given, should be a positive number
            if (input_data["TemperatureData"].contains("InitUndercooling"))
                if (input_data["TemperatureData"]["InitUndercooling"] < 0)
                    throw std::runtime_error("Error: optional temperature data argument InitUndercooling should be "
                                             "greater than or equal to zero");
            if (simulation_type == "SingleGrain")
                temperature.init_undercooling = input_data["TemperatureData"]["InitUndercooling"];
            else if ((simulation_type == "Directional") && (input_data["TemperatureData"].contains("InitUndercooling")))
                temperature.init_undercooling = input_data["TemperatureData"]["InitUndercooling"];
            if ((temperature.G > 0) && (Kokkos::fabs(temperature.R) < 0.000001)) {
                // Throw error for edge case where the cooling rate is 0, but cells in the domain would be initialized
                // above the liquidus temperature (i.e., cells that would never solidify)
                int location_init_undercooling;
                if (simulation_type == "Directional")
                    location_init_undercooling = 0;
                else
                    location_init_undercooling = Kokkos::floorf(static_cast<float>(domain.nz) / 2.0);
                int location_liquidus_isotherm =
                    location_init_undercooling +
                    Kokkos::round(temperature.init_undercooling / (temperature.G * domain.deltax));
                if ((temperature.R == 0) && (location_liquidus_isotherm <= domain.nz - 1))
                    throw std::runtime_error("Error: Domain will not fully solidify based on the size in Z, initial "
                                             "undercooling, thermal gradient, and cooling rate in the input file");
            }
        }

        // Substrate inputs:
        if (simulation_type == "Directional") {
            // Must contain inputs corresponding to one of the three modes
            // fract_surface_sites_active only used for mode (i)
            // surface_site_density only used for mode (ii) - given in grains/mm^2
            // grain_locations_x/grain_locations_y/grain_ids only used for mode (iii)
            // Determine surface initialization mode
            bool mode_i_input = input_data["Substrate"].contains("SurfaceSiteFraction");
            bool mode_ii_input = input_data["Substrate"].contains("SurfaceSiteDensity");
            bool mode_iii_input = ((input_data["Substrate"].contains("GrainLocationsX")) &&
                                   (input_data["Substrate"].contains("GrainLocationsY")) &&
                                   (input_data["Substrate"].contains("GrainIDs")));
            // Check for deprecated input FractionSurfaceSitesActive - which used the logic now associated with the
            // SurfaceSiteDensity input
            bool mode_deprecated_input = input_data["Substrate"].contains("FractionSurfaceSitesActive");
            if ((mode_i_input && mode_ii_input) || (mode_i_input && mode_iii_input) ||
                (mode_ii_input && mode_iii_input) || (mode_i_input && mode_ii_input && mode_iii_input))
                throw std::runtime_error(
                    "Error: Conflicting substrate input parameters detected. See examples/README for details");
            if (mode_i_input) {
                substrate.surface_init_mode = "SurfaceSiteFraction";
                substrate.fract_surface_sites_active = input_data["Substrate"]["SurfaceSiteFraction"];
            }
            else if (mode_ii_input) {
                substrate.surface_init_mode = "SurfaceSiteDensity";
                substrate.surface_site_density = input_data["Substrate"]["SurfaceSiteDensity"];
            }
            else if (mode_iii_input) {
                substrate.surface_init_mode = "Custom";
                // Ensure the number of grains is consistent
                std::size_t num_grains = input_data["Substrate"]["GrainLocationsX"].size();
                if ((input_data["Substrate"]["GrainLocationsY"].size() != num_grains) ||
                    (input_data["Substrate"]["GrainIDs"].size() != num_grains))
                    throw std::runtime_error(
                        "Error: GrainLocationsX, GrainLocationsY, and GrainIDs must be the same size");
                else {
                    for (std::size_t graincount = 0; graincount < num_grains; graincount++) {
                        substrate.grain_locations_x.push_back(input_data["Substrate"]["GrainLocationsX"][graincount]);
                        substrate.grain_locations_y.push_back(input_data["Substrate"]["GrainLocationsY"][graincount]);
                        substrate.grain_ids.push_back(input_data["Substrate"]["GrainIDs"][graincount]);
                    }
                }
            }
            else if (mode_deprecated_input) {
                if (id == 0)
                    std::cout
                        << "Warning: FractionSurfaceSitesActive is a deprecated input for problem type C and support "
                           "will be removed in a future release. See README for updated initialization options"
                        << std::endl;
                // The deprecated input uses the logic currently associated with the SurfaceSiteDensity initialization
                // mode, so the fraction should be transformed into a density
                substrate.surface_init_mode = "SurfaceSiteDensity";
                double fract_surface_active_old = input_data["Substrate"]["FractionSurfaceSitesActive"];
                substrate.surface_site_density =
                    fract_surface_active_old / (domain.deltax * domain.deltax * pow(10, 12));
            }
            else
                throw std::runtime_error("Error: either SurfaceSiteFraction, SurfaceSiteDensity, or lists of "
                                         "GrainLocationsX/GrainLocationsY/GrainIDs are "
                                         "required to initialize this problem type");
            if (input_data["Substrate"].contains("FillBottomSurface"))
                substrate.fill_bottom_surface = input_data["Substrate"]["FillBottomSurface"];
        }
        else if (simulation_type == "SingleGrain") {
            // Orientation of the single grain at the domain center
            substrate.single_grain_orientation = input_data["Substrate"]["GrainOrientation"];
        }
        else {
            // Substrate data - should data come from an initial size or a file?
            if ((input_data["Substrate"].contains("SubstrateFilename")) &&
                (input_data["Substrate"].contains("MeanSize")))
                throw std::runtime_error(
                    "Error: only one of substrate grain size and substrate structure filename should "
                    "be provided in the input file");
            else if (input_data["Substrate"].contains("SubstrateFilename")) {
                substrate.substrate_filename = input_data["Substrate"]["SubstrateFilename"];
                substrate.use_substrate_file = true;
            }
            else if (input_data["Substrate"].contains("MeanSize")) {
                substrate.substrate_grain_spacing = input_data["Substrate"]["MeanSize"];
                substrate.use_substrate_file = false;
            }
            if (simulation_type == "Spot") {
                // No powder for this problem type, baseplate top at Z = 0
                substrate.baseplate_top_z = 0.0;
                substrate.baseplate_through_powder = false;
            }
            else {
                // Should the baseplate microstructure be extended through the powder layers? Default is false
                if (input_data["Substrate"].contains("ExtendSubstrateThroughPower"))
                    substrate.baseplate_through_powder = input_data["Substrate"]["ExtendSubstrateThroughPower"];
                // defaults to a unique grain at each site in the powder layers if not given
                if (input_data["Substrate"].contains("PowderDensity")) {
                    // powder density is given as a density per unit volume, normalized by 10^12 m^-3 --> convert this
                    // into a density of sites active on the CA grid (0 to 1)
                    substrate.powder_active_fraction = input_data["Substrate"]["PowderDensity"];
                    substrate.powder_active_fraction =
                        substrate.powder_active_fraction * pow(10, 12) * pow(domain.deltax, 3);
                    if ((substrate.powder_active_fraction < 0.0) || (substrate.powder_active_fraction > 1.0))
                        throw std::runtime_error(
                            "Error: Density of powder surface sites active must be larger than 0 and less "
                            "than 1/(CA cell volume)");
                }
                if ((input_data["Substrate"].contains("PowderFirstLayer")) && (id == 0))
                    std::cout
                        << "Warning: PowderFirstLayer input is no longer used, the top of the first layer must be "
                           "specified using BaseplateTopZ (which will otherwise default to Z = 0)"
                        << std::endl;
                // The top of the baseplate is designated using BaseplateTopZ (assumed to be Z = 0 if not given in input
                // file)
                if (input_data["Substrate"].contains("BaseplateTopZ"))
                    substrate.baseplate_top_z = input_data["Substrate"]["BaseplateTopZ"];
                if ((substrate.baseplate_through_powder) && (input_data["Substrate"].contains("PowderDensity")))
                    throw std::runtime_error(
                        "Error: if the option to extend the baseplate through the powder layers is "
                        "toggled, a powder layer density cannot be given");
            }
            // Optional input for initial size of octhedra when a cell begins solidification (in units of CA cells)
            if (input_data["Substrate"].contains("InitOctahedronSize")) {
                substrate.init_oct_size = input_data["Substrate"]["InitOctahedronSize"];
                if ((substrate.init_oct_size < 0.0) || (substrate.init_oct_size >= 1.0))
                    throw std::runtime_error("Error: InitOctahedronSize should be at least 0, and less than 1");
            }
        }

        // Printing inputs
        getPrintDataFromInputFile(input_data, id);

        // Print information to console about the input file data read
        if (id == 0) {
            std::cout << "Material simulated is " << material_filename << std::endl;
            std::cout << "CA cell size is " << domain.deltax * pow(10, 6) << " microns" << std::endl;
            std::cout << "Nucleation density is " << nucleation.n_max << " per m^3" << std::endl;
            std::cout << "Mean nucleation undercooling is " << nucleation.dtn
                      << " K, standard deviation of distribution is " << nucleation.dtsigma << "K" << std::endl;
            if (simulation_type == "Directional") {
                std::cout << "CA Simulation using a unidirectional, fixed thermal gradient of " << temperature.G
                          << " K/m and a cooling rate of " << temperature.R << " K/s" << std::endl;
                std::cout << "The time step is " << domain.deltat * pow(10, 6) << " microseconds" << std::endl;
                if (substrate.surface_init_mode == "Custom")
                    std::cout << "Input grain locations and ID values will be used at the bottom surface" << std::endl;
                else if (substrate.surface_init_mode == "SurfaceSiteFraction")
                    std::cout << "The fraction of CA cells at the bottom surface that are active is "
                              << substrate.fract_surface_sites_active << std::endl;
                else if (substrate.surface_init_mode == "SurfaceSiteDensity")
                    std::cout << "The density of active CA cells at the bottom surface is "
                              << substrate.surface_site_density << " grains per µm2" << std::endl;
            }
            else if (simulation_type == "Spot") {
                std::cout << "CA Simulation using a radial, fixed thermal gradient of " << temperature.G
                          << " K/m across a hemispherical spot of radius " << domain.spot_radius
                          << " cells, with a cooling rate of " << temperature.R << " K/s" << std::endl;
                std::cout << "The time step is " << domain.deltat * pow(10, 6) << " microseconds" << std::endl;
            }
            else if (simulation_type == "FromFile") {
                std::cout << "CA Simulation using temperature data from file(s)" << std::endl;
                std::cout << "The time step is " << domain.deltat << " seconds" << std::endl;
                std::cout << "The first temperature data file to be read is " << temperature.temp_paths[0]
                          << ", and there are " << temperature.temp_files_in_series << " in the series" << std::endl;
                std::cout << "A total of " << domain.number_of_layers << " layers of solidification offset by "
                          << domain.layer_height << " CA cells will be simulated" << std::endl;
            }
        }
    }

    // Updated version, where fields are organized into intralayer and interlayer
    std::vector<bool> getPrintFieldValues(nlohmann::json input_data, const std::string fieldtype,
                                          const std::vector<std::string> fieldnames_key) {
        int num_fields_key = fieldnames_key.size();
        int num_fields_given = input_data["Printing"][fieldtype]["Fields"].size();
        std::vector<bool> print_fields_given(num_fields_key, false);
        // Check each given field against each possible input field name
        for (int field_given = 0; field_given < num_fields_given; field_given++) {
            for (int field_key = 0; field_key < num_fields_key; field_key++) {
                if (input_data["Printing"][fieldtype]["Fields"][field_given] == fieldnames_key[field_key])
                    print_fields_given[field_key] = true;
            }
        }
        return print_fields_given;
    }

    // Read the input data file and initialize appropriate variables to non-default values if necessary
    void getPrintDataFromInputFile(nlohmann::json input_data, const int id) {
        // Path to output data
        print.path_to_output = input_data["Printing"]["PathToOutput"];
        // Name of output data
        print.base_filename = input_data["Printing"]["OutputFile"];
        if (simulation_type == "Directional")
            if (input_data["Printing"].contains("PrintFrontUndercooling"))
                print.print_front_undercooling = input_data["Printing"]["PrintFrontUndercooling"];
        // Should ASCII or binary be used to print vtk data? Defaults to ASCII if not given
        if (input_data["Printing"].contains("PrintBinary"))
            print.print_binary = input_data["Printing"]["PrintBinary"];
        // Should default ExaConstit output be printed after the simulation? If so, what size RVE?
        // If a size of 0 is given, this is set to false
        if (input_data["Printing"].contains("PrintExaConstitSize")) {
            print.rve_size = input_data["Printing"]["PrintExaConstitSize"];
            if (print.rve_size != 0)
                print.print_default_rve = true;
        }

        if ((input_data["Printing"].contains("PrintFieldsInit")) ||
            (input_data["Printing"].contains("PrintIntermediateOutput")) ||
            (input_data["Printing"].contains("PrintIntermediateOutput"))) {
            throw std::runtime_error("Error: The old print input format is no longer compatible with ExaCA. See "
                                     "examples/README for updated format");
        }
        else {
            if (input_data["Printing"].contains("Intralayer")) {
                // Fields to be printed during a simulation or during a layer of a simulation - if increment
                if (input_data["Printing"]["Intralayer"]["Increment"] == 0) {
                    // Files only printed for the initial state of a given layer
                    print.intralayer_increment = INT_MAX;
                }
                else {
                    // Files to be printed at some interval during each layer
                    print.intralayer_increment = input_data["Printing"]["Intralayer"]["Increment"];
                    print.intralayer_idle_frames = input_data["Printing"]["Intralayer"]["PrintIdleFrames"];
                }
                // Which fields should be printed during the layers?
                std::vector<bool> print_fields_intralayer =
                    getPrintFieldValues(input_data, "Intralayer", print.fieldnames_key);
                if (print_fields_intralayer[0])
                    print.intralayer_grain_id = true;
                if (print_fields_intralayer[1])
                    print.intralayer_layer_id = true;
                if (print_fields_intralayer[2])
                    print.intralayer_grain_misorientation = true;
                if (print_fields_intralayer[3])
                    print.intralayer_undercooling_current = true;
                if (print_fields_intralayer[4])
                    print.intralayer_undercooling_solidification_start = true;
                if (print_fields_intralayer[5])
                    print.intralayer_melt_time_step = true;
                if (print_fields_intralayer[6])
                    print.intralayer_crit_time_step = true;
                if (print_fields_intralayer[7])
                    print.intralayer_undercooling_change = true;
                if (print_fields_intralayer[8])
                    print.intralayer_cell_type = true;
                if (print_fields_intralayer[9])
                    print.intralayer_diagonal_length = true;
                if (print_fields_intralayer[10])
                    print.intralayer_solidification_event_counter = true;
                if (print_fields_intralayer[11])
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
            if (input_data["Printing"]["Interlayer"].contains("Layers")) {
                if ((id == 0) && (input_data["Printing"]["Interlayer"].contains("Increment")))
                    std::cout << "Warning: A list of layers to print and a layer increment were both present in the "
                                 "input file print options, the layer increment will be ignored"
                              << std::endl;
                int num_print_layers = input_data["Printing"]["Interlayer"]["Layers"].size();
                for (int n = 0; n < num_print_layers; n++) {
                    int print_layer_val = input_data["Printing"]["Interlayer"]["Layers"][n];
                    if (print_layer_val >= domain.number_of_layers) {
                        if (id == 0)
                            std::cout << "Note: adjusting layer value of " << print_layer_val << " to "
                                      << domain.number_of_layers - 1
                                      << " as the simulation only contains layers 0 through " << domain.number_of_layers
                                      << std::endl;
                        print_layer_val = domain.number_of_layers - 1;
                    }
                    print.print_layer_number.push_back(print_layer_val);
                }
                // Make sure files print after last layer, even if it wasn't listed
                if (print.print_layer_number[num_print_layers - 1] != domain.number_of_layers - 1)
                    print.print_layer_number.push_back(domain.number_of_layers - 1);
            }
            else if (input_data["Printing"]["Interlayer"].contains("Increment")) {
                // Print layer numbers starting at 0 and at interlayer_increment, always including the last layer
                int interlayer_increment = input_data["Printing"]["Interlayer"]["Increment"];
                for (int n = 0; n < domain.number_of_layers - 1; n += interlayer_increment)
                    print.print_layer_number.push_back(n);
                print.print_layer_number.push_back(domain.number_of_layers - 1);
            }
            else
                print.print_layer_number.push_back(domain.number_of_layers - 1);

            // Which fields should be printed during the layers?
            std::vector<bool> print_fields_interlayer =
                getPrintFieldValues(input_data, "Interlayer", print.fieldnames_key);
            if (print_fields_interlayer[0])
                print.interlayer_grain_id = true;
            if (print_fields_interlayer[1])
                print.interlayer_layer_id = true;
            if (print_fields_interlayer[2])
                print.interlayer_grain_misorientation = true;
            if (print_fields_interlayer[3])
                print.interlayer_undercooling_current = true;
            if (print_fields_interlayer[4])
                print.interlayer_undercooling_solidification_start = true;
            if (print_fields_interlayer[5])
                print.interlayer_melt_time_step = true;
            if (print_fields_interlayer[6])
                print.interlayer_crit_time_step = true;
            if (print_fields_interlayer[7])
                print.interlayer_undercooling_change = true;
            if (print_fields_interlayer[8])
                print.interlayer_cell_type = true;
            if (print_fields_interlayer[9])
                print.interlayer_diagonal_length = true;
            if (print_fields_interlayer[10])
                print.interlayer_solidification_event_counter = true;
            if (print_fields_interlayer[11])
                print.interlayer_number_of_solidification_events = true;
            if ((print.interlayer_grain_id) || (print.interlayer_layer_id) || (print.interlayer_undercooling_current) ||
                (print.interlayer_undercooling_solidification_start))
                print.interlayer_full = true;
            // First 5 inputs are full domain inputs - check if any of the others were toggled
            int num_interlayer_current_inputs = print_fields_interlayer.size();
            for (int n = 5; n < num_interlayer_current_inputs; n++) {
                if (print_fields_interlayer[n])
                    print.interlayer_current = true;
            }
            // Should starting undercooling for solidification be stored?
            if ((print.intralayer_undercooling_solidification_start) ||
                (print.interlayer_undercooling_solidification_start) || (print.print_front_undercooling))
                print.store_solidification_start = true;
        }
        if (id == 0)
            std::cout << "Successfully parsed data printing options from input file" << std::endl;
    }

    // Ensure that input powder layer init options are compatible with this domain size, if needed for this problem type
    // TODO: Expand check that inputs are valid for the problem type
    void checkPowderOverflow(const int nx, const int ny, const int layer_height, const int number_of_layers) {
        // Check to make sure powder grain density is compatible with the number of powder sites
        // If this problem type includes a powder layer of some grain density, ensure that integer overflow won't occur
        // when assigning powder layer GrainIDs
        if (!(substrate.baseplate_through_powder)) {
            long int num_cells_powder_layers = static_cast<long int>(nx) * static_cast<long int>(ny) *
                                               static_cast<long int>(layer_height) *
                                               static_cast<long int>(number_of_layers - 1);
            long int num_assigned_cells_powder_layers = std::lround(
                Kokkos::round(static_cast<double>(num_cells_powder_layers) * substrate.powder_active_fraction));
            if (num_assigned_cells_powder_layers > INT_MAX)
                throw std::runtime_error(
                    "Error: A smaller value for powder density is required to avoid potential integer "
                    "overflow when assigning powder layer GrainID");
        }
    }

    // Print a log file for this ExaCA run in json file format, containing information about the run parameters used
    // from the input file as well as the decomposition scheme
    // Note: Passing external values for inputs like deltax that will later be stored in the grid class, with the grid
    // class passed to this function
    void printExaCALog(const int id, const int np, const int ny_local, const int y_offset, const double deltax,
                       const int number_of_layers, const int layer_height, const int nx, const int ny, const int nz,
                       Timers timers, const int cycle, const double x_min, const double x_max, const double y_min,
                       const double y_max, const double z_min, const double z_max, const float vol_fraction_nucleated) {

        int *ny_local_allranks = new int[np];
        int *y_offset_allranks = new int[np];
        MPI_Gather(&ny_local, 1, MPI_INT, ny_local_allranks, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(&y_offset, 1, MPI_INT, y_offset_allranks, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (id == 0) {
            std::string FName = print.path_to_output + print.base_filename + ".json";
            std::cout << "Printing ExaCA log file" << std::endl;
            std::ofstream exaca_log;
            exaca_log.open(FName);
            exaca_log << "{" << std::endl;
            exaca_log << "   \"ExaCAVersion\": \"" << version() << "\", " << std::endl;
            exaca_log << "   \"ExaCACommitHash\": \"" << gitCommitHash() << "\", " << std::endl;
            exaca_log << "   \"KokkosVersion\": \"" << kokkosVersion() << "\", " << std::endl;
            exaca_log << "   \"InputFile\": \"" << file_name << "\", " << std::endl;
            exaca_log << "   \"TimeStepOfOutput\": " << cycle << "," << std::endl;
            exaca_log << "   \"SimulationType\": \"" << simulation_type << "\"," << std::endl;
            exaca_log << "   \"GrainOrientationFile\": \"" << grain_orientation_file << "\"," << std::endl;
            exaca_log << "   \"Domain\": {" << std::endl;
            exaca_log << "      \"Nx\": " << nx << "," << std::endl;
            exaca_log << "      \"Ny\": " << ny << "," << std::endl;
            exaca_log << "      \"Nz\": " << nz << "," << std::endl;
            exaca_log << "      \"CellSize\": " << deltax << "," << std::endl;
            exaca_log << "      \"TimeStep\": " << domain.deltat << "," << std::endl;
            exaca_log << "      \"XBounds\": [" << x_min << "," << x_max << "]," << std::endl;
            exaca_log << "      \"YBounds\": [" << y_min << "," << y_max << "]," << std::endl;
            exaca_log << "      \"ZBounds\": [" << z_min << "," << z_max << "]";
            if (simulation_type == "FromFile") {
                exaca_log << "," << std::endl;
                exaca_log << "      \"NumberOfLayers\": " << number_of_layers << "," << std::endl;
                exaca_log << "      \"LayerOffset\": " << layer_height;
            }
            else if (simulation_type == "Spot") {
                exaca_log << "," << std::endl;
                exaca_log << "      \"SpotRadius\": " << domain.spot_radius;
            }
            exaca_log << std::endl;
            exaca_log << "   }," << std::endl;
            exaca_log << "   \"Nucleation\": {" << std::endl;
            exaca_log << "      \"Density\": " << nucleation.n_max << "," << std::endl;
            exaca_log << "      \"MeanUndercooling\": " << nucleation.dtn << "," << std::endl;
            exaca_log << "      \"StDevUndercooling\": " << nucleation.dtsigma << "," << std::endl;
            exaca_log << "      \"VolFractionNucleated\": " << vol_fraction_nucleated << std::endl;
            exaca_log << "   }," << std::endl;
            exaca_log << "   \"TemperatureData\": {" << std::endl;
            if (simulation_type == "FromFile") {
                exaca_log << "       \"TemperatureFiles\": [";
                for (int i = 0; i < temperature.temp_files_in_series - 1; i++) {
                    exaca_log << "\"" << temperature.temp_paths[i] << "\", ";
                }
                exaca_log << "\"" << temperature.temp_paths[temperature.temp_files_in_series - 1] << "\"]" << std::endl;
            }
            else {
                exaca_log << "      \"G\": " << temperature.G << "," << std::endl;
                exaca_log << "      \"R\": " << temperature.R << std::endl;
            }
            exaca_log << "   }," << std::endl;
            exaca_log << "   \"Substrate\": {" << std::endl;
            if (simulation_type == "Directional")
                exaca_log << "       \"SurfaceSiteFraction\": " << substrate.fract_surface_sites_active << std::endl;
            else if (simulation_type == "SingleGrain")
                exaca_log << "       \"GrainOrientation\": " << substrate.single_grain_orientation << std::endl;
            else {
                if (substrate.use_substrate_file)
                    exaca_log << "       \"SubstrateFilename\": " << substrate.substrate_filename << std::endl;
                else
                    exaca_log << "       \"MeanSize\": " << substrate.substrate_grain_spacing << std::endl;
            }
            exaca_log << "   }," << std::endl;
            exaca_log << "   \"InterfacialResponse\": {" << std::endl;
            exaca_log << "       \"Function\": "
                      << "\"" << irf.function << "\"," << std::endl;
            exaca_log << "       \"A\": " << (irf.A) << "," << std::endl;
            exaca_log << "       \"B\": " << (irf.B) << "," << std::endl;
            exaca_log << "       \"C\": " << (irf.C) << "," << std::endl;
            if (irf.function == irf.cubic)
                exaca_log << "       \"D\": " << (irf.D) << "," << std::endl;
            exaca_log << "       \"FreezingRange\": " << (irf.freezing_range) << std::endl;
            exaca_log << "   },";
            exaca_log << "   \"NumberMPIRanks\": " << np << "," << std::endl;
            exaca_log << "   \"Decomposition\": {" << std::endl;
            exaca_log << "       \"SubdomainYSize\": [";
            for (int i = 0; i < np - 1; i++)
                exaca_log << ny_local_allranks[i] << ",";
            exaca_log << ny_local_allranks[np - 1] << "]," << std::endl;
            exaca_log << "       \"SubdomainYOffset\": [";
            for (int i = 0; i < np - 1; i++)
                exaca_log << y_offset_allranks[i] << ",";
            exaca_log << y_offset_allranks[np - 1] << "]" << std::endl;
            exaca_log << "   }," << std::endl;
            exaca_log << timers.printLog() << std::endl;
            exaca_log << "}" << std::endl;
            exaca_log.close();
        }
    }

    void parseIRF(const int id) {
        if (id == 0)
            std::cout << "Parsing material file using json input format" << std::endl;
        std::ifstream material_data(material_filename);
        nlohmann::json data = nlohmann::json::parse(material_data);
        irf.A = data["coefficients"]["A"];
        irf.B = data["coefficients"]["B"];
        irf.C = data["coefficients"]["C"];
        std::string functionform = data["function"];
        if (functionform == "cubic") {
            irf.D = data["coefficients"]["D"];
            irf.function = irf.cubic;
        }
        else if ((functionform == "quadratic") || (functionform == "power")) {
            // D should not have been given, this functional form only takes 3 input fitting parameters
            if (data["coefficients"]["D"] != nullptr) {
                std::string error = "Error: functional form of this type takes only A, B, and C as inputs";
                throw std::runtime_error(error);
            }
            if (functionform == "quadratic")
                irf.function = irf.quadratic;
            else if (functionform == "power")
                irf.function = irf.power;
        }
        else
            throw std::runtime_error("Error: Unrecognized functional form for interfacial response function, currently "
                                     "supported options are quadratic, cubic, and exponential");
        irf.freezing_range = data["freezing_range"];
    }
};

#endif
