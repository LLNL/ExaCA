// Copyright Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_INPUTDATA_HPP
#define EXACA_INPUTDATA_HPP

#include "CAinfo.hpp"
#include "CAtimers.hpp"

#include "mpi.h"

#include <nlohmann/json.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Error if this is not a valid simulation type.
inline void validSimulationType(std::string simulation_type) {
    if (simulation_type != "Directional" && simulation_type != "Spot" && simulation_type != "SingleGrain" &&
        simulation_type != "FromFile" && simulation_type != "FromFinch")
        throw std::runtime_error("Error: unknown problem type \"" + simulation_type + "\".");
}

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

#endif
