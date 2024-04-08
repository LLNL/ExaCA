// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_TEMPS_HPP
#define EXACA_TEMPS_HPP

#include "CAgrid.hpp"
#include "CAinputs.hpp"
#include "CAparsefiles.hpp"
#include "mpi.h"

#include <Kokkos_Core.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Reduced form of the time-temperature history and temperature field variables used by ExaCA
template <typename MemorySpace>
struct Temperature {

    using memory_space = MemorySpace;
    using view_type_int = Kokkos::View<int *, memory_space>;
    using view_type_float = Kokkos::View<float *, memory_space>;
    using view_type_double = Kokkos::View<double *, memory_space>;
    using view_type_double_2d = Kokkos::View<double **, memory_space>;
    using view_type_float_2d = Kokkos::View<float **, memory_space>;
    using view_type_float_3d = Kokkos::View<float ***, memory_space>;
    using view_type_int_host = typename view_type_int::HostMirror;
    using view_type_float_host = typename view_type_float::HostMirror;
    using view_type_double_host = typename view_type_double::HostMirror;
    using view_type_double_2d_host = typename view_type_double_2d::HostMirror;
    using view_type_float_2d_host = typename view_type_float_2d::HostMirror;
    using view_type_float_3d_host = typename view_type_float_3d::HostMirror;

    // Using the default exec space for this memory space.
    using execution_space = typename memory_space::execution_space;

    // Maximum number of times each CA cell in a given layer undergoes solidification
    view_type_int max_solidification_events;
    // For each cell in the current layer (index 1), and each time solidification happens (index 2), hold the values
    // that will be used for MeltTimeStep, CritTimeStep, and UndercoolingChange (index 3)
    view_type_float_3d layer_time_temp_history;
    // The number of times that each CA cell will undergo solidification during this layer
    view_type_int number_of_solidification_events;
    // A counter for the number of times each CA cell has undergone solidification so far this layer
    view_type_int solidification_event_counter;
    // The current undercooling of a CA cell (if superheated liquid or hasn't undergone solidification yet, equals 0)
    // Also maintained for the full multilayer domain
    view_type_float undercooling_current_all_layers, undercooling_current;
    // Data structure for storing raw temperature data from file(s)
    // Store data as double - needed for small time steps to resolve local differences in solidification conditions
    // Each data point has 6 values (X, Y, Z coordinates, melting time, liquidus time, and cooling rate)
    const int num_temperature_components = 6;
    view_type_double_2d_host raw_temperature_data;
    // These contain "number_of_layers" values corresponding to the location within "raw_temperature_data" of the first
    // data element in each temperature file, if used
    view_type_int_host first_value, last_value;
    // Undercooling when a solidification first began in a cell (optional to store based on selected inputs)
    bool _store_solidification_start;
    view_type_float undercooling_solidification_start_all_layers, undercooling_solidification_start;
    // Temperature field inputs from file
    TemperatureInputs _inputs;

    // Constructor creates views with size based on the grid inputs - each cell assumed to solidify once by default,
    // layer_time_temp_history modified to account for multiple events if needed undercooling_current and
    // solidification_event_counter are default initialized to zeros
    Temperature(const Grid &grid, TemperatureInputs inputs, const bool store_solidification_start = false,
                const int est_num_temperature_data_points = 100000)
        : max_solidification_events(
              view_type_int(Kokkos::ViewAllocateWithoutInitializing("number_of_layers"), grid.number_of_layers))
        , layer_time_temp_history(view_type_float_3d(Kokkos::ViewAllocateWithoutInitializing("layer_time_temp_history"),
                                                     grid.domain_size, 1, 3))
        , number_of_solidification_events(view_type_int(
              Kokkos::ViewAllocateWithoutInitializing("number_of_solidification_events"), grid.domain_size))
        , solidification_event_counter(view_type_int("solidification_event_counter", grid.domain_size))
        , undercooling_current_all_layers(view_type_float("undercooling_current", grid.domain_size_all_layers))
        , raw_temperature_data(view_type_double_2d_host(Kokkos::ViewAllocateWithoutInitializing("raw_temperature_data"),
                                                        est_num_temperature_data_points, num_temperature_components))
        , first_value(view_type_int_host(Kokkos::ViewAllocateWithoutInitializing("first_value"), grid.number_of_layers))
        , last_value(view_type_int_host(Kokkos::ViewAllocateWithoutInitializing("last_value"), grid.number_of_layers))
        , _store_solidification_start(store_solidification_start)
        , _inputs(inputs) {

        if (_store_solidification_start) {
            // Default init starting undercooling in cells to zero
            undercooling_solidification_start_all_layers =
                view_type_float("undercooling_solidification_start", grid.domain_size_all_layers);
            getCurrentLayerStartingUndercooling(grid.layer_range);
        }
        getCurrentLayerUndercooling(grid.layer_range);
    }

    // Read and parse the temperature file (double precision values in a comma-separated, ASCII format with a header
    // line - or a binary string of double precision values), storing the x, y, z, tm, tl, cr values in the RawData
    // vector. Each rank only contains the points corresponding to cells within the associated Y bounds.
    // number_of_temperature_data_points is incremented on each rank as data is added to RawData
    void parseTemperatureData(const std::string tempfile_thislayer, const double y_min, const double deltax,
                              const int lower_y_bound, const int upper_y_bound, int &number_of_temperature_data_points,
                              const bool binary_input_data, const int temperature_buffer_increment = 100000) {

        std::ifstream temperature_filestream;
        temperature_filestream.open(tempfile_thislayer);
        if (binary_input_data) {
            while (!temperature_filestream.eof()) {
                double x_temperature_point = readBinaryData<double>(temperature_filestream);
                double y_temperature_point = readBinaryData<double>(temperature_filestream);
                // If no data was extracted from the stream, the end of the file was reached
                if (!(temperature_filestream))
                    break;
                // Check the y value from parsed_line, to check if this point is stored on this rank
                // Check the CA grid positions of the data point to see which rank(s) should store it
                int y_int = Kokkos::round((y_temperature_point - y_min) / deltax);
                if ((y_int >= lower_y_bound) && (y_int <= upper_y_bound)) {
                    // This data point is inside the bounds of interest for this MPI rank
                    // Store the x and y values in RawData
                    raw_temperature_data(number_of_temperature_data_points, 0) = x_temperature_point;
                    raw_temperature_data(number_of_temperature_data_points, 1) = y_temperature_point;
                    // Parse the remaining 4 components (z, tm, tl, cr) from the line and store in RawData
                    for (int component = 2; component < num_temperature_components; component++) {
                        raw_temperature_data(number_of_temperature_data_points, component) =
                            readBinaryData<double>(temperature_filestream);
                    }
                    number_of_temperature_data_points++;
                    int raw_temperature_data_extent = raw_temperature_data.extent(0);
                    // Adjust size of RawData if it is near full
                    if (number_of_temperature_data_points >= raw_temperature_data_extent) {
                        Kokkos::resize(raw_temperature_data, raw_temperature_data_extent + temperature_buffer_increment,
                                       num_temperature_components);
                    }
                }
                else {
                    // This data point is inside the bounds of interest for this MPI rank
                    // ignore the z, tm, tl, cr values associated with it
                    unsigned char temp[4 * sizeof(double)];
                    temperature_filestream.read(reinterpret_cast<char *>(temp), 4 * sizeof(double));
                }
            }
        }
        else {
            // Get number of columns in this temperature file
            std::string header_line;
            getline(temperature_filestream, header_line);
            int vals_per_line = checkForHeaderValues(header_line);
            while (!temperature_filestream.eof()) {
                std::vector<std::string> parsed_line(
                    num_temperature_components); // Each line has an x, y, z, tm, tl, cr
                std::string read_line;
                if (!getline(temperature_filestream, read_line))
                    break;
                // Only parse the first 6 columns of the temperature data
                splitString(read_line, parsed_line, vals_per_line);
                // Check the y value from parsed_line, to check if this point is stored on this rank
                double y_temperature_point = getInputDouble(parsed_line[1]);
                // Check the CA grid positions of the data point to see which rank(s) should store it
                int y_int = Kokkos::round((y_temperature_point - y_min) / deltax);
                if ((y_int >= lower_y_bound) && (y_int <= upper_y_bound)) {
                    // This data point is inside the bounds of interest for this MPI rank: Store the x, z, tm, tl, and
                    // cr vals inside of RawData, incrementing with each value added
                    for (int component = 0; component < num_temperature_components; component++) {
                        raw_temperature_data(number_of_temperature_data_points, component) =
                            getInputDouble(parsed_line[component]);
                    }
                    number_of_temperature_data_points++;
                    // Adjust size of RawData if it is near full
                    int raw_temperature_data_extent = raw_temperature_data.extent(0);
                    // Adjust size of RawData if it is near full
                    if (number_of_temperature_data_points >= raw_temperature_data_extent) {
                        Kokkos::resize(raw_temperature_data, raw_temperature_data_extent + temperature_buffer_increment,
                                       num_temperature_components);
                    }
                }
            }
        }
    }

    // Read in temperature data from files, stored in the host view "RawData", with the appropriate MPI ranks storing
    // the appropriate data
    void readTemperatureData(int id, const Grid &grid, int layernumber) {

        // Y coordinates of this rank's data, inclusive and including ghost nodes
        int lower_y_bound = grid.y_offset;
        int upper_y_bound = grid.y_offset + grid.ny_local - 1;

        std::cout << "On MPI rank " << id << ", the Y bounds (in cells) are [" << lower_y_bound << "," << upper_y_bound
                  << "]" << std::endl;
        // Store raw data relevant to each rank in the vector structure RawData
        // Two passes through reading temperature data files- this is the second pass, reading the actual X/Y/Z/liquidus
        // time/cooling rate data and each rank stores the data relevant to itself in "RawData". With remelting
        // (simulation_type == "RM"), this is the same except that some X/Y/Z coordinates may be repeated in a file, and
        // a "melting time" value is stored in addition to liquidus time and cooling rate
        int number_of_temperature_data_points = 0;
        // Second pass through the files - ignore header line
        int first_layer_to_read, last_layer_to_read;
        if (_inputs.layerwise_temp_read) {
            first_layer_to_read = layernumber;
            last_layer_to_read = layernumber;
        }
        else {
            first_layer_to_read = 0;
            last_layer_to_read = std::min(grid.number_of_layers, _inputs.temp_files_in_series) - 1;
        }
        // Which temperature files should be read? Just the one file for layer "layernumber", or all of them?
        for (int layer_read_count = first_layer_to_read; layer_read_count <= last_layer_to_read; layer_read_count++) {

            std::string tempfile_thislayer;
            if (_inputs.layerwise_temp_read) {
                int LayerInSeries = layernumber % _inputs.temp_files_in_series;
                tempfile_thislayer = _inputs.temp_paths[LayerInSeries];
            }
            else
                tempfile_thislayer = _inputs.temp_paths[layer_read_count];

            first_value(layer_read_count) = number_of_temperature_data_points;
            // Read and parse temperature file for either binary or ASCII, storing the appropriate values on each MPI
            // rank within RawData and incrementing number_of_temperature_data_points appropriately
            bool binary_input_data = checkTemperatureFileFormat(tempfile_thislayer);
            parseTemperatureData(tempfile_thislayer, grid.y_min, grid.deltax, lower_y_bound, upper_y_bound,
                                 number_of_temperature_data_points, binary_input_data);
            last_value(layer_read_count) = number_of_temperature_data_points;
        } // End loop over all files read for all layers
        Kokkos::resize(raw_temperature_data, number_of_temperature_data_points, num_temperature_components);
        // Determine start values for each layer's data within "RawData", if all layers were read
        if (!(_inputs.layerwise_temp_read)) {
            if (grid.number_of_layers > _inputs.temp_files_in_series) {
                for (int layer_read_count = _inputs.temp_files_in_series; layer_read_count < grid.number_of_layers;
                     layer_read_count++) {
                    if (_inputs.temp_files_in_series == 1) {
                        // Since all layers have the same temperature data, each layer's "ZMinLayer" is just
                        // translated from that of the first layer
                        first_value(layer_read_count) = first_value(layer_read_count - 1);
                        last_value(layer_read_count) = last_value(layer_read_count - 1);
                    }
                    else {
                        // All layers have different temperature data but in a repeating pattern
                        int repeated_file = (layer_read_count) % _inputs.temp_files_in_series;
                        first_value(layer_read_count) = first_value(repeated_file);
                        last_value(layer_read_count) = last_value(repeated_file);
                    }
                }
            }
        }
    }

    // Initialize temperature data with a fixed thermal gradient in Z (can also be zero) for constrained/single grain
    // problem types
    void initialize(const int id, const std::string simulation_type, const Grid &grid, const double deltat) {

        // Initialize temperature field in Z direction with thermal gradient G set in input file
        // Liquidus front (InitUndercooling = 0) is at domain bottom for directional solidification, is at domain center
        // (with custom InitUndercooling value) for single grain solidification
        int location_init_undercooling, location_liquidus_isotherm;
        if (simulation_type == "C")
            location_init_undercooling = 0;
        else
            location_init_undercooling = Kokkos::floorf(static_cast<float>(grid.nz) / 2.0);

        // If thermal gradient is 0, liquidus isotherm does not exist - initialize to nz to avoid divide by zero error
        // and ensure all cells are initialized as undercooled (i.e., at Z coordinates less than nz)
        if (_inputs.G == 0)
            location_liquidus_isotherm = grid.nz;
        else
            location_liquidus_isotherm =
                location_init_undercooling + Kokkos::round(_inputs.init_undercooling / (_inputs.G * grid.deltax));

        // Local copies for lambda capture.
        auto layer_time_temp_history_local = layer_time_temp_history;
        auto max_solidification_events_local = max_solidification_events;
        auto number_solidification_events_local = number_of_solidification_events;
        auto undercooling_current_local = undercooling_current;
        double init_undercooling_local = _inputs.init_undercooling;
        double G_local = _inputs.G;
        double R_local = _inputs.R;
        auto policy = Kokkos::RangePolicy<execution_space>(0, grid.domain_size);
        Kokkos::parallel_for(
            "TempInitG", policy, KOKKOS_LAMBDA(const int &index) {
                // All cells past melting time step
                layer_time_temp_history_local(index, 0, 0) = -1;
                // Negative dist_from_liquidus and dist_from_init_undercooling values for cells below the liquidus
                // isotherm
                int coord_z = grid.getCoordZ(index);
                int dist_from_liquidus = coord_z - location_liquidus_isotherm;
                // Cells reach liquidus at a time dependent on their Z coordinate
                // Cells with negative liquidus time values are already undercooled, should have positive undercooling
                // and negative liquidus time step
                if (dist_from_liquidus < 0) {
                    layer_time_temp_history_local(index, 0, 1) = -1;
                    int dist_from_init_undercooling = coord_z - location_init_undercooling;
                    undercooling_current_local(index) =
                        init_undercooling_local - dist_from_init_undercooling * G_local * grid.deltax;
                }
                else {
                    // Cells with positive liquidus time values are not yet tracked - leave current undercooling as
                    // default zeros and set liquidus time step. R_local will never be zero here, as all cells at a
                    // fixed undercooling must be below the liquidus (distFromLiquidus < 0)
                    layer_time_temp_history_local(index, 0, 1) =
                        dist_from_liquidus * G_local * grid.deltax / (R_local * deltat);
                }
                // Cells cool at a constant rate
                layer_time_temp_history_local(index, 0, 2) = R_local * deltat;
                // All cells solidify once
                max_solidification_events_local(0) = 1;
                number_solidification_events_local(index) = 1;
            });
        if (id == 0) {
            std::cout << "Temperature field initialized for unidirectional solidification with G = " << G_local
                      << " K/m, initial undercooling at Z = " << location_init_undercooling << " of "
                      << _inputs.init_undercooling << " K below the liquidus" << std::endl;
            std::cout << "Done with temperature field initialization" << std::endl;
        }
    }

    // Initialize temperature data for a hemispherical spot melt (single layer simulation)
    void initialize(const int id, const Grid &grid, const double freezing_range, double deltat, double spot_radius) {

        // Each cell solidifies at most one time
        Kokkos::deep_copy(max_solidification_events, 1);

        // Outer edges of spots are initialized at the liquidus temperature
        // Spots cool at constant rate R, spot thermal gradient = G
        float isotherm_velocity = (_inputs.R / _inputs.G) * deltat / grid.deltax; // in cells per time step
        int spot_time_est = spot_radius / isotherm_velocity + (freezing_range / _inputs.R) / deltat; // in time steps
        // Spot center location - center of domain in X and Y, top of domain in Z
        float spot_center_x = spot_radius + 1;
        float spot_center_y = spot_radius + 1;
        float spot_center_z = grid.nz - 0.5;

        if (id == 0)
            std::cout << "Initializing temperature field for the hemispherical spot, which will take approximately "
                      << spot_time_est << " time steps to fully solidify" << std::endl;

        // Initialize layer_time_temp_history data values for this spot
        auto layer_time_temp_history_local = layer_time_temp_history;
        auto number_solidification_events_local = number_of_solidification_events;
        double R_local = _inputs.R;
        auto md_policy =
            Kokkos::MDRangePolicy<execution_space, Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>(
                {0, 0, 0}, {grid.nz, grid.nx, grid.ny_local});
        Kokkos::parallel_for(
            "SpotTemperatureInit", md_policy, KOKKOS_LAMBDA(const int coord_z, const int coord_x, const int coord_y) {
                // 1D cell index
                int index = grid.get1DIndex(coord_x, coord_y, coord_z);
                // Distance of this cell from the spot center
                float dist_z = spot_center_z - coord_z;
                float dist_x = spot_center_x - coord_x;
                int coord_y_global = coord_y + grid.y_offset;
                float dist_y = spot_center_y - coord_y_global;
                float tot_dist = Kokkos::hypot(dist_x, dist_y, dist_z);
                if (tot_dist <= spot_radius) {
                    // Melt time step (first time step)
                    layer_time_temp_history_local(index, 0, 0) = 1;
                    // Liquidus time step (related to distance to spot edge)
                    layer_time_temp_history_local(index, 0, 1) =
                        Kokkos::round((spot_radius - tot_dist) / isotherm_velocity) + 1;
                    // Cooling rate per time step
                    layer_time_temp_history_local(index, 0, 2) = R_local * deltat;
                    number_solidification_events_local(index) = 1;
                }
                else {
                    // Dummy layer_time_temp_history data
                    for (int l = 0; l < 3; l++) {
                        layer_time_temp_history_local(index, 0, l) = -1.0;
                    }
                    number_solidification_events_local(index) = 0;
                }
            });
        MPI_Barrier(MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Spot melt temperature field initialized" << std::endl;
    }

    // Calculate the number of times that a cell in layer "layernumber" undergoes melting/solidification, and store in
    // max_solidification_events_host
    void calcMaxSolidificationEvents(const int id, const int layernumber,
                                     const view_type_int_host max_solidification_events_host, const int start_range,
                                     const int end_range, const Grid &grid) {

        if (layernumber > _inputs.temp_files_in_series) {
            // Use the value from a previously checked layer, since the time-temperature history is reused
            if (_inputs.temp_files_in_series == 1) {
                // All layers have the same temperature data, max_solidification_events for this layer is the same as
                // the last
                max_solidification_events_host(layernumber) = max_solidification_events_host(layernumber - 1);
            }
            else {
                // All layers have different temperature data but in a repeating pattern
                int repeated_file = layernumber % _inputs.temp_files_in_series;
                max_solidification_events_host(layernumber) = max_solidification_events_host(repeated_file);
            }
        }
        else {
            // Need to calculate max_solidification_events(layernumber) from the values in RawData
            // Init to 0
            view_type_int_host temp_melt_count("temp_melt_count", grid.domain_size);

            for (int i = start_range; i < end_range; i++) {

                // Get the integer X, Y, Z coordinates associated with this data point, with the Y coordinate based on
                // local MPI rank's grid
                int coord_x = getTempCoordX(i, grid.x_min, grid.deltax);
                int coord_y = getTempCoordY(i, grid.y_min, grid.deltax, grid.y_offset);
                int coord_z = getTempCoordZ(i, grid.deltax, grid.layer_height, layernumber, grid.z_min_layer);
                // Convert to 1D coordinate in the current layer's domain
                int index = grid.get1DIndex(coord_x, coord_y, coord_z);
                temp_melt_count(index)++;
            }
            int max_count = 0;
            for (int i = 0; i < grid.domain_size; i++) {
                if (temp_melt_count(i) > max_count)
                    max_count = temp_melt_count(i);
            }
            int max_count_global;
            MPI_Allreduce(&max_count, &max_count_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            max_solidification_events_host(layernumber) = max_count_global;
        }
        if (id == 0)
            std::cout << "The maximum number of melting/solidification events during layer " << layernumber << " is "
                      << max_solidification_events_host(layernumber) << std::endl;
    }

    // Read data from storage, and calculate the normalized x value of the data point
    int getTempCoordX(const int i, const double x_min, const double deltax) {
        int x_coord = Kokkos::round((raw_temperature_data(i, 0) - x_min) / deltax);
        return x_coord;
    }
    // Read data from storage, and calculate the normalized y value of the data point. If the optional offset argument
    // is given, the return value is calculated relative to the edge of the MPI rank's local simulation domain (which is
    // offset by y_offset cells from the global domain edge)
    int getTempCoordY(const int i, const double y_min, const double deltax, const int y_offset = 0) {
        int y_coord = Kokkos::round((raw_temperature_data(i, 1) - y_min) / deltax) - y_offset;
        return y_coord;
    }
    // Read data from storage, and calculate the normalized z value of the data point
    int getTempCoordZ(const int i, const double deltax, const int layer_height, const int layer_counter,
                      const view_type_double_host z_min_layer) {
        int z_coord = Kokkos::round(
            (raw_temperature_data(i, 2) + deltax * layer_height * layer_counter - z_min_layer[layer_counter]) / deltax);
        return z_coord;
    }
    // Read data from storage, obtain melting time
    double getTempCoordTM(const int i) {
        double TMelting = raw_temperature_data(i, 3);
        return TMelting;
    }
    // Read data from storage, obtain liquidus time
    double getTempCoordTL(const int i) {
        double TLiquidus = raw_temperature_data(i, 4);
        return TLiquidus;
    }
    // Read data from storage, obtain cooling rate
    double getTempCoordCR(const int i) {
        double cooling_rate = raw_temperature_data(i, 5);
        return cooling_rate;
    }

    // Initialize temperature fields for layer "layernumber" in case where temperature data comes from file(s)
    // TODO: This can be performed on the device as the dirS problem is
    void initialize(const int layernumber, const int id, const Grid &grid, const double freezing_range,
                    const double deltat) {

        // Data was already read into the "raw_temperature_data" data structure
        // Determine which section of "raw_temperature_data" is relevant for this layer of the overall domain
        int start_range = first_value[layernumber];
        int end_range = last_value[layernumber];

        // Temporary host view for the maximum number of times a cell in a given layer will solidify (different for each
        // layer)
        view_type_int_host max_solidification_events_host =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), max_solidification_events);
        calcMaxSolidificationEvents(id, layernumber, max_solidification_events_host, start_range, end_range, grid);
        int max_num_solidification_events = max_solidification_events_host(layernumber);

        // Resize layer_time_temp_history now that the max number of solidification events is known for this layer
        Kokkos::resize(layer_time_temp_history, grid.domain_size, max_num_solidification_events, 3);

        // These views are initialized to zeros on the host, filled with data, and then copied to the device for layer
        // "layernumber"
        view_type_float_3d_host layer_time_temp_history_host("TimeTempHistory_H", grid.domain_size,
                                                             max_num_solidification_events, 3);
        view_type_int_host number_of_solidification_events_host("NumSEvents_H", grid.domain_size);

        double largest_time = 0;
        double largest_time_global = 0;
        if (id == 0)
            std::cout << "Range of raw data for layer " << layernumber << " on rank 0 is " << start_range << " to "
                      << end_range << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        for (int i = start_range; i < end_range; i++) {

            // Get the integer X, Y, Z coordinates associated with this data point, along with the associated TM, TL, CR
            // values
            // coord_y is relative to this MPI rank's grid, while coord_y_global is relative to the overall simulation
            // domain
            int coord_x = getTempCoordX(i, grid.x_min, grid.deltax);
            int coord_y = getTempCoordY(i, grid.y_min, grid.deltax, grid.y_offset);
            int coord_z = getTempCoordZ(i, grid.deltax, grid.layer_height, layernumber, grid.z_min_layer);
            double t_melting = getTempCoordTM(i);
            double t_liquidus = getTempCoordTL(i);
            double cooling_rate = getTempCoordCR(i);

            // 1D cell coordinate on this MPI rank's domain
            int index = grid.get1DIndex(coord_x, coord_y, coord_z);
            // Store TM, TL, CR values for this solidification event in layer_time_temp_history
            layer_time_temp_history_host(index, number_of_solidification_events_host(index), 0) =
                Kokkos::round(t_melting / deltat) + 1;
            layer_time_temp_history_host(index, number_of_solidification_events_host(index), 1) =
                Kokkos::round(t_liquidus / deltat) + 1;
            layer_time_temp_history_host(index, number_of_solidification_events_host(index), 2) =
                std::abs(cooling_rate) * deltat;
            // Increment number of solidification events for this cell
            number_of_solidification_events_host(index)++;
            // Estimate of the time step where the last possible solidification is expected to occur
            double solidus_time = t_liquidus + freezing_range / cooling_rate;
            if (solidus_time > largest_time)
                largest_time = solidus_time;
        }
        MPI_Allreduce(&largest_time, &largest_time_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (id == 0)
            std::cout << " Layer " << layernumber
                      << " time step where all cells are cooled below solidus for the final time is "
                      << Kokkos::round((largest_time_global) / deltat) << " (or " << largest_time_global << " seconds)"
                      << std::endl;
        if (id == 0)
            std::cout << "Layer " << layernumber << " temperatures read" << std::endl;

        // Reorder solidification events in layer_time_temp_history(location,event number,component) so that they are in
        // order based on the melting time values (component = 0)
        for (int index = 0; index < grid.domain_size; index++) {
            int n_solidification_events_cell = number_of_solidification_events_host(index);
            if (n_solidification_events_cell > 0) {
                for (int i = 0; i < n_solidification_events_cell - 1; i++) {
                    for (int j = (i + 1); j < n_solidification_events_cell; j++) {
                        if (layer_time_temp_history_host(index, i, 0) > layer_time_temp_history_host(index, j, 0)) {
                            // Swap these two points - melting event "j" happens before event "i"
                            float old_melt_val = layer_time_temp_history_host(index, i, 0);
                            float old_liq_val = layer_time_temp_history_host(index, i, 1);
                            float old_cr_val = layer_time_temp_history_host(index, i, 2);
                            layer_time_temp_history_host(index, i, 0) = layer_time_temp_history_host(index, j, 0);
                            layer_time_temp_history_host(index, i, 1) = layer_time_temp_history_host(index, j, 1);
                            layer_time_temp_history_host(index, i, 2) = layer_time_temp_history_host(index, j, 2);
                            layer_time_temp_history_host(index, j, 0) = old_melt_val;
                            layer_time_temp_history_host(index, j, 1) = old_liq_val;
                            layer_time_temp_history_host(index, j, 2) = old_cr_val;
                        }
                    }
                }
            }
        }
        // If a cell melts twice before reaching the liquidus temperature, this is a double counted solidification
        // event and should be removed
        for (int index = 0; index < grid.domain_size; index++) {
            int n_solidification_events_cell = number_of_solidification_events_host(index);
            if (n_solidification_events_cell > 1) {
                for (int i = 0; i < n_solidification_events_cell - 1; i++) {
                    if (layer_time_temp_history_host(index, i + 1, 0) < layer_time_temp_history_host(index, i, 1)) {
                        std::cout << "Cell " << index << " removing anomalous event " << i + 1 << " out of "
                                  << n_solidification_events_cell - 1 << std::endl;
                        // Keep whichever event has the larger liquidus time
                        if (layer_time_temp_history_host(index, i + 1, 1) > layer_time_temp_history_host(index, i, 1)) {
                            layer_time_temp_history_host(index, i, 0) = layer_time_temp_history_host(index, i + 1, 0);
                            layer_time_temp_history_host(index, i, 1) = layer_time_temp_history_host(index, i + 1, 1);
                            layer_time_temp_history_host(index, i, 2) = layer_time_temp_history_host(index, i + 1, 2);
                        }
                        layer_time_temp_history_host(index, i + 1, 0) = 0.0;
                        layer_time_temp_history_host(index, i + 1, 1) = 0.0;
                        layer_time_temp_history_host(index, i + 1, 2) = 0.0;
                        // Reshuffle other solidification events over if needed
                        for (int ii = (i + 1); ii < n_solidification_events_cell - 1; ii++) {
                            layer_time_temp_history_host(index, ii, 0) = layer_time_temp_history_host(index, ii + 1, 0);
                            layer_time_temp_history_host(index, ii, 1) = layer_time_temp_history_host(index, ii + 1, 1);
                            layer_time_temp_history_host(index, ii, 2) = layer_time_temp_history_host(index, ii + 1, 2);
                        }
                        number_of_solidification_events_host(index)--;
                    }
                }
            }
        }

        // Copy host view data back to device
        max_solidification_events = Kokkos::create_mirror_view_and_copy(memory_space(), max_solidification_events_host);
        layer_time_temp_history = Kokkos::create_mirror_view_and_copy(memory_space(), layer_time_temp_history_host);
        number_of_solidification_events =
            Kokkos::create_mirror_view_and_copy(memory_space(), number_of_solidification_events_host);

        if (id == 0) {
            std::cout << "Layer " << layernumber << " temperature field is from Z = " << grid.z_layer_bottom
                      << " through " << grid.nz_layer + grid.z_layer_bottom - 1 << " of the global domain" << std::endl;
            std::cout << "Done with temperature field initialization" << std::endl;
        }
    }

    // Get the subview associated with the undercooling of cells in the current layer. Do not reset the undercooling of
    // cells from the prior layer to zero as this information will be stored for a potential print (and a cell that
    // remelts in the current layer will have its undercooling reset to 0 and recalculated)
    void getCurrentLayerUndercooling(std::pair<int, int> layer_range) {
        undercooling_current = Kokkos::subview(undercooling_current_all_layers, layer_range);
    }

    // (Optional based on selected inputs) Get the subview associated with the initial undercooling of cells during
    // solidification start in the current layer. Do not reset the undercooling of cells from the prior layer to zero as
    // this information will be stored for a potential print (and a cell that remelts in the current layer will have its
    // undercooling reset to 0 and recalculated)
    void getCurrentLayerStartingUndercooling(std::pair<int, int> layer_range) {
        undercooling_solidification_start = Kokkos::subview(undercooling_solidification_start_all_layers, layer_range);
    }

    // For each Z coordinate, find the smallest undercooling at which solidification started and finished, writing this
    // data to an output file
    view_type_float_2d_host getFrontUndercoolingStartFinish(const int id, const Grid &grid) {
        view_type_float start_solidification_z(Kokkos::ViewAllocateWithoutInitializing("start_solidification_z"),
                                               grid.nz);
        view_type_float end_solidification_z(Kokkos::ViewAllocateWithoutInitializing("end_solidification_z"), grid.nz);
        auto _undercooling_solidification_start = undercooling_solidification_start;
        auto _undercooling_current = undercooling_current;
        auto policy = Kokkos::RangePolicy<execution_space>(0, grid.nz);
        Kokkos::parallel_for(
            "GetMinUndercooling", policy, KOKKOS_LAMBDA(const int &coord_z) {
                float min_start_undercooling = Kokkos::Experimental::finite_max_v<float>;
                float min_end_undercooling = Kokkos::Experimental::finite_max_v<float>;
                for (int coord_x = 0; coord_x < grid.nx; coord_x++) {
                    for (int coord_y = 1; coord_y < grid.ny_local - 1; coord_y++) {
                        int index = grid.get1DIndex(coord_x, coord_y, coord_z);
                        if (_undercooling_solidification_start(index) < min_start_undercooling)
                            min_start_undercooling = _undercooling_solidification_start(index);
                        if (_undercooling_current(index) < min_end_undercooling)
                            min_end_undercooling = _undercooling_current(index);
                    }
                }
                start_solidification_z(coord_z) = min_start_undercooling;
                end_solidification_z(coord_z) = min_end_undercooling;
            });

        // Rank 0 - collect min values from all ranks to get the global start/end solidification undercoolings
        view_type_float start_solidification_z_reduced(
            Kokkos::ViewAllocateWithoutInitializing("start_solidification_z_red"), grid.nz);
        view_type_float end_solidification_z_reduced(
            Kokkos::ViewAllocateWithoutInitializing("end_solidification_z_red"), grid.nz);
        // Could potentially perform these reductions on the 2D view storing both start and end values, but depends on
        // 2D view data layout
        MPI_Reduce(start_solidification_z.data(), start_solidification_z_reduced.data(), grid.nz, MPI_FLOAT, MPI_MIN, 0,
                   MPI_COMM_WORLD);
        MPI_Reduce(end_solidification_z.data(), end_solidification_z_reduced.data(), grid.nz, MPI_FLOAT, MPI_MIN, 0,
                   MPI_COMM_WORLD);

        // Rank 0 - copy to host view
        view_type_float_2d_host start_end_solidification_z_host(
            Kokkos::ViewAllocateWithoutInitializing("start_end_solidification_z_host"), grid.nz, 2);
        if (id == 0) {
            view_type_float_host start_solidification_z_host =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), start_solidification_z_reduced);
            view_type_float_host end_solidification_z_host =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), end_solidification_z_reduced);
            for (int coord_z = 0; coord_z < grid.nz; coord_z++) {
                start_end_solidification_z_host(coord_z, 0) = start_solidification_z_host(coord_z);
                start_end_solidification_z_host(coord_z, 1) = end_solidification_z_host(coord_z);
            }
        }
        return start_end_solidification_z_host;
    }

    // Reset local cell undercooling (and if needed, the cell's starting undercooling) to 0
    KOKKOS_INLINE_FUNCTION
    void resetUndercooling(const int index) const {
        if (_store_solidification_start)
            undercooling_solidification_start(index) = 0.0;
        undercooling_current(index) = 0.0;
    }

    // Update local cell undercooling for the current melt-resolidification event
    KOKKOS_INLINE_FUNCTION
    void updateUndercooling(const int index) const {
        undercooling_current(index) += layer_time_temp_history(index, solidification_event_counter(index), 2);
    }

    // (Optional based on inputs) Set the starting undercooling in the cell for the solidification event that just
    // started
    KOKKOS_INLINE_FUNCTION
    void setStartingUndercooling(const int index) const {
        if (_store_solidification_start)
            undercooling_solidification_start(index) = undercooling_current(index);
    }

    // Update the solidification event counter for a cell that has not finished the previous solidification event (i.e.,
    // solidification is not complete in this cell as it is not tempsolid or solid type)
    KOKKOS_INLINE_FUNCTION
    void updateSolidificationCounter(const int index) const { solidification_event_counter(index)++; }

    // Update the solidification event counter for the cell (which is either tempsolid or solid type) and return whether
    // all solidification events have completed in the cell
    KOKKOS_INLINE_FUNCTION
    bool updateCheckSolidificationCounter(const int index) const {
        bool solidification_complete;
        solidification_event_counter(index)++;
        if (solidification_event_counter(index) == number_of_solidification_events(index))
            solidification_complete = true;
        else
            solidification_complete = false;
        return solidification_complete;
    }

    // Reset solidification event counter and get the subview associated with the undercooling field for the next layer
    void resetLayerEventsUndercooling(const Grid &grid) {
        getCurrentLayerUndercooling(grid.layer_range);
        if (_store_solidification_start)
            getCurrentLayerStartingUndercooling(grid.layer_range);
        Kokkos::realloc(solidification_event_counter, grid.domain_size);
        Kokkos::deep_copy(solidification_event_counter, 0);
    }

    // Extract the next time that this point undergoes melting
    KOKKOS_INLINE_FUNCTION
    int getMeltTimeStep(const int cycle, const int index) const {
        int melt_time_step;
        int solidification_event_counter_cell = solidification_event_counter(index);
        melt_time_step = static_cast<int>(layer_time_temp_history(index, solidification_event_counter_cell, 0));
        if (cycle > melt_time_step) {
            // If the cell has already exceeded the melt time step for the current melt-solidification event, get the
            // melt time step associated with the next solidification event - or, if there is no next
            // melt-solidification event, return the max possible int as the cell will not melt again during this layer
            // of the multilayer problem
            if (solidification_event_counter_cell < (number_of_solidification_events(index) - 1))
                melt_time_step =
                    static_cast<int>(layer_time_temp_history(index, solidification_event_counter_cell + 1, 0));
            else
                melt_time_step = INT_MAX;
        }
        return melt_time_step;
    }

    // Extract the next time that this point cools below the liquidus
    // Uses the current value of the solidification event counter
    KOKKOS_INLINE_FUNCTION
    int getCritTimeStep(const int index) const {
        int solidification_event_counter_cell = solidification_event_counter(index);
        int crit_time_step = static_cast<int>(layer_time_temp_history(index, solidification_event_counter_cell, 1));
        return crit_time_step;
    }
    // Uses a specified solidification event
    KOKKOS_INLINE_FUNCTION
    int getCritTimeStep(const int index, const int solidification_event_counter_cell) const {
        int crit_time_step = static_cast<int>(layer_time_temp_history(index, solidification_event_counter_cell, 1));
        return crit_time_step;
    }

    // Extract the cooling rate associated with a specified solidificaiton event
    KOKKOS_INLINE_FUNCTION
    float getUndercoolingChange(const int index, const int solidification_event_counter_cell) const {
        float undercooling_change = layer_time_temp_history(index, solidification_event_counter_cell, 2);
        return undercooling_change;
    }

    // Extract either the last time step that all points undergo melting in the layer, the last time they cools below
    // the liquidus, or the rate at which they cools from the liquidus from layer_time_temp_history (corresponds to
    // solidification event number `NumSolidificationEvents-1` for the cell) (can't just use subview here since
    // NumSolidificationEvents is different for each cell) If the cell does not undergo solidification, either print -1
    // or the specified default value
    template <typename extracted_view_data_type>
    extracted_view_data_type extractTmTlCrData(const int extracted_val, const int domain_size,
                                               const int default_val = -1) {
        extracted_view_data_type extracted_data(Kokkos::ViewAllocateWithoutInitializing("extracted_data"), domain_size);
        using extracted_value_type = typename extracted_view_data_type::value_type;

        // Local copies for lambda capture.
        auto layer_time_temp_history_local = layer_time_temp_history;
        auto number_of_solidification_events_local = number_of_solidification_events;

        auto policy = Kokkos::RangePolicy<execution_space>(0, domain_size);
        Kokkos::parallel_for(
            "Extract_tm_tl_cr_data", policy, KOKKOS_LAMBDA(const int &index) {
                int num_solidification_events_this_cell = number_of_solidification_events_local(index);
                // If this cell doesn't undergo solidification at all, print -1
                if (num_solidification_events_this_cell == 0)
                    extracted_data(index) = static_cast<extracted_value_type>(default_val);
                else
                    extracted_data(index) = static_cast<extracted_value_type>(
                        layer_time_temp_history_local(index, num_solidification_events_this_cell - 1, extracted_val));
            });
        return extracted_data;
    }
};

#endif
