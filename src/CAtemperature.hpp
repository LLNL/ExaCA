// Copyright Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_TEMPS_HPP
#define EXACA_TEMPS_HPP

#include "CAgrid.hpp"
#include "CAinputdata.hpp"
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
    using view_type_int_3d = Kokkos::View<int ***, memory_space>;
    using view_type_double = Kokkos::View<double *, memory_space>;
    using view_type_double_2d = Kokkos::View<double **, memory_space>;
    using view_type_float = Kokkos::View<float *, memory_space>;
    using view_type_float_2d = Kokkos::View<float **, memory_space>;
    using view_type_int_host = typename view_type_int::HostMirror;
    using view_type_int_3d_host = typename view_type_int_3d::HostMirror;
    using view_type_double_host = typename view_type_double::HostMirror;
    using view_type_double_2d_host = typename view_type_double_2d::HostMirror;
    using view_type_float_host = typename view_type_float::HostMirror;
    using view_type_float_2d_host = typename view_type_float_2d::HostMirror;
    using view_type_coupled = Kokkos::View<double **, Kokkos::LayoutLeft, Kokkos::HostSpace>;

    // Using the default exec space for this memory space.
    using execution_space = typename memory_space::execution_space;

    // Maximum number of times each CA cell in a given layer undergoes solidification
    view_type_int max_solidification_events;
    // Each time that a cell (index 0) goes above and below the liquidus time (index 2) for each melt-solidification
    // event (index 1)
    view_type_int_3d liquidus_time;
    // Cooling rate for each cell (index 1) for each solidification event (index 2)
    view_type_float_2d cooling_rate;
    // The number of times that each CA cell will traverse the liquidus during this layer
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

    // Constructor creates views with size based on the grid inputs - liquidus_time and cooling_rate are later modified
    // to account for multiple events if needed, undercooling_current and solidification_event_counter are default
    // initialized to zeros
    Temperature(const Grid &grid, TemperatureInputs inputs, const bool store_solidification_start = false,
                const int est_num_temperature_data_points = 100000)
        : max_solidification_events(
              view_type_int(Kokkos::ViewAllocateWithoutInitializing("number_of_layers"), grid.number_of_layers))
        , liquidus_time(
              view_type_int_3d(Kokkos::ViewAllocateWithoutInitializing("liquidus_time"), grid.domain_size, 1, 2))
        , cooling_rate(view_type_float_2d(Kokkos::ViewAllocateWithoutInitializing("cooling_rate"), grid.domain_size, 1))
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

    // Constructor using in-memory temperature data from external source (input_temperature_data)
    Temperature(const int id, const int np, const Grid &grid, TemperatureInputs inputs,
                view_type_coupled input_temperature_data, const bool store_solidification_start = false)
        : max_solidification_events(
              view_type_int(Kokkos::ViewAllocateWithoutInitializing("number_of_layers"), grid.number_of_layers))
        , liquidus_time(
              view_type_int_3d(Kokkos::ViewAllocateWithoutInitializing("liquidus_time"), grid.domain_size, 1, 2))
        , cooling_rate(view_type_float_2d(Kokkos::ViewAllocateWithoutInitializing("cooling_rate"), grid.domain_size, 1))
        , number_of_solidification_events(view_type_int(
              Kokkos::ViewAllocateWithoutInitializing("number_of_solidification_events"), grid.domain_size))
        , solidification_event_counter(view_type_int("solidification_event_counter", grid.domain_size))
        , undercooling_current_all_layers(view_type_float("undercooling_current", grid.domain_size_all_layers))
        , raw_temperature_data(view_type_double_2d_host(Kokkos::ViewAllocateWithoutInitializing("raw_temperature_data"),
                                                        input_temperature_data.extent(0) * (inputs.number_of_copies),
                                                        num_temperature_components))
        , first_value(view_type_int_host(Kokkos::ViewAllocateWithoutInitializing("first_value"), grid.number_of_layers))
        , last_value(view_type_int_host(Kokkos::ViewAllocateWithoutInitializing("last_value"), grid.number_of_layers))
        , _store_solidification_start(store_solidification_start)
        , _inputs(inputs) {

        copyTemperatureData(id, np, grid, input_temperature_data);
        if (_store_solidification_start) {
            // Default init starting undercooling in cells to zero
            undercooling_solidification_start_all_layers =
                view_type_float("undercooling_solidification_start", grid.domain_size_all_layers);
            getCurrentLayerStartingUndercooling(grid.layer_range);
        }
        getCurrentLayerUndercooling(grid.layer_range);
    }

    // If using a problem type that involves translating or mirroring temperature data, output a new X or Y value based
    // on value from the file/from finch
    double getTranslatedXY(const double input_xy_value, const double xy_min, const double xy_max,
                           const double xy_offset, const int line_number, const bool mirror) {
        double xy_translated = input_xy_value + line_number * xy_offset;
        if ((mirror) && (line_number % 2 == 1))
            xy_translated = xy_max - xy_translated + xy_min;
        return xy_translated;
    }

    // Based on the X coordinate compared to the global domain bounds and the Y coordinate compared to the MPI rank Y
    // bounds (integers), return whether or not this translated coordinate is in bounds on this rank
    bool translatedXYInBounds(const double x_translated, const double y_translated, const Grid grid) {
        const int coord_x = Kokkos::round((x_translated - grid.x_min) / grid.deltax);
        const int coord_y_global = Kokkos::round((y_translated - grid.y_min) / grid.deltax);
        bool xy_in_bounds;
        if ((coord_x >= 0) && (coord_x < grid.nx) && (coord_y_global >= grid.y_offset) &&
            (coord_y_global < grid.y_offset + grid.ny_local))
            xy_in_bounds = true;
        else
            xy_in_bounds = false;
        return xy_in_bounds;
    }

    // Copy data from external source onto the appropriate MPI ranks
    void copyTemperatureData(const int id, const int np, const Grid &grid, view_type_coupled input_temperature_data,
                             const int resize_padding = 100) {

        // Take first num_temperature_components columns of input_temperature_data
        int finch_data_size = input_temperature_data.extent(0);
        int finch_temp_components = input_temperature_data.extent(1);
        std::cout << "Rank " << id << " has " << finch_data_size << " events from the Finch simulation" << std::endl;

        // First, store data with Y coordinates in bounds for this rank in raw_temperature_data
        int temperature_point_counter = 0;
        extractTemperatureData(input_temperature_data, raw_temperature_data, temperature_point_counter, grid);

        // Communication pattern - sending to right, receiving from left
        int left, right;
        if (id == 0)
            left = np - 1;
        else
            left = id - 1;
        if (id == np - 1)
            right = 0;
        else
            right = id + 1;

        // Send and recieve data so each rank parses all finch data points
        int send_data_size = finch_data_size;
        for (int i = 0; i < np - 1; i++) {
            // Get size for sending/receiving
            int recv_data_size;
            MPI_Request send_request_size, recv_request_size;
            MPI_Isend(&send_data_size, 1, MPI_INT, right, 0, MPI_COMM_WORLD, &send_request_size);
            MPI_Irecv(&recv_data_size, 1, MPI_INT, left, 0, MPI_COMM_WORLD, &recv_request_size);
            MPI_Wait(&send_request_size, MPI_STATUS_IGNORE);
            MPI_Wait(&recv_request_size, MPI_STATUS_IGNORE);
            // Allocate view for received data
            view_type_coupled finch_data_recv(Kokkos::ViewAllocateWithoutInitializing("finch_data_recv"),
                                              recv_data_size, finch_temp_components);

            // Send data to the right, recieve data from the left - if needed, increase size of data stored on this rank
            // to accomodate received data
            MPI_Request send_request_data, recv_request_data;
            MPI_Isend(input_temperature_data.data(), send_data_size * finch_temp_components, MPI_DOUBLE, right, 1,
                      MPI_COMM_WORLD, &send_request_data);
            MPI_Irecv(finch_data_recv.data(), recv_data_size * finch_temp_components, MPI_DOUBLE, left, 1,
                      MPI_COMM_WORLD, &recv_request_data);
            int current_size = raw_temperature_data.extent(0);
            if (temperature_point_counter + _inputs.number_of_copies * recv_data_size >= current_size)
                Kokkos::resize(raw_temperature_data,
                               temperature_point_counter + _inputs.number_of_copies * recv_data_size + resize_padding,
                               num_temperature_components);
            MPI_Wait(&send_request_data, MPI_STATUS_IGNORE);
            MPI_Wait(&recv_request_data, MPI_STATUS_IGNORE);

            // Unpack the appropriate received data into raw_temperature_data
            extractTemperatureData(finch_data_recv, raw_temperature_data, temperature_point_counter, grid);

            // Replace send buffer with the received data
            input_temperature_data = finch_data_recv;
            send_data_size = recv_data_size;
        }
        // Resize with the number of temperature data points on this rank now known
        Kokkos::resize(raw_temperature_data, temperature_point_counter, num_temperature_components);
        for (int n = 0; n < grid.number_of_layers; n++) {
            first_value(n) = 0;
            last_value(n) = temperature_point_counter;
        }
        std::cout << "Rank " << id << " has " << temperature_point_counter << " events to simulate with ExaCA"
                  << std::endl;
    }

    // Get temperature data from Finch and store as raw temperatures.
    template <typename SrcViewType, typename DstViewType>
    void extractTemperatureData(const SrcViewType &temp_src, DstViewType &temp_dst, int &temp_count, const Grid grid) {
        int data_size = temp_src.extent(0);
        for (int n = 0; n < data_size; n++) {
            for (int l = 0; l < _inputs.number_of_copies; l++) {
                // Consider each translated point in either X or Y, optionally mirroring the point in the non-translated
                // direction
                // If using a problem type that involves translating or mirroring temperature data, output a new X value
                // based on value from the file/from finch
                const double x_translated =
                    getTranslatedXY(temp_src(n, 0), grid.x_min, grid.x_max, _inputs.x_offset, l, _inputs.mirror_x);
                const double y_translated =
                    getTranslatedXY(temp_src(n, 1), grid.y_min, grid.y_max, _inputs.y_offset, l, _inputs.mirror_y);
                const bool xy_in_bounds = translatedXYInBounds(x_translated, y_translated, grid);
                // Only store if the point is in-bounds in X and Y
                if (xy_in_bounds) {
                    temp_dst(temp_count, 0) = x_translated;
                    temp_dst(temp_count, 1) = y_translated;
                    temp_dst(temp_count, 2) = temp_src(n, 2);
                    // Offset for melting and liquidus times
                    double time_offset = l * _inputs.temporal_offset;
                    temp_dst(temp_count, 3) = temp_src(n, 3) + time_offset;
                    temp_dst(temp_count, 4) = temp_src(n, 4) + time_offset;
                    temp_dst(temp_count, 5) = temp_src(n, 5);
                    // Increment counter for each point stored on this rank
                    temp_count++;
                }
            }
        }
    }

    // Read and parse the temperature file (double precision values in a comma-separated, ASCII format with a header
    // line - or a binary string of double precision values), storing the x, y, z, tm, tl, cr values in the
    // raw_temperature_data view. Each rank only contains the points corresponding to cells within the associated Y
    // bounds. number_of_temperature_data_points is incremented on each rank as data is added to raw_temperature_data
    void parseTemperatureData(const std::string tempfile_thislayer, const Grid grid,
                              int &number_of_temperature_data_points, const bool binary_input_data,
                              const int temperature_buffer_increment = 100000) {

        std::ifstream temperature_filestream;
        temperature_filestream.open(tempfile_thislayer);
        if (binary_input_data) {
            while (!temperature_filestream.eof()) {
                std::vector<double> parsed_data(num_temperature_components);
                // Get x, y, z, tm, tl, cr - additional temperature components (Gx, Gy, Gz) have so far have been unused
                // in ExaCA - read and ignore these
                for (int component = 0; component < num_temperature_components; component++)
                    parsed_data[component] = readBinaryData<double>(temperature_filestream);
                // If no data was extracted from the stream, the end of the file was reached
                if (!(temperature_filestream))
                    break;
                // If translating/rotating line scan data, do this in a loop
                for (int l = 0; l < _inputs.number_of_copies; l++) {
                    // Consider each translated point in either X or Y, optionally mirroring the point in the
                    // non-translated direction If using a problem type that involves translating or mirroring
                    // temperature data, output a new X value based on value from the file/from finch
                    const double x_translated =
                        getTranslatedXY(parsed_data[0], grid.x_min, grid.x_max, _inputs.x_offset, l, _inputs.mirror_x);
                    const double y_translated =
                        getTranslatedXY(parsed_data[1], grid.y_min, grid.y_max, _inputs.y_offset, l, _inputs.mirror_y);
                    const bool xy_in_bounds = translatedXYInBounds(x_translated, y_translated, grid);
                    if (xy_in_bounds) {
                        // This data point is inside the bounds of interest for this MPI rank
                        // Store the x and y values in RawData
                        raw_temperature_data(number_of_temperature_data_points, 0) = x_translated;
                        raw_temperature_data(number_of_temperature_data_points, 1) = y_translated;
                        // z coordinate is stored as-is
                        raw_temperature_data(number_of_temperature_data_points, 2) = parsed_data[2];
                        // tm and tl may require an offset in time if l > 0
                        const double time_offset = l * _inputs.temporal_offset;
                        raw_temperature_data(number_of_temperature_data_points, 3) = parsed_data[3] + time_offset;
                        raw_temperature_data(number_of_temperature_data_points, 4) = parsed_data[4] + time_offset;
                        // cr is stored as-is
                        raw_temperature_data(number_of_temperature_data_points, 5) = parsed_data[5];
                        // Increment number of temperature data points stored
                        number_of_temperature_data_points++;
                        int raw_temperature_data_extent = raw_temperature_data.extent(0);
                        // Adjust size of raw_temperature_data if it is near full
                        if (number_of_temperature_data_points >= raw_temperature_data_extent)
                            Kokkos::resize(raw_temperature_data,
                                           raw_temperature_data_extent + temperature_buffer_increment,
                                           num_temperature_components);
                    }
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
                // If translating/rotating line scan data, do this in a loop
                for (int l = 0; l < _inputs.number_of_copies; l++) {
                    // Consider each translated point in either X or Y, optionally mirroring the point in the
                    // non-translated direction If using a problem type that involves translating or mirroring
                    // temperature data, output a new X value based on value from the file/from finch
                    const double x_temperature_point = getInputDouble(parsed_line[0]);
                    const double y_temperature_point = getInputDouble(parsed_line[1]);
                    const double x_translated = getTranslatedXY(x_temperature_point, grid.x_min, grid.x_max,
                                                                _inputs.x_offset, l, _inputs.mirror_x);
                    const double y_translated = getTranslatedXY(y_temperature_point, grid.y_min, grid.y_max,
                                                                _inputs.y_offset, l, _inputs.mirror_y);
                    const bool xy_in_bounds = translatedXYInBounds(x_translated, y_translated, grid);
                    // Only store if the point is in-bounds in X and Y
                    if (xy_in_bounds) {
                        // Store the translated/mirrored x, y, and offset tm and tl values in raw_temperature_data along
                        // with z and cr value (unchanged during any translation/mirroring of data)
                        raw_temperature_data(number_of_temperature_data_points, 0) = x_translated;
                        raw_temperature_data(number_of_temperature_data_points, 1) = y_translated;
                        raw_temperature_data(number_of_temperature_data_points, 2) = getInputDouble(parsed_line[2]);
                        const double time_offset = l * _inputs.temporal_offset;
                        raw_temperature_data(number_of_temperature_data_points, 3) =
                            getInputDouble(parsed_line[3]) + time_offset;
                        raw_temperature_data(number_of_temperature_data_points, 4) =
                            getInputDouble(parsed_line[4]) + time_offset;
                        raw_temperature_data(number_of_temperature_data_points, 5) = getInputDouble(parsed_line[5]);
                        number_of_temperature_data_points++;
                        // Adjust size of raw_temperature_data if it is near full
                        int raw_temperature_data_extent = raw_temperature_data.extent(0);
                        // Adjust size of raw_temperature_data if it is near full
                        if (number_of_temperature_data_points >= raw_temperature_data_extent) {
                            Kokkos::resize(raw_temperature_data,
                                           raw_temperature_data_extent + temperature_buffer_increment,
                                           num_temperature_components);
                        }
                    }
                }
            }
        }
    }

    // Read in temperature data from files, stored in the host view raw_temperature_data, with the appropriate MPI ranks
    // storing the appropriate data
    void readTemperatureData(int id, const Grid &grid, int layernumber) {

        // Y coordinates of this rank's data, inclusive and including ghost nodes
        int lower_y_bound = grid.y_offset;
        int upper_y_bound = grid.y_offset + grid.ny_local - 1;

        std::cout << "On MPI rank " << id << ", the Y bounds (in cells) are [" << lower_y_bound << "," << upper_y_bound
                  << "]" << std::endl;
        int number_of_temperature_data_points = 0;
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
            // rank within raw_temperature_data and incrementing number_of_temperature_data_points appropriately
            bool binary_input_data = checkTemperatureFileFormat(tempfile_thislayer);
            parseTemperatureData(tempfile_thislayer, grid, number_of_temperature_data_points, binary_input_data);
            last_value(layer_read_count) = number_of_temperature_data_points;
        } // End loop over all files read for all layers
        Kokkos::resize(raw_temperature_data, number_of_temperature_data_points, num_temperature_components);
        // Determine start values for each layer's data within raw_temperature_data, if all layers were read
        if (!(_inputs.layerwise_temp_read)) {
            if (grid.number_of_layers > _inputs.temp_files_in_series) {
                for (int layer_read_count = _inputs.temp_files_in_series; layer_read_count < grid.number_of_layers;
                     layer_read_count++) {
                    if (_inputs.temp_files_in_series == 1) {
                        // All layers have the same temperature data
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

    // Initialize temperature data with a fixed thermal gradient in Z (can also be zero) for Directional/SingleGrain
    // problem types
    void initialize(const int id, const std::string simulation_type, const Grid &grid, const double deltat) {

        // Initialize temperature field in Z direction with thermal gradient G set in input file
        // Liquidus front (InitUndercooling = 0) is at domain bottom for directional solidification, is at domain center
        // (with custom InitUndercooling value) for single grain solidification
        int location_init_undercooling, location_liquidus_isotherm;
        if (simulation_type == "Directional")
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
        auto liquidus_time_local = liquidus_time;
        auto cooling_rate_local = cooling_rate;
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
                liquidus_time_local(index, 0, 0) = -1;
                // Negative dist_from_liquidus and dist_from_init_undercooling values for cells below the liquidus
                // isotherm
                int coord_z = grid.getCoordZ(index);
                int dist_from_liquidus = coord_z - location_liquidus_isotherm;
                // Cells reach liquidus at a time dependent on their Z coordinate
                // Cells with negative liquidus time values are already undercooled, should have positive undercooling
                // and negative liquidus time step
                if (dist_from_liquidus < 0) {
                    liquidus_time_local(index, 0, 1) = -1;
                    int dist_from_init_undercooling = coord_z - location_init_undercooling;
                    undercooling_current_local(index) =
                        init_undercooling_local - dist_from_init_undercooling * G_local * grid.deltax;
                }
                else {
                    // Cells with positive liquidus time values are not yet tracked - leave current undercooling as
                    // default zeros and set liquidus time step. R_local will never be zero here, as all cells at a
                    // fixed undercooling must be below the liquidus (distFromLiquidus < 0)
                    liquidus_time_local(index, 0, 1) = dist_from_liquidus * G_local * grid.deltax / (R_local * deltat);
                }
                // Cells cool at a constant rate
                cooling_rate_local(index, 0) = R_local * deltat;
                // All cells solidify once
                max_solidification_events_local(0) = 1;
                number_solidification_events_local(index) = 1;
            });
        MPI_Barrier(MPI_COMM_WORLD);
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

        // Initialize liquidus_time and cooling_rate values for this spot
        auto liquidus_time_local = liquidus_time;
        auto cooling_rate_local = cooling_rate;
        auto number_of_solidification_events_local = number_of_solidification_events;
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
                    liquidus_time_local(index, 0, 0) = 1;
                    // Liquidus time step (related to distance to spot edge). Must be at least 1 time step after melt
                    // time step
                    liquidus_time_local(index, 0, 1) = Kokkos::round((spot_radius - tot_dist) / isotherm_velocity) + 2;
                    // Cooling rate per time step
                    cooling_rate_local(index, 0) = R_local * deltat;
                    number_of_solidification_events_local(index) = 1;
                }
                else {
                    // Dummy layer_time_temp_history data
                    liquidus_time_local(index, 0, 0) = -1;
                    liquidus_time_local(index, 0, 1) = -1;
                    cooling_rate_local(index, 0) = -1.0;
                    number_of_solidification_events_local(index) = 0;
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
                                     const int end_range, const Grid &grid, const std::string simulation_type) {

        bool calc_remelting_events;
        if (simulation_type == "FromFile") {
            if (layernumber > _inputs.temp_files_in_series)
                calc_remelting_events = false;
            else
                calc_remelting_events = true;
        }
        else {
            if (layernumber > 0)
                calc_remelting_events = false;
            else
                calc_remelting_events = true;
        }
        if (!calc_remelting_events) {
            // Use the value from a previously checked layer, since the time-temperature history is reused
            if ((simulation_type == "FromFinch") || (_inputs.temp_files_in_series == 1)) {
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
            // Need to calculate max_solidification_events(layernumber)
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
        MPI_Barrier(MPI_COMM_WORLD);
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
        double t_melting = raw_temperature_data(i, 3);
        return t_melting;
    }
    // Read data from storage, obtain liquidus time
    double getTempCoordTL(const int i) {
        double t_liquidus = raw_temperature_data(i, 4);
        return t_liquidus;
    }
    // Read data from storage, obtain cooling rate
    double getTempCoordCR(const int i) {
        double cooling_rate = raw_temperature_data(i, 5);
        return cooling_rate;
    }

    // Initialize temperature fields for layer "layernumber" in case where temperature data comes from file(s)
    // TODO: This can be performed on the device as the dirS problem is
    void initialize(const int layernumber, const int id, const Grid &grid, const double freezing_range,
                    const double deltat, const std::string simulation_type) {

        // Data was already read into the "raw_temperature_data" data structure
        // Determine which section of "raw_temperature_data" is relevant for this layer of the overall domain
        int start_range = first_value[layernumber];
        int end_range = last_value[layernumber];

        // Temporary host view for the maximum number of times a cell in a given layer will solidify (different for each
        // layer)
        view_type_int_host max_solidification_events_host =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), max_solidification_events);
        calcMaxSolidificationEvents(id, layernumber, max_solidification_events_host, start_range, end_range, grid,
                                    simulation_type);
        int max_num_solidification_events = max_solidification_events_host(layernumber);

        // Resize liquidus_time now that the max number of solidification events is known for this layer
        Kokkos::resize(liquidus_time, grid.domain_size, max_num_solidification_events, 2);

        // These views are initialized to zeros on the host, filled with data, and then copied to the device for layer
        // "layernumber"
        view_type_int_3d_host liquidus_time_host("LiquidusTime_H", grid.domain_size, max_num_solidification_events, 2);
        view_type_float_2d_host cooling_rate_host("CoolingRate_H", grid.domain_size, max_num_solidification_events);
        view_type_int_host number_of_solidification_events_host("NumSEvents_H", grid.domain_size);

        double largest_time = 0;
        double largest_time_global = 0;
        if (id == 0)
            std::cout << "Range of raw data for layer " << layernumber << " on rank 0 is " << start_range << " to "
                      << end_range << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        for (int i = start_range; i < end_range; i++) {
            // Get the integer X, Y, Z coordinates associated with this data point, along with the associated TM,
            // TL, CR values coord_y is relative to this MPI rank's grid, while coord_y_global is relative to the
            // overall simulation domain
            int coord_x = getTempCoordX(i, grid.x_min, grid.deltax);
            int coord_y = getTempCoordY(i, grid.y_min, grid.deltax, grid.y_offset);
            int coord_z = getTempCoordZ(i, grid.deltax, grid.layer_height, layernumber, grid.z_min_layer);
            double t_melting = getTempCoordTM(i);
            double t_liquidus = getTempCoordTL(i);
            double cooling_rate_cell = getTempCoordCR(i);

            // 1D cell coordinate on this MPI rank's domain
            int index = grid.get1DIndex(coord_x, coord_y, coord_z);
            // Store TM, TL, CR values for this solidification event in liquidus_time_host/cooling_rate_host
            const int n_solidification_events_cell = number_of_solidification_events_host(index);
            liquidus_time_host(index, n_solidification_events_cell, 0) = Kokkos::round(t_melting / deltat) + 1;
            liquidus_time_host(index, n_solidification_events_cell, 1) = Kokkos::round(t_liquidus / deltat) + 1;
            // Cannot go above/below liquidus on same time step - must offset by at least 1
            if (liquidus_time_host(index, n_solidification_events_cell, 0) ==
                liquidus_time_host(index, n_solidification_events_cell, 1))
                liquidus_time_host(index, n_solidification_events_cell, 1)++;
            cooling_rate_host(index, n_solidification_events_cell) = std::abs(cooling_rate_cell) * deltat;
            // Increment number of solidification events for this cell
            number_of_solidification_events_host(index)++;
            // Estimate of the time step where the last possible solidification is expected to occur
            double solidus_time = t_liquidus + freezing_range / cooling_rate_cell;
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

        // Reorder solidification events in liquidus_time_host(location,event number,component) and
        // cooling_rate_host(location,event number) so that they are in order based on the melting time values
        // (component = 0 in liquidus_time_host)
        for (int index = 0; index < grid.domain_size; index++) {
            int n_solidification_events_cell = number_of_solidification_events_host(index);
            if (n_solidification_events_cell > 0) {
                for (int i = 0; i < n_solidification_events_cell - 1; i++) {
                    for (int j = (i + 1); j < n_solidification_events_cell; j++) {
                        if (liquidus_time_host(index, i, 0) > liquidus_time_host(index, j, 0)) {
                            // Swap these two points - melting event "j" happens before event "i"
                            float old_melt_val = liquidus_time_host(index, i, 0);
                            float old_liq_val = liquidus_time_host(index, i, 1);
                            float old_cr_val = cooling_rate_host(index, i);
                            liquidus_time_host(index, i, 0) = liquidus_time_host(index, j, 0);
                            liquidus_time_host(index, i, 1) = liquidus_time_host(index, j, 1);
                            cooling_rate_host(index, i) = cooling_rate_host(index, j);
                            liquidus_time_host(index, j, 0) = old_melt_val;
                            liquidus_time_host(index, j, 1) = old_liq_val;
                            cooling_rate_host(index, j) = old_cr_val;
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
                    if (liquidus_time_host(index, i + 1, 0) < liquidus_time_host(index, i, 1)) {
                        std::cout << "Cell " << index << " removing anomalous event " << i + 1 << " out of "
                                  << n_solidification_events_cell - 1 << std::endl;
                        // Keep whichever event has the larger liquidus time
                        if (liquidus_time_host(index, i + 1, 1) > liquidus_time_host(index, i, 1)) {
                            liquidus_time_host(index, i, 0) = liquidus_time_host(index, i + 1, 0);
                            liquidus_time_host(index, i, 1) = liquidus_time_host(index, i + 1, 1);
                            cooling_rate_host(index, i) = cooling_rate_host(index, i + 1);
                        }
                        liquidus_time_host(index, i + 1, 0) = 0.0;
                        liquidus_time_host(index, i + 1, 1) = 0.0;
                        cooling_rate_host(index, i + 1) = 0.0;
                        // Reshuffle other solidification events over if needed
                        for (int ii = (i + 1); ii < n_solidification_events_cell - 1; ii++) {
                            liquidus_time_host(index, ii, 0) = liquidus_time_host(index, ii + 1, 0);
                            liquidus_time_host(index, ii, 1) = liquidus_time_host(index, ii + 1, 1);
                            cooling_rate_host(index, ii) = cooling_rate_host(index, ii + 1);
                        }
                        number_of_solidification_events_host(index)--;
                    }
                }
            }
        }

        // Copy host view data back to device
        max_solidification_events = Kokkos::create_mirror_view_and_copy(memory_space(), max_solidification_events_host);
        liquidus_time = Kokkos::create_mirror_view_and_copy(memory_space(), liquidus_time_host);
        cooling_rate = Kokkos::create_mirror_view_and_copy(memory_space(), cooling_rate_host);
        number_of_solidification_events =
            Kokkos::create_mirror_view_and_copy(memory_space(), number_of_solidification_events_host);
        MPI_Barrier(MPI_COMM_WORLD);
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
        undercooling_current(index) += cooling_rate(index, solidification_event_counter(index));
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
        melt_time_step = liquidus_time(index, solidification_event_counter_cell, 0);
        if (cycle > melt_time_step) {
            // If the cell has already exceeded the melt time step for the current melt-solidification event, get the
            // melt time step associated with the next solidification event - or, if there is no next
            // melt-solidification event, return the max possible int as the cell will not melt again during this layer
            // of the multilayer problem
            if (solidification_event_counter_cell < (number_of_solidification_events(index) - 1))
                melt_time_step = liquidus_time(index, solidification_event_counter_cell + 1, 0);
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
        int crit_time_step = liquidus_time(index, solidification_event_counter_cell, 1);
        return crit_time_step;
    }
    // Uses a specified solidification event
    KOKKOS_INLINE_FUNCTION
    int getCritTimeStep(const int index, const int solidification_event_counter_cell) const {
        int crit_time_step = liquidus_time(index, solidification_event_counter_cell, 1);
        return crit_time_step;
    }

    // Extract the cooling rate associated with a specified solidificaiton event
    KOKKOS_INLINE_FUNCTION
    float getUndercoolingChange(const int index, const int solidification_event_counter_cell) const {
        float undercooling_change = cooling_rate(index, solidification_event_counter_cell);
        return undercooling_change;
    }

    // Extract either the last time step that all points undergo melting in the layer or the last time they cools below
    // the liquidus from liquidus_time (corresponds to solidification event number_of_solidification_events-1
    // for the cell) (can't just use subview here since number_of_solidification_events is different for each cell) If
    // the cell does not undergo solidification, either print -1 or the specified default value
    template <typename extracted_view_data_type>
    extracted_view_data_type extractTmTlData(const int extracted_val, const int domain_size,
                                             const int default_val = -1) {
        extracted_view_data_type extracted_data(Kokkos::ViewAllocateWithoutInitializing("extracted_data"), domain_size);

        // Local copies for lambda capture.
        auto liquidus_time_local = liquidus_time;
        auto number_of_solidification_events_local = number_of_solidification_events;

        auto policy = Kokkos::RangePolicy<execution_space>(0, domain_size);
        Kokkos::parallel_for(
            "Extract_tm_tl_data", policy, KOKKOS_LAMBDA(const int &index) {
                // If this cell doesn't undergo solidification at all, print -1
                const int num_solidification_events_this_cell = number_of_solidification_events_local(index);
                if (num_solidification_events_this_cell == 0)
                    extracted_data(index) = default_val;
                else
                    extracted_data(index) =
                        liquidus_time_local(index, num_solidification_events_this_cell - 1, extracted_val);
            });
        return extracted_data;
    }

    // Extract the cooling rate corresponding to solidification event number `number_of_solidification_events-1` for the
    // cell (can't just use subview here since number_of_solidification_events is different for each cell) If the cell
    // does not undergo solidification, either print -1 or the specified default value
    template <typename extracted_view_data_type>
    extracted_view_data_type extractCrData(const int domain_size, const int default_val = -1) {
        extracted_view_data_type extracted_data(Kokkos::ViewAllocateWithoutInitializing("extracted_data"), domain_size);

        // Local copies for lambda capture.
        auto cooling_rate_local = cooling_rate;
        auto number_of_solidification_events_local = number_of_solidification_events;

        auto policy = Kokkos::RangePolicy<execution_space>(0, domain_size);
        Kokkos::parallel_for(
            "Extract_tm_tl_data", policy, KOKKOS_LAMBDA(const int &index) {
                // If this cell doesn't undergo solidification at all, print -1
                const int num_solidification_events_this_cell = number_of_solidification_events_local(index);
                if (num_solidification_events_this_cell == 0)
                    extracted_data(index) = default_val;
                else
                    extracted_data(index) = cooling_rate_local(index, num_solidification_events_this_cell - 1);
            });
        return extracted_data;
    }
};

#endif
