// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_GRID_HPP
#define EXACA_GRID_HPP

#include "CAinputs.hpp"
#include "CAparsefiles.hpp"

#include "mpi.h"

#include <Kokkos_Core.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Data for local (individual MPI ranks) and global (across all ranks) grids
struct Grid {

    // Global domain bounds, all layers
    int nx, ny, nz, domain_size_all_layers;
    double deltax, x_min, x_max, y_min, y_max, z_min, z_max;

    // Variables characterizing local processor grids relative to global domain
    // 1D decomposition in Y: Each MPI rank has a subset consisting of of ny_local cells, out of ny cells in Y
    // Each MPI rank's subdomain is offset by y_offset cells from the lower bound of the domain in Y
    int ny_local, y_offset;
    // Variables characterizing process IDs of neighboring MPI ranks on the grid
    // Positive Y/NegativeY directions are North/South
    int neighbor_rank_north, neighbor_rank_south;
    // Variables denoting whether or not each MPI rank's grid is at a global domain boundary
    bool at_north_boundary, at_south_boundary;

    // Multilayer problem information
    using view_type_double_host = Kokkos::View<double *, Kokkos::HostSpace>;
    view_type_double_host z_min_layer, z_max_layer;
    int number_of_layers, layer_height;
    double HT_deltax;

    // Current layer information for multilayer problems
    int z_layer_bottom, z_layer_top, nz_layer, domain_size, bottom_of_current_layer, top_of_current_layer;
    std::pair<int, int> layer_range;

    // Domain inputs from file
    DomainInputs _inputs;
    // Temperature inputs from file
    TemperatureInputs _t_inputs;

    // Creates grid struct with uninitialized values, used in unit tests
    Grid(int number_of_layers_temp = 1)
        : z_min_layer(
              view_type_double_host(Kokkos::ViewAllocateWithoutInitializing("z_min_layer"), number_of_layers_temp))
        , z_max_layer(
              view_type_double_host(Kokkos::ViewAllocateWithoutInitializing("z_max_layer"), number_of_layers_temp)) {
        number_of_layers = number_of_layers_temp;
    };

    // Constructor for grid used in ExaCA
    Grid(const std::string SimulationType, const int id, const int np, const int number_of_layers_temp,
         DomainInputs inputs, TemperatureInputs t_inputs)
        : z_min_layer(
              view_type_double_host(Kokkos::ViewAllocateWithoutInitializing("z_min_layer"), number_of_layers_temp))
        , z_max_layer(
              view_type_double_host(Kokkos::ViewAllocateWithoutInitializing("z_max_layer"), number_of_layers_temp))
        , _inputs(inputs)
        , _t_inputs(t_inputs) {

        // Copy from inputs structs
        deltax = _inputs.deltax;
        layer_height = _inputs.LayerHeight;
        layer_height = _inputs.LayerHeight;
        number_of_layers = number_of_layers_temp;

        // Obtain global domain bounds
        // For problem type R (reading data from file), need to parse temperature data files to obtain domain bounds for
        // each layer and for the multilayer domain
        if (SimulationType == "R") {
            // For simulations using input temperature data with remelting: even if only LayerwiseTempRead is true, all
            // files need to be read to determine the domain bounds
            find_xyz_bounds(id);
        }
        else {
            // Copy inputs from inputs struct into grid struct
            nx = _inputs.nx;
            ny = _inputs.ny;
            nz = _inputs.nz;
            // Domain origin is placed at 0
            x_min = 0.0;
            y_min = 0.0;
            z_min = 0.0;
            x_max = nx * deltax;
            y_max = ny * deltax;
            z_max = nz * deltax;
            if (SimulationType == "S") {
                // If this is a spot melt problem, also set the z_min/z_max for each layer
                for (int n = 0; n < number_of_layers; n++) {
                    z_min_layer(n) = deltax * (layer_height * n);
                    z_max_layer(n) = deltax * (_inputs.SpotRadius + layer_height * n);
                }
            }
        }
        if (id == 0) {
            std::cout << "Domain size: " << nx << " by " << ny << " by " << nz << std::endl;
            std::cout << "X Limits of domain: " << x_min << " and " << x_max << std::endl;
            std::cout << "Y Limits of domain: " << y_min << " and " << y_max << std::endl;
            std::cout << "Z Limits of domain: " << z_min << " and " << z_max << std::endl;
            std::cout << "================================================================" << std::endl;
        }

        // Decompose the domain into subdomains on each MPI rank: Calculate ny_local and y_offset for each rank, where
        // each subdomain contains "ny_local" in Y, offset from the full domain origin by "y_offset" cells in Y
        // (previously "DomainDecomposition" in CAinitialize.cpp) First, compare total MPI ranks to total Y cells.
        if (np > ny)
            throw std::runtime_error(
                "Error: Cannot run with more MPI ranks than cells in Y (decomposition direction).");

        // The following was previously performed in "DomainDecomposition" in CAinitialize.cpp:
        // Domain size across all ranks and all layers
        domain_size_all_layers = get_domain_size_all_layers();

        // Get neighboring ranks to each direction and set boundary indicator variables
        neighbor_rank_north = get_neighbor_rank_north(id, np);
        neighbor_rank_south = get_neighbor_rank_south(id, np);
        at_north_boundary = get_at_north_boundary();
        at_south_boundary = get_at_south_boundary();

        // Determine, for each MPI process id, the local grid size in y (and the offset in y relative to the overall
        // simulation domain)
        y_offset = get_y_offset(id, np);
        ny_local = get_ny_local(id, np);

        // Add halo regions with a width of 1 in +/- Y if this MPI rank is not as a domain boundary in said direction
        add_halo();
        std::cout << "Rank " << id << " spans Y = " << y_offset << " through " << y_offset + ny_local - 1 << std::endl;

        // Bounds of layer 0: Z coordinates span z_layer_bottom-z_layer_top, inclusive (functions previously in
        // CAinitialize.cpp and CAcelldata.hpp)
        z_layer_bottom = calc_z_layer_bottom(SimulationType, 0);
        z_layer_top = calc_z_layer_top(SimulationType, _inputs.SpotRadius, 0);
        nz_layer = calc_nz_layer(id, 0);
        domain_size = calc_domain_size(); // Number of cells in the current layer on this MPI rank
        bottom_of_current_layer = get_bottom_of_current_layer();
        top_of_current_layer = get_top_of_current_layer();
        layer_range = std::make_pair(bottom_of_current_layer, top_of_current_layer);
        MPI_Barrier(MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Mesh initialized: initial domain size is " << nz_layer << " out of " << nz
                      << " total cells in the Z direction" << std::endl;
    };

    // For simulation type R, obtain the physical XYZ bounds of the domain by reading temperature data files and parsing
    // the coordinates. Previously in CAinitialize.cpp
    void find_xyz_bounds(const int id) {

        // Two passes through reading temperature data files- the first pass only reads the headers to
        // determine units and X/Y/Z bounds of the simulaton domain. Using the X/Y/Z bounds of the simulation domain,
        // nx, ny, and nz can be calculated and the domain decomposed among MPI processes. The maximum number of
        // remelting events in the simulation can also be calculated. The second pass reads the actual X/Y/Z/liquidus
        // time/cooling rate data and each rank stores the data relevant to itself in "RawData" - this is done in the
        // subroutine "ReadTemperatureData"
        x_min = std::numeric_limits<double>::max();
        y_min = std::numeric_limits<double>::max();
        z_min = std::numeric_limits<double>::max();
        x_max = std::numeric_limits<double>::lowest();
        y_max = std::numeric_limits<double>::lowest();
        z_max = std::numeric_limits<double>::lowest();

        // Read the first temperature file
        std::ifstream FirstTemperatureFile;
        FirstTemperatureFile.open(_t_inputs.temp_paths[0]);
        std::string FirstLineFirstFile;
        // Header line
        getline(FirstTemperatureFile, FirstLineFirstFile);

        // Read all data files to determine the domain bounds, max number of remelting events
        // for simulations with remelting
        int layers_to_read = std::min(number_of_layers, _t_inputs.TempFilesInSeries); // was given in input file
        for (int layer_read_count = 1; layer_read_count <= layers_to_read; layer_read_count++) {

            std::string tempfile_thislayer = _t_inputs.temp_paths[layer_read_count - 1];
            // Get min and max x coordinates in this file, which can be a binary or ASCII input file
            // binary file type uses extension .catemp, all other file types assumed to be comma-separated ASCII input
            bool binary_input_data = checkTemperatureFileFormat(tempfile_thislayer);
            // { x_min, x_max, y_min, y_max, z_min, z_max }
            std::array<double, 6> xyz_min_max_this_layer =
                parseTemperatureCoordinateMinMax(tempfile_thislayer, binary_input_data);

            // Based on the input file's layer offset, adjust z_min/z_max from the temperature data coordinate
            // system to the multilayer CA coordinate system Check to see in the XYZ bounds for this layer are
            // also limiting for the entire multilayer CA coordinate system
            xyz_min_max_this_layer[4] += deltax * layer_height * (layer_read_count - 1);
            xyz_min_max_this_layer[5] += deltax * layer_height * (layer_read_count - 1);
            if (xyz_min_max_this_layer[0] < x_min)
                x_min = xyz_min_max_this_layer[0];
            if (xyz_min_max_this_layer[1] > x_max)
                x_max = xyz_min_max_this_layer[1];
            if (xyz_min_max_this_layer[2] < y_min)
                y_min = xyz_min_max_this_layer[2];
            if (xyz_min_max_this_layer[3] > y_max)
                y_max = xyz_min_max_this_layer[3];
            if (xyz_min_max_this_layer[4] < z_min)
                z_min = xyz_min_max_this_layer[4];
            if (xyz_min_max_this_layer[5] > z_max)
                z_max = xyz_min_max_this_layer[5];
            z_min_layer[layer_read_count - 1] = xyz_min_max_this_layer[4];
            z_max_layer[layer_read_count - 1] = xyz_min_max_this_layer[5];
            if (id == 0)
                std::cout << "Layer = " << layer_read_count << " Z Bounds are " << xyz_min_max_this_layer[4] << " "
                          << xyz_min_max_this_layer[5] << std::endl;
        }
        // Extend domain in Z (build) direction if the number of layers are simulated is greater than the number
        // of temperature files read
        if (number_of_layers > _t_inputs.TempFilesInSeries) {
            for (int layer_read_count = _t_inputs.TempFilesInSeries; layer_read_count < number_of_layers;
                 layer_read_count++) {
                if (_t_inputs.TempFilesInSeries == 1) {
                    // Only one temperature file was read, so the upper Z bound should account for an additional
                    // "number_of_layers-1" worth of data Since all layers have the same temperature data, each
                    // layer's "z_min_layer" is just translated from that of the first layer
                    z_min_layer[layer_read_count] = z_min_layer[layer_read_count - 1] + deltax * layer_height;
                    z_max_layer[layer_read_count] = z_max_layer[layer_read_count - 1] + deltax * layer_height;
                    z_max += deltax * layer_height;
                }
                else {
                    // "TempFilesInSeries" temperature files was read, so the upper Z bound should account for
                    // an additional "number_of_layers-TempFilesInSeries" worth of data
                    int RepeatedFile = (layer_read_count) % _t_inputs.TempFilesInSeries;
                    int RepeatUnit = layer_read_count / _t_inputs.TempFilesInSeries;
                    z_min_layer[layer_read_count] =
                        z_min_layer[RepeatedFile] + RepeatUnit * _t_inputs.TempFilesInSeries * deltax * layer_height;
                    z_max_layer[layer_read_count] =
                        z_max_layer[RepeatedFile] + RepeatUnit * _t_inputs.TempFilesInSeries * deltax * layer_height;
                    z_max += deltax * layer_height;
                }
            }
        }

        // Now at the conclusion of "Loop 0", the decomposition can be performed as the domain bounds are known
        // (all header lines from all files have been read)
        // CA cells in each direction span from the lower to the higher bound of the temperature data - without wall
        // cells or padding around the simulation edges
        nx = round((x_max - x_min) / deltax) + 1;
        ny = round((y_max - y_min) / deltax) + 1;
        nz = round((z_max - z_min) / deltax) + 1;
    }

    int get_domain_size_all_layers() {
        int domain_size_all_layers_local = nx * ny * nz;
        return domain_size_all_layers_local;
    }

    // Determine the MPI rank of the processor in the +Y direction
    int get_neighbor_rank_north(const int id, const int np) {
        int neighbor_rank_north_local;
        if (np > 1) {
            neighbor_rank_north_local = id + 1;
            if (id == np - 1)
                neighbor_rank_north_local = MPI_PROC_NULL;
        }
        else {
            // No MPI communication
            neighbor_rank_north_local = MPI_PROC_NULL;
        }
        return neighbor_rank_north_local;
    }

    // Determine the MPI rank of the processor in the -Y direction
    int get_neighbor_rank_south(const int id, const int np) {
        int neighbor_rank_south_local;
        if (np > 1) {
            neighbor_rank_south_local = id - 1;
            if (id == 0)
                neighbor_rank_south_local = MPI_PROC_NULL;
        }
        else {
            // No MPI communication
            neighbor_rank_south_local = MPI_PROC_NULL;
        }
        return neighbor_rank_south_local;
    }

    // Determine if this MPI rank is at the domain boundary in the +Y direction
    bool get_at_north_boundary() {
        bool at_north_boundary_local;
        if (neighbor_rank_north == MPI_PROC_NULL)
            at_north_boundary_local = true;
        else
            at_north_boundary_local = false;
        return at_north_boundary_local;
    }

    // Determine if this MPI rank is at the domain boundary in the -Y direction
    bool get_at_south_boundary() {
        bool at_south_boundary_local;
        if (neighbor_rank_south == MPI_PROC_NULL)
            at_south_boundary_local = true;
        else
            at_south_boundary_local = false;
        return at_south_boundary_local;
    }

    // Subdivide ny into ny_local across np ranks as evenly as possible, return ny_local on rank id
    int get_ny_local(const int id, const int np) const {
        int ny_local_local = 0;
        int ny_local_est = ny / np;
        int y_remainder = ny % np;
        if (y_remainder == 0) {
            ny_local_local = ny_local_est;
        }
        else {
            if (y_remainder > id) {
                ny_local_local = ny_local_est + 1;
            }
            else {
                ny_local_local = ny_local_est;
            }
        }
        return ny_local_local;
    }

    // Get the offset in the Y direction of the local grid on rank id from the global grid of np ranks, which has ny
    // total cells in Y
    int get_y_offset(const int id, const int np) const {
        int y_offset_local = 0;
        int y_offset_est = ny / np;
        int y_remainder = ny % np;
        if (y_remainder == 0) {
            y_offset_local = id * y_offset_est;
        }
        else {
            if (y_remainder > id) {
                y_offset_local = id * (y_offset_est + 1);
            }
            else {
                y_offset_local =
                    (y_offset_est + 1) * (y_remainder - 1) + y_offset_est + 1 + (id - y_remainder) * y_offset_est;
            }
        }
        return y_offset_local;
    }

    // Add ghost nodes to the appropriate subdomains (added where the subdomains overlap, but not at edges of physical
    // domain)
    void add_halo() {

        // Add halo regions in Y direction if this subdomain borders subdomains on other processors
        // If only 1 rank in the y direction, no halo regions - subdomain is coincident with overall simulation domain
        // If multiple ranks in the y direction, either 1 halo region (borders another rank's subdomain in either the +y
        // or -y direction) or 2 halo regions (if it borders other rank's subdomains in both the +y and -y directions)
        if (neighbor_rank_north != MPI_PROC_NULL)
            ny_local++;
        if (neighbor_rank_south != MPI_PROC_NULL) {
            ny_local++;
            // Also adjust subdomain offset, as these ghost nodes were added on the -y side of the subdomain
            y_offset--;
        }
    }

    // Get the Z coordinate of the lower bound of iteration for layer layernumber
    int calc_z_layer_bottom(const std::string SimulationType, const int layernumber) {

        int z_layer_bottom_local = -1; // assign dummy initial value
        if ((SimulationType == "C") || (SimulationType == "SingleGrain")) {
            // Not a multilayer problem, top of "layer" is the top of the overall simulation domain
            z_layer_bottom_local = 0;
        }
        else if (SimulationType == "S") {
            // lower bound of domain is an integer multiple of the layer spacing, since the temperature field is the
            // same for every layer
            z_layer_bottom_local = layer_height * layernumber;
        }
        else if (SimulationType == "R") {
            // lower bound of domain is based on the data read from the file(s)
            z_layer_bottom_local = round((z_min_layer[layernumber] - z_min) / deltax);
        }
        if (z_layer_bottom_local == -1)
            throw std::runtime_error("Error: ZBound_Low went uninitialized, problem type must be C, S, or R");
        return z_layer_bottom_local;
    }

    // Get the Z coordinate of the upper bound of iteration for layer layernumber
    int calc_z_layer_top(const std::string SimulationType, const int SpotRadius, const int layernumber) {

        int z_layer_top_local = -1; // assign dummy initial value
        if ((SimulationType == "C") || (SimulationType == "SingleGrain")) {
            // Not a multilayer problem, top of "layer" is the top of the overall simulation domain
            z_layer_top_local = nz - 1;
        }
        else if (SimulationType == "S") {
            // Top of layer is equal to the spot radius for a problem of hemispherical spot solidification, plus an
            // offset depending on the layer number
            z_layer_top_local = SpotRadius + layer_height * layernumber;
        }
        else if (SimulationType == "R") {
            // Top of layer comes from the layer's file data (implicitly assumes bottom of layer 0 is the bottom of the
            // overall domain - this should be fixed in the future for edge cases where this isn't true)
            z_layer_top_local = round((z_max_layer[layernumber] - z_min) / deltax);
        }
        if (z_layer_top_local == -1)
            throw std::runtime_error(
                "Error: ZBound_High went uninitialized, problem type must be SingleGrain, C, S, or R");
        return z_layer_top_local;
    }

    // Calculate the size of the active domain in Z for layer layernumber
    int calc_nz_layer(const int id, const int layernumber) {
        int nz_layer_local = z_layer_top - z_layer_bottom + 1;
        if (id == 0)
            std::cout << "Layer " << layernumber << "'s active domain is from Z = " << z_layer_bottom << " through "
                      << z_layer_top << " (" << nz_layer_local << ") cells" << std::endl;
        return nz_layer_local;
    }

    // Calculate the size of the domain, as a number of cells, for this layer of a multilayer problem
    int calc_domain_size() {
        int domain_size_local = nx * ny_local * nz_layer;
        return domain_size_local;
    }

    // Get 1D coordinate corresponding to the first cell in this layer of a multilayer problem
    int get_bottom_of_current_layer() {
        int bottom_of_current_layer_local = z_layer_bottom * nx * ny_local;
        return bottom_of_current_layer_local;
    }

    // Get 1D coordinate corresponding to the cell after the last cell in this layer of a multilayer problem
    int get_top_of_current_layer() {
        int top_of_current_layer_local = z_layer_bottom * nx * ny_local + domain_size;
        return top_of_current_layer_local;
    }

    // Determine new active cell domain size and offset from bottom of global domain
    void init_next_layer(const int id, const std::string SimulationType, const int next_layer_number,
                         const int SpotRadius) {
        z_layer_bottom = calc_z_layer_bottom(SimulationType, next_layer_number);
        z_layer_top = calc_z_layer_top(SimulationType, SpotRadius, next_layer_number);
        nz_layer = calc_nz_layer(id, next_layer_number);
        domain_size = calc_domain_size();
        bottom_of_current_layer = get_bottom_of_current_layer();
        top_of_current_layer = get_top_of_current_layer();
        layer_range = std::make_pair(bottom_of_current_layer, top_of_current_layer);
    }

    // Get the 1D cell coordinate from the x, y, and z cell positions
    KOKKOS_INLINE_FUNCTION
    int get_1D_index(const int coord_x, const int coord_y_local, const int coord_z) const {
        int index = coord_z * nx * ny_local + coord_x * ny_local + coord_y_local;
        return index;
    }
    // TODO: There is probably some creative way to combine these functions and return one object containing the x, y,
    // and z positions Get the z cell position of the cell from the 1D cell coordinate
    KOKKOS_INLINE_FUNCTION
    int get_coord_Z(const int index) const {
        int coord_z = index / (nx * ny_local);
        return coord_z;
    }
    // Get the y cell position of the cell from the 1D cell coordinate
    KOKKOS_INLINE_FUNCTION
    int get_coord_Y(const int index) const {
        int rem = index % (nx * ny_local);
        int coord_y = rem % ny_local;
        return coord_y;
    }
    // Get the x cell position of the cell from the 1D cell coordinate
    KOKKOS_INLINE_FUNCTION
    int get_coord_X(const int index) const {
        int rem = index % (nx * ny_local);
        int coord_x = rem / ny_local;
        return coord_x;
    }
    // Get the z cell position of the cell from the 1D cell coordinate with respect to the overall simulation domain
    // (all MPI ranks)
    KOKKOS_INLINE_FUNCTION
    int get_coord_Z_global(const int index) const {
        int coord_z = index / (nx * ny);
        return coord_z;
    }
    // Get the y cell position of the cell from the 1D cell coordinate with respect to the overall simulation domain
    // (all MPI ranks)
    KOKKOS_INLINE_FUNCTION
    int get_coord_Y_global(const int index) const {
        int rem = index % (nx * ny);
        int coord_y = rem % ny;
        return coord_y;
    }
    // Get the x cell position of the cell from the 1D cell coordinate with respect to the overall simulation domain
    // (all MPI ranks)
    KOKKOS_INLINE_FUNCTION
    int get_coord_X_global(const int index) const {
        int rem = index % (nx * ny);
        int coord_x = rem / ny;
        return coord_x;
    }
};

#endif
