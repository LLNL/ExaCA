// Copyright Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_CELLDATA_HPP
#define EXACA_CELLDATA_HPP

#include "CAconfig.hpp"
#include "CAgrid.hpp"
#include "CAinputdata.hpp"
#include "CAparsefiles.hpp"
#include "CAtypes.hpp"
#include "mpi.h"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

// Data that belongs to individual cells governing the current state and various ID values
template <typename MemorySpace>
struct CellData {

    using memory_space = MemorySpace;
    using view_type_int = Kokkos::View<int *, memory_space>;
    using view_type_int_host = typename view_type_int::HostMirror;
    using device_layout = typename view_type_int::array_layout;
    using view_type_int_2d_host = Kokkos::View<int **, device_layout, Kokkos::HostSpace>;
    using view_type_int_unmanaged = Kokkos::View<int *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using view_type_bool = Kokkos::View<bool *, memory_space>;
    using view_type_short = Kokkos::View<short *, memory_space>;
    using view_type_float = Kokkos::View<float *, memory_space>;

    // Using the default exec space for this memory space.
    using execution_space = typename memory_space::execution_space;

    int next_layer_first_epitaxial_grain_id;
    view_type_int grain_id_all_layers, cell_type;
    view_type_short layer_id_all_layers;
    view_type_bool melt_edge, melt_edge_all_layers;
    // Substrate inputs from file
    SubstrateInputs _inputs;
    // Storing of whether or not a cell is at an edge of a melt pool (FromFile and FromFinch problem types only)
    bool _store_melt_pool_edge;

    // Constructor for views and view bounds for current layer
    // GrainID is initialized to zeros, while others are not initialized
    // cell_type only exists for the current layer of a multilayer problem
    CellData(const Grid &grid, SubstrateInputs inputs, const bool store_melt_pool_edge = false)
        : grain_id_all_layers(view_type_int("GrainID", grid.domain_size_all_layers))
        , cell_type(view_type_int(Kokkos::ViewAllocateWithoutInitializing("cell_type"), grid.domain_size))
        , layer_id_all_layers(
              view_type_short(Kokkos::ViewAllocateWithoutInitializing("LayerID"), grid.domain_size_all_layers))
        , _inputs(inputs)
        , _store_melt_pool_edge(store_melt_pool_edge) {
        if (_store_melt_pool_edge) {
            // Default init to zero
            melt_edge_all_layers = view_type_bool("melt_edge", grid.domain_size_all_layers);
            // Current layer
            getCurrentLayerMeltEdge(grid.layer_range);
        }
    }

    // Get the subview associated with the edge indicator of cells in the current layer
    void getCurrentLayerMeltEdge(std::pair<int, int> layer_range) {
        melt_edge = Kokkos::subview(melt_edge_all_layers, layer_range);
    }

    // Initializes the single active cell and associated active cell data structures for the single grain at the domain
    // center
    void initSubstrate(const int id, const Grid &grid) {

        // Location of the single grain
        int grain_location_x = Kokkos::floorf(static_cast<float>(grid.nx) / 2.0);
        int grain_location_y = Kokkos::floorf(static_cast<float>(grid.ny) / 2.0);
        int grain_location_z = Kokkos::floorf(static_cast<float>(grid.nz) / 2.0);

        // Local copies for lambda capture.
        auto cell_type_local = cell_type;
        auto grain_id_all_layers_local = grain_id_all_layers;
        int single_grain_orientation_local = _inputs.single_grain_orientation;

        auto policy = Kokkos::RangePolicy<execution_space>(0, grid.domain_size);
        Kokkos::parallel_for(
            "SingleGrainInit", policy, KOKKOS_LAMBDA(const int &index) {
                int coord_x = grid.getCoordX(index);
                int coord_y = grid.getCoordY(index);
                int coord_z = grid.getCoordZ(index);
                int coord_y_global = coord_y + grid.y_offset;
                if ((coord_x == grain_location_x) && (coord_y_global == grain_location_y) &&
                    (coord_z == grain_location_z)) {
                    cell_type_local(index) = FutureActive;
                    grain_id_all_layers_local(index) = single_grain_orientation_local + 1;
                }
                else
                    cell_type_local(index) = Liquid;
            });
        if ((grain_location_y >= grid.y_offset) && (grain_location_y < grid.y_offset + grid.ny_local))
            std::cout << "Rank " << id << " initialized a grain with orientation " << single_grain_orientation_local
                      << " initialized at X = " << grain_location_x << ", Y = " << grain_location_y
                      << ", Z = " << grain_location_z << std::endl;

        MPI_Barrier(MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Grain struct initialized" << std::endl;
    }

    // Get the X, Y coordinates and grain ID values for grains at the bottom surface for problem type Directional
    view_type_int_2d_host getSurfaceActiveCellData(int &substrate_act_cells, const Grid &grid,
                                                   const unsigned long rng_seed) {

        // Number of cells at the bottom surface that could potentially be assigned GrainID values
        const int bottom_surface_size = grid.nx * grid.ny;

        // First get number of substrate grains for each initialization condition
        if (_inputs.surface_init_mode == "SurfaceSiteFraction")
            substrate_act_cells = Kokkos::round(_inputs.fract_surface_sites_active * bottom_surface_size);
        else if (_inputs.surface_init_mode == "SurfaceSiteDensity")
            substrate_act_cells =
                Kokkos::round(_inputs.surface_site_density * static_cast<double>(bottom_surface_size) * grid.deltax *
                              grid.deltax * pow(10, 12));
        if (_inputs.surface_init_mode == "Custom")
            substrate_act_cells = _inputs.grain_locations_x.size();

        // View for storing surface grain locations and IDs
        view_type_int_2d_host act_cell_data_host(Kokkos::ViewAllocateWithoutInitializing("ActCellData_Host"),
                                                 substrate_act_cells, 3);
        if (_inputs.surface_init_mode == "SurfaceSiteFraction") {
            // Fraction of surface sites active was given - ensure the appropriate number of cells are assigned
            // GrainIDs. Note that the physical locations of these sites (x,y) will vary based on the domain size/cell
            // size Create list of grain IDs and shuffle - leave 0s for cells without substrate grains
            std::vector<int> grain_locations_1d(bottom_surface_size, 0);
            for (int n = 0; n < substrate_act_cells; n++) {
                grain_locations_1d[n] = n + 1; // grain ID for epitaxial seeds must be > 0
            }
            std::mt19937_64 gen(rng_seed);
            std::shuffle(grain_locations_1d.begin(), grain_locations_1d.end(), gen);
            // Fill act_cell_data_host with the non-zero grain IDs and their associated X and Y based on the position in
            // the 1D vector
            int act_cell_count = 0;
            for (int n = 0; n < bottom_surface_size; n++) {
                if (grain_locations_1d[n] != 0) {
                    act_cell_data_host(act_cell_count, 0) = grid.getCoordXGlobal(n);
                    act_cell_data_host(act_cell_count, 1) = grid.getCoordYGlobal(n);
                    act_cell_data_host(act_cell_count, 2) = grain_locations_1d[n];
                    act_cell_count++;
                }
            }
        }
        else if (_inputs.surface_init_mode == "SurfaceSiteDensity") {
            // Calls to Xdist(gen) and Y dist(gen) return random locations for grain seeds
            // Since X = 0 and X = nx-1 are the cell centers of the last cells in X, locations are evenly scattered
            // between X = -0.49999 and X = nx - 0.5, as the cells have a half width of 0.5. Note that if the number of
            // grains is large compared to the number of cells, multiple grain IDs may be assigned to one cell and the
            // total density will be underestimated on the given grid
            std::mt19937_64 gen(rng_seed);
            std::uniform_real_distribution<double> x_dist(-0.49999, grid.nx - 0.5);
            std::uniform_real_distribution<double> y_dist(-0.49999, grid.ny - 0.5);
            // Randomly locate substrate grain seeds for cells in the interior of this subdomain (at the k = 0 bottom
            // surface)
            for (int n = 0; n < substrate_act_cells; n++) {
                double x_location = x_dist(gen);
                double y_location = y_dist(gen);
                // Randomly select integer coordinates between 0 and nx-1 or ny-1
                act_cell_data_host(n, 0) = Kokkos::round(x_location);
                act_cell_data_host(n, 1) = Kokkos::round(y_location);
                act_cell_data_host(n, 2) = n + 1; // grain ID for epitaxial seeds must be > 0
            }
        }
        else if (_inputs.surface_init_mode == "Custom") {
            // Values were already given in the inputs struct, copy these to the view
            for (int n = 0; n < substrate_act_cells; n++) {
                act_cell_data_host(n, 0) = _inputs.grain_locations_x[n];
                act_cell_data_host(n, 1) = _inputs.grain_locations_y[n];
                act_cell_data_host(n, 2) = _inputs.grain_ids[n];
            }
        }
        return act_cell_data_host;
    }

    // Initializes cell types and epitaxial Grain ID values where substrate grains are future active cells on the bottom
    // surface of the constrained domain
    void initSubstrate(const int id, const Grid &grid, const unsigned long rng_seed) {

        // Fill the view of cell X, Y, and ID values, updating the number of substrate active cells appropriately
        // TODO: Could generate random numbers on GPU, instead of using host view and copying over - but would also need
        // inputs struct to store device data for grain locations in X, Y, and GrainIDs
        int substrate_act_cells;
        view_type_int_2d_host act_cell_data_host = getSurfaceActiveCellData(substrate_act_cells, grid, rng_seed);
        // Copy views of substrate grain locations and IDs back to the device
        auto act_cell_data = Kokkos::create_mirror_view_and_copy(memory_space(), act_cell_data_host);

        // Start with all cells as liquid prior to locating substrate grain seeds
        // All cells have LayerID = 0 as this is not a multilayer problem
        Kokkos::deep_copy(cell_type, Liquid);
        Kokkos::deep_copy(layer_id_all_layers, 0);

        // Local copies for lambda capture.
        auto cell_type_local = cell_type;
        auto grain_id_all_layers_local = grain_id_all_layers;

        // Determine which grains/active cells belong to which MPI ranks
        auto policy = Kokkos::RangePolicy<execution_space>(0, substrate_act_cells);
        Kokkos::parallel_for(
            "ConstrainedGrainInit", policy, KOKKOS_LAMBDA(const int &n) {
                // What are the X and Y coordinates of this active cell relative to the X and Y bounds of this rank?
                if ((act_cell_data(n, 1) >= grid.y_offset) && (act_cell_data(n, 1) < grid.y_offset + grid.ny_local)) {
                    // Convert X and Y coordinates to values relative to this MPI rank's grid (Z = 0 for these active
                    // cells, at bottom surface) GrainIDs come from the position on the list of substrate active cells
                    // to avoid reusing the same value
                    int coord_x = act_cell_data(n, 0);
                    int coord_y = act_cell_data(n, 1) - grid.y_offset;
                    int coord_z = 0;
                    int index = grid.get1DIndex(coord_x, coord_y, coord_z);
                    cell_type_local(index) = FutureActive;
                    grain_id_all_layers_local(index) = act_cell_data(n, 2); // assign GrainID > 0 to epitaxial seeds
                }
            });
        // Option to fill empty sites at bottom surface with the grain ID of the nearest grain
        if (_inputs.fill_bottom_surface) {
            auto md_policy =
                Kokkos::MDRangePolicy<execution_space, Kokkos::Rank<2, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>(
                    {0, 0}, {grid.nx, grid.ny_local});
            // For cells that are not associated with grain centers, optionally assign them the GrainID of the nearest
            // grain center
            Kokkos::parallel_for(
                "BaseplateGen", md_policy, KOKKOS_LAMBDA(const int coord_x, const int coord_y) {
                    int index_all_layers = grid.get1DIndex(coord_x, coord_y, 0);
                    if (grain_id_all_layers_local(index_all_layers) == 0) {
                        // This cell needs to be assigned a GrainID value
                        // Check each possible baseplate grain center to find the closest one
                        float min_distance_to_this_grain = grid.nx * grid.ny;
                        int min_distance_to_this_grain_grain_id = 0;
                        for (int n = 0; n < substrate_act_cells; n++) {
                            // Substrate grain center at coord_x_grain, coord_y_grain - how far is the cell at i,
                            // j+y_offset?
                            int coord_y_grain_global = act_cell_data(n, 1);
                            int coord_x_grain = act_cell_data(n, 0);
                            int coord_y_global = coord_y + grid.y_offset;
                            float distance_to_this_grain_x = coord_x - coord_x_grain;
                            float distance_to_this_grain_y = coord_y_global - coord_y_grain_global;
                            float distance_to_this_grain =
                                Kokkos::sqrt(distance_to_this_grain_x * distance_to_this_grain_x +
                                             distance_to_this_grain_y * distance_to_this_grain_y);
                            if (distance_to_this_grain < min_distance_to_this_grain) {
                                // This is the closest grain center to cell at "CAGridLocation" - update values
                                min_distance_to_this_grain = distance_to_this_grain;
                                min_distance_to_this_grain_grain_id = act_cell_data(n, 2);
                            }
                        }
                        // GrainID associated with the closest baseplate grain center
                        grain_id_all_layers_local(index_all_layers) = min_distance_to_this_grain_grain_id;
                        // These cells are also future active cells
                        // For directional solidification, only one "layer" so index for all layers can be used here to
                        // index cell_type
                        cell_type_local(index_all_layers) = FutureActive;
                    }
                });
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (id == 0) {
            std::cout << "Number of substrate active cells across all ranks: " << substrate_act_cells << std::endl;
            std::cout << "Grain struct initialized" << std::endl;
        }
    }

    // Initializes substrate grain structure using either a Voronoi assignment of grain ID values to cells or reading
    // the data from a file
    void initSubstrate(const int id, const Grid &grid, const unsigned long rng_seed,
                       view_type_int number_of_solidification_events) {

        // Determine the number of cells in the Z direction that are part of the baseplate
        int baseplate_size_z = getBaseplateSizeZ(id, grid);
        // Generate the baseplate microstructure, or read it from a file, to initialize the grain ID values from Z = 0
        // up to but not including Z = BaseplateTopZ
        if (_inputs.use_substrate_file)
            initBaseplateGrainID(id, grid, baseplate_size_z);
        else
            initBaseplateGrainID(id, grid, rng_seed, baseplate_size_z);

        // Powder layer extends from Z = powder_bottom_z up to but not including Z = powder_top_z
        // Bottom of layer is the next coordinate up from the baseplate
        int powder_bottom_z = Kokkos::round((_inputs.baseplate_top_z - grid.z_min) / grid.deltax) + 1;
        int powder_top_z = Kokkos::round((grid.z_max_layer[0] - grid.z_min) / grid.deltax) + 1;
        // Generate powder grain structure grain IDs for top of layer 0 if needed (i.e, if the powder layer height is
        // more than zero cells)
        if (powder_top_z > powder_bottom_z)
            initPowderGrainID(0, id, rng_seed, grid, powder_bottom_z, powder_top_z);

        // LayerID starts at -1 for all cells
        Kokkos::deep_copy(layer_id_all_layers, -1);

        // Initialize cell types and layer IDs based on whether cells will solidify in layer 0 or not
        initCellTypeLayerID(0, id, grid, number_of_solidification_events);

        MPI_Barrier(MPI_COMM_WORLD);
        if (id == 0) {
            std::cout << "Grain struct initialized" << std::endl;
        }
    }

    // Determine the height of the baseplate, in CA cells
    // The baseplate always starts at the simulation bottom (Z coordinate corresponding to z_min, Z index = 0),
    // regardless of whether the first layer melts the cells at the bottom or not. If baseplate_through_powder is true,
    // the baseplate microstructure extends through the entire simulation domain in Z (size nz). If
    // baseplate_through_powder is false, the baseplate top from the input file is used, or it is assumed that the top
    // of the baseplate is at Z = 0 microns
    int getBaseplateSizeZ(const int id, const Grid &grid) {
        int baseplate_size_z;
        if (_inputs.baseplate_through_powder)
            baseplate_size_z = grid.nz;
        else {
            baseplate_size_z = Kokkos::round((_inputs.baseplate_top_z - grid.z_min) / grid.deltax) + 1;
            int max_baseplate_size_z = Kokkos::round((grid.z_max_layer[0] - grid.z_min) / grid.deltax) + 1;
            if (baseplate_size_z > max_baseplate_size_z) {
                baseplate_size_z = max_baseplate_size_z;
                if (id == 0)
                    std::cout << "Warning: The specified location of the baseplate top is located above the "
                                 "temperature data for the first layer of the problem; it will be truncated to stop at "
                                 "the top Z coordinate of the first layer's temperature data"
                              << std::endl;
            }
        }
        if ((id == 0) && (baseplate_size_z < 1))
            std::cout << "Warning: no baseplate microstructure will be used, as the temperature data does not overlap "
                         "with the input value for the baseplate top"
                      << std::endl;
        return baseplate_size_z;
    }

    // Check that the substrate bounds span the simulation domain in the requested dimension, within some tolerance
    void checkSubstrateBound(const double substrate_bound_low, const int substrate_num_points,
                             const double substrate_cell_spacing, const double simulation_bound_low,
                             const double simulation_bound_high, std::string dim, const double tolerance) {
        if (simulation_bound_low < substrate_bound_low - tolerance) {
            std::string error = "Error: lower bound of simulation in " + dim + " ( " +
                                std::to_string(simulation_bound_low) +
                                ") extends beyond the lower bound of the substrate grain structure from the file (" +
                                std::to_string(substrate_bound_low) + ")";
            throw std::runtime_error(error);
        }
        const double substrate_bound_high =
            substrate_bound_low + static_cast<double>(substrate_num_points - 1) * substrate_cell_spacing;
        if (simulation_bound_high > substrate_bound_high + tolerance) {
            std::string error = "Error: upper bound of simulation in " + dim + " ( " +
                                std::to_string(simulation_bound_high) +
                                ") extends beyond the upper bound of the substrate grain structure from the file (" +
                                std::to_string(substrate_bound_high) + ")";
            throw std::runtime_error(error);
        }
    }

    // Initializes Grain ID values where the substrate comes from a file
    void initBaseplateGrainID(const int id, const Grid &grid, const int baseplate_size_z) {

        // Parse substrate input file
        std::ifstream substrate;
        substrate.open(_inputs.substrate_filename);
        // Ignore first two header lines
        skipLines(substrate, 2);

        // Is this data binary or ASCII
        std::string read_line;
        getline(substrate, read_line);
        const bool binary_s = (read_line.find("BINARY") != std::string::npos);

        // Ignore line
        skipLines(substrate, 1);

        // Get nx_s, ny_s, and nz_s
        std::vector<std::string> dims_read = splitString(substrate);
        const int nx_s = getInputInt(dims_read[1]);
        const int ny_s = getInputInt(dims_read[2]);
        const int nz_s = getInputInt(dims_read[3]);

        // Get origin location
        std::vector<std::string> org_read = splitString(substrate);
        const double x_min_s = getInputDouble(org_read[1]);
        const double y_min_s = getInputDouble(org_read[2]);
        const double z_min_s = getInputDouble(org_read[3]);

        // Get voxel spacing
        std::vector<std::string> vox_spacing_read = splitString(substrate);
        const double deltax_s = getInputDouble(vox_spacing_read[1]);
        if ((deltax_s != getInputDouble(vox_spacing_read[2])) || (deltax_s != getInputDouble(vox_spacing_read[3])))
            throw std::runtime_error("Error: substrate data must have same spacing in all directions");

        // Ensure substrate dimensions can sufficiently cover the solidification domain
        if (id == 0)
            std::cout << "Substrate dimensions from file are " << nx_s << " by " << ny_s << " by " << nz_s
                      << ", voxel spacing is " << deltax_s << std::endl;
        checkSubstrateBound(x_min_s, nx_s, deltax_s, grid.x_min, grid.x_max, "X", grid.deltax);
        checkSubstrateBound(y_min_s, ny_s, deltax_s, grid.y_min, grid.y_max, "Y", grid.deltax);
        checkSubstrateBound(z_min_s, nz_s, deltax_s, grid.z_min, _inputs.baseplate_top_z, "Z", grid.deltax);

        // Ignore line
        skipLines(substrate, 1);

        // Ensure data is of type integer
        getline(substrate, read_line);
        if (!read_line.find("int"))
            throw std::runtime_error("Error: substrate grain ID data must be type int");

        // Ignore line
        skipLines(substrate, 1);

        // Read grain ID data from file into the host view
        Kokkos::View<int ***, Kokkos::HostSpace> grain_id_s_host(Kokkos::ViewAllocateWithoutInitializing("GrainID_S"),
                                                                 nz_s, nx_s, ny_s);
        if (binary_s)
            grain_id_s_host =
                readBinaryField<Kokkos::View<int ***, Kokkos::HostSpace>, int>(substrate, nx_s, ny_s, nz_s, "GrainID");
        else
            grain_id_s_host =
                readASCIIField<Kokkos::View<int ***, Kokkos::HostSpace>>(substrate, nx_s, ny_s, nz_s, "GrainID");
        if (id == 0)
            std::cout << "Successfully read substrate GrainID data from the file" << std::endl;
        substrate.close();
        // Copy host view data to device
        auto grain_id_s = Kokkos::create_mirror_view_and_copy(memory_space(), grain_id_s_host);

        // Assign each CA cell the GrainID of the nearest voxel in the substrate, returning the max GrainID
        const int baseplate_top_coord_1D = grid.nx * grid.ny_local * baseplate_size_z;
        int max_grain_id_local = 0;
        auto grain_id_all_layers_local = grain_id_all_layers;
        auto policy = Kokkos::RangePolicy<execution_space>(0, baseplate_top_coord_1D);
        Kokkos::parallel_reduce(
            "BaseplateInit", policy,
            KOKKOS_LAMBDA(const int &index_all_layers, int &update) {
                // x, y, z associated with this 1D cell location, with respect to the global simulation bounds
                const int coord_x_global = grid.getCoordX(index_all_layers);
                const int coord_y_global = grid.getCoordY(index_all_layers) + grid.y_offset;
                const int coord_z_global = grid.getCoordZ(index_all_layers);
                const double x_location = grid.x_min + coord_x_global * grid.deltax;
                const double y_location = grid.y_min + coord_y_global * grid.deltax;
                const double z_location = grid.z_min + coord_z_global * grid.deltax;
                // What voxel does this correspond to in the substrate?
                const int coord_x_s = Kokkos::round((x_location - x_min_s) / deltax_s);
                const int coord_y_s = Kokkos::round((y_location - y_min_s) / deltax_s);
                const int coord_z_s = Kokkos::round((z_location - z_min_s) / deltax_s);
                grain_id_all_layers_local(index_all_layers) = grain_id_s(coord_z_s, coord_x_s, coord_y_s);
                if (grain_id_all_layers_local(index_all_layers) > update)
                    update = grain_id_all_layers_local(index_all_layers);
            },
            Kokkos::Max<int>(max_grain_id_local));
        Kokkos::fence();
        // Avoid reusing GrainID in next layer's powder grain structure
        // TODO: Also account for negative grain ids in the baseplate, in the case where this comes from previous ExaCA
        // data during future restart file reads
        max_grain_id_local++;
        MPI_Allreduce(&max_grain_id_local, &next_layer_first_epitaxial_grain_id, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Substrate file read complete" << std::endl;
    }

    // Initializes Grain ID values where the baseplate is generated using an input grain spacing and a Voronoi
    // Tessellation
    void initBaseplateGrainID(const int id, const Grid &grid, const unsigned long rng_seed,
                              const int baseplate_size_z) {

        std::mt19937_64 gen(rng_seed);

        // Based on the baseplate volume (convert to cubic microns to match units) and the substrate grain spacing,
        // determine the number of baseplate grains
        int baseplate_volume = grid.nx * grid.ny * baseplate_size_z;
        double baseplate_volume_microns = baseplate_volume * Kokkos::pow(grid.deltax, 3) * Kokkos::pow(10, 18);
        double substrate_mean_grain_volume_microns = Kokkos::pow(_inputs.substrate_grain_spacing, 3);
        int number_of_baseplate_grains = Kokkos::round(baseplate_volume_microns / substrate_mean_grain_volume_microns);
        // Need at least 1 baseplate grain, cannot have more baseplate grains than cells in the baseplate
        number_of_baseplate_grains = std::max(number_of_baseplate_grains, 1);
        number_of_baseplate_grains = std::min(number_of_baseplate_grains, grid.nx * grid.ny * baseplate_size_z);
        // TODO: Use device RNG to generate baseplate grain locations, instead of host with copy
        // List of potential grain IDs (starting at 1) - index corresponds to the associated CA cell location
        // Assign positive values for indices 0 through number_of_baseplate_grains-1, assign zeros to the remaining
        // indices
        std::vector<int> baseplate_grain_locations(baseplate_volume);
        std::vector<int> baseplate_grain_ids(baseplate_volume);
        for (int i = 0; i < baseplate_volume; i++) {
            baseplate_grain_locations[i] = i;
            if (i < number_of_baseplate_grains)
                baseplate_grain_ids[i] = i + 1;
            else
                baseplate_grain_ids[i] = 0;
        }
        // Shuffle list of grain IDs
        std::shuffle(baseplate_grain_ids.begin(), baseplate_grain_ids.end(), gen);

        // Create views of baseplate grain IDs and locations - copying baseplate_grain_ids and baseplate_grain_locations
        // values only for indices where baseplate_grain_ids(i) != 0 (cells with baseplate_grain_ids(i) = 0 were not
        // assigned baseplate center locations)
        view_type_int_host baseplate_grain_locations_host(
            Kokkos::ViewAllocateWithoutInitializing("baseplate_grain_locations_Host"), number_of_baseplate_grains);
        view_type_int_host baseplate_grain_ids_host(Kokkos::ViewAllocateWithoutInitializing("baseplate_grain_ids_Host"),
                                                    number_of_baseplate_grains);
        int count = 0;
        for (int i = 0; i < baseplate_volume; i++) {
            if (baseplate_grain_ids[i] != 0) {
                baseplate_grain_locations_host(count) = baseplate_grain_locations[i];
                baseplate_grain_ids_host(count) = baseplate_grain_ids[i];
                count++;
            }
        }

        // Copy baseplate views to the device
        auto baseplate_grain_ids_device = Kokkos::create_mirror_view_and_copy(memory_space(), baseplate_grain_ids_host);
        auto baseplate_grain_locations_device =
            Kokkos::create_mirror_view_and_copy(memory_space(), baseplate_grain_locations_host);
        if (id == 0) {
            std::cout << "Baseplate spanning domain coordinates Z = 0 through " << baseplate_size_z - 1 << std::endl;
            std::cout << "Number of baseplate grains: " << number_of_baseplate_grains << std::endl;
        }

        // Local copies for lambda capture.
        auto grain_id_all_layers_local = grain_id_all_layers;

        // First, assign cells that are associated with grain centers the appropriate non-zero GrainID values (assumes
        // GrainID values were initialized to zeros)
        auto policy = Kokkos::RangePolicy<execution_space>(0, number_of_baseplate_grains);
        Kokkos::parallel_for(
            "BaseplateInit", policy, KOKKOS_LAMBDA(const int &n) {
                int baseplate_grain_loc = baseplate_grain_locations_device(n);
                // x, y, z associated with baseplate grain "n", at 1D coordinate "BaseplateGrainLoc"
                int coord_z_all_layers = grid.getCoordZGlobal(baseplate_grain_loc);
                int coord_y_global = grid.getCoordYGlobal(baseplate_grain_loc);
                int coord_x = grid.getCoordXGlobal(baseplate_grain_loc);
                if ((coord_y_global >= grid.y_offset) && (coord_y_global < grid.y_offset + grid.ny_local)) {
                    // This grain is associated with a cell on this MPI rank
                    int coord_y = coord_y_global - grid.y_offset;
                    int index_all_layers = grid.get1DIndex(coord_x, coord_y, coord_z_all_layers);
                    grain_id_all_layers_local(index_all_layers) = baseplate_grain_ids_device(n);
                }
            });
        Kokkos::fence();

        auto md_policy =
            Kokkos::MDRangePolicy<execution_space, Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>(
                {0, 0, 0}, {baseplate_size_z, grid.nx, grid.ny_local});

        // For cells that are not associated with grain centers, assign them the GrainID of the nearest grain center
        Kokkos::parallel_for(
            "BaseplateGen", md_policy,
            KOKKOS_LAMBDA(const int coord_z_all_layers, const int coord_x, const int coord_y) {
                int index_all_layers = grid.get1DIndex(coord_x, coord_y, coord_z_all_layers);
                if (grain_id_all_layers_local(index_all_layers) == 0) {
                    // This cell needs to be assigned a GrainID value
                    // Check each possible baseplate grain center to find the closest one
                    float min_distance_to_this_grain = grid.nx * grid.ny * baseplate_size_z;
                    int min_distance_to_this_grain_grain_id = 0;
                    for (int n = 0; n < number_of_baseplate_grains; n++) {
                        // Baseplate grain center at x_n, y_n, z_n - how far is the cell at i, j+y_offset, k?
                        int coord_z_grain_all_layers = grid.getCoordZGlobal(baseplate_grain_locations_device(n));
                        int coord_y_grain_global = grid.getCoordYGlobal(baseplate_grain_locations_device(n));
                        int coord_x_grain = grid.getCoordXGlobal(baseplate_grain_locations_device(n));
                        int coord_y_global = coord_y + grid.y_offset;
                        float distance_to_this_grain_x = coord_x - coord_x_grain;
                        float distance_to_this_grain_y = coord_y_global - coord_y_grain_global;
                        float distance_to_this_grain_z = coord_z_grain_all_layers - coord_z_all_layers;
                        float distance_to_this_grain =
                            Kokkos::hypot(distance_to_this_grain_x, distance_to_this_grain_y, distance_to_this_grain_z);
                        if (distance_to_this_grain < min_distance_to_this_grain) {
                            // This is the closest grain center to cell at "CAGridLocation" - update values
                            min_distance_to_this_grain = distance_to_this_grain;
                            min_distance_to_this_grain_grain_id = baseplate_grain_ids_device(n);
                        }
                    }
                    // GrainID associated with the closest baseplate grain center
                    grain_id_all_layers_local(index_all_layers) = min_distance_to_this_grain_grain_id;
                }
            });

        next_layer_first_epitaxial_grain_id =
            number_of_baseplate_grains + 1; // avoid reusing GrainID in next layer's powder grain structure
        if (id == 0)
            std::cout << "Baseplate grain structure initialized" << std::endl;
    }

    // Each layer's top Z coordinates are seeded with CA-cell sized substrate grains (emulating bulk nucleation
    // alongside the edges of partially melted powder particles). These Z coordinates span powder_bottom_z up to but not
    // including powder_top_z
    void initPowderGrainID(const int layernumber, const int id, const unsigned long rng_seed, const Grid &grid,
                           const int powder_bottom_z, const int powder_top_z) {

        // On all ranks, generate list of powder grain IDs (starting with next_layer_first_epitaxial_grain_id, and
        // shuffle them so that their locations aren't sequential and depend on the rng_seed (different for each layer)
        std::mt19937_64 gen(rng_seed + static_cast<unsigned long>(layernumber));
        std::uniform_real_distribution<double> dis(0.0, 1.0);

        // TODO: This should be performed on the device, rather than the host
        int powder_layer_height = powder_top_z - powder_bottom_z;
        int powder_layer_cells = grid.nx * grid.ny * powder_layer_height;
        int powder_layer_assigned_cells =
            Kokkos::round(static_cast<double>(powder_layer_cells) * _inputs.powder_active_fraction);
        MPI_Barrier(MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Initializing powder layer for Z = " << powder_bottom_z << " through " << powder_top_z - 1
                      << " (" << grid.nx * grid.ny * powder_layer_height << " cells): powder layer has "
                      << powder_layer_assigned_cells << " cells assigned new grain ID values" << std::endl;

        // Associate powder grain IDs with CA cells in the powder layer, if nonzero number of powder cells
        if (powder_layer_assigned_cells > 0) {

            std::vector<int> powder_grain_ids(powder_layer_cells, 0);
            for (int n = 0; n < powder_layer_assigned_cells; n++) {
                powder_grain_ids[n] = n + next_layer_first_epitaxial_grain_id; // assigned a nonzero GrainID
            }
            std::shuffle(powder_grain_ids.begin(), powder_grain_ids.end(), gen);
            // Wrap powder layer GrainIDs into an unmanaged view, then copy to the device
            view_type_int_unmanaged powder_grain_ids_host(powder_grain_ids.data(), powder_layer_cells);
            auto powder_grain_ids_device = Kokkos::create_mirror_view_and_copy(memory_space(), powder_grain_ids_host);

            int powder_start = grid.nx * grid.ny * powder_bottom_z;
            if (id == 0)
                std::cout << "Powder layer grain ID values range from " << next_layer_first_epitaxial_grain_id
                          << " through " << next_layer_first_epitaxial_grain_id + powder_layer_assigned_cells - 1
                          << std::endl;

            // Iterate over all cells in the powder layer, on each rank loading the powder grain ID data for local cell
            // locations
            auto grain_id_all_layers_local = grain_id_all_layers;
            auto powder_policy =
                Kokkos::MDRangePolicy<execution_space, Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>(
                    {powder_bottom_z, 0, 0}, {powder_top_z, grid.nx, grid.ny});
            Kokkos::parallel_for(
                "PowderGrainInit", powder_policy,
                KOKKOS_LAMBDA(const int coord_z_all_layers, const int coord_x, const int coord_y_global) {
                    // Is this powder coordinate in X and Y in bounds for this rank? Is the grain id of this site
                    // unassigned (wasn't captured during solidification of the previous layer)?
                    if ((coord_y_global >= grid.y_offset) && (coord_y_global < grid.y_offset + grid.ny_local)) {
                        int coord_y = coord_y_global - grid.y_offset;
                        int index_all_layers = grid.get1DIndex(coord_x, coord_y, coord_z_all_layers);
                        if (grain_id_all_layers_local(index_all_layers) == 0) {
                            int index_all_ranks_all_layers =
                                coord_z_all_layers * grid.nx * grid.ny + coord_x * grid.ny + coord_y_global;
                            grain_id_all_layers_local(index_all_layers) =
                                powder_grain_ids_device(index_all_ranks_all_layers - powder_start);
                        }
                    }
                });
            Kokkos::fence();

            // Update next_layer_first_epitaxial_grain_id for next layer
            next_layer_first_epitaxial_grain_id += powder_layer_assigned_cells;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Initialized powder grain structure for layer " << layernumber << std::endl;
    }

    // Sets up views, powder layer (if necessary), and cell types for the next layer of a multilayer problem
    void initNextLayer(const int nextlayernumber, const int id, const Grid &grid, const unsigned long rng_seed,
                       view_type_int number_of_solidification_events) {

        // Subviews for the next layer's grain id, layer id, cell type are constructed based on updated layer bound
        // z_layer_bottom
        // Powder layer extends from Z = powder_bottom_z (1 cell above the top of the previous layer) up to but not
        // including Z = powder_top_z
        int powder_bottom_z = Kokkos::round((grid.z_max_layer[nextlayernumber - 1] - grid.z_min) / grid.deltax) + 1;
        int powder_top_z = Kokkos::round((grid.z_max_layer[nextlayernumber] - grid.z_min) / grid.deltax) + 1;
        if (!(_inputs.baseplate_through_powder))
            initPowderGrainID(nextlayernumber, id, rng_seed, grid, powder_bottom_z, powder_top_z);

        // Initialize active cell data structures and nuclei locations for the next layer "layernumber + 1"
        initCellTypeLayerID(nextlayernumber, id, grid, number_of_solidification_events);
        // Current layer melt pool edge indicator, if needed
        if (_store_melt_pool_edge)
            getCurrentLayerMeltEdge(grid.layer_range);
    }

    // Initializes cells for the current layer as either solid (don't resolidify) or tempsolid (will melt and
    // resolidify)
    void initCellTypeLayerID(const int layernumber, const int id, const Grid &grid,
                             view_type_int number_of_solidification_events) {

        int melt_pool_cell_count;
        // Realloc celltype to the domain size of the next layer
        Kokkos::realloc(cell_type, grid.domain_size);
        // Local copies for lambda capture.
        auto cell_type_local = cell_type;
        auto layer_id_all_layers_local = layer_id_all_layers;

        auto policy = Kokkos::RangePolicy<execution_space>(0, grid.domain_size);

        Kokkos::parallel_reduce(
            "cell_typeInitSolidRM", policy,
            KOKKOS_LAMBDA(const int &index, int &local_count) {
                int index_all_layers = index + grid.z_layer_bottom * grid.nx * grid.ny_local;
                if (number_of_solidification_events(index) > 0) {
                    cell_type_local(index) = TempSolid;
                    layer_id_all_layers_local(index_all_layers) = layernumber;
                    local_count++;
                }
                else
                    cell_type_local(index) = Solid;
            },
            melt_pool_cell_count);
        int total_melt_pool_cell_count;
        MPI_Reduce(&melt_pool_cell_count, &total_melt_pool_cell_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Number of cells across all ranks to undergo solidification at least once: "
                      << total_melt_pool_cell_count << std::endl;
    }

    // Stores/returns the volume fraction of nucleated grains to the console
    // Moved from CAfunctions.hpp
    float calcVolFractionNucleated(const int id, const Grid &grid) {

        // For interior cells, add the number of cells that underwent melting/solidification and the number of cells
        // with sub-zero grain IDs
        int melted_cells_local = 0;
        int nucleated_grain_cells_local = 0;
        // Local copies for lambda capture.
        auto grain_id_all_layers_local = grain_id_all_layers;
        auto layer_id_all_layers_local = layer_id_all_layers;
        Kokkos::parallel_reduce(
            "NumSolidifiedCells", grid.domain_size,
            KOKKOS_LAMBDA(const int index, int &update_meltcount, int &update_nucleatecount) {
                int coord_y = grid.getCoordY(index);
                // Is this Y coordinate in the halo region? If so, do not increment counter
                bool in_halo_region = false;
                if (((coord_y == 0) && (!grid.at_south_boundary)) ||
                    ((coord_y == grid.ny_local - 1) && (!grid.at_north_boundary)))
                    in_halo_region = true;
                if ((grain_id_all_layers_local(index) < 0) && (!in_halo_region))
                    update_nucleatecount++;
                if ((layer_id_all_layers_local(index) != -1) && (!in_halo_region))
                    update_meltcount++;
            },
            melted_cells_local, nucleated_grain_cells_local);

        // Reduce the values by summing over all ranks
        int melted_cells_global, nucleated_grain_cells_global;
        MPI_Allreduce(&melted_cells_local, &melted_cells_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&nucleated_grain_cells_local, &nucleated_grain_cells_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        // Calculate nucleated grain fraction
        float vol_fraction_nucleated =
            static_cast<float>(nucleated_grain_cells_global) / static_cast<float>(melted_cells_global);
        if (id == 0)
            std::cout << "The fraction of the solidified material consisting of nucleated grains is "
                      << vol_fraction_nucleated << std::endl;
        return vol_fraction_nucleated;
    }

    // If storing whether or not a cell is at a melt pool edge, update the value
    KOKKOS_INLINE_FUNCTION
    void setMeltEdge(const int index, const bool updated_state) const {
        if (_store_melt_pool_edge)
            melt_edge(index) = updated_state;
    }

    // Take a view consisting of data for all layers, and return a subview of the same type consisting of just the cells
    // corresponding to the current layer of a multilayer problem
    auto getGrainIDSubview(const Grid &grid) const { return Kokkos::subview(grain_id_all_layers, grid.layer_range); }
    auto getLayerIDSubview(const Grid &grid) const { return Kokkos::subview(layer_id_all_layers, grid.layer_range); }
};

#endif
