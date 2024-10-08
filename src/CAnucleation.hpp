// Copyright Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_NUCLEATION_HPP
#define EXACA_NUCLEATION_HPP

#include "CAcelldata.hpp"
#include "CAgrid.hpp"
#include "CAinputs.hpp"
#include "CAinterface.hpp"
#include "CAtemperature.hpp"

#include "mpi.h"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <random>
#include <string>
#include <vector>

// Data regarding nucleation events in the domain
template <typename MemorySpace>
struct Nucleation {

    using memory_space = MemorySpace;
    using view_type_int = Kokkos::View<int *, memory_space>;
    using view_type_int_host = typename view_type_int::HostMirror;

    // Using the default exec space for this memory space.
    using execution_space = typename memory_space::execution_space;

    // Four counters tracked here:
    // 1. nuclei_whole_domain - tracks all nuclei (regardless of whether an event would be possible based on the layer
    // ID and cell type), used for Grain ID assignment to ensure that no Grain ID get reused - same on all MPI ranks
    // 2. possible_nuclei - the subset of Nuclei_ThisLayer that are located within the bounds of a
    // given MPI rank that may possibly occur (nuclei locations are associated with a liquid cell with a layer ID that
    // matches this layer number). Starts at 0 each layer
    // 3. nucleation_counter - the number of nucleation events that have actually either failed or succeeded on a given
    // MPI rank. Starts at 0 each layer
    // 4. successful_nucleation_counter - the number of nucleation events that have successfully occurred at a given
    // point in the simulation of a layer. Starts at 0 each layer.
    int nuclei_whole_domain, possible_nuclei, nucleation_counter, successful_nucleation_counter;

    // The time steps at which nucleation events will occur in the given layer, on the host
    view_type_int_host nucleation_times_host;
    // The locations and grain IDs of the potential nuclei in the layer
    view_type_int nuclei_locations, nuclei_grain_id;
    // Nucleation inputs from file
    NucleationInputs _inputs;

    // Constructor - initialize CA views using input for initial guess at number of possible events
    // Default is that no nucleation has occurred to this point, optional input argument to start with the counter at a
    // specific value
    Nucleation(const int estimated_nuclei_this_rank_this_layer, const NucleationInputs inputs,
               const int num_prior_nuclei = 0)
        : nucleation_times_host(view_type_int_host(Kokkos::ViewAllocateWithoutInitializing("NucleationTimes_Host"),
                                                   estimated_nuclei_this_rank_this_layer))
        , nuclei_locations(view_type_int(Kokkos::ViewAllocateWithoutInitializing("NucleiLocations"),
                                         estimated_nuclei_this_rank_this_layer))
        , nuclei_grain_id(view_type_int(Kokkos::ViewAllocateWithoutInitializing("NucleiGrainID"),
                                        estimated_nuclei_this_rank_this_layer))
        , _inputs(inputs) {

        nuclei_whole_domain = num_prior_nuclei;
        resetNucleiCounters(); // start counters at 0
    }

    // Reset the nuclei counters prior to nuclei initialization for a layer
    void resetNucleiCounters() {
        // Init appropriate counters for the layer to 0 - possible nuclei for this layer will be calculated in this
        // function
        possible_nuclei = 0;
        nucleation_counter = 0;
        successful_nucleation_counter = 0;
    }

    // Initialize nucleation site locations, GrainID values, and time at which nucleation events will potentially occur,
    // accounting for multiple possible nucleation events in cells that melt and solidify multiple times
    template <class... Params>
    void placeNuclei(const Temperature<memory_space> &temperature, const unsigned long rng_seed, const int layernumber,
                     const Grid &grid, const int id) {

        // TODO: convert this subroutine into kokkos kernels, rather than copying data back to the host, and nucleation
        // data back to the device again. This is currently performed on the device due to heavy usage of standard
        // library algorithm functions Copy temperature data into temporary host views for this subroutine
        auto max_solidification_events_host =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), temperature.max_solidification_events);
        auto number_of_solidification_events_host =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), temperature.number_of_solidification_events);
        auto liquidus_time_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), temperature.liquidus_time);
        auto cooling_rate_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), temperature.cooling_rate);

        // Use new RNG seed for each layer
        std::mt19937_64 generator(rng_seed + static_cast<unsigned long>(layernumber));
        // Uniform distributions for nuclei location assignment - associate each nucleation event with an XYZ coordinate
        // in meters
        std::uniform_real_distribution<double> nucleation_site_dist_x(grid.x_min, grid.x_max);
        std::uniform_real_distribution<double> nucleation_site_dist_y(grid.y_min, grid.y_max);
        std::uniform_real_distribution<double> nucleation_site_dist_z(grid.z_min_layer(layernumber),
                                                                      grid.z_max_layer(layernumber));
        // Gaussian distribution of nucleation undercooling
        std::normal_distribution<double> nucleation_undercooling_dist(_inputs.dtn, _inputs.dtsigma);

        // Max number of nucleated grains in this layer
        // Use long int in intermediate steps calculating the number of nucleated grains, though the number should be
        // small enough to be stored as an int
        const double domain_volume = (grid.x_max - grid.x_min) * (grid.y_max - grid.y_min) *
                                     (grid.z_max_layer(layernumber) - grid.z_min_layer(layernumber));
        // If each cell underwent solidification 1x, the number of potential nuclei in the layer
        const long int nuclei_this_layer_single_long = std::lround(_inputs.n_max * domain_volume);
        // Multiplier for the number of nucleation events per layer, based on the max number of solidification events
        long int nuclei_multiplier_long = static_cast<long int>(max_solidification_events_host(layernumber));
        long int nuclei_this_layer_long = nuclei_this_layer_single_long * nuclei_multiplier_long;
        // Vector and view sizes should be type int for indexing and grain ID assignment purposes -
        // Nuclei_ThisLayer_long should be less than INT_MAX
        if (nuclei_this_layer_long > INT_MAX)
            throw std::runtime_error("Error: Number of potential nucleation sites in the system exceeds the number of "
                                     "valid GrainID; either nucleation density, the number of melt-solidification "
                                     "events in the temperature data, or the domain size should be reduced");
        int nuclei_this_layer_single = static_cast<int>(nuclei_this_layer_single_long);
        int nuclei_this_layer = static_cast<int>(nuclei_this_layer_long);
        int nuclei_multiplier = static_cast<int>(nuclei_multiplier_long);

        // Nuclei Grain ID are assigned to avoid reusing values from previous layers
        std::vector<int> nuclei_grain_id_whole_domain_v(nuclei_this_layer);
        std::vector<double> nuclei_undercooling_whole_domain_v(nuclei_this_layer);
        // Views for storing potential nucleated grain coordinates
        view_type_int_host nuclei_x(Kokkos::ViewAllocateWithoutInitializing("NucleiX"), nuclei_this_layer);
        view_type_int_host nuclei_y(Kokkos::ViewAllocateWithoutInitializing("NucleiY"), nuclei_this_layer);
        view_type_int_host nuclei_z(Kokkos::ViewAllocateWithoutInitializing("NucleiZ"), nuclei_this_layer);

        for (int meltevent = 0; meltevent < nuclei_multiplier; meltevent++) {
            for (int n = 0; n < nuclei_this_layer_single; n++) {
                int n_event = meltevent * nuclei_this_layer_single + n;
                // Generate possible nuclei locations
                double nuclei_x_unrounded = nucleation_site_dist_x(generator);
                double nuclei_y_unrounded = nucleation_site_dist_y(generator);
                double nuclei_z_unrounded = nucleation_site_dist_z(generator);
                // Associate these locations with a specific cell on the grid
                nuclei_x(n_event) = Kokkos::round((nuclei_x_unrounded - grid.x_min) / grid.deltax);
                nuclei_y(n_event) = Kokkos::round((nuclei_y_unrounded - grid.y_min) / grid.deltax);
                nuclei_z(n_event) = Kokkos::round((nuclei_z_unrounded - grid.z_min_layer[layernumber]) / grid.deltax);
                // Assign each nuclei a Grain ID (negative values used for nucleated grains) and an undercooling
                nuclei_grain_id_whole_domain_v[n_event] =
                    -(nuclei_whole_domain + n_event + 1); // avoid using grain ID 0
                nuclei_undercooling_whole_domain_v[n_event] = nucleation_undercooling_dist(generator);
            }
        }

        // Shuffle these vectors to make sure the same grain IDs and undercooling don't end up in the same spots each
        // layer
        std::shuffle(nuclei_grain_id_whole_domain_v.begin(), nuclei_grain_id_whole_domain_v.end(), generator);
        std::shuffle(nuclei_undercooling_whole_domain_v.begin(), nuclei_undercooling_whole_domain_v.end(), generator);
        if ((id == 0) && (nuclei_this_layer > 0))
            std::cout << "Range of Grain IDs from which layer " << layernumber << " nucleation events were selected: -"
                      << nuclei_whole_domain + 1 << " through -" << nuclei_whole_domain + nuclei_this_layer
                      << std::endl;
        // Update number of nuclei counter for whole domain based on the number of nuclei in this layer
        nuclei_whole_domain += nuclei_this_layer;

        // Loop through nuclei for this layer - each MPI rank storing the nucleation events that are possible (i.e,
        // nucleation event is associated with a CA cell on that MPI rank's subdomain, the cell is liquid type, and the
        // cell is associated with the current layer of the multilayer problem) Don't put nuclei in "ghost" cells -
        // those nucleation events occur on other ranks and the existing halo exchange functionality will handle this
        std::vector<int> nuclei_grain_id_myrank_v(nuclei_this_layer), nuclei_location_myrank_v(nuclei_this_layer);
        // Store nucleation times as doubles to correctly order events in times, then convert to time steps later
        std::vector<double> nucleation_times_myrank_v(nuclei_this_layer);
        for (int meltevent = 0; meltevent < nuclei_multiplier; meltevent++) {
            for (int n = 0; n < nuclei_this_layer_single; n++) {
                int n_event = meltevent * nuclei_this_layer_single + n;
                if (((nuclei_y(n_event) > grid.y_offset) || (grid.at_south_boundary)) &&
                    ((nuclei_y(n_event) < grid.y_offset + grid.ny_local - 1) || (grid.at_north_boundary))) {
                    // Convert 3D location (using global X and Y coordinates) into a 1D location (using local X and Y
                    // coordinates) for the possible nucleation event, both as relative to the bottom of this layer
                    int nuclei_location_this_layer =
                        grid.get1DIndex(nuclei_x(n_event), nuclei_y(n_event) - grid.y_offset, nuclei_z(n_event));
                    // Criteria for placing a nucleus - whether or not this nuclei is associated with a solidification
                    // event
                    if (meltevent < number_of_solidification_events_host(nuclei_location_this_layer)) {
                        // Nucleation event is possible - cell undergoes solidification at least once, this nucleation
                        // event is associated with one of the time periods during which the associated cell undergoes
                        // solidification
                        nuclei_location_myrank_v[possible_nuclei] = nuclei_location_this_layer;

                        double liq_time_this_event =
                            static_cast<double>(liquidus_time_host(nuclei_location_this_layer, meltevent, 1));
                        double cooling_rate_this_event =
                            static_cast<double>(cooling_rate_host(nuclei_location_this_layer, meltevent));
                        double time_to_nuc_und =
                            liq_time_this_event + nuclei_undercooling_whole_domain_v[n_event] / cooling_rate_this_event;
                        if (liq_time_this_event > time_to_nuc_und)
                            time_to_nuc_und = liq_time_this_event;
                        nucleation_times_myrank_v[possible_nuclei] = time_to_nuc_und;
                        // Assign this cell the potential nucleated grain ID
                        nuclei_grain_id_myrank_v[possible_nuclei] = nuclei_grain_id_whole_domain_v[n_event];
                        // Increment counter on this MPI rank
                        possible_nuclei++;
                    }
                }
            }
        }

        // How many nucleation events are actually possible (associated with a cell in this layer that will undergo
        // solidification)?
        int possible_nuclei_all_ranks_this_layer;
        MPI_Reduce(&possible_nuclei, &possible_nuclei_all_ranks_this_layer, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Number of potential nucleation events in layer " << layernumber << " : "
                      << possible_nuclei_all_ranks_this_layer << std::endl;

        // Now that the number of nucleation events on each rank is known, resize these vectors
        nuclei_location_myrank_v.resize(possible_nuclei);
        nuclei_grain_id_myrank_v.resize(possible_nuclei);
        nucleation_times_myrank_v.resize(possible_nuclei);

        // Sort the list of time steps at which nucleation occurs, keeping the time steps paired with the corresponding
        // locations for nucleation events and grain IDs
        std::vector<std::tuple<double, int, int>> nucleation_time_loc_id;
        nucleation_time_loc_id.reserve(possible_nuclei);
        for (int n = 0; n < possible_nuclei; n++) {
            nucleation_time_loc_id.push_back(std::make_tuple(nucleation_times_myrank_v[n], nuclei_location_myrank_v[n],
                                                             nuclei_grain_id_myrank_v[n]));
        }
        // Sorting from low to high
        std::sort(nucleation_time_loc_id.begin(), nucleation_time_loc_id.end());

        // With possible_nuclei for this rank and this layer now known, resize views appropriately
        Kokkos::resize(nucleation_times_host, possible_nuclei);
        Kokkos::resize(nuclei_locations, possible_nuclei);
        Kokkos::resize(nuclei_grain_id, possible_nuclei);

        // Create temporary view to store nucleation locations, grain ID data initialized on the host
        // nucleation_times_host are stored using a host view that is passed to the nucleateGrain subroutine and later
        // used
        view_type_int_host nuclei_locations_host(Kokkos::ViewAllocateWithoutInitializing("NucleiLocations_Host"),
                                                 possible_nuclei);
        view_type_int_host nuclei_grain_id_host(Kokkos::ViewAllocateWithoutInitializing("NucleiGrainID_Host"),
                                                possible_nuclei);
        for (int n = 0; n < possible_nuclei; n++) {
            double nucleation_time = std::get<0>(nucleation_time_loc_id[n]);
            // Convert nucleation time to a time step
            nucleation_times_host(n) = Kokkos::round(nucleation_time);
            nuclei_locations_host(n) = std::get<1>(nucleation_time_loc_id[n]);
            nuclei_grain_id_host(n) = std::get<2>(nucleation_time_loc_id[n]);
        }
        // Copy nucleation data to the device
        nuclei_locations = Kokkos::create_mirror_view_and_copy(memory_space(), nuclei_locations_host);
        nuclei_grain_id = Kokkos::create_mirror_view_and_copy(memory_space(), nuclei_grain_id_host);
        MPI_Barrier(MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Nuclei initialized" << std::endl;
    }

    // Check for nucleation events on this time step, updating the corresponding cell appropriately for any successful
    // nucleation event
    void nucleateGrain(const int cycle, const Grid &grid, CellData<memory_space> &celldata,
                       Interface<memory_space> interface) {

        auto grain_id = celldata.getGrainIDSubview(grid);

        // Is there nucleation left in this layer to check?
        if (nucleation_counter < possible_nuclei) {
            // Is there at least one potential nucleation event on this rank, at this time step?
            if (cycle == nucleation_times_host(nucleation_counter)) {
                bool nucleation_check = true;
                int first_event = nucleation_counter; // first potential nucleation event to check
                // Are there any other nucleation events this time step to check?
                while (nucleation_check) {
                    nucleation_counter++;
                    // If the previous nucleation event was the last one for this layer of the simulation, exit loop
                    if (nucleation_counter == possible_nuclei)
                        break;
                    // If the next nucleation event corresponds to a future time step, finish check
                    if (cycle != nucleation_times_host(nucleation_counter))
                        nucleation_check = false;
                }
                int last_event = nucleation_counter;
                // parallel_reduce checks each potential nucleation event this time step (first_event, up to but not
                // including last_event)
                int nucleation_this_dt = 0; // return number of successful event from parallel_reduce
                auto nuclei_locations_local = nuclei_locations;
                auto nuclei_grain_id_local = nuclei_grain_id;

                // Launch kokkos kernel - check if the corresponding CA cell location is liquid
                auto policy = Kokkos::RangePolicy<execution_space>(first_event, last_event);
                Kokkos::parallel_reduce(
                    "NucleiUpdateLoop", policy,
                    KOKKOS_LAMBDA(const int nucleation_counter_device, int &update) {
                        int nucleation_event_location = nuclei_locations_local(nucleation_counter_device);
                        int update_val =
                            FutureActive; // added to steering vector to become a new active cell as part of cellCapture
                        int old_val = Liquid;
                        int old_cell_type_value = Kokkos::atomic_compare_exchange(
                            &celldata.cell_type(nucleation_event_location), old_val, update_val);
                        if (old_cell_type_value == Liquid) {
                            // Successful nucleation event - atomic update of cell type, proceeded if the atomic
                            // exchange is successful (cell was liquid) Add future active cell location to steering
                            // vector and change cell type, assign new Grain ID
                            grain_id(nucleation_event_location) = nuclei_grain_id_local(nucleation_counter_device);
                            interface.steering_vector(Kokkos::atomic_fetch_add(&interface.num_steer(0), 1)) =
                                nucleation_event_location;
                            // This cell was not at the edge of the temperature field - set indicator to false if this
                            // is being tracked
                            celldata.setMeltEdge(nucleation_event_location, false);
                            update++;
                        }
                    },
                    nucleation_this_dt);
                // Update the number of successful nuclei counter with the number of successful nucleation events from
                // this time step (nucleation_this_dt)
                successful_nucleation_counter += nucleation_this_dt;
            }
        }
    }
};

#endif
