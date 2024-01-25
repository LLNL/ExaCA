// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "ExaCA.hpp"

#include "mpi.h"

#include <string>
#include <vector>

void runExaCA(int id, int np, std::string input_file) {

    // Run on the default space.
    using memory_space = Kokkos::DefaultExecutionSpace::memory_space;

    double nucl_time = 0.0, create_sv_time = 0.0, capture_time = 0.0, ghost_time = 0.0;
    double start_nucl_time, start_create_sv_time, start_capture_time, start_ghost_time;
    double start_init_time = MPI_Wtime();

    // Read input file
    Inputs inputs(id, input_file);
    std::string simulation_type = inputs.simulation_type;

    // Setup local and global grids, decomposing domain
    Grid grid(simulation_type, id, np, inputs.domain.number_of_layers, inputs.domain, inputs.temperature);

    // Material response function
    InterfacialResponseFunction irf(id, inputs.material_filename, inputs.domain.deltat, grid.deltax);

    // Ensure that input powder layer init options are compatible with this domain size, if needed for this problem type
    if ((simulation_type == "R") || (simulation_type == "S"))
        inputs.checkPowderOverflow(grid.nx, grid.ny, grid.layer_height, grid.number_of_layers);

    // Temperature fields characterized by data in this structure
    Temperature<memory_space> temperature(grid, inputs.temperature);
    // Read temperature data if necessary
    if (simulation_type == "R")
        temperature.readTemperatureData(id, grid, 0);
    // Initialize the temperature fields for the simualtion type of interest
    if ((simulation_type == "C") || (simulation_type == "SingleGrain"))
        temperature.initialize(id, simulation_type, grid, inputs.domain.deltat);
    else if (simulation_type == "S")
        temperature.initialize(id, grid, irf.freezing_range, inputs);
    else if (simulation_type == "R")
        temperature.initialize(0, id, grid, irf.freezing_range, inputs.domain.deltat);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Done with temperature field initialization" << std::endl;

    // Initialize grain orientations
    Orientation<memory_space> orientation(inputs.grain_orientation_file, false);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Done with orientation initialization " << std::endl;

    // Initialize cell types, grain IDs, and layer IDs
    CellData<memory_space> celldata(grid.domain_size, grid.domain_size_all_layers, inputs.substrate);
    if (simulation_type == "C")
        celldata.initSubstrate(id, grid, inputs.rng_seed);
    else if (simulation_type == "SingleGrain")
        celldata.initSubstrate(id, grid);
    else
        celldata.initSubstrate(id, grid, inputs.rng_seed, temperature.number_of_solidification_events);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Grain struct initialized" << std::endl;

    // Variables characterizing the active cell region within each rank's grid, including buffers for ghost node data
    // (fixed size) and the steering vector/steering vector size on host/device
    Interface<memory_space> interface(grid.domain_size);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Done with interface struct initialization " << std::endl;

    // Nucleation data structure, containing views of nuclei locations, time steps, and ids, and nucleation event
    // counters - initialized with an estimate on the number of nuclei in the layer Without knowing
    // PossibleNuclei_ThisRankThisLayer yet, initialize nucleation data structures to estimated sizes, resize inside of
    // NucleiInit when the number of nuclei per rank is known
    int estimated_nuclei_this_rank_this_layer = inputs.nucleation.n_max * pow(grid.deltax, 3) * grid.domain_size;
    Nucleation<memory_space> nucleation(estimated_nuclei_this_rank_this_layer, grid.deltax, inputs.nucleation);
    // Fill in nucleation data structures, and assign nucleation undercooling values to potential nucleation events
    // Potential nucleation grains are only associated with liquid cells in layer 0 - they will be initialized for each
    // successive layer when layer 0 in complete
    nucleation.placeNuclei(temperature, inputs.rng_seed, 0, grid, id);

    // Initialize printing struct from inputs
    Print print(grid, np, inputs.print);

    // If specified, print initial values in some views for debugging purposes
    double init_time = MPI_Wtime() - start_init_time;
    if (id == 0)
        std::cout << "Data initialized: Time spent: " << init_time << " s" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    int cycle = 0;
    double start_run_time = MPI_Wtime();
    for (int layernumber = 0; layernumber < grid.number_of_layers; layernumber++) {

        int x_switch = 0;
        double layer_time_1 = MPI_Wtime();

        // Loop continues until all liquid cells claimed by solid grains
        do {

            // Start of time step - optional print current state of ExaCA simulation (up to and including the current
            // layer's data)
            print.printIntralayer(id, np, layernumber, inputs.domain.deltat, cycle, grid, celldata, temperature,
                                  interface, orientation);
            cycle++;

            // Update cells on GPU - undercooling and diagonal length updates, nucleation
            // Cells with a successful nucleation event are marked and added to a steering vector, later dealt with in
            // CellCapture
            start_nucl_time = MPI_Wtime();
            nucleation.nucleateGrain(cycle, grid, celldata, interface);
            nucl_time += MPI_Wtime() - start_nucl_time;

            // Update cells on GPU - new active cells, solidification of old active cells
            // Cell capture performed in two steps - first, adding cells of interest to a steering vector (different
            // subroutine called with versus without remelting), and second, iterating over the steering vector to
            // perform active cell creation and capture operations
            // Constrained/directional solidification problem cells only undergo 1 solidificaton event and are
            // initialized as liquid - the steering vector operation for this problem can be constructed using
            // FillSteeringVector_NoRemelt (a simplified version of FillSteeringVector_Remelt
            start_create_sv_time = MPI_Wtime();

            if ((simulation_type == "C") || (simulation_type == "SingleGrain"))
                fillSteeringVector_NoRemelt(cycle, grid, celldata, temperature, interface);
            else
                fillSteeringVector_Remelt(cycle, grid, celldata, temperature, interface);
            create_sv_time += MPI_Wtime() - start_create_sv_time;

            start_capture_time = MPI_Wtime();
            // Cell capture and checking of the MPI buffers to ensure that all appropriate interface updates in the halo
            // regions were recorded
            cellCapture(cycle, np, grid, irf, celldata, temperature, interface, orientation);
            checkBuffers(id, cycle, grid, celldata, interface, orientation.n_grain_orientations);
            capture_time += MPI_Wtime() - start_capture_time;

            if (np > 1) {
                // Update ghost nodes
                start_ghost_time = MPI_Wtime();
                haloUpdate(cycle, id, grid, celldata, interface, orientation);
                ghost_time += MPI_Wtime() - start_ghost_time;
            }

            if ((cycle % 1000 == 0) && (simulation_type != "SingleGrain")) {
                intermediateOutputAndCheck(id, np, cycle, grid, nucleation.successful_nucleation_counter, x_switch,
                                           celldata, temperature, inputs.simulation_type, layernumber, orientation,
                                           print, inputs.domain.deltat, interface);
            }
            else if (simulation_type == "SingleGrain") {
                intermediateOutputAndCheck(id, cycle, grid, x_switch, celldata.cell_type);
            }

        } while (x_switch == 0);
        if (layernumber != grid.number_of_layers - 1) {
            MPI_Barrier(MPI_COMM_WORLD);

            // Optional print current state of ExaCA
            print.printInterlayer(id, np, layernumber, grid, celldata, temperature, interface, orientation);
            if (id == 0)
                std::cout << "Layer " << layernumber << " finished solidification; initializing layer "
                          << layernumber + 1 << std::endl;

            // Determine new active cell domain size and offset from bottom of global domain
            grid.initNextLayer(id, simulation_type, layernumber + 1, inputs.domain.spot_radius);

            // For simulation type R, need to initialize new temperature field data for layer "layernumber + 1"
            // TODO: reorganize these temperature functions calls into a temperature.init_next_layer as done with the
            // substrate
            if (simulation_type == "R") {
                // If the next layer's temperature data isn't already stored, it should be read
                if (inputs.temperature.layerwise_temp_read)
                    temperature.readTemperatureData(id, grid, layernumber + 1);
                // Initialize next layer's temperature data
                temperature.initialize(layernumber + 1, id, grid, irf.freezing_range, inputs.domain.deltat);
            }

            // Reset solidification event counter of all cells to zeros for the next layer, resizing to number of cells
            // associated with the next layer, and get the subview for undercooling
            temperature.resetLayerEventsUndercooling(grid);
            // Resize and zero all view data relating to the active region from the last layer, in preparation for the
            // next layer
            interface.initNextLayer(grid.domain_size);
            MPI_Barrier(MPI_COMM_WORLD);

            // Sets up views, powder layer (if necessary), and cell types for the next layer of a multilayer problem
            celldata.initNextLayer(layernumber + 1, id, grid, inputs.rng_seed,
                                   temperature.number_of_solidification_events);

            // Initialize potential nucleation event data for next layer "layernumber + 1"
            // Views containing nucleation data will be resized to the possible number of nuclei on a given MPI rank for
            // the next layer
            nucleation.resetNucleiCounters(); // start counters at 0
            nucleation.placeNuclei(temperature, inputs.rng_seed, layernumber + 1, grid, id);

            x_switch = 0;
            MPI_Barrier(MPI_COMM_WORLD);
            double layer_time_2 = MPI_Wtime();
            cycle = 0;
            if (id == 0)
                std::cout << "Time for layer number " << layernumber << " was " << layer_time_2 - layer_time_1
                          << " s, starting layer " << layernumber + 1 << std::endl;
        }
        else {
            MPI_Barrier(MPI_COMM_WORLD);
            double layer_time_2 = MPI_Wtime();
            if (id == 0)
                std::cout << "Time for final layer was " << layer_time_2 - layer_time_1 << " s" << std::endl;
        }
    }

    double run_time = MPI_Wtime() - start_run_time;
    double start_out_time = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);
    // Collect and print specified final fields to output files
    print.printInterlayer(id, np, grid.number_of_layers - 1, grid, celldata, temperature, interface, orientation);

    // Calculate volume fraction of solidified domain consisting of nucleated grains
    float vol_fraction_nucleated = celldata.calcVolFractionNucleated(id, grid);

    double out_time = MPI_Wtime() - start_out_time;
    double init_max_time, init_min_time, out_max_time, out_min_time = 0.0;
    double nucl_max_time, nucl_min_time, create_sv_min_time, create_sv_max_time, capture_max_time, capture_min_time,
        ghost_max_time, ghost_min_time = 0.0;
    MPI_Allreduce(&init_time, &init_max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&init_time, &init_min_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&nucl_time, &nucl_max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&nucl_time, &nucl_min_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&create_sv_time, &create_sv_max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&create_sv_time, &create_sv_min_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&capture_time, &capture_max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&capture_time, &capture_min_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&ghost_time, &ghost_max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&ghost_time, &ghost_min_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&out_time, &out_max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&out_time, &out_min_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    // Print the log file with JSON format, timing information to the console
    inputs.printExaCALog(id, np, input_file, grid.ny_local, grid.y_offset, irf, grid.deltax, grid.number_of_layers,
                         grid.layer_height, grid.nx, grid.ny, grid.nz, init_time, run_time, out_time, cycle,
                         init_max_time, init_min_time, nucl_max_time, nucl_min_time, create_sv_min_time,
                         create_sv_max_time, capture_max_time, capture_min_time, ghost_max_time, ghost_min_time,
                         out_max_time, out_min_time, grid.x_min, grid.x_max, grid.y_min, grid.y_max, grid.z_min,
                         grid.z_max, vol_fraction_nucleated);
    if (id == 0)
        printExaCATiming(np, init_time, run_time, out_time, cycle, init_max_time, init_min_time, nucl_max_time,
                         nucl_min_time, create_sv_min_time, create_sv_max_time, capture_max_time, capture_min_time,
                         ghost_max_time, ghost_min_time, out_max_time, out_min_time);
}
