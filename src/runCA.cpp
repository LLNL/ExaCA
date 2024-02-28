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

    // Create timers
    Timers timers(id);
    timers.startInit();

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
    Temperature<memory_space> temperature(grid, inputs.temperature, inputs.print.store_solidification_start);
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

    // Initialize grain orientations
    Orientation<memory_space> orientation(id, inputs.grain_orientation_file, false);
    MPI_Barrier(MPI_COMM_WORLD);

    // Initialize cell types, grain IDs, and layer IDs
    CellData<memory_space> celldata(grid.domain_size, grid.domain_size_all_layers, inputs.substrate);
    if (simulation_type == "C")
        celldata.initSubstrate(id, grid, inputs.rng_seed);
    else if (simulation_type == "SingleGrain")
        celldata.initSubstrate(id, grid);
    else
        celldata.initSubstrate(id, grid, inputs.rng_seed, temperature.number_of_solidification_events);
    MPI_Barrier(MPI_COMM_WORLD);

    // Variables characterizing the active cell region within each rank's grid, including buffers for ghost node data
    // (fixed size) and the steering vector/steering vector size on host/device
    Interface<memory_space> interface(id, grid.domain_size);
    MPI_Barrier(MPI_COMM_WORLD);

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

    // End of initialization
    timers.stopInit();
    MPI_Barrier(MPI_COMM_WORLD);

    int cycle = 0;
    timers.startRun();
    for (int layernumber = 0; layernumber < grid.number_of_layers; layernumber++) {

        int x_switch = 0;
        timers.startLayer();

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
            timers.startNucleation();
            nucleation.nucleateGrain(cycle, grid, celldata, interface);
            timers.stopNucleation();

            // Update cells on GPU - new active cells, solidification of old active cells
            // Cell capture performed in two steps - first, adding cells of interest to a steering vector (different
            // subroutine called with versus without remelting), and second, iterating over the steering vector to
            // perform active cell creation and capture operations
            // Constrained/directional solidification problem cells only undergo 1 solidificaton event and are
            // initialized as liquid - the steering vector operation for this problem can be constructed using
            // FillSteeringVector_NoRemelt (a simplified version of FillSteeringVector_Remelt
            timers.startSV();
            if ((simulation_type == "C") || (simulation_type == "SingleGrain"))
                fillSteeringVector_NoRemelt(cycle, grid, celldata, temperature, interface);
            else
                fillSteeringVector_Remelt(cycle, grid, celldata, temperature, interface);
            timers.stopSV();

            timers.startCapture();
            // Cell capture and checking of the MPI buffers to ensure that all appropriate interface updates in the halo
            // regions were recorded
            cellCapture(cycle, np, grid, irf, celldata, temperature, interface, orientation);
            checkBuffers(id, cycle, grid, celldata, interface, orientation.n_grain_orientations);
            timers.stopCapture();

            if (np > 1) {
                // Update ghost nodes
                timers.startComm();
                haloUpdate(cycle, id, grid, celldata, interface, orientation);
                timers.stopComm();
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
            cycle = 0;
            timers.stopLayer(layernumber);
        }
        else {
            MPI_Barrier(MPI_COMM_WORLD);
            timers.stopLayer();
        }
    }
    timers.stopRun();
    MPI_Barrier(MPI_COMM_WORLD);
    timers.startOutput();

    // Collect and print specified final fields to output files
    print.printInterlayer(id, np, grid.number_of_layers - 1, grid, celldata, temperature, interface, orientation);

    // Calculate volume fraction of solidified domain consisting of nucleated grains
    float vol_fraction_nucleated = celldata.calcVolFractionNucleated(id, grid);
    // Get MPI timing statisticss
    timers.stopOutput();
    timers.reduceMPI();

    // Print the log file with JSON format
    inputs.printExaCALog(id, np, input_file, grid.ny_local, grid.y_offset, irf, grid.deltax, grid.number_of_layers,
                         grid.layer_height, grid.nx, grid.ny, grid.nz, timers, cycle, grid.x_min, grid.x_max,
                         grid.y_min, grid.y_max, grid.z_min, grid.z_max, vol_fraction_nucleated);

    // Print timing information to the console
    timers.printFinal(np, cycle);
}
