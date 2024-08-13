// Copyright Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef RUNCA_HPP
#define RUNCA_HPP

#include "ExaCA.hpp"

#include "mpi.h"

#include <string>
#include <vector>

template <typename MemorySpace>
void runExaCA(int id, int np, Inputs inputs, Timers timers, Grid grid, Temperature<MemorySpace> temperature) {

    // Run on the default space.
    using memory_space = MemorySpace;

    std::string simulation_type = inputs.simulation_type;

    // Material response function
    InterfacialResponseFunction irf(inputs.domain.deltat, grid.deltax, inputs.irf);

    // Read temperature data if necessary
    if (simulation_type == "FromFile")
        temperature.readTemperatureData(id, grid, 0);
    // Initialize the temperature fields for the simulation type of interest
    if ((simulation_type == "Directional") || (simulation_type == "SingleGrain"))
        temperature.initialize(id, simulation_type, grid, inputs.domain.deltat);
    else if (simulation_type == "Spot")
        temperature.initialize(id, grid, irf.freezingRange(), inputs.domain.deltat, inputs.domain.spot_radius);
    else if ((simulation_type == "FromFile") || (simulation_type == "FromFinch"))
        temperature.initialize(0, id, grid, irf.freezingRange(), inputs.domain.deltat, simulation_type);
    MPI_Barrier(MPI_COMM_WORLD);

    // Initialize grain orientations
    Orientation<memory_space> orientation(id, inputs.grain_orientation_file, false);
    MPI_Barrier(MPI_COMM_WORLD);

    // Initialize cell types, grain IDs, and layer IDs
    CellData<memory_space> celldata(grid, inputs.substrate, inputs.print.store_melt_pool_edge);
    if (simulation_type == "Directional")
        celldata.initSubstrate(id, grid, inputs.rng_seed);
    else if (simulation_type == "SingleGrain")
        celldata.initSubstrate(id, grid);
    else
        celldata.initSubstrate(id, grid, inputs.rng_seed, temperature.number_of_solidification_events);
    MPI_Barrier(MPI_COMM_WORLD);

    // Variables characterizing the active cell region within each rank's grid, including buffers for ghost node data
    // (fixed size) and the steering vector/steering vector size on host/device
    Interface<memory_space> interface(id, grid.domain_size, inputs.substrate.init_oct_size);
    MPI_Barrier(MPI_COMM_WORLD);

    // Nucleation data structure, containing views of nuclei locations, time steps, and ids, and nucleation event
    // counters - initialized with an estimate on the number of nuclei in the layer Without knowing
    // estimated_nuclei_this_rank_this_layer yet, initialize nucleation data structures to estimated sizes, resize
    // inside of placeNuclei when the number of nuclei per rank is known
    int estimated_nuclei_this_rank_this_layer = inputs.nucleation.n_max * pow(grid.deltax, 3) * grid.domain_size;
    Nucleation<memory_space> nucleation(estimated_nuclei_this_rank_this_layer, inputs.nucleation);
    // Fill in nucleation data structures, and assign nucleation undercooling values to potential nucleation events
    // Potential nucleation grains are only associated with liquid cells in layer 0 - they will be initialized for each
    // successive layer when layer 0 is complete
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

        // Loop continues until all liquid cells claimed by solid grains, and no solid cells undergo remelting
        do {

            // Start of time step - optional print current state of ExaCA simulation (up to and including the current
            // layer's data)
            print.printIntralayer(id, np, layernumber, inputs.domain.deltat, cycle, grid, celldata, temperature,
                                  interface, orientation);
            cycle++;

            // Cells with a successful nucleation event are marked and added to a steering vector, later dealt with in
            // cellCapture
            timers.startNucleation();
            nucleation.nucleateGrain(cycle, grid, celldata, interface);
            timers.stopNucleation();

            // Cells that have a successful nucleation event, and other cells that are at the solid-liquid interface are
            // added to a steering vector. Logic in fillSteeringVector_NoRemelt is a simpified version of
            // fillSteeringVector_Remelt
            timers.startSV();
            if ((simulation_type == "Directional") || (simulation_type == "SingleGrain"))
                fillSteeringVector_NoRemelt(cycle, grid, celldata, temperature, interface);
            else
                fillSteeringVector_Remelt(cycle, grid, celldata, temperature, interface);
            timers.stopSV();

            // Iterate over the steering vector to perform active cell creation and capture operations, and if needed,
            // melting of cells that have gone above the liquidus. Also places halo cell data into send buffers, later
            // checking the MPI buffers to ensure that all appropriate interface updates in the halo regions were
            // recorded
            timers.startCapture();
            cellCapture(cycle, np, grid, irf, celldata, temperature, interface, orientation);
            checkBuffers(id, cycle, grid, celldata, interface, orientation.n_grain_orientations);
            timers.stopCapture();

            if (np > 1) {
                // Update ghost nodes
                timers.startComm();
                haloUpdate(cycle, id, grid, celldata, interface, orientation);
                timers.stopComm();
            }

            // Check on progress of solidification simulation of the current layer, setting x_switch = 1 if complete
            if ((cycle % 1000 == 0) && (simulation_type != "SingleGrain")) {
                intermediateOutputAndCheck(id, np, cycle, grid, nucleation.successful_nucleation_counter, x_switch,
                                           celldata, temperature, inputs.simulation_type, layernumber, orientation,
                                           print, inputs.domain.deltat, interface);
            }
            else if (simulation_type == "SingleGrain") {
                intermediateOutputAndCheck(id, cycle, grid, x_switch, celldata.cell_type);
            }

        } while (x_switch == 0);

        // Reset intralayer print counter and print time series file for previous layer's intralayer data (if needed)
        print.resetIntralayer(id, layernumber);

        if (layernumber != grid.number_of_layers - 1) {
            MPI_Barrier(MPI_COMM_WORLD);

            // Optional print current state of ExaCA
            print.printInterlayer(id, np, layernumber, grid, celldata, temperature, interface, orientation);

            // Determine new active cell domain size and offset from bottom of global domain
            grid.initNextLayer(id, simulation_type, layernumber + 1);

            // Initialize new temperature field data for layer "layernumber + 1"
            // TODO: reorganize these temperature functions calls into a temperature.init_next_layer as done with the
            // substrate
            // If the next layer's temperature data isn't already stored, it should be read
            if ((simulation_type == "FromFile") && (inputs.temperature.layerwise_temp_read))
                temperature.readTemperatureData(id, grid, layernumber + 1);
            MPI_Barrier(MPI_COMM_WORLD);
            // Initialize next layer's temperature data
            temperature.initialize(layernumber + 1, id, grid, irf.freezingRange(), inputs.domain.deltat,
                                   simulation_type);

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
    inputs.printExaCALog(id, np, cycle, grid, timers, vol_fraction_nucleated);

    // Print timing information to the console
    timers.printFinal(np, cycle);
}

#endif
