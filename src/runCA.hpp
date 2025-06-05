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

// Run ExaCA using thermal data for a given problem type
template <typename MemorySpace>
void runExaCALayer(int id, int np, int layernumber, int &cycle, Inputs inputs, Timers &timers, Grid &grid,
                   Temperature<MemorySpace> temperature, InterfacialResponseFunction irf,
                   Orientation<MemorySpace> &orientation, CellData<MemorySpace> &celldata,
                   Interface<MemorySpace> &interface, Nucleation<MemorySpace> &nucleation, Print &print,
                   std::string simulation_type) {

    int x_switch = 0;

    // Loop continues until all liquid cells claimed by solid grains, and no solid cells undergo remelting
    do {

        // Start of time step - optional print current state of ExaCA simulation (up to and including the current
        // layer's data)
        print.printIntralayer(id, np, layernumber, inputs.domain.deltat, cycle, grid, celldata, temperature, interface,
                              orientation);
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
                                       celldata, temperature, inputs.simulation_type, layernumber, orientation, print,
                                       inputs.domain.deltat, interface);
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
        print.printInterlayer(id, np, layernumber, inputs.domain.deltat, cycle, grid, celldata, temperature, interface,
                              orientation);

        // If more layers to simulate, determine new active cell domain size and offset from bottom of global domain
        grid.initNextLayer(id, simulation_type, layernumber + 1);
    }
}

// Initialize structs for the next layer of a multilayer ExaCA simulation
template <typename MemorySpace>
void initExaCALayer(int id, int layernumber, int &cycle, Inputs inputs, Grid &grid,
                    Temperature<MemorySpace> temperature, CellData<MemorySpace> &celldata,
                    Interface<MemorySpace> &interface, Nucleation<MemorySpace> &nucleation) {

    if (layernumber != grid.number_of_layers - 1) {
        // Resize and zero all view data relating to the active region from the last layer, in preparation for the
        // next layer
        interface.initNextLayer(grid.domain_size);
        MPI_Barrier(MPI_COMM_WORLD);

        // Sets up views, powder layer (if necessary), and cell types for the next layer of a multilayer problem
        celldata.initNextLayer(layernumber + 1, id, grid, inputs.rng_seed, temperature.number_of_solidification_events);

        // Initialize potential nucleation event data for next layer "layernumber + 1"
        // Views containing nucleation data will be resized to the possible number of nuclei on a given MPI rank for
        // the next layer
        nucleation.resetNucleiCounters(); // start counters at 0
        nucleation.placeNuclei(temperature, inputs.rng_seed, layernumber + 1, grid, id);

        // If not on the last layer, reset the time step counter to 0 (if on the last layer, the time step is stored and
        // printed in the log file)
        cycle = 0;
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

// Finalize timer values, print end-of-run output to the console/files, and print the log file
template <typename MemorySpace>
void finalizeExaCA(int id, int np, int cycle, Inputs inputs, Timers timers, Grid &grid,
                   Temperature<MemorySpace> temperature, Orientation<MemorySpace> &orientation,
                   CellData<MemorySpace> &celldata, Interface<MemorySpace> &interface, Print &print) {

    timers.startOutput();

    // Collect and print specified final fields to output files
    print.printInterlayer(id, np, grid.number_of_layers - 1, inputs.domain.deltat, cycle, grid, celldata, temperature,
                          interface, orientation);

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
