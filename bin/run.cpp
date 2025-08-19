// Copyright Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "ExaCA.hpp"

#include <Kokkos_Core.hpp>

#include "mpi.h"

#include <stdexcept>
#include <string>

int main(int argc, char *argv[]) {
    // Initialize MPI
    int id, np;
    MPI_Init(&argc, &argv);
    // Initialize Kokkos
    Kokkos::initialize(argc, argv);
    {
        using memory_space = Kokkos::DefaultExecutionSpace::memory_space;

        // Get number of processes
        MPI_Comm_size(MPI_COMM_WORLD, &np);
        // Get individual process ID
        MPI_Comm_rank(MPI_COMM_WORLD, &id);

        if (id == 0) {
            std::cout << "ExaCA version: " << version() << " \nExaCA commit:  " << gitCommitHash()
                      << "\nKokkos version: " << kokkosVersion() << std::endl;
            Kokkos::DefaultExecutionSpace().print_configuration(std::cout);
            std::cout << "Number of MPI ranks = " << np << std::endl;
        }
        if (argc < 2) {
            throw std::runtime_error("Error: Must provide path to input file on the command line.");
        }
        else {

            // Create timers
            Timers timers(id);
            timers.startInit();

            // Run CA code using reduced temperature data format
            std::string input_file = argv[1];
            // Read input file
            Inputs inputs(id, input_file);
            std::string simulation_type = inputs.simulation_type;

            // Setup local and global grids, decomposing domain (needed to construct temperature)
            Grid grid(inputs.simulation_type, id, np, inputs.domain.number_of_layers, inputs.domain, inputs.substrate,
                      inputs.temperature);
            // Temperature fields characterized by data in this structure
            Temperature<memory_space> temperature(grid, inputs.temperature, inputs.print.store_solidification_start);

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

            // Variables characterizing the active cell region within each rank's grid, including buffers for ghost node
            // data (fixed size) and the steering vector/steering vector size on host/device
            Interface<memory_space> interface(id, grid.domain_size, inputs.substrate.init_oct_size);
            MPI_Barrier(MPI_COMM_WORLD);

            // Nucleation data structure, containing views of nuclei locations, time steps, and ids, and nucleation
            // event counters - initialized with an estimate on the number of nuclei in the layer Without knowing
            // estimated_nuclei_this_rank_this_layer yet, initialize nucleation data structures to estimated sizes,
            // resize inside of placeNuclei when the number of nuclei per rank is known
            int estimated_nuclei_this_rank_this_layer =
                inputs.nucleation.n_max * pow(grid.deltax, 3) * grid.domain_size;
            Nucleation<memory_space> nucleation(estimated_nuclei_this_rank_this_layer, inputs.nucleation);
            // Fill in nucleation data structures, and assign nucleation undercooling values to potential nucleation
            // events Potential nucleation grains are only associated with liquid cells in layer 0 - they will be
            // initialized for each successive layer when layer 0 is complete
            nucleation.placeNuclei(temperature, inputs.rng_seed, 0, grid, id);

            // Initialize printing struct from inputs
            Print print(grid, np, inputs.print);

            // End of initialization
            timers.stopInit();
            MPI_Barrier(MPI_COMM_WORLD);

            int cycle = 0;
            timers.startRun();

            // Run ExaCA to model solidification of each layer
            for (int layernumber = 0; layernumber < grid.number_of_layers; layernumber++) {
                timers.startLayer();
                runExaCALayer(id, np, layernumber, cycle, inputs, timers, grid, temperature, irf, orientation, celldata,
                              interface, nucleation, print, simulation_type);

                if (layernumber != grid.number_of_layers - 1) {
                    // Initialize new temperature field data for layer "layernumber + 1"
                    // TODO: reorganize these temperature functions calls into a temperature.init_next_layer as done
                    // with the substrate If the next layer's temperature data isn't already stored, it should be read
                    if ((simulation_type == "FromFile") && (inputs.temperature.layerwise_temp_read))
                        temperature.readTemperatureData(id, grid, layernumber + 1);
                    MPI_Barrier(MPI_COMM_WORLD);
                    // Initialize next layer's temperature data
                    temperature.initialize(layernumber + 1, id, grid, irf.freezingRange(), inputs.domain.deltat,
                                           simulation_type);
                    // Reset solidification event counter of all cells to zeros for the next layer, resizing to number
                    // of cells associated with the next layer, and get the subview for undercooling
                    temperature.resetLayerEventsUndercooling(grid);

                    // Initialize next layer of the simulation
                    initExaCALayer(id, layernumber, cycle, inputs, grid, temperature, celldata, interface, nucleation);
                    timers.stopLayer(layernumber);
                }
                else {
                    MPI_Barrier(MPI_COMM_WORLD);
                    timers.stopLayer();
                }
            }
            timers.stopRun();
            MPI_Barrier(MPI_COMM_WORLD);

            // Print ExaCA end-of-run data
            finalizeExaCA(id, np, cycle, inputs, timers, grid, temperature, orientation, celldata, interface, print);
        }
    }
    // Finalize Kokkos
    Kokkos::finalize();
    // Finalize MPI
    MPI_Finalize();
    return 0;
}
