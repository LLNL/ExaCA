// Copyright Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "ExaCA.hpp"
#include "Finch_Core.hpp"

#include <Kokkos_Core.hpp>

#include "mpi.h"

#include <cctype>
#include <stdexcept>
#include <string>

// Run Finch and store solidification event data for use by ExaCA. Either called for a specific layer, or for a series
// of layers where the return view stores time-temperature histories for layers l=0,1,2... where first_value[l] is the
// index of the first event and last_value[l]-1 is the index of the last event associated with layer l
Kokkos::View<double **, Kokkos::LayoutLeft, Kokkos::HostSpace>
getFinchData(const int, const int, const int first_finch_simulation, const int num_finch_simulations,
             const int number_of_layers, Inputs exaca_inputs, std::array<double, 3> &exaca_low_corner,
             std::array<double, 3> &exaca_high_corner, std::vector<int> &first_value_finch,
             std::vector<int> &last_value_finch) {

    using exec_space = Kokkos::DefaultExecutionSpace;
    using memory_space = exec_space::memory_space;
    Kokkos::View<double **, Kokkos::LayoutLeft, Kokkos::HostSpace> input_temperature_data(
        Kokkos::ViewAllocateWithoutInitializing("finch_temperature_data"), 0);
    // Run Finch, storing solidification data, domain bounds, and the locations in the view where each layer's data
    // starts and ends
    const int last_finch_simulation = first_finch_simulation + num_finch_simulations;
    for (int finch_input_file_num = first_finch_simulation; finch_input_file_num < last_finch_simulation;
         finch_input_file_num++) {
        // Setup Finch simulation
        Finch::Inputs finch_inputs(MPI_COMM_WORLD, exaca_inputs.temperature.temp_paths[finch_input_file_num]);

        // initialize a moving beam
        Finch::MovingBeam beam(finch_inputs.source.scan_path_file);

        // Define boundary condition details.
        std::array<std::string, 6> bc_types = {"adiabatic", "adiabatic", "adiabatic",
                                               "adiabatic", "adiabatic", "adiabatic"};

        // create the global mesh
        Finch::Grid<memory_space> finch_grid(MPI_COMM_WORLD, finch_inputs.space.cell_size,
                                             finch_inputs.space.global_low_corner,
                                             finch_inputs.space.global_high_corner, finch_inputs.space.ranks_per_dim,
                                             bc_types, finch_inputs.space.initial_temperature);

        // Run Finch heat transfer
        Finch::Layer app(finch_inputs, finch_grid);
        auto fd = Finch::createSolver(finch_inputs, finch_grid);
        app.run(exec_space(), finch_inputs, finch_grid, beam, fd);

        // Bounds of Finch domain for this layer stored for use by ExaCA - overwrite with bounds of Finch events if trim
        // options were specified
        std::array<double, 3> exaca_low_corner_layer = finch_inputs.space.global_low_corner;
        std::array<double, 3> exaca_high_corner_layer = finch_inputs.space.global_high_corner;
        std::array<double, 3> finch_low_corner_layer = app.getLowerSolidificationDataBounds();
        std::array<double, 3> finch_high_corner_layer = app.getUpperSolidificationDataBounds();
        if (exaca_inputs.temperature.trim_unmelted_region_xy) {
            exaca_low_corner_layer[0] = finch_low_corner_layer[0];
            exaca_high_corner_layer[0] = finch_high_corner_layer[0];
            exaca_low_corner_layer[1] = finch_low_corner_layer[1];
            exaca_high_corner_layer[1] = finch_high_corner_layer[1];
        }
        if (exaca_inputs.temperature.trim_unmelted_region_z) {
            exaca_low_corner_layer[2] = finch_low_corner_layer[2];
            exaca_high_corner_layer[2] = finch_high_corner_layer[2];
        }

        // Extend overall Finch bounds if necessary
        if (finch_input_file_num == 0) {
            for (int i = 0; i < 3; i++) {
                exaca_low_corner[i] = exaca_low_corner_layer[i];
                exaca_high_corner[i] = exaca_high_corner_layer[i];
            }
        }
        else {
            for (int i = 0; i < 3; i++) {
                exaca_low_corner[i] = std::min(exaca_low_corner_layer[i], exaca_low_corner[i]);
                exaca_high_corner[i] = std::min(exaca_high_corner_layer[i], exaca_high_corner[i]);
            }
        }

        // Append this layer's solidification data to input_temperature_data
        app.appendSolidificationData(input_temperature_data, first_value_finch, last_value_finch, finch_input_file_num,
                                     num_finch_simulations);
        // If performing multiple finch simulations during initialization, fill first/last_value_finch with values for
        // repeated data
        if (num_finch_simulations > 1) {
            for (int repeated_layer = num_finch_simulations; repeated_layer < number_of_layers; repeated_layer++) {
                const int repeated_file = repeated_layer % num_finch_simulations;
                first_value_finch[repeated_layer] = first_value_finch[repeated_file];
                last_value_finch[repeated_layer] = last_value_finch[repeated_file];
            }
        }
    }
    return input_temperature_data;
}

int main(int argc, char *argv[]) {
    // Initialize Kokkos & MPI
    int id, np;
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);
    {
        // Get number of processes
        MPI_Comm_size(MPI_COMM_WORLD, &np);
        // Get individual process ID
        MPI_Comm_rank(MPI_COMM_WORLD, &id);
        using exec_space = Kokkos::DefaultExecutionSpace;
        using memory_space = exec_space::memory_space;

        if (id == 0) {
            std::cout << "ExaCA version: " << version() << " \nExaCA commit:  " << gitCommitHash()
                      << "\nKokkos version: " << kokkosVersion() << std::endl;
            Kokkos::DefaultExecutionSpace().print_configuration(std::cout);
            std::cout << "Number of MPI ranks = " << np << std::endl;
        }
        std::string exaca_input_file;
        if (argc < 2)
            throw std::runtime_error("Error: Must provide path to input file on the command line.");
        else if (argc == 2)
            exaca_input_file = argv[1];
        else
            exaca_input_file = argv[2];

        // Setup ExaCA simulation details
        Inputs inputs(id, exaca_input_file);

        // Without Finch input file on the command line, ensure that at least 1 Finch input file is listed in the ExaCA
        // input file. Otherwise store Finch input filename given on the command line within the ExaCA inputs struct
        if ((argc == 2) && (inputs.temperature.temp_files_in_series == 0))
            throw std::runtime_error(
                "Error: Finch input file must be given on the command line if absent from the ExaCA input file");
        else if (argc > 2) {
            std::string finch_input_file = argv[1];
            inputs.temperature.temp_files_in_series = 1;
            inputs.temperature.temp_paths.push_back(finch_input_file);
            if (id == 0)
                std::cout << "Warning: Ability to specify Finch input file on the command line is deprecated and will "
                             "be removed in a future release, please use the FinchInputFile input list in the "
                             "TemperatureData portion of the ExaCA input file to specify desired Finch input file(s)"
                          << std::endl;
        }

        // Create timers
        Timers timers(id);

        // ExaCA either uses the same bounds as Finch, or a bounding box that only includes the regions that underwent
        // melting and solidification
        std::array<double, 3> exaca_low_corner, exaca_high_corner;
        // Should Finch simulations for all layers be performed and time-temperature history data stored prior to
        // running ExaCA, or should each Finch simulation be run one at a time between ExaCA-simulated layers?
        int num_finch_simulations;
        if (inputs.temperature.layerwise_temp_read)
            num_finch_simulations = 1;
        else
            num_finch_simulations = inputs.temperature.temp_files_in_series;
        // Store values denoting which portions of input_temperature_data belong to which simulated layer in Finch
        std::vector<int> first_value_finch(inputs.domain.number_of_layers),
            last_value_finch(inputs.domain.number_of_layers);
        // View for time-temperature history data from Finch - run either just the first layer, or run all Finch
        // simulations and store results
        timers.startHeatTransfer();
        Kokkos::View<double **, Kokkos::LayoutLeft, Kokkos::HostSpace> input_temperature_data =
            getFinchData(id, np, 0, num_finch_simulations, inputs.domain.number_of_layers, inputs, exaca_low_corner,
                         exaca_high_corner, first_value_finch, last_value_finch);
        MPI_Barrier(MPI_COMM_WORLD);
        timers.stopHeatTransfer();

        // Start of ExaCA init time.
        timers.startInit();
        // Initialize ExaCA grid, using bounds of input Finch data
        Grid grid(id, np, inputs.domain, inputs.substrate, inputs.temperature, inputs.domain.deltax, exaca_low_corner,
                  exaca_high_corner);

        // Material response function
        InterfacialResponseFunction irf(inputs.domain.deltat, grid.deltax, inputs.irf);

        // Temperature fields characterized by data in this structure - rearrange data in input_temperature_data to
        // match ExaCA's domain decomposition. Store which temperature data is associated with which layers
        Temperature<memory_space> temperature(id, np, grid, inputs.temperature, input_temperature_data,
                                              first_value_finch, last_value_finch,
                                              inputs.print.store_solidification_start);
        temperature.initialize(0, id, grid, irf.freezingRange(), inputs.domain.deltat, "FromFinch");
        MPI_Barrier(MPI_COMM_WORLD);

        // Initialize grain orientations
        Orientation<memory_space> orientation(id, inputs.grain_orientation_file, false);
        MPI_Barrier(MPI_COMM_WORLD);

        // Initialize cell types, grain IDs, and layer IDs
        CellData<memory_space> celldata(grid, inputs.substrate, inputs.print.store_melt_pool_edge);
        celldata.initSubstrate(id, grid, inputs.rng_seed, temperature.number_of_solidification_events);
        MPI_Barrier(MPI_COMM_WORLD);

        // Variables characterizing the active cell region within each rank's grid, including buffers for ghost node
        // data (fixed size) and the steering vector/steering vector size on host/device
        Interface<memory_space> interface(id, grid.domain_size, inputs.substrate.init_oct_size);
        MPI_Barrier(MPI_COMM_WORLD);

        // Nucleation data structure, containing views of nuclei locations, time steps, and ids, and nucleation event
        // counters - initialized with an estimate on the number of nuclei in the layer Without knowing
        // estimated_nuclei_this_rank_this_layer yet, initialize nucleation data structures to estimated sizes, resize
        // inside of placeNuclei when the number of nuclei per rank is known
        int estimated_nuclei_this_rank_this_layer = inputs.nucleation.n_max * pow(grid.deltax, 3) * grid.domain_size;
        Nucleation<memory_space> nucleation(estimated_nuclei_this_rank_this_layer, inputs.nucleation);
        // Fill in nucleation data structures, and assign nucleation undercooling values to potential nucleation events
        // Potential nucleation grains are only associated with liquid cells in layer 0 - they will be initialized for
        // each successive layer when layer 0 is complete
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
                          interface, nucleation, print, "FromFinch");

            if (layernumber != grid.number_of_layers - 1) {

                // Initialize new temperature field data for layer "layernumber + 1" - this is either stored in the
                // temperature struct, or comes from a new Finch simulation
                if (inputs.temperature.layerwise_temp_read) {
                    timers.startHeatTransfer();
                    input_temperature_data =
                        getFinchData(id, np, layernumber + 1, num_finch_simulations, inputs.domain.number_of_layers,
                                     inputs, exaca_low_corner, exaca_high_corner, first_value_finch, last_value_finch);
                    timers.stopHeatTransfer();
                    temperature.copyTemperatureData(id, np, layernumber + 1, grid, input_temperature_data,
                                                    first_value_finch, last_value_finch);
                }
                MPI_Barrier(MPI_COMM_WORLD);
                // Initialize next layer's temperature data
                temperature.initialize(layernumber + 1, id, grid, irf.freezingRange(), inputs.domain.deltat,
                                       "FromFinch");

                // Reset solidification event counter of all cells to zeros for the next layer, resizing to number of
                // cells associated with the next layer, and get the subview for undercooling (TODO: part of
                // temperature.initialize in the future)
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

    // Finalize Kokkos & MPI
    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}
