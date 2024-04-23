// Copyright Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_PRINT_HPP
#define EXACA_PRINT_HPP

#include "CAcelldata.hpp"
#include "CAgrid.hpp"
#include "CAinputdata.hpp"
#include "CAinterface.hpp"
#include "CAinterfacialresponse.hpp"
#include "CAorientation.hpp"
#include "CAparsefiles.hpp"
#include "CAtemperature.hpp"
#include "CAtypes.hpp"
#include "mpi.h"

#include <Kokkos_Core.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Write data of type PrintType as ascii or binary, with option to convert between big and small endian binary
template <typename PrintType>
void writeData(std::ofstream &outstream, PrintType print_value, bool print_binary, bool swap_endian_yn = false) {
    if (print_binary) {
        if (swap_endian_yn)
            swapEndian(print_value);
        int var_size = sizeof(PrintType);
        outstream.write((char *)&print_value, var_size);
    }
    else
        outstream << print_value << " ";
}

// Struct to hold data printing options and functions
struct Print {

    // If printing data at the end of layers, the counter for the number of intermediate files that have been printed
    int interlayer_file_count = 0;

    // Message sizes and data offsets for data send/received to/from other ranks- message size different for different
    // ranks
    using view_type_int_host = Kokkos::View<int *, Kokkos::HostSpace>;
    using view_type_float_host = Kokkos::View<float *, Kokkos::HostSpace>;
    view_type_int_host recv_y_offset, recv_ny_local, recv_buf_size;
    // Y coordinates for a given rank's data being send/loaded into the view of all domain data on rank 0=
    int send_buf_start_y, send_buf_end_y, send_buf_size;
    // Holds print options from input file
    PrintInputs _inputs;
    // Combined path/file prefix for output files
    std::string path_base_filename;

    // Default constructor - options are set in getPrintDataFromFile and copied into this struct
    Print(const Grid &grid, const int np, PrintInputs inputs)
        : recv_y_offset(view_type_int_host(Kokkos::ViewAllocateWithoutInitializing("Recv_y_offset"), np))
        , recv_ny_local(view_type_int_host(Kokkos::ViewAllocateWithoutInitializing("Recv_ny_local"), np))
        , recv_buf_size(view_type_int_host(Kokkos::ViewAllocateWithoutInitializing("RBufSize"), np))
        , _inputs(inputs) {

        // Buffers for sending/receiving data across ranks
        for (int recvrank = 0; recvrank < np; recvrank++) {
            recv_y_offset(recvrank) = grid.getYOffset(recvrank, np);
            recv_ny_local(recvrank) = grid.getNyLocal(recvrank, np);
            recv_buf_size(recvrank) = grid.nx * recv_ny_local(recvrank) * grid.nz;
        }

        // Y coordinates for a given rank's data being send/loaded into the view of all domain data on rank 0
        if (grid.y_offset == 0)
            send_buf_start_y = 0;
        else
            send_buf_start_y = 1;
        if (grid.ny_local + grid.y_offset == grid.ny)
            send_buf_end_y = grid.ny_local;
        else
            send_buf_end_y = grid.ny_local - 1;
        send_buf_size = grid.nx * (send_buf_end_y - send_buf_start_y) * grid.nz;

        path_base_filename = _inputs.path_to_output + _inputs.base_filename;
    }

    // Called on rank 0 to collect view data from other ranks, or on other ranks to send data to rank 0
    // MPI datatype corresponding to the view data
    template <typename Collect1DViewTypeDevice>
    auto collectViewData(int id, int np, const Grid &grid, const bool current_layer_only, MPI_Datatype msg_type,
                         Collect1DViewTypeDevice view_data_this_rank_device) {

        // Either collect data from all Z coordinates, or just the ones associated with the current layer of a problem
        int z_print_size;
        if (current_layer_only)
            z_print_size = grid.nz_layer;
        else
            z_print_size = grid.nz;

        // View (int or float extracted from 1D View) for 3D data (no initial size, only given size/filled on rank 0)
        using value_type = typename Collect1DViewTypeDevice::value_type;
        Kokkos::View<value_type ***, Kokkos::HostSpace> view_data_whole_domain(
            Kokkos::ViewAllocateWithoutInitializing(view_data_this_rank_device.label() + "_WholeDomain"), 0, 0, 0);

        // Get the host view type
        using host_view_type = typename Collect1DViewTypeDevice::HostMirror;

        // Copy view data of interest host
        auto view_data_this_rank = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), view_data_this_rank_device);
        if (id == 0) {
            // View (short, int, or float) for 3D data given a size
            Kokkos::resize(view_data_whole_domain, z_print_size, grid.nx, grid.ny);
            // Place rank 0 data into view for whole domain
            for (int coord_z = 0; coord_z < z_print_size; coord_z++) {
                for (int coord_x = 0; coord_x < grid.nx; coord_x++) {
                    for (int coord_y_local = 0; coord_y_local < grid.ny_local; coord_y_local++) {
                        int index = grid.get1DIndex(coord_x, coord_y_local, coord_z);
                        view_data_whole_domain(coord_z, coord_x, coord_y_local) = view_data_this_rank(index);
                    }
                }
            }

            // Receive values from other ranks - message size different for different ranks
            for (int recvrank = 1; recvrank < np; recvrank++) {
                int recv_buf_size_this_rank = recv_buf_size(recvrank);
                host_view_type recv_buf(Kokkos::ViewAllocateWithoutInitializing("RecvBufData"),
                                        recv_buf_size_this_rank);
                MPI_Recv(recv_buf.data(), recv_buf_size_this_rank, msg_type, recvrank, 0, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
                int data_counter = 0;
                for (int coord_z = 0; coord_z < z_print_size; coord_z++) {
                    for (int coord_x = 0; coord_x < grid.nx; coord_x++) {
                        for (int coord_y_local = 0; coord_y_local < recv_ny_local(recvrank); coord_y_local++) {
                            int coord_y_global = coord_y_local + recv_y_offset(recvrank);
                            view_data_whole_domain(coord_z, coord_x, coord_y_global) = recv_buf(data_counter);
                            data_counter++;
                        }
                    }
                }
            }
        }
        else {
            // Send non-ghost node data to rank 0
            int data_counter = 0;
            host_view_type send_buf(Kokkos::ViewAllocateWithoutInitializing("SendBuf"), send_buf_size);
            for (int coord_z = 0; coord_z < z_print_size; coord_z++) {
                for (int coord_x = 0; coord_x < grid.nx; coord_x++) {
                    for (int coord_y_local = send_buf_start_y; coord_y_local < send_buf_end_y; coord_y_local++) {
                        int index = grid.get1DIndex(coord_x, coord_y_local, coord_z);
                        send_buf(data_counter) = view_data_this_rank(index);
                        data_counter++;
                    }
                }
            }
            MPI_Send(send_buf.data(), send_buf_size, msg_type, 0, 0, MPI_COMM_WORLD);
        }
        return view_data_whole_domain;
    }

    // Called on rank 0, prints current values of selected data structures to vtk files
    // For intralayer print (current_layer_only = true):
    // Creates a file "Basefilename_layer[layernumber]_[time].vtk" containing the portions of data structures
    // corresponding to the current layer of a multilayer simulation The exception to this is for print
    // GrainMisorientation, where data for all layer up to the current layer are printed to a separate file
    // "Basefilename_layer[layernumber]_[time]_Misorientations.vtk"
    template <typename MemorySpace>
    void printIntralayer(const int id, const int np, const int layernumber, const double deltat, const int cycle,
                         const Grid &grid, CellData<MemorySpace> &celldata, Temperature<MemorySpace> &temperature,
                         Interface<MemorySpace> &interface, Orientation<MemorySpace> &orientation) {

        using view_type_float = Kokkos::View<float *, MemorySpace>;
        using view_type_int = Kokkos::View<int *, MemorySpace>;
        if ((_inputs.intralayer) && (cycle % _inputs.intralayer_increment == 0)) {
            // Current time in microseconds
            double current_time = static_cast<double>(cycle) * deltat * Kokkos::pow(10, 6);
            std::string vtk_filename_current_layer_base =
                path_base_filename + "_layer" + std::to_string(layernumber) + "_" + std::to_string(current_time);
            std::string vtk_filename_current_layer = vtk_filename_current_layer_base + ".vtk";
            std::ofstream intralayer_ofstream;
            if (id == 0) {
                std::cout << "Printing data structures for current layer of the simulation to file "
                          << vtk_filename_current_layer << std::endl;
                // Write header for printing current layer's data to file
                writeHeader(intralayer_ofstream, vtk_filename_current_layer, grid, true);
            }
            // Print current layer vtk data (Z coordinates of overall domain spanning layer_range_z)
            if (_inputs.intralayer_grain_id) {
                auto grain_id_current_layer = celldata.getGrainIDSubview(grid);
                auto grain_id_whole_domain = collectViewData(id, np, grid, true, MPI_INT, grain_id_current_layer);
                printViewData(id, intralayer_ofstream, grid, true, "int", "GrainID", grain_id_whole_domain);
            }
            if (_inputs.intralayer_layer_id) {
                auto layer_id_current_layer = celldata.getLayerIDSubview(grid);
                auto layer_id_whole_domain = collectViewData(id, np, grid, true, MPI_SHORT, layer_id_current_layer);
                printViewData(id, intralayer_ofstream, grid, true, "short", "LayerID", layer_id_whole_domain);
            }
            if (_inputs.intralayer_undercooling_current) {
                // TODO: remove this, the subview undercooling_current_layer is already stored in the temperature struct
                auto undercooling_current_layer =
                    Kokkos::subview(temperature.undercooling_current_all_layers, grid.layer_range);
                auto undercooling_whole_domain =
                    collectViewData(id, np, grid, true, MPI_FLOAT, undercooling_current_layer);
                printViewData(id, intralayer_ofstream, grid, true, "float", "UndercoolingCurrent",
                              undercooling_whole_domain);
            }
            if (_inputs.intralayer_undercooling_solidification_start) {
                auto undercooling_start_whole_domain =
                    collectViewData(id, np, grid, true, MPI_FLOAT, temperature.undercooling_solidification_start);
                printViewData(id, intralayer_ofstream, grid, true, "float", "UndercoolingStart",
                              undercooling_start_whole_domain);
            }
            if (_inputs.intralayer_melt_time_step) {
                auto melt_time_step = temperature.template extractTmTlCrData<view_type_int>(0, grid.domain_size);
                auto melt_time_step_whole_domain = collectViewData(id, np, grid, true, MPI_INT, melt_time_step);
                printViewData(id, intralayer_ofstream, grid, true, "int", "MeltTimeStep", melt_time_step_whole_domain);
            }
            if (_inputs.intralayer_crit_time_step) {
                auto crit_time_step = temperature.template extractTmTlCrData<view_type_int>(1, grid.domain_size);
                auto crit_time_step_whole_domain = collectViewData(id, np, grid, true, MPI_INT, crit_time_step);
                printViewData(id, intralayer_ofstream, grid, true, "int", "CritTimeStep", crit_time_step_whole_domain);
            }
            if (_inputs.intralayer_undercooling_change) {
                auto undercooling_change =
                    temperature.template extractTmTlCrData<view_type_float>(2, grid.domain_size, 0);
                auto undercooling_change_whole_domain =
                    collectViewData(id, np, grid, true, MPI_FLOAT, undercooling_change);
                printViewData(id, intralayer_ofstream, grid, true, "float", "UndercoolingChange",
                              undercooling_change_whole_domain);
            }
            if (_inputs.intralayer_cell_type) {
                auto cell_type_whole_domain = collectViewData(id, np, grid, true, MPI_INT, celldata.cell_type);
                printViewData(id, intralayer_ofstream, grid, true, "int", "CellType", cell_type_whole_domain);
            }
            if (_inputs.intralayer_diagonal_length) {
                auto diagonal_length_whole_domain =
                    collectViewData(id, np, grid, true, MPI_FLOAT, interface.diagonal_length);
                printViewData(id, intralayer_ofstream, grid, true, "float", "DiagonalLength",
                              diagonal_length_whole_domain);
            }
            if (_inputs.intralayer_solidification_event_counter) {
                auto solidification_event_counter_whole_domain =
                    collectViewData(id, np, grid, true, MPI_INT, temperature.solidification_event_counter);
                printViewData(id, intralayer_ofstream, grid, true, "int", "SolidificationEventCounter",
                              solidification_event_counter_whole_domain);
            }
            if (_inputs.intralayer_number_of_solidification_events) {
                auto number_of_solidification_events_whole_domain =
                    collectViewData(id, np, grid, true, MPI_INT, temperature.number_of_solidification_events);
                printViewData(id, intralayer_ofstream, grid, true, "int", "NumberOfSolidificationEvents",
                              number_of_solidification_events_whole_domain);
            }
            if (id == 0)
                intralayer_ofstream.close();

            // If necessary, collect/print GrainMisorientation field to separate file
            if (_inputs.intralayer_grain_misorientation) {
                // Get GrainID data for all layers
                auto grain_id_all_layers_whole_domain =
                    collectViewData(id, np, grid, false, MPI_INT, celldata.grain_id_all_layers);
                // Get CellType data for current layer
                auto cell_type_whole_domain = collectViewData(id, np, grid, true, MPI_INT, celldata.cell_type);
                if (id == 0) {
                    // Print GrainMisorientation for all layers up to the current layer
                    std::string misorientations_filename = vtk_filename_current_layer_base + "_Misorientations.vtk";
                    std::cout << "Printing file of grain misorientations " << misorientations_filename << std::endl;
                    printGrainMisorientations(misorientations_filename, grid, grain_id_all_layers_whole_domain,
                                              cell_type_whole_domain, orientation);
                }
            }
        }
        else
            return;
    }

    // Called on rank 0, prints current values of selected data structures to vtk files
    // Creates a file "Basefilename_layer[layernumber].vtk" (the exception being the final layer, named
    // "Basefilename.vtk") containing data for the whole simulation domain. Creates a second file
    // "Basefilename_layer[layernumber]_layeronly.vtk" (the exception being the final layer, named
    // "Basefilename_layeronly.vtk") containing data from the layer of the simulation domain that just completed. The
    // grain misorientations, if desired, will print to an additional file "Basefilename_Misorientations.vtk" for the
    // full simulation domain
    template <typename MemorySpace>
    void printInterlayer(const int id, const int np, const int layernumber, const Grid &grid,
                         CellData<MemorySpace> &celldata, Temperature<MemorySpace> &temperature,
                         Interface<MemorySpace> &interface, Orientation<MemorySpace> &orientation) {
        if (id == 0)
            std::cout << "Layer " << layernumber << " finished solidification" << std::endl;

        using view_type_float = Kokkos::View<float *, MemorySpace>;
        using view_type_int = Kokkos::View<int *, MemorySpace>;
        if (layernumber == _inputs.print_layer_number[interlayer_file_count]) {
            std::string vtk_filename_base;
            if (layernumber != grid.number_of_layers - 1)
                vtk_filename_base = path_base_filename + "_layer" + std::to_string(layernumber);
            else
                vtk_filename_base = path_base_filename;

            // Collect GrainID data for whole domain (nearly always needed for vtk file print of final output of a
            // layer)
            auto grain_id_all_layers_whole_domain =
                collectViewData(id, np, grid, false, MPI_INT, celldata.grain_id_all_layers);

            // Views where data should be printed for all layers up to and including the current one: GrainID, LayerID,
            // UndercoolingCurrent (and GrainMisorientation, but that is printed to a separate file
            if (_inputs.interlayer_full) {
                std::string vtk_filename_all_layers = vtk_filename_base + ".vtk";
                std::ofstream interlayer_all_layers_ofstream;
                if (id == 0) {
                    std::cout << "Printing data structures for all layers of the simulation to file "
                              << vtk_filename_all_layers << std::endl;
                    // Write header for printing full domain data to file
                    writeHeader(interlayer_all_layers_ofstream, vtk_filename_all_layers, grid, false);
                }
                // Get GrainID data for all layers
                if (_inputs.interlayer_grain_id)
                    printViewData(id, interlayer_all_layers_ofstream, grid, false, "int", "GrainID",
                                  grain_id_all_layers_whole_domain);
                if (_inputs.interlayer_layer_id) {
                    auto layer_id_all_layers_whole_domain =
                        collectViewData(id, np, grid, false, MPI_SHORT, celldata.layer_id_all_layers);
                    printViewData(id, interlayer_all_layers_ofstream, grid, false, "short", "LayerID",
                                  layer_id_all_layers_whole_domain);
                }
                if (_inputs.interlayer_undercooling_current) {
                    auto undercooling_all_layers_whole_domain =
                        collectViewData(id, np, grid, false, MPI_FLOAT, temperature.undercooling_current_all_layers);
                    printViewData(id, interlayer_all_layers_ofstream, grid, false, "float", "UndercoolingCurrent",
                                  undercooling_all_layers_whole_domain);
                }
                if (_inputs.interlayer_undercooling_solidification_start) {
                    auto undercooling_start_all_layers_whole_domain = collectViewData(
                        id, np, grid, false, MPI_FLOAT, temperature.undercooling_solidification_start_all_layers);
                    printViewData(id, interlayer_all_layers_ofstream, grid, false, "float", "UndercoolingStart",
                                  undercooling_start_all_layers_whole_domain);
                }
                interlayer_all_layers_ofstream.close();
            }
            // File of front undercooling at interface
            if (_inputs.print_front_undercooling) {
                auto start_end_solidification_z = temperature.getFrontUndercoolingStartFinish(id, grid);
                if (id == 0)
                    printSolidificationFrontUndercooling(grid.deltax, grid.nz, start_end_solidification_z);
            }
            // Views where data should be printed only for the layer of the problem that just finished
            if (_inputs.interlayer_current) {
                std::string vtk_filename_current_layer = vtk_filename_base + "_layeronly.vtk";
                std::ofstream currentlayer_ofstream;
                if (id == 0) {
                    std::cout << "Printing data structures for current layer of the simulation to file "
                              << vtk_filename_current_layer << std::endl;
                    // Write header for printing current layer's data to file
                    writeHeader(currentlayer_ofstream, vtk_filename_current_layer, grid, true);
                }
                if (_inputs.interlayer_melt_time_step) {
                    auto melt_time_step = temperature.template extractTmTlCrData<view_type_int>(0, grid.domain_size);
                    auto melt_time_step_whole_domain = collectViewData(id, np, grid, true, MPI_INT, melt_time_step);
                    printViewData(id, currentlayer_ofstream, grid, true, "int", "MeltTimeStep",
                                  melt_time_step_whole_domain);
                }
                if (_inputs.interlayer_crit_time_step) {
                    auto crit_time_step = temperature.template extractTmTlCrData<view_type_int>(1, grid.domain_size);
                    auto crit_time_step_whole_domain = collectViewData(id, np, grid, true, MPI_INT, crit_time_step);
                    printViewData(id, currentlayer_ofstream, grid, true, "int", "CritTimeStep",
                                  crit_time_step_whole_domain);
                }
                if (_inputs.interlayer_undercooling_change) {
                    auto undercooling_change =
                        temperature.template extractTmTlCrData<view_type_float>(2, grid.domain_size, 0);
                    auto undercooling_change_whole_domain =
                        collectViewData(id, np, grid, true, MPI_FLOAT, undercooling_change);
                    printViewData(id, currentlayer_ofstream, grid, true, "float", "UndercoolingChange",
                                  undercooling_change_whole_domain);
                }
                if (_inputs.interlayer_cell_type) {
                    auto cell_type_whole_domain = collectViewData(id, np, grid, true, MPI_INT, celldata.cell_type);
                    printViewData(id, currentlayer_ofstream, grid, true, "int", "CellType", cell_type_whole_domain);
                }
                if (_inputs.interlayer_diagonal_length) {
                    auto diagonal_length_whole_domain =
                        collectViewData(id, np, grid, true, MPI_FLOAT, interface.diagonal_length);
                    printViewData(id, currentlayer_ofstream, grid, true, "float", "DiagonalLength",
                                  diagonal_length_whole_domain);
                }
                if (_inputs.interlayer_solidification_event_counter) {
                    auto solidification_event_counter_whole_domain =
                        collectViewData(id, np, grid, true, MPI_INT, temperature.solidification_event_counter);
                    printViewData(id, currentlayer_ofstream, grid, true, "int", "SolidificationEventCounter",
                                  solidification_event_counter_whole_domain);
                }
                if (_inputs.interlayer_number_of_solidification_events) {
                    auto number_of_solidification_events_whole_domain =
                        collectViewData(id, np, grid, true, MPI_INT, temperature.number_of_solidification_events);
                    printViewData(id, currentlayer_ofstream, grid, true, "int", "NumberOfSolidificationEvents",
                                  number_of_solidification_events_whole_domain);
                }
                if (id == 0)
                    currentlayer_ofstream.close();
            }

            // If necessary, collect/print GrainMisorientation field to separate file
            if (_inputs.interlayer_grain_misorientation) {
                // Get CellType data for current layer
                auto cell_type_whole_domain = collectViewData(id, np, grid, true, MPI_INT, celldata.cell_type);
                if (id == 0) {
                    // Print GrainMisorientation for all layers up to the current layer
                    std::string misorientations_filename = vtk_filename_base + "_Misorientations.vtk";
                    std::cout << "Printing file of grain misorientations " << misorientations_filename << std::endl;
                    printGrainMisorientations(misorientations_filename, grid, grain_id_all_layers_whole_domain,
                                              cell_type_whole_domain, orientation);
                }
            }

            // If necessary, print RVE data for ExaConstit input
            if (_inputs.print_default_rve) {
                auto layer_id_all_layers_whole_domain =
                    collectViewData(id, np, grid, false, MPI_SHORT, celldata.layer_id_all_layers);
                if (id == 0)
                    printExaConstitDefaultRVE(grid, layer_id_all_layers_whole_domain, grain_id_all_layers_whole_domain);
            }
            // Update file counter
            interlayer_file_count++;
        }
        else
            return;
    }

    // Check if intermediate output is to be printed during a series of "skipped" time steps (i.e., time steps where no
    // melting or solidification occurs). If it should be printed on the time step, call
    // printIntermediateGrainMisorientation
    template <typename MemorySpace>
    void printIdleIntralayer(const int id, const int np, const int layernumber, const int deltat, const int cycle,
                             const Grid &grid, CellData<MemorySpace> &celldata, Temperature<MemorySpace> &temperature,
                             Interface<MemorySpace> &interface, Orientation<MemorySpace> &orientation,
                             const int global_next_melt_time_step) {

        if (_inputs.intralayer_idle_frames) {
            // Print current state of ExaCA simulation (up to and including the current layer's data) during the skipped
            // time steps, if intermediate output is toggled
            for (int cycle_jump = cycle + 1; cycle_jump < global_next_melt_time_step; cycle_jump++) {
                printIntralayer(id, np, layernumber, deltat, cycle_jump, grid, celldata, temperature, interface,
                                orientation);
            }
        }
    }

    // Print vtk file header data on rank 0, for either binary or ASCII output
    // current_layer_only = true: Printing only data for the Z coordinates corresponding to the current layer of a
    // multilayer problem current_layer_only = false: Printing data up to and including the top of the current layer of
    // a multilayer problem. This will be the top of the overall simulation domain if this was the final layer
    void writeHeader(std::ofstream &output_fstream, std::string filename, const Grid &grid,
                     const bool current_layer_only) {

        int z_print_size;
        float z_print_origin;
        if (current_layer_only) {
            z_print_size = grid.nz_layer;
            z_print_origin = grid.z_min + grid.z_layer_bottom * grid.deltax;
        }
        else {
            z_print_size = grid.z_layer_bottom + grid.nz_layer;
            z_print_origin = grid.z_min;
        }

        if (_inputs.print_binary)
            output_fstream.open(filename, std::ios::out | std::ios::binary);
        else
            output_fstream.open(filename);
        output_fstream << "# vtk DataFile Version 3.0" << std::endl;
        output_fstream << "vtk output" << std::endl;
        if (_inputs.print_binary)
            output_fstream << "BINARY" << std::endl;
        else
            output_fstream << "ASCII" << std::endl;
        output_fstream << "DATASET STRUCTURED_POINTS" << std::endl;
        output_fstream << "DIMENSIONS " << grid.nx << " " << grid.ny << " " << z_print_size << std::endl;
        output_fstream << "ORIGIN " << grid.x_min << " " << grid.y_min << " " << z_print_origin << std::endl;
        output_fstream << "SPACING " << grid.deltax << " " << grid.deltax << " " << grid.deltax << std::endl;
        output_fstream << std::fixed << "POINT_DATA " << grid.nx * grid.ny * z_print_size << std::endl;
    }

    // Called on rank 0 to write view data to the vtk file
    template <typename Print3DViewType>
    void printViewData(const int id, std::ofstream &output_fstream, const Grid &grid, const bool current_layer_only,
                       std::string data_label, std::string var_name_label, Print3DViewType view_data_whole_domain) {
        if (id != 0)
            return;

        // Printing the z coordinates spanning the current layer, or from the overall simulation bottom (k = 0) through
        // the top of the current layer
        int z_start = 0;
        int z_end;
        if (current_layer_only)
            z_end = grid.nz_layer;
        else
            z_end = grid.z_layer_bottom + grid.nz_layer;
        // Print data to the vtk file - casting to the appropriate type if necessary
        output_fstream << "SCALARS " << var_name_label << " " << data_label << " 1" << std::endl;
        output_fstream << "LOOKUP_TABLE default" << std::endl;
        for (int coord_z = z_start; coord_z < z_end; coord_z++) {
            for (int coord_y_global = 0; coord_y_global < grid.ny; coord_y_global++) {
                for (int coord_x = 0; coord_x < grid.nx; coord_x++) {
                    if (data_label == "int") {
                        int writeval = static_cast<int>(view_data_whole_domain(coord_z, coord_x, coord_y_global));
                        writeData(output_fstream, writeval, _inputs.print_binary, true);
                    }
                    else if (data_label == "short") {
                        short writeval = static_cast<short>(view_data_whole_domain(coord_z, coord_x, coord_y_global));
                        writeData(output_fstream, writeval, _inputs.print_binary, true);
                    }
                    else if (data_label == "float") {
                        float writeval = static_cast<float>(view_data_whole_domain(coord_z, coord_x, coord_y_global));
                        writeData(output_fstream, writeval, _inputs.print_binary, true);
                    }
                }
            }
            // Do not insert newline character if using binary writing, as this will break the binary data read by
            // adding a blank line
            if (!(_inputs.print_binary))
                output_fstream << std::endl;
        }
    }

    // Called on rank 0 to write solidification undercooling data to a file
    template <typename Print2DViewType>
    void printSolidificationFrontUndercooling(const double deltax, const int nz, Print2DViewType front_undercooling) {

        std::string front_undercooling_filename = path_base_filename + "_FrontUndercooling.csv";
        std::ofstream und_ofstream;
        und_ofstream.open(front_undercooling_filename);
        und_ofstream << "Z (micrometers), Initial Undercooling, Final Undercooling" << std::endl;
        for (int coord_z = 0; coord_z < nz - 1; coord_z++) {
            und_ofstream << static_cast<double>(coord_z) * deltax * pow(10, 6) << "," << front_undercooling(coord_z, 0)
                         << "," << front_undercooling(coord_z, 1) << std::endl;
        }
        und_ofstream << static_cast<double>(nz - 1) * deltax * pow(10, 6) << "," << front_undercooling(nz - 1, 0) << ","
                     << front_undercooling(nz - 1, 1);
        und_ofstream.close();
    }

    // On rank 0, print grain misorientation, 0-62 for epitaxial grains and 100-162 for nucleated grains, to a vtk
    // file. If printing an intermediate state, print -1 for cells that are liquid
    template <typename ViewTypeInt3DHost, typename OrientationMemory>
    void printGrainMisorientations(const std::string misorientations_filename, const Grid &grid,
                                   ViewTypeInt3DHost grain_id_whole_domain, ViewTypeInt3DHost cell_type_whole_domain,
                                   Orientation<OrientationMemory> &orientation) {

        // Print grain orientations to file - either all layers (print_region = 2), or if in an intermediate state, the
        // layers up to the current one (print_region = 1) z_end will equal grid.nz if this is the final layer
        int z_end = grid.z_layer_bottom + grid.nz_layer;
        std::ofstream misorientations_ofstream;
        writeHeader(misorientations_ofstream, misorientations_filename, grid, false);
        misorientations_ofstream << "SCALARS Angle_z short 1" << std::endl;
        misorientations_ofstream << "LOOKUP_TABLE default" << std::endl;

        // Get grain <100> misorientation relative to the Z direction for each orientation
        auto grain_misorientation = orientation.misorientationCalc(2);
        // For cells that are currently liquid (possible only for intermediate state print, as the final state will only
        // have solid cells), -1 is printed as the misorienatation. Misorientations for grains from the baseplate or
        // powder layer are between 0-62 (degrees, rounded to nearest integer), and cells that are associated with
        // nucleated grains are assigned values between 100-162 to differentiate them. Additionally, 200 is printed as
        // the misorientation for cells in the powder layer that have not been assigned a grain ID.
        // For prior layers, cell type check is unnecessary as these regions have all solidified
        for (int k = 0; k < grid.z_layer_bottom; k++) {
            for (int j = 0; j < grid.ny; j++) {
                for (int i = 0; i < grid.nx; i++) {
                    short int_print_val;
                    if (grain_id_whole_domain(k, i, j) == 0)
                        int_print_val = 200;
                    else {
                        int my_orientation =
                            getGrainOrientation(grain_id_whole_domain(k, i, j), orientation.n_grain_orientations);
                        if (grain_id_whole_domain(k, i, j) < 0)
                            int_print_val = static_cast<short>(std::round(grain_misorientation(my_orientation)) + 100);
                        else
                            int_print_val = static_cast<short>(std::round(grain_misorientation(my_orientation)));
                    }
                    writeData(misorientations_ofstream, int_print_val, _inputs.print_binary, true);
                }
            }
            if (!(_inputs.print_binary))
                misorientations_ofstream << std::endl;
        }
        // For current layer, check cell types to see if -1 should be printed (if this is a print following a layer, all
        // cells will be solid and no -1s should be written)
        for (int k = grid.z_layer_bottom; k < z_end; k++) {
            for (int j = 0; j < grid.ny; j++) {
                for (int i = 0; i < grid.nx; i++) {
                    short int_print_val;
                    if (grain_id_whole_domain(k, i, j) == 0)
                        int_print_val = 200;
                    else {
                        int my_orientation =
                            getGrainOrientation(grain_id_whole_domain(k, i, j), orientation.n_grain_orientations);
                        if (grain_id_whole_domain(k, i, j) < 0)
                            int_print_val = static_cast<short>(std::round(grain_misorientation(my_orientation)) + 100);
                        else
                            int_print_val = static_cast<short>(std::round(grain_misorientation(my_orientation)));
                    }
                    // Offset in Z as cell type values only exist for current layer
                    if (cell_type_whole_domain(k - grid.z_layer_bottom, i, j) == Liquid)
                        int_print_val = -1;
                    writeData(misorientations_ofstream, int_print_val, _inputs.print_binary, true);
                }
            }
            if (!(_inputs.print_binary))
                misorientations_ofstream << std::endl;
        }
        misorientations_ofstream.close();
    }

    // Print a representative volume element (RVE) from this multilayer simulation from the "default" location in the
    // domain. The default location is as close to the center of the domain in X and Y as possible, and as close to the
    // top of the domain while not including the final layer's microstructure. If an RVE size was not specified in the
    // input file, the default size is 0.5 by 0.5 by 0.5 mm
    template <typename ViewTypeShort3DHost, typename ViewTypeInt3DHost>
    void printExaConstitDefaultRVE(const Grid &grid, ViewTypeShort3DHost layer_id_whole_domain,
                                   ViewTypeInt3DHost grain_id_whole_domain) {

        // Determine the lower and upper Y bounds of the RVE
        int rve_xlow = std::floor(grid.nx / 2) - std::floor(_inputs.rve_size / 2);
        int rve_xhigh = rve_xlow + _inputs.rve_size - 1;
        int rve_ylow = std::floor(grid.ny / 2) - std::floor(_inputs.rve_size / 2);
        int rve_yhigh = rve_ylow + _inputs.rve_size - 1;

        // Make sure the RVE fits in the simulation domain in X and Y
        if ((rve_xlow < 0) || (rve_xhigh > grid.nx - 1) || (rve_ylow < 0) || (rve_yhigh > grid.ny - 1)) {
            std::cout << "WARNING: Simulation domain is too small to obtain default RVE data (should be at least "
                      << _inputs.rve_size << " cells in the X and Y directions" << std::endl;
            if (rve_xlow < 0)
                rve_xlow = 0;
            if (rve_xhigh > grid.nx - 1)
                rve_xhigh = grid.nx - 1;
            if (rve_ylow < 0)
                rve_ylow = 0;
            if (rve_yhigh > grid.ny - 1)
                rve_yhigh = grid.ny - 1;
        }

        // Determine the upper Z bound of the RVE - largest Z for which all cells in the RVE do not contain LayerID
        // values of the last layer
        int rve_zhigh = grid.nz - 1;
        for (int k = grid.nz - 1; k >= 0; k--) {
            [&] {
                for (int j = rve_ylow; j <= rve_yhigh; j++) {
                    for (int i = rve_xlow; i <= rve_xhigh; i++) {
                        if (layer_id_whole_domain(k, i, j) == (grid.number_of_layers - 1))
                            return; // check next lowest value for k
                    }
                }
                rve_zhigh = k;
                k = 0; // leave loop
            }();
        }

        // Determine the lower Z bound of the RVE, and make sure the RVE fits in the simulation domain in X and Y
        int rve_zlow = rve_zhigh - _inputs.rve_size + 1;
        if (rve_zlow < 0) {
            std::cout << "WARNING: Simulation domain is too small to obtain default RVE data (should be at least "
                      << _inputs.rve_size << " cells in the Z direction, more layers are required" << std::endl;
            rve_zlow = 0;
            rve_zhigh = grid.nz - 1;
        }

        // Print RVE data to file
        std::string filename = path_base_filename + "_ExaConstit.csv";
        std::cout << "Default size RVE with X coordinates " << rve_xlow << "," << rve_xhigh << "; Y coordinates "
                  << rve_ylow << "," << rve_yhigh << "; Z coordinates " << rve_zlow << "," << rve_zhigh
                  << " being printed to file " << filename << " for ExaConstit" << std::endl;
        std::ofstream exaconstit_ofstream;
        exaconstit_ofstream.open(filename);
        exaconstit_ofstream << "Coordinates are in CA units, 1 cell = " << grid.deltax
                            << " m. Data is cell-centered. Origin at " << rve_xlow << "," << rve_ylow << "," << rve_zlow
                            << " , domain size is " << rve_xhigh - rve_xlow + 1 << " by " << rve_yhigh - rve_ylow + 1
                            << " by " << rve_zhigh - rve_zlow + 1 << " cells" << std::endl;
        exaconstit_ofstream << "X coord, Y coord, Z coord, Grain ID" << std::endl;
        for (int k = rve_zlow; k <= rve_zhigh; k++) {
            for (int j = rve_ylow; j <= rve_yhigh; j++) {
                for (int i = rve_xlow; i <= rve_xhigh; i++) {
                    exaconstit_ofstream << i << "," << j << "," << k << "," << grain_id_whole_domain(k, i, j)
                                        << std::endl;
                }
            }
        }
        exaconstit_ofstream.close();
    }
};

#endif
