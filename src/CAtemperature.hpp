// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_TEMPS_HPP
#define EXACA_TEMPS_HPP

#include "CAgrid.hpp"
#include "CAinputs.hpp"
#include "CAparsefiles.hpp"
#include "mpi.h"

#include <Kokkos_Core.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Reduced form of the time-temperature history and temperature field variables used by ExaCA
template <typename MemorySpace>
struct Temperature {

    using memory_space = MemorySpace;
    using view_type_int = Kokkos::View<int *, memory_space>;
    using view_type_float = Kokkos::View<float *, memory_space>;
    using view_type_double = Kokkos::View<double *, memory_space>;
    using view_type_float_3d = Kokkos::View<float ***, memory_space>;
    using view_type_int_host = typename view_type_int::HostMirror;
    using view_type_double_host = typename view_type_double::HostMirror;
    using view_type_float_3d_host = typename view_type_float_3d::HostMirror;

    // Using the default exec space for this memory space.
    using execution_space = typename memory_space::execution_space;

    // Maximum number of times each CA cell in a given layer undergoes solidification
    view_type_int MaxSolidificationEvents;
    // For each cell in the current layer (index 1), and each time solidification happens (index 2), hold the values
    // that will be used for MeltTimeStep, CritTimeStep, and UndercoolingChange (index 3)
    view_type_float_3d LayerTimeTempHistory;
    // The number of times that each CA cell will undergo solidification during this layer
    view_type_int NumberOfSolidificationEvents;
    // A counter for the number of times each CA cell has undergone solidification so far this layer
    view_type_int SolidificationEventCounter;
    // The current undercooling of a CA cell (if superheated liquid or hasn't undergone solidification yet, equals 0)
    // Also maintained for the full multilayer domain
    view_type_float UndercoolingCurrent_AllLayers, UndercoolingCurrent;
    // Data structure for storing raw temperature data from file(s)
    // Store data as double - needed for small time steps to resolve local differences in solidification conditions
    // Each data point has 6 values (X, Y, Z coordinates, melting time, liquidus time, and cooling rate)
    view_type_double_host RawTemperatureData;
    // These contain "number_of_layers" values corresponding to the location within "RawTemperatureData" of the first
    // data element in each temperature file, if used
    view_type_int_host FirstValue, LastValue;
    // Temperature field inputs from file
    TemperatureInputs _inputs;

    // Constructor creates views with size based on the grid inputs - each cell assumed to solidify once by default,
    // LayerTimeTempHistory modified to account for multiple events if needed UndercoolingCurrent and
    // SolidificationEventCounter are default initialized to zeros
    Temperature(const Grid &grid, TemperatureInputs inputs, const int EstNumTemperatureDataPoints = 1000000)
        : MaxSolidificationEvents(
              view_type_int(Kokkos::ViewAllocateWithoutInitializing("number_of_layers"), grid.number_of_layers))
        , LayerTimeTempHistory(view_type_float_3d(Kokkos::ViewAllocateWithoutInitializing("LayerTimeTempHistory"),
                                                  grid.domain_size, 1, 3))
        , NumberOfSolidificationEvents(
              view_type_int(Kokkos::ViewAllocateWithoutInitializing("NumberOfSolidificationEvents"), grid.domain_size))
        , SolidificationEventCounter(view_type_int("SolidificationEventCounter", grid.domain_size))
        , UndercoolingCurrent_AllLayers(view_type_float("UndercoolingCurrent", grid.domain_size_all_layers))
        , RawTemperatureData(view_type_double_host(Kokkos::ViewAllocateWithoutInitializing("RawTemperatureData"),
                                                   EstNumTemperatureDataPoints))
        , FirstValue(view_type_int_host(Kokkos::ViewAllocateWithoutInitializing("FirstValue"), grid.number_of_layers))
        , LastValue(view_type_int_host(Kokkos::ViewAllocateWithoutInitializing("LastValue"), grid.number_of_layers))
        , _inputs(inputs) {

        get_current_layer_undercooling(grid.layer_range);
    }

    // Read in temperature data from files, stored in the host view "RawData", with the appropriate MPI ranks storing
    // the appropriate data
    void readTemperatureData(int id, const Grid &grid, int layernumber) {

        // Y coordinates of this rank's data, inclusive and including ghost nodes
        int LowerYBound = grid.y_offset;
        int UpperYBound = grid.y_offset + grid.ny_local - 1;

        std::cout << "On MPI rank " << id << ", the Y bounds (in cells) are [" << LowerYBound << "," << UpperYBound
                  << "]" << std::endl;
        // Store raw data relevant to each rank in the vector structure RawData
        // Two passes through reading temperature data files- this is the second pass, reading the actual X/Y/Z/liquidus
        // time/cooling rate data and each rank stores the data relevant to itself in "RawData". With remelting
        // (SimulationType == "RM"), this is the same except that some X/Y/Z coordinates may be repeated in a file, and
        // a "melting time" value is stored in addition to liquidus time and cooling rate
        int NumberOfTemperatureDataPoints = 0;
        // Second pass through the files - ignore header line
        int FirstLayerToRead, LastLayerToRead;
        if (_inputs.LayerwiseTempRead) {
            FirstLayerToRead = layernumber;
            LastLayerToRead = layernumber;
        }
        else {
            FirstLayerToRead = 0;
            LastLayerToRead = std::min(grid.number_of_layers, _inputs.TempFilesInSeries) - 1;
        }
        // Which temperature files should be read? Just the one file for layer "layernumber", or all of them?
        for (int LayerReadCount = FirstLayerToRead; LayerReadCount <= LastLayerToRead; LayerReadCount++) {

            std::string tempfile_thislayer;
            if (_inputs.LayerwiseTempRead) {
                int LayerInSeries = layernumber % _inputs.TempFilesInSeries;
                tempfile_thislayer = _inputs.temp_paths[LayerInSeries];
            }
            else
                tempfile_thislayer = _inputs.temp_paths[LayerReadCount];

            FirstValue(LayerReadCount) = NumberOfTemperatureDataPoints;
            // Read and parse temperature file for either binary or ASCII, storing the appropriate values on each MPI
            // rank within RawData and incrementing NumberOfTemperatureDataPoints appropriately
            bool BinaryInputData = checkTemperatureFileFormat(tempfile_thislayer);
            parseTemperatureData(tempfile_thislayer, grid.y_min, grid.deltax, LowerYBound, UpperYBound,
                                 NumberOfTemperatureDataPoints, BinaryInputData, RawTemperatureData);
            LastValue(LayerReadCount) = NumberOfTemperatureDataPoints;
        } // End loop over all files read for all layers
        Kokkos::resize(RawTemperatureData, NumberOfTemperatureDataPoints);
        // Determine start values for each layer's data within "RawData", if all layers were read
        if (!(_inputs.LayerwiseTempRead)) {
            if (grid.number_of_layers > _inputs.TempFilesInSeries) {
                for (int LayerReadCount = _inputs.TempFilesInSeries; LayerReadCount < grid.number_of_layers;
                     LayerReadCount++) {
                    if (_inputs.TempFilesInSeries == 1) {
                        // Since all layers have the same temperature data, each layer's "ZMinLayer" is just
                        // translated from that of the first layer
                        FirstValue(LayerReadCount) = FirstValue(LayerReadCount - 1);
                        LastValue(LayerReadCount) = LastValue(LayerReadCount - 1);
                    }
                    else {
                        // All layers have different temperature data but in a repeating pattern
                        int RepeatedFile = (LayerReadCount) % _inputs.TempFilesInSeries;
                        FirstValue(LayerReadCount) = FirstValue(RepeatedFile);
                        LastValue(LayerReadCount) = LastValue(RepeatedFile);
                    }
                }
            }
        }
    }

    // Initialize temperature data without a thermal gradient for constrained/single grain problem types
    void initialize(int id, const Grid &grid, double deltat) {

        // Initialize temperature field in Z direction with thermal gradient G set in input file
        // Liquidus front (InitUndercooling = 0) is at domain bottom for directional solidification, is at domain center
        // (with custom InitUndercooling value) for single grain solidification

        // Local copies for lambda capture.
        auto LayerTimeTempHistory_local = LayerTimeTempHistory;
        auto MaxSolidificationEvents_local = MaxSolidificationEvents;
        auto NumberOfSolidificationEvents_local = NumberOfSolidificationEvents;
        auto UndercoolingCurrent_local = UndercoolingCurrent;
        double initUndercooling_local = _inputs.initUndercooling;
        double R_local = _inputs.R;
        // Uniform undercooling field
        auto policy = Kokkos::RangePolicy<execution_space>(0, grid.domain_size);
        Kokkos::parallel_for(
            "TempInitUniform", policy, KOKKOS_LAMBDA(const int &index) {
                // All cells past melting time step and liquidus time step
                LayerTimeTempHistory_local(index, 0, 0) = -1;
                // Cells reach liquidus at a time dependent on their Z coordinate
                LayerTimeTempHistory_local(index, 0, 1) = -1;
                // Cells cool at a constant rate
                LayerTimeTempHistory_local(index, 0, 2) = R_local * deltat;
                // All cells solidify once
                MaxSolidificationEvents_local(0) = 1;
                NumberOfSolidificationEvents_local(index) = 1;
                // All cells at init undercooling
                UndercoolingCurrent_local(index) = initUndercooling_local;
            });
        if (id == 0)
            std::cout << "Undercooling field initialized to = " << _inputs.initUndercooling << " K for all cells"
                      << std::endl;
    }

    // Initialize temperature data with a thermal gradient in Z for constrained/single grain problem types
    void initialize(int id, std::string SimulationType, const Grid &grid, double deltat) {

        // Initialize temperature field in Z direction with thermal gradient G set in input file
        // Liquidus front (InitUndercooling = 0) is at domain bottom for directional solidification, is at domain center
        // (with custom InitUndercooling value) for single grain solidification
        int locationOfInitUndercooling;
        if (SimulationType == "C")
            locationOfInitUndercooling = 0;
        else
            locationOfInitUndercooling = Kokkos::floorf(static_cast<float>(grid.nz) / 2.0);
        int locationOfLiquidus =
            locationOfInitUndercooling + Kokkos::round(_inputs.initUndercooling / (_inputs.G * grid.deltax));
        // Local copies for lambda capture.
        auto LayerTimeTempHistory_local = LayerTimeTempHistory;
        auto MaxSolidificationEvents_local = MaxSolidificationEvents;
        auto NumberOfSolidificationEvents_local = NumberOfSolidificationEvents;
        auto UndercoolingCurrent_local = UndercoolingCurrent;
        double initUndercooling_local = _inputs.initUndercooling;
        double G_local = _inputs.G;
        double R_local = _inputs.R;

        auto policy = Kokkos::RangePolicy<execution_space>(0, grid.domain_size);
        Kokkos::parallel_for(
            "TempInitG", policy, KOKKOS_LAMBDA(const int &index) {
                int coord_z = grid.get_coord_Z(index);
                // Negative distFromLiquidus values for cells below the liquidus
                int distFromLiquidus = coord_z - locationOfLiquidus;
                // All cells past melting time step
                LayerTimeTempHistory_local(index, 0, 0) = -1;
                // Cells reach liquidus at a time dependent on their Z coordinate
                // Cells with negative liquidus time values are already undercooled, should have positive undercooling
                // and negative liquidus time step Cells with positive liquidus time values are not yet tracked - leave
                // current undercooling as default zeros and set liquidus time step
                if (distFromLiquidus < 0) {
                    LayerTimeTempHistory_local(index, 0, 1) = -1;
                    UndercoolingCurrent_local(index) =
                        initUndercooling_local - distFromLiquidus * (G_local * grid.deltax);
                }
                else
                    LayerTimeTempHistory_local(index, 0, 1) =
                        distFromLiquidus * G_local * grid.deltax / (R_local * deltat);
                // Cells cool at a constant rate
                LayerTimeTempHistory_local(index, 0, 2) = R_local * deltat;
                // All cells solidify once
                MaxSolidificationEvents_local(0) = 1;
                NumberOfSolidificationEvents_local(index) = 1;
            });
        if (id == 0)
            std::cout << "Temperature field initialized for unidirectional solidification with G = " << G_local
                      << " K/m" << std::endl;
    }

    // For an overlapping spot melt pattern, determine max number of times a cell will melt/solidify as part of a layer
    int calcMaxSolidificationEvents(const Grid &grid, int NumberOfSpots, int NSpotsX, int SpotRadius, int SpotOffset) {

        Kokkos::View<int **, Kokkos::HostSpace> MaxSolidificationEvents_Temp("SEvents_Temp", grid.nx, grid.ny_local);
        for (int n = 0; n < NumberOfSpots; n++) {
            int XSpotPos = SpotRadius + (n % NSpotsX) * SpotOffset;
            int YSpotPos = SpotRadius + (n / NSpotsX) * SpotOffset;
            for (int i = 0; i < grid.nx; i++) {
                float DistX = (float)(XSpotPos - i);
                for (int j = 0; j < grid.ny_local; j++) {
                    int YGlobal = j + grid.y_offset;
                    float DistY = (float)(YSpotPos - YGlobal);
                    float TotDist = Kokkos::sqrt(DistX * DistX + DistY * DistY);
                    if (TotDist <= SpotRadius) {
                        MaxSolidificationEvents_Temp(i, j)++;
                    }
                }
            }
        }
        int TempMax = 0;
        for (int i = 0; i < grid.nx; i++) {
            for (int j = 0; j < grid.ny_local; j++) {
                if (MaxSolidificationEvents_Temp(i, j) > TempMax) {
                    TempMax = MaxSolidificationEvents_Temp(i, j);
                }
            }
        }
        // Max solidification events should be the same on each rank (and each layer, since the pattern for all layers
        // is identical)
        int GlobalMaxSEvents;
        MPI_Allreduce(&TempMax, &GlobalMaxSEvents, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        return GlobalMaxSEvents;
    }

    // Initialize temperature data for an array of overlapping spot melts. As every layer is the same, this only needs
    // to be done at the start of the simulation
    // TODO: This can be performed on the device as the dirS problem is
    void initialize(int id, const Grid &grid, double FreezingRange, Inputs &inputs) {

        int NumberOfSpots = inputs.domain.NSpotsX * inputs.domain.NSpotsY;

        // Temporary host view for the maximum number of times a cell in a given layer will solidify (same for every
        // layer)
        view_type_int_host MaxSolidificationEvents_Host =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), MaxSolidificationEvents);
        int MaxNumSolidificationEvents = calcMaxSolidificationEvents(
            grid, NumberOfSpots, inputs.domain.NSpotsX, inputs.domain.SpotRadius, inputs.domain.SpotOffset);
        for (int layernumber = 0; layernumber < grid.number_of_layers; layernumber++)
            MaxSolidificationEvents_Host(layernumber) = MaxNumSolidificationEvents;

        // Resize LayerTimeTempHistory now that the max number of solidification events is known
        Kokkos::resize(LayerTimeTempHistory, grid.domain_size, MaxNumSolidificationEvents, 3);

        // Temporary host views filled with data and copied to the device. Initialize to zeros
        view_type_float_3d_host LayerTimeTempHistory_Host("TimeTempHistory_H", grid.domain_size,
                                                          MaxNumSolidificationEvents, 3);
        view_type_int_host NumberOfSolidificationEvents_Host("NumSEvents_H", grid.domain_size);

        // Outer edges of spots are initialized at the liquidus temperature
        // Spots cool at constant rate R, spot thermal gradient = G
        // Time between "start" of next spot is the time it takes for the previous spot
        // to have entirely gone below the solidus temperature
        float IsothermVelocity = (_inputs.R / _inputs.G) * inputs.domain.deltat / grid.deltax; // in cells per time step
        int TimeBetweenSpots = inputs.domain.SpotRadius / IsothermVelocity +
                               (FreezingRange / _inputs.R) / inputs.domain.deltat; // in time steps

        if (id == 0)
            std::cout << "Initializing temperature field for " << NumberOfSpots
                      << ", each of which takes approximately " << TimeBetweenSpots << " time steps to fully solidify"
                      << std::endl;

        for (int n = 0; n < NumberOfSpots; n++) {
            if (id == 0)
                std::cout << "Initializing spot " << n << std::endl;
            // Initialize LayerTimeTempHistory data values for this spot/this layer - relative to the layer bottom
            int XSpotPos = inputs.domain.SpotRadius + (n % inputs.domain.NSpotsX) * inputs.domain.SpotOffset;
            int YSpotPos = inputs.domain.SpotRadius + (n / inputs.domain.NSpotsX) * inputs.domain.SpotOffset;
            for (int coord_z = 0; coord_z <= inputs.domain.SpotRadius; coord_z++) {
                // Distance of this cell from the spot center
                float DistZ = (float)(inputs.domain.SpotRadius - coord_z);
                for (int coord_x = 0; coord_x < grid.nx; coord_x++) {
                    float DistX = (float)(XSpotPos - coord_x);
                    for (int coord_y = 0; coord_y < grid.ny_local; coord_y++) {
                        int coord_y_global = coord_y + grid.y_offset;
                        float DistY = (float)(YSpotPos - coord_y_global);
                        float TotDist = Kokkos::hypot(DistX, DistY, DistZ);
                        if (TotDist <= inputs.domain.SpotRadius) {
                            int index = grid.get_1D_index(coord_x, coord_y, coord_z);
                            // Melt time
                            LayerTimeTempHistory_Host(index, NumberOfSolidificationEvents_Host(index), 0) =
                                1 + TimeBetweenSpots * n;
                            // Liquidus time
                            int LiquidusTime = Kokkos::round((static_cast<float>(inputs.domain.SpotRadius) - TotDist) /
                                                             IsothermVelocity) +
                                               TimeBetweenSpots * n;
                            LayerTimeTempHistory_Host(index, NumberOfSolidificationEvents_Host(index), 1) =
                                1 + LiquidusTime;
                            // Cooling rate
                            LayerTimeTempHistory_Host(index, NumberOfSolidificationEvents_Host(index), 2) =
                                _inputs.R * inputs.domain.deltat;
                            NumberOfSolidificationEvents_Host(index)++;
                        }
                    }
                }
            }
        }

        // Copy host view data back to device
        MaxSolidificationEvents = Kokkos::create_mirror_view_and_copy(memory_space(), MaxSolidificationEvents_Host);
        LayerTimeTempHistory = Kokkos::create_mirror_view_and_copy(memory_space(), LayerTimeTempHistory_Host);
        NumberOfSolidificationEvents =
            Kokkos::create_mirror_view_and_copy(memory_space(), NumberOfSolidificationEvents_Host);
        MPI_Barrier(MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Spot melt temperature field initialized; each cell will solidify up to "
                      << MaxNumSolidificationEvents << " times" << std::endl;
    }

    // Calculate the number of times that a cell in layer "layernumber" undergoes melting/solidification, and store in
    // MaxSolidificationEvents_Host
    void calcMaxSolidificationEvents(int id, int layernumber, view_type_int_host MaxSolidificationEvents_Host,
                                     int StartRange, int EndRange, const Grid &grid) {

        if (layernumber > _inputs.TempFilesInSeries) {
            // Use the value from a previously checked layer, since the time-temperature history is reused
            if (_inputs.TempFilesInSeries == 1) {
                // All layers have the same temperature data, MaxSolidificationEvents for this layer is the same as the
                // last
                MaxSolidificationEvents_Host(layernumber) = MaxSolidificationEvents_Host(layernumber - 1);
            }
            else {
                // All layers have different temperature data but in a repeating pattern
                int RepeatedFile = layernumber % _inputs.TempFilesInSeries;
                MaxSolidificationEvents_Host(layernumber) = MaxSolidificationEvents_Host(RepeatedFile);
            }
        }
        else {
            // Need to calculate MaxSolidificationEvents(layernumber) from the values in RawData
            // Init to 0
            view_type_int_host TempMeltCount("TempMeltCount", grid.domain_size);

            for (int i = StartRange; i < EndRange; i += 6) {

                // Get the integer X, Y, Z coordinates associated with this data point, with the Y coordinate based on
                // local MPI rank's grid
                int coord_x = getTempCoordX(i, grid.x_min, grid.deltax);
                int coord_y = getTempCoordY(i, grid.y_min, grid.deltax, grid.y_offset);
                int coord_z = getTempCoordZ(i, grid.deltax, grid.layer_height, layernumber, grid.z_min_layer);
                // Convert to 1D coordinate in the current layer's domain
                int index = grid.get_1D_index(coord_x, coord_y, coord_z);
                TempMeltCount(index)++;
            }
            int MaxCount = 0;
            for (int i = 0; i < grid.domain_size; i++) {
                if (TempMeltCount(i) > MaxCount)
                    MaxCount = TempMeltCount(i);
            }
            int MaxCountGlobal;
            MPI_Allreduce(&MaxCount, &MaxCountGlobal, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            MaxSolidificationEvents_Host(layernumber) = MaxCountGlobal;
        }
        if (id == 0)
            std::cout << "The maximum number of melting/solidification events during layer " << layernumber << " is "
                      << MaxSolidificationEvents_Host(layernumber) << std::endl;
    }

    // Read data from storage, and calculate the normalized x value of the data point
    int getTempCoordX(int i, double XMin, double deltax) {
        int x_coord = Kokkos::round((RawTemperatureData(i) - XMin) / deltax);
        return x_coord;
    }
    // Read data from storage, and calculate the normalized y value of the data point. If the optional offset argument
    // is given, the return value is calculated relative to the edge of the MPI rank's local simulation domain (which is
    // offset by y_offset cells from the global domain edge)
    int getTempCoordY(int i, double YMin, double deltax, int y_offset = 0) {
        int y_coord = Kokkos::round((RawTemperatureData(i + 1) - YMin) / deltax) - y_offset;
        return y_coord;
    }
    // Read data from storage, and calculate the normalized z value of the data point
    int getTempCoordZ(int i, double deltax, int LayerHeight, int LayerCounter, view_type_double_host ZMinLayer) {
        int z_coord = Kokkos::round(
            (RawTemperatureData(i + 2) + deltax * LayerHeight * LayerCounter - ZMinLayer[LayerCounter]) / deltax);
        return z_coord;
    }
    // Read data from storage, obtain melting time
    double getTempCoordTM(int i) {
        double TMelting = RawTemperatureData(i + 3);
        return TMelting;
    }
    // Read data from storage, obtain liquidus time
    double getTempCoordTL(int i) {
        double TLiquidus = RawTemperatureData(i + 4);
        return TLiquidus;
    }
    // Read data from storage, obtain cooling rate
    double getTempCoordCR(int i) {
        double CoolingRate = RawTemperatureData(i + 5);
        return CoolingRate;
    }

    // Initialize temperature fields for layer "layernumber" in case where temperature data comes from file(s)
    // TODO: This can be performed on the device as the dirS problem is
    void initialize(int layernumber, int id, const Grid &grid, double FreezingRange, double deltat) {

        // Data was already read into the "RawTemperatureData" data structure
        // Determine which section of "RawTemperatureData" is relevant for this layer of the overall domain
        int StartRange = FirstValue[layernumber];
        int EndRange = LastValue[layernumber];

        // Temporary host view for the maximum number of times a cell in a given layer will solidify (different for each
        // layer)
        view_type_int_host MaxSolidificationEvents_Host =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), MaxSolidificationEvents);
        calcMaxSolidificationEvents(id, layernumber, MaxSolidificationEvents_Host, StartRange, EndRange, grid);
        int MaxNumSolidificationEvents = MaxSolidificationEvents_Host(layernumber);

        // Resize LayerTimeTempHistory now that the max number of solidification events is known for this layer
        Kokkos::resize(LayerTimeTempHistory, grid.domain_size, MaxNumSolidificationEvents, 3);

        // These views are initialized to zeros on the host, filled with data, and then copied to the device for layer
        // "layernumber"
        view_type_float_3d_host LayerTimeTempHistory_Host("TimeTempHistory_H", grid.domain_size,
                                                          MaxNumSolidificationEvents, 3);
        view_type_int_host NumberOfSolidificationEvents_Host("NumSEvents_H", grid.domain_size);

        double LargestTime = 0;
        double LargestTime_Global = 0;
        if (id == 0)
            std::cout << "Range of raw data for layer " << layernumber << " on rank 0 is " << StartRange << " to "
                      << EndRange << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        for (int i = StartRange; i < EndRange; i += 6) {

            // Get the integer X, Y, Z coordinates associated with this data point, along with the associated TM, TL, CR
            // values
            // coord_y is relative to ths MPI rank's grid, while coord_y_global is relative to the overall simulation
            // domain
            int coord_x = getTempCoordX(i, grid.x_min, grid.deltax);
            int coord_y = getTempCoordY(i, grid.y_min, grid.deltax, grid.y_offset);
            int coord_z = getTempCoordZ(i, grid.deltax, grid.layer_height, layernumber, grid.z_min_layer);
            double TMelting = getTempCoordTM(i);
            double TLiquidus = getTempCoordTL(i);
            double CoolingRate = getTempCoordCR(i);

            // 1D cell coordinate on this MPI rank's domain
            int index = grid.get_1D_index(coord_x, coord_y, coord_z);
            // Store TM, TL, CR values for this solidification event in LayerTimeTempHistory
            LayerTimeTempHistory_Host(index, NumberOfSolidificationEvents_Host(index), 0) =
                Kokkos::round(TMelting / deltat) + 1;
            LayerTimeTempHistory_Host(index, NumberOfSolidificationEvents_Host(index), 1) =
                Kokkos::round(TLiquidus / deltat) + 1;
            LayerTimeTempHistory_Host(index, NumberOfSolidificationEvents_Host(index), 2) =
                std::abs(CoolingRate) * deltat;
            // Increment number of solidification events for this cell
            NumberOfSolidificationEvents_Host(index)++;
            // Estimate of the time step where the last possible solidification is expected to occur
            double SolidusTime = TLiquidus + FreezingRange / CoolingRate;
            if (SolidusTime > LargestTime)
                LargestTime = SolidusTime;
        }
        MPI_Allreduce(&LargestTime, &LargestTime_Global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (id == 0)
            std::cout << " Layer " << layernumber
                      << " time step where all cells are cooled below solidus for the final time is "
                      << Kokkos::round((LargestTime_Global) / deltat) << " (or " << LargestTime_Global << " seconds)"
                      << std::endl;
        if (id == 0)
            std::cout << "Layer " << layernumber << " temperatures read" << std::endl;

        // Reorder solidification events in LayerTimeTempHistory(location,event number,component) so that they are in
        // order based on the melting time values (component = 0)
        for (int index = 0; index < grid.domain_size; index++) {
            int NSolidificationEvents_cell = NumberOfSolidificationEvents_Host(index);
            if (NSolidificationEvents_cell > 0) {
                for (int i = 0; i < NSolidificationEvents_cell - 1; i++) {
                    for (int j = (i + 1); j < NSolidificationEvents_cell; j++) {
                        if (LayerTimeTempHistory_Host(index, i, 0) > LayerTimeTempHistory_Host(index, j, 0)) {
                            // Swap these two points - melting event "j" happens before event "i"
                            float OldMeltVal = LayerTimeTempHistory_Host(index, i, 0);
                            float OldLiqVal = LayerTimeTempHistory_Host(index, i, 1);
                            float OldCRVal = LayerTimeTempHistory_Host(index, i, 2);
                            LayerTimeTempHistory_Host(index, i, 0) = LayerTimeTempHistory_Host(index, j, 0);
                            LayerTimeTempHistory_Host(index, i, 1) = LayerTimeTempHistory_Host(index, j, 1);
                            LayerTimeTempHistory_Host(index, i, 2) = LayerTimeTempHistory_Host(index, j, 2);
                            LayerTimeTempHistory_Host(index, j, 0) = OldMeltVal;
                            LayerTimeTempHistory_Host(index, j, 1) = OldLiqVal;
                            LayerTimeTempHistory_Host(index, j, 2) = OldCRVal;
                        }
                    }
                }
            }
        }
        // If a cell melts twice before reaching the liquidus temperature, this is a double counted solidification
        // event and should be removed
        for (int index = 0; index < grid.domain_size; index++) {
            int NSolidificationEvents_cell = NumberOfSolidificationEvents_Host(index);
            if (NSolidificationEvents_cell > 1) {
                for (int i = 0; i < NSolidificationEvents_cell - 1; i++) {
                    if (LayerTimeTempHistory_Host(index, i + 1, 0) < LayerTimeTempHistory_Host(index, i, 1)) {
                        std::cout << "Cell " << index << " removing anomalous event " << i + 1 << " out of "
                                  << NSolidificationEvents_cell - 1 << std::endl;
                        // Keep whichever event has the larger liquidus time
                        if (LayerTimeTempHistory_Host(index, i + 1, 1) > LayerTimeTempHistory_Host(index, i, 1)) {
                            LayerTimeTempHistory_Host(index, i, 0) = LayerTimeTempHistory_Host(index, i + 1, 0);
                            LayerTimeTempHistory_Host(index, i, 1) = LayerTimeTempHistory_Host(index, i + 1, 1);
                            LayerTimeTempHistory_Host(index, i, 2) = LayerTimeTempHistory_Host(index, i + 1, 2);
                        }
                        LayerTimeTempHistory_Host(index, i + 1, 0) = 0.0;
                        LayerTimeTempHistory_Host(index, i + 1, 1) = 0.0;
                        LayerTimeTempHistory_Host(index, i + 1, 2) = 0.0;
                        // Reshuffle other solidification events over if needed
                        for (int ii = (i + 1); ii < NSolidificationEvents_cell - 1; ii++) {
                            LayerTimeTempHistory_Host(index, ii, 0) = LayerTimeTempHistory_Host(index, ii + 1, 0);
                            LayerTimeTempHistory_Host(index, ii, 1) = LayerTimeTempHistory_Host(index, ii + 1, 1);
                            LayerTimeTempHistory_Host(index, ii, 2) = LayerTimeTempHistory_Host(index, ii + 1, 2);
                        }
                        NumberOfSolidificationEvents_Host(index)--;
                    }
                }
            }
        }

        // Copy host view data back to device
        MaxSolidificationEvents = Kokkos::create_mirror_view_and_copy(memory_space(), MaxSolidificationEvents_Host);
        LayerTimeTempHistory = Kokkos::create_mirror_view_and_copy(memory_space(), LayerTimeTempHistory_Host);
        NumberOfSolidificationEvents =
            Kokkos::create_mirror_view_and_copy(memory_space(), NumberOfSolidificationEvents_Host);

        if (id == 0)
            std::cout << "Layer " << layernumber << " temperature field is from Z = " << grid.z_layer_bottom
                      << " through " << grid.nz_layer + grid.z_layer_bottom - 1 << " of the global domain" << std::endl;
    }

    // Get the subview associated with the undercooling of cells in the current layer. Do not reset the undercooling of
    // cells from the prior layer to zero as this information will be stored for a potential print (and a cell that
    // remelts in the current layer will have its undercooling reset to 0 and recalculated)
    void get_current_layer_undercooling(std::pair<int, int> LayerRange) {
        UndercoolingCurrent = Kokkos::subview(UndercoolingCurrent_AllLayers, LayerRange);
    }

    // Reset local cell undercooling to 0
    KOKKOS_INLINE_FUNCTION
    void reset_undercooling(const int index) const { UndercoolingCurrent(index) = 0.0; }

    // Update local cell undercooling for the current melt-resolidification event
    KOKKOS_INLINE_FUNCTION
    void update_undercooling(const int index) const {
        UndercoolingCurrent(index) += LayerTimeTempHistory(index, SolidificationEventCounter(index), 2);
    }

    // Update the solidification event counter for a cell that has not finished the previous solidification event (i.e.,
    // solidification is not complete in this cell as it is not tempsolid or solid type)
    KOKKOS_INLINE_FUNCTION
    void update_solidification_counter(const int index) const { SolidificationEventCounter(index)++; }

    // Update the solidification event counter for the cell (which is either tempsolid or solid type) and return whether
    // all solidification events have completed in the cell
    KOKKOS_INLINE_FUNCTION
    bool update_check_solidification_counter(const int index) const {
        bool solidification_complete;
        SolidificationEventCounter(index)++;
        if (SolidificationEventCounter(index) == NumberOfSolidificationEvents(index))
            solidification_complete = true;
        else
            solidification_complete = false;
        return solidification_complete;
    }

    // Reset solidification event counter and get the subview associated with the undercooling field for the next layer
    void reset_layer_events_undercooling(const Grid &grid) {
        get_current_layer_undercooling(grid.layer_range);
        Kokkos::realloc(SolidificationEventCounter, grid.domain_size);
        Kokkos::deep_copy(SolidificationEventCounter, 0);
    }

    // Extract the next time that this point undergoes melting
    KOKKOS_INLINE_FUNCTION
    int getMeltTimeStep(const int cycle, const int index) const {
        double MeltTimeStep;
        int SolidificationEventCounter_cell = SolidificationEventCounter(index);
        MeltTimeStep = static_cast<int>(LayerTimeTempHistory(index, SolidificationEventCounter_cell, 0));
        if (cycle > MeltTimeStep) {
            // If the cell has already exceeded the melt time step for the current melt-solidification event, get the
            // melt time step associated with the next solidification event - or, if there is no next
            // melt-solidification event, return the max possible int as the cell will not melt again during this layer
            // of the multilayer problem
            if (SolidificationEventCounter_cell < (NumberOfSolidificationEvents(index) - 1))
                MeltTimeStep = static_cast<int>(LayerTimeTempHistory(index, SolidificationEventCounter_cell + 1, 0));
            else
                MeltTimeStep = INT_MAX;
        }
        return MeltTimeStep;
    }

    // Extract the next time that this point cools below the liquidus
    // Uses the current value of the solidification event counter
    KOKKOS_INLINE_FUNCTION
    int getCritTimeStep(const int index) const {
        int SolidificationEventCounter_cell = SolidificationEventCounter(index);
        int CritTimeStep = static_cast<int>(LayerTimeTempHistory(index, SolidificationEventCounter_cell, 1));
        return CritTimeStep;
    }
    // Uses a specified solidification event
    KOKKOS_INLINE_FUNCTION
    int getCritTimeStep(const int index, const int SolidificationEventCounter_cell) const {
        int CritTimeStep = static_cast<int>(LayerTimeTempHistory(index, SolidificationEventCounter_cell, 1));
        return CritTimeStep;
    }

    // Extract the cooling rate associated with a specified solidificaiton event
    KOKKOS_INLINE_FUNCTION
    float getUndercoolingChange(const int index, const int SolidificationEventCounter_cell) const {
        float UndercoolingChange = LayerTimeTempHistory(index, SolidificationEventCounter_cell, 2);
        return UndercoolingChange;
    }

    // Extract either the last time step that all points undergo melting in the layer, the last time they cools below
    // the liquidus, or the rate at which they cools from the liquidus from LayerTimeTempHistory (corresponds to
    // solidification event number `NumSolidificationEvents-1` for the cell) (can't just use subview here since
    // NumSolidificationEvents is different for each cell) If the cell does not undergo solidification, either print -1
    // or the specified default value
    template <typename ExtractedViewDataType>
    ExtractedViewDataType extract_tm_tl_cr_data(const int extracted_val, const int domain_size,
                                                const int DefaultVal = -1) {
        ExtractedViewDataType ExtractedData(Kokkos::ViewAllocateWithoutInitializing("ExtractedData"), domain_size);
        using extracted_value_type = typename ExtractedViewDataType::value_type;

        // Local copies for lambda capture.
        auto LayerTimeTempHistory_local = LayerTimeTempHistory;
        auto NumberOfSolidificationEvents_local = NumberOfSolidificationEvents;

        auto policy = Kokkos::RangePolicy<execution_space>(0, domain_size);
        Kokkos::parallel_for(
            "Extract_tm_tl_cr_data", policy, KOKKOS_LAMBDA(const int &index) {
                int NumSolidificationEvents_ThisCell = NumberOfSolidificationEvents_local(index);
                // If this cell doesn't undergo solidification at all, print -1
                if (NumSolidificationEvents_ThisCell == 0)
                    ExtractedData(index) = static_cast<extracted_value_type>(DefaultVal);
                else
                    ExtractedData(index) = static_cast<extracted_value_type>(
                        LayerTimeTempHistory_local(index, NumSolidificationEvents_ThisCell - 1, extracted_val));
            });
        return ExtractedData;
    }
};

#endif
