// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_TEMPS_HPP
#define EXACA_TEMPS_HPP

#include "CAconfig.hpp"
#include "CAparsefiles.hpp"
#include "CAtypes.hpp"
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
    view_type_float UndercoolingCurrent;
    // Data structure for storing raw temperature data from file(s)
    // Store data as double - needed for small time steps to resolve local differences in solidification conditions
    // Each data point has 6 values (X, Y, Z coordinates, melting time, liquidus time, and cooling rate)
    view_type_double_host RawTemperatureData;
    // These contain "NumberOfLayers" values corresponding to the location within "RawTemperatureData" of the first data
    // element in each temperature file, if used
    view_type_int_host FirstValue, LastValue;

    // Constructor creates views with size based on the grid inputs - each cell assumed to solidify once by default,
    // LayerTimeTempHistory modified to account for multiple events if needed UndercoolingCurrent and
    // SolidificationEventCounter are default initialized to zeros
    Temperature(const int DomainSize, const int NumberOfLayers, const int EstNumTemperatureDataPoints = 1000000)
        : MaxSolidificationEvents(
              view_type_int(Kokkos::ViewAllocateWithoutInitializing("NumberOfLayers"), NumberOfLayers))
        , LayerTimeTempHistory(
              view_type_float_3d(Kokkos::ViewAllocateWithoutInitializing("LayerTimeTempHistory"), DomainSize, 1, 3))
        , NumberOfSolidificationEvents(
              view_type_int(Kokkos::ViewAllocateWithoutInitializing("NumberOfSolidificationEvents"), DomainSize))
        , SolidificationEventCounter(view_type_int("SolidificationEventCounter", DomainSize))
        , UndercoolingCurrent(view_type_float("UndercoolingCurrent", DomainSize))
        , RawTemperatureData(view_type_double_host(Kokkos::ViewAllocateWithoutInitializing("RawTemperatureData"),
                                                   EstNumTemperatureDataPoints))
        , FirstValue(view_type_int_host(Kokkos::ViewAllocateWithoutInitializing("FirstValue"), NumberOfLayers))
        , LastValue(view_type_int_host(Kokkos::ViewAllocateWithoutInitializing("LastValue"), NumberOfLayers)) {}

    // Check if the temperature data is in ASCII or binary format
    bool checkTemperatureFileFormat(std::string tempfile_thislayer) {
        bool BinaryInputData;
        std::size_t found = tempfile_thislayer.find(".catemp");
        if (found == std::string::npos)
            BinaryInputData = false;
        else
            BinaryInputData = true;
        return BinaryInputData;
    };

    // Read in temperature data from files, stored in the host view "RawData", with the appropriate MPI ranks storing
    // the appropriate data
    void readTemperatureData(int id, double &deltax, double HT_deltax, int &HTtoCAratio, int y_offset, int ny_local,
                             double YMin, std::vector<std::string> &temp_paths, int NumberOfLayers,
                             int TempFilesInSeries, bool LayerwiseTempRead, int layernumber) {

        double HTtoCAratio_unrounded = HT_deltax / deltax;
        double HTtoCAratio_floor = floor(HTtoCAratio_unrounded);
        if (((HTtoCAratio_unrounded - HTtoCAratio_floor) > 0.0005) && (id == 0)) {
            std::string error = "Error: Temperature data point spacing not evenly divisible by CA cell size";
            throw std::runtime_error(error);
        }
        else if (((HTtoCAratio_unrounded - HTtoCAratio_floor) > 0.000001) && (id == 0)) {
            std::cout << "Note: Adjusting cell size from " << deltax << " to " << HT_deltax / HTtoCAratio_floor
                      << " to "
                         "ensure even divisibility of CA cell size into temperature data spacing"
                      << std::endl;
        }
        // Adjust deltax to exact value based on temperature data spacing and ratio between heat transport/CA cell sizes
        deltax = HT_deltax / HTtoCAratio_floor;
        HTtoCAratio = round(HT_deltax / deltax); // OpenFOAM/CA cell size ratio
        // If HTtoCAratio > 1, an interpolation of input temperature data is needed
        // The Y bounds are the region (for this MPI rank) of the physical domain that needs to be
        // read extends past the actual spatial extent of the local domain for purposes of interpolating
        // from HT_deltax to deltax
        int LowerYBound = y_offset - (y_offset % HTtoCAratio);
        int UpperYBound;
        if (HTtoCAratio == 1)
            UpperYBound = y_offset + ny_local - 1;
        else
            UpperYBound = y_offset + ny_local - 1 + HTtoCAratio - (y_offset + ny_local - 1) % HTtoCAratio;

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
        if (LayerwiseTempRead) {
            FirstLayerToRead = layernumber;
            LastLayerToRead = layernumber;
        }
        else {
            FirstLayerToRead = 0;
            LastLayerToRead = std::min(NumberOfLayers, TempFilesInSeries) - 1;
        }
        // Which temperature files should be read? Just the one file for layer "layernumber", or all of them?
        for (int LayerReadCount = FirstLayerToRead; LayerReadCount <= LastLayerToRead; LayerReadCount++) {

            std::string tempfile_thislayer;
            if (LayerwiseTempRead) {
                int LayerInSeries = layernumber % TempFilesInSeries;
                tempfile_thislayer = temp_paths[LayerInSeries];
            }
            else
                tempfile_thislayer = temp_paths[LayerReadCount];

            FirstValue(LayerReadCount) = NumberOfTemperatureDataPoints;
            // Read and parse temperature file for either binary or ASCII, storing the appropriate values on each MPI
            // rank within RawData and incrementing NumberOfTemperatureDataPoints appropriately
            bool BinaryInputData = checkTemperatureFileFormat(tempfile_thislayer);
            parseTemperatureData(tempfile_thislayer, YMin, deltax, LowerYBound, UpperYBound,
                                 NumberOfTemperatureDataPoints, BinaryInputData, RawTemperatureData);
            LastValue(LayerReadCount) = NumberOfTemperatureDataPoints;
        } // End loop over all files read for all layers
        Kokkos::resize(RawTemperatureData, NumberOfTemperatureDataPoints);
        // Determine start values for each layer's data within "RawData", if all layers were read
        if (!(LayerwiseTempRead)) {
            if (NumberOfLayers > TempFilesInSeries) {
                for (int LayerReadCount = TempFilesInSeries; LayerReadCount < NumberOfLayers; LayerReadCount++) {
                    if (TempFilesInSeries == 1) {
                        // Since all layers have the same temperature data, each layer's "ZMinLayer" is just
                        // translated from that of the first layer
                        FirstValue(LayerReadCount) = FirstValue(LayerReadCount - 1);
                        LastValue(LayerReadCount) = LastValue(LayerReadCount - 1);
                    }
                    else {
                        // All layers have different temperature data but in a repeating pattern
                        int RepeatedFile = (LayerReadCount) % TempFilesInSeries;
                        FirstValue(LayerReadCount) = FirstValue(RepeatedFile);
                        LastValue(LayerReadCount) = LastValue(RepeatedFile);
                    }
                }
            }
        }
    }

    // Initialize temperature data without a thermal gradient for constrained/single grain problem types
    void initialize(double R, int id, double deltat, int DomainSize, double initUndercooling) {

        // Initialize temperature field in Z direction with thermal gradient G set in input file
        // Liquidus front (InitUndercooling = 0) is at domain bottom for directional solidification, is at domain center
        // (with custom InitUndercooling value) for single grain solidification
        auto LayerTimeTempHistory_local = LayerTimeTempHistory;
        auto MaxSolidificationEvents_local = MaxSolidificationEvents;
        auto NumberOfSolidificationEvents_local = NumberOfSolidificationEvents;
        auto UndercoolingCurrent_local = UndercoolingCurrent;

        // Uniform undercooling field
        Kokkos::parallel_for(
            "TempInitUniform", DomainSize, KOKKOS_LAMBDA(const int &index) {
                // All cells past melting time step and liquidus time step
                LayerTimeTempHistory_local(index, 0, 0) = -1;
                // Cells reach liquidus at a time dependent on their Z coordinate
                LayerTimeTempHistory_local(index, 0, 1) = -1;
                // Cells cool at a constant rate
                LayerTimeTempHistory_local(index, 0, 2) = R * deltat;
                // All cells solidify once
                MaxSolidificationEvents_local(0) = 1;
                NumberOfSolidificationEvents_local(index) = 1;
                // All cells at init undercooling
                UndercoolingCurrent(index) = initUndercooling;
            });
        if (id == 0)
            std::cout << "Undercooling field initialized to = " << initUndercooling << " K for all cells" << std::endl;
    }

    // Initialize temperature data with a thermal gradient in Z for constrained/single grain problem types
    void initialize(std::string SimulationType, double G, double R, int id, int nx, int ny_local, int nz, double deltax,
                    double deltat, int DomainSize, double initUndercooling) {

        // Initialize temperature field in Z direction with thermal gradient G set in input file
        // Liquidus front (InitUndercooling = 0) is at domain bottom for directional solidification, is at domain center
        // (with custom InitUndercooling value) for single grain solidification
        int locationOfInitUndercooling;
        if (SimulationType == "C")
            locationOfInitUndercooling = 0;
        else
            locationOfInitUndercooling = floorf(static_cast<float>(nz) / 2.0);
        int locationOfLiquidus = locationOfInitUndercooling + round(initUndercooling / (G * deltax));

        auto LayerTimeTempHistory_local = LayerTimeTempHistory;
        auto MaxSolidificationEvents_local = MaxSolidificationEvents;
        auto NumberOfSolidificationEvents_local = NumberOfSolidificationEvents;
        auto UndercoolingCurrent_local = UndercoolingCurrent;
        Kokkos::parallel_for(
            "TempInitG", DomainSize, KOKKOS_LAMBDA(const int &index) {
                int coord_z = getCoordZ(index, nx, ny_local);
                // Negative distFromLiquidus values for cells below the liquidus
                int distFromLiquidus = coord_z - locationOfLiquidus;
                // All cells past melting time step
                LayerTimeTempHistory_local(index, 0, 0) = -1;
                // Cells reach liquidus at a time dependent on their Z coordinate
                float liquidusTime = distFromLiquidus * G * deltax / (R * deltat);
                // Cells with negative liquidus time values are already undercooled, should have positive undercooling
                // Cells with positive liquidus time values are not yet tracked
                // Leave current undercooling as default zeros
                if (liquidusTime < 0)
                    LayerTimeTempHistory_local(index, 0, 1) = -1;
                else
                    LayerTimeTempHistory_local(index, 0, 1) = liquidusTime;
                // Cells cool at a constant rate
                LayerTimeTempHistory_local(index, 0, 2) = R * deltat;
                // All cells solidify once
                MaxSolidificationEvents_local(0) = 1;
                NumberOfSolidificationEvents_local(index) = 1;
            });
        if (id == 0)
            std::cout << "Temperature field initialized for unidirectional solidification with G = " << G << " K/m"
                      << std::endl;
    }

    // For an overlapping spot melt pattern, determine max number of times a cell will melt/solidify as part of a layer
    int calcMaxSolidificationEvents(int nx, int ny_local, int NumberOfSpots, int NSpotsX, int SpotRadius,
                                    int SpotOffset, int y_offset) {

        ViewI2D_H MaxSolidificationEvents_Temp("SEvents_Temp", nx, ny_local);
        for (int n = 0; n < NumberOfSpots; n++) {
            int XSpotPos = SpotRadius + (n % NSpotsX) * SpotOffset;
            int YSpotPos = SpotRadius + (n / NSpotsX) * SpotOffset;
            for (int i = 0; i < nx; i++) {
                float DistX = (float)(XSpotPos - i);
                for (int j = 0; j < ny_local; j++) {
                    int YGlobal = j + y_offset;
                    float DistY = (float)(YSpotPos - YGlobal);
                    float TotDist = sqrt(DistX * DistX + DistY * DistY);
                    if (TotDist <= SpotRadius) {
                        MaxSolidificationEvents_Temp(i, j)++;
                    }
                }
            }
        }
        int TempMax = 0;
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny_local; j++) {
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
    void initialize(double G, double R, int id, int nx, int ny_local, int y_offset, double deltax, double deltat,
                    int DomainSize, double FreezingRange, int NSpotsX, int NSpotsY, int SpotRadius, int SpotOffset,
                    int NumberOfLayers) {

        int NumberOfSpots = NSpotsX * NSpotsY;

        // Temporary host view for the maximum number of times a cell in a given layer will solidify (same for every
        // layer)
        view_type_int_host MaxSolidificationEvents_Host =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), MaxSolidificationEvents);
        int MaxNumSolidificationEvents =
            calcMaxSolidificationEvents(nx, ny_local, NumberOfSpots, NSpotsX, SpotRadius, SpotOffset, y_offset);
        for (int layernumber = 0; layernumber < NumberOfLayers; layernumber++)
            MaxSolidificationEvents_Host(layernumber) = MaxNumSolidificationEvents;

        // Resize LayerTimeTempHistory now that the max number of solidification events is known
        Kokkos::resize(LayerTimeTempHistory, DomainSize, MaxNumSolidificationEvents, 3);

        // Temporary host views filled with data and copied to the device. Initialize to zeros
        view_type_float_3d_host LayerTimeTempHistory_Host("TimeTempHistory_H", DomainSize, MaxNumSolidificationEvents,
                                                          3);
        view_type_int_host NumberOfSolidificationEvents_Host("NumSEvents_H", DomainSize);

        // Outer edges of spots are initialized at the liquidus temperature
        // Spots cool at constant rate R, spot thermal gradient = G
        // Time between "start" of next spot is the time it takes for the previous spot
        // to have entirely gone below the solidus temperature
        float IsothermVelocity = (R / G) * deltat / deltax;                                  // in cells per time step
        int TimeBetweenSpots = SpotRadius / IsothermVelocity + (FreezingRange / R) / deltat; // in time steps

        if (id == 0)
            std::cout << "Initializing temperature field for " << NumberOfSpots
                      << ", each of which takes approximately " << TimeBetweenSpots << " time steps to fully solidify"
                      << std::endl;

        for (int n = 0; n < NumberOfSpots; n++) {
            if (id == 0)
                std::cout << "Initializing spot " << n << std::endl;
            // Initialize LayerTimeTempHistory data values for this spot/this layer - relative to the layer bottom
            int XSpotPos = SpotRadius + (n % NSpotsX) * SpotOffset;
            int YSpotPos = SpotRadius + (n / NSpotsX) * SpotOffset;
            for (int coord_z = 0; coord_z <= SpotRadius; coord_z++) {
                // Distance of this cell from the spot center
                float DistZ = (float)(SpotRadius - coord_z);
                for (int coord_x = 0; coord_x < nx; coord_x++) {
                    float DistX = (float)(XSpotPos - coord_x);
                    for (int coord_y = 0; coord_y < ny_local; coord_y++) {
                        int coord_y_global = coord_y + y_offset;
                        float DistY = (float)(YSpotPos - coord_y_global);
                        float TotDist = sqrt(DistX * DistX + DistY * DistY + DistZ * DistZ);
                        if (TotDist <= SpotRadius) {
                            int index = get1Dindex(coord_x, coord_y, coord_z, nx, ny_local);
                            // Melt time
                            LayerTimeTempHistory_Host(index, NumberOfSolidificationEvents_Host(index), 0) =
                                1 + TimeBetweenSpots * n;
                            // Liquidus time
                            int LiquidusTime = round((static_cast<float>(SpotRadius) - TotDist) / IsothermVelocity) +
                                               TimeBetweenSpots * n;
                            LayerTimeTempHistory_Host(index, NumberOfSolidificationEvents_Host(index), 1) =
                                1 + LiquidusTime;
                            // Cooling rate
                            LayerTimeTempHistory_Host(index, NumberOfSolidificationEvents_Host(index), 2) = R * deltat;
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
    void calcMaxSolidificationEvents(int id, int layernumber, int TempFilesInSeries,
                                     ViewI_H MaxSolidificationEvents_Host, int StartRange, int EndRange, double XMin,
                                     double YMin, double deltax, double *ZMinLayer, int LayerHeight, int nx,
                                     int ny_local, int y_offset, int DomainSize) {

        if (layernumber > TempFilesInSeries) {
            // Use the value from a previously checked layer, since the time-temperature history is reused
            if (TempFilesInSeries == 1) {
                // All layers have the same temperature data, MaxSolidificationEvents for this layer is the same as the
                // last
                MaxSolidificationEvents_Host(layernumber) = MaxSolidificationEvents_Host(layernumber - 1);
            }
            else {
                // All layers have different temperature data but in a repeating pattern
                int RepeatedFile = layernumber % TempFilesInSeries;
                MaxSolidificationEvents_Host(layernumber) = MaxSolidificationEvents_Host(RepeatedFile);
            }
        }
        else {
            // Need to calculate MaxSolidificationEvents(layernumber) from the values in RawData
            // Init to 0
            view_type_int_host TempMeltCount("TempMeltCount", DomainSize);

            for (int i = StartRange; i < EndRange; i += 6) {

                // Get the integer X, Y, Z coordinates associated with this data point, with the Y coordinate based on
                // local MPI rank's grid
                int coord_x = getTempCoordX(i, XMin, deltax);
                int coord_y = getTempCoordY(i, YMin, deltax, y_offset);
                int coord_z = getTempCoordZ(i, deltax, LayerHeight, layernumber, ZMinLayer);
                // Convert to 1D coordinate in the current layer's domain
                int index = get1Dindex(coord_x, coord_y, coord_z, nx, ny_local);
                TempMeltCount(index)++;
            }
            int MaxCount = 0;
            for (int i = 0; i < DomainSize; i++) {
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
        int x_coord = round((RawTemperatureData(i) - XMin) / deltax);
        return x_coord;
    }
    // Read data from storage, and calculate the normalized y value of the data point. If the optional offset argument
    // is given, the return value is calculated relative to the edge of the MPI rank's local simulation domain (which is
    // offset by y_offset cells from the global domain edge)
    int getTempCoordY(int i, double YMin, double deltax, int y_offset = 0) {
        int y_coord = round((RawTemperatureData(i + 1) - YMin) / deltax) - y_offset;
        return y_coord;
    }
    // Read data from storage, and calculate the normalized z value of the data point
    int getTempCoordZ(int i, double deltax, int LayerHeight, int LayerCounter, double *ZMinLayer) {
        int z_coord =
            round((RawTemperatureData(i + 2) + deltax * LayerHeight * LayerCounter - ZMinLayer[LayerCounter]) / deltax);
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
    void initialize(int layernumber, int id, int nx, int ny_local, int DomainSize, int y_offset, double deltax,
                    double deltat, double FreezingRange, double XMin, double YMin, double *ZMinLayer, int LayerHeight,
                    int nz_layer, int z_layer_bottom, int *FinishTimeStep, int TempFilesInSeries) {

        // Data was already read into the "RawTemperatureData" data structure
        // Determine which section of "RawTemperatureData" is relevant for this layer of the overall domain
        int StartRange = FirstValue[layernumber];
        int EndRange = LastValue[layernumber];

        // Temporary host view for the maximum number of times a cell in a given layer will solidify (different for each
        // layer)
        view_type_int_host MaxSolidificationEvents_Host =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), MaxSolidificationEvents);
        calcMaxSolidificationEvents(id, layernumber, TempFilesInSeries, MaxSolidificationEvents_Host, StartRange,
                                    EndRange, XMin, YMin, deltax, ZMinLayer, LayerHeight, nx, ny_local, y_offset,
                                    DomainSize);
        int MaxNumSolidificationEvents = MaxSolidificationEvents_Host(0);

        // Resize LayerTimeTempHistory now that the max number of solidification events is known for this layer
        Kokkos::resize(LayerTimeTempHistory, DomainSize, MaxNumSolidificationEvents, 3);

        // These views are initialized to zeros on the host, filled with data, and then copied to the device for layer
        // "layernumber"
        view_type_float_3d_host LayerTimeTempHistory_Host("TimeTempHistory_H", DomainSize, MaxNumSolidificationEvents,
                                                          3);
        view_type_int_host NumberOfSolidificationEvents_Host("NumSEvents_H", DomainSize);

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
            int coord_x = getTempCoordX(i, XMin, deltax);
            int coord_y = getTempCoordY(i, YMin, deltax, y_offset);
            int coord_z = getTempCoordZ(i, deltax, LayerHeight, layernumber, ZMinLayer);
            double TMelting = getTempCoordTM(i);
            double TLiquidus = getTempCoordTL(i);
            double CoolingRate = getTempCoordCR(i);

            // 1D cell coordinate on this MPI rank's domain
            int index = get1Dindex(coord_x, coord_y, coord_z, nx, ny_local);
            // Store TM, TL, CR values for this solidification event in LayerTimeTempHistory
            LayerTimeTempHistory_Host(index, NumberOfSolidificationEvents_Host(index), 0) =
                round(TMelting / deltat) + 1;
            LayerTimeTempHistory_Host(index, NumberOfSolidificationEvents_Host(index), 1) =
                round(TLiquidus / deltat) + 1;
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
            std::cout << "Largest time globally for layer " << layernumber << " is " << LargestTime_Global << std::endl;
        FinishTimeStep[layernumber] = round((LargestTime_Global) / deltat);
        if (id == 0)
            std::cout << " Layer " << layernumber << " FINISH TIME STEP IS " << FinishTimeStep[layernumber]
                      << std::endl;
        if (id == 0)
            std::cout << "Layer " << layernumber << " temperatures read" << std::endl;

        // Reorder solidification events in LayerTimeTempHistory(location,event number,component) so that they are in
        // order based on the melting time values (component = 0)
        for (int index = 0; index < DomainSize; index++) {
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
        for (int index = 0; index < DomainSize; index++) {
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
            std::cout << "Layer " << layernumber << " temperature field is from Z = " << z_layer_bottom << " through "
                      << nz_layer + z_layer_bottom - 1 << " of the global domain" << std::endl;
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

    // Reset solidification event counter and the undercooling to zero for all cells, resizing for the number of cells
    // associated with the next layer's domain
    void reset_layer_events_undercooling(const int DomainSize) {
        Kokkos::realloc(UndercoolingCurrent, DomainSize);
        Kokkos::realloc(SolidificationEventCounter, DomainSize);
        Kokkos::deep_copy(UndercoolingCurrent, 0.0);
        Kokkos::deep_copy(SolidificationEventCounter, 0);
    }

    // Extract the next time that this point undergoes melting
    KOKKOS_INLINE_FUNCTION
    int getMeltTimeStep(const int index) const {
        int SolidificationEventCounter_cell = SolidificationEventCounter(index);
        int MeltTimeStep = static_cast<int>(LayerTimeTempHistory(index, SolidificationEventCounter_cell, 0));
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
    ExtractedViewDataType extract_tm_tl_cr_data(const int extracted_val, const int DomainSize,
                                                const int DefaultVal = -1) {
        ExtractedViewDataType ExtractedData(Kokkos::ViewAllocateWithoutInitializing("ExtractedData"), DomainSize);
        using extracted_value_type = typename ExtractedViewDataType::value_type;
        auto LayerTimeTempHistory_local = LayerTimeTempHistory;
        auto NumberOfSolidificationEvents_local = NumberOfSolidificationEvents;
        Kokkos::parallel_for(
            "Extract_tm_tl_cr_data", DomainSize, KOKKOS_LAMBDA(const int &index) {
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
