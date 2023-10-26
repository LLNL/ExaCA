// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_GRID_HPP
#define EXACA_GRID_HPP

#include "CAinputs.hpp"
#include "CAparsefiles.hpp"
#include "CAtypes.hpp"
#include "mpi.h"

#include <Kokkos_Core.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Data for local (individual MPI ranks) and global (across all ranks) grids
struct Grid {

    // Global domain bounds, all layers
    int nx, ny, nz, DomainSize_AllLayers;
    double deltax, XMin, XMax, YMin, YMax, ZMin, ZMax;

    // Variables characterizing local processor grids relative to global domain
    // 1D decomposition in Y: Each MPI rank has a subset consisting of of ny_local cells, out of ny cells in Y
    // Each MPI rank's subdomain is offset by y_offset cells from the lower bound of the domain in Y
    int ny_local, y_offset;
    // Variables characterizing process IDs of neighboring MPI ranks on the grid
    // Positive Y/NegativeY directions are North/South
    int NeighborRank_North, NeighborRank_South;
    // Variables denoting whether or not each MPI rank's grid is at a global domain boundary
    bool AtNorthBoundary, AtSouthBoundary;

    // Multilayer problem information
    ViewD_H ZMinLayer, ZMaxLayer;
    int NumberOfLayers, LayerHeight;
    double HT_deltax;

    // Current layer information for multilayer problems
    int z_layer_bottom, z_layer_top, nz_layer, DomainSize, BottomOfCurrentLayer, TopOfCurrentLayer;
    std::pair<int, int> LayerRange;

    // Neighbor lists
    NList NeighborX, NeighborY, NeighborZ;

    // TODO: No longer used, should be removed
    int HTtoCAratio;

    // Domain inputs from file
    DomainInputs _inputs;
    // Temperature inputs from file
    TemperatureInputs _t_inputs;

    Grid(std::string SimulationType, const int id, const int np, const int NumberOfLayers_temp, DomainInputs inputs,
         TemperatureInputs t_inputs)
        : ZMinLayer(ViewD_H(Kokkos::ViewAllocateWithoutInitializing("ZMinLayer"), NumberOfLayers_temp))
        , ZMaxLayer(ViewD_H(Kokkos::ViewAllocateWithoutInitializing("ZMaxLayer"), NumberOfLayers_temp))
        , _inputs(inputs)
        , _t_inputs(t_inputs) {

        // Copy from inputs structs
        deltax = _inputs.deltax;
        LayerHeight = _inputs.LayerHeight;
        HT_deltax = _t_inputs.HT_deltax;
        LayerHeight = _inputs.LayerHeight;
        NumberOfLayers = NumberOfLayers_temp;

        // Obtain global domain bounds
        // For problem type R (reading data from file), need to parse temperature data files to obtain domain bounds for
        // each layer and for the multilayer domain
        if (SimulationType == "R") {
            // For simulations using input temperature data with remelting: even if only LayerwiseTempRead is true, all
            // files need to be read to determine the domain bounds
            find_xyz_bounds(id);
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
            // Adjust deltax to exact value based on temperature data spacing and ratio between heat transport/CA cell
            // sizes
            deltax = HT_deltax / HTtoCAratio_floor;
            HTtoCAratio = round(HT_deltax / deltax); // OpenFOAM/CA cell size ratio
        }
        else {
            // Copy inputs from inputs struct into grid struct
            nx = _inputs.nx;
            ny = _inputs.ny;
            nz = _inputs.nz;
            // Domain origin is placed at 0
            XMin = 0.0;
            YMin = 0.0;
            ZMin = 0.0;
            XMax = nx * deltax;
            YMax = ny * deltax;
            ZMax = nz * deltax;
            if (SimulationType == "S") {
                // If this is a spot melt problem, also set the ZMin/ZMax for each layer
                for (int n = 0; n < NumberOfLayers; n++) {
                    ZMinLayer(n) = deltax * (LayerHeight * n);
                    ZMaxLayer(n) = deltax * (_inputs.SpotRadius + LayerHeight * n);
                }
            }
        }
        if (id == 0) {
            std::cout << "Domain size: " << nx << " by " << ny << " by " << nz << std::endl;
            std::cout << "X Limits of domain: " << XMin << " and " << XMax << std::endl;
            std::cout << "Y Limits of domain: " << YMin << " and " << YMax << std::endl;
            std::cout << "Z Limits of domain: " << ZMin << " and " << ZMax << std::endl;
            std::cout << "================================================================" << std::endl;
        }

        // Set cell X, Y, and Z direction neighbors
        neighbor_list_init();

        // Decompose the domain into subdomains on each MPI rank: Calculate ny_local and y_offset for each rank, where
        // each subdomain contains "ny_local" in Y, offset from the full domain origin by "y_offset" cells in Y
        // (previously "DomainDecomposition" in CAinitialize.cpp) First, compare total MPI ranks to total Y cells.
        if (np > ny)
            throw std::runtime_error(
                "Error: Cannot run with more MPI ranks than cells in Y (decomposition direction).");

        // The following was previously performed in "DomainDecomposition" in CAinitialize.cpp:
        // Domain size across all ranks and all layers
        DomainSize_AllLayers = getDomainSize_AllLayers();

        // Get neighboring ranks to each direction and set boundary indicator variables
        NeighborRank_North = getNeighborRank_North(id, np);
        NeighborRank_South = getNeighborRank_South(id, np);
        AtNorthBoundary = getAtNorthBoundary();
        AtSouthBoundary = getAtSouthBoundary();

        // Determine, for each MPI process id, the local grid size in y (and the offset in y relative to the overall
        // simulation domain)
        y_offset = get_yoffset(id, np);
        ny_local = get_nylocal(id, np);

        // Add halo regions with a width of 1 in +/- Y if this MPI rank is not as a domain boundary in said direction
        add_halo();
            std::cout << "Rank " << id << " spans Y = " << y_offset << " through " << y_offset + ny_local - 1 << std::endl;

        // Bounds of layer 0: Z coordinates span z_layer_bottom-z_layer_top, inclusive (functions previously in
        // CAinitialize.cpp and CAcelldata.hpp)
        z_layer_bottom = calc_z_layer_bottom(SimulationType, 0);
        z_layer_top = calc_z_layer_top(SimulationType, _inputs.SpotRadius, 0);
        nz_layer = calc_nz_layer(id, 0);
        DomainSize = calcDomainSize(); // Number of cells in the current layer on this MPI rank
        BottomOfCurrentLayer = getBottomOfCurrentLayer();
        TopOfCurrentLayer = getTopOfCurrentLayer();
        LayerRange = std::make_pair(BottomOfCurrentLayer, TopOfCurrentLayer);
        MPI_Barrier(MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Mesh initialized: initial domain size is " << nz_layer << " out of " << nz
                      << " total cells in the Z direction" << std::endl;
    };

    // For simulation type R, obtain the physical XYZ bounds of the domain by reading temperature data files and parsing
    // the coordinates. Previously in CAinitialize.cpp
    void find_xyz_bounds(int id) {

        // Two passes through reading temperature data files- the first pass only reads the headers to
        // determine units and X/Y/Z bounds of the simulaton domain. Using the X/Y/Z bounds of the simulation domain,
        // nx, ny, and nz can be calculated and the domain decomposed among MPI processes. The maximum number of
        // remelting events in the simulation can also be calculated. The second pass reads the actual X/Y/Z/liquidus
        // time/cooling rate data and each rank stores the data relevant to itself in "RawData" - this is done in the
        // subroutine "ReadTemperatureData"
        XMin = std::numeric_limits<double>::max();
        YMin = std::numeric_limits<double>::max();
        ZMin = std::numeric_limits<double>::max();
        XMax = std::numeric_limits<double>::lowest();
        YMax = std::numeric_limits<double>::lowest();
        ZMax = std::numeric_limits<double>::lowest();

        // Read the first temperature file
        std::ifstream FirstTemperatureFile;
        FirstTemperatureFile.open(_t_inputs.temp_paths[0]);
        std::string FirstLineFirstFile;
        // Header line
        getline(FirstTemperatureFile, FirstLineFirstFile);

        // Read all data files to determine the domain bounds, max number of remelting events
        // for simulations with remelting
        int LayersToRead = std::min(NumberOfLayers, _t_inputs.TempFilesInSeries); // was given in input file
        for (int LayerReadCount = 1; LayerReadCount <= LayersToRead; LayerReadCount++) {

            std::string tempfile_thislayer = _t_inputs.temp_paths[LayerReadCount - 1];
            // Get min and max x coordinates in this file, which can be a binary or ASCII input file
            // binary file type uses extension .catemp, all other file types assumed to be comma-separated ASCII input
            bool BinaryInputData = checkTemperatureFileFormat(tempfile_thislayer);
            // { Xmin, Xmax, Ymin, Ymax, Zmin, Zmax }
            std::array<double, 6> XYZMinMax_ThisLayer =
                parseTemperatureCoordinateMinMax(tempfile_thislayer, BinaryInputData);

            // Based on the input file's layer offset, adjust ZMin/ZMax from the temperature data coordinate
            // system to the multilayer CA coordinate system Check to see in the XYZ bounds for this layer are
            // also limiting for the entire multilayer CA coordinate system
            XYZMinMax_ThisLayer[4] += deltax * LayerHeight * (LayerReadCount - 1);
            XYZMinMax_ThisLayer[5] += deltax * LayerHeight * (LayerReadCount - 1);
            if (XYZMinMax_ThisLayer[0] < XMin)
                XMin = XYZMinMax_ThisLayer[0];
            if (XYZMinMax_ThisLayer[1] > XMax)
                XMax = XYZMinMax_ThisLayer[1];
            if (XYZMinMax_ThisLayer[2] < YMin)
                YMin = XYZMinMax_ThisLayer[2];
            if (XYZMinMax_ThisLayer[3] > YMax)
                YMax = XYZMinMax_ThisLayer[3];
            if (XYZMinMax_ThisLayer[4] < ZMin)
                ZMin = XYZMinMax_ThisLayer[4];
            if (XYZMinMax_ThisLayer[5] > ZMax)
                ZMax = XYZMinMax_ThisLayer[5];
            ZMinLayer[LayerReadCount - 1] = XYZMinMax_ThisLayer[4];
            ZMaxLayer[LayerReadCount - 1] = XYZMinMax_ThisLayer[5];
            if (id == 0)
                std::cout << "Layer = " << LayerReadCount << " Z Bounds are " << XYZMinMax_ThisLayer[4] << " "
                          << XYZMinMax_ThisLayer[5] << std::endl;
        }
        // Extend domain in Z (build) direction if the number of layers are simulated is greater than the number
        // of temperature files read
        if (NumberOfLayers > _t_inputs.TempFilesInSeries) {
            for (int LayerReadCount = _t_inputs.TempFilesInSeries; LayerReadCount < NumberOfLayers; LayerReadCount++) {
                if (_t_inputs.TempFilesInSeries == 1) {
                    // Only one temperature file was read, so the upper Z bound should account for an additional
                    // "NumberOfLayers-1" worth of data Since all layers have the same temperature data, each
                    // layer's "ZMinLayer" is just translated from that of the first layer
                    ZMinLayer[LayerReadCount] = ZMinLayer[LayerReadCount - 1] + deltax * LayerHeight;
                    ZMaxLayer[LayerReadCount] = ZMaxLayer[LayerReadCount - 1] + deltax * LayerHeight;
                    ZMax += deltax * LayerHeight;
                }
                else {
                    // "TempFilesInSeries" temperature files was read, so the upper Z bound should account for
                    // an additional "NumberOfLayers-TempFilesInSeries" worth of data
                    int RepeatedFile = (LayerReadCount) % _t_inputs.TempFilesInSeries;
                    int RepeatUnit = LayerReadCount / _t_inputs.TempFilesInSeries;
                    ZMinLayer[LayerReadCount] =
                        ZMinLayer[RepeatedFile] + RepeatUnit * _t_inputs.TempFilesInSeries * deltax * LayerHeight;
                    ZMaxLayer[LayerReadCount] =
                        ZMaxLayer[RepeatedFile] + RepeatUnit * _t_inputs.TempFilesInSeries * deltax * LayerHeight;
                    ZMax += deltax * LayerHeight;
                }
            }
        }

        // Now at the conclusion of "Loop 0", the decomposition can be performed as the domain bounds are known
        // (all header lines from all files have been read)
        // CA cells in each direction span from the lower to the higher bound of the temperature data - without wall
        // cells or padding around the simulation edges
        nx = round((XMax - XMin) / deltax) + 1;
        ny = round((YMax - YMin) / deltax) + 1;
        nz = round((ZMax - ZMin) / deltax) + 1;
    }

    // Intialize neighbor list structures (NeighborX, NeighborY, NeighborZ)
    void neighbor_list_init() {

        // Assignment of neighbors around a cell "X" is as follows (in order of closest to furthest from cell "X")
        // Neighbors 0 through 8 are in the -Y direction
        // Neighbors 9 through 16 are in the XY plane with cell X
        // Neighbors 17 through 25 are in the +Y direction
        NeighborX = {0, 1, -1, 0, 0, -1, 1, -1, 1, 0, 0, 1, -1, 1, -1, 1, -1, 0, 1, -1, 0, 0, 1, -1, 1, -1};
        NeighborY = {-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        NeighborZ = {0, 0, 0, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 0, 0, 0, 0, 0, 1, -1, 1, 1, -1, -1};
    }

    int getDomainSize_AllLayers() {
        int DomainSize_AllLayers_local = nx * ny * nz;
        return DomainSize_AllLayers_local;
    }

    // Determine the MPI rank of the processor in the +Y direction
    int getNeighborRank_North(const int id, const int np) {
        int NeighborRank_North_local;
        if (np > 1) {
            NeighborRank_North_local = id + 1;
            if (id == np - 1)
                NeighborRank_North_local = MPI_PROC_NULL;
        }
        else {
            // No MPI communication
            NeighborRank_North_local = MPI_PROC_NULL;
        }
        return NeighborRank_North_local;
    }

    // Determine the MPI rank of the processor in the -Y direction
    int getNeighborRank_South(const int id, const int np) {
        int NeighborRank_South_local;
        if (np > 1) {
            NeighborRank_South_local = id - 1;
            if (id == 0)
                NeighborRank_South_local = MPI_PROC_NULL;
        }
        else {
            // No MPI communication
            NeighborRank_South_local = MPI_PROC_NULL;
        }
        return NeighborRank_South_local;
    }

    // Determine if this MPI rank is at the domain boundary in the +Y direction
    bool getAtNorthBoundary() {
        bool AtNorthBoundary_local;
        if (NeighborRank_North == MPI_PROC_NULL)
            AtNorthBoundary_local = true;
        else
            AtNorthBoundary_local = false;
        return AtNorthBoundary_local;
    }

    // Determine if this MPI rank is at the domain boundary in the -Y direction
    bool getAtSouthBoundary() {
        bool AtSouthBoundary_local;
        if (NeighborRank_South == MPI_PROC_NULL)
            AtSouthBoundary_local = true;
        else
            AtSouthBoundary_local = false;
        return AtSouthBoundary_local;
    }

    // Subdivide ny into ny_local across np ranks as evenly as possible, return ny_local on rank id
    int get_nylocal(const int id, const int np) {
        int ny_local_local = 0;
        int ny_local_est = ny / np;
        int YRemainder = ny % np;
        if (YRemainder == 0) {
            ny_local_local = ny_local_est;
        }
        else {
            if (YRemainder > id) {
                ny_local_local = ny_local_est + 1;
            }
            else {
                ny_local_local = ny_local_est;
            }
        }
        return ny_local_local;
    }

    // Get the offset in the Y direction of the local grid on rank id from the global grid of np ranks, which has ny
    // total cells in Y
    int get_yoffset(const int id, const int np) {
        int y_offset_local = 0;
        int y_offset_est = ny / np;
        int YRemainder = ny % np;
        if (YRemainder == 0) {
            y_offset_local = id * y_offset_est;
        }
        else {
            if (YRemainder > id) {
                y_offset_local = id * (y_offset_est + 1);
            }
            else {
                y_offset_local =
                    (y_offset_est + 1) * (YRemainder - 1) + y_offset_est + 1 + (id - YRemainder) * y_offset_est;
            }
        }
        return y_offset_local;
    }

    // Add ghost nodes to the appropriate subdomains (added where the subdomains overlap, but not at edges of physical
    // domain)
    void add_halo() {

        // Add halo regions in Y direction if this subdomain borders subdomains on other processors
        // If only 1 rank in the y direction, no halo regions - subdomain is coincident with overall simulation domain
        // If multiple ranks in the y direction, either 1 halo region (borders another rank's subdomain in either the +y
        // or -y direction) or 2 halo regions (if it borders other rank's subdomains in both the +y and -y directions)
        if (NeighborRank_North != MPI_PROC_NULL)
            ny_local++;
        if (NeighborRank_South != MPI_PROC_NULL) {
            ny_local++;
            // Also adjust subdomain offset, as these ghost nodes were added on the -y side of the subdomain
            y_offset--;
        }
    }

    // Get the Z coordinate of the lower bound of iteration for layer layernumber
    int calc_z_layer_bottom(std::string SimulationType, int layernumber) {

        int z_layer_bottom_local = -1; // assign dummy initial value
        if ((SimulationType == "C") || (SimulationType == "SingleGrain")) {
            // Not a multilayer problem, top of "layer" is the top of the overall simulation domain
            z_layer_bottom_local = 0;
        }
        else if (SimulationType == "S") {
            // lower bound of domain is an integer multiple of the layer spacing, since the temperature field is the
            // same for every layer
            z_layer_bottom_local = LayerHeight * layernumber;
        }
        else if (SimulationType == "R") {
            // lower bound of domain is based on the data read from the file(s)
            z_layer_bottom_local = round((ZMinLayer[layernumber] - ZMin) / deltax);
        }
        if (z_layer_bottom == -1)
            throw std::runtime_error("Error: ZBound_Low went uninitialized, problem type must be C, S, or R");
        return z_layer_bottom_local;
    }

    // Get the Z coordinate of the upper bound of iteration for layer layernumber
    int calc_z_layer_top(std::string SimulationType, int SpotRadius, int layernumber) {

        int z_layer_top_local = -1; // assign dummy initial value
        if ((SimulationType == "C") || (SimulationType == "SingleGrain")) {
            // Not a multilayer problem, top of "layer" is the top of the overall simulation domain
            z_layer_top_local = nz - 1;
        }
        else if (SimulationType == "S") {
            // Top of layer is equal to the spot radius for a problem of hemispherical spot solidification, plus an
            // offset depending on the layer number
            z_layer_top_local = SpotRadius + LayerHeight * layernumber;
        }
        else if (SimulationType == "R") {
            // Top of layer comes from the layer's file data (implicitly assumes bottom of layer 0 is the bottom of the
            // overall domain - this should be fixed in the future for edge cases where this isn't true)
            z_layer_top_local = round((ZMaxLayer[layernumber] - ZMin) / deltax);
        }
        if (z_layer_top == -1)
            throw std::runtime_error(
                "Error: ZBound_High went uninitialized, problem type must be SingleGrain, C, S, or R");
        return z_layer_top_local;
    }

    // Calculate the size of the active domain in Z for layer layernumber
    int calc_nz_layer(int id, int layernumber) {
        int nz_layer_local = z_layer_top - z_layer_bottom + 1;
        if (id == 0)
            std::cout << "Layer " << layernumber << "'s active domain is from Z = " << z_layer_bottom << " through "
                      << z_layer_top << " (" << nz_layer_local << ") cells" << std::endl;
        return nz_layer_local;
    }

    // Calculate the size of the domain, as a number of cells, for this layer of a multilayer problem
    int calcDomainSize() {
        int DomainSize_local = nx * ny_local * nz_layer;
        return DomainSize_local;
    }

    // Get 1D coordinate corresponding to the first cell in this layer of a multilayer problem
    int getBottomOfCurrentLayer() {
        int BottomOfCurrentLayer_local = z_layer_bottom * nx * ny_local;
        return BottomOfCurrentLayer_local;
    }

    // Get 1D coordinate corresponding to the cell after the last cell in this layer of a multilayer problem
    int getTopOfCurrentLayer() {
        int TopOfCurrentLayer_local = z_layer_bottom * nx * ny_local + DomainSize;
        return TopOfCurrentLayer_local;
    }

    // Determine new active cell domain size and offset from bottom of global domain
    void init_next_layer(int id, std::string SimulationType, int nextlayernumber, int SpotRadius) {
        z_layer_bottom = calc_z_layer_bottom(SimulationType, nextlayernumber);
        z_layer_top = calc_z_layer_top(SimulationType, SpotRadius, nextlayernumber);
        nz_layer = calc_nz_layer(id, nextlayernumber);
        DomainSize = calcDomainSize();
        BottomOfCurrentLayer = getBottomOfCurrentLayer();
        TopOfCurrentLayer = getTopOfCurrentLayer();
        LayerRange = std::make_pair(BottomOfCurrentLayer, TopOfCurrentLayer);
    }

    // Get the 1D cell coordinate from the x, y, and z cell positions
    KOKKOS_INLINE_FUNCTION
    int get1Dindex(const int coord_x, const int coord_y_local, const int coord_z) const {
        int index = coord_z * nx * ny_local + coord_x * ny_local + coord_y_local;
        return index;
    }
    // TODO: There is probably some creative way to combine these functions and return one object containing the x, y,
    // and z positions Get the z cell position of the cell from the 1D cell coordinate
    KOKKOS_INLINE_FUNCTION
    int getCoordZ(const int index) const {
        int coord_z = index / (nx * ny_local);
        return coord_z;
    }
    // Get the y cell position of the cell from the 1D cell coordinate
    KOKKOS_INLINE_FUNCTION
    int getCoordY(const int index) const {
        int Rem = index % (nx * ny_local);
        int coord_y = Rem % ny_local;
        return coord_y;
    }
    // Get the x cell position of the cell from the 1D cell coordinate
    KOKKOS_INLINE_FUNCTION
    int getCoordX(const int index) const {
        int Rem = index % (nx * ny_local);
        int coord_x = Rem / ny_local;
        return coord_x;
    }
    // Get the z cell position of the cell from the 1D cell coordinate with respect to the overall simulation domain (all MPI ranks)
    KOKKOS_INLINE_FUNCTION
    int getCoordZ_global(const int index) const {
        int coord_z = index / (nx * ny);
        return coord_z;
    }
    // Get the y cell position of the cell from the 1D cell coordinate with respect to the overall simulation domain
    // (all MPI ranks)
    KOKKOS_INLINE_FUNCTION
    int getCoordY_global(const int index) const {
        int Rem = index % (nx * ny);
        int coord_y = Rem % ny;
        return coord_y;
    }
    // Get the x cell position of the cell from the 1D cell coordinate with respect to the overall simulation domain
    // (all MPI ranks)
    KOKKOS_INLINE_FUNCTION
    int getCoordX_global(const int index) const {
        int Rem = index % (nx * ny);
        int coord_x = Rem / ny;
        return coord_x;
    }
};

#endif
