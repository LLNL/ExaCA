// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAinitialize.hpp"

#include "mpi.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <regex>

void checkPowderOverflow(int nx, int ny, int LayerHeight, int NumberOfLayers, Inputs &inputs) {

    // Check to make sure powder grain density is compatible with the number of powder sites
    // If this problem type includes a powder layer of some grain density, ensure that integer overflow won't occur when
    // assigning powder layer GrainIDs
    if (!(inputs.substrateInputs.BaseplateThroughPowder)) {
        long int NumCellsPowderLayers =
            (long int)(nx) * (long int)(ny) * (long int)(LayerHeight) * (long int)(NumberOfLayers - 1);
        long int NumAssignedCellsPowderLayers =
            std::lround(round(static_cast<double>(NumCellsPowderLayers) * inputs.substrateInputs.PowderActiveFraction));
        if (NumAssignedCellsPowderLayers > INT_MAX)
            throw std::runtime_error("Error: A smaller value for powder density is required to avoid potential integer "
                                     "overflow when assigning powder layer GrainID");
    }
}

//*****************************************************************************/
// Intialize neighbor list structures (NeighborX, NeighborY, NeighborZ)
void NeighborListInit(NList &NeighborX, NList &NeighborY, NList &NeighborZ) {

    // Assignment of neighbors around a cell "X" is as follows (in order of closest to furthest from cell "X")
    // Neighbors 0 through 8 are in the -Y direction
    // Neighbors 9 through 16 are in the XY plane with cell X
    // Neighbors 17 through 25 are in the +Y direction
    NeighborX = {0, 1, -1, 0, 0, -1, 1, -1, 1, 0, 0, 1, -1, 1, -1, 1, -1, 0, 1, -1, 0, 0, 1, -1, 1, -1};
    NeighborY = {-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    NeighborZ = {0, 0, 0, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 0, 0, 0, 0, 0, 1, -1, 1, 1, -1, -1};
}

// Check if the temperature data is in ASCII or binary format
bool checkTemperatureFileFormat(std::string tempfile_thislayer) {
    bool BinaryInputData;
    std::size_t found = tempfile_thislayer.find(".catemp");
    if (found == std::string::npos)
        BinaryInputData = false;
    else
        BinaryInputData = true;
    return BinaryInputData;
}

// Obtain the physical XYZ bounds of the domain, using either domain size from the input file, or reading temperature
// data files and parsing the coordinates
void FindXYZBounds(int id, double &deltax, int &nx, int &ny, int &nz, double &XMin, double &XMax, double &YMin,
                   double &YMax, double &ZMin, double &ZMax, double *ZMinLayer, double *ZMaxLayer, int NumberOfLayers,
                   int LayerHeight, Inputs<device_memory_space> &inputs) {

    if (inputs.SimulationType == "R") {
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

        // Read the first temperature file, first line to determine if the "new" OpenFOAM output format (with a 1 line
        // header) is used, or whether the "old" OpenFOAM header (which contains information like the X/Y/Z bounds of
        // the simulation domain) is
        std::ifstream FirstTemperatureFile;
        FirstTemperatureFile.open(inputs.temperatureInputs.temp_paths[0]);
        std::string FirstLineFirstFile;
        getline(FirstTemperatureFile, FirstLineFirstFile);
        std::size_t found = FirstLineFirstFile.find("Number of temperature data points");
        if (found != std::string::npos) {
            // Old temperature data format detected - no longer supported by ExaCA
            std::string error = "Error: Old header and temperature file format no longer supported";
            throw std::runtime_error(error);
        }

        // Read all data files to determine the domain bounds, max number of remelting events
        // for simulations with remelting
        int LayersToRead =
            std::min(NumberOfLayers, inputs.temperatureInputs.TempFilesInSeries); // was given in input file
        for (int LayerReadCount = 1; LayerReadCount <= LayersToRead; LayerReadCount++) {

            std::string tempfile_thislayer = inputs.temperatureInputs.temp_paths[LayerReadCount - 1];
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
        if (NumberOfLayers > inputs.temperatureInputs.TempFilesInSeries) {
            for (int LayerReadCount = inputs.temperatureInputs.TempFilesInSeries; LayerReadCount < NumberOfLayers;
                 LayerReadCount++) {
                if (inputs.temperatureInputs.TempFilesInSeries == 1) {
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
                    int RepeatedFile = (LayerReadCount) % inputs.temperatureInputs.TempFilesInSeries;
                    int RepeatUnit = LayerReadCount / inputs.temperatureInputs.TempFilesInSeries;
                    ZMinLayer[LayerReadCount] =
                        ZMinLayer[RepeatedFile] +
                        RepeatUnit * inputs.temperatureInputs.TempFilesInSeries * deltax * LayerHeight;
                    ZMaxLayer[LayerReadCount] =
                        ZMaxLayer[RepeatedFile] +
                        RepeatUnit * inputs.temperatureInputs.TempFilesInSeries * deltax * LayerHeight;
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
    else {
        // Using fixed G/R values to set temperature field - no temperature data to read
        // Setting physical domain bounds for consistency with problems that use input temperature data
        // Let XMin/YMin/ZMin be equal to 0, XMax/YMax/ZMax equal to nx*deltax, ny*deltax, nz*deltax
        XMin = 0.0;
        YMin = 0.0;
        ZMin = 0.0;
        XMax = nx * deltax;
        YMax = ny * deltax;
        ZMax = nz * deltax;
        // If this is a spot melt problem, also set the ZMin/ZMax for each layer
        for (int n = 0; n < NumberOfLayers; n++) {
            ZMinLayer[n] = deltax * (LayerHeight * n);
            ZMaxLayer[n] = deltax * (inputs.domainInputs.SpotRadius + LayerHeight * n);
        }
    }
    if (id == 0) {
        std::cout << "Domain size: " << nx << " by " << ny << " by " << nz << std::endl;
        std::cout << "X Limits of domain: " << XMin << " and " << XMax << std::endl;
        std::cout << "Y Limits of domain: " << YMin << " and " << YMax << std::endl;
        std::cout << "Z Limits of domain: " << ZMin << " and " << ZMax << std::endl;
        std::cout << "================================================================" << std::endl;
    }
}

// Decompose the domain into subdomains on each MPI rank: Calculate MyYSlices and MyYOffset for each rank, where each
// subdomain contains "MyYSlices" in Y, offset from the full domain origin by "MyYOffset" cells in Y
void DomainDecomposition(int id, int np, int &ny_local, int &y_offset, int &NeighborRank_North, int &NeighborRank_South,
                         int &nx, int &ny, int &nz, int &DomainSize_AllLayers, bool &AtNorthBoundary,
                         bool &AtSouthBoundary) {

    // Compare total MPI ranks to total Y cells.
    if (np > ny)
        throw std::runtime_error("Error: Cannot run with more MPI ranks than cells in Y (decomposition direction).");

    // Determine which subdomains are at which locations on the grid relative to the others
    InitialDecomposition(id, np, NeighborRank_North, NeighborRank_South, AtNorthBoundary, AtSouthBoundary);
    // Determine, for each MPI process id, the local grid size in x and y (and the offsets in x and y relative to the
    // overall simulation domain)
    y_offset = get_yoffset(id, ny, np);
    ny_local = get_nylocal(id, ny, np);

    // Add ghost nodes at subdomain overlaps
    AddGhostNodes(NeighborRank_North, NeighborRank_South, ny_local, y_offset);

    DomainSize_AllLayers = nx * ny_local * nz; // Number of cells on this MPI rank across all layers of the problem
}

//*****************************************************************************/
// Get the Z coordinate of the lower bound of iteration
int calc_z_layer_bottom(std::string SimulationType, int LayerHeight, int layernumber, double *ZMinLayer, double ZMin,
                        double deltax) {

    int z_layer_bottom = -1; // assign dummy initial value
    if ((SimulationType == "C") || (SimulationType == "SingleGrain")) {
        // Not a multilayer problem, top of "layer" is the top of the overall simulation domain
        z_layer_bottom = 0;
    }
    else if (SimulationType == "S") {
        // lower bound of domain is an integer multiple of the layer spacing, since the temperature field is the
        // same for every layer
        z_layer_bottom = LayerHeight * layernumber;
    }
    else if (SimulationType == "R") {
        // lower bound of domain is based on the data read from the file(s)
        z_layer_bottom = round((ZMinLayer[layernumber] - ZMin) / deltax);
    }
    if (z_layer_bottom == -1)
        throw std::runtime_error("Error: ZBound_Low went uninitialized, problem type must be C, S, or R");
    return z_layer_bottom;
}
//*****************************************************************************/
// Get the Z coordinate of the upper bound of iteration
int calc_z_layer_top(std::string SimulationType, int SpotRadius, int LayerHeight, int layernumber, double ZMin,
                     double deltax, int nz, double *ZMaxLayer) {

    int z_layer_top = -1; // assign dummy initial value
    if ((SimulationType == "C") || (SimulationType == "SingleGrain")) {
        // Not a multilayer problem, top of "layer" is the top of the overall simulation domain
        z_layer_top = nz - 1;
    }
    else if (SimulationType == "S") {
        // Top of layer is equal to the spot radius for a problem of hemispherical spot solidification, plus an offset
        // depending on the layer number
        z_layer_top = SpotRadius + LayerHeight * layernumber;
    }
    else if (SimulationType == "R") {
        // Top of layer comes from the layer's file data (implicitly assumes bottom of layer 0 is the bottom of the
        // overall domain - this should be fixed in the future for edge cases where this isn't true)
        z_layer_top = round((ZMaxLayer[layernumber] - ZMin) / deltax);
    }
    if (z_layer_top == -1)
        throw std::runtime_error("Error: ZBound_High went uninitialized, problem type must be C, S, or R");
    return z_layer_top;
}
//*****************************************************************************/
// Calculate the size of the active domain in Z
int calc_nz_layer(int z_layer_bottom, int z_layer_top, int id, int layernumber) {
    int nz_layer = z_layer_top - z_layer_bottom + 1;
    if (id == 0)
        std::cout << "Layer " << layernumber << "'s active domain is from Z = " << z_layer_bottom << " through "
                  << z_layer_top << " (" << nz_layer << ") cells" << std::endl;
    return nz_layer;
}
//*****************************************************************************/
// Calculate the size of the domain, as a number of cells
int calcLayerDomainSize(int nx, int ny_local, int nz_layer) {
    int DomainSize = nx * ny_local * nz_layer;
    return DomainSize;
}
//*****************************************************************************/
void ZeroResetViews(int LocalActiveDomainSize, ViewF &DiagonalLength, ViewF &CritDiagonalLength, ViewF &DOCenter,
                    ViewI &SteeringVector) {

    // Realloc steering vector as LocalActiveDomainSize may have changed (old values aren't needed)
    Kokkos::realloc(SteeringVector, LocalActiveDomainSize);

    // Realloc active cell data structure and halo regions on device (old values not needed)
    Kokkos::realloc(DiagonalLength, LocalActiveDomainSize);
    Kokkos::realloc(DOCenter, 3 * LocalActiveDomainSize);
    Kokkos::realloc(CritDiagonalLength, 26 * LocalActiveDomainSize);

    // Reset active cell data structures on device
    Kokkos::deep_copy(DiagonalLength, 0);
    Kokkos::deep_copy(DOCenter, 0);
    Kokkos::deep_copy(CritDiagonalLength, 0);
}
