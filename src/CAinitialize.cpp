// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAinitialize.hpp"

#include "CAconfig.hpp"
#include "CAfunctions.hpp"
#include "CAghostnodes.hpp"
#include "CAparsefiles.hpp"
#include "CAupdate.hpp"

#include "mpi.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <regex>

// Initializes input parameters, mesh, temperature field, and grain structures for CA simulations

// Read ExaCA input file (JSON format)
void InputReadFromFile(int id, std::string InputFile, std::string &SimulationType, double &deltax, double &NMax,
                       double &dTN, double &dTsigma, std::string &GrainOrientationFile, int &TempFilesInSeries,
                       std::vector<std::string> &temp_paths, double &HT_deltax, double &deltat, int &NumberOfLayers,
                       int &LayerHeight, std::string &MaterialFileName, std::string &SubstrateFileName,
                       float &SubstrateGrainSpacing, bool &UseSubstrateFile, double &G, double &R, int &nx, int &ny,
                       int &nz, double &FractSurfaceSitesActive, int &NSpotsX, int &NSpotsY, int &SpotOffset,
                       int &SpotRadius, double &RNGSeed, bool &BaseplateThroughPowder, double &PowderActiveFraction,
                       bool &LayerwiseTempRead, bool &PowderFirstLayer, Print &print, double &initUndercooling,
                       int &singleGrainOrientation) {

    std::ifstream InputData(InputFile);
    nlohmann::json inputdata = nlohmann::json::parse(InputData);

    // General inputs
    SimulationType = inputdata["SimulationType"];
    // "C": constrained (directional) solidification
    // "S": array of overlapping hemispherical spots
    // "R": time-temperature history comes from external files
    // Check if simulation type includes remelting ("M" suffix to input problem type) - all simulations now use
    // remelting, so in the absence of this suffix, print warning that the problem will use remelting DirSoldification
    // problems now include remelting logic
    if (SimulationType == "RM") {
        SimulationType = "R";
    }
    else if (SimulationType == "SM") {
        SimulationType = "S";
    }
    else if ((SimulationType == "S") || (SimulationType == "R")) {
        if (id == 0) {
            if (SimulationType == "S")
                std::cout << "Warning: While the specified problem type did not include remelting, all simulations now "
                             "include remelting"
                          << std::endl;
        }
    }
    // Input files that should be present for all problem types
    std::string MaterialFileName_Read = inputdata["MaterialFileName"];
    std::string GrainOrientationFile_Read = inputdata["GrainOrientationFile"];
    // Path to file of materials constants based on install/source location
    MaterialFileName = checkFileInstalled(MaterialFileName_Read, id);
    checkFileNotEmpty(MaterialFileName);
    // Path to file of grain orientations based on install/source location
    GrainOrientationFile = checkFileInstalled(GrainOrientationFile_Read, id);
    checkFileNotEmpty(GrainOrientationFile);
    // Seed for random number generator (defaults to 0 if not given)
    if (inputdata.contains("RandomSeed"))
        RNGSeed = inputdata["RandomSeed"];
    else
        RNGSeed = 0.0;

    // Domain inputs:
    // Cell size - given in meters, stored in micrometers
    deltax = inputdata["Domain"]["CellSize"];
    deltax = deltax * pow(10, -6);
    // Time step - given in seconds, stored in microseconds
    deltat = inputdata["Domain"]["TimeStep"];
    deltat = deltat * pow(10, -6);
    if ((SimulationType == "C") || (SimulationType == "SingleGrain")) {
        // Domain size, in cells
        nx = inputdata["Domain"]["Nx"];
        ny = inputdata["Domain"]["Ny"];
        nz = inputdata["Domain"]["Nz"];
        NumberOfLayers = 1;
        LayerHeight = nz;
    }
    else {
        // Number of layers, layer height are needed for problem types S and R
        NumberOfLayers = inputdata["Domain"]["NumberOfLayers"];
        LayerHeight = inputdata["Domain"]["LayerOffset"];
        // Type S needs spot information, which is then used to compute the domain bounds
        if (SimulationType == "S") {
            NSpotsX = inputdata["Domain"]["NSpotsX"];
            NSpotsY = inputdata["Domain"]["NSpotsY"];
            // Radius and offset are given in micrometers, convert to cells
            SpotRadius = inputdata["Domain"]["RSpots"];
            SpotRadius = SpotRadius * pow(10, -6) / deltax;
            SpotOffset = inputdata["Domain"]["SpotOffset"];
            SpotOffset = SpotOffset * pow(10, -6) / deltax;
            // Calculate nx, ny, and nz based on spot array pattern and number of layers
            nz = SpotRadius + 1 + (NumberOfLayers - 1) * LayerHeight;
            nx = 2 * SpotRadius + 1 + SpotOffset * (NSpotsX - 1);
            ny = 2 * SpotRadius + 1 + SpotOffset * (NSpotsY - 1);
        }
    }

    // Nucleation inputs:
    // Nucleation density (normalized by 10^12 m^-3), mean nucleation undercooling/st dev undercooling(K)
    // SingleGrain problem type does not have nucleation, just a grain of a single orientation in the domain center
    if (SimulationType != "SingleGrain") {
        NMax = inputdata["Nucleation"]["Density"];
        NMax = NMax * pow(10, 12);
        dTN = inputdata["Nucleation"]["MeanUndercooling"];
        dTsigma = inputdata["Nucleation"]["StDev"];
    }
    else {
        NMax = 0.0;
        dTN = 0.0;
        dTsigma = 0.0;
    }

    // Temperature inputs:
    if (SimulationType == "R") {
        // Temperature data resolution - default to using CA cell size if the assumed temperature data resolution if not
        // given
        if (inputdata["TemperatureData"].contains("HeatTransferCellSize")) {
            HT_deltax = inputdata["TemperatureData"]["HeatTransferCellSize"];
            // Value is given in micrometers, convert to meters
            HT_deltax = HT_deltax * pow(10, -6);
        }
        else
            HT_deltax = deltax;
        if (HT_deltax != deltax)
            throw std::runtime_error("Error: CA cell size and input temperature data resolution must be equivalent");
        // Read all temperature files at once (default), or one at a time?
        if (inputdata["TemperatureData"].contains("LayerwiseTempRead")) {
            LayerwiseTempRead = inputdata["TemperatureData"]["LayerwiseTempRead"];
        }
        else
            LayerwiseTempRead = false;
        // Get the paths/number of/names of the temperature data files used
        TempFilesInSeries = inputdata["TemperatureData"]["TemperatureFiles"].size();
        if (TempFilesInSeries == 0)
            throw std::runtime_error("Error: No temperature files listed in the temperature instructions file");
        else {
            for (int filename = 0; filename < TempFilesInSeries; filename++)
                temp_paths.push_back(inputdata["TemperatureData"]["TemperatureFiles"][filename]);
        }
    }
    else {
        // Temperature data uses fixed thermal gradient (K/m) and cooling rate (K/s)
        G = inputdata["TemperatureData"]["G"];
        R = inputdata["TemperatureData"]["R"];
        if (SimulationType == "SingleGrain")
            initUndercooling = inputdata["TemperatureData"]["InitUndercooling"];
        else
            initUndercooling = 0.0;
    }

    // Substrate inputs:
    if (SimulationType == "C") {
        // Fraction of sites at bottom surface active
        FractSurfaceSitesActive = inputdata["Substrate"]["FractionSurfaceSitesActive"];
    }
    else if (SimulationType == "SingleGrain") {
        // Orientation of the single grain at the domain center
        singleGrainOrientation = inputdata["Substrate"]["GrainOrientation"];
    }
    else {
        // Substrate data - should data come from an initial size or a file?
        if ((inputdata["Substrate"].contains("SubstrateFilename")) && (inputdata["Substrate"].contains("MeanSize")))
            throw std::runtime_error("Error: only one of substrate grain size and substrate structure filename should "
                                     "be provided in the input file");
        else if (inputdata["Substrate"].contains("SubstrateFilename")) {
            SubstrateFileName = inputdata["Substrate"]["SubstrateFilename"];
            UseSubstrateFile = true;
        }
        else if (inputdata["Substrate"].contains("MeanSize")) {
            SubstrateGrainSpacing = inputdata["Substrate"]["MeanSize"];
            UseSubstrateFile = false;
        }
        // Should the baseplate microstructure be extended through the powder layers? Default is false
        if (inputdata["Substrate"].contains("ExtendSubstrateThroughPower"))
            BaseplateThroughPowder = inputdata["Substrate"]["ExtendSubstrateThroughPower"];
        else
            BaseplateThroughPowder = false;
        if (inputdata["Substrate"].contains("PowderDensity")) {
            // powder density is given as a density per unit volume, normalized by 10^12 m^-3 --> convert this into a
            // density of sites active on the CA grid (0 to 1)
            PowderActiveFraction = inputdata["Substrate"]["PowderDensity"];
            PowderActiveFraction = PowderActiveFraction * pow(10, 12) * pow(deltax, 3);
            if ((PowderActiveFraction < 0.0) || (PowderActiveFraction > 1.0))
                throw std::runtime_error("Error: Density of powder surface sites active must be larger than 0 and less "
                                         "than 1/(CA cell volume)");
        }
        else
            PowderActiveFraction = 1.0; // defaults to a unique grain at each site in the powder layers
        // Should a powder layer be initialized at the top of the baseplate for the first layer? (Defaults to false,
        // where the baseplate spans all of layer 0, and only starting with layer 1 is a powder layer present)
        if (inputdata["Substrate"].contains("PowderFirstLayer"))
            PowderFirstLayer = inputdata["Substrate"]["PowderFirstLayer"];
        else
            PowderFirstLayer = false;
        if ((BaseplateThroughPowder) && ((PowderFirstLayer) || (inputdata["Substrate"].contains("PowderDensity"))))
            throw std::runtime_error(
                "Error: if the option to extend the baseplate through the powder layers is toggled, options regarding "
                "the powder layer (PowderFirstLayer/PowderDensity cannot be given");
    }

    // Printing inputs:
    print.getPrintDataFromInputFile(inputdata, id, deltat);

    // Print information to console about the input file data read
    if (id == 0) {
        std::cout << "Material simulated is " << MaterialFileName << std::endl;
        std::cout << "CA cell size is " << deltax * pow(10, 6) << " microns" << std::endl;
        std::cout << "Nucleation density is " << NMax << " per m^3" << std::endl;
        std::cout << "Mean nucleation undercooling is " << dTN << " K, standard deviation of distribution is "
                  << dTsigma << "K" << std::endl;
        if (SimulationType == "C") {
            std::cout << "CA Simulation using a unidirectional, fixed thermal gradient of " << G
                      << " K/m and a cooling rate of " << R << " K/s" << std::endl;
            std::cout << "The time step is " << deltat * pow(10, 6) << " microseconds" << std::endl;
            std::cout << "The fraction of CA cells at the bottom surface that are active is " << FractSurfaceSitesActive
                      << std::endl;
        }
        else if (SimulationType == "S") {
            std::cout << "CA Simulation using a radial, fixed thermal gradient of " << G
                      << " K/m as a series of hemispherical spots, and a cooling rate of " << R << " K/s" << std::endl;
            std::cout << "A total of " << NumberOfLayers << " spots per layer, with layers offset by " << LayerHeight
                      << " CA cells will be simulated" << std::endl;
            std::cout << "The time step is " << deltat * pow(10, 6) << " microseconds" << std::endl;
        }
        else if (SimulationType == "R") {
            std::cout << "CA Simulation using temperature data from file(s)" << std::endl;
            std::cout << "The time step is " << deltat << " seconds" << std::endl;
            std::cout << "The first temperature data file to be read is " << temp_paths[0] << ", and there are "
                      << TempFilesInSeries << " in the series" << std::endl;
            std::cout << "A total of " << NumberOfLayers << " layers of solidification offset by " << LayerHeight
                      << " CA cells will be simulated" << std::endl;
        }
    }
}

void checkPowderOverflow(int nx, int ny, int LayerHeight, int NumberOfLayers, bool BaseplateThroughPowder,
                         double PowderActiveFraction) {

    // Check to make sure powder grain density is compatible with the number of powder sites
    // If this problem type includes a powder layer of some grain density, ensure that integer overflow won't occur when
    // assigning powder layer GrainIDs
    if (!(BaseplateThroughPowder)) {
        long int NumCellsPowderLayers =
            (long int)(nx) * (long int)(ny) * (long int)(LayerHeight) * (long int)(NumberOfLayers - 1);
        long int NumAssignedCellsPowderLayers =
            std::lround(round(static_cast<double>(NumCellsPowderLayers) * PowderActiveFraction));
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
void FindXYZBounds(std::string SimulationType, int id, double &deltax, int &nx, int &ny, int &nz,
                   std::vector<std::string> &temp_paths, double &XMin, double &XMax, double &YMin, double &YMax,
                   double &ZMin, double &ZMax, int &LayerHeight, int NumberOfLayers, int TempFilesInSeries,
                   double *ZMinLayer, double *ZMaxLayer, int SpotRadius) {

    if (SimulationType == "R") {
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
        FirstTemperatureFile.open(temp_paths[0]);
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
        int LayersToRead = std::min(NumberOfLayers, TempFilesInSeries); // was given in input file
        for (int LayerReadCount = 1; LayerReadCount <= LayersToRead; LayerReadCount++) {

            std::string tempfile_thislayer = temp_paths[LayerReadCount - 1];
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
        if (NumberOfLayers > TempFilesInSeries) {
            for (int LayerReadCount = TempFilesInSeries; LayerReadCount < NumberOfLayers; LayerReadCount++) {
                if (TempFilesInSeries == 1) {
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
                    int RepeatedFile = (LayerReadCount) % TempFilesInSeries;
                    int RepeatUnit = LayerReadCount / TempFilesInSeries;
                    ZMinLayer[LayerReadCount] =
                        ZMinLayer[RepeatedFile] + RepeatUnit * TempFilesInSeries * deltax * LayerHeight;
                    ZMaxLayer[LayerReadCount] =
                        ZMaxLayer[RepeatedFile] + RepeatUnit * TempFilesInSeries * deltax * LayerHeight;
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
            ZMaxLayer[n] = deltax * (SpotRadius + LayerHeight * n);
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
