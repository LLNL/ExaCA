// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_INPUTS_HPP
#define EXACA_INPUTS_HPP

#include "CAinfo.hpp"
#include "CAinterfacialresponse.hpp"
#include "CAtypes.hpp"

#include "mpi.h"

#include <nlohmann/json.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Structs to organize data within inputs struct
struct DomainInputs {
    double deltax, deltat;
    // number of CA cells in each direction only initialized here for problem types C, S, and SingleGrain
    int nx, ny, nz;
    // multilayer problems only
    int NumberOfLayers, LayerHeight;
    // problem type S only
    int NSpotsX, NSpotsY, SpotRadius, SpotOffset;
};

struct NucleationInputs {
    // unused for single grain problem type, zero by default
    double NMax = 0.0;
    double dTN = 0.0;
    double dTsigma = 0.0;
};

struct TemperatureInputs {
    // Used for problem type R (by default, all temperature files read during init)
    bool LayerwiseTempRead = false;
    int TempFilesInSeries;
    double HT_deltax;
    std::vector<std::string> temp_paths;
    // Use for problem types other than R (no temperature files to read) - default to no initial undercooling at
    // solidification front
    double G, R;
    double initUndercooling = 0.0;
};

struct SubstrateInputs {
    // problem type C only
    double FractSurfaceSitesActive;
    // problem type SingleGrain only
    double singleGrainOrientation;
    // problem types S and R only
    bool UseSubstrateFile = false;
    bool BaseplateThroughPowder = false;
    std::string SubstrateFileName;
    double SubstrateGrainSpacing;
    // defaults to all sites in powder layer initialized with a new grain
    double PowderActiveFraction = 1.0;
    // Top of baseplate assumed at Z = 0 if not otherwise given
    double BaseplateTopZ;
};

struct PrintInputs {
    // Base name of CA output
    std::string BaseFileName;
    // Path to CA output
    std::string PathToOutput = "";
    // Fields to be printed at start of run: GrainID, LayerID, MeltTimeStep, CritTimeStep, UndercoolingChange (all
    // for first layer)
    bool PrintInitGrainID = false;
    bool PrintInitLayerID = false;
    bool PrintInitMeltTimeStep = false;
    bool PrintInitCritTimeStep = false;
    bool PrintInitUndercoolingChange = false;
    // Fields to be printed at end of run: GrainID, LayerID, GrainMisorientation, UndercoolingCurrent (for whole
    // domain) and MeltTimeStep, CritTimeStep, UndercoolingChange, CellType (for the final layer)
    bool PrintFinalGrainID = false;
    bool PrintFinalLayerID = false;
    bool PrintFinalMisorientation = false;
    bool PrintFinalUndercoolingCurrent = false;
    bool PrintFinalMeltTimeStep = false;
    bool PrintFinalCritTimeStep = false;
    bool PrintFinalUndercoolingChange = false;
    bool PrintFinalCellType = false;

    // Print intermediate output of grain misorientation in time
    bool PrintTimeSeries = false;
    // If printing the time series, should microstructure be printed if unchanged from the previous frame
    bool PrintIdleTimeSeriesFrames = false;
    // If printing the time series, the increment, in time steps, between printing of intermediate output
    int TimeSeriesInc = 1;
    // If printing the time series, the counter for the number of intermediate files that have been printed for a
    // given layer of a multilayer problem
    int IntermediateFileCounter = 0;

    // Should binary be used for printing vtk data?
    bool PrintBinary = false;

    // Should the default RVE data for ExaConstit be printed? If so, with what size?
    bool PrintDefaultRVE = false;
    int RVESize = 0;
};

template <typename MemorySpace>
struct Inputs {

    using memory_space = MemorySpace;

    std::string SimulationType, MaterialFileName, GrainOrientationFile;
    double RNGSeed = 0.0;
    DomainInputs domainInputs;
    NucleationInputs nucleationInputs;
    TemperatureInputs temperatureInputs;
    SubstrateInputs substrateInputs;
    PrintInputs printInputs;

    // Creates input struct with uninitialized/default values, used in unit tests
    Inputs(){};

    Inputs(int id, std::string InputFile) {

        // Open and read JSON input file
        std::ifstream InputData(InputFile);
        nlohmann::json inputdata = nlohmann::json::parse(InputData);

        // General inputs
        SimulationType = inputdata["SimulationType"];
        // "C": constrained (directional) solidification
        // "S": array of overlapping hemispherical spots
        // "R": time-temperature history comes from external files
        // Check if simulation type includes remelting ("M" suffix to input problem type) - all simulations now use
        // remelting, so in the absence of this suffix, print warning that the problem will use remelting
        // DirSoldification problems now include remelting logic
        if (SimulationType == "RM") {
            SimulationType = "R";
        }
        else if (SimulationType == "SM") {
            SimulationType = "S";
        }
        else if ((SimulationType == "S") || (SimulationType == "R")) {
            if (id == 0) {
                if (SimulationType == "S")
                    std::cout
                        << "Warning: While the specified problem type did not include remelting, all simulations now "
                           "include remelting"
                        << std::endl;
            }
        }
        if ((SimulationType == "S") && (id == 0))
            std::cout << "Warning: The spot melt array simulation type (Problem type S) is now deprecated and will be "
                         "removed in a future release"
                      << std::endl;

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

        // Domain inputs:
        // Cell size - given in meters, stored in micrometers
        domainInputs.deltax = inputdata["Domain"]["CellSize"];
        domainInputs.deltax = domainInputs.deltax * pow(10, -6);
        // Time step - given in seconds, stored in microseconds
        domainInputs.deltat = inputdata["Domain"]["TimeStep"];
        domainInputs.deltat = domainInputs.deltat * pow(10, -6);
        if ((SimulationType == "C") || (SimulationType == "SingleGrain")) {
            // Domain size, in cells
            domainInputs.nx = inputdata["Domain"]["Nx"];
            domainInputs.ny = inputdata["Domain"]["Ny"];
            domainInputs.nz = inputdata["Domain"]["Nz"];
            domainInputs.NumberOfLayers = 1;
            domainInputs.LayerHeight = domainInputs.nz;
        }
        else {
            // Number of layers, layer height are needed for problem types S and R
            domainInputs.NumberOfLayers = inputdata["Domain"]["NumberOfLayers"];
            domainInputs.LayerHeight = inputdata["Domain"]["LayerOffset"];
            // Type S needs spot information, which is then used to compute the domain bounds
            if (SimulationType == "S") {
                domainInputs.NSpotsX = inputdata["Domain"]["NSpotsX"];
                domainInputs.NSpotsY = inputdata["Domain"]["NSpotsY"];
                // Radius and offset are given in micrometers, convert to cells
                domainInputs.SpotRadius = inputdata["Domain"]["RSpots"];
                domainInputs.SpotRadius = domainInputs.SpotRadius * pow(10, -6) / domainInputs.deltax;
                domainInputs.SpotOffset = inputdata["Domain"]["SpotOffset"];
                domainInputs.SpotOffset = domainInputs.SpotOffset * pow(10, -6) / domainInputs.deltax;
                // Calculate nx, ny, and nz based on spot array pattern and number of layers
                domainInputs.nz =
                    domainInputs.SpotRadius + 1 + (domainInputs.NumberOfLayers - 1) * domainInputs.LayerHeight;
                domainInputs.nx =
                    2 * domainInputs.SpotRadius + 1 + domainInputs.SpotOffset * (domainInputs.NSpotsX - 1);
                domainInputs.ny =
                    2 * domainInputs.SpotRadius + 1 + domainInputs.SpotOffset * (domainInputs.NSpotsY - 1);
            }
        }

        // Nucleation inputs:
        // Nucleation density (normalized by 10^12 m^-3), mean nucleation undercooling/st dev undercooling(K)
        // SingleGrain problem type does not have nucleation, just a grain of a single orientation in the domain center
        if (SimulationType != "SingleGrain") {
            nucleationInputs.NMax = inputdata["Nucleation"]["Density"];
            nucleationInputs.NMax = nucleationInputs.NMax * pow(10, 12);
            nucleationInputs.dTN = inputdata["Nucleation"]["MeanUndercooling"];
            nucleationInputs.dTsigma = inputdata["Nucleation"]["StDev"];
        }

        // Temperature inputs:
        if (SimulationType == "R") {
            // Temperature data resolution - default to using CA cell size if the assumed temperature data resolution if
            // not given
            if (inputdata["TemperatureData"].contains("HeatTransferCellSize")) {
                temperatureInputs.HT_deltax = inputdata["TemperatureData"]["HeatTransferCellSize"];
                // Value is given in micrometers, convert to meters
                temperatureInputs.HT_deltax = temperatureInputs.HT_deltax * pow(10, -6);
            }
            else
                temperatureInputs.HT_deltax = domainInputs.deltax;
            if (temperatureInputs.HT_deltax != domainInputs.deltax)
                throw std::runtime_error(
                    "Error: CA cell size and input temperature data resolution must be equivalent");
            // Read all temperature files at once (default), or one at a time?
            if (inputdata["TemperatureData"].contains("LayerwiseTempRead")) {
                temperatureInputs.LayerwiseTempRead = inputdata["TemperatureData"]["LayerwiseTempRead"];
            }
            // Get the paths/number of/names of the temperature data files used
            temperatureInputs.TempFilesInSeries = inputdata["TemperatureData"]["TemperatureFiles"].size();
            if (temperatureInputs.TempFilesInSeries == 0)
                throw std::runtime_error("Error: No temperature files listed in the temperature instructions file");
            else {
                for (int filename = 0; filename < temperatureInputs.TempFilesInSeries; filename++)
                    temperatureInputs.temp_paths.push_back(inputdata["TemperatureData"]["TemperatureFiles"][filename]);
            }
        }
        else {
            // Temperature data uses fixed thermal gradient (K/m) and cooling rate (K/s)
            temperatureInputs.G = inputdata["TemperatureData"]["G"];
            temperatureInputs.R = inputdata["TemperatureData"]["R"];
            if (SimulationType == "SingleGrain")
                temperatureInputs.initUndercooling = inputdata["TemperatureData"]["InitUndercooling"];
        }

        // Substrate inputs:
        if (SimulationType == "C") {
            // Fraction of sites at bottom surface active
            substrateInputs.FractSurfaceSitesActive = inputdata["Substrate"]["FractionSurfaceSitesActive"];
        }
        else if (SimulationType == "SingleGrain") {
            // Orientation of the single grain at the domain center
            substrateInputs.singleGrainOrientation = inputdata["Substrate"]["GrainOrientation"];
        }
        else {
            // Substrate data - should data come from an initial size or a file?
            if ((inputdata["Substrate"].contains("SubstrateFilename")) && (inputdata["Substrate"].contains("MeanSize")))
                throw std::runtime_error(
                    "Error: only one of substrate grain size and substrate structure filename should "
                    "be provided in the input file");
            else if (inputdata["Substrate"].contains("SubstrateFilename")) {
                substrateInputs.SubstrateFileName = inputdata["Substrate"]["SubstrateFilename"];
                substrateInputs.UseSubstrateFile = true;
            }
            else if (inputdata["Substrate"].contains("MeanSize")) {
                substrateInputs.SubstrateGrainSpacing = inputdata["Substrate"]["MeanSize"];
                substrateInputs.UseSubstrateFile = false;
            }
            // Should the baseplate microstructure be extended through the powder layers? Default is false
            if (inputdata["Substrate"].contains("ExtendSubstrateThroughPower"))
                substrateInputs.BaseplateThroughPowder = inputdata["Substrate"]["ExtendSubstrateThroughPower"];
            // defaults to a unique grain at each site in the powder layers if not given
            if (inputdata["Substrate"].contains("PowderDensity")) {
                // powder density is given as a density per unit volume, normalized by 10^12 m^-3 --> convert this into
                // a density of sites active on the CA grid (0 to 1)
                substrateInputs.PowderActiveFraction = inputdata["Substrate"]["PowderDensity"];
                substrateInputs.PowderActiveFraction =
                    substrateInputs.PowderActiveFraction * pow(10, 12) * pow(domainInputs.deltax, 3);
                if ((substrateInputs.PowderActiveFraction < 0.0) || (substrateInputs.PowderActiveFraction > 1.0))
                    throw std::runtime_error(
                        "Error: Density of powder surface sites active must be larger than 0 and less "
                        "than 1/(CA cell volume)");
            }
            if ((inputdata["Substrate"].contains("PowderFirstLayer")) && (id == 0))
                std::cout << "Warning: PowderFirstLayer input is no longer used, the top of the first layer must be "
                             "specified using BaseplateTopZ (which will otherwise default to Z = 0)"
                          << std::endl;
            // The top of the baseplate is designated using BaseplateTopZ (assumed to be Z = 0 if not given in input
            // file)
            if (inputdata["Substrate"].contains("BaseplateTopZ"))
                substrateInputs.BaseplateTopZ = inputdata["Substrate"]["BaseplateTopZ"];
            else
                substrateInputs.BaseplateTopZ = 0.0;
            if ((substrateInputs.BaseplateThroughPowder) && (inputdata["Substrate"].contains("PowderDensity")))
                throw std::runtime_error("Error: if the option to extend the baseplate through the powder layers is "
                                         "toggled, a powder layer density cannot be given");
        }

        // Printing inputs
        getPrintDataFromInputFile(inputdata, id, domainInputs.deltat);

        // Print information to console about the input file data read
        if (id == 0) {
            std::cout << "Material simulated is " << MaterialFileName << std::endl;
            std::cout << "CA cell size is " << domainInputs.deltax * pow(10, 6) << " microns" << std::endl;
            std::cout << "Nucleation density is " << nucleationInputs.NMax << " per m^3" << std::endl;
            std::cout << "Mean nucleation undercooling is " << nucleationInputs.dTN
                      << " K, standard deviation of distribution is " << nucleationInputs.dTsigma << "K" << std::endl;
            if (SimulationType == "C") {
                std::cout << "CA Simulation using a unidirectional, fixed thermal gradient of " << temperatureInputs.G
                          << " K/m and a cooling rate of " << temperatureInputs.R << " K/s" << std::endl;
                std::cout << "The time step is " << domainInputs.deltat * pow(10, 6) << " microseconds" << std::endl;
                std::cout << "The fraction of CA cells at the bottom surface that are active is "
                          << substrateInputs.FractSurfaceSitesActive << std::endl;
            }
            else if (SimulationType == "S") {
                std::cout << "CA Simulation using a radial, fixed thermal gradient of " << temperatureInputs.G
                          << " K/m as a series of hemispherical spots, and a cooling rate of " << temperatureInputs.R
                          << " K/s" << std::endl;
                std::cout << "A total of " << domainInputs.NumberOfLayers << " spots per layer, with layers offset by "
                          << domainInputs.LayerHeight << " CA cells will be simulated" << std::endl;
                std::cout << "The time step is " << domainInputs.deltat * pow(10, 6) << " microseconds" << std::endl;
            }
            else if (SimulationType == "R") {
                std::cout << "CA Simulation using temperature data from file(s)" << std::endl;
                std::cout << "The time step is " << domainInputs.deltat << " seconds" << std::endl;
                std::cout << "The first temperature data file to be read is " << temperatureInputs.temp_paths[0]
                          << ", and there are " << temperatureInputs.TempFilesInSeries << " in the series" << std::endl;
                std::cout << "A total of " << domainInputs.NumberOfLayers << " layers of solidification offset by "
                          << domainInputs.LayerHeight << " CA cells will be simulated" << std::endl;
            }
        }
    }

    // Read the input data file and initialize appropriate variables to non-default values if necessary
    void getPrintDataFromInputFile(nlohmann::json inputdata, int id, double deltat) {
        // Path to output data
        printInputs.PathToOutput = inputdata["Printing"]["PathToOutput"];
        // Name of output data
        printInputs.BaseFileName = inputdata["Printing"]["OutputFile"];
        // Should ASCII or binary be used to print vtk data? Defaults to ASCII if not given
        if (inputdata["Printing"].contains("PrintBinary"))
            printInputs.PrintBinary = inputdata["Printing"]["PrintBinary"];
        // Should default ExaConstit output be printed after the simulation? If so, what size RVE?
        // If a size of 0 is given, this is set to false
        if (inputdata["Printing"].contains("PrintExaConstitSize")) {
            printInputs.RVESize = inputdata["Printing"]["PrintExaConstitSize"];
            if (printInputs.RVESize != 0)
                printInputs.PrintDefaultRVE = true;
        }

        // Which fields should be printed at the start and end of the simulation?
        std::vector<std::string> InitFieldnames_key = {"GrainID", "LayerID", "MeltTimeStep", "CritTimeStep",
                                                       "UndercoolingChange"};
        std::vector<bool> PrintFieldsInit = getPrintFieldValues(inputdata, "PrintFieldsInit", InitFieldnames_key);
        if (PrintFieldsInit[0])
            printInputs.PrintInitGrainID = true;
        if (PrintFieldsInit[1])
            printInputs.PrintInitLayerID = true;
        if (PrintFieldsInit[2])
            printInputs.PrintInitMeltTimeStep = true;
        if (PrintFieldsInit[3])
            printInputs.PrintInitCritTimeStep = true;
        if (PrintFieldsInit[4])
            printInputs.PrintInitUndercoolingChange = true;

        std::vector<std::string> FinalFieldnames_key = {
            "GrainID",      "LayerID",      "GrainMisorientation", "UndercoolingCurrent",
            "MeltTimeStep", "CritTimeStep", "UndercoolingChange",  "CellType"};
        std::vector<bool> PrintFieldsFinal = getPrintFieldValues(inputdata, "PrintFieldsFinal", FinalFieldnames_key);
        if (PrintFieldsFinal[0])
            printInputs.PrintFinalGrainID = true;
        if (PrintFieldsFinal[1])
            printInputs.PrintFinalLayerID = true;
        if (PrintFieldsFinal[2])
            printInputs.PrintFinalMisorientation = true;
        if (PrintFieldsFinal[3])
            printInputs.PrintFinalUndercoolingCurrent = true;
        if (PrintFieldsFinal[4])
            printInputs.PrintFinalMeltTimeStep = true;
        if (PrintFieldsFinal[5])
            printInputs.PrintFinalCritTimeStep = true;
        if (PrintFieldsFinal[6])
            printInputs.PrintFinalUndercoolingChange = true;
        if (PrintFieldsFinal[7])
            printInputs.PrintFinalCellType = true;

        // Should intermediate output be printed?
        if (inputdata["Printing"].contains("PrintIntermediateOutput")) {
            // An increment of 0 will set the intermediate file printing to false
            double TimeSeriesFrameInc_time = inputdata["Printing"]["PrintIntermediateOutput"]["Frequency"];
            if (TimeSeriesFrameInc_time != 0) {
                printInputs.PrintTimeSeries = true;
                // Increment is given in microseconds, convert to seconds
                TimeSeriesFrameInc_time = TimeSeriesFrameInc_time * pow(10, -6);
                printInputs.TimeSeriesInc = round(TimeSeriesFrameInc_time / deltat);
                // Should the intermediate output be printed even if the simulation was unchanged from the previous
                // output step?
                printInputs.PrintIdleTimeSeriesFrames =
                    inputdata["Printing"]["PrintIntermediateOutput"]["PrintIdleFrames"];
                if (id == 0)
                    std::cout << "Intermediate output for movie frames will be printed every "
                              << printInputs.TimeSeriesInc << " time steps (or every "
                              << printInputs.TimeSeriesInc * deltat << " microseconds)" << std::endl;
            }
        }
        if (id == 0)
            std::cout << "Successfully parsed data printing options from input file" << std::endl;
    }

    // Print a log file for this ExaCA run in json file format, containing information about the run parameters used
    // from the input file as well as the decomposition scheme
    // Note: Passing external values for inputs like deltax that will later be stored in the grid class, with the grid
    // class passed to this function
    void PrintExaCALog(int id, int np, std::string InputFile, int ny_local, int y_offset,
                       InterfacialResponseFunction irf, double deltax, double HT_deltax, int NumberOfLayers,
                       int LayerHeight, int nx, int ny, int nz, double InitTime, double RunTime, double OutTime,
                       int cycle, double InitMaxTime, double InitMinTime, double NuclMaxTime, double NuclMinTime,
                       double CreateSVMinTime, double CreateSVMaxTime, double CaptureMaxTime, double CaptureMinTime,
                       double GhostMaxTime, double GhostMinTime, double OutMaxTime, double OutMinTime, double XMin,
                       double XMax, double YMin, double YMax, double ZMin, double ZMax, float VolFractionNucleated) {

        int *ny_local_allranks = new int[np];
        int *y_offset_allranks = new int[np];
        MPI_Gather(&ny_local, 1, MPI_INT, ny_local_allranks, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(&y_offset, 1, MPI_INT, y_offset_allranks, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (id == 0) {
            std::string FName = printInputs.PathToOutput + printInputs.BaseFileName + ".json";
            std::cout << "Printing ExaCA log file" << std::endl;
            std::ofstream ExaCALog;
            ExaCALog.open(FName);
            ExaCALog << "{" << std::endl;
            ExaCALog << "   \"ExaCAVersion\": \"" << version() << "\", " << std::endl;
            ExaCALog << "   \"ExaCACommitHash\": \"" << gitCommitHash() << "\", " << std::endl;
            ExaCALog << "   \"KokkosVersion\": \"" << kokkosVersion() << "\", " << std::endl;
            ExaCALog << "   \"InputFile\": \"" << InputFile << "\", " << std::endl;
            ExaCALog << "   \"TimeStepOfOutput\": " << cycle << "," << std::endl;
            ExaCALog << "   \"SimulationType\": \"" << SimulationType << "\"," << std::endl;
            ExaCALog << "   \"GrainOrientationFile\": \"" << GrainOrientationFile << "\"," << std::endl;
            ExaCALog << "   \"Domain\": {" << std::endl;
            ExaCALog << "      \"Nx\": " << nx << "," << std::endl;
            ExaCALog << "      \"Ny\": " << ny << "," << std::endl;
            ExaCALog << "      \"Nz\": " << nz << "," << std::endl;
            ExaCALog << "      \"CellSize\": " << deltax << "," << std::endl;
            ExaCALog << "      \"TimeStep\": " << domainInputs.deltat << "," << std::endl;
            ExaCALog << "      \"XBounds\": [" << XMin << "," << XMax << "]," << std::endl;
            ExaCALog << "      \"YBounds\": [" << YMin << "," << YMax << "]," << std::endl;
            ExaCALog << "      \"ZBounds\": [" << ZMin << "," << ZMax << "]";
            if (SimulationType != "C") {
                ExaCALog << "," << std::endl;
                ExaCALog << "      \"NumberOfLayers\": " << NumberOfLayers << "," << std::endl;
                ExaCALog << "      \"LayerOffset\": " << LayerHeight;
                if (SimulationType == "S") {
                    ExaCALog << "," << std::endl;
                    ExaCALog << "      \"NSpotsX\": " << domainInputs.NSpotsX << "," << std::endl;
                    ExaCALog << "      \"NSpotsY\": " << domainInputs.NSpotsY << "," << std::endl;
                    ExaCALog << "      \"RSpots\": " << domainInputs.SpotRadius << "," << std::endl;
                    ExaCALog << "      \"SpotOffset\": " << domainInputs.SpotOffset << std::endl;
                }
                else
                    ExaCALog << std::endl;
            }
            else
                ExaCALog << std::endl;
            ExaCALog << "   }," << std::endl;
            ExaCALog << "   \"Nucleation\": {" << std::endl;
            ExaCALog << "      \"Density\": " << nucleationInputs.NMax << "," << std::endl;
            ExaCALog << "      \"MeanUndercooling\": " << nucleationInputs.dTN << "," << std::endl;
            ExaCALog << "      \"StDevUndercooling\": " << nucleationInputs.dTsigma << "," << std::endl;
            ExaCALog << "      \"VolFractionNucleated\": " << VolFractionNucleated << std::endl;
            ExaCALog << "   }," << std::endl;
            ExaCALog << "   \"TemperatureData\": {" << std::endl;
            if (SimulationType == "R") {
                ExaCALog << "       \"TemperatureFiles\": [";
                for (int i = 0; i < temperatureInputs.TempFilesInSeries - 1; i++) {
                    ExaCALog << "\"" << temperatureInputs.temp_paths[i] << "\", ";
                }
                ExaCALog << "\"" << temperatureInputs.temp_paths[temperatureInputs.TempFilesInSeries - 1] << "\"],"
                         << std::endl;
                ExaCALog << "       \"HeatTransferCellSize\": " << HT_deltax << std::endl;
            }
            else {
                ExaCALog << "      \"G\": " << temperatureInputs.G << "," << std::endl;
                ExaCALog << "      \"R\": " << temperatureInputs.R << std::endl;
            }
            ExaCALog << "   }," << std::endl;
            ExaCALog << "   \"Substrate\": {" << std::endl;
            if (SimulationType == "C")
                ExaCALog << "       \"FractionSurfaceSitesActive\": " << substrateInputs.FractSurfaceSitesActive
                         << std::endl;
            else if (SimulationType == "SingleGrain")
                ExaCALog << "       \"GrainOrientation\": " << substrateInputs.singleGrainOrientation << std::endl;
            else {
                if (substrateInputs.UseSubstrateFile)
                    ExaCALog << "       \"SubstrateFilename\": " << substrateInputs.SubstrateFileName << std::endl;
                else
                    ExaCALog << "       \"MeanSize\": " << substrateInputs.SubstrateGrainSpacing << std::endl;
            }
            ExaCALog << "   }," << std::endl;
            ExaCALog << irf.print() << std::endl;
            ExaCALog << "   \"NumberMPIRanks\": " << np << "," << std::endl;
            ExaCALog << "   \"Decomposition\": {" << std::endl;
            ExaCALog << "       \"SubdomainYSize\": [";
            for (int i = 0; i < np - 1; i++)
                ExaCALog << ny_local_allranks[i] << ",";
            ExaCALog << ny_local_allranks[np - 1] << "]," << std::endl;
            ExaCALog << "       \"SubdomainYOffset\": [";
            for (int i = 0; i < np - 1; i++)
                ExaCALog << y_offset_allranks[i] << ",";
            ExaCALog << y_offset_allranks[np - 1] << "]" << std::endl;
            ExaCALog << "   }," << std::endl;
            ExaCALog << "   \"Timing\": {" << std::endl;
            ExaCALog << "       \"Runtime\": " << InitTime + RunTime + OutTime << "," << std::endl;
            ExaCALog << "       \"InitRunOutputBreakdown\": [" << InitTime << "," << RunTime << "," << OutTime << "],"
                     << std::endl;
            ExaCALog << "       \"MaxMinInitTime\": [" << InitMaxTime << "," << InitMinTime << "]," << std::endl;
            ExaCALog << "       \"MaxMinNucleationTime\": [" << NuclMaxTime << "," << NuclMinTime << "]," << std::endl;
            ExaCALog << "       \"MaxMinSteeringVectorCreationTime\": [" << CreateSVMaxTime << "," << CreateSVMinTime
                     << "]," << std::endl;
            ExaCALog << "       \"MaxMinCellCaptureTime\": [" << CaptureMaxTime << "," << CaptureMinTime << "],"
                     << std::endl;
            ExaCALog << "       \"MaxMinGhostExchangeTime\": [" << GhostMaxTime << "," << GhostMinTime << "],"
                     << std::endl;
            ExaCALog << "       \"MaxMinOutputTime\": [" << OutMaxTime << "," << OutMinTime << "]" << std::endl;
            ExaCALog << "   }" << std::endl;
            ExaCALog << "}" << std::endl;
            ExaCALog.close();
        }
    }
};

#endif
