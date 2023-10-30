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
    double deltax = 0.0, deltat = 0.0;
    // number of CA cells in each direction only initialized here for problem types C, S, and SingleGrain
    int nx = 0, ny = 0, nz = 0;
    // multilayer problems only
    int NumberOfLayers = 1, LayerHeight = 0;
    // problem type S only
    int NSpotsX = 0, NSpotsY = 0, SpotRadius = 0, SpotOffset = 0;
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
    int TempFilesInSeries = 0;
    double HT_deltax = 0.0;
    std::vector<std::string> temp_paths;
    // Use for problem types other than R (no temperature files to read) - default to no initial undercooling at
    // solidification front
    double G = 0, R = 0;
    double initUndercooling = 0.0;
};

struct SubstrateInputs {
    // problem type C only
    double FractSurfaceSitesActive = 0.0;
    // problem type SingleGrain only
    int singleGrainOrientation = 0;
    // problem types S and R only
    bool UseSubstrateFile = false;
    bool BaseplateThroughPowder = false;
    std::string SubstrateFileName = "";
    double SubstrateGrainSpacing = 0.0;
    // defaults to all sites in powder layer initialized with a new grain
    double PowderActiveFraction = 1.0;
    // Top of baseplate assumed at Z = 0 if not otherwise given
    double BaseplateTopZ = 0.0;
};

struct PrintInputs {
    // Base name of CA output
    std::string BaseFileName = "";
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

    // Should binary be used for printing vtk data?
    bool PrintBinary = false;

    // Should the default RVE data for ExaConstit be printed? If so, with what size?
    bool PrintDefaultRVE = false;
    int RVESize = 0;
};

struct Inputs {

    std::string SimulationType = "", MaterialFileName = "", GrainOrientationFile = "";
    double RNGSeed = 0.0;
    DomainInputs domain;
    NucleationInputs nucleation;
    TemperatureInputs temperature;
    SubstrateInputs substrate;
    PrintInputs print;

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
        domain.deltax = inputdata["Domain"]["CellSize"];
        domain.deltax = domain.deltax * pow(10, -6);
        // Time step - given in seconds, stored in microseconds
        domain.deltat = inputdata["Domain"]["TimeStep"];
        domain.deltat = domain.deltat * pow(10, -6);
        if ((SimulationType == "C") || (SimulationType == "SingleGrain")) {
            // Domain size, in cells
            domain.nx = inputdata["Domain"]["Nx"];
            domain.ny = inputdata["Domain"]["Ny"];
            domain.nz = inputdata["Domain"]["Nz"];
        }
        else {
            // Number of layers, layer height are needed for problem types S and R
            domain.NumberOfLayers = inputdata["Domain"]["NumberOfLayers"];
            domain.LayerHeight = inputdata["Domain"]["LayerOffset"];
            // Type S needs spot information, which is then used to compute the domain bounds
            if (SimulationType == "S") {
                domain.NSpotsX = inputdata["Domain"]["NSpotsX"];
                domain.NSpotsY = inputdata["Domain"]["NSpotsY"];
                // Radius and offset are given in micrometers, convert to cells
                domain.SpotRadius = inputdata["Domain"]["RSpots"];
                domain.SpotRadius = domain.SpotRadius * pow(10, -6) / domain.deltax;
                domain.SpotOffset = inputdata["Domain"]["SpotOffset"];
                domain.SpotOffset = domain.SpotOffset * pow(10, -6) / domain.deltax;
                // Calculate nx, ny, and nz based on spot array pattern and number of layers
                domain.nz = domain.SpotRadius + 1 + (domain.NumberOfLayers - 1) * domain.LayerHeight;
                domain.nx = 2 * domain.SpotRadius + 1 + domain.SpotOffset * (domain.NSpotsX - 1);
                domain.ny = 2 * domain.SpotRadius + 1 + domain.SpotOffset * (domain.NSpotsY - 1);
            }
        }

        // Nucleation inputs:
        // Nucleation density (normalized by 10^12 m^-3), mean nucleation undercooling/st dev undercooling(K)
        // SingleGrain problem type does not have nucleation, just a grain of a single orientation in the domain center
        if (SimulationType != "SingleGrain") {
            nucleation.NMax = inputdata["Nucleation"]["Density"];
            nucleation.NMax = nucleation.NMax * pow(10, 12);
            nucleation.dTN = inputdata["Nucleation"]["MeanUndercooling"];
            nucleation.dTsigma = inputdata["Nucleation"]["StDev"];
        }

        // Temperature inputs:
        if (SimulationType == "R") {
            // Temperature data resolution - default to using CA cell size if the assumed temperature data resolution if
            // not given
            if (inputdata["TemperatureData"].contains("HeatTransferCellSize")) {
                temperature.HT_deltax = inputdata["TemperatureData"]["HeatTransferCellSize"];
                // Value is given in micrometers, convert to meters
                temperature.HT_deltax = temperature.HT_deltax * pow(10, -6);
            }
            else
                temperature.HT_deltax = domain.deltax;
            if (temperature.HT_deltax != domain.deltax)
                throw std::runtime_error(
                    "Error: CA cell size and input temperature data resolution must be equivalent");
            // Read all temperature files at once (default), or one at a time?
            if (inputdata["TemperatureData"].contains("LayerwiseTempRead")) {
                temperature.LayerwiseTempRead = inputdata["TemperatureData"]["LayerwiseTempRead"];
            }
            // Get the paths/number of/names of the temperature data files used
            temperature.TempFilesInSeries = inputdata["TemperatureData"]["TemperatureFiles"].size();
            if (temperature.TempFilesInSeries == 0)
                throw std::runtime_error("Error: No temperature files listed in the temperature instructions file");
            else {
                for (int filename = 0; filename < temperature.TempFilesInSeries; filename++)
                    temperature.temp_paths.push_back(inputdata["TemperatureData"]["TemperatureFiles"][filename]);
            }
        }
        else {
            // Temperature data uses fixed thermal gradient (K/m) and cooling rate (K/s)
            temperature.G = inputdata["TemperatureData"]["G"];
            temperature.R = inputdata["TemperatureData"]["R"];
            if (SimulationType == "SingleGrain")
                temperature.initUndercooling = inputdata["TemperatureData"]["InitUndercooling"];
        }

        // Substrate inputs:
        if (SimulationType == "C") {
            // Fraction of sites at bottom surface active
            substrate.FractSurfaceSitesActive = inputdata["Substrate"]["FractionSurfaceSitesActive"];
        }
        else if (SimulationType == "SingleGrain") {
            // Orientation of the single grain at the domain center
            substrate.singleGrainOrientation = inputdata["Substrate"]["GrainOrientation"];
        }
        else {
            // Substrate data - should data come from an initial size or a file?
            if ((inputdata["Substrate"].contains("SubstrateFilename")) && (inputdata["Substrate"].contains("MeanSize")))
                throw std::runtime_error(
                    "Error: only one of substrate grain size and substrate structure filename should "
                    "be provided in the input file");
            else if (inputdata["Substrate"].contains("SubstrateFilename")) {
                substrate.SubstrateFileName = inputdata["Substrate"]["SubstrateFilename"];
                substrate.UseSubstrateFile = true;
            }
            else if (inputdata["Substrate"].contains("MeanSize")) {
                substrate.SubstrateGrainSpacing = inputdata["Substrate"]["MeanSize"];
                substrate.UseSubstrateFile = false;
            }
            // Should the baseplate microstructure be extended through the powder layers? Default is false
            if (inputdata["Substrate"].contains("ExtendSubstrateThroughPower"))
                substrate.BaseplateThroughPowder = inputdata["Substrate"]["ExtendSubstrateThroughPower"];
            // defaults to a unique grain at each site in the powder layers if not given
            if (inputdata["Substrate"].contains("PowderDensity")) {
                // powder density is given as a density per unit volume, normalized by 10^12 m^-3 --> convert this into
                // a density of sites active on the CA grid (0 to 1)
                substrate.PowderActiveFraction = inputdata["Substrate"]["PowderDensity"];
                substrate.PowderActiveFraction = substrate.PowderActiveFraction * pow(10, 12) * pow(domain.deltax, 3);
                if ((substrate.PowderActiveFraction < 0.0) || (substrate.PowderActiveFraction > 1.0))
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
                substrate.BaseplateTopZ = inputdata["Substrate"]["BaseplateTopZ"];
            else if (SimulationType == "S") {
                // Baseplate top if not otherwise given for spot array problem is at the top of the first layer
                substrate.BaseplateTopZ = domain.deltax * domain.SpotRadius;
            }
            if ((substrate.BaseplateThroughPowder) && (inputdata["Substrate"].contains("PowderDensity")))
                throw std::runtime_error("Error: if the option to extend the baseplate through the powder layers is "
                                         "toggled, a powder layer density cannot be given");
        }

        // Printing inputs
        getPrintDataFromInputFile(inputdata, id, domain.deltat);

        // Print information to console about the input file data read
        if (id == 0) {
            std::cout << "Material simulated is " << MaterialFileName << std::endl;
            std::cout << "CA cell size is " << domain.deltax * pow(10, 6) << " microns" << std::endl;
            std::cout << "Nucleation density is " << nucleation.NMax << " per m^3" << std::endl;
            std::cout << "Mean nucleation undercooling is " << nucleation.dTN
                      << " K, standard deviation of distribution is " << nucleation.dTsigma << "K" << std::endl;
            if (SimulationType == "C") {
                std::cout << "CA Simulation using a unidirectional, fixed thermal gradient of " << temperature.G
                          << " K/m and a cooling rate of " << temperature.R << " K/s" << std::endl;
                std::cout << "The time step is " << domain.deltat * pow(10, 6) << " microseconds" << std::endl;
                std::cout << "The fraction of CA cells at the bottom surface that are active is "
                          << substrate.FractSurfaceSitesActive << std::endl;
            }
            else if (SimulationType == "S") {
                std::cout << "CA Simulation using a radial, fixed thermal gradient of " << temperature.G
                          << " K/m as a series of hemispherical spots, and a cooling rate of " << temperature.R
                          << " K/s" << std::endl;
                std::cout << "A total of " << domain.NumberOfLayers << " spots per layer, with layers offset by "
                          << domain.LayerHeight << " CA cells will be simulated" << std::endl;
                std::cout << "The time step is " << domain.deltat * pow(10, 6) << " microseconds" << std::endl;
            }
            else if (SimulationType == "R") {
                std::cout << "CA Simulation using temperature data from file(s)" << std::endl;
                std::cout << "The time step is " << domain.deltat << " seconds" << std::endl;
                std::cout << "The first temperature data file to be read is " << temperature.temp_paths[0]
                          << ", and there are " << temperature.TempFilesInSeries << " in the series" << std::endl;
                std::cout << "A total of " << domain.NumberOfLayers << " layers of solidification offset by "
                          << domain.LayerHeight << " CA cells will be simulated" << std::endl;
            }
        }
    }

    // Read the input data file and initialize appropriate variables to non-default values if necessary
    void getPrintDataFromInputFile(nlohmann::json inputdata, int id, double deltat) {
        // Path to output data
        print.PathToOutput = inputdata["Printing"]["PathToOutput"];
        // Name of output data
        print.BaseFileName = inputdata["Printing"]["OutputFile"];
        // Should ASCII or binary be used to print vtk data? Defaults to ASCII if not given
        if (inputdata["Printing"].contains("PrintBinary"))
            print.PrintBinary = inputdata["Printing"]["PrintBinary"];
        // Should default ExaConstit output be printed after the simulation? If so, what size RVE?
        // If a size of 0 is given, this is set to false
        if (inputdata["Printing"].contains("PrintExaConstitSize")) {
            print.RVESize = inputdata["Printing"]["PrintExaConstitSize"];
            if (print.RVESize != 0)
                print.PrintDefaultRVE = true;
        }

        // Which fields should be printed at the start and end of the simulation?
        std::vector<std::string> InitFieldnames_key = {"GrainID", "LayerID", "MeltTimeStep", "CritTimeStep",
                                                       "UndercoolingChange"};
        std::vector<bool> PrintFieldsInit = getPrintFieldValues(inputdata, "PrintFieldsInit", InitFieldnames_key);
        if (PrintFieldsInit[0])
            print.PrintInitGrainID = true;
        if (PrintFieldsInit[1])
            print.PrintInitLayerID = true;
        if (PrintFieldsInit[2])
            print.PrintInitMeltTimeStep = true;
        if (PrintFieldsInit[3])
            print.PrintInitCritTimeStep = true;
        if (PrintFieldsInit[4])
            print.PrintInitUndercoolingChange = true;

        std::vector<std::string> FinalFieldnames_key = {
            "GrainID",      "LayerID",      "GrainMisorientation", "UndercoolingCurrent",
            "MeltTimeStep", "CritTimeStep", "UndercoolingChange",  "CellType"};
        std::vector<bool> PrintFieldsFinal = getPrintFieldValues(inputdata, "PrintFieldsFinal", FinalFieldnames_key);
        if (PrintFieldsFinal[0])
            print.PrintFinalGrainID = true;
        if (PrintFieldsFinal[1])
            print.PrintFinalLayerID = true;
        if (PrintFieldsFinal[2])
            print.PrintFinalMisorientation = true;
        if (PrintFieldsFinal[3])
            print.PrintFinalUndercoolingCurrent = true;
        if (PrintFieldsFinal[4])
            print.PrintFinalMeltTimeStep = true;
        if (PrintFieldsFinal[5])
            print.PrintFinalCritTimeStep = true;
        if (PrintFieldsFinal[6])
            print.PrintFinalUndercoolingChange = true;
        if (PrintFieldsFinal[7])
            print.PrintFinalCellType = true;

        // Should intermediate output be printed?
        if (inputdata["Printing"].contains("PrintIntermediateOutput")) {
            // An increment of 0 will set the intermediate file printing to false
            double TimeSeriesFrameInc_time = inputdata["Printing"]["PrintIntermediateOutput"]["Frequency"];
            if (TimeSeriesFrameInc_time != 0) {
                print.PrintTimeSeries = true;
                // Increment is given in microseconds, convert to seconds
                TimeSeriesFrameInc_time = TimeSeriesFrameInc_time * pow(10, -6);
                print.TimeSeriesInc = round(TimeSeriesFrameInc_time / deltat);
                // Should the intermediate output be printed even if the simulation was unchanged from the previous
                // output step?
                print.PrintIdleTimeSeriesFrames = inputdata["Printing"]["PrintIntermediateOutput"]["PrintIdleFrames"];
                if (id == 0)
                    std::cout << "Intermediate output for movie frames will be printed every " << print.TimeSeriesInc
                              << " time steps (or every " << print.TimeSeriesInc * deltat << " microseconds)"
                              << std::endl;
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
            std::string FName = print.PathToOutput + print.BaseFileName + ".json";
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
            ExaCALog << "      \"TimeStep\": " << domain.deltat << "," << std::endl;
            ExaCALog << "      \"XBounds\": [" << XMin << "," << XMax << "]," << std::endl;
            ExaCALog << "      \"YBounds\": [" << YMin << "," << YMax << "]," << std::endl;
            ExaCALog << "      \"ZBounds\": [" << ZMin << "," << ZMax << "]";
            if (SimulationType != "C") {
                ExaCALog << "," << std::endl;
                ExaCALog << "      \"NumberOfLayers\": " << NumberOfLayers << "," << std::endl;
                ExaCALog << "      \"LayerOffset\": " << LayerHeight;
                if (SimulationType == "S") {
                    ExaCALog << "," << std::endl;
                    ExaCALog << "      \"NSpotsX\": " << domain.NSpotsX << "," << std::endl;
                    ExaCALog << "      \"NSpotsY\": " << domain.NSpotsY << "," << std::endl;
                    ExaCALog << "      \"RSpots\": " << domain.SpotRadius << "," << std::endl;
                    ExaCALog << "      \"SpotOffset\": " << domain.SpotOffset << std::endl;
                }
                else
                    ExaCALog << std::endl;
            }
            else
                ExaCALog << std::endl;
            ExaCALog << "   }," << std::endl;
            ExaCALog << "   \"Nucleation\": {" << std::endl;
            ExaCALog << "      \"Density\": " << nucleation.NMax << "," << std::endl;
            ExaCALog << "      \"MeanUndercooling\": " << nucleation.dTN << "," << std::endl;
            ExaCALog << "      \"StDevUndercooling\": " << nucleation.dTsigma << "," << std::endl;
            ExaCALog << "      \"VolFractionNucleated\": " << VolFractionNucleated << std::endl;
            ExaCALog << "   }," << std::endl;
            ExaCALog << "   \"TemperatureData\": {" << std::endl;
            if (SimulationType == "R") {
                ExaCALog << "       \"TemperatureFiles\": [";
                for (int i = 0; i < temperature.TempFilesInSeries - 1; i++) {
                    ExaCALog << "\"" << temperature.temp_paths[i] << "\", ";
                }
                ExaCALog << "\"" << temperature.temp_paths[temperature.TempFilesInSeries - 1] << "\"]," << std::endl;
                ExaCALog << "       \"HeatTransferCellSize\": " << HT_deltax << std::endl;
            }
            else {
                ExaCALog << "      \"G\": " << temperature.G << "," << std::endl;
                ExaCALog << "      \"R\": " << temperature.R << std::endl;
            }
            ExaCALog << "   }," << std::endl;
            ExaCALog << "   \"Substrate\": {" << std::endl;
            if (SimulationType == "C")
                ExaCALog << "       \"FractionSurfaceSitesActive\": " << substrate.FractSurfaceSitesActive << std::endl;
            else if (SimulationType == "SingleGrain")
                ExaCALog << "       \"GrainOrientation\": " << substrate.singleGrainOrientation << std::endl;
            else {
                if (substrate.UseSubstrateFile)
                    ExaCALog << "       \"SubstrateFilename\": " << substrate.SubstrateFileName << std::endl;
                else
                    ExaCALog << "       \"MeanSize\": " << substrate.SubstrateGrainSpacing << std::endl;
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
