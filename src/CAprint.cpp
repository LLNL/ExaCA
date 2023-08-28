// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT
#include "CAprint.hpp"
#include "CAfunctions.hpp"
#include "CAinterfacialresponse.hpp"
#include "CAtypes.hpp"
#include "mpi.h"

#include <Kokkos_Core.hpp>

#include <nlohmann/json.hpp>

std::string version() { return ExaCA_VERSION; }

std::string gitCommitHash() { return ExaCA_GIT_COMMIT_HASH; }

std::string kokkosVersion() { return ExaCA_Kokkos_VERSION_STRING; }

// Print a log file for this ExaCA run in json file format, containing information about the run parameters used
// from the input file as well as the decomposition scheme
void PrintExaCALog(int id, int np, std::string InputFile, std::string PathToOutput, std::string BaseFileName,
                   std::string SimulationType, int ny_local, int y_offset, InterfacialResponseFunction irf,
                   double deltax, double NMax, double dTN, double dTsigma, std::vector<std::string> temp_paths,
                   int TempFilesInSeries, double HT_deltax, double deltat, int NumberOfLayers, int LayerHeight,
                   std::string SubstrateFileName, double SubstrateGrainSpacing, bool SubstrateFile, double G, double R,
                   int nx, int ny, int nz, double FractSurfaceSitesActive, int NSpotsX, int NSpotsY, int SpotOffset,
                   int SpotRadius, double InitTime, double RunTime, double OutTime, int cycle, double InitMaxTime,
                   double InitMinTime, double NuclMaxTime, double NuclMinTime, double CreateSVMinTime,
                   double CreateSVMaxTime, double CaptureMaxTime, double CaptureMinTime, double GhostMaxTime,
                   double GhostMinTime, double OutMaxTime, double OutMinTime, double XMin, double XMax, double YMin,
                   double YMax, double ZMin, double ZMax, std::string GrainOrientationFile, float VolFractionNucleated,
                   int singleGrainOrientation) {

    int *ny_local_allranks = new int[np];
    int *y_offset_allranks = new int[np];
    MPI_Gather(&ny_local, 1, MPI_INT, ny_local_allranks, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&y_offset, 1, MPI_INT, y_offset_allranks, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (id == 0) {
        std::string FName = PathToOutput + BaseFileName + ".json";
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
        ExaCALog << "      \"TimeStep\": " << deltat << "," << std::endl;
        ExaCALog << "      \"XBounds\": [" << XMin << "," << XMax << "]," << std::endl;
        ExaCALog << "      \"YBounds\": [" << YMin << "," << YMax << "]," << std::endl;
        ExaCALog << "      \"ZBounds\": [" << ZMin << "," << ZMax << "]";
        if (SimulationType != "C") {
            ExaCALog << "," << std::endl;
            ExaCALog << "      \"NumberOfLayers\": " << NumberOfLayers << "," << std::endl;
            ExaCALog << "      \"LayerOffset\": " << LayerHeight;
            if (SimulationType == "S") {
                ExaCALog << "," << std::endl;
                ExaCALog << "      \"NSpotsX\": " << NSpotsX << "," << std::endl;
                ExaCALog << "      \"NSpotsY\": " << NSpotsY << "," << std::endl;
                ExaCALog << "      \"RSpots\": " << SpotRadius << "," << std::endl;
                ExaCALog << "      \"SpotOffset\": " << SpotOffset << std::endl;
            }
            else
                ExaCALog << std::endl;
        }
        else
            ExaCALog << std::endl;
        ExaCALog << "   }," << std::endl;
        ExaCALog << "   \"Nucleation\": {" << std::endl;
        ExaCALog << "      \"Density\": " << NMax << "," << std::endl;
        ExaCALog << "      \"MeanUndercooling\": " << dTN << "," << std::endl;
        ExaCALog << "      \"StDevUndercooling\": " << dTsigma << "," << std::endl;
        ExaCALog << "      \"VolFractionNucleated\": " << VolFractionNucleated << std::endl;
        ExaCALog << "   }," << std::endl;
        ExaCALog << "   \"TemperatureData\": {" << std::endl;
        if (SimulationType == "R") {
            ExaCALog << "       \"TemperatureFiles\": [";
            for (int i = 0; i < TempFilesInSeries - 1; i++) {
                ExaCALog << "\"" << temp_paths[i] << "\", ";
            }
            ExaCALog << "\"" << temp_paths[TempFilesInSeries - 1] << "\"]," << std::endl;
            ExaCALog << "       \"HeatTransferCellSize\": " << HT_deltax << std::endl;
        }
        else {
            ExaCALog << "      \"G\": " << G << "," << std::endl;
            ExaCALog << "      \"R\": " << R << std::endl;
        }
        ExaCALog << "   }," << std::endl;
        ExaCALog << "   \"Substrate\": {" << std::endl;
        if (SimulationType == "C")
            ExaCALog << "       \"FractionSurfaceSitesActive\": " << FractSurfaceSitesActive << std::endl;
        else if (SimulationType == "SingleGrain")
            ExaCALog << "       \"GrainOrientation\": " << singleGrainOrientation << std::endl;
        else {
            if (SubstrateFile)
                ExaCALog << "       \"SubstrateFilename\": " << SubstrateFileName << std::endl;
            else
                ExaCALog << "       \"MeanSize\": " << SubstrateGrainSpacing << std::endl;
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
        ExaCALog << "       \"MaxMinGhostExchangeTime\": [" << GhostMaxTime << "," << GhostMinTime << "]," << std::endl;
        ExaCALog << "       \"MaxMinOutputTime\": [" << OutMaxTime << "," << OutMinTime << "]" << std::endl;
        ExaCALog << "   }" << std::endl;
        ExaCALog << "}" << std::endl;
        ExaCALog.close();
    }
}

// Print timing info to console
void PrintExaCATiming(int np, double InitTime, double RunTime, double OutTime, int cycle, double InitMaxTime,
                      double InitMinTime, double NuclMaxTime, double NuclMinTime, double CreateSVMinTime,
                      double CreateSVMaxTime, double CaptureMaxTime, double CaptureMinTime, double GhostMaxTime,
                      double GhostMinTime, double OutMaxTime, double OutMinTime) {

    std::cout << "===================================================================================" << std::endl;
    std::cout << "Having run with = " << np << " processors" << std::endl;
    std::cout << "Output written at cycle = " << cycle << std::endl;
    std::cout << "Total time = " << InitTime + RunTime + OutTime << std::endl;
    std::cout << "Time spent initializing data = " << InitTime << " s" << std::endl;
    std::cout << "Time spent performing CA calculations = " << RunTime << " s" << std::endl;
    std::cout << "Time spent collecting and printing output data = " << OutTime << " s\n" << std::endl;

    std::cout << "Max/min rank time initializing data  = " << InitMaxTime << " / " << InitMinTime << " s" << std::endl;
    std::cout << "Max/min rank time in CA nucleation   = " << NuclMaxTime << " / " << NuclMinTime << " s" << std::endl;
    std::cout << "Max/min rank time in CA steering vector creation = " << CreateSVMaxTime << " / " << CreateSVMinTime
              << " s" << std::endl;
    std::cout << "Max/min rank time in CA cell capture = " << CaptureMaxTime << " / " << CaptureMinTime << " s"
              << std::endl;
    std::cout << "Max/min rank time in CA ghosting     = " << GhostMaxTime << " / " << GhostMinTime << " s"
              << std::endl;
    std::cout << "Max/min rank time exporting data     = " << OutMaxTime << " / " << OutMinTime << " s\n" << std::endl;

    std::cout << "===================================================================================" << std::endl;
}
