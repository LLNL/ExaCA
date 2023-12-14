// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAparsefiles.hpp"
#include "GArepresentativeregion.hpp"
#include "GAutils.hpp"

#include "ExaCA.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// The path to/name of the ExaCA analysis input file and the path to/ name of the base filename for the ExaCA
// microstructure data (without extension) are given on the command line
int main(int argc, char *argv[]) {

    // Initialize Kokkos
    Kokkos::initialize();
    {
        using memory_space = typename Kokkos::DefaultExecutionSpace::memory_space;
        // Read command line input to obtain name of analysis file
        std::string AnalysisFile;
        if (argc < 2) {
            throw std::runtime_error("Error: Full path to and name of analysis file must be given on the command line");
        }
        AnalysisFile = argv[1];
        std::string BaseFileName = argv[2];
        std::string LogFile = BaseFileName + ".json";
        std::string MicrostructureFile = BaseFileName + ".vtk";
        std::cout << "Performing analysis of " << MicrostructureFile << " , using the log file " << LogFile
                  << " and the options specified in " << AnalysisFile << std::endl;

        std::string RotationFilename, EulerAnglesFilename, RGBFilename;
        double deltax;
        int nx, ny, nz, NumberOfLayers;
        std::vector<double> XYZBounds(6);
        ParseLogFile(LogFile, nx, ny, nz, deltax, NumberOfLayers, XYZBounds, RotationFilename, EulerAnglesFilename,
                     RGBFilename, false);

        // Allocate memory blocks for GrainID and LayerID data
        Kokkos::View<int ***, Kokkos::HostSpace> GrainID(Kokkos::ViewAllocateWithoutInitializing("GrainID"), nz, nx,
                                                         ny);
        Kokkos::View<short ***, Kokkos::HostSpace> LayerID(Kokkos::ViewAllocateWithoutInitializing("LayerID"), nz, nx,
                                                           ny);

        // Fill arrays with data from paraview file
        InitializeData(MicrostructureFile, nx, ny, nz, GrainID, LayerID);

        // Grain unit vectors, grain euler angles, RGB colors for IPF-Z coloring
        // (9*NumberOfOrientations,  3*NumberOfOrientations, and 3*NumberOfOrientations in size, respectively)
        int NumberOfOrientations = 0;
        using view_float = Kokkos::View<float *, memory_space>;
        view_float GrainUnitVector(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"),
                                   9 * NumberOfOrientations);
        view_float GrainEulerAngles(Kokkos::ViewAllocateWithoutInitializing("GrainEulerAngles"),
                                    3 * NumberOfOrientations);
        view_float GrainRGBValues(Kokkos::ViewAllocateWithoutInitializing("GrainRGBValues"), 3 * NumberOfOrientations);

        // Initialize, then copy back to host
        OrientationInit(0, NumberOfOrientations, GrainUnitVector, RotationFilename, 9);
        auto GrainUnitVector_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainUnitVector);
        OrientationInit(0, NumberOfOrientations, GrainEulerAngles, EulerAnglesFilename, 3);
        auto GrainEulerAngles_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainEulerAngles);
        OrientationInit(0, NumberOfOrientations, GrainRGBValues, RGBFilename, 3);
        auto GrainRGBValues_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainRGBValues);

        // Representative region creation
        std::ifstream AnalysisDataStream(AnalysisFile);
        nlohmann::json AnalysisData = nlohmann::json::parse(AnalysisDataStream);
        nlohmann::json RegionsData = AnalysisData["Regions"];
        int NumberOfRegions = AnalysisData["Regions"].size();
        std::cout << "There are " << NumberOfRegions << " regions to analyze" << std::endl;
        for (auto it = RegionsData.begin(); it != RegionsData.end(); it++) {
            // Create region
            std::string RegionName = it.key();
            RepresentativeRegion representativeRegion(AnalysisData, RegionName, nx, ny, nz, deltax, XYZBounds, GrainID);
            std::string BaseFileNameThisRegion = BaseFileName + "_" + RegionName;

            // Output file stream for quantities of interest
            std::ofstream QoIs;
            std::string QoIs_fname = BaseFileNameThisRegion + "_QoIs.txt";
            QoIs.open(QoIs_fname);
            // Header data for QoIs file
            representativeRegion.printAnalysisHeader(QoIs);

            // Fraction of region consisting of nucleated grains, unmelted material
            if (representativeRegion.AnalysisOptions_StatsYN[0])
                representativeRegion.printGrainTypeFractions(QoIs, GrainID, LayerID);

            // Calculate and if specified, print misorientation data
            std::vector<float> GrainMisorientationXVector =
                representativeRegion.getGrainMisorientation("X", GrainUnitVector_Host, NumberOfOrientations);
            std::vector<float> GrainMisorientationYVector =
                representativeRegion.getGrainMisorientation("Y", GrainUnitVector_Host, NumberOfOrientations);
            std::vector<float> GrainMisorientationZVector =
                representativeRegion.getGrainMisorientation("Z", GrainUnitVector_Host, NumberOfOrientations);
            if (representativeRegion.AnalysisOptions_StatsYN[1])
                representativeRegion.printMeanMisorientations(QoIs, GrainMisorientationXVector,
                                                              GrainMisorientationYVector, GrainMisorientationZVector);

            // Print mean size data if specified
            if (representativeRegion.AnalysisOptions_StatsYN[2])
                representativeRegion.printMeanSize(QoIs);

            // If XExtent, YExtent, ZExtent, or BuildTransAspectRatio/Extent are toggled for general stats printing
            // or per grain printing, calculate grain extents for the necessary direction(s) (otherwise don't, since
            // it can be slow for large volumes) Extents are calculated in microns
            // If options StatsYN[3] or PerGrainStatsYN[5] is toggled, grain extents are needed
            // If StatsYN[4] or PerGrainStatsYN[2] is toggled, X extents are needed
            // If StatsYN[5] or PerGrainStatsYN[3] is toggled, Y extents are needed
            // If StatsYN[6] or PerGrainStatsYN[4] is toggled, Z extents are needed
            representativeRegion.calcNecessaryGrainExtents(GrainID, deltax);
            std::vector<float> BuildTransAspectRatio(representativeRegion.NumberOfGrains);
            if ((representativeRegion.AnalysisOptions_StatsYN[3]) ||
                (representativeRegion.AnalysisOptions_PerGrainStatsYN[5]))
                representativeRegion.calcBuildTransAspectRatio(BuildTransAspectRatio);
            if (representativeRegion.AnalysisOptions_StatsYN[3])
                representativeRegion.printMeanBuildTransAspectRatio(QoIs);

            if (representativeRegion.AnalysisOptions_StatsYN[4])
                representativeRegion.printMeanExtent(QoIs, "X");
            if (representativeRegion.AnalysisOptions_StatsYN[5])
                representativeRegion.printMeanExtent(QoIs, "Y");
            if (representativeRegion.AnalysisOptions_StatsYN[6])
                representativeRegion.printMeanExtent(QoIs, "Z");

            // Determine IPF-Z color of each grain relative to each direction: 0 (red), 1 (green), 2 (blue)
            std::vector<float> GrainRed =
                representativeRegion.getIPFZColor(0, NumberOfOrientations, GrainRGBValues_Host);
            std::vector<float> GrainGreen =
                representativeRegion.getIPFZColor(1, NumberOfOrientations, GrainRGBValues_Host);
            std::vector<float> GrainBlue =
                representativeRegion.getIPFZColor(1, NumberOfOrientations, GrainRGBValues_Host);

            // Write grain area data as a function of Z location in the representative volume if the options are
            // toggled, writing to files
            // "[BaseFileNameThisRegion]_WeightedGrainAreas.csv" and "[BaseFileNameThisRegion]_GrainAreas.csv",
            // respectively
            if ((representativeRegion.AnalysisOptions_LayerwiseStatsYN[0]) ||
                (representativeRegion.AnalysisOptions_LayerwiseStatsYN[1]))
                representativeRegion.writeAreaSeries(BaseFileNameThisRegion, deltax, GrainID);
            QoIs.close();

            // Write per-grain stats for the analysis types specified to the file
            // "[BaseFileNameThisRegion]_grains.csv"
            if (representativeRegion.PrintPerGrainStatsYN)
                representativeRegion.writePerGrainStats(BaseFileNameThisRegion, GrainMisorientationXVector,
                                                        GrainMisorientationYVector, GrainMisorientationZVector,
                                                        BuildTransAspectRatio, GrainRed, GrainGreen, GrainBlue);

            // ExaConstit print a file named "[BaseFileNameThisRegion]_ExaConstit.csv"
            if (representativeRegion.PrintExaConstitYN) {
                representativeRegion.writeExaConstitRVE(BaseFileNameThisRegion, deltax, GrainID);
            }

            // Pole figure print a file named "[BaseFileNameThisRegion]_PoleFigureData.txt"
            if (representativeRegion.PrintPoleFigureYN) {
                auto GOHistogram = representativeRegion.getOrientationHistogram(NumberOfOrientations, GrainID, LayerID);
                representativeRegion.writePoleFigure(BaseFileNameThisRegion, NumberOfOrientations,
                                                     GrainEulerAngles_Host, GOHistogram);
            }

            // IPF map for area print a file named "[BaseFileNameThisRegion]_IPFCrossSectionData.txt"
            if (representativeRegion.PrintInversePoleFigureMapYN) {
                representativeRegion.writeIPFColoredCrossSection(BaseFileNameThisRegion, GrainID, GrainEulerAngles_Host,
                                                                 deltax, NumberOfOrientations);
            }
            std::cout << "Finished analysis for region " << RegionName << std::endl;
        } // end loop over all representative regions in analysis file
    }     // end scope for kokkos
    // Finalize kokkos and end program
    Kokkos::finalize();
    return 0;
}
