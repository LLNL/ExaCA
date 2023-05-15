// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAparsefiles.hpp"
#include "GAprint.hpp"
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
        ViewI3D_H GrainID(Kokkos::ViewAllocateWithoutInitializing("GrainID"), nz, nx, ny);
        ViewI3D_H LayerID(Kokkos::ViewAllocateWithoutInitializing("LayerID"), nz, nx, ny);

        // Fill arrays with data from paraview file
        InitializeData(MicrostructureFile, nx, ny, nz, GrainID, LayerID);

        // Grain unit vectors, grain euler angles, RGB colors for IPF-Z coloring
        // (9*NumberOfOrientations,  3*NumberOfOrientations, and 3*NumberOfOrientations in size, respectively)
        int NumberOfOrientations = 0;
        ViewF GrainUnitVector(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), 9 * NumberOfOrientations);
        ViewF GrainEulerAngles(Kokkos::ViewAllocateWithoutInitializing("GrainEulerAngles"), 3 * NumberOfOrientations);
        ViewF GrainRGBValues(Kokkos::ViewAllocateWithoutInitializing("GrainRGBValues"), 3 * NumberOfOrientations);

        // Initialize, then copy back to host
        OrientationInit(0, NumberOfOrientations, GrainUnitVector, RotationFilename, 9);
        ViewF_H GrainUnitVector_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainUnitVector);
        OrientationInit(0, NumberOfOrientations, GrainEulerAngles, EulerAnglesFilename, 3);
        ViewF_H GrainEulerAngles_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainEulerAngles);
        OrientationInit(0, NumberOfOrientations, GrainRGBValues, RGBFilename, 3);
        ViewF_H GrainRGBValues_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainRGBValues);

        // Representative region creation
        std::ifstream AnalysisDataStream(AnalysisFile);
        nlohmann::json AnalysisData = nlohmann::json::parse(AnalysisDataStream);
        nlohmann::json RegionsData = AnalysisData["Regions"];
        int NumberOfRegions = AnalysisData["Regions"].size();
        std::cout << "There are " << NumberOfRegions << " regions to analyze" << std::endl;
        for (auto it = RegionsData.begin(); it != RegionsData.end(); it++) {
            // Create region
            std::string RegionName = it.key();
            std::cout << "Parsing data for region " << RegionName << std::endl;
            nlohmann::json RegionData = AnalysisData["Regions"][RegionName];
            RepresentativeRegion representativeRegion(RegionData, nx, ny, nz, deltax, XYZBounds);
            std::cout << "Loaded analysis options for region " << RegionName << std::endl;
            std::string BaseFileNameThisRegion = BaseFileName + "_" + RegionName;
            // List of grain ID values in the representative region
            std::vector<int> GrainIDVector = representativeRegion.getGrainIDVector(GrainID);
            // List of unique grain IDs associated with the region
            std::vector<int> UniqueGrainIDVector = representativeRegion.getUniqueGrainIDVector(GrainIDVector);
            int NumberOfGrains = UniqueGrainIDVector.size();

            // Get the size (in units of length, area, or volume) associated with each unique grain ID value
            // TODO: Once Json is required, move all options calculating aspects of the representative region into
            // the struct and out of GAutils, deleting the redundant GAutils subroutines
            std::vector<float> GrainSizeVector_Microns =
                representativeRegion.getGrainSizeVector(GrainIDVector, UniqueGrainIDVector, NumberOfGrains, deltax);

            // Output file stream for quantities of interest
            std::ofstream QoIs;
            std::string QoIs_fname = BaseFileNameThisRegion + "_QoIs.txt";
            QoIs.open(QoIs_fname);
            // Header data for QoIs file (probably just pass the whole class to the subroutines in a future release
            // once Json is required)
            // TODO: After removing non-JSON input option, pass entire object to printAnaylsisHeader routines
            if (representativeRegion.regionType == "area")
                printAnalysisHeader_Area(QoIs, representativeRegion.xBounds_Cells[0],
                                         representativeRegion.xBounds_Cells[1], representativeRegion.yBounds_Cells[0],
                                         representativeRegion.yBounds_Cells[1], representativeRegion.zBounds_Cells[0],
                                         representativeRegion.zBounds_Cells[1], representativeRegion.xBounds_Meters[0],
                                         representativeRegion.xBounds_Meters[1], representativeRegion.yBounds_Meters[0],
                                         representativeRegion.yBounds_Meters[1], representativeRegion.zBounds_Meters[0],
                                         representativeRegion.zBounds_Meters[1], RegionName,
                                         representativeRegion.regionOrientation);
            else if (representativeRegion.regionType == "volume")
                representativeRegion.printAnalysisHeader_Volume(QoIs, RegionName);

            // TODO: After removing non-JSON input option, pass entire object to printGrainTypeFractions routine
            // Fraction of region consisting of nucleated grains, unmelted material
            if (representativeRegion.AnalysisOptions_StatsYN[0]) {
                printGrainTypeFractions(QoIs, representativeRegion.xBounds_Cells[0],
                                        representativeRegion.xBounds_Cells[1], representativeRegion.yBounds_Cells[0],
                                        representativeRegion.yBounds_Cells[1], representativeRegion.zBounds_Cells[0],
                                        representativeRegion.zBounds_Cells[1], GrainID, LayerID,
                                        representativeRegion.regionSize_Cells);
            }

            // Calculate and if specified, print misorientation data
            std::vector<float> GrainMisorientationXVector = getGrainMisorientation(
                "X", GrainUnitVector_Host, UniqueGrainIDVector, NumberOfOrientations, NumberOfGrains);
            std::vector<float> GrainMisorientationYVector = getGrainMisorientation(
                "Y", GrainUnitVector_Host, UniqueGrainIDVector, NumberOfOrientations, NumberOfGrains);
            std::vector<float> GrainMisorientationZVector = getGrainMisorientation(
                "Z", GrainUnitVector_Host, UniqueGrainIDVector, NumberOfOrientations, NumberOfGrains);
            if (representativeRegion.AnalysisOptions_StatsYN[1]) {
                printMeanMisorientations(QoIs, NumberOfGrains, GrainMisorientationXVector, GrainMisorientationYVector,
                                         GrainMisorientationZVector, GrainSizeVector_Microns,
                                         representativeRegion.regionSize_Meters);
            }

            // Print mean size data if specified
            if (representativeRegion.AnalysisOptions_StatsYN[2])
                printMeanSize(QoIs, NumberOfGrains, representativeRegion.regionSize_Microns,
                              representativeRegion.regionType, representativeRegion.units_dimension);

            // If XExtent, YExtent, ZExtent, or BuildTransAspectRatio/Extent are toggled for general stats printing
            // or per grain printing, calculate grain extents for the necessary direction(s) (otherwise don't, since
            // it can be slow for large volumes) Extents are calculated in microns
            // TODO: After removing non-JSON input option, pass entire object to calcGrainExtent routine
            std::vector<float> GrainExtentX(NumberOfGrains);
            std::vector<float> GrainExtentY(NumberOfGrains);
            std::vector<float> GrainExtentZ(NumberOfGrains);
            std::vector<float> BuildTransAspectRatio(NumberOfGrains);
            if ((representativeRegion.AnalysisOptions_StatsYN[3]) ||
                (representativeRegion.AnalysisOptions_StatsYN[4]) ||
                (representativeRegion.AnalysisOptions_PerGrainStatsYN[2]) ||
                (representativeRegion.AnalysisOptions_PerGrainStatsYN[5]))
                calcGrainExtent(GrainExtentX, GrainID, UniqueGrainIDVector, GrainSizeVector_Microns, NumberOfGrains,
                                representativeRegion.xBounds_Cells[0], representativeRegion.xBounds_Cells[1],
                                representativeRegion.yBounds_Cells[0], representativeRegion.yBounds_Cells[1],
                                representativeRegion.zBounds_Cells[0], representativeRegion.zBounds_Cells[1], "X",
                                deltax, representativeRegion.regionType);
            if ((representativeRegion.AnalysisOptions_StatsYN[3]) ||
                (representativeRegion.AnalysisOptions_StatsYN[5]) ||
                (representativeRegion.AnalysisOptions_PerGrainStatsYN[3]) ||
                (representativeRegion.AnalysisOptions_PerGrainStatsYN[5]))
                calcGrainExtent(GrainExtentY, GrainID, UniqueGrainIDVector, GrainSizeVector_Microns, NumberOfGrains,
                                representativeRegion.xBounds_Cells[0], representativeRegion.xBounds_Cells[1],
                                representativeRegion.yBounds_Cells[0], representativeRegion.yBounds_Cells[1],
                                representativeRegion.zBounds_Cells[0], representativeRegion.zBounds_Cells[1], "Y",
                                deltax, representativeRegion.regionType);
            if ((representativeRegion.AnalysisOptions_StatsYN[3]) ||
                (representativeRegion.AnalysisOptions_StatsYN[6]) ||
                (representativeRegion.AnalysisOptions_PerGrainStatsYN[4]) ||
                (representativeRegion.AnalysisOptions_PerGrainStatsYN[5]))
                calcGrainExtent(GrainExtentZ, GrainID, UniqueGrainIDVector, GrainSizeVector_Microns, NumberOfGrains,
                                representativeRegion.xBounds_Cells[0], representativeRegion.xBounds_Cells[1],
                                representativeRegion.yBounds_Cells[0], representativeRegion.yBounds_Cells[1],
                                representativeRegion.zBounds_Cells[0], representativeRegion.zBounds_Cells[1], "Z",
                                deltax, representativeRegion.regionType);
            if ((representativeRegion.AnalysisOptions_StatsYN[3]) ||
                (representativeRegion.AnalysisOptions_PerGrainStatsYN[5]))
                calcBuildTransAspectRatio(BuildTransAspectRatio, GrainExtentX, GrainExtentY, GrainExtentZ,
                                          NumberOfGrains);
            if (representativeRegion.AnalysisOptions_StatsYN[3])
                printMeanBuildTransAspectRatio(QoIs, GrainExtentX, GrainExtentY, GrainExtentZ, GrainSizeVector_Microns,
                                               representativeRegion.regionSize_Microns, NumberOfGrains);
            if (representativeRegion.AnalysisOptions_StatsYN[4]) {
                printMeanExtent(QoIs, GrainExtentX, "X", NumberOfGrains);
            }
            if (representativeRegion.AnalysisOptions_StatsYN[5])
                printMeanExtent(QoIs, GrainExtentY, "Y", NumberOfGrains);
            if (representativeRegion.AnalysisOptions_StatsYN[6])
                printMeanExtent(QoIs, GrainExtentZ, "Z", NumberOfGrains);

            // Determine IPF-Z color of each grain relative to each direction: 0 (red), 1 (green), 2 (blue)
            std::vector<float> GrainRed =
                getIPFZColor(0, UniqueGrainIDVector, NumberOfOrientations, GrainRGBValues_Host, NumberOfGrains);
            std::vector<float> GrainGreen =
                getIPFZColor(1, UniqueGrainIDVector, NumberOfOrientations, GrainRGBValues_Host, NumberOfGrains);
            std::vector<float> GrainBlue =
                getIPFZColor(1, UniqueGrainIDVector, NumberOfOrientations, GrainRGBValues_Host, NumberOfGrains);

            // Write grain area data as a function of Z location in the representative volume if the options are
            // toggled, writing to files
            // "[BaseFileNameThisRegion]_WeightedGrainAreas.csv" and "[BaseFileNameThisRegion]_GrainAreas.csv",
            // respectively
            // TODO: After removing non-JSON input option, pass entire object to writeAreaSeries routine
            if ((representativeRegion.AnalysisOptions_LayerwiseStatsYN[0]) ||
                (representativeRegion.AnalysisOptions_LayerwiseStatsYN[1]))
                writeAreaSeries(representativeRegion.AnalysisOptions_LayerwiseStatsYN[1],
                                representativeRegion.AnalysisOptions_LayerwiseStatsYN[0], BaseFileNameThisRegion,
                                deltax, representativeRegion.xBounds_Cells[0], representativeRegion.xBounds_Cells[1],
                                representativeRegion.yBounds_Cells[0], representativeRegion.yBounds_Cells[1],
                                representativeRegion.zBounds_Cells[0], representativeRegion.zBounds_Cells[1], GrainID,
                                representativeRegion.zBounds_Meters[0]);
            QoIs.close();

            // Write per-grain stats for the analysis types specified to the file
            // "[BaseFileNameThisRegion]_grains.csv"
            if (representativeRegion.PrintPerGrainStatsYN)
                representativeRegion.writePerGrainStats(
                    BaseFileNameThisRegion, UniqueGrainIDVector, GrainMisorientationXVector, GrainMisorientationYVector,
                    GrainMisorientationZVector, GrainSizeVector_Microns, GrainExtentX, GrainExtentY, GrainExtentZ,
                    BuildTransAspectRatio, NumberOfGrains, GrainRed, GrainGreen, GrainBlue);

            // ExaConstit print a file named "[BaseFileNameThisRegion]_ExaConstit.csv"
            if (representativeRegion.PrintExaConstitYN) {
                representativeRegion.writeExaConstitRVE(BaseFileNameThisRegion, deltax, GrainID);
            }

            // Pole figure print a file named "[BaseFileNameThisRegion]_PoleFigureData.txt"
            if (representativeRegion.PrintPoleFigureYN) {
                ViewI_H GOHistogram = getOrientationHistogram(
                    NumberOfOrientations, GrainID, LayerID, representativeRegion.xBounds_Cells[0],
                    representativeRegion.xBounds_Cells[1], representativeRegion.yBounds_Cells[0],
                    representativeRegion.yBounds_Cells[1], representativeRegion.zBounds_Cells[0],
                    representativeRegion.zBounds_Cells[1]);
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
