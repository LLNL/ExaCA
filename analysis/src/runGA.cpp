// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

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

// The name of the ExaCA output analysis file (.txt, located in ExaCA/analysis), the ExaCA microstructure file (.vtk),
// and the ExaCA log file (.log) is given on the command line - these files should all have the same name, other than
// the file extension
int main(int argc, char *argv[]) {

    // Read command line input to obtain name of analysis file
    bool OrientationFilesInInput;
    std::string AnalysisFile, LogFile, MicrostructureFile, RotationFilename, EulerAnglesFilename, OutputFileName,
        RGBFilename;
    double deltax;
    if (argc < 2) {
        throw std::runtime_error("Error: Full path to and name of analysis file must be given on the command line");
    }
    else {
        AnalysisFile = argv[1];
        // Get path to/name of all files of interest by reading the analysis file
        // If the file of orientations is given in both rotation matrix and Euler angle form, it is assumed that the new
        // output format is used for cross-section orientation data Otherwise, the old format is used
        ParseFilenames(AnalysisFile, LogFile, MicrostructureFile, RotationFilename, OutputFileName, EulerAnglesFilename,
                       RGBFilename, OrientationFilesInInput);
    }
    std::cout << "Performing analysis of " << MicrostructureFile << " , using the log file " << LogFile
              << " and the options specified in " << AnalysisFile << std::endl;

    // Initialize Kokkos
    Kokkos::initialize();
    {

        // Given the input file name, parse the paraview file for the cell size, x, y, and z dimensions, number of
        // layers
        int nx, ny, nz, NumberOfLayers;
        std::vector<double> XYZBounds(6);
        int LogFormat = checkLogFormat(LogFile);
        if (LogFormat == 0)
            ParseLogFile_OldNoColon(LogFile, nx, ny, nz, deltax, NumberOfLayers, true, XYZBounds);
        else if (LogFormat == 1)
            ParseLogFile_Old(LogFile, nx, ny, nz, deltax, NumberOfLayers, XYZBounds, RotationFilename,
                             EulerAnglesFilename, RGBFilename, OrientationFilesInInput);
        else {
#ifdef ExaCA_ENABLE_JSON
            ParseLogFile(LogFile, nx, ny, nz, deltax, NumberOfLayers, XYZBounds, RotationFilename, EulerAnglesFilename,
                         RGBFilename, OrientationFilesInInput);
#endif
        }

        // Allocate memory blocks for GrainID and LayerID data
        ViewI3D_H GrainID(Kokkos::ViewAllocateWithoutInitializing("GrainID"), nz, nx, ny);
        ViewI3D_H LayerID(Kokkos::ViewAllocateWithoutInitializing("LayerID"), nz, nx, ny);

        // Fill arrays with data from paraview file
        InitializeData(MicrostructureFile, nx, ny, nz, GrainID, LayerID);

        // Read analysis file ("ExaCA/examples/Outputs.txt") to determine which analysis should be done
        // There are three parts to this analysis:
        // Part 1: Determine the number of, and bounds of, specified RVEs for ExaConstits
        int NumberOfRVEs = 0;
        std::vector<int> XLow_RVE, XHigh_RVE, YLow_RVE, YHigh_RVE, ZLow_RVE,
            ZHigh_RVE; // Contains data on each individual RVE bounds

        // Part 2: Determine the number of, and planes for, cross-sections
        int NumberOfCrossSections = 0;
        std::vector<std::string>
            CrossSectionPlane;             // Contains type of cross-section/location of each cross-section as a string
        std::vector<std::string> CSLabels; // Labels for each cross-section when printed to files
        std::vector<int> CrossSectionLocation; // Contains location of each cross-section as an integer
        std::vector<bool> PrintSectionPF, PrintSectionIPF,
            BimodalAnalysis; // Whether or not to print pole figure or IPF-colormap data
                             // for each individual cross section, whether or not to separate grain areas into two
                             // distributions or not
        // Part 3: Analysis options on a representative region in x, y, and z (can be different than the x, y, and z for
        // ExaConstit)
        bool AnalysisTypes[8] = {0};            // which analysis modes other than the defaults should be considered?
        int XMin, XMax, YMin, YMax, ZMin, ZMax; // bounds of analysis region

        // Analysis also requires reading orientation files of specified names
        int NumberOfOrientations;
        ParseAnalysisFile(AnalysisFile, RotationFilename, NumberOfOrientations, AnalysisTypes, XLow_RVE, XHigh_RVE,
                          YLow_RVE, YHigh_RVE, ZLow_RVE, ZHigh_RVE, NumberOfRVEs, CrossSectionPlane,
                          CrossSectionLocation, NumberOfCrossSections, XMin, XMax, YMin, YMax, ZMin, ZMax, nx, ny, nz,
                          LayerID, NumberOfLayers, PrintSectionPF, PrintSectionIPF, BimodalAnalysis, CSLabels);

        // Allocate memory for grain unit vectors, grain euler angles, RGB colors for IPF-Z coloring
        // (9*NumberOfOrientations,  3*NumberOfOrientations, and 3*NumberOfOrientations in size, respectively)
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

        // Analysis routines
        // Part 1: ExaConstit-specific RVE(s)
        writeExaConstitRVE(NumberOfRVEs, OutputFileName, nx, ny, nz, deltax, GrainID, XLow_RVE, XHigh_RVE, YLow_RVE,
                           YHigh_RVE, ZLow_RVE, ZHigh_RVE); // if > 0 RVEs, print file(s) "n_ExaConstit.csv"

        // Part 2: Cross-sections for area statistics, pole figures, inverse pole figure coloring maps of microstructure
        printCrossSectionData(NumberOfCrossSections, OutputFileName, CrossSectionPlane, CrossSectionLocation, nx, ny,
                              nz, NumberOfOrientations, GrainID, PrintSectionPF, PrintSectionIPF, BimodalAnalysis,
                              deltax, GrainUnitVector_Host, GrainEulerAngles_Host, GrainRGBValues_Host, CSLabels);

        // Part 3: Representative volume grain statistics
        // List of grain ID values in the representative region
        int RepresentativeRegionSize_Cells = (XMax - XMin + 1) * (YMax - YMin + 1) * (ZMax - ZMin + 1);
        double RepresentativeRegionSize_Microns = RepresentativeRegionSize_Cells * convertToMicrons(deltax, "volume");

        std::vector<int> GrainIDVector = getRepresentativeRegionGrainIDs(GrainID, XMin, XMax, YMin, YMax, ZMin, ZMax,
                                                                         RepresentativeRegionSize_Cells);
        // Get the number of grains in the representative region and a list of the unique grain IDs
        int NumberOfGrains;
        std::vector<int> UniqueGrainIDVector = getUniqueGrains(GrainIDVector, NumberOfGrains);
        // Get the size (in units of length, area, or volume) associated with each unique grain ID value
        std::vector<float> GrainSizeVector =
            getGrainSizes(GrainIDVector, UniqueGrainIDVector, NumberOfGrains, deltax, "volume");

        // Output file stream for quantities of interest
        std::ofstream QoIs;
        // Header data for QoIs file
        printAnalysisHeader(OutputFileName, QoIs, XMin, XMax, YMin, YMax, ZMin, ZMax, XYZBounds);
        // Fraction of region consisting of nucleated grains, unmelted material
        printGrainTypeFractions(QoIs, XMin, XMax, YMin, YMax, ZMin, ZMax, GrainID, LayerID,
                                RepresentativeRegionSize_Cells);

        // Calculate misorientation data
        std::vector<float> GrainMisorientationXVector =
            getGrainMisorientation("X", GrainUnitVector, UniqueGrainIDVector, NumberOfOrientations, NumberOfGrains);
        std::vector<float> GrainMisorientationYVector =
            getGrainMisorientation("Y", GrainUnitVector, UniqueGrainIDVector, NumberOfOrientations, NumberOfGrains);
        std::vector<float> GrainMisorientationZVector =
            getGrainMisorientation("Z", GrainUnitVector, UniqueGrainIDVector, NumberOfOrientations, NumberOfGrains);
        // Print mean misorientation data to stats file
        printMeanMisorientations(QoIs, NumberOfGrains, GrainMisorientationXVector, GrainMisorientationYVector,
                                 GrainMisorientationZVector, GrainSizeVector, RepresentativeRegionSize_Microns);
        // Also print older misorientation stats data (to be removed in a future release)
        printMisorientationDataOld(XMin, XMax, YMin, YMax, ZMin, ZMax, LayerID, GrainUnitVector_Host, GrainID,
                                   NumberOfOrientations);
        // Print mean size data
        printMeanSize(QoIs, NumberOfGrains, RepresentativeRegionSize_Microns, "volume");

        // Need grain extents in x and y if analysis options 2 or 6 are toggled
        // Need grain extents in z if analysis options 2 or 5 are toggled
        // Need aspect ratios if analysis option 2 is toggled (after obtaining grain extents)
        // Extents are calculated in microns
        std::vector<float> GrainExtentX(NumberOfGrains);
        std::vector<float> GrainExtentY(NumberOfGrains);
        std::vector<float> GrainExtentZ(NumberOfGrains);
        std::vector<float> BuildTransAspectRatio(NumberOfGrains);
        if ((AnalysisTypes[2]) || (AnalysisTypes[6])) {
            calcGrainExtent(GrainExtentX, GrainID, UniqueGrainIDVector, GrainSizeVector, NumberOfGrains, XMin, XMax,
                            YMin, YMax, ZMin, ZMax, "X", deltax, "volume");
            calcGrainExtent(GrainExtentY, GrainID, UniqueGrainIDVector, GrainSizeVector, NumberOfGrains, XMin, XMax,
                            YMin, YMax, ZMin, ZMax, "Y", deltax, "volume");
            printMeanExtent(QoIs, GrainExtentX, "X", NumberOfGrains);
            printMeanExtent(QoIs, GrainExtentY, "Y", NumberOfGrains);
        }
        if ((AnalysisTypes[2]) || (AnalysisTypes[5])) {
            calcGrainExtent(GrainExtentZ, GrainID, UniqueGrainIDVector, GrainSizeVector, NumberOfGrains, XMin, XMax,
                            YMin, YMax, ZMin, ZMax, "Z", deltax, "volume");
            printMeanExtent(QoIs, GrainExtentZ, "Z", NumberOfGrains);
        }
        if (AnalysisTypes[2]) {
            calcBuildTransAspectRatio(BuildTransAspectRatio, GrainExtentX, GrainExtentY, GrainExtentZ, NumberOfGrains);
            printMeanBuildTransAspectRatio(QoIs, GrainExtentX, GrainExtentY, GrainExtentZ, GrainSizeVector,
                                           RepresentativeRegionSize_Microns, NumberOfGrains);
        }

        // If analysis option 6 is toggled, also print grain width data (x and y direction, and old grain width options)
        if (AnalysisTypes[6]) {
            printSizeOld(OutputFileName, NumberOfGrains, GrainExtentX, GrainExtentY, XMin, XMax, YMin, YMax, ZMax,
                         deltax, GrainID);
        }

        // Write grain area data as a function of Z location in the representative volume if option 3 (weighted mean
        // grain area) or option 4 (unweighted mean grain area) are toggled, writing to files
        // "[OutputFileName]_WeightedGrainAreas.csv" and "[OutputFileName]_GrainAreas.csv", respectively
        if ((AnalysisTypes[3]) || (AnalysisTypes[4]))
            writeAreaSeries(AnalysisTypes[3], AnalysisTypes[4], OutputFileName, deltax, XMin, XMax, YMin, YMax, ZMin,
                            ZMax, GrainID, XYZBounds[4]);
        QoIs.close();

        // Determine IPF-Z color of each grain relative to each direction: 0 (red), 1 (green), 2 (blue)
        bool PrintIPFRGB = false; // not yet used outside of grain_analysis_amb
        std::vector<float> GrainRed =
            getIPFZColor(0, UniqueGrainIDVector, NumberOfOrientations, GrainRGBValues_Host, NumberOfGrains);
        std::vector<float> GrainGreen =
            getIPFZColor(1, UniqueGrainIDVector, NumberOfOrientations, GrainRGBValues_Host, NumberOfGrains);
        std::vector<float> GrainBlue =
            getIPFZColor(1, UniqueGrainIDVector, NumberOfOrientations, GrainRGBValues_Host, NumberOfGrains);

        // Write per-grain stats for the analysis types specified to the file "[OutputFileName]_grains.csv"
        writePerGrainStats(OutputFileName, "volume", UniqueGrainIDVector, GrainMisorientationXVector,
                           GrainMisorientationYVector, GrainMisorientationZVector, GrainSizeVector, GrainExtentX,
                           GrainExtentY, GrainExtentZ, BuildTransAspectRatio, AnalysisTypes, NumberOfGrains,
                           PrintIPFRGB, GrainRed, GrainGreen, GrainBlue);
        // If analysis option 7 is toggled, print orientation data to file
        // "[OutputFileName]_RegionLabel_PoleFigureData.txt"
        if (AnalysisTypes[7]) {
            ViewI_H GOHistogram =
                getOrientationHistogram(NumberOfOrientations, GrainID, LayerID, XMin, XMax, YMin, YMax, ZMin, ZMax);
            std::string RegionLabel = "VolumeX" + std::to_string(XMin) + "-" + std::to_string(XMax) + "Y" +
                                      std::to_string(YMin) + "-" + std::to_string(YMax) + "Z" + std::to_string(ZMin) +
                                      "-" + std::to_string(ZMax);
            writePoleFigure(OutputFileName, RegionLabel, NumberOfOrientations, GrainEulerAngles, GOHistogram);
        }
    }
    // Finalize kokkos and end program
    Kokkos::finalize();

    return 0;
}
