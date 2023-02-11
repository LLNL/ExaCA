// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAinitialize.hpp"
#include "CAparsefiles.hpp"
#include "CAtypes.hpp"
#include "GAprint.hpp"
#include "GAutils.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// TODO: This file will be removed in a future release, and these hardcoded analysis options will be run using
// grain_analysis and a default analysis input file The name of the ExaCA output analysis file (.txt, located in
// ExaCA/analysis), the ExaCA microstructure file (.vtk), and the ExaCA log file (.log) is given on the command line -
// these files should all have the same name, other than the file extension
int main(int argc, char *argv[]) {

    // Read command line input to obtain name of analysis file
    std::string BaseFileName, LogFile, MicrostructureFile, RotationFilename, EulerAnglesFilename, RGBFilename,
        OutputFileName;
    double deltax;
    if (argc < 2) {
        throw std::runtime_error("Error: Full path to/prefix for vtk and log data must be given");
    }
    else {
        BaseFileName = argv[1];
#ifdef ExaCA_ENABLE_JSON
        // If json is enabled, first check if log is in json format
        LogFile = BaseFileName + ".json";
        bool Json_found = checkFileExists(LogFile, 0, false);
        // If json is enabled and the .json file does not exist, try the .log file
        if (!(Json_found))
            LogFile = BaseFileName + ".log";
#else
        LogFile = BaseFileName + ".log";
#endif
        MicrostructureFile = BaseFileName + ".vtk";
    }
    std::cout << "Performing analysis of " << MicrostructureFile << " , using the log file " << LogFile << std::endl;

    // Initialize Kokkos
    Kokkos::initialize();
    {
        // Given the input file name, parse the paraview file for the cell size, x, y, and z dimensions, number of
        // layers, X, Y, Z lower and upper bounds
        int nx, ny, nz, NumberOfLayers;
        std::vector<double> XYZBounds(6);
        int LogFormat = checkLogFormat(LogFile);
        if (LogFormat == 0) {
            ParseLogFile_OldNoColon(LogFile, nx, ny, nz, deltax, NumberOfLayers, true, XYZBounds);
            // Use default file names as they were not in the log file, and there is no analysis input file
            CheckInputFiles(LogFile, MicrostructureFile, RotationFilename, EulerAnglesFilename, RGBFilename);
        }
        else if (LogFormat == 1) {
            ParseLogFile_Old(LogFile, nx, ny, nz, deltax, NumberOfLayers, XYZBounds, RotationFilename,
                             EulerAnglesFilename, RGBFilename, false);
        }
        else {
#ifdef ExaCA_ENABLE_JSON
            ParseLogFile(LogFile, nx, ny, nz, deltax, NumberOfLayers, XYZBounds, RotationFilename, EulerAnglesFilename,
                         RGBFilename, false);
#endif
        }

        // Allocate memory blocks for GrainID and LayerID data
        ViewI3D_H GrainID(Kokkos::ViewAllocateWithoutInitializing("GrainID"), nz, nx, ny);
        ViewI3D_H LayerID(Kokkos::ViewAllocateWithoutInitializing("LayerID"), nz, nx, ny);

        // Fill arrays with data from paraview file
        InitializeData(MicrostructureFile, nx, ny, nz, GrainID, LayerID);

        // Allocate memory for grain unit vectors, grain euler angles, RGB colors for IPF-Z coloring
        // (9*NumberOfOrientations,  3*NumberOfOrientations, and 3*NumberOfOrientations in size, respectively)
        // No initialize size yet, will be resized in OrientationInit when NumberOfOrientations is known
        int NumberOfOrientations;
        ViewF GrainUnitVector(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), 0);
        ViewF GrainEulerAngles(Kokkos::ViewAllocateWithoutInitializing("GrainEulerAngles"), 0);
        ViewF GrainRGBValues(Kokkos::ViewAllocateWithoutInitializing("GrainRGBValues"), 0);

        // Initialize, then copy back to host
        OrientationInit(0, NumberOfOrientations, GrainUnitVector, RotationFilename, 9);
        ViewF_H GrainUnitVector_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainUnitVector);
        OrientationInit(0, NumberOfOrientations, GrainEulerAngles, EulerAnglesFilename, 3);
        ViewF_H GrainEulerAngles_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainEulerAngles);
        OrientationInit(0, NumberOfOrientations, GrainRGBValues, RGBFilename, 3);
        ViewF_H GrainRGBValues_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainRGBValues);

        // Print information about 2 specific cross-sections: Z = 1 mm above the base plate, and X = midway through the
        // volume of the microstructure
        // TODO: create class "AnalysisRegion" to store this information, run analysis on each region as part of a for
        // loop Class will contain: regionUnits (cells or microns, converting to whichever wasn't given)
        // regionBoundsXLower (stored in cells and in microns)
        // regionBoundsXUpper (stored in cells and in microns)
        // regionBoundsYLower (stored in cells and in microns)
        // regionBoundsYUpper (stored in cells and in microns)
        // regionBoundsZLower (stored in cells and in microns)
        // regionBoundsZUpper (stored in cells and in microns)
        // regionSizeX, regionSizeY, regionSizeZ (stored in cells)
        // regionSize (number of cells in the region)
        // regionType ("length", "area", or "volume", depending on the regionBounds supplied)
        // regionOrientation ("X", "Y", or "Z" for length, "XY", "YZ", or "XY" for area, not used for volume)
        // regionName (string label for output files)
        // true/false analysis options:
        // Name                               // Region      // Print location
        // ================================== // =========== // ====================================
        // Grain Type Fractions               // All         // std::out/QoIs file
        // Misorientations                    // All         // std::out/QoIs file and perGrain file
        // Size                               // All         // std::out/QoIs file and perGrain file
        // BuildTrans Aspect Ratio            // Volume      // std::out/QoIs file and perGrain file
        // Extent in X, Y, or Z               // All         // std::out/QoIs file and perGrain file
        // Area/weighted area vs build height // Volume      // separate file(s)
        // IPF-Z coloring in X, Y, or Z       // Area        // perGrain file
        // Pole figure data for MTEX          // Area/Volume // separate file
        // Inverse pole figure data for MTEX  // Area        // separate file
        // Per grain stats                    // All         // perGrain file (true if at least one perGrain option is
        // true, false otherwise)
        // ExaConstit RVE                     // Volume      // separate file

        int NumberOfRegions = 2;
        // First region is an XY cross-section located at ZLocCrossSection
        // Second region is a YZ cross-section located at XLocCrossSection
        int ZLocCrossSection = std::round((0.001 - XYZBounds[2]) / deltax);
        int XLocCrossSection = std::round(nx / 2);
        std::vector<int> regionBoundsXLower_cells = {0, XLocCrossSection};
        std::vector<int> regionBoundsXUpper_cells = {nx - 1, XLocCrossSection};
        std::vector<int> regionSizeX = {nx - 1, 1};
        std::vector<int> regionBoundsYLower_cells = {0, 0};
        std::vector<int> regionBoundsYUpper_cells = {ny - 1, ny - 1};
        std::vector<int> regionSizeY = {ny - 1, ny - 1};
        std::vector<int> regionBoundsZLower_cells = {ZLocCrossSection, 0};
        std::vector<int> regionBoundsZUpper_cells = {ZLocCrossSection, nz - 1};
        std::vector<int> regionSizeZ = {1, nz - 1};
        std::vector<int> regionSize = {nx * ny, ny * nz};
        std::vector<double> regionBoundsXLower_microns(2), regionBoundsXUpper_microns(2), regionBoundsYLower_microns(2),
            regionBoundsYUpper_microns(2), regionBoundsZLower_microns(2), regionBoundsZUpper_microns(2);
        for (int n = 0; n < 2; n++) {
            regionBoundsXLower_microns[n] = regionBoundsXLower_cells[n] * deltax;
            regionBoundsXUpper_microns[n] = regionBoundsXLower_microns[n] + regionSizeX[n] * deltax;
            regionBoundsYLower_microns[n] = regionBoundsYLower_cells[n] * deltax;
            regionBoundsYUpper_microns[n] = regionBoundsYLower_microns[n] + regionSizeX[n] * deltax;
            regionBoundsZLower_microns[n] = regionBoundsZLower_cells[n] * deltax;
            regionBoundsZUpper_microns[n] = regionBoundsZLower_microns[n] + regionSizeX[n] * deltax;
        }
        std::vector<std::string> regionType = {"plane", "plane"};
        std::vector<std::string> regionOrientation = {"XY", "YZ"};
        std::vector<std::string> regionName = {"XY", "YZ"};
        std::vector<bool> printGrainTypeFractionsYN = {true, true};
        std::vector<bool> printSizeYN = {true, false};
        std::vector<bool> printIPFZColorYN = {true, false};
        std::vector<bool> printPerGrainStatsYN = {true, false};
        std::vector<bool> printPoleFigureYN = {true, true};
        std::vector<bool> printInversePoleFigureYN = {true, true};

        std::ofstream QoIs;
        std::string QoIs_fname = BaseFileName + "_QoIs.txt";
        QoIs.open(QoIs_fname);

        for (int n = 0; n < NumberOfRegions; n++) {
            printAnalysisHeader(QoIs, regionBoundsXLower_cells[n], regionBoundsXUpper_cells[n],
                                regionBoundsYLower_cells[n], regionBoundsYUpper_cells[n], regionBoundsZLower_cells[n],
                                regionBoundsZUpper_cells[n], regionBoundsXLower_microns[n],
                                regionBoundsXUpper_microns[n], regionBoundsYLower_microns[n],
                                regionBoundsYUpper_microns[n], regionBoundsZLower_microns[n],
                                regionBoundsZUpper_microns[n], regionName[n], regionOrientation[n]);
            if (printGrainTypeFractionsYN[n])
                printGrainTypeFractions(QoIs, regionBoundsXLower_cells[n], regionBoundsXUpper_cells[n],
                                        regionBoundsYLower_cells[n], regionBoundsYUpper_cells[n],
                                        regionBoundsZLower_cells[n], regionBoundsZUpper_cells[n], GrainID, LayerID,
                                        regionSize[n]);

            // TODO: Each of these vectors should also belong to the object of the AnalysisRegion class
            // List of grain IDs in this region
            std::vector<int> GrainIDVector = getRepresentativeRegionGrainIDs(
                GrainID, regionBoundsXLower_cells[n], regionBoundsXUpper_cells[n], regionBoundsYLower_cells[n],
                regionBoundsYUpper_cells[n], regionBoundsZLower_cells[n], regionBoundsZUpper_cells[n], regionSize[n]);
            // Get the number of grains in the representative region and a list of the unique grain IDs
            int NumberOfGrains;
            std::vector<int> UniqueGrainIDVector = getUniqueGrains(GrainIDVector, NumberOfGrains);
            // Get the size (in units of length, area, or volume) associated with each unique grain ID value
            std::vector<float> GrainSizeVector =
                getGrainSizes(GrainIDVector, UniqueGrainIDVector, NumberOfGrains, deltax, "area");
            // Get the IPF-Z (R,G,B) values corresponding to the IPF mapping of this grain from MTEX
            std::vector<float> GrainRed =
                getIPFZColor(0, UniqueGrainIDVector, NumberOfOrientations, GrainRGBValues_Host, NumberOfGrains);
            std::vector<float> GrainGreen =
                getIPFZColor(1, UniqueGrainIDVector, NumberOfOrientations, GrainRGBValues_Host, NumberOfGrains);
            std::vector<float> GrainBlue =
                getIPFZColor(2, UniqueGrainIDVector, NumberOfOrientations, GrainRGBValues_Host, NumberOfGrains);

            if (printPerGrainStatsYN[n]) {
                // TODO: Combine this and writePerGrainStats
                std::string stats_fname = BaseFileName + "_grains.csv";
                std::ofstream GrainStats;
                GrainStats.open(stats_fname);
                GrainStats << "GrainID,area,IPFZ_r,IPFZ_g,IPFZ_b" << std::endl;
                for (int i = 0; i < NumberOfGrains; i++) {
                    GrainStats << UniqueGrainIDVector[i] << "," << GrainSizeVector[i] << "," << GrainRed[i] << ","
                               << GrainGreen[i] << "," << GrainBlue[i] << std::endl;
                }
            }

            if (printPoleFigureYN[n]) {
                ViewI_H GOHistogram = getOrientationHistogram(NumberOfOrientations, GrainID, LayerID,
                                                              regionBoundsXLower_cells[n], regionBoundsXUpper_cells[n],
                                                              regionBoundsYLower_cells[n], regionBoundsYUpper_cells[n],
                                                              regionBoundsZLower_cells[n], regionBoundsZUpper_cells[n]);
                writePoleFigure(BaseFileName, regionName[n], NumberOfOrientations, GrainEulerAngles_Host, GOHistogram);
            }
            if (printInversePoleFigureYN[n]) {
                int Index1Low, Index1High, Index2Low, Index2High, CrossSectionOutOfPlaneLocation;
                if (regionOrientation[n] == "XY") {
                    Index1Low = regionBoundsXLower_cells[n];
                    Index1High = regionBoundsXUpper_cells[n];
                    Index2Low = regionBoundsYLower_cells[n];
                    Index2High = regionBoundsYUpper_cells[n];
                    CrossSectionOutOfPlaneLocation = regionBoundsZLower_cells[n];
                }
                else if (regionOrientation[n] == "XZ") {
                    Index1Low = regionBoundsXLower_cells[n];
                    Index1High = regionBoundsXUpper_cells[n];
                    Index2Low = regionBoundsZLower_cells[n];
                    Index2High = regionBoundsZUpper_cells[n];
                    CrossSectionOutOfPlaneLocation = regionBoundsYLower_cells[n];
                }
                else if (regionOrientation[n] == "YZ") {
                    Index1Low = regionBoundsYLower_cells[n];
                    Index1High = regionBoundsYUpper_cells[n];
                    Index2Low = regionBoundsZLower_cells[n];
                    Index2High = regionBoundsZUpper_cells[n];
                    CrossSectionOutOfPlaneLocation = regionBoundsXLower_cells[n];
                }
                else
                    throw std::runtime_error("Invalid region type");
                writeIPFColoredCrossSection(BaseFileName, regionName[n], regionOrientation[n], Index1Low, Index1High,
                                            Index2Low, Index2High, CrossSectionOutOfPlaneLocation, GrainID,
                                            GrainEulerAngles_Host, deltax, NumberOfOrientations);
            }
        }
    }
    // Finalize kokkos and end program
    Kokkos::finalize();

    return 0;
}
