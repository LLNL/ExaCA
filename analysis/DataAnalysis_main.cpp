// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "header.hpp"

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
    std::string BaseFileName, AnalysisFile, LogFile, MicrostructureFile, RotationFilename;
    double deltax;
    if (argc < 2) {
        throw std::runtime_error("Error: Analysis file name must be given on the command line");
    }
    else {
        AnalysisFile = argv[1];
        // Get path to/name of all files of interest by reading the analysis file
        ParseFilenames(AnalysisFile, LogFile, MicrostructureFile, RotationFilename, BaseFileName);
    }
    std::cout << "Performing analysis of " << MicrostructureFile << " , using the log file " << LogFile
              << " and the options specified in " << AnalysisFile << std::endl;

    // Initialize Kokkos
    Kokkos::initialize();
    {

        // Given the input file name, parse the paraview file for the cell size, x, y, and z dimensions, number of
        // layers
        int nx, ny, nz, NumberOfLayers;
        ParseLogFile(LogFile, nx, ny, nz, deltax, NumberOfLayers);

        // Allocate memory blocks for GrainID, LayerID, and Melted data
        ViewI3D_H GrainID(Kokkos::ViewAllocateWithoutInitializing("GrainID"), nz, ny, ny);
        ViewI3D_H LayerID(Kokkos::ViewAllocateWithoutInitializing("LayerID"), nz, ny, ny);
        ViewI3D_H Melted(Kokkos::ViewAllocateWithoutInitializing("Melted"), nz, ny, ny);

        // Fill arrays with data from paraview file
        InitializeData(MicrostructureFile, nx, ny, nz, GrainID, LayerID, Melted);

        // Read analysis file ("ExaCA/examples/Outputs.txt") to determine which analysis should be done
        // There are three parts to this analysis:
        // Part 1: Determine the number of, and bounds of, specified RVEs for ExaConstits
        int NumberOfRVEs = 0;
        std::vector<int> XLow_RVE, XHigh_RVE, YLow_RVE, YHigh_RVE, ZLow_RVE,
            ZHigh_RVE; // Contains data on each individual RVE bounds

        // Part 2: Determine the number of, and planes for, cross-sections to be analyzed by pyEBSD/MTEX
        int NumberOfCrossSections = 0;
        std::vector<int> CrossSectionPlane, CrossSectionLocation; // Contains data on each individual ANG CrossSection

        // Part 3: Analysis options on a representative region in x, y, and z (can be different than the x, y, and z for
        // ExaConstit)
        bool AnalysisTypes[8] = {0};            // which analysis modes other than the defaults should be considered?
        int XMin, XMax, YMin, YMax, ZMin, ZMax; // bounds of analysis region

        // Analysis also requires reading orientation files of specified names
        int NumberOfOrientations;
        ParseAnalysisFile(AnalysisFile, RotationFilename, NumberOfOrientations, AnalysisTypes, XLow_RVE, XHigh_RVE,
                          YLow_RVE, YHigh_RVE, ZLow_RVE, ZHigh_RVE, NumberOfRVEs, CrossSectionPlane,
                          CrossSectionLocation, NumberOfCrossSections, XMin, XMax, YMin, YMax, ZMin, ZMax, nx, ny, nz,
                          LayerID, Melted, NumberOfLayers);

        // Allocate memory for grain unit vectors (NumberOfOrientations by 3 by 3)
        ViewF3D_H GrainUnitVector(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), NumberOfOrientations, 3,
                                  3);

        ParseGrainOrientationFiles(RotationFilename, NumberOfOrientations, GrainUnitVector);

        // Analysis routines
        // Part 1: ExaConstit-specific RVE(s)
        PrintExaConstitRVEData(NumberOfRVEs, BaseFileName, nx, ny, nz, deltax, GrainID, XLow_RVE, XHigh_RVE, YLow_RVE,
                               YHigh_RVE, ZLow_RVE, ZHigh_RVE); // if > 0 RVEs, print file(s) "n_ExaConstit.csv"

        // Part 2: Cross-sections for inverse pole figures
        PrintInversePoleFigureCrossSections(NumberOfCrossSections, BaseFileName, CrossSectionPlane,
                                            CrossSectionLocation, nx, ny, nz, NumberOfOrientations, GrainID);

        // Part 3: Representative volume grain statistics
        // Print data to std::out and if analysis option 0 is toggled, to the file
        // "[BaseFileName]_MisorientationFrequency.csv"
        PrintMisorientationData(AnalysisTypes, BaseFileName, XMin, XMax, YMin, YMax, ZMin, ZMax, Melted,
                                GrainUnitVector, GrainID, NumberOfOrientations);
        // Print data to std::out and if analysis options 1, 2, or 5 are toggled, print data to files
        // "[BaseFileName]_VolumeFrequency.csv", "[BaseFileName]_AspectRatioFrequency.csv", and
        // [BaseFileName]_GrainHeightDistribution.csv", respectively
        PrintSizeData(AnalysisTypes, BaseFileName, XMin, XMax, YMin, YMax, ZMin, ZMax, nx, ny, nz, Melted, GrainID,
                      deltax);
        // Print data to std::out and if analysis options 3, 4, or 6 are toggled, print data to files
        // "[BaseFileName]_GrainAreas.csv", "[BaseFileName]_WeightedGrainAreas.csv", and
        // "[BaseFileName]_GrainWidthDistribution.csv", respectively
        PrintGrainAreaData(AnalysisTypes, BaseFileName, deltax, XMin, XMax, YMin, YMax, ZMin, ZMax, GrainID);
        // If analysis option 7 is toggled, print orientation data to files "[BaseFileName]_pyEBSDOrientations.csv" and
        // "[BaseFileName]_MTEXOrientations.csv"
        PrintPoleFigureData(AnalysisTypes, BaseFileName, NumberOfOrientations, XMin, XMax, YMin, YMax, ZMin, ZMax,
                            GrainID, Melted);
    }
    // Finalize kokkos and end program
    Kokkos::finalize();

    return 0;
}
