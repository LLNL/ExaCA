// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "runGA.hpp"
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
    bool NewOrientationFormatYN;
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
                       NewOrientationFormatYN, RGBFilename);
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
        ParseLogFile(LogFile, nx, ny, nz, deltax, NumberOfLayers, false, XYZBounds);

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
                          LayerID, NumberOfLayers, PrintSectionPF, PrintSectionIPF, BimodalAnalysis,
                          NewOrientationFormatYN, CSLabels);

        // Allocate memory for grain unit vectors, grain euler angles, RGB colors for IPF-Z coloring
        // (9*NumberOfOrientations,  3*NumberOfOrientations, and 3*NumberOfOrientations in size, respectively)
        ViewF GrainUnitVector(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), 9 * NumberOfOrientations);
        ViewF GrainEulerAngles(Kokkos::ViewAllocateWithoutInitializing("GrainEulerAngles"), 3 * NumberOfOrientations);
        ViewF GrainRGBValues(Kokkos::ViewAllocateWithoutInitializing("GrainRGBValues"), 3 * NumberOfOrientations);

        // Initialize, then copy back to host. GrainEulerAngles only used with new input file format
        OrientationInit(0, NumberOfOrientations, GrainUnitVector, RotationFilename, 9);
        ViewF_H GrainUnitVector_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainUnitVector);
        if (NewOrientationFormatYN) {
            OrientationInit(0, NumberOfOrientations, GrainEulerAngles, EulerAnglesFilename, 3);
        }
        ViewF_H GrainEulerAngles_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainEulerAngles);
        OrientationInit(0, NumberOfOrientations, GrainRGBValues, RGBFilename, 3);
        ViewF_H GrainRGBValues_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainRGBValues);

        // Analysis routines
        // Part 1: ExaConstit-specific RVE(s)
        PrintExaConstitRVEData(NumberOfRVEs, OutputFileName, nx, ny, nz, deltax, GrainID, XLow_RVE, XHigh_RVE, YLow_RVE,
                               YHigh_RVE, ZLow_RVE, ZHigh_RVE); // if > 0 RVEs, print file(s) "n_ExaConstit.csv"

        // Part 2: Cross-sections for area statistics, pole figures, inverse pole figure coloring maps of microstructure
        PrintCrossSectionData(NumberOfCrossSections, OutputFileName, CrossSectionPlane, CrossSectionLocation, nx, ny,
                              nz, NumberOfOrientations, GrainID, PrintSectionPF, PrintSectionIPF, BimodalAnalysis,
                              NewOrientationFormatYN, deltax, GrainUnitVector_Host, GrainEulerAngles_Host,
                              GrainRGBValues_Host, CSLabels);

        // Part 3: Representative volume grain statistics
        // Print data to std::out and if analysis option 0 is toggled, to the file
        // "[OutputFileName]_MisorientationFrequency.csv"
        PrintMisorientationData(AnalysisTypes, OutputFileName, XMin, XMax, YMin, YMax, ZMin, ZMax, LayerID,
                                GrainUnitVector_Host, GrainID, NumberOfOrientations);
        // Print data to std::out and if analysis options 1, 2, or 5 are toggled, print data to files
        // "[OutputFileName]_VolumeFrequency.csv", "[OutputFileName]_AspectRatioFrequency.csv", and
        // [OutputFileName]_GrainHeightDistribution.csv", respectively
        PrintSizeData(AnalysisTypes, OutputFileName, XMin, XMax, YMin, YMax, ZMin, ZMax, nx, ny, nz, LayerID, GrainID,
                      deltax);
        // Print data to std::out and if analysis options 3, 4, or 6 are toggled, print data to files
        // "[OutputFileName]_GrainAreas.csv", "[OutputFileName]_WeightedGrainAreas.csv", and
        // "[OutputFileName]_GrainWidthDistribution.csv", respectively
        PrintGrainAreaData(AnalysisTypes, OutputFileName, deltax, XMin, XMax, YMin, YMax, ZMin, ZMax, GrainID);
        // If analysis option 7 is toggled, print orientation data to files "[OutputFileName]_pyEBSDOrientations.csv"
        // and "[OutputFileName]_MTEXOrientations.csv"
        PrintPoleFigureData(AnalysisTypes, OutputFileName, NumberOfOrientations, XMin, XMax, YMin, YMax, ZMin, ZMax,
                            GrainID, LayerID, NewOrientationFormatYN, GrainEulerAngles_Host);
    }
    // Finalize kokkos and end program
    Kokkos::finalize();

    return 0;
}
