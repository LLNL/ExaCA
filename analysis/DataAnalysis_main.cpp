// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "header.h"

#include <stdexcept>
#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

// Required inputs on the command line: first an ExaCA microstructure (.vtk format), then an ExaCA output analysis file
// There should also exist a log file with the same name as the ExaCA microstructure (with .log instead of .vtk), in the same directory as the ExaCA microstructure file
int main(int argc, char *argv[]) {
    
    // Read command line input to obtain ExaCA microstructure file name - throw error if no file is given
    std::string InputFile, OutputAnalysisFile, BaseFileName;
    double deltax;
    if (argc < 2) {
        throw std::runtime_error("Error: Neither ExaCA microstructure nor analysis file were given on the command line");
    }
    else if (argc != 3) {
        throw std::runtime_error("Error: Either one of the two required inputs (ExaCA microstructure file, output analysis file) was not given on the command line, or an excessive number of inputs was given");
    }
    else {
        InputFile = argv[1];
        OutputAnalysisFile = argv[2];
        std::size_t StrLength = InputFile.length();
        BaseFileName = InputFile.substr(0,StrLength-4);
    }
    
    // Given the input file name, parse the paraview file for the cell size, x, y, and z dimensions, number of layers
    int nx, ny, nz, NumberOfLayers;
    ParseLogFile(BaseFileName, nx, ny, nz, deltax, NumberOfLayers);
    
    // Allocate memory blocks for GrainID, LayerID, and Melted data
    // Alocate them as 2D arrays of size nx * ny x i.e., no of 2D Arrays
    int*** GrainID = new int**[nz];
    int*** LayerID = new int**[nz];
    int*** Melted = new int**[nz];
    for (int k = 0; k < nz; k++) {
        // Allocate memory blocks for rows of each 2D array
        GrainID[k] = new int*[nx];
        LayerID[k] = new int*[nx];
        Melted[k] = new int*[nx];
        for (int i = 0; i < nx; i++) {
            // Allocate memory blocks for columns of each 2D array
            GrainID[k][i] = new int[ny];
            LayerID[k][i] = new int[ny];
            Melted[k][i] = new int[ny];
        }
    }
    
    // Fill arrays with data from paraview file
    InitializeData(InputFile, nx, ny, nz, GrainID, LayerID, Melted);
    
    // Read analysis file ("ExaCA/examples/Outputs.txt") to determine which analysis should be done
    // There are three parts to this analysis:
    // Part 1: Determine the number of, and bounds of, specified RVEs for ExaConstits
    int NumberOfRVEs = 0;
    std::vector<int> XLow_RVE, XHigh_RVE, YLow_RVE, YHigh_RVE, ZLow_RVE, ZHigh_RVE; // Contains data on each individual RVE bounds
    
    // Part 2: Determine the number of, and planes for, cross-sections to be analyzed by pyEBSD/MTEX
    int NumberOfCrossSections = 0;
    std::vector<int> CrossSectionPlane, CrossSectionLocation; // Contains data on each individual ANG CrossSection
    
    // Part 3: Analysis options on a representative region in x, y, and z (can be different than the x, y, and z for ExaConstit)
    bool AnalysisTypes[8] = {0}; // which analysis modes other than the defaults should be considered?
    int XMin, XMax, YMin, YMax, ZMin, ZMax; // bounds of analysis region
    
    // Analysis also requires reading orientation files of specified names
    int NumberOfOrientations;
    std::string RotationFilename, EulerFilename; // Names of orientation files specified by the analysis file
    ParseAnalysisFile(OutputAnalysisFile, RotationFilename, EulerFilename, NumberOfOrientations, AnalysisTypes, XLow_RVE, XHigh_RVE, YLow_RVE, YHigh_RVE, ZLow_RVE, ZHigh_RVE, NumberOfRVEs, CrossSectionPlane, CrossSectionLocation, NumberOfCrossSections, XMin, XMax, YMin, YMax, ZMin, ZMax, nx, ny, nz, LayerID, Melted, NumberOfLayers);
    
    // Allocate memory blocks for grain unit vectors (NumberOfOrientations by 3 by 3)
    // grain euler angles (NumberOfOrientations by 3) will be used in the future for pyEBSD analysis/data printing
    float*** GrainUnitVector = new float**[NumberOfOrientations];
    float** GrainEulerAngles = new float*[NumberOfOrientations];
    for (int n = 0; n < NumberOfOrientations; n++) {
        // Allocate memory blocks for rows of the 1D array (GrainEulerAngles not yet used/allocated - part of future work with pyEBSD)
        // GrainEulerAngles[n] = new float[3];
        // Allocate memory blocks for rows of the 2D array
        GrainUnitVector[n] = new float*[3];
        for (int i = 0; i < 3; i++) {
            // Allocate memory blocks for columns of the 2D array
            GrainUnitVector[n][i] = new float[3];
        }
    }
    ParseGrainOrientationFiles(RotationFilename, EulerFilename, NumberOfOrientations, GrainUnitVector, GrainEulerAngles);
    
    // Analysis routines
    // Part 1: ExaConstit-specific RVE(s)
    PrintExaConstitRVEData(NumberOfRVEs, BaseFileName, nx, ny, nz, deltax, GrainID, XLow_RVE, XHigh_RVE, YLow_RVE, YHigh_RVE, ZLow_RVE, ZHigh_RVE); // if > 0 RVEs, print file(s) "n_ExaConstit.csv"
    
    // Part 2: Cross-sections for inverse pole figures
    PrintInversePoleFigureCrossSections(NumberOfCrossSections, BaseFileName, CrossSectionPlane, CrossSectionLocation, nx, ny, nz, NumberOfOrientations, GrainID, GrainEulerAngles);

    // Part 3: Representative volume grain statistics
    // Print data to std::out and if analysis option 0 is toggled, to the file "[BaseFileName]_MisorientationFrequency.csv"
    PrintMisorientationData(AnalysisTypes, BaseFileName, XMin, XMax, YMin, YMax, ZMin, ZMax, Melted, GrainUnitVector, GrainID, NumberOfOrientations);
    // Print data to std::out and if analysis options 1, 2, or 5 are toggled, print data to files "[BaseFileName]_VolumeFrequency.csv", "[BaseFileName]_AspectRatioFrequency.csv", and [BaseFileName]_GrainHeightDistribution.csv", respectively
    PrintSizeData(AnalysisTypes, BaseFileName, XMin, XMax, YMin, YMax, ZMin, ZMax, nx, ny, nz, Melted, GrainID, deltax);
    // Print data to std::out and if analysis options 3, 4, or 6 are toggled, print data to files "[BaseFileName]_GrainAreas.csv", "[BaseFileName]_WeightedGrainAreas.csv", and "[BaseFileName]_GrainWidthDistribution.csv", respectively
    PrintGrainAreaData(AnalysisTypes, BaseFileName, deltax, XMin, XMax, YMin, YMax, ZMin, ZMax, GrainID);
    // If analysis option 7 is toggled, print orientation data to files "[BaseFileName]_pyEBSDOrientations.csv" and "[BaseFileName]_MTEXOrientations.csv"
    PrintPoleFigureData(AnalysisTypes, BaseFileName, NumberOfOrientations, GrainEulerAngles, XMin, XMax, YMin, YMax, ZMin, ZMax, GrainID, Melted);
    
    return 0;
}
