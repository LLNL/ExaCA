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
    if (argc == 0) {
        throw std::runtime_error("Error: Neither ExaCA microstructure nor analysis file were given on the command line");
    }
    else if (argc == 1) {
        throw std::runtime_error("Error: One of the two required inptus (ExaCA microstructure file, output analysis file) was not given on the command line");
    }
    else {
        InputFile = argv[1];
        OutputAnalysisFile = argv[2];
        std::size_t StrLength = InputFile.length();
        BaseFileName = InputFile.substr(0,StrLength-3);
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
    // Determine the number of, and bounds of, specified RVEs and ANG (pyEBSD) cross-sections
    bool AnalysisTypes[8] = {0};
    int NumberOfRVEs = 0;
    int NumberOfANGCrossSections = 0;
    std::vector<int> CenterX_RVE, CenterY_RVE, CenterZ_RVE, Size_RVE, ANGCrossSectionPlane, ANGCrossSectionLocation; // Contains data on each individual RVE/ANGCrossSection
    std::vector<bool> ExaConstitPrint, pyEBSDPrint; // whether or not each RVE is printed to files for ExaConstit or pyEBSD
    std::string RotationFilename, EulerFilename; // Names of orientation files specified by the analysis file
    int XMin, XMax, YMin, YMax; // bounds of analysis cross-section
    int NumberOfOrientations, NumberOfMeltedCells;
    ParseAnalysisFile(OutputAnalysisFile, RotationFilename, EulerFilename, NumberOfOrientations, AnalysisTypes, CenterX_RVE, CenterY_RVE, CenterZ_RVE, Size_RVE, ExaConstitPrint, pyEBSDPrint, NumberOfRVEs, ANGCrossSectionPlane, ANGCrossSectionLocation, NumberOfANGCrossSections, XMin, XMax, YMin, YMax, nx, ny, nz, LayerID, Melted, NumberOfLayers);
    // Allocate memory blocks for grain unit vectors (NumberOfOrientations by 3 by 3) and grain euler angles (NumberOfOrientations by 3)
    float*** GrainUnitVector = new float**[NumberOfOrientations];
    float** GrainEulerAngles = new float*[NumberOfOrientations];
    for (int n = 0; n < NumberOfOrientations; n++) {
        // Allocate memory blocks for rows of the 1D array
        GrainEulerAngles[n] = new float[3];
        // Allocate memory blocks for rows of the 2D array
        GrainUnitVector[n] = new float*[3];
        for (int i = 0; i < 3; i++) {
            // Allocate memory blocks for columns of the 2D array
            GrainUnitVector[n][i] = new float[3];
        }
    }
    ParseGrainOrientationFiles(RotationFilename, EulerFilename, NumberOfOrientations, GrainUnitVector, GrainEulerAngles);
    
    // Analysis routines
    PrintRVEData(NumberOfRVEs, BaseFileName, nx, ny, nz, deltax, GrainID, CenterX_RVE, CenterY_RVE, CenterZ_RVE, Size_RVE, ExaConstitPrint, pyEBSDPrint); // if > 0 RVEs, print file(s) "n_ExaConstit.csv"
    // PrintANGCrossSections();
    PrintMisorientationData(AnalysisTypes, BaseFileName, XMin, XMax, YMin, YMax, nz, Melted, GrainUnitVector, GrainID, NumberOfOrientations, NumberOfMeltedCells); // if analysis options 0 or 1 are toggled, print data to std::out and/or file "_MisorientationFrequency.txt"
    PrintVolumeData(AnalysisTypes, XMin, XMax, YMin, YMax, nz, NumberOfMeltedCells, Melted, GrainID); // Print volume fraction of nucleated grains, options for mean grain volume and volume distribution
    PrintGrainAreaData(AnalysisTypes, BaseFileName, deltax, XMin, XMax, YMin, YMax, nz, GrainID);
    
    return 0;
}
