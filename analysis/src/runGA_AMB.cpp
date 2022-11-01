// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "runGA_AMB.hpp"
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

// The name of the ExaCA output analysis file (.txt, located in ExaCA/analysis), the ExaCA microstructure file (.vtk),
// and the ExaCA log file (.log) is given on the command line - these files should all have the same name, other than
// the file extension
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
        CheckInputFiles_AMB(BaseFileName, LogFile, MicrostructureFile, RotationFilename, EulerAnglesFilename,
                            RGBFilename);
    }
    std::cout << "Performing analysis of " << MicrostructureFile << " , using the log file " << LogFile << std::endl;

    // Initialize Kokkos
    Kokkos::initialize();
    {

        // Given the input file name, parse the paraview file for the cell size, x, y, and z dimensions, number of
        // layers, X, Y, Z lower and upper bounds
        int nx, ny, nz, NumberOfLayers;
        std::vector<double> XYZBounds(6);
        ParseLogFile(LogFile, nx, ny, nz, deltax, NumberOfLayers, true, XYZBounds);

        // Allocate memory blocks for GrainID and LayerID data
        ViewI3D_H GrainID(Kokkos::ViewAllocateWithoutInitializing("GrainID"), nz, nx, ny);
        ViewI3D_H LayerID(Kokkos::ViewAllocateWithoutInitializing("LayerID"), nz, nx, ny);

        // Fill arrays with data from paraview file
        InitializeData(MicrostructureFile, nx, ny, nz, GrainID, LayerID);

        // Allocate memory for grain unit vectors, grain euler angles, RGB colors for IPF-Z coloring
        // (9*NumberOfOrientations,  3*NumberOfOrientations, and 3*NumberOfOrientations in size, respectively)
        int NumberOfOrientations = 10000;
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

        // Print information about 2 specific cross-sections: Z = 1 mm above the base plate, and X = midway through the
        // volume of the microstructure
        int NumberOfCrossSections = 2;
        std::vector<std::string> CrossSectionPlane = {"XY", "YZ"};
        int ZLocCrossSection = std::round((0.001 - XYZBounds[2]) / deltax);
        int XLocCrossSection = std::round(nx / 2);
        std::vector<int> CrossSectionLocation = {ZLocCrossSection, XLocCrossSection};
        std::vector<bool> PrintSectionPF = {true, true};
        std::vector<bool> PrintSectionIPF = {true, true};
        std::vector<bool> BimodalAnalysis = {true, false};
        std::string CSNameXY =
            "XY Cross-section located approx 1 mm above the baseplate (Z = " + std::to_string(ZLocCrossSection) + "):";
        std::string CSNameYZ =
            "YZ cross-section located through simulation center (X = " + std::to_string(XLocCrossSection) + "):";
        std::vector<std::string> CSLabels = {CSNameXY, CSNameYZ};
        PrintCrossSectionData(NumberOfCrossSections, BaseFileName, CrossSectionPlane, CrossSectionLocation, nx, ny, nz,
                              NumberOfOrientations, GrainID, PrintSectionPF, PrintSectionIPF, BimodalAnalysis, true,
                              deltax, GrainUnitVector_Host, GrainEulerAngles_Host, GrainRGBValues_Host, CSLabels);
    }
    // Finalize kokkos and end program
    Kokkos::finalize();

    return 0;
}
