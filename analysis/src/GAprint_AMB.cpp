// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT
#include "GAprint_AMB.hpp"
#include "CAfunctions.hpp"
#include "GAutils.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Helper function for PrintDefaultCrossSections_AMB and PrintCrossSectionSeries_AMB
void WriteEulerAngleData_AMB(std::ofstream &Grainplot, ViewI3D_H &GrainID, int i, int j, int k,
                             int NumberOfOrientations, double deltax, ViewF_H &GrainEulerAngles,
                             int &NucleatedGrainCells) {
    int GOVal = (abs(GrainID(k, i, j)) - 1) % NumberOfOrientations;
    // The grain structure is phase "1" - any unindexed points (which are possible from regions
    // that didn't undergo melting) are assigned phase "0"
    if (GOVal == -1)
        Grainplot << "0 0 0 0 " << i * deltax * pow(10, 6) << " " << j * deltax * pow(10, 6) << std::endl;
    else
        Grainplot << GrainEulerAngles(3 * GOVal) << " " << GrainEulerAngles(3 * GOVal + 1) << " "
                  << GrainEulerAngles(3 * GOVal + 2) << " 1 " << i * deltax * pow(10, 6) << " "
                  << j * deltax * pow(10, 6) << std::endl;
    // Count number of cells in this cross-section have GrainID < 0 (grains formed via nucleation)
    if (GrainID(k, i, j) < 0)
        NucleatedGrainCells++;
}

// For every 4 cells, print XY cross-section data to be plotted with the inverse pole figure colormaps or as pole
// figures
void PrintXYCrossSectionSeries_AMB(std::string BaseFileName, int nx, int ny, int nz, int NumberOfOrientations,
                                   ViewI3D_H &GrainID, double deltax, ViewF_H &GrainEulerAngles, double ZMin) {

    std::string OutputFilenameNGF = BaseFileName + "_FractNucGrains.csv";
    std::ofstream NGF;
    NGF.open(OutputFilenameNGF);
    NGF << "Z, FractNucleatedGrains" << std::endl;
    for (int k = 0; k < nz; k += 4) {
        std::cout << "Printing cross-section at k = " << k << std::endl;
        std::string OutputFilename = BaseFileName + "-" + std::to_string(k) + "_XYCrossSection.txt";
        std::ofstream Grainplot;
        Grainplot.open(OutputFilename);
        int NucleatedGrainCells = 0;
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                WriteEulerAngleData_AMB(Grainplot, GrainID, i, j, k, NumberOfOrientations, deltax, GrainEulerAngles,
                                        NucleatedGrainCells);
            }
        }
        double FractNucGrains = static_cast<double>(NucleatedGrainCells) / static_cast<double>(nx * ny);
        NGF << ZMin + deltax * k << ", " << FractNucGrains << std::endl;
        std::cout << "The fraction of the cross-sectional area at Z = " << ZMin + deltax * k
                  << " consisting of nucleated grains is: " << FractNucGrains << std::endl;
        Grainplot.close();
    }
    NGF.close();
}

// Print XY cross-section data corresponding to Z = 1 mm in the build, print
// - Area fraction consisting of nucleated grains
// - Area fraction of grains < 1500 sq microns ("small" grains)
// - Average area of small, large grains
// - Area of each small, large grain
// - Average misorientation (relative to X, Y, and Z) of small, large grains
// - Misorientation of each small, large grain

// For the centermost YZ cross-section of the build, print
// - Fraction of cross-section consisting of nucleated grains
// - Area of each grain
// - Average grain area

// Grains consisting of 2 or fewer pixels will be filtered out
void PrintDefaultCrossSectionData_AMB(std::string BaseFileName, int nx, int ny, int nz, int NumberOfOrientations,
                                      ViewI3D_H GrainID, double deltax, ViewF_H GrainEulerAngles,
                                      ViewF_H GrainUnitVector, double ZMin) {

    // Open file that will contain data used for ExaAM calibration/experimental comparison
    std::ofstream AMB_QoIs;
    std::string OutputFilenameAMBQoIs = BaseFileName + "_AMBQoI.txt";
    AMB_QoIs.open(OutputFilenameAMBQoIs);

    // What Z coordinate corresponds to 1 mm above the baseplate?
    int ZLocCrossSection = std::round((0.001 - ZMin) / deltax);
    std::cout << "XY cross-section at approx 1 mm above the baseplate is located at Z = " << ZLocCrossSection
              << std::endl;
    AMB_QoIs << "For the XY cross-section at approx 1 mm above the baseplate..." << std::endl;
    std::cout << "For the XY cross-section at approx 1 mm above the baseplate..." << std::endl;

    // Print MTEX data for the the XY cross-section, and store list of GrainIDs in this cross-section
    std::string XYMTEXFile = BaseFileName + "_XYMTEX.txt";
    std::ofstream XYMTEX;
    XYMTEX.open(XYMTEXFile);
    int NucleatedGrainCellsXY = 0;
    std::vector<int> XYGrainIDs(nx * ny);
    int Counter = 0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            WriteEulerAngleData_AMB(XYMTEX, GrainID, i, j, ZLocCrossSection, NumberOfOrientations, deltax,
                                    GrainEulerAngles, NucleatedGrainCellsXY);
            XYGrainIDs[Counter] = GrainID(ZLocCrossSection, i, j);
            Counter++;
        }
    }
    double FractNucGrainsXY = static_cast<double>(NucleatedGrainCellsXY) / static_cast<double>(nx * ny);
    std::cout << "Area fraction of nucleated grains is: " << FractNucGrainsXY << std::endl;
    AMB_QoIs << "Area fraction of nucleated grains is: " << FractNucGrainsXY << std::endl;
    XYMTEX.close();

    // Make list of unique grains in this XY cross-section
    std::vector<int> XYUniqueGrainIDs = FindUniqueGrains(XYGrainIDs);
    int XYNumberOfGrains = XYUniqueGrainIDs.size();

    // Make a list of grain areas in this XY cross-section
    std::vector<int> XYGrainAreas(XYNumberOfGrains, 0);
    for (int n = 0; n < XYNumberOfGrains; n++) {
        for (int nn = 0; nn < nx * ny; nn++) {
            if (XYUniqueGrainIDs[n] == XYGrainIDs[nn])
                XYGrainAreas[n]++;
        }
    }

    // Make list of grains < 5 cells ("too small to include in statistics"), larger than 1500 sq microns ("large"), and
    // in between ("small")
    int MinGrainSize = 5; // in cells
    int SmallLargeCutoff = 1500 / (deltax * deltax * pow(10, 12));
    int NumTooSmallGrains = 0, NumSmallGrains = 0, NumLargeGrains = 0;
    int AreaTooSmallGrains = 0, AreaSmallGrains = 0, AreaLargeGrains = 0;
    for (int n = 0; n < XYNumberOfGrains; n++) {
        if (XYGrainAreas[n] < MinGrainSize) {
            AreaTooSmallGrains += XYGrainAreas[n];
            NumTooSmallGrains++;
        }
        else {
            if (XYGrainAreas[n] < SmallLargeCutoff) {
                AreaSmallGrains += XYGrainAreas[n];
                NumSmallGrains++;
            }
            else {
                AreaLargeGrains += XYGrainAreas[n];
                NumLargeGrains++;
            }
        }
    }

    // What fraction of the cross-sectional area consists of each type of grain?
    double FractAreaTooSmall = static_cast<double>(AreaTooSmallGrains) / static_cast<double>(nx * ny);
    double FractAreaSmall = static_cast<double>(AreaSmallGrains) / static_cast<double>(nx * ny);
    double FractAreaLarge = static_cast<double>(AreaLargeGrains) / static_cast<double>(nx * ny);
    std::cout << "Area fraction of grains too small to include in statistics ( < " << MinGrainSize
              << " cells in area or less): " << FractAreaTooSmall << std::endl;
    std::cout << "Area fraction of grains smaller than 1500 sq microns: " << FractAreaSmall << std::endl;
    AMB_QoIs << "Area fraction of grains smaller than 1500 sq microns: " << FractAreaSmall << std::endl;
    std::cout << "Area fraction of grains greater than or equal to 1500 sq microns: " << FractAreaLarge << std::endl;

    // What's the average area for small and large grains?
    double AvgAreaSmall = static_cast<double>(nx * ny * FractAreaSmall) / static_cast<double>(NumSmallGrains);
    double AvgAreaLarge = static_cast<double>(nx * ny * FractAreaLarge) / static_cast<double>(NumLargeGrains);
    std::cout << "Average area (in square microns) for small grains: " << AvgAreaSmall * deltax * deltax * pow(10, 12)
              << std::endl;
    std::cout << "Average area (in square microns) for large grains: " << AvgAreaLarge * deltax * deltax * pow(10, 12)
              << std::endl;
    AMB_QoIs << "Average area (in square microns) for small grains: " << AvgAreaSmall * deltax * deltax * pow(10, 12)
             << std::endl;
    AMB_QoIs << "Average area (in square microns) for large grains: " << AvgAreaLarge * deltax * deltax * pow(10, 12)
             << std::endl;

    // Print area of each small grain to a file, and area of each large grain to a file
    std::string XYSmallAreasFile = BaseFileName + "_XYSmallAreas.txt";
    std::string XYLargeAreasFile = BaseFileName + "_XYLargeAreas.txt";
    std::ofstream XYSmallAreas, XYLargeAreas;
    XYSmallAreas.open(XYSmallAreasFile);
    XYLargeAreas.open(XYLargeAreasFile);
    for (int n = 0; n < XYNumberOfGrains; n++) {
        if (XYGrainAreas[n] >= MinGrainSize) {
            if (XYGrainAreas[n] < SmallLargeCutoff)
                XYSmallAreas << static_cast<double>(XYGrainAreas[n] * deltax * deltax * pow(10, 12)) << std::endl;
            else
                XYLargeAreas << static_cast<double>(XYGrainAreas[n] * deltax * deltax * pow(10, 12)) << std::endl;
        }
    }
    XYSmallAreas.close();
    XYLargeAreas.close();

    // Print misorientation relative to the X, Y, and Z directions for small and large grains, weighted by grain area
    std::string MisorientationFileSmallX = BaseFileName + "_XMisorientationSmallAreas.txt";
    std::string MisorientationFileSmallY = BaseFileName + "_YMisorientationSmallAreas.txt";
    std::string MisorientationFileSmallZ = BaseFileName + "_ZMisorientationSmallAreas.txt";
    std::string MisorientationFileLargeX = BaseFileName + "_XMisorientationLargeAreas.txt";
    std::string MisorientationFileLargeY = BaseFileName + "_YMisorientationLargeAreas.txt";
    std::string MisorientationFileLargeZ = BaseFileName + "_ZMisorientationLargeAreas.txt";
    std::ofstream MisorientationSmallX, MisorientationSmallY, MisorientationSmallZ, MisorientationLargeX,
        MisorientationLargeY, MisorientationLargeZ;
    MisorientationSmallX.open(MisorientationFileSmallX);
    MisorientationSmallY.open(MisorientationFileSmallY);
    MisorientationSmallZ.open(MisorientationFileSmallZ);
    MisorientationLargeX.open(MisorientationFileLargeX);
    MisorientationLargeY.open(MisorientationFileLargeY);
    MisorientationLargeZ.open(MisorientationFileLargeZ);

    ViewF_H GrainMisorientationX = MisorientationCalc(NumberOfOrientations, GrainUnitVector, 0);
    ViewF_H GrainMisorientationY = MisorientationCalc(NumberOfOrientations, GrainUnitVector, 1);
    ViewF_H GrainMisorientationZ = MisorientationCalc(NumberOfOrientations, GrainUnitVector, 2);
    float GrainMisorientation_small_sumx = 0, GrainMisorientation_small_sumy = 0, GrainMisorientation_small_sumz = 0;
    float GrainMisorientation_large_sumx = 0, GrainMisorientation_large_sumy = 0, GrainMisorientation_large_sumz = 0;
    for (int n = 0; n < XYNumberOfGrains; n++) {
        if (XYGrainAreas[n] >= MinGrainSize) {

            int ThisGrainOrientation = getGrainOrientation(XYUniqueGrainIDs[n], NumberOfOrientations);
            if (ThisGrainOrientation < 0)
                throw std::runtime_error(
                    "Error analyzing grain misorientations: GrainID = 0 grain is present in the XY cross-section");
            float ThisGrainMisorientationX = GrainMisorientationX(ThisGrainOrientation);
            float ThisGrainMisorientationY = GrainMisorientationY(ThisGrainOrientation);
            float ThisGrainMisorientationZ = GrainMisorientationZ(ThisGrainOrientation);

            if (XYGrainAreas[n] < SmallLargeCutoff) {
                GrainMisorientation_small_sumx += XYGrainAreas[n] * ThisGrainMisorientationX;
                GrainMisorientation_small_sumy += XYGrainAreas[n] * ThisGrainMisorientationY;
                GrainMisorientation_small_sumz += XYGrainAreas[n] * ThisGrainMisorientationZ;
                for (int nn = 0; nn < XYGrainAreas[n]; nn++) {
                    MisorientationSmallX << ThisGrainMisorientationX << std::endl;
                    MisorientationSmallY << ThisGrainMisorientationY << std::endl;
                    MisorientationSmallZ << ThisGrainMisorientationZ << std::endl;
                }
            }
            else {
                GrainMisorientation_large_sumx += XYGrainAreas[n] * ThisGrainMisorientationX;
                GrainMisorientation_large_sumy += XYGrainAreas[n] * ThisGrainMisorientationY;
                GrainMisorientation_large_sumz += XYGrainAreas[n] * ThisGrainMisorientationZ;
                for (int nn = 0; nn < XYGrainAreas[n]; nn++) {
                    MisorientationLargeX << ThisGrainMisorientationX << std::endl;
                    MisorientationLargeY << ThisGrainMisorientationY << std::endl;
                    MisorientationLargeZ << ThisGrainMisorientationZ << std::endl;
                }
            }
        }
    }
    MisorientationSmallX.close();
    MisorientationSmallY.close();
    MisorientationSmallZ.close();
    MisorientationLargeX.close();
    MisorientationLargeY.close();
    MisorientationLargeZ.close();
    float AvgMisorientation_Small_X = GrainMisorientation_small_sumx / static_cast<float>(AreaSmallGrains);
    float AvgMisorientation_Small_Y = GrainMisorientation_small_sumy / static_cast<float>(AreaSmallGrains);
    float AvgMisorientation_Small_Z = GrainMisorientation_small_sumz / static_cast<float>(AreaSmallGrains);
    float AvgMisorientation_Large_X = GrainMisorientation_large_sumx / static_cast<float>(AreaLargeGrains);
    float AvgMisorientation_Large_Y = GrainMisorientation_large_sumy / static_cast<float>(AreaLargeGrains);
    float AvgMisorientation_Large_Z = GrainMisorientation_large_sumz / static_cast<float>(AreaLargeGrains);
    std::cout << "Average misorientation for small grains relative to the X direction: " << AvgMisorientation_Small_X
              << std::endl;
    std::cout << "Average misorientation for small grains relative to the Y direction: " << AvgMisorientation_Small_Y
              << std::endl;
    std::cout << "Average misorientation for small grains relative to the Z direction: " << AvgMisorientation_Small_Z
              << std::endl;
    std::cout << "Average misorientation for large grains relative to the X direction: " << AvgMisorientation_Large_X
              << std::endl;
    std::cout << "Average misorientation for large grains relative to the Y direction: " << AvgMisorientation_Large_Y
              << std::endl;
    std::cout << "Average misorientation for large grains relative to the Z direction: " << AvgMisorientation_Large_Z
              << std::endl;
    AMB_QoIs << "Average misorientation for small grains relative to the X direction: " << AvgMisorientation_Small_X
             << std::endl;
    AMB_QoIs << "Average misorientation for small grains relative to the Y direction: " << AvgMisorientation_Small_Y
             << std::endl;
    AMB_QoIs << "Average misorientation for small grains relative to the Z direction: " << AvgMisorientation_Small_Z
             << std::endl;
    AMB_QoIs << "Average misorientation for large grains relative to the X direction: " << AvgMisorientation_Large_X
             << std::endl;
    AMB_QoIs << "Average misorientation for large grains relative to the Y direction: " << AvgMisorientation_Large_Y
             << std::endl;
    AMB_QoIs << "Average misorientation for large grains relative to the Z direction: " << AvgMisorientation_Large_Z
             << std::endl;

    // YZ cross-section analysis
    // What X coordinate corresponds to the centerline of the simulation?
    int XLocCrossSection = std::round(nx / 2);
    std::cout << "YZ cross-section at center of simulation is located at X = " << XLocCrossSection << std::endl;
    AMB_QoIs << "For the YZ cross-section bisecting the simulation volume..." << std::endl;
    std::cout << "For the YZ cross-section bisecting the simulation volume..." << std::endl;

    // Print MTEX data for the the YZ cross-section
    std::string YZMTEXFile = BaseFileName + "_YZMTEX.txt";
    std::ofstream YZMTEX;
    YZMTEX.open(YZMTEXFile);
    int NucleatedGrainCellsYZ = 0;
    std::vector<int> YZGrainIDs(ny * nz);
    Counter = 0;
    for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
            WriteEulerAngleData_AMB(YZMTEX, GrainID, XLocCrossSection, j, k, NumberOfOrientations, deltax,
                                    GrainEulerAngles, NucleatedGrainCellsYZ);
            YZGrainIDs[Counter] = GrainID(k, XLocCrossSection, j);
            Counter++;
        }
    }
    double FractNucGrainsYZ = static_cast<double>(NucleatedGrainCellsYZ) / static_cast<double>(ny * nz);
    std::cout << "Area fraction of nucleated grains is: " << FractNucGrainsYZ << std::endl;
    AMB_QoIs << "Area fraction of nucleated grains is: " << FractNucGrainsYZ << std::endl;
    YZMTEX.close();

    // Make list of unique grains in this YZ cross-section
    std::vector<int> YZUniqueGrainIDs = FindUniqueGrains(YZGrainIDs);
    int YZNumberOfGrains = YZUniqueGrainIDs.size();

    // Make a list of grain areas in this YZ cross-section
    std::vector<int> YZGrainAreas(YZNumberOfGrains, 0);
    for (int n = 0; n < YZNumberOfGrains; n++) {
        for (int nn = 0; nn < ny * nz; nn++) {
            if (YZUniqueGrainIDs[n] == YZGrainIDs[nn])
                YZGrainAreas[n]++;
        }
    }

    // Count grain areas large enough to be counted
    int NumGrainsAboveThreshold = 0, AreaGrainsAboveThreshold = 0;
    for (int n = 0; n < YZNumberOfGrains; n++) {
        if (YZGrainAreas[n] >= MinGrainSize) {
            AreaGrainsAboveThreshold += YZGrainAreas[n];
            NumGrainsAboveThreshold++;
        }
    }

    // What fraction of the cross-sectional area consists of grains large enough to be counted?
    double FractAreaAboveThreshold = static_cast<double>(AreaGrainsAboveThreshold) / static_cast<double>(ny * nz);
    std::cout << "YZ area fraction of grains too small to include in statistics ( < " << MinGrainSize
              << " cells in area or less): " << 1 - FractAreaAboveThreshold << std::endl;

    // What's the average area for the grains large enough to be counted?
    double AvgAreaAboveThreshold =
        static_cast<double>(ny * nz * FractAreaAboveThreshold) / static_cast<double>(NumGrainsAboveThreshold);
    std::cout << "YZ average grain area (in square microns): " << AvgAreaAboveThreshold * deltax * deltax * pow(10, 12)
              << std::endl;
    AMB_QoIs << "YZ average grain area (in square microns): " << AvgAreaAboveThreshold * deltax * deltax * pow(10, 12)
             << std::endl;

    // Print area of each grain to a file, as a fraction of the total area filled by grains large enough to be
    // considered
    std::string YZAreasFile = BaseFileName + "_YZAreas.txt";
    std::ofstream YZAreas;
    YZAreas.open(YZAreasFile);
    for (int n = 0; n < YZNumberOfGrains; n++) {
        if (YZGrainAreas[n] > MinGrainSize) {
            YZAreas << static_cast<double>(YZGrainAreas[n]) / static_cast<double>(AreaGrainsAboveThreshold)
                    << std::endl;
        }
    }
    YZAreas.close();

    AMB_QoIs.close();
}
