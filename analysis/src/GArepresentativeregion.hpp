// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef GA_REPREGION_HPP
#define GA_REPREGION_HPP

#include "CAconfig.hpp"
#include "CAtypes.hpp"
#include "GAutils.hpp"

#include <Kokkos_Core.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifdef ExaCA_ENABLE_JSON
#include <nlohmann/json.hpp>
// Structure that holds the details of a representative region of the larger microstructure, along which the specific
// analyses to be performed
struct RepresentativeRegion {

    std::string regionType;
    std::string regionOrientation;
    std::string units;
    std::string units_dimension;
    std::vector<double> xBounds_Meters = std::vector<double>(2);
    std::vector<double> yBounds_Meters = std::vector<double>(2);
    std::vector<double> zBounds_Meters = std::vector<double>(2);
    double regionSize_Meters;
    double regionSize_Microns;
    std::vector<int> xBounds_Cells = std::vector<int>(2);
    std::vector<int> yBounds_Cells = std::vector<int>(2);
    std::vector<int> zBounds_Cells = std::vector<int>(2);
    int regionSize_Cells;

    // Analysis options for printing stats to the screen (_Stats), to the file of grainwise statistics (_PerGrainStats),
    // or to files of layerwise statistics (_LayerwiseStats)
    std::vector<std::string> AnalysisOptions_Stats_key = {
        "GrainTypeFractions",    // all regions
        "Misorientation",        // all regions
        "Size",                  // all regions
        "BuildTransAspectRatio", // volumes only
        "XExtent",               // all regions
        "YExtent",               // all regions
        "ZExtent",               // all regions
    };
    std::vector<bool> AnalysisOptions_StatsYN = std::vector<bool>(7, false);
    std::vector<std::string> AnalysisOptions_PerGrainStats_key = {
        "Misorientation",        // all regions
        "Size",                  // all regions - volume, area, or length
        "XExtent",               // all regions
        "YExtent",               // all regions
        "ZExtent",               // all regions
        "BuildTransAspectRatio", // volumes only
        "IPFZ-RGB"               // all regions
    };
    std::vector<bool> AnalysisOptions_PerGrainStatsYN = std::vector<bool>(5, false);
    std::vector<std::string> AnalysisOptions_LayerwiseStats_key = {
        "MeanGrainArea",        // volume only
        "MeanWeightedGrainArea" // volume only
    };
    bool PrintPerGrainStatsYN;
    std::vector<bool> AnalysisOptions_LayerwiseStatsYN = std::vector<bool>(2, false);

    // Analysis options that print separate files
    bool PrintExaConstitYN, PrintPoleFigureYN, PrintInversePoleFigureMapYN;

    // Constructor
    RepresentativeRegion(nlohmann::json RegionData, int nx, int ny, int nz, double deltax,
                         std::vector<double> XYZBounds) {

        // Are the bounds given in cells, or in microns? Store both representations
        setUnits(RegionData);
        // Obtain the bounds of the region in x, y, and z, in both cells and microns
        ConvertBounds(RegionData, deltax, XYZBounds, nx, ny, nz);
        // Deduce region type/orientation/size units from the bounds given
        setRegionTypeOrientation();
        setRegionSize(deltax);
        setUnitDimension();

        // Check which overall stats and per grain stats should be printed for this region
        ReadAnalysisOptionsFromList(RegionData, "printStats", AnalysisOptions_Stats_key, AnalysisOptions_StatsYN);
        // PrintPerGrainStatsYN = true if any one of the options are toggled
        ReadAnalysisOptionsFromList(RegionData, "printPerGrainStats", AnalysisOptions_PerGrainStats_key,
                                    AnalysisOptions_PerGrainStatsYN);
        int NumAnalysisOptions_PerGrainStats = AnalysisOptions_PerGrainStatsYN.size();
        for (int n = 0; n < NumAnalysisOptions_PerGrainStats; n++) {
            if (AnalysisOptions_PerGrainStatsYN[n])
                PrintPerGrainStatsYN = true;
        }
        std::cout << "Read analysis options" << std::endl;

        // Layerwise stats are for volumes only
        if (regionType == "volume")
            ReadAnalysisOptionsFromList(RegionData, "printLayerwiseData", AnalysisOptions_LayerwiseStats_key,
                                        AnalysisOptions_LayerwiseStatsYN);
        // Check other Y/N options (false by default or if not an allowed option for the region type)
        ReadSeparateFileAnalysisOptions(RegionData);
    }

    void ReadSeparateFileAnalysisOptions(nlohmann::json RegionData) {

        // ExaConstit print/area series prints for volumes only
        if (regionType == "volume") {
            if (RegionData.contains("printExaConstit")) {
                PrintExaConstitYN = RegionData["printExaConstit"];
            }
            else
                PrintExaConstitYN = false;
        }
        else
            PrintExaConstitYN = false;

        // Pole figure data print for volumes and areas only
        if ((regionType == "volume") || (regionType == "area")) {
            if (RegionData.contains("printPoleFigure")) {
                PrintPoleFigureYN = RegionData["printPoleFigure"];
            }
            else
                PrintPoleFigureYN = false;
        }
        else
            PrintPoleFigureYN = false;

        // IPF map data for areas only
        if (regionType == "area") {
            if (RegionData.contains("printInversePoleFigureMap")) {
                PrintInversePoleFigureMapYN = RegionData["printInversePoleFigureMap"];
            }
            else
                PrintInversePoleFigureMapYN = false;
        }
        else
            PrintInversePoleFigureMapYN = false;
    }

    void setRegionTypeOrientation() {
        bool FlatX = (xBounds_Cells[0] == xBounds_Cells[1]);
        bool FlatY = (yBounds_Cells[0] == yBounds_Cells[1]);
        bool FlatZ = (zBounds_Cells[0] == zBounds_Cells[1]);
        if ((!FlatX) && (!FlatY) && (!FlatZ)) {
            regionType = "volume";
            regionOrientation = "XYZ";
        }
        else if ((FlatX) && (FlatY) && (FlatZ))
            throw std::runtime_error("Error: Cannot analyze a single cell");
        else {
            // Either an area or a length
            if ((FlatX) && (FlatY)) {
                regionType = "length";
                regionOrientation = "Z";
            }
            else if ((FlatX) && (FlatZ)) {
                regionType = "length";
                regionOrientation = "Y";
            }
            else if ((FlatY) && (FlatZ)) {
                regionType = "length";
                regionOrientation = "X";
            }
            else {
                regionType = "area";
                if (FlatX)
                    regionOrientation = "YZ";
                else if (FlatY)
                    regionOrientation = "XZ";
                else if (FlatZ)
                    regionOrientation = "XY";
            }
        }
    }

    void ReadAnalysisOptionsFromList(nlohmann::json RegionData, std::string AnalysisType,
                                     std::vector<std::string> AnalysisKey, std::vector<bool> &AnalysisYN) {
        int NumOptionsRead = RegionData[AnalysisType].size();
        int NumOptionsTotal = AnalysisYN.size();
        for (int i = 0; i < NumOptionsRead; i++) {
            std::string ValueRead = RegionData[AnalysisType][i];
            for (int j = 0; j < NumOptionsTotal; j++) {
                std::string ValueKey = AnalysisKey[j];
                if (ValueRead == AnalysisKey[j]) {
                    AnalysisYN[j] = true;
                }
            }
        }
    }

    void setUnits(nlohmann::json RegionData) {
        // Check that the units provided are Meters or Cells, throwing an error otherwise
        units = RegionData["units"];
        if (units != "Meters" && units != "meters" && units != "Cells" && units != "cells")
            throw std::runtime_error("Error: Invalid units, should be Meters or Cells");
    }

    // Get the lower and upper bounds for x, y, or z, in microns, converting from whichever format (cells or microns)
    // was present in the analysis input file
    void ConvertBounds(nlohmann::json RegionData, double deltax, std::vector<double> XYZBounds, int nx, int ny,
                       int nz) {

        std::vector<std::string> BoundDirections = {"x", "y", "z"};
        for (int dir = 0; dir < 3; dir++) {

            // Initialize bounds in both Meters and Cells
            std::vector<int> Bounds_Cells(2);
            std::vector<double> Bounds_Meters(2);

            // Get the minimium coordinate along the specified axis
            double CoordinateMinVal = XYZBounds[dir];

            // Two bounds or no bounds should be given
            // If no bounds are given, defaults to the size of the domain along that direction
            std::string twobounds = BoundDirections[dir] + "Bounds";
            if (!(RegionData.contains(twobounds))) {
                Bounds_Cells[0] = 0;
                Bounds_Meters[0] = CoordinateMinVal;
                if (dir == 0) {
                    Bounds_Cells[1] = nx - 1;
                    Bounds_Meters[1] = XYZBounds[3];
                }
                else if (dir == 1) {
                    Bounds_Cells[1] = ny - 1;
                    Bounds_Meters[1] = XYZBounds[4];
                }
                else if (dir == 2) {
                    Bounds_Cells[1] = nz - 1;
                    Bounds_Meters[1] = XYZBounds[5];
                }
            }
            else {
                std::vector<double> Bounds_Read(2);
                Bounds_Read[0] = RegionData[twobounds][0];
                Bounds_Read[1] = RegionData[twobounds][1];
                if ((units == "Meters") || (units == "meters")) {
                    // Bounds provided were in meters, convert to cells
                    Bounds_Meters[0] = Bounds_Read[0];
                    Bounds_Meters[1] = Bounds_Read[1];
                    Bounds_Cells[0] = std::round((Bounds_Read[0] - CoordinateMinVal) / deltax);
                    Bounds_Cells[1] = std::round((Bounds_Read[1] - CoordinateMinVal) / deltax);
                }
                else if ((units == "Cells") || (units == "cells")) {
                    // Bounds provided were in cells, convert to microns
                    Bounds_Cells[0] = static_cast<int>(Bounds_Read[0]);
                    Bounds_Cells[1] = static_cast<int>(Bounds_Read[1]);
                    Bounds_Meters[0] = CoordinateMinVal + deltax * Bounds_Read[0];
                    Bounds_Meters[1] = CoordinateMinVal + deltax * Bounds_Read[1];
                }
            }

            // Store bounds for the given direction
            for (int i = 0; i < 2; i++) {
                if (dir == 0) {
                    xBounds_Meters[i] = Bounds_Meters[i];
                    xBounds_Cells[i] = Bounds_Cells[i];
                }
                else if (dir == 1) {
                    yBounds_Meters[i] = Bounds_Meters[i];
                    yBounds_Cells[i] = Bounds_Cells[i];
                }
                else if (dir == 2) {
                    zBounds_Meters[i] = Bounds_Meters[i];
                    zBounds_Cells[i] = Bounds_Cells[i];
                }
            }
        }
    }

    // Get the size of the region in cells, and in meters/microns (either 1D, 2D, or 3D units depending on region type)
    void setRegionSize(double deltax) {
        regionSize_Cells = (xBounds_Cells[1] - xBounds_Cells[0] + 1) * (yBounds_Cells[1] - yBounds_Cells[0] + 1) *
                           (zBounds_Cells[1] - zBounds_Cells[0] + 1);
        if (regionType == "length") {
            regionSize_Meters = regionSize_Cells * deltax;
            regionSize_Microns = regionSize_Meters * pow(10, 6);
        }
        else if (regionType == "area") {
            regionSize_Meters = regionSize_Cells * pow(deltax, 2);
            regionSize_Microns = regionSize_Meters * pow(10, 12);
        }
        else if (regionType == "volume") {
            regionSize_Meters = regionSize_Cells * pow(deltax, 3);
            regionSize_Microns = regionSize_Meters * pow(10, 18);
        }
        else
            throw std::runtime_error("Error: Invalid region type in GetRegionSize during region construction");
    }

    // Set the appropriate units for the region based on the dimensions (in microns)
    void setUnitDimension() {
        regionSize_Cells = (xBounds_Cells[1] - xBounds_Cells[0] + 1) * (yBounds_Cells[1] - yBounds_Cells[0] + 1) *
                           (zBounds_Cells[1] - zBounds_Cells[0] + 1);
        if (regionType == "length")
            units_dimension = "microns";
        else if (regionType == "area")
            units_dimension = "square microns";
        else if (regionType == "volume")
            units_dimension = "cubic microns";
        else
            throw std::runtime_error("Error: Invalid region type in setRegionSize during region construction");
    }

    // Get the list of Grain IDs associated with the representative region
    std::vector<int> getGrainIDVector(ViewI3D_H GrainID) {

        std::vector<int> GrainIDVector(regionSize_Cells);
        int count = 0;
        for (int k = zBounds_Cells[0]; k <= zBounds_Cells[1]; k++) {
            for (int i = xBounds_Cells[0]; i <= xBounds_Cells[1]; i++) {
                for (int j = yBounds_Cells[0]; j <= yBounds_Cells[1]; j++) {
                    GrainIDVector[count] = GrainID(k, i, j);
                    count++;
                }
            }
        }
        return GrainIDVector;
    }

    // Get the list of unique Grain IDs associated with the representative region
    std::vector<int> getUniqueGrainIDVector(std::vector<int> GrainIDVector) {

        std::vector<int> UniqueGrainIDVector = GrainIDVector;
        std::sort(UniqueGrainIDVector.begin(), UniqueGrainIDVector.end());
        std::vector<int>::iterator it;
        it = std::unique(UniqueGrainIDVector.begin(), UniqueGrainIDVector.end());
        UniqueGrainIDVector.resize(std::distance(UniqueGrainIDVector.begin(), it));
        return UniqueGrainIDVector;
    }

    // Given an input vector of integer Grain ID values "GrainIDVector", and an input vector of the
    // unique Grain ID values "UniqueGrainIDVector" (of size "NumberOfGrains"), return a third vector "GrainSizeVector"
    // listing the size of each of the "NumberOfGrains" grains (scaled by the cell size deltax and depending on whether
    // the region is 1D, 2D, or 3D)
    std::vector<float> getGrainSizeVector(const std::vector<int> GrainIDVector,
                                          const std::vector<int> UniqueGrainIDVector, const int NumberOfGrains,
                                          double deltax) {

        std::vector<float> GrainSizeVector(NumberOfGrains);
        double conv = convertToMicrons(deltax, regionType);
        for (int n = 0; n < NumberOfGrains; n++) {
            int GrainSizeCells = std::count(GrainIDVector.begin(), GrainIDVector.end(), UniqueGrainIDVector[n]);
            // convert to either microns, square microns, or cubic microns
            GrainSizeVector[n] = conv * GrainSizeCells;
        }
        return GrainSizeVector;
    }

    // TODO: Consolidate volume/area header print routines from GAprint when non-JSON option is removed
    // Print header for grain statistics file for a representative volume
    void printAnalysisHeader_Volume(std::ofstream &QoIs, std::string RegionName) {

        std::cout << "Stats for " << RegionName << " volume:" << std::endl;
        QoIs << "Stats for " << RegionName << " volume:" << std::endl;
        std::cout << "The representative volume specified is bounded by X = [" << xBounds_Meters[0] << ","
                  << xBounds_Meters[1] << "], Y = [" << yBounds_Meters[0] << "," << yBounds_Meters[1] << "], and Z = ["
                  << zBounds_Meters[0] << "," << zBounds_Meters[1] << "] m" << std::endl;
        std::cout << "The representative volume specified is bounded by cells spanning X = [" << xBounds_Cells[0] << ","
                  << xBounds_Cells[1] << "], Y = [" << yBounds_Cells[0] << "," << yBounds_Cells[1] << "], and Z = ["
                  << zBounds_Cells[0] << "," << zBounds_Cells[1] << "]" << std::endl;
        QoIs << "The representative volume specified is bounded by X = [" << xBounds_Meters[0] << ","
             << xBounds_Meters[1] << "], Y = [" << yBounds_Meters[0] << "," << yBounds_Meters[1] << "], and Z = ["
             << zBounds_Meters[0] << "," << zBounds_Meters[1] << "] m" << std::endl;
        QoIs << "The representative volume specified is bounded by cells spanning X = [" << xBounds_Cells[0] << ","
             << xBounds_Cells[1] << "], Y = [" << yBounds_Cells[0] << "," << yBounds_Cells[1] << "], and Z = ["
             << zBounds_Cells[0] << "," << zBounds_Cells[1] << "]" << std::endl;
    }

    // Write pole figure data for this region to a file to be read by MTEX
    void writePoleFigure(std::string BaseFileNameThisRegion, int NumberOfOrientations, ViewF_H GrainEulerAngles,
                         ViewI_H GOHistogram) {

        // Using new format, write pole figure data to "Filename"
        std::string Filename = BaseFileNameThisRegion + "_PoleFigureData.txt";
        std::ofstream GrainplotPF;
        GrainplotPF.open(Filename);
        GrainplotPF << "% MTEX ODF" << std::endl;
        GrainplotPF << "% crystal symmetry: \"m-3m\"" << std::endl;
        GrainplotPF << "% specimen symmetry: \"43\"" << std::endl;
        GrainplotPF << "% phi1    Phi     phi2    value" << std::endl;
        GrainplotPF << std::fixed << std::setprecision(6);
        for (int i = 0; i < NumberOfOrientations; i++) {
            GrainplotPF << GrainEulerAngles(3 * i) << " " << GrainEulerAngles(3 * i + 1) << " "
                        << GrainEulerAngles(3 * i + 2) << " " << (float)(GOHistogram(i)) << std::endl;
        }
        GrainplotPF.close();
    }

    // Print data to be read by MTEX to plot the cross-section using the inverse pole figure colormap
    void writeIPFColoredCrossSection(std::string BaseFileNameThisRegion, ViewI3D_H GrainID, ViewF_H GrainEulerAngles,
                                     double deltax, int NumberOfOrientations) {

        // What portion of the area is to be printed?
        int Index1Low, Index1High, Index2Low, Index2High, CrossSectionOutOfPlaneLocation;
        if (regionOrientation == "XY") {
            Index1Low = xBounds_Cells[0];
            Index1High = xBounds_Cells[1];
            Index2Low = yBounds_Cells[0];
            Index2High = yBounds_Cells[1];
            CrossSectionOutOfPlaneLocation = zBounds_Cells[0];
        }
        else if (regionOrientation == "XZ") {
            Index1Low = xBounds_Cells[0];
            Index1High = xBounds_Cells[1];
            Index2Low = zBounds_Cells[0];
            Index2High = zBounds_Cells[1];
            CrossSectionOutOfPlaneLocation = yBounds_Cells[0];
        }
        else if (regionOrientation == "YZ") {
            Index1Low = yBounds_Cells[0];
            Index1High = yBounds_Cells[1];
            Index2Low = zBounds_Cells[0];
            Index2High = zBounds_Cells[1];
            CrossSectionOutOfPlaneLocation = xBounds_Cells[0];
        }
        else
            throw std::runtime_error(
                "Error: Invalid plane type in writeIPFColoredCrossSection: should be XY, XZ, or YZ");

        std::string FNameIPF = BaseFileNameThisRegion + "_IPFCrossSectionData.txt";
        std::ofstream GrainplotIPF;
        GrainplotIPF.open(FNameIPF);
        GrainplotIPF << std::fixed << std::setprecision(6);
        for (int Index1 = Index1Low; Index1 < Index1High; Index1++) {
            for (int Index2 = Index2Low; Index2 < Index2High; Index2++) {
                // How do Index1, Index2, and out of plane location correspond to GrainID(Z loc, X loc, Yloc)?
                int ZLoc, XLoc, YLoc;
                if (regionOrientation == "XY") {
                    XLoc = Index1;
                    YLoc = Index2;
                    ZLoc = CrossSectionOutOfPlaneLocation;
                }
                else if (regionOrientation == "YZ") {
                    XLoc = CrossSectionOutOfPlaneLocation;
                    YLoc = Index1;
                    ZLoc = Index2;
                }
                else if (regionOrientation == "XZ") {
                    XLoc = Index1;
                    YLoc = CrossSectionOutOfPlaneLocation;
                    ZLoc = Index2;
                }
                else
                    throw std::runtime_error("Error: Unknown regionOrientation input for WriteIPFColoredCrossSection: "
                                             "should be XY, YZ, or XZ");
                // What orientation does this grain id correspond to? Should be between 0 and NumberOfOrientations-1
                int GOVal = (abs(GrainID(ZLoc, XLoc, YLoc)) - 1) % NumberOfOrientations;
                // The grain structure is phase "1" - any unindexed points with GOVal = -1 (which are possible from
                // regions that didn't undergo melting) are assigned phase "0"
                if (GOVal == -1)
                    GrainplotIPF << "0 0 0 0 " << Index1 * deltax * pow(10, 6) << " " << Index2 * deltax * pow(10, 6)
                                 << std::endl;
                else
                    GrainplotIPF << GrainEulerAngles(3 * GOVal) << " " << GrainEulerAngles(3 * GOVal + 1) << " "
                                 << GrainEulerAngles(3 * GOVal + 2) << " 1 " << Index1 * deltax * pow(10, 6) << " "
                                 << Index2 * deltax * pow(10, 6) << std::endl;
            }
        }
        GrainplotIPF.close();
    }

    // Write data for the specified RVE for ExaConstit
    void writeExaConstitRVE(std::string BaseFileNameThisRegion, double deltax, ViewI3D_H GrainID) {

        int XLow = xBounds_Cells[0];
        int XHigh = xBounds_Cells[1];
        int YLow = yBounds_Cells[0];
        int YHigh = yBounds_Cells[1];
        int ZLow = zBounds_Cells[0];
        int ZHigh = zBounds_Cells[1];
        std::string FName = BaseFileNameThisRegion + "_ExaConstit.csv";
        std::ofstream Grainplot;
        Grainplot.open(FName);
        Grainplot << "Coordinates are in CA units, 1 cell = " << deltax << " m. Data is cell-centered. Origin at "
                  << XLow << "," << YLow << "," << ZLow << " , domain size is " << XHigh - XLow + 1 << " by "
                  << YHigh - YLow + 1 << " by " << ZHigh - ZLow + 1 << " cells" << std::endl;
        Grainplot << "X coord, Y coord, Z coord, Grain ID" << std::endl;
        int NucleatedGrainCells = 0;
        for (int k = ZLow; k <= ZHigh; k++) {
            for (int j = YLow; j <= YHigh; j++) {
                for (int i = XLow; i <= XHigh; i++) {
                    Grainplot << i << "," << j << "," << k << "," << GrainID(k, i, j) << std::endl;
                    if (GrainID(k, i, j) < 0)
                        NucleatedGrainCells++;
                }
            }
        }
        Grainplot.close();
    }

    // Write a csv file of stats for each grain
    void writePerGrainStats(std::string OutputFileName, std::vector<int> UniqueGrainIDVector,
                            std::vector<float> GrainMisorientationXVector,
                            std::vector<float> GrainMisorientationYVector,
                            std::vector<float> GrainMisorientationZVector, std::vector<float> GrainSizeVector,
                            std::vector<float> GrainExtentX, std::vector<float> GrainExtentY,
                            std::vector<float> GrainExtentZ, std::vector<float> BuildTransAspectRatio,
                            int NumberOfGrains, std::vector<float> GrainRed, std::vector<float> GrainGreen,
                            std::vector<float> GrainBlue) {

        // Which quantities should be printed?
        std::string stats_fname = OutputFileName + "_grains.csv";
        std::ofstream GrainStats;
        GrainStats.open(stats_fname);
        GrainStats << "GrainID";
        // Which quantities should be printed for each grain?
        if (AnalysisOptions_PerGrainStatsYN[0])
            GrainStats << ",misorientationX,misorientationY,misorientationZ";
        if (AnalysisOptions_PerGrainStatsYN[1])
            GrainStats << "," << regionType;
        if (AnalysisOptions_PerGrainStatsYN[2])
            GrainStats << ",extentX";
        if (AnalysisOptions_PerGrainStatsYN[3])
            GrainStats << ",extentY";
        if (AnalysisOptions_PerGrainStatsYN[4])
            GrainStats << ",extentZ";
        if (AnalysisOptions_PerGrainStatsYN[5])
            GrainStats << ",buildtransAR";
        if (AnalysisOptions_PerGrainStatsYN[6])
            GrainStats << ",IPFZ_r,IPFZ_g,IPFZ_b";
        GrainStats << std::endl;

        // Print the specified quantities to the csv file
        for (int n = 0; n < NumberOfGrains; n++) {
            GrainStats << UniqueGrainIDVector[n];
            if (AnalysisOptions_PerGrainStatsYN[0])
                GrainStats << "," << GrainMisorientationXVector[n] << "," << GrainMisorientationYVector[n] << ","
                           << GrainMisorientationZVector[n];
            if (AnalysisOptions_PerGrainStatsYN[1])
                GrainStats << "," << GrainSizeVector[n];
            if (AnalysisOptions_PerGrainStatsYN[2])
                GrainStats << "," << GrainExtentX[n];
            if (AnalysisOptions_PerGrainStatsYN[3])
                GrainStats << "," << GrainExtentY[n];
            if (AnalysisOptions_PerGrainStatsYN[4])
                GrainStats << "," << GrainExtentZ[n];
            if (AnalysisOptions_PerGrainStatsYN[5])
                GrainStats << "," << BuildTransAspectRatio[n];
            if (AnalysisOptions_PerGrainStatsYN[6])
                GrainStats << "," << GrainRed[n] << "," << GrainGreen[n] << "," << GrainBlue[n];
            GrainStats << std::endl;
        }
        GrainStats.close();
    }
};

#endif

#endif
