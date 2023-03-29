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
        std::string Units = RegionData["units"];
        // Obtain the bounds of the region in x, y, and z, in both cells and microns
        ConvertBounds(RegionData, Units, deltax, XYZBounds, nx, ny, nz);
        // Deduce region type/orientation from the bounds given
        GetRegionTypeOrientation();
        GetRegionSize(deltax);

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

    void GetRegionTypeOrientation() {
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

    // Get the lower and upper bounds for x, y, or z, in microns, converting from whichever format (cells or microns)
    // was present in the analysis input file
    void ConvertBounds(nlohmann::json RegionData, std::string Units, double deltax, std::vector<double> XYZBounds,
                       int nx, int ny, int nz) {

        std::vector<std::string> BoundDirections = {"x", "y", "z"};
        for (int dir = 0; dir < 3; dir++) {
            // Get the minimium coordinate along the specified axis
            double CoordinateMinVal = XYZBounds[dir];

            // Next, check if just one bound is given, or two bounds
            // If no bounds are given, defaults to the size of the domain along that direction
            std::string boundtype = BoundDirections[dir];
            std::vector<double> Bounds_Read(2);
            std::string onebound = boundtype + "Bound";
            std::string twobounds = boundtype + "Bounds";
            if ((RegionData.contains(onebound)) && (RegionData.contains(twobounds))) {
                std::string error =
                    "Error: region for analysis cannot have values for both " + onebound + " and " + twobounds;
                throw std::runtime_error(error);
            }
            else if ((!RegionData.contains(onebound)) && (!RegionData.contains(twobounds))) {
                Bounds_Read[0] = 0;
                if (Units == "Meters") {
                    if (dir == 0)
                        Bounds_Read[1] = XYZBounds[3];
                    else if (dir == 1)
                        Bounds_Read[1] = XYZBounds[4];
                    else if (dir == 2)
                        Bounds_Read[1] = XYZBounds[5];
                }
                else if (Units == "Cells") {
                    if (dir == 0)
                        Bounds_Read[1] = nx - 1;
                    else if (dir == 1)
                        Bounds_Read[1] = ny - 1;
                    else if (dir == 2)
                        Bounds_Read[1] = nz - 1;
                }
                else
                    throw std::runtime_error("Error: Invalid units in ConvertBounds, should be Meters or Cells");
            }
            if (RegionData.contains(onebound)) {
                Bounds_Read[0] = RegionData[onebound];
                Bounds_Read[1] = Bounds_Read[0];
            }
            else if (RegionData.contains(twobounds)) {
                Bounds_Read[0] = RegionData[twobounds][0];
                Bounds_Read[1] = RegionData[twobounds][1];
            }
            // Obtain bounds in both cell units and coordinates
            std::vector<int> Bounds_Cells(2);
            std::vector<double> Bounds_Meters(2);
            if (Units == "Cells") {
                // Bounds provided were in cells, convert to microns
                Bounds_Cells[0] = static_cast<int>(Bounds_Read[0]);
                Bounds_Cells[1] = static_cast<int>(Bounds_Read[1]);
                Bounds_Meters[0] = CoordinateMinVal + deltax * Bounds_Read[0];
                Bounds_Meters[1] = CoordinateMinVal + deltax * Bounds_Read[1];
            }
            else if (Units == "Meters") {
                // Bounds provided were in microns, convert to cells
                Bounds_Meters[0] = Bounds_Read[0];
                Bounds_Meters[1] = Bounds_Read[1];
                Bounds_Cells[0] = std::round((Bounds_Read[0] - CoordinateMinVal) / deltax);
                Bounds_Cells[1] = std::round((Bounds_Read[1] - CoordinateMinVal) / deltax);
            }
            else
                throw std::runtime_error(
                    "Error: Units for a given region in the analysis input file should be Cells or Meters");
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
    void GetRegionSize(double deltax) {
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
};

#endif

#endif
