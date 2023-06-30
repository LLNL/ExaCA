// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
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

#include <nlohmann/json.hpp>
// Structure that holds the details of a representative region of the larger microstructure, along which the specific
// analyses to be performed
struct RepresentativeRegion {

    std::string regionName;
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
    std::vector<float> GrainExtentX, GrainExtentY, GrainExtentZ;

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
    std::vector<bool> AnalysisOptions_PerGrainStatsYN = std::vector<bool>(7, false);
    std::vector<std::string> AnalysisOptions_LayerwiseStats_key = {
        "MeanGrainArea",        // volume only
        "MeanWeightedGrainArea" // volume only
    };
    bool PrintPerGrainStatsYN;
    std::vector<bool> AnalysisOptions_LayerwiseStatsYN = std::vector<bool>(2, false);

    // Analysis options that print separate files
    bool PrintExaConstitYN, PrintPoleFigureYN, PrintInversePoleFigureMapYN;

    // List of grain ID values in the representative region
    std::vector<int> GrainIDVector;
    // List of unique grain IDs associated with the region and the number of grains
    int NumberOfGrains;
    std::vector<int> UniqueGrainIDVector;
    // Size (in units of length, area, or volume, depending on regionType) associated with each grain
    std::vector<float> GrainSizeVector_Microns;

    // Constructor
    RepresentativeRegion(nlohmann::json AnalysisData, std::string RegionName, int nx, int ny, int nz, double deltax,
                         std::vector<double> XYZBounds, ViewI3D_H GrainID) {

        // Data for the specific region of interest
        regionName = RegionName;
        std::cout << "Parsing data for region " << RegionName << std::endl;
        nlohmann::json RegionData = AnalysisData["Regions"][RegionName];
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

        // List of grain ID values in the representative region
        GrainIDVector = getGrainIDVector(GrainID);
        // List of unique grain IDs associated with the region, also initialize the number of grains
        UniqueGrainIDVector = getUniqueGrains();
        // Size (in units of length, area, or volume, depending on regionType) associated with each grain
        GrainSizeVector_Microns = getGrainSizeVector(deltax);
        std::cout << "Loaded analysis options for region " << RegionName << std::endl;
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
            if (RegionData.contains("printPoleFigureData")) {
                PrintPoleFigureYN = RegionData["printPoleFigureData"];
            }
            else
                PrintPoleFigureYN = false;
        }
        else
            PrintPoleFigureYN = false;

        // IPF map data for areas only
        if (regionType == "area") {
            if (RegionData.contains("printInversePoleFigureData")) {
                PrintInversePoleFigureMapYN = RegionData["printInversePoleFigureData"];
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
        if (regionType == "length")
            units_dimension = "microns";
        else if (regionType == "area")
            units_dimension = "square microns";
        else if (regionType == "volume")
            units_dimension = "cubic microns";
        else
            throw std::runtime_error("Error: Invalid region type in setRegionSize during region construction");
    }

    // Helper functions used to analyze a region of microstructure data
    // Subroutines starting with "get" return the data specified
    // Subroutines starting with "calc" calculate the quantity specified but do not return it

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

    // Given an input vector of integer Grain ID values, return an output vector consisting of the unique Grain ID
    // values, sorted from lowest to highest. Store the number of grains
    std::vector<int> getUniqueGrains() {
        std::vector<int> UniqueGrainIDVector = GrainIDVector;
        std::sort(UniqueGrainIDVector.begin(), UniqueGrainIDVector.end());
        std::vector<int>::iterator it;
        it = std::unique(UniqueGrainIDVector.begin(), UniqueGrainIDVector.end());
        UniqueGrainIDVector.resize(std::distance(UniqueGrainIDVector.begin(), it));
        NumberOfGrains = UniqueGrainIDVector.size();
        return UniqueGrainIDVector;
    }

    // Given an input vector of integer Grain ID values "GrainIDVector", and an input vector of the
    // unique Grain ID values "UniqueGrainIDVector" (of size "NumberOfGrains"), return a third vector "GrainSizeVector"
    // listing the size of each of the "NumberOfGrains" grains (scaled by the cell size deltax and depending on whether
    // the region is 1D, 2D, or 3D)
    std::vector<float> getGrainSizeVector(double deltax) {

        std::vector<float> GrainSizeVector(NumberOfGrains);
        double conv = convertToMicrons(deltax, regionType);
        for (int n = 0; n < NumberOfGrains; n++) {
            int GrainSizeCells = std::count(GrainIDVector.begin(), GrainIDVector.end(), UniqueGrainIDVector[n]);
            // convert to either microns, square microns, or cubic microns
            GrainSizeVector[n] = conv * GrainSizeCells;
        }
        return GrainSizeVector;
    }

    // Given the 3D grain structure "GrainID", determine the extent in the direction specified of each of the
    // "NumberOfGrains" unique grains from the volume bounded by [XLow, XHigh], [YLow, YHigh], [ZLow, ZHigh]. Extent is
    // calculated in microns
    void calcGrainExtent(std::vector<float> &GrainExtent, ViewI3D_H GrainID, std::string Direction, double deltax) {

        for (int n = 0; n < NumberOfGrains; n++) {
            int ThisGrainID = UniqueGrainIDVector[n];
            int ThisGrainSize = std::round(GrainSizeVector_Microns[n] * convertToCells(deltax, regionType));
            std::vector<int> GrainCoordinate(ThisGrainSize);
            int count = 0;
            for (int k = zBounds_Cells[0]; k <= zBounds_Cells[1]; k++) {
                for (int i = xBounds_Cells[0]; i <= xBounds_Cells[1]; i++) {
                    for (int j = yBounds_Cells[0]; j <= yBounds_Cells[1]; j++) {
                        if (GrainID(k, i, j) == ThisGrainID) {
                            if (Direction == "X")
                                GrainCoordinate[count] = i;
                            else if (Direction == "Y")
                                GrainCoordinate[count] = j;
                            else if (Direction == "Z")
                                GrainCoordinate[count] = k;
                            count++;
                        }
                    }
                }
            }
            int MinCoord = *std::min_element(GrainCoordinate.begin(), GrainCoordinate.end());
            int MaxCoord = *std::max_element(GrainCoordinate.begin(), GrainCoordinate.end());
            GrainExtent[n] = (MaxCoord - MinCoord + 1) * convertToMicrons(deltax, "length");
        }
    }

    void calcNecessaryGrainExtents(ViewI3D_H GrainID, double deltax) {
        // For the analysis options:
        // If StatsYN[3] or PerGrainStatsYN[5] is toggled, all grain extents are needed. If StatsYN[4] or
        // PerGrainStatsYN[2] is toggled, X extents are needed. If StatsYN[5] or PerGrainStatsYN[3] is toggled, Y
        // extents are needed. If StatsYN[6] or PerGrainStatsYN[4] is toggled, Z extents are needed
        bool calcExtentX = false;
        bool calcExtentY = false;
        bool calcExtentZ = false;
        if ((AnalysisOptions_StatsYN[3]) || (AnalysisOptions_StatsYN[4])) {
            calcExtentX = true;
            calcExtentY = true;
            calcExtentZ = true;
        }
        else {
            if ((AnalysisOptions_StatsYN[4]) || (AnalysisOptions_PerGrainStatsYN[2]))
                calcExtentX = true;
            if ((AnalysisOptions_StatsYN[5]) || (AnalysisOptions_PerGrainStatsYN[3]))
                calcExtentY = true;
            if ((AnalysisOptions_StatsYN[6]) || (AnalysisOptions_PerGrainStatsYN[4]))
                calcExtentZ = true;
        }
        if (calcExtentX) {
            GrainExtentX.resize(NumberOfGrains);
            calcGrainExtent(GrainExtentX, GrainID, "X", deltax);
        }
        if (calcExtentY) {
            GrainExtentY.resize(NumberOfGrains);
            calcGrainExtent(GrainExtentY, GrainID, "Y", deltax);
        }
        if (calcExtentZ) {
            GrainExtentZ.resize(NumberOfGrains);
            calcGrainExtent(GrainExtentZ, GrainID, "Z", deltax);
        }
    }

    // Calculate misorientation relative to the specified cardinal direction for each grain and store in the vector
    std::vector<float> getGrainMisorientation(std::string Direction, ViewF_H GrainUnitVector,
                                              int NumberOfOrientations) {

        std::vector<float> GrainMisorientationVector(NumberOfGrains);
        int direction;
        if (Direction == "X")
            direction = 0;
        else if (Direction == "Y")
            direction = 1;
        else if (Direction == "Z")
            direction = 2;
        else
            throw std::runtime_error(
                "Error: invalid direction specified in calcGrainMisorientation: should be X, Y, or Z");
        ViewF_H GrainMisorientation = MisorientationCalc(NumberOfOrientations, GrainUnitVector, direction);
        for (int n = 0; n < NumberOfGrains; n++) {
            int MyOrientation = getGrainOrientation(UniqueGrainIDVector[n], NumberOfOrientations);
            float MyMisorientation = GrainMisorientation(MyOrientation);
            GrainMisorientationVector[n] = MyMisorientation;
        }
        return GrainMisorientationVector;
    }

    // Create a histogram of orientations for texture determination, using the GrainID values in the volume bounded by
    // [XMin,XMax], [YMin,YMax], [ZMin,ZMax] and excluding and cells that did not undergo melting (GrainID = -1)
    ViewI_H getOrientationHistogram(int NumberOfOrientations, ViewI3D_H GrainID, ViewI3D_H LayerID) {

        // Init histogram values to zero
        ViewI_H GOHistogram("GOHistogram", NumberOfOrientations);
        for (int k = zBounds_Cells[0]; k <= zBounds_Cells[1]; k++) {
            for (int j = yBounds_Cells[0]; j <= yBounds_Cells[1]; j++) {
                for (int i = xBounds_Cells[0]; i <= xBounds_Cells[1]; i++) {
                    if (LayerID(k, i, j) != -1) {
                        int GOVal = getGrainOrientation(GrainID(k, i, j), NumberOfOrientations);
                        GOHistogram(GOVal)++;
                    }
                }
            }
        }
        return GOHistogram;
    }

    //*****************************************************************************/
    // "print" routines print average quantities to the screen (and to the file stream specified by QoIs)
    // "write" routines write data to specific files
    //*****************************************************************************/
    // Print information about a representative area or volume to the console/QoIs file
    void printAnalysisHeader(std::ofstream &QoIs) {

        std::string temp;
        temp = "Stats for " + regionType + " " + regionName + "\n";
        if (regionOrientation == "XY") {
            temp += "The representative area is located at Z = " + std::to_string(zBounds_Meters[0]) + " m\n";
            temp += "(in CA units, Z = " + std::to_string(zBounds_Cells[0]) + ")\n";
            temp += "The representative area is bounded by the region spanning X = [" +
                    std::to_string(xBounds_Meters[0]) + "," + std::to_string(xBounds_Meters[1]) + "], Y = [ " +
                    std::to_string(yBounds_Meters[0]) + "," + std::to_string(yBounds_Meters[1]) + "] m\n";
            temp += "(in CA units X = [" + std::to_string(xBounds_Cells[0]) + "," + std::to_string(xBounds_Cells[1]) +
                    "], Y = [" + std::to_string(yBounds_Cells[0]) + "," + std::to_string(yBounds_Cells[1]) + "])\n";
        }
        else if (regionOrientation == "XZ") {
            temp += "The representative area is located at Y = " + std::to_string(yBounds_Meters[0]) + " m\n";
            temp += "(in CA units, Y = " + std::to_string(yBounds_Cells[0]) + ")\n";
            temp += "The representative area is bounded by the region spanning X = [" +
                    std::to_string(xBounds_Meters[0]) + "," + std::to_string(xBounds_Meters[1]) + "], Z = [" +
                    std::to_string(zBounds_Meters[0]) + "," + std::to_string(zBounds_Meters[1]) + "] m\n";
            temp += "(in CA units X = [" + std::to_string(xBounds_Cells[0]) + "," + std::to_string(xBounds_Cells[1]) +
                    "], Z = [" + std::to_string(zBounds_Cells[0]) + "," + std::to_string(zBounds_Cells[1]) + "])\n";
        }
        else if (regionOrientation == "YZ") {
            temp += "The representative area is located at X = " + std::to_string(xBounds_Meters[0]) + " m\n";
            temp += "(in CA units, X = " + std::to_string(xBounds_Cells[0]) + ")\n";
            temp += "The representative area is bounded by the region spanning Y = [" +
                    std::to_string(yBounds_Meters[0]) + "," + std::to_string(yBounds_Meters[1]) + "], Z = [" +
                    std::to_string(zBounds_Meters[0]) + "," + std::to_string(zBounds_Meters[1]) + "] m\n";
            temp += "(in CA units Y = [" + std::to_string(yBounds_Cells[0]) + "," + std::to_string(yBounds_Cells[1]) +
                    "], Z = [" + std::to_string(zBounds_Cells[0]) + "," + std::to_string(zBounds_Cells[1]) + "])\n";
        }
        else if (regionOrientation == "XYZ") {
            temp += "The representative volume specified is bounded by X = [" + std::to_string(xBounds_Meters[0]) +
                    "," + std::to_string(xBounds_Meters[1]) + "], Y = [" + std::to_string(yBounds_Meters[0]) + "," +
                    std::to_string(yBounds_Meters[1]) + "], and Z = [" + std::to_string(zBounds_Meters[0]) + "," +
                    std::to_string(zBounds_Meters[1]) + "] m\n";
            temp += "The representative volume specified is bounded by cells spanning X = [" +
                    std::to_string(xBounds_Cells[0]) + "," + std::to_string(xBounds_Cells[1]) + "], Y = [" +
                    std::to_string(yBounds_Cells[0]) + "," + std::to_string(yBounds_Cells[1]) + "], and Z = [" +
                    std::to_string(zBounds_Cells[0]) + "," + std::to_string(zBounds_Cells[1]) + "]\n";
        }
        else
            throw std::runtime_error("Error: invalid region type");
        dual_print(temp, std::cout, QoIs);
    }

    // Print number of cells in the representative region that did not undergo melting, fraction consisting of nucleated
    // grains to the console/QoIs file
    void printGrainTypeFractions(std::ofstream &QoIs, ViewI3D_H GrainID, ViewI3D_H LayerID) {

        int NumberOfUnmeltedCells = 0;
        int NumberOfNucleatedGrainCells = 0;
        for (int k = zBounds_Cells[0]; k <= zBounds_Cells[1]; k++) {
            for (int i = xBounds_Cells[0]; i <= xBounds_Cells[1]; i++) {
                for (int j = yBounds_Cells[0]; j <= yBounds_Cells[1]; j++) {
                    if (LayerID(k, i, j) == -1)
                        NumberOfUnmeltedCells++;
                    if (GrainID(k, i, j) < 0)
                        NumberOfNucleatedGrainCells++;
                }
            }
        }
        std::string temp;
        temp = "-- The representative region consists of " + std::to_string(regionSize_Cells) + " cells\n";
        temp += "-- The number of cells in the region that did not undergo melting is " +
                std::to_string(NumberOfUnmeltedCells) + "\n";
        float VolFractNucGrains = DivideCast<float>(NumberOfNucleatedGrainCells, regionSize_Cells);
        temp += "-- The volume fraction consisting of nucleated grains is " + std::to_string(VolFractNucGrains) + "\n";
        dual_print(temp, std::cout, QoIs);
    }

    // Print average misorientation for the region relative to each direction
    void printMeanMisorientations(std::ofstream &QoIs, std::vector<float> GrainMisorientationXVector,
                                  std::vector<float> GrainMisorientationYVector,
                                  std::vector<float> GrainMisorientationZVector) {

        std::vector<std::string> MisorientationDirectionLabels = {"misorientationX", "misorientationY",
                                                                  "misorientationZ"};
        std::vector<std::string> MisorientationDirectionLabelsShort = {"+X", "+Y", "+Z"};
        float GrainMisorientation_sum_x = 0.0;
        float GrainMisorientation_sum_y = 0.0;
        float GrainMisorientation_sum_z = 0.0;
        for (int n = 0; n < NumberOfGrains; n++) {
            GrainMisorientation_sum_x += GrainMisorientationXVector[n] * GrainSizeVector_Microns[n];
            GrainMisorientation_sum_y += GrainMisorientationYVector[n] * GrainSizeVector_Microns[n];
            GrainMisorientation_sum_z += GrainMisorientationZVector[n] * GrainSizeVector_Microns[n];
        }
        float AvgMisorientationX = DivideCast<float>(GrainMisorientation_sum_x, regionSize_Meters);
        float AvgMisorientationY = DivideCast<float>(GrainMisorientation_sum_y, regionSize_Meters);
        float AvgMisorientationZ = DivideCast<float>(GrainMisorientation_sum_z, regionSize_Meters);
        std::string temp;
        temp = "-- Average misorientation (weighted by size) for grains relative to the X direction (in degrees): " +
               std::to_string(AvgMisorientationX) + "\n";
        temp += "-- Average misorientation (weighted by size) for grains relative to the Y direction (in degrees): " +
                std::to_string(AvgMisorientationY) + "\n";
        temp += "-- Average misorientation (weighted by size) for grains relative to the Z direction (in degrees): " +
                std::to_string(AvgMisorientationZ) + "\n";
        dual_print(temp, std::cout, QoIs);
    }

    // Print the average grain size and the number of grains in the region
    void printMeanSize(std::ofstream &QoIs) {

        double AvgSizePerGrain = DivideCast<double>(regionSize_Microns, NumberOfGrains);
        std::string temp = "-- There are " + std::to_string(NumberOfGrains) + " grains in this " + regionType +
                           " , and the mean grain " + regionType + " is " + std::to_string(AvgSizePerGrain) + " " +
                           units_dimension + "\n";
        dual_print(temp, std::cout, QoIs);
    }

    // Print average grain extent in the specified direction
    void printMeanExtent(std::ofstream &QoIs, std::string Direction) {

        float GrainExtentSum = 0.0;
        std::vector<float> GrainExtent_Microns;
        if (Direction == "X")
            GrainExtent_Microns = GrainExtentX;
        else if (Direction == "Y")
            GrainExtent_Microns = GrainExtentY;
        else if (Direction == "Z")
            GrainExtent_Microns = GrainExtentZ;
        for (int n = 0; n < NumberOfGrains; n++)
            GrainExtentSum += GrainExtent_Microns[n];
        float AvgGrainExtent = DivideCast<float>(GrainExtentSum, NumberOfGrains);
        std::string temp = "-- The mean grain extent in the " + Direction + " direction is " +
                           std::to_string(AvgGrainExtent) + " microns\n";
        dual_print(temp, std::cout, QoIs);
    }

    // Print average aspect ratio in the build to the average of the transverse directions
    void printMeanBuildTransAspectRatio(std::ofstream &QoIs) {

        std::vector<float> GrainAspectRatios(NumberOfGrains);
        float ARSum = 0.0;
        float VolWtARSum = 0.0;
        for (int n = 0; n < NumberOfGrains; n++) {
            float AR_XZ = GrainExtentZ[n] / GrainExtentX[n];
            float AR_YZ = GrainExtentZ[n] / GrainExtentY[n];
            GrainAspectRatios[n] = 0.5 * (AR_XZ + AR_YZ);
            ARSum += GrainAspectRatios[n];
            VolWtARSum += GrainAspectRatios[n] * GrainSizeVector_Microns[n];
        }
        double RepresentativeRegionSize_Microns = regionSize_Meters * pow(10, 18);
        std::string temp;
        temp = "-- The mean grain aspect ratio (Z direction to transverse) is " +
               std::to_string(DivideCast<float>(ARSum, NumberOfGrains)) + "\n";
        temp += "-- The mean volume-weighted grain aspect ratio (Z direction to transverse) is " +
                std::to_string(DivideCast<float>(VolWtARSum, RepresentativeRegionSize_Microns)) + "\n";
        dual_print(temp, std::cout, QoIs);
    }

    // Write unweighted and/or weighted grain areas as a function of build height to file(s)
    // TODO: Remove weighted grain area calculations from the next release
    void writeAreaSeries(std::string BaseFileName, double deltax, ViewI3D_H GrainID) {

        bool PrintUnweightedAreas = AnalysisOptions_LayerwiseStatsYN[0];
        bool PrintWeightedAreas = AnalysisOptions_LayerwiseStatsYN[1];
        std::string FName1 = BaseFileName + "_GrainAreas.csv";
        std::string FName2 = BaseFileName + "_WeightedGrainAreas.csv";
        std::ofstream Grainplot1, Grainplot2;
        if (PrintUnweightedAreas) {
            std::cout << "Printing file " << FName1 << " of grain area values (in square microns) for all Z coordinates"
                      << std::endl;
            Grainplot1.open(FName1);
            Grainplot1 << "Zcoordinate(µm),MeanArea(µm2)" << std::endl;
        }
        if (PrintWeightedAreas) {
            std::cout << "Printing file " << FName2
                      << " of weighted grain area values (in square microns) for every 5th Z coordinate" << std::endl;
            std::cout << "Note: Option to print weighted grain area data will be removed in a future release"
                      << std::endl;
            Grainplot2.open(FName2);
            Grainplot2 << "Zcoordinate(µm),WeightedMeanArea(µm2)" << std::endl;
        }

        int LayerArea = (xBounds_Cells[1] - xBounds_Cells[0] + 1) * (yBounds_Cells[1] - yBounds_Cells[0] + 1);
        for (int k = zBounds_Cells[0]; k <= zBounds_Cells[1]; k++) {
            // All GrainID values in the representative area at Z = k
            std::vector<int> GrainIDVector_Area(regionSize_Cells);
            int count = 0;
            for (int i = xBounds_Cells[0]; i <= xBounds_Cells[1]; i++) {
                for (int j = yBounds_Cells[0]; j <= yBounds_Cells[1]; j++) {
                    GrainIDVector_Area[count] = GrainID(k, i, j);
                    count++;
                }
            }
            // Get the number of grains in the representative area and a list of the unique grain IDs
            std::vector<int> UniqueGrainIDVector_Area = GrainIDVector_Area;
            std::sort(UniqueGrainIDVector_Area.begin(), UniqueGrainIDVector_Area.end());
            std::vector<int>::iterator it;
            it = std::unique(UniqueGrainIDVector_Area.begin(), UniqueGrainIDVector_Area.end());
            UniqueGrainIDVector_Area.resize(std::distance(UniqueGrainIDVector_Area.begin(), it));
            int NumberOfGrains_Area = UniqueGrainIDVector_Area.size();
            double MeanGrainAreaThisLayer = DivideCast<double>(LayerArea, NumberOfGrains_Area);
            if (PrintUnweightedAreas)
                Grainplot1 << zBounds_Meters[0] * pow(10, 6) + convertToMicrons(deltax, "length") << ","
                           << MeanGrainAreaThisLayer * convertToMicrons(deltax, "area") << std::endl;
            if ((PrintWeightedAreas) && (k % 5 == 0)) {
                std::vector<float> GrainSizeVector_Microns_Area(NumberOfGrains_Area);
                double conv = convertToMicrons(deltax, "area");
                for (int n = 0; n < NumberOfGrains_Area; n++) {
                    int GrainSizeCells =
                        std::count(GrainIDVector_Area.begin(), GrainIDVector_Area.end(), UniqueGrainIDVector_Area[n]);
                    GrainSizeVector_Microns_Area[n] = conv * GrainSizeCells;
                }
                float AreaXArea = 0.0;
                for (int n = 0; n < NumberOfGrains_Area; n++)
                    AreaXArea += GrainSizeVector_Microns_Area[n] * GrainSizeVector_Microns_Area[n];
                double WeightedArea = DivideCast<double>(AreaXArea, LayerArea);
                Grainplot2 << zBounds_Meters[0] * pow(10, 6) + convertToMicrons(deltax, "length") << "," << WeightedArea
                           << std::endl;
                if (k == zBounds_Cells[1])
                    std::cout
                        << "[Note: this will no longer be printed in a future release] The mean weighted grain area "
                           "at the representative region top (Z coordinate = "
                        << zBounds_Cells[1] << ") is " << WeightedArea << " square microns" << std::endl;
            }
            if (k == zBounds_Cells[1])
                std::cout << "[Note: this will no longer be printed in a future release] The mean grain area at the "
                             "representative region top (Z coordinate = "
                          << zBounds_Cells[1] << ") is " << MeanGrainAreaThisLayer << " square microns" << std::endl;
        }
        if (PrintUnweightedAreas)
            Grainplot1.close();
        if (PrintWeightedAreas)
            Grainplot2.close();
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

    // From the grain extents in x, y, and z, calcualte the aspect ratio for each grain in the build to the average of
    // the transverse directions
    void calcBuildTransAspectRatio(std::vector<float> &BuildTransAspectRatio) {

        for (int n = 0; n < NumberOfGrains; n++) {
            float AR_XZ = GrainExtentZ[n] / GrainExtentX[n];
            float AR_YZ = GrainExtentZ[n] / GrainExtentY[n];
            BuildTransAspectRatio[n] = 0.5 * (AR_XZ + AR_YZ);
        }
    }

    // Determine the R, G, or B value for this grain orientation for the IPF-Z inverse pole figure colormap (Color = 0
    // for R, 1 for G, 2 for B)
    std::vector<float> getIPFZColor(int Color, int NumberOfOrientations, ViewF_H GrainRGBValues) {
        std::vector<float> IPFZColor(NumberOfGrains);
        for (int n = 0; n < NumberOfGrains; n++) {
            int MyOrientation = getGrainOrientation(UniqueGrainIDVector[n], NumberOfOrientations);
            IPFZColor[n] = GrainRGBValues(3 * MyOrientation + Color);
        }
        return IPFZColor;
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
        for (int k = ZLow; k <= ZHigh; k++) {
            for (int j = YLow; j <= YHigh; j++) {
                for (int i = XLow; i <= XHigh; i++) {
                    Grainplot << i << "," << j << "," << k << "," << GrainID(k, i, j) << std::endl;
                }
            }
        }
        Grainplot.close();
    }

    // Write a csv file of stats for each grain
    void writePerGrainStats(std::string OutputFileName, std::vector<float> GrainMisorientationXVector,
                            std::vector<float> GrainMisorientationYVector,
                            std::vector<float> GrainMisorientationZVector, std::vector<float> BuildTransAspectRatio,
                            std::vector<float> GrainRed, std::vector<float> GrainGreen, std::vector<float> GrainBlue) {

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
                GrainStats << "," << GrainSizeVector_Microns[n];
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
