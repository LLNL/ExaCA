// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "GAprint.hpp"
#include "GAutils.hpp"

//*****************************************************************************/
// "print" routines print average quantities to the screen (and to the file stream specified by QoIs)
// "write" routines write data to specific files
// Print routines endining in "Old" print quantities that will no longer be computer/printed to files in a future
// release
//*****************************************************************************/
// Write header of unique grain ID values to stats file, and print information about the representative region to the
// console/QoIs file
void printAnalysisHeader(std::ofstream &QoIs, const int XLow, const int XHigh, const int YLow, const int YHigh,
                         const int ZLow, const int ZHigh, std::vector<double> XYZBounds) {

    // TODO: Remove redundant code by making a print routine that prints a string to both std::cout and the QoIs file,
    // since the same information is printed to both
    std::cout << "The representative volume specified is bounded by X = [" << XYZBounds[0] << "," << XYZBounds[1]
              << "], Y = [" << XYZBounds[2] << "," << XYZBounds[3] << "], and Z = [" << XYZBounds[4] << ","
              << XYZBounds[5] << "] m" << std::endl;
    std::cout << "The representative volume specified is bounded by cells spanning X = [" << XLow << "," << XHigh
              << "], Y = [ " << YLow << "," << YHigh << "], and Z = [" << ZLow << "," << ZHigh << "] m" << std::endl;
    QoIs << "The representative volume specified is bounded by X = [" << XYZBounds[0] << "," << XYZBounds[1]
         << "], Y = [ " << XYZBounds[2] << "," << XYZBounds[3] << "], and Z = [" << XYZBounds[4] << "," << XYZBounds[5]
         << "] m" << std::endl;
    QoIs << "The representative volume specified is bounded by cells spanning X = [" << XLow << "," << XHigh
         << "], Y = [ " << YLow << "," << YHigh << "], and Z = [" << ZLow << "," << ZHigh << "] m" << std::endl;
}

// Print information about a representative area to the console/QoIs file
void printAnalysisHeader(std::ofstream &QoIs, const int XLow_cells, const int XHigh_cells, const int YLow_cells,
                         const int YHigh_cells, const int ZLow_cells, const int ZHigh_cells, const double XLow_microns,
                         const double XHigh_microns, const double YLow_microns, const double YHigh_microns,
                         const double ZLow_microns, const double ZHigh_microns, const std::string regionName,
                         const std::string regionOrientation) {

    std::cout << "Stats for " << regionName << " plane:" << std::endl;
    QoIs << "Stats for " << regionName << " plane:" << std::endl;
    if (regionOrientation == "XY") {
        std::cout << "The representative area is located at Z = " << ZLow_microns << " microns" << std::endl;
        std::cout << "(in CA units, Z = " << ZLow_cells << ")" << std::endl;
        std::cout << "The representative area is bounded by the region spanning X = [" << XLow_microns << ","
                  << XHigh_microns << "], Y = [ " << YLow_microns << "," << YHigh_microns << "] microns" << std::endl;
        std::cout << "(in CA units X = [" << XLow_cells << "," << XHigh_cells << "], Y = [ " << YLow_cells << ","
                  << YHigh_cells << "])" << std::endl;
        QoIs << "The representative area is located at Z = " << ZLow_microns << " microns" << std::endl;
        QoIs << "(in CA units, Z = " << ZLow_cells << ")" << std::endl;
        QoIs << "The representative area is bounded by the region spanning X = [" << XLow_microns << ","
             << XHigh_microns << "], Y = [ " << YLow_microns << "," << YHigh_microns << "] microns" << std::endl;
        QoIs << "(in CA units X = [" << XLow_cells << "," << XHigh_cells << "], Y = [ " << YLow_cells << ","
             << YHigh_cells << "])" << std::endl;
    }
    else if (regionOrientation == "XZ") {
        std::cout << "The representative area is located at Y = " << YLow_microns << " microns" << std::endl;
        std::cout << "(in CA units, Y = " << YLow_cells << ")" << std::endl;
        std::cout << "The representative area is bounded by the region spanning X = [" << XLow_microns << ","
                  << XHigh_microns << "], Z = [ " << ZLow_microns << "," << ZHigh_microns << "] microns" << std::endl;
        std::cout << "(in CA units X = [" << XLow_cells << "," << XHigh_cells << "], Z = [ " << ZLow_cells << ","
                  << ZHigh_cells << "])" << std::endl;
        QoIs << "The representative area is located at Y = " << YLow_microns << " microns" << std::endl;
        QoIs << "(in CA units, Y = " << YLow_cells << ")" << std::endl;
        QoIs << "The representative area is bounded by the region spanning X = [" << XLow_microns << ","
             << XHigh_microns << "], Z = [ " << ZLow_microns << "," << ZHigh_microns << "] microns" << std::endl;
        QoIs << "(in CA units X = [" << XLow_cells << "," << XHigh_cells << "], Z = [ " << ZLow_cells << ","
             << ZHigh_cells << "])" << std::endl;
    }
    else if (regionOrientation == "YZ") {
        std::cout << "The representative area is located at X = " << XLow_microns << " microns" << std::endl;
        std::cout << "(in CA units, X = " << XLow_cells << ")" << std::endl;
        std::cout << "The representative area is bounded by the region spanning Y = [" << YLow_microns << ","
                  << YHigh_microns << "], Z = [ " << ZLow_microns << "," << ZHigh_microns << "] microns" << std::endl;
        std::cout << "(in CA units Y = [" << YLow_cells << "," << YHigh_cells << "], Z = [ " << ZLow_cells << ","
                  << ZHigh_cells << "])" << std::endl;
        QoIs << "The representative area is located at X = " << XLow_microns << " microns" << std::endl;
        QoIs << "(in CA units, X = " << XLow_cells << ")" << std::endl;
        QoIs << "The representative area is bounded by the region spanning Y = [" << YLow_microns << ","
             << YHigh_microns << "], Z = [ " << ZLow_microns << "," << ZHigh_microns << "] microns" << std::endl;
        QoIs << "(in CA units Y = [" << YLow_cells << "," << YHigh_cells << "], Z = [ " << ZLow_cells << ","
             << ZHigh_cells << "])" << std::endl;
    }
    else
        throw std::runtime_error("Error: invalid region type");
}

// Print number of cells in the representative region that did not undergo melting, fraction consisting of nucleated
// grains to the console/QoIs file
void printGrainTypeFractions(std::ofstream &QoIs, const int XLow, const int XHigh, const int YLow, const int YHigh,
                             const int ZLow, const int ZHigh, ViewI3D_H GrainID, ViewI3D_H LayerID,
                             int RepresentativeRegionSize) {

    int NumberOfUnmeltedCells = 0;
    int NumberOfNucleatedGrainCells = 0;
    for (int k = ZLow; k <= ZHigh; k++) {
        for (int i = XLow; i <= XHigh; i++) {
            for (int j = YLow; j <= YHigh; j++) {
                if (LayerID(k, i, j) == -1)
                    NumberOfUnmeltedCells++;
                if (GrainID(k, i, j) < 0)
                    NumberOfNucleatedGrainCells++;
            }
        }
    }
    std::cout << "-- The representative region consists of " << RepresentativeRegionSize << " cells" << std::endl;
    QoIs << "-- The representative region consists of " << RepresentativeRegionSize << " cells" << std::endl;
    std::cout << "-- The number of cells in the region that did not undergo melting is " << NumberOfUnmeltedCells
              << std::endl;
    QoIs << "-- The number of cells in the region that did not undergo melting is " << NumberOfUnmeltedCells
         << std::endl;
    float VolFractNucGrains = DivideCast<float>(NumberOfNucleatedGrainCells, RepresentativeRegionSize);
    std::cout << "-- The volume fraction consisting of nucleated grains is " << VolFractNucGrains << std::endl;
    QoIs << "-- The volume fraction consisting of nucleated grains is " << VolFractNucGrains << std::endl;
}

// Print average misorientation for the region relative to each direction
void printMeanMisorientations(std::ofstream &QoIs, int NumberOfGrains, std::vector<float> GrainMisorientationXVector,
                              std::vector<float> GrainMisorientationYVector,
                              std::vector<float> GrainMisorientationZVector, std::vector<float> GrainSizeVector,
                              double RepresentativeRegionSize_Microns) {

    std::vector<std::string> MisorientationDirectionLabels = {"misorientationX", "misorientationY", "misorientationZ"};
    std::vector<std::string> MisorientationDirectionLabelsShort = {"+X", "+Y", "+Z"};
    float GrainMisorientation_sum_x = 0.0;
    float GrainMisorientation_sum_y = 0.0;
    float GrainMisorientation_sum_z = 0.0;
    for (int n = 0; n < NumberOfGrains; n++) {
        GrainMisorientation_sum_x += GrainMisorientationXVector[n] * GrainSizeVector[n];
        GrainMisorientation_sum_y += GrainMisorientationYVector[n] * GrainSizeVector[n];
        GrainMisorientation_sum_z += GrainMisorientationZVector[n] * GrainSizeVector[n];
    }
    float AvgMisorientationX = DivideCast<float>(GrainMisorientation_sum_x, RepresentativeRegionSize_Microns);
    float AvgMisorientationY = DivideCast<float>(GrainMisorientation_sum_y, RepresentativeRegionSize_Microns);
    float AvgMisorientationZ = DivideCast<float>(GrainMisorientation_sum_z, RepresentativeRegionSize_Microns);
    std::cout << "-- Average misorientation (weighted by size) for grains relative to the X direction (in degrees): "
              << AvgMisorientationX << std::endl;
    QoIs << "-- Average misorientation (weighted by size) relative to the X direction (in degrees): "
         << AvgMisorientationX << std::endl;
    std::cout << "-- Average misorientation (weighted by size) for grains relative to the Y direction (in degrees): "
              << AvgMisorientationY << std::endl;
    QoIs << "-- Average misorientation (weighted by size) relative to the Y direction (in degrees): "
         << AvgMisorientationY << std::endl;
    std::cout << "-- Average misorientation (weighted by size) for grains relative to the Z direction (in degrees): "
              << AvgMisorientationZ << std::endl;
    QoIs << "-- Average misorientation (weighted by size) relative to the Z direction (in degrees): "
         << AvgMisorientationZ << std::endl;
}

// Print mean misorientation of portion of the volume that underwent melting, and the mean misorientation at the volume
// top to the screen
void printMisorientationDataOld(int XMin, int XMax, int YMin, int YMax, int ZMin, int ZMax, ViewI3D_H LayerID,
                                ViewF_H GrainUnitVector, ViewI3D_H GrainID, int NumberOfOrientations) {

    int NumberOfMeltedCells = 0;
    long double MisorientationSum = 0.0;
    int NumberOfMeltedCellsTop = 0;
    long double MisorientationSumTop = 0.0;
    ViewF_H GrainMisorientation = MisorientationCalc(NumberOfOrientations, GrainUnitVector, 2);
    for (int k = ZMin; k <= ZMax; k++) {
        for (int i = XMin; i <= XMax; i++) {
            for (int j = YMin; j <= YMax; j++) {
                // Only take data from cells in the representative area that underwent melting (LayerID >= 0)
                if (LayerID(k, i, j) != -1) {
                    int MyOrientation = getGrainOrientation(GrainID(k, i, j), NumberOfOrientations);
                    float MyMisorientation = GrainMisorientation(MyOrientation);
                    MisorientationSum += MyMisorientation;
                    if (k == ZMax) {
                        MisorientationSumTop += MyMisorientation;
                        NumberOfMeltedCellsTop++;
                    }
                    NumberOfMeltedCells++;
                }
            }
        }
    }
    std::cout << "[Note: this will no longer be printed in a future release] The mean misorientation relative to the "
                 "+Z direction, excluding regions that didn't undergo "
                 "melting, is "
              << DivideCast<float>(MisorientationSum, NumberOfMeltedCells) << " degrees" << std::endl;
    std::cout << "[Note: this will no longer be printed in a future release] The mean misorientation relative to the "
                 "+Z direction at the representative region top (Z = "
              << ZMax << "), excluding regions that didn't undergo melting, is "
              << MisorientationSumTop / ((double)(NumberOfMeltedCellsTop)) << " degrees" << std::endl;
}

// Print the average grain size and the number of grains in the region
void printMeanSize(std::ofstream &QoIs, int NumberOfGrains, double RepresentativeRegionSize_Microns,
                   std::string RegionType) {

    std::string Units;
    if (RegionType == "length")
        Units = "microns";
    else if (RegionType == "area")
        Units = "square microns";
    else if (RegionType == "volume")
        Units = "cubic microns";
    else
        throw std::runtime_error("Error: unknown region type in printMeanSize");
    double AvgSizePerGrain = DivideCast<double>(RepresentativeRegionSize_Microns, NumberOfGrains);
    std::cout << "-- There are " << NumberOfGrains << " grains in this " << RegionType << " , and the mean grain "
              << RegionType << " is " << AvgSizePerGrain << " " << Units << std::endl;
    QoIs << "-- There are " << NumberOfGrains << " grains in this " << RegionType << " , and the mean grain "
         << RegionType << " is " << AvgSizePerGrain << " " << Units << std::endl;
}

// Print average aspect ratio in the build to the average of the transverse directions
void printMeanBuildTransAspectRatio(std::ofstream &QoIs, std::vector<float> GrainExtentX,
                                    std::vector<float> GrainExtentY, std::vector<float> GrainExtentZ,
                                    std::vector<float> GrainSizeVector, double RepresentativeRegionSize_Microns,
                                    int NumberOfGrains) {

    std::vector<float> GrainAspectRatios(NumberOfGrains);
    float ARSum = 0.0;
    float VolWtARSum = 0.0;
    for (int n = 0; n < NumberOfGrains; n++) {
        float AR_XZ = GrainExtentZ[n] / GrainExtentX[n];
        float AR_YZ = GrainExtentZ[n] / GrainExtentY[n];
        GrainAspectRatios[n] = 0.5 * (AR_XZ + AR_YZ);
        ARSum += GrainAspectRatios[n];
        VolWtARSum += GrainAspectRatios[n] * GrainSizeVector[n];
    }
    std::cout << "-- The mean grain aspect ratio (Z direction to transverse) is "
              << DivideCast<float>(ARSum, NumberOfGrains) << std::endl;
    QoIs << "-- The mean grain aspect ratio (Z direction to transverse) is " << DivideCast<float>(ARSum, NumberOfGrains)
         << std::endl;
    std::cout << "-- The mean volume-weighted grain aspect ratio (Z direction to transverse) is "
              << DivideCast<float>(VolWtARSum, RepresentativeRegionSize_Microns) << std::endl;
    QoIs << "-- The mean volume-weighted grain aspect ratio (Z direction to transverse) is "
         << DivideCast<float>(VolWtARSum, RepresentativeRegionSize_Microns) << std::endl;
}

// Print average grain extent in the specified direction
void printMeanExtent(std::ofstream &QoIs, std::vector<float> GrainExtent, std::string Direction, int NumberOfGrains) {

    float GrainExtentSum = 0.0;
    for (int n = 0; n < NumberOfGrains; n++)
        GrainExtentSum += GrainExtent[n];
    float AvgGrainExtent = DivideCast<float>(GrainExtentSum, NumberOfGrains);
    std::cout << "-- The mean grain extent in the " << Direction << " direction is " << AvgGrainExtent << " microns"
              << std::endl;
    QoIs << "-- The mean grain extent in the " << Direction << " direction is " << AvgGrainExtent << " microns"
         << std::endl;
}

// Print mean grain width (average of extents in x and y), and print the grain width distribution at the top of the
// representative region
void printSizeOld(std::string BaseFileName, int NumberOfGrains, std::vector<float> GrainExtentX,
                  std::vector<float> GrainExtentY, const int XMin, const int XMax, const int YMin, const int YMax,
                  const int ZMax, double deltax, ViewI3D_H GrainID) {

    float GrainExtentSumXY = 0.0;
    for (int n = 0; n < NumberOfGrains; n++) {
        GrainExtentSumXY += GrainExtentX[n];
        GrainExtentSumXY += GrainExtentY[n];
    }
    float AvgGrainExtentXY = DivideCast<float>(GrainExtentSumXY, 2 * NumberOfGrains);
    std::cout << "[Note: this will no longer be printed in a future release] The mean grain width is "
              << AvgGrainExtentXY << " microns" << std::endl;
    std::cout << "Note: option to print grain width distribution at the top of the representative volume as frequency "
                 "files will be removed in a future release"
              << std::endl;
    std::string FNameX = BaseFileName + "_GrainWidthDistributionX.csv";
    std::string FNameY = BaseFileName + "_GrainWidthDistributionY.csv";
    std::ofstream GrainplotX;
    std::ofstream GrainplotY;
    std::cout << "Printing files " << FNameX << " and " << FNameY
              << " of grain width distributions in x and y (in microns) at Z = " << ZMax << std::endl;
    GrainplotX.open(FNameX);
    GrainplotY.open(FNameY);
    // Get list of unique grain IDs within [XMin, XMax], [YMin, YMax] at Z = ZMax
    int RepresentativeAreaSize = (XMax - XMin + 1) * (YMax - YMin + 1);
    std::vector<int> GrainIDVector_Area =
        getRepresentativeRegionGrainIDs(GrainID, XMin, XMax, YMin, YMax, ZMax, ZMax, RepresentativeAreaSize);
    // Get the number of grains in the representative area and a list of the unique grain IDs
    int NumberOfGrains_Area;
    std::vector<int> UniqueGrainIDVector_Area = getUniqueGrains(GrainIDVector_Area, NumberOfGrains_Area);
    std::cout << "[Note: this will no longer be printed in a future release] Number of grains at the representative "
                 "region top (Z = "
              << ZMax << "): " << NumberOfGrains_Area << std::endl;

    // Get the number of cells associated with each unique grain ID value in this area
    std::vector<float> GrainSizeVector_Area =
        getGrainSizes(GrainIDVector_Area, UniqueGrainIDVector_Area, NumberOfGrains_Area, deltax, "area");
    std::vector<float> GrainExtentX_Area(NumberOfGrains_Area);
    std::vector<float> GrainExtentY_Area(NumberOfGrains_Area);
    calcGrainExtent(GrainExtentX_Area, GrainID, UniqueGrainIDVector_Area, GrainSizeVector_Area, NumberOfGrains_Area,
                    XMin, XMax, YMin, YMax, ZMax, ZMax, "X", deltax, "area");
    calcGrainExtent(GrainExtentY_Area, GrainID, UniqueGrainIDVector_Area, GrainSizeVector_Area, NumberOfGrains_Area,
                    XMin, XMax, YMin, YMax, ZMax, ZMax, "Y", deltax, "area");

    // Print grain widths in x and y for this area
    for (int n = 0; n < NumberOfGrains_Area; n++) {
        GrainplotX << GrainExtentX_Area[n] << std::endl;
        GrainplotY << GrainExtentY_Area[n] << std::endl;
    }
    GrainplotX.close();
    GrainplotY.close();
}

// Write unweighted and/or weighted grain areas as a function of build height to file(s)
void writeAreaSeries(bool PrintWeightedAreas, bool PrintUnweightedAreas, std::string BaseFileName, double deltax,
                     int XMin, int XMax, int YMin, int YMax, int ZMin, int ZMax, ViewI3D_H GrainID,
                     double ZMin_Coordinate) {

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
        Grainplot2.open(FName2);
        Grainplot2 << "Zcoordinate(µm),WeightedMeanArea(µm2)" << std::endl;
    }

    int LayerArea = (XMax - XMin + 1) * (YMax - YMin + 1);
    for (int k = ZMin; k <= ZMax; k++) {
        // All GrainID values in the representative area at Z = k
        std::vector<int> GrainIDVector_Area =
            getRepresentativeRegionGrainIDs(GrainID, XMin, XMax, YMin, YMax, k, k, LayerArea);
        // Get the number of grains in the representative area and a list of the unique grain IDs
        int NumberOfGrains_Area;
        std::vector<int> UniqueGrainIDVector_Area = getUniqueGrains(GrainIDVector_Area, NumberOfGrains_Area);
        double MeanGrainAreaThisLayer = DivideCast<double>(LayerArea, NumberOfGrains_Area);
        if (PrintUnweightedAreas)
            Grainplot1 << ZMin_Coordinate + convertToMicrons(deltax, "length") << ","
                       << MeanGrainAreaThisLayer * convertToMicrons(deltax, "area") << std::endl;
        if ((PrintWeightedAreas) && (k % 5 == 0)) {
            std::vector<float> GrainSizeVector_Area =
                getGrainSizes(GrainIDVector_Area, UniqueGrainIDVector_Area, NumberOfGrains_Area, deltax, "area");
            float AreaXArea = 0.0;
            for (int n = 0; n < NumberOfGrains_Area; n++)
                AreaXArea += GrainSizeVector_Area[n] * GrainSizeVector_Area[n];
            double WeightedArea = DivideCast<double>(AreaXArea, LayerArea);
            Grainplot2 << ZMin_Coordinate + convertToMicrons(deltax, "length") << "," << WeightedArea << std::endl;
            if (k == ZMax)
                std::cout << "[Note: this will no longer be printed in a future release] The mean weighted grain area "
                             "at the representative region top (Z coordinate = "
                          << ZMax << ") is " << WeightedArea << " square microns" << std::endl;
        }
        if (k == ZMax)
            std::cout << "[Note: this will no longer be printed in a future release] The mean grain area at the "
                         "representative region top (Z coordinate = "
                      << ZMax << ") is " << MeanGrainAreaThisLayer << " square microns" << std::endl;
    }
    if (PrintUnweightedAreas)
        Grainplot1.close();
    if (PrintWeightedAreas)
        Grainplot2.close();
}

// Write a csv file of stats for each grain
void writePerGrainStats(std::string OutputFileName, std::string RegionType, std::vector<int> UniqueGrainIDVector,
                        std::vector<float> GrainMisorientationXVector, std::vector<float> GrainMisorientationYVector,
                        std::vector<float> GrainMisorientationZVector, std::vector<float> GrainSizeVector,
                        std::vector<float> GrainExtentX, std::vector<float> GrainExtentY,
                        std::vector<float> GrainExtentZ, std::vector<float> BuildTransAspectRatio, bool *AnalysisTypes,
                        int NumberOfGrains, bool PrintIPFRGB, std::vector<float> GrainRed,
                        std::vector<float> GrainGreen, std::vector<float> GrainBlue) {

    // Which quantities should be printed?
    bool PrintMisorientation = AnalysisTypes[0];
    bool PrintSize = AnalysisTypes[1];
    bool PrintBuildTransAspectRatio = AnalysisTypes[2];
    bool PrintExtentBuild = AnalysisTypes[5];
    bool PrintExtentTrans = AnalysisTypes[6];
    std::string stats_fname = OutputFileName + "_grains.csv";
    std::ofstream GrainStats;
    GrainStats.open(stats_fname);
    GrainStats << "GrainID";
    // Which quantities should be printed for each grain?
    if (PrintMisorientation)
        GrainStats << ",misorientationX,misorientationY,misorientationZ";
    if (PrintSize)
        GrainStats << "," << RegionType;
    if (PrintExtentTrans)
        GrainStats << ",extentX,extentY";
    if (PrintExtentBuild)
        GrainStats << ",extentZ";
    if (PrintBuildTransAspectRatio)
        GrainStats << ",buildtransAR";
    if (PrintIPFRGB)
        GrainStats << ",IPFZ_r,IPFZ_g,IPFZ_b";
    GrainStats << std::endl;

    // Print the specified quantities to the csv file
    for (int n = 0; n < NumberOfGrains; n++) {
        GrainStats << UniqueGrainIDVector[n];
        if (PrintMisorientation)
            GrainStats << "," << GrainMisorientationXVector[n] << "," << GrainMisorientationYVector[n] << ","
                       << GrainMisorientationZVector[n];
        if (PrintSize)
            GrainStats << "," << GrainSizeVector[n];
        if (PrintExtentTrans)
            GrainStats << "," << GrainExtentX[n] << "," << GrainExtentY[n];
        if (PrintExtentBuild)
            GrainStats << "," << GrainExtentZ[n];
        if (PrintBuildTransAspectRatio)
            GrainStats << "," << BuildTransAspectRatio[n];
        if (PrintIPFRGB)
            GrainStats << "," << GrainRed[n] << "," << GrainGreen[n] << "," << GrainBlue[n];
        GrainStats << std::endl;
    }
    GrainStats.close();
}
//*****************************************************************************/
// Helper function for unimodal analysis of the grains in the specified cross-section
void AnalyzeCrossSection_Unimodal(std::ofstream &QoIs, std::string BaseFileName, std::string ThisCrossSectionPlane,
                                  double deltax, int NumberOfGrains, int CrossSectionSize, std::vector<int> GrainAreas,
                                  float MinGrainSize_microns = 7.8125) {

    // Count grain areas large enough to be counted
    int NumGrainsAboveThreshold = 0, AreaGrainsAboveThreshold = 0;
    int MinGrainSize =
        std::round(MinGrainSize_microns / (deltax * deltax * pow(10, 12))); // convert area to value in cells
    for (int n = 0; n < NumberOfGrains; n++) {
        if (GrainAreas[n] >= MinGrainSize) {
            AreaGrainsAboveThreshold += GrainAreas[n];
            NumGrainsAboveThreshold++;
        }
    }

    // What fraction of the cross-sectional area consists of grains large enough to be counted?
    double FractAreaAboveThreshold = DivideCast<double>(AreaGrainsAboveThreshold, CrossSectionSize);
    std::cout << "Area fraction of grains too small to include in statistics ( < " << MinGrainSize
              << " cells in area or less): " << 1 - FractAreaAboveThreshold << std::endl;

    // What's the average area for the grains large enough to be counted?
    double AvgAreaAboveThreshold =
        DivideCast<double>(CrossSectionSize * FractAreaAboveThreshold, NumGrainsAboveThreshold);
    std::cout << "Average grain area (in square microns): " << AvgAreaAboveThreshold * deltax * deltax * pow(10, 12)
              << std::endl;
    QoIs << "Average grain area (in square microns): " << AvgAreaAboveThreshold * deltax * deltax * pow(10, 12)
         << std::endl;

    // Print area of each grain to a file, as a fraction of the total area filled by grains large enough to be
    // considered
    std::string AreasFile = BaseFileName + "_" + ThisCrossSectionPlane + "_Areas.txt";
    std::ofstream Areas;
    Areas.open(AreasFile);
    for (int n = 0; n < NumberOfGrains; n++) {
        if (GrainAreas[n] >= MinGrainSize) {
            double ThisGrainArea = DivideCast<double>(GrainAreas[n], AreaGrainsAboveThreshold);
            Areas << ThisGrainArea << std::endl;
        }
    }
    Areas.close();
}

// Helper function for bimodal analysis of the grains in the specified cross-section
void AnalyzeCrossSection_Bimodal(std::ofstream &QoIs, std::string BaseFileName, std::string ThisCrossSectionPlane,
                                 double deltax, int NumberOfGrains, int CrossSectionSize,
                                 std::vector<int> UniqueGrainIDs, std::vector<int> GrainAreas, int NumberOfOrientations,
                                 ViewF_H GrainUnitVector, ViewF_H GrainRGBValues, float MinGrainSize_microns = 7.8125,
                                 float SmallLargeCutoff_microns = 1500) {

    // Make list of grains < "MinGrainSize" cells ("too small to include in statistics"), larger than
    // SmallLargeCutoff_microns sq microns ("large"), and in between ("small")
    int MinGrainSize =
        std::round(MinGrainSize_microns / (deltax * deltax * pow(10, 12))); // convert area to value in cells
    int SmallLargeCutoff =
        std::round(SmallLargeCutoff_microns / (deltax * deltax * pow(10, 12))); // convert area to value in cells
    int NumSmallGrains = 0, NumLargeGrains = 0;
    int AreaTooSmallGrains = 0, AreaSmallGrains = 0, AreaLargeGrains = 0;
    for (int n = 0; n < NumberOfGrains; n++) {
        if (GrainAreas[n] < MinGrainSize) {
            AreaTooSmallGrains += GrainAreas[n];
        }
        else {
            if (GrainAreas[n] < SmallLargeCutoff) {
                AreaSmallGrains += GrainAreas[n];
                NumSmallGrains++;
            }
            else {
                AreaLargeGrains += GrainAreas[n];
                NumLargeGrains++;
            }
        }
    }

    // What fraction of the cross-sectional area consists of each type of grain?
    double FractAreaTooSmall = DivideCast<double>(AreaTooSmallGrains, CrossSectionSize);
    double FractAreaSmall = DivideCast<double>(AreaSmallGrains, CrossSectionSize);
    double FractAreaLarge = DivideCast<double>(AreaLargeGrains, CrossSectionSize);
    std::cout << "Area fraction of grains too small to include in statistics ( < " << MinGrainSize
              << " cells in area or less): " << FractAreaTooSmall << std::endl;
    std::cout << "Area fraction of grains smaller than 1500 sq microns: " << FractAreaSmall << std::endl;
    QoIs << "Area fraction of grains smaller than 1500 sq microns: " << FractAreaSmall << std::endl;
    std::cout << "Area fraction of grains greater than or equal to 1500 sq microns: " << FractAreaLarge << std::endl;

    // What's the average area for small and large grains?
    double AvgAreaSmall = DivideCast<double>(CrossSectionSize * FractAreaSmall, NumSmallGrains);
    double AvgAreaLarge = DivideCast<double>(CrossSectionSize * FractAreaLarge, NumLargeGrains);
    std::cout << "Average area (in square microns) for small grains: " << AvgAreaSmall * deltax * deltax * pow(10, 12)
              << std::endl;
    std::cout << "Average area (in square microns) for large grains: " << AvgAreaLarge * deltax * deltax * pow(10, 12)
              << std::endl;
    QoIs << "Average area (in square microns) for small grains: " << AvgAreaSmall * deltax * deltax * pow(10, 12)
         << std::endl;
    QoIs << "Average area (in square microns) for large grains: " << AvgAreaLarge * deltax * deltax * pow(10, 12)
         << std::endl;

    // Print area of each small grain to a file, and area of each large grain to a file
    std::string SmallAreasFile = BaseFileName + "_" + ThisCrossSectionPlane + "_SmallAreas.txt";
    std::string LargeAreasFile = BaseFileName + "_" + ThisCrossSectionPlane + "_LargeAreas.txt";
    std::ofstream SmallAreas, LargeAreas;
    SmallAreas.open(SmallAreasFile);
    LargeAreas.open(LargeAreasFile);
    for (int n = 0; n < NumberOfGrains; n++) {
        if (GrainAreas[n] >= MinGrainSize) {
            if (GrainAreas[n] < SmallLargeCutoff)
                SmallAreas << static_cast<double>(GrainAreas[n] * deltax * deltax * pow(10, 12)) << std::endl;
            else
                LargeAreas << static_cast<double>(GrainAreas[n] * deltax * deltax * pow(10, 12)) << std::endl;
        }
    }
    SmallAreas.close();
    LargeAreas.close();

    // Print misorientation relative to the X, Y, and Z directions for small and large grains, weighted by grain area
    // For large grains, calculate the ratio of green (101 orientations) to blue (111 orientations)
    std::string MisorientationFileSmallX =
        BaseFileName + "_" + ThisCrossSectionPlane + "_XMisorientationSmallAreas.txt";
    std::string MisorientationFileSmallY =
        BaseFileName + "_" + ThisCrossSectionPlane + "_YMisorientationSmallAreas.txt";
    std::string MisorientationFileSmallZ =
        BaseFileName + "_" + ThisCrossSectionPlane + "_ZMisorientationSmallAreas.txt";
    std::string MisorientationFileLargeX =
        BaseFileName + "_" + ThisCrossSectionPlane + "_XMisorientationLargeAreas.txt";
    std::string MisorientationFileLargeY =
        BaseFileName + "_" + ThisCrossSectionPlane + "_YMisorientationLargeAreas.txt";
    std::string MisorientationFileLargeZ =
        BaseFileName + "_" + ThisCrossSectionPlane + "_ZMisorientationLargeAreas.txt";
    std::string GreenBlueFreqFile = BaseFileName + "_" + ThisCrossSectionPlane + "_GreenBlueRatio.txt";
    std::ofstream MisorientationSmallX, MisorientationSmallY, MisorientationSmallZ, MisorientationLargeX,
        MisorientationLargeY, MisorientationLargeZ, GreenBlueFreq;
    MisorientationSmallX.open(MisorientationFileSmallX);
    MisorientationSmallY.open(MisorientationFileSmallY);
    MisorientationSmallZ.open(MisorientationFileSmallZ);
    MisorientationLargeX.open(MisorientationFileLargeX);
    MisorientationLargeY.open(MisorientationFileLargeY);
    MisorientationLargeZ.open(MisorientationFileLargeZ);
    GreenBlueFreq.open(GreenBlueFreqFile);

    ViewF_H GrainMisorientationX = MisorientationCalc(NumberOfOrientations, GrainUnitVector, 0);
    ViewF_H GrainMisorientationY = MisorientationCalc(NumberOfOrientations, GrainUnitVector, 1);
    ViewF_H GrainMisorientationZ = MisorientationCalc(NumberOfOrientations, GrainUnitVector, 2);
    float GrainMisorientation_small_sumx = 0, GrainMisorientation_small_sumy = 0, GrainMisorientation_small_sumz = 0;
    float GrainMisorientation_large_sumx = 0, GrainMisorientation_large_sumy = 0, GrainMisorientation_large_sumz = 0;
    float GreenBlueRatioSum = 0.0;
    for (int n = 0; n < NumberOfGrains; n++) {
        if (GrainAreas[n] >= MinGrainSize) {

            int ThisGrainOrientation = getGrainOrientation(UniqueGrainIDs[n], NumberOfOrientations);
            if (ThisGrainOrientation < 0)
                throw std::runtime_error(
                    "Error analyzing grain misorientations: GrainID = 0 grain is present in the XY cross-section");
            float ThisGrainMisorientationX = GrainMisorientationX(ThisGrainOrientation);
            float ThisGrainMisorientationY = GrainMisorientationY(ThisGrainOrientation);
            float ThisGrainMisorientationZ = GrainMisorientationZ(ThisGrainOrientation);

            if (GrainAreas[n] < SmallLargeCutoff) {
                GrainMisorientation_small_sumx += GrainAreas[n] * ThisGrainMisorientationX;
                GrainMisorientation_small_sumy += GrainAreas[n] * ThisGrainMisorientationY;
                GrainMisorientation_small_sumz += GrainAreas[n] * ThisGrainMisorientationZ;
                for (int nn = 0; nn < GrainAreas[n]; nn++) {
                    MisorientationSmallX << ThisGrainMisorientationX << std::endl;
                    MisorientationSmallY << ThisGrainMisorientationY << std::endl;
                    MisorientationSmallZ << ThisGrainMisorientationZ << std::endl;
                }
            }
            else {
                GrainMisorientation_large_sumx += GrainAreas[n] * ThisGrainMisorientationX;
                GrainMisorientation_large_sumy += GrainAreas[n] * ThisGrainMisorientationY;
                GrainMisorientation_large_sumz += GrainAreas[n] * ThisGrainMisorientationZ;
                float Green = GrainRGBValues(3 * ThisGrainOrientation + 1);
                float Blue = GrainRGBValues(3 * ThisGrainOrientation + 2);
                float GreenBlueRatio_ThisGrain = Green / (Green + Blue);
                GreenBlueFreq << GreenBlueRatio_ThisGrain << std::endl;
                GreenBlueRatioSum += GreenBlueRatio_ThisGrain;
                for (int nn = 0; nn < GrainAreas[n]; nn++) {
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
    GreenBlueFreq.close();
    float AvgMisorientation_Small_X = DivideCast<float>(GrainMisorientation_small_sumx, AreaSmallGrains);
    float AvgMisorientation_Small_Y = DivideCast<float>(GrainMisorientation_small_sumy, AreaSmallGrains);
    float AvgMisorientation_Small_Z = DivideCast<float>(GrainMisorientation_small_sumz, AreaSmallGrains);
    float AvgMisorientation_Large_X = DivideCast<float>(GrainMisorientation_large_sumx, AreaLargeGrains);
    float AvgMisorientation_Large_Y = DivideCast<float>(GrainMisorientation_large_sumy, AreaLargeGrains);
    float AvgMisorientation_Large_Z = DivideCast<float>(GrainMisorientation_large_sumz, AreaLargeGrains);
    float AvgGreenBlueRatio = GreenBlueRatioSum / NumLargeGrains;
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
    std::cout << "Average 101:111 ratio for large grains relative to the Z direction: " << AvgGreenBlueRatio
              << std::endl;
    QoIs << "Average misorientation for small grains relative to the X direction: " << AvgMisorientation_Small_X
         << std::endl;
    QoIs << "Average misorientation for small grains relative to the Y direction: " << AvgMisorientation_Small_Y
         << std::endl;
    QoIs << "Average misorientation for small grains relative to the Z direction: " << AvgMisorientation_Small_Z
         << std::endl;
    QoIs << "Average misorientation for large grains relative to the X direction: " << AvgMisorientation_Large_X
         << std::endl;
    QoIs << "Average misorientation for large grains relative to the Y direction: " << AvgMisorientation_Large_Y
         << std::endl;
    QoIs << "Average misorientation for large grains relative to the Z direction: " << AvgMisorientation_Large_Z
         << std::endl;
    QoIs << "Average 101:111 ratio for large grains relative to the Z direction: " << AvgGreenBlueRatio << std::endl;
}

// Analysis of data for the speified cross-section(s)
// FIXME: This subroutine should be removed in the future, as the specific analysis modes on a given volume or
// cross-section are called directly from main
void printCrossSectionData(int NumberOfCrossSections, std::string BaseFileName,
                           std::vector<std::string> CrossSectionPlane, std::vector<int> CrossSectionLocation, int nx,
                           int ny, int nz, int NumberOfOrientations, ViewI3D_H GrainID,
                           std::vector<bool> PrintSectionPF, std::vector<bool> PrintSectionIPF,
                           std::vector<bool> BimodalAnalysis, double deltax, ViewF_H GrainUnitVector,
                           ViewF_H GrainEulerAngles, ViewF_H GrainRGBValues, std::vector<std::string> CSLabels) {

    // Open file of cross-section quantities of interest
    std::ofstream QoIs;
    std::string OutputFilenameQoIs = BaseFileName + "_QoI.txt";
    QoIs.open(OutputFilenameQoIs);

    // Loop over each cross-section specified in the analysis file (grain_analysis) or over the two default
    // cross-sections (grain_analysis_amb)
    for (int n = 0; n < NumberOfCrossSections; n++) {

        // Print cross-section label to file
        QoIs << CSLabels[n] << std::endl;

        // Print data for pyEBSD/MTEX
        std::string ThisCrossSectionPlane = CrossSectionPlane[n]; // Which kind of cross-section is this?
        int Index1Low = 0;
        int Index2Low = 0;
        int Index1High, Index2High; // Values depend on the cross-section axes: nx, ny, or nz
        std::string Plane;
        if (ThisCrossSectionPlane.find("XZ") != std::string::npos) {
            Index1High = nx;
            Index2High = nz;
            Plane = "XZ";
        }
        else if (ThisCrossSectionPlane.find("YZ") != std::string::npos) {
            Index1High = ny;
            Index2High = nz;
            Plane = "YZ";
        }
        else if (ThisCrossSectionPlane.find("XY") != std::string::npos) {
            Index1High = nx;
            Index2High = ny;
            Plane = "XY";
        }
        else
            throw std::runtime_error("Error: cross-section for analysis must be XZ, YZ, or XY");
        int CrossSectionOutOfPlaneLocation = CrossSectionLocation[n];

        std::cout << "Printing cross-section data for cross-section " << ThisCrossSectionPlane << std::endl;
        int NucleatedGrainCells = 0;
        int UnmeltedCells = 0;
        int CrossSectionSize = (Index1High - Index1Low) * (Index2High - Index2Low);
        std::vector<int> CrossSectionGrainIDs(CrossSectionSize);
        int Counter = 0;
        for (int Index1 = Index1Low; Index1 < Index1High; Index1++) {
            for (int Index2 = Index2Low; Index2 < Index2High; Index2++) {
                // How do Index1, Index2, CrossSectionOutOfPlaneLocation correspond to GrainID(Z loc, X loc, Yloc)?
                int ZLoc, XLoc, YLoc;
                if (Plane == "XY") {
                    XLoc = Index1;
                    YLoc = Index2;
                    ZLoc = CrossSectionOutOfPlaneLocation;
                }
                else if (Plane == "YZ") {
                    XLoc = CrossSectionOutOfPlaneLocation;
                    YLoc = Index1;
                    ZLoc = Index2;
                }
                else {
                    XLoc = Index1;
                    YLoc = CrossSectionOutOfPlaneLocation;
                    ZLoc = Index2;
                }
                // Count number of cells in this cross-section have GrainID < 0 (grains formed via nucleation)
                if (GrainID(ZLoc, XLoc, YLoc) < 0)
                    NucleatedGrainCells++;
                // Add this GrainID to the vector of GrainIDs for this cross-section only if it is not equal to 0
                if (GrainID(ZLoc, XLoc, YLoc) == 0)
                    UnmeltedCells++;
                else {
                    CrossSectionGrainIDs[Counter] = GrainID(ZLoc, XLoc, YLoc);
                    Counter++;
                }
            }
        }
        double AreaFractNucleatedGrains = DivideCast<double>(NucleatedGrainCells, CrossSectionSize);
        double AreaFractUnmelted = DivideCast<double>(UnmeltedCells, CrossSectionSize);
        std::cout << "The fraction of the cross-section that went unmelted (not assigned a GrainID) is "
                  << AreaFractUnmelted << std::endl;
        // Resize cross-section to exclude unmelted cells
        CrossSectionSize = Counter;
        CrossSectionGrainIDs.resize(CrossSectionSize);
        QoIs << "The fraction of grains in this cross-section formed via nucleation events is "
             << AreaFractNucleatedGrains << std::endl;
        // Collect grain euler angles for the given plane to write to a file to be read by MTEX/plotted as inverse pole
        // figure-colored cross-sections
        if (PrintSectionIPF[n])
            writeIPFColoredCrossSection(BaseFileName, ThisCrossSectionPlane, Plane, Index1Low, Index1High, Index2Low,
                                        Index2High, CrossSectionOutOfPlaneLocation, GrainID, GrainEulerAngles, deltax,
                                        NumberOfOrientations);
        // Collect grain orientation frequency data and write to a file to be read by MTEX/plotted as pole figures
        if (PrintSectionPF[n]) {
            ViewI_H GOHistogram = getOrientationHistogram(NumberOfOrientations, CrossSectionGrainIDs, CrossSectionSize);
            writePoleFigure(BaseFileName, ThisCrossSectionPlane, NumberOfOrientations, GrainEulerAngles, GOHistogram);
        }
        // Make list of unique grains and corresponding grain areas
        int NumberOfGrains;
        std::vector<int> UniqueGrainIDs = getUniqueGrains(CrossSectionGrainIDs, NumberOfGrains);
        std::cout << "The number of grains in this cross-section is " << NumberOfGrains << std::endl;
        std::vector<int> GrainAreas(NumberOfGrains, 0);
        for (int i = 0; i < NumberOfGrains; i++) {
            for (int j = 0; j < CrossSectionSize; j++) {
                if (UniqueGrainIDs[i] == CrossSectionGrainIDs[j])
                    GrainAreas[i]++;
            }
        }

        // Should these grains be analyzed as a single distribution of grain areas, or a bimodal distribution
        // (bimodal distribution option also includes printing of misorientation data)
        if (BimodalAnalysis[n])
            AnalyzeCrossSection_Bimodal(QoIs, BaseFileName, ThisCrossSectionPlane, deltax, NumberOfGrains,
                                        CrossSectionSize, UniqueGrainIDs, GrainAreas, NumberOfOrientations,
                                        GrainUnitVector, GrainRGBValues);
        else
            AnalyzeCrossSection_Unimodal(QoIs, BaseFileName, ThisCrossSectionPlane, deltax, NumberOfGrains,
                                         CrossSectionSize, GrainAreas);
    }
    QoIs.close();
}

//*****************************************************************************/
void writePoleFigure(std::string BaseFileName, std::string RegionLabel, int NumberOfOrientations,
                     ViewF_H GrainEulerAngles, ViewI_H GOHistogram) {

    // Using new format, write pole figure data to "Filename"
    std::string Filename = BaseFileName + "_" + RegionLabel + "_PoleFigureData.txt";
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

//*****************************************************************************/
// For the region bounded by [Index1Low,Index1High] and [Index2Low,Index2High], at out of plane location given by
// CrossSectionOutOfPlaneLocation, print data to be read by MTEX to plot the cross-section using the inverse pole figure
// colormap. Identities of the in plane and out of plane indices depend on the value for "Plane"
void writeIPFColoredCrossSection(std::string BaseFileName, std::string CrossSectionLabel, std::string Plane,
                                 int Index1Low, int Index1High, int Index2Low, int Index2High,
                                 int CrossSectionOutOfPlaneLocation, ViewI3D_H GrainID, ViewF_H GrainEulerAngles,
                                 double deltax, int NumberOfOrientations) {

    std::string FNameIPF = BaseFileName + "_" + CrossSectionLabel + "_IPFCrossSectionData.txt";
    std::ofstream GrainplotIPF;
    GrainplotIPF.open(FNameIPF);
    GrainplotIPF << std::fixed << std::setprecision(6);
    for (int Index1 = Index1Low; Index1 < Index1High; Index1++) {
        for (int Index2 = Index2Low; Index2 < Index2High; Index2++) {
            int Index3 = CrossSectionOutOfPlaneLocation;
            // How do Index1, Index2, Index3 correspond to GrainID(Z loc, X loc, Yloc)?
            int ZLoc, XLoc, YLoc;
            if (Plane == "XY") {
                XLoc = Index1;
                YLoc = Index2;
                ZLoc = Index3;
            }
            else if (Plane == "YZ") {
                XLoc = Index3;
                YLoc = Index1;
                ZLoc = Index2;
            }
            else if (Plane == "XZ") {
                XLoc = Index1;
                YLoc = Index3;
                ZLoc = Index2;
            }
            else
                throw std::runtime_error(
                    "Error: Unknown plane input for WriteIPFColoredCrossSection: should be XY, YZ, or XZ");
            // What orientation does this grain id correspond to? Should be between 0 and NumberOfOrientations-1
            int GOVal = (abs(GrainID(ZLoc, XLoc, YLoc)) - 1) % NumberOfOrientations;
            // The grain structure is phase "1" - any unindexed points with GOVal = -1 (which are possible from regions
            // that didn't undergo melting) are assigned phase "0"
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

//*****************************************************************************/
void writeExaConstitRVE(int NumberOfRVEs, std::string BaseFileName, int, int, int, double deltax, ViewI3D_H GrainID,
                        std::vector<int> XLow_RVE, std::vector<int> XHigh_RVE, std::vector<int> YLow_RVE,
                        std::vector<int> YHigh_RVE, std::vector<int> ZLow_RVE, std::vector<int> ZHigh_RVE) {

    // Loop over each RVE specified in the file "AnalysisOutputs.txt"
    for (int n = 0; n < NumberOfRVEs; n++) {
        // Print data for ExaConstit
        std::string FNameE = BaseFileName + std::to_string(n) + "_ExaConstit.csv";
        std::cout << "RVE number " << n + 1 << " with X coordinates " << XLow_RVE[n] << "," << XHigh_RVE[n]
                  << "; Y coordinates " << YLow_RVE[n] << "," << YHigh_RVE[n] << "; Z coordinates " << ZLow_RVE[n]
                  << "," << ZHigh_RVE[n] << " being printed to file " << FNameE << " for ExaConstit" << std::endl;
        std::ofstream GrainplotE;
        GrainplotE.open(FNameE);
        GrainplotE << "Coordinates are in CA units, 1 cell = " << deltax << " m. Data is cell-centered. Origin at "
                   << XLow_RVE[n] << "," << YLow_RVE[n] << "," << ZLow_RVE[n] << " , domain size is "
                   << XHigh_RVE[n] - XLow_RVE[n] + 1 << " by " << YHigh_RVE[n] - YLow_RVE[n] + 1 << " by "
                   << ZHigh_RVE[n] - ZLow_RVE[n] + 1 << " cells" << std::endl;
        GrainplotE << "X coord, Y coord, Z coord, Grain ID" << std::endl;
        int NucleatedGrainCells = 0;
        for (int k = ZLow_RVE[n]; k <= ZHigh_RVE[n]; k++) {
            for (int j = YLow_RVE[n]; j <= YHigh_RVE[n]; j++) {
                for (int i = XLow_RVE[n]; i <= XHigh_RVE[n]; i++) {
                    GrainplotE << i << "," << j << "," << k << "," << GrainID(k, i, j) << std::endl;
                    if (GrainID(k, i, j) < 0)
                        NucleatedGrainCells++;
                }
            }
        }
        GrainplotE.close();
        int RVESize = (XHigh_RVE[n] - XLow_RVE[n]) * (YHigh_RVE[n] - YLow_RVE[n]) * (ZHigh_RVE[n] - ZLow_RVE[n]);
        double RVEFractNucleatedGrains = DivideCast<double>(NucleatedGrainCells, RVESize);
        std::cout << "The fraction of grains formed via nucleation events in this RVE is " << RVEFractNucleatedGrains
                  << std::endl;
    }
}
