// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "GAprint.hpp"
#include "GAutils.hpp"

//*****************************************************************************/
// "print" routines print average quantities to the screen (and to the file stream specified by QoIs)
// "write" routines write data to specific files
//*****************************************************************************/
// Print information about a representative area to the console/QoIs file
void printAnalysisHeader_Area(std::ofstream &QoIs, const int XLow_cells, const int XHigh_cells, const int YLow_cells,
                              const int YHigh_cells, const int ZLow_cells, const int ZHigh_cells,
                              const double XLow_Meters, const double XHigh_Meters, const double YLow_Meters,
                              const double YHigh_Meters, const double ZLow_Meters, const double ZHigh_Meters,
                              const std::string regionName, const std::string regionOrientation) {

    std::cout << "Stats for " << regionName << " plane:" << std::endl;
    QoIs << "Stats for " << regionName << " plane:" << std::endl;
    if (regionOrientation == "XY") {
        std::cout << "The representative area is located at Z = " << ZLow_Meters << " m" << std::endl;
        std::cout << "(in CA units, Z = " << ZLow_cells << ")" << std::endl;
        std::cout << "The representative area is bounded by the region spanning X = [" << XLow_Meters << ","
                  << XHigh_Meters << "], Y = [ " << YLow_Meters << "," << YHigh_Meters << "] m" << std::endl;
        std::cout << "(in CA units X = [" << XLow_cells << "," << XHigh_cells << "], Y = [" << YLow_cells << ","
                  << YHigh_cells << "])" << std::endl;
        QoIs << "The representative area is located at Z = " << ZLow_Meters << " m" << std::endl;
        QoIs << "(in CA units, Z = " << ZLow_cells << ")" << std::endl;
        QoIs << "The representative area is bounded by the region spanning X = [" << XLow_Meters << "," << XHigh_Meters
             << "], Y = [" << YLow_Meters << "," << YHigh_Meters << "] m" << std::endl;
        QoIs << "(in CA units X = [" << XLow_cells << "," << XHigh_cells << "], Y = [" << YLow_cells << ","
             << YHigh_cells << "])" << std::endl;
    }
    else if (regionOrientation == "XZ") {
        std::cout << "The representative area is located at Y = " << YLow_Meters << " m" << std::endl;
        std::cout << "(in CA units, Y = " << YLow_cells << ")" << std::endl;
        std::cout << "The representative area is bounded by the region spanning X = [" << XLow_Meters << ","
                  << XHigh_Meters << "], Z = [" << ZLow_Meters << "," << ZHigh_Meters << "] m" << std::endl;
        std::cout << "(in CA units X = [" << XLow_cells << "," << XHigh_cells << "], Z = [" << ZLow_cells << ","
                  << ZHigh_cells << "])" << std::endl;
        QoIs << "The representative area is located at Y = " << YLow_Meters << " m" << std::endl;
        QoIs << "(in CA units, Y = " << YLow_cells << ")" << std::endl;
        QoIs << "The representative area is bounded by the region spanning X = [" << XLow_Meters << "," << XHigh_Meters
             << "], Z = [" << ZLow_Meters << "," << ZHigh_Meters << "] m" << std::endl;
        QoIs << "(in CA units X = [" << XLow_cells << "," << XHigh_cells << "], Z = [" << ZLow_cells << ","
             << ZHigh_cells << "])" << std::endl;
    }
    else if (regionOrientation == "YZ") {
        std::cout << "The representative area is located at X = " << XLow_Meters << " m" << std::endl;
        std::cout << "(in CA units, X = " << XLow_cells << ")" << std::endl;
        std::cout << "The representative area is bounded by the region spanning Y = [" << YLow_Meters << ","
                  << YHigh_Meters << "], Z = [" << ZLow_Meters << "," << ZHigh_Meters << "] m" << std::endl;
        std::cout << "(in CA units Y = [" << YLow_cells << "," << YHigh_cells << "], Z = [" << ZLow_cells << ","
                  << ZHigh_cells << "])" << std::endl;
        QoIs << "The representative area is located at X = " << XLow_Meters << " m" << std::endl;
        QoIs << "(in CA units, X = " << XLow_cells << ")" << std::endl;
        QoIs << "The representative area is bounded by the region spanning Y = [" << YLow_Meters << "," << YHigh_Meters
             << "], Z = [" << ZLow_Meters << "," << ZHigh_Meters << "] m" << std::endl;
        QoIs << "(in CA units Y = [" << YLow_cells << "," << YHigh_cells << "], Z = [" << ZLow_cells << ","
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
                              double RepresentativeRegionSize_Meters) {

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
    float AvgMisorientationX = DivideCast<float>(GrainMisorientation_sum_x, RepresentativeRegionSize_Meters);
    float AvgMisorientationY = DivideCast<float>(GrainMisorientation_sum_y, RepresentativeRegionSize_Meters);
    float AvgMisorientationZ = DivideCast<float>(GrainMisorientation_sum_z, RepresentativeRegionSize_Meters);
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

// Print the average grain size and the number of grains in the region
void printMeanSize(std::ofstream &QoIs, int NumberOfGrains, double RepresentativeRegionSize_Microns,
                   std::string RegionType, std::string Units) {

    double AvgSizePerGrain = DivideCast<double>(RepresentativeRegionSize_Microns, NumberOfGrains);
    std::cout << "-- There are " << NumberOfGrains << " grains in this " << RegionType << " , and the mean grain "
              << RegionType << " is " << AvgSizePerGrain << " " << Units << std::endl;
    QoIs << "-- There are " << NumberOfGrains << " grains in this " << RegionType << " , and the mean grain "
         << RegionType << " is " << AvgSizePerGrain << " " << Units << std::endl;
}

// Print average aspect ratio in the build to the average of the transverse directions
void printMeanBuildTransAspectRatio(std::ofstream &QoIs, std::vector<float> GrainExtentX,
                                    std::vector<float> GrainExtentY, std::vector<float> GrainExtentZ,
                                    std::vector<float> GrainSizeVector, double RepresentativeRegionSize_Meters,
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
    double RepresentativeRegionSize_Microns = RepresentativeRegionSize_Meters * pow(10, 18);
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
void printMeanExtent(std::ofstream &QoIs, std::vector<float> GrainExtent_Microns, std::string Direction,
                     int NumberOfGrains) {

    float GrainExtentSum = 0.0;
    for (int n = 0; n < NumberOfGrains; n++)
        GrainExtentSum += GrainExtent_Microns[n];
    float AvgGrainExtent = DivideCast<float>(GrainExtentSum, NumberOfGrains);
    std::cout << "-- The mean grain extent in the " << Direction << " direction is " << AvgGrainExtent << " microns"
              << std::endl;
    QoIs << "-- The mean grain extent in the " << Direction << " direction is " << AvgGrainExtent << " microns"
         << std::endl;
}

// Write unweighted and/or weighted grain areas as a function of build height to file(s)
void writeAreaSeries(bool PrintWeightedAreas, bool PrintUnweightedAreas, std::string BaseFileName, double deltax,
                     int XMin, int XMax, int YMin, int YMax, int ZMin, int ZMax, ViewI3D_H GrainID,
                     double ZMin_Meters) {

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
        std::cout << "Note: Option to print weighted grain area data will be removed in a future release" << std::endl;
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
            Grainplot1 << ZMin_Meters * pow(10, 6) + convertToMicrons(deltax, "length") << ","
                       << MeanGrainAreaThisLayer * convertToMicrons(deltax, "area") << std::endl;
        if ((PrintWeightedAreas) && (k % 5 == 0)) {
            std::vector<float> GrainSizeVector_Area =
                getGrainSizes(GrainIDVector_Area, UniqueGrainIDVector_Area, NumberOfGrains_Area, deltax, "area");
            float AreaXArea = 0.0;
            for (int n = 0; n < NumberOfGrains_Area; n++)
                AreaXArea += GrainSizeVector_Area[n] * GrainSizeVector_Area[n];
            double WeightedArea = DivideCast<double>(AreaXArea, LayerArea);
            Grainplot2 << ZMin_Meters * pow(10, 6) + convertToMicrons(deltax, "length") << "," << WeightedArea
                       << std::endl;
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
