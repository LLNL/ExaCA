// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT
#include "GAprint.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

//*****************************************************************************/
void PrintMisorientationData(bool *AnalysisTypes, std::string BaseFileName, int XMin, int XMax, int YMin, int YMax,
                             int ZMin, int ZMax, ViewI3D_H LayerID, ViewF_H GrainUnitVector, ViewI3D_H GrainID,
                             int NumberOfOrientations) {

    // Frequency of misorientations in the selected region
    std::ofstream MisorientationPlot;
    std::string FNameM = BaseFileName + "_MisorientationFrequency.csv";
    if (AnalysisTypes[0]) {
        MisorientationPlot.open(FNameM);
        std::cout << "Printing file " << FNameM
                  << " of grain misorientations relative to the +Z direction, for all cells in selected volume (not "
                     "including substrate)"
                  << std::endl;
    }
    int NumberOfMeltedCells = 0;
    long double MisorientationSum = 0.0;
    int NumberOfMeltedCellsTop = 0;
    long double MisorientationSumTop = 0.0;
    ViewF_H GrainMisorientation(Kokkos::ViewAllocateWithoutInitializing("GrainMisorientation"), NumberOfOrientations);
    for (int n = 0; n < NumberOfOrientations; n++) {
        double AngleZmin = 62.7;
        for (int ll = 0; ll < 3; ll++) {
            double AngleZ = std::abs((180 / M_PI) * acos(GrainUnitVector(9 * n + 3 * ll + 2)));
            if (AngleZ < AngleZmin) {
                AngleZmin = AngleZ;
            }
        }
        GrainMisorientation(n) = AngleZmin;
    }
    for (int k = ZMin; k <= ZMax; k++) {
        for (int i = XMin; i <= XMax; i++) {
            for (int j = YMin; j <= YMax; j++) {
                // Only take data from cells in the representative area that underwent melting (LayerID >= 0)
                if (LayerID(k, i, j) != -1) {
                    int MyOrientation = ((abs(GrainID(k, i, j)) - 1) % NumberOfOrientations);
                    float MyMisorientation = GrainMisorientation(MyOrientation);
                    if (AnalysisTypes[0])
                        MisorientationPlot << MyMisorientation << std::endl;
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
    if (AnalysisTypes[0])
        MisorientationPlot.close();
    std::cout << "Within the representative region consisting of "
              << (XMax - XMin + 1) * (YMax - YMin + 1) * (ZMax - ZMin + 1) << " cells, " << NumberOfMeltedCells
              << " underwent melting and:" << std::endl;
    std::cout << "-- The mean misorientation relative to the +Z direction is "
              << MisorientationSum / ((double)(NumberOfMeltedCells)) << " degrees" << std::endl;
    std::cout << "-- The mean misorientation relative to the +Z direction at the representative region top (Z = "
              << ZMax << ") is " << MisorientationSumTop / ((double)(NumberOfMeltedCellsTop)) << " degrees"
              << std::endl;
}

//*****************************************************************************/

void PrintSizeData(bool *AnalysisTypes, std::string BaseFileName, int XMin, int XMax, int YMin, int YMax, int ZMin,
                   int ZMax, int nx, int ny, int nz, ViewI3D_H, ViewI3D_H GrainID, double deltax) {

    // Get vector of unique GrainIDs
    int Counter = 0;
    int NucleatedGrainCells = 0;
    int DomainVol = (ZMax - ZMin + 1) * (YMax - YMin + 1) * (XMax - XMin + 1);
    std::vector<int> UniqueGrainIDs(DomainVol);
    for (int k = ZMin; k <= ZMax; k++) {
        for (int i = XMin; i <= XMax; i++) {
            for (int j = YMin; j <= YMax; j++) {
                UniqueGrainIDs[Counter] = GrainID(k, i, j);
                Counter++;
                if (GrainID(k, i, j) < 0)
                    NucleatedGrainCells++;
            }
        }
    }
    std::sort(UniqueGrainIDs.begin(), UniqueGrainIDs.end());
    std::vector<int>::iterator it;
    it = std::unique(UniqueGrainIDs.begin(), UniqueGrainIDs.end());
    UniqueGrainIDs.resize(std::distance(UniqueGrainIDs.begin(), it));
    int NumberOfGrains = UniqueGrainIDs.size();
    std::cout << "-- There are " << NumberOfGrains << " grains in this volume, and the mean grain volume is "
              << deltax * deltax * pow(10, 12) * ((double)(DomainVol) / (double)(NumberOfGrains)) << " cubic microns"
              << std::endl;
    std::cout << "-- The volume fraction consisting of nucleated grains is "
              << ((double)(NucleatedGrainCells)) / ((double)(DomainVol)) << std::endl;

    float *AspectRatio = new float[NumberOfGrains];
    int *VolGrain = new int[NumberOfGrains];
    float *GrainHeight = new float[NumberOfGrains];
    float ARSum = 0.0;
    float VolWtARSum = 0.0;
    float GrainHeightSum = 0.0;
    float GrainWidthSum = 0.0;
    for (int n = 0; n < NumberOfGrains; n++) {
        AspectRatio[n] = 0;
        VolGrain[n] = 0;
        int TempTopX = 0;
        int TempTopY = 0;
        int TempTopZ = 0;
        int TempBottomX = nx - 1;
        int TempBottomY = ny - 1;
        int TempBottomZ = nz - 1;
        int ThisGrainID = UniqueGrainIDs[n];
        for (int k = ZMin; k <= ZMax; k++) {
            for (int i = XMin; i <= XMax; i++) {
                for (int j = YMin; j <= YMax; j++) {
                    // Only take data from cells in the representative area that underwent melting
                    if (GrainID(k, i, j) == ThisGrainID) {
                        VolGrain[n]++;
                        if (i > TempTopX)
                            TempTopX = i;
                        if (i < TempBottomX)
                            TempBottomX = i;
                        if (j > TempTopY)
                            TempTopY = j;
                        if (j < TempBottomY)
                            TempBottomY = j;
                        if (k > TempTopZ)
                            TempTopZ = k;
                        if (k < TempBottomZ)
                            TempBottomZ = k;
                    }
                }
            }
        }
        GrainHeight[n] = TempTopZ - TempBottomZ + 1;
        float GrainWidthX = TempTopX - TempBottomX + 1;
        float GrainWidthY = TempTopY - TempBottomY + 1;
        GrainHeightSum += GrainHeight[n];
        GrainWidthSum += 0.5 * (GrainWidthX + GrainWidthY);
        float AR_XZ = (float)(TempTopZ - TempBottomZ + 1) / (float)(TempTopX - TempBottomX + 1);
        float AR_YZ = (float)(TempTopZ - TempBottomZ + 1) / (float)(TempTopY - TempBottomY + 1);
        AspectRatio[n] = 0.5 * (AR_XZ + AR_YZ);
        ARSum += AspectRatio[n];
        VolWtARSum += AspectRatio[n] * (float)(VolGrain[n]);
    }
    std::cout << "-- The mean grain aspect ratio (Z direction to transverse) is " << ARSum / (float)(NumberOfGrains)
              << std::endl;
    std::cout << "-- The mean volume-weighted grain aspect ratio (Z direction to transverse) is "
              << VolWtARSum / ((float)(DomainVol)) << std::endl;
    std::cout << "-- The mean grain height is " << GrainHeightSum / (float)(NumberOfGrains) << " microns" << std::endl;
    std::cout << "-- The mean grain width is " << GrainWidthSum / (float)(NumberOfGrains) << " microns" << std::endl;

    if (AnalysisTypes[1]) {
        std::ofstream VolumePlot;
        std::string FNameV = BaseFileName + "_VolumeFrequency.csv";
        VolumePlot.open(FNameV);
        std::cout << "Printing file " << FNameV << " of grain volumes (in cubic microns) in selected volume"
                  << std::endl;
        for (int n = 0; n < NumberOfGrains; n++) {
            VolumePlot << VolGrain[n] * deltax * deltax * pow(10, 12) << std::endl;
        }
        VolumePlot.close();
    }
    if (AnalysisTypes[2]) {
        std::ofstream ARPlot;
        std::string FNameAR = BaseFileName + "_AspectRatioFrequency.csv";
        ARPlot.open(FNameAR);
        std::cout << "Printing file " << FNameAR << " of grain aspect ratios in selected volume" << std::endl;
        for (int n = 0; n < NumberOfGrains; n++) {
            ARPlot << AspectRatio[n] << std::endl;
        }
        ARPlot.close();
    }
    if (AnalysisTypes[5]) {
        std::ofstream GrainHeightPlot;
        std::string FNameGH = BaseFileName + "_GrainHeightDistribution.csv";
        std::cout << "Printing file " << FNameGH << " of mean height distribution (in microns) for the selected volume"
                  << std::endl;
        GrainHeightPlot.open(FNameGH);
        for (int n = 0; n < NumberOfGrains; n++) {
            GrainHeightPlot << GrainHeight[n] << std::endl;
        }
        GrainHeightPlot.close();
    }
}

//*****************************************************************************/
void PrintGrainAreaData(bool *AnalysisTypes, std::string BaseFileName, double deltax, int XMin, int XMax, int YMin,
                        int YMax, int ZMin, int ZMax, ViewI3D_H GrainID) {

    std::string FName1 = BaseFileName + "_GrainAreas.csv";
    std::string FName2 = BaseFileName + "_WeightedGrainAreas.csv";
    std::ofstream Grainplot1, Grainplot2;
    if (AnalysisTypes[3]) {
        std::cout << "Printing file " << FName1 << " of grain area values (in square microns) for all Z coordinates"
                  << std::endl;
        Grainplot1.open(FName1);
    }
    if (AnalysisTypes[4]) {
        std::cout << "Printing file " << FName2
                  << " of weighted grain area values (in square microns) for every 5th Z coordinate" << std::endl;
        Grainplot2.open(FName2);
    }

    if (!(AnalysisTypes[2]) && (!(AnalysisTypes[3])))
        ZMin = ZMax; // only print grain area/weighted grain area to screen at this one Z coordinate

    int LayerArea = (XMax - XMin + 1) * (YMax - YMin + 1);
    for (int k = ZMin; k <= ZMax; k++) {
        std::vector<int> GIDAllVals_ThisLayer(LayerArea);
        int Counter = 0;
        for (int i = XMin; i <= XMax; i++) {
            for (int j = YMin; j <= YMax; j++) {
                GIDAllVals_ThisLayer[Counter] = GrainID(k, i, j);
                Counter++;
            }
        }
        std::vector<int> GIDVals_ThisLayer;
        GIDVals_ThisLayer = GIDAllVals_ThisLayer;
        std::vector<int>::iterator it;
        std::sort(GIDVals_ThisLayer.begin(), GIDVals_ThisLayer.end());
        it = std::unique(GIDVals_ThisLayer.begin(), GIDVals_ThisLayer.end());
        GIDVals_ThisLayer.resize(std::distance(GIDVals_ThisLayer.begin(), it));
        int GrainsThisLayer = GIDVals_ThisLayer.size();
        double MeanGrainAreaThisLayer = (double)(LayerArea) / (double)(GrainsThisLayer);
        if (k == ZMax) {
            std::cout << "Number of grains at the representative region top (Z = " << ZMax << "): " << GrainsThisLayer
                      << std::endl;
            if (AnalysisTypes[6]) {
                std::string FName3 = BaseFileName + "_GrainWidthDistributionX.csv";
                std::string FName4 = BaseFileName + "_GrainWidthDistributionY.csv";
                std::ofstream Grainplot3;
                std::ofstream Grainplot4;
                std::cout << "Printing files " << FName3 << " and " << FName4
                          << " of grain width distributions in x and y (in microns) at Z = " << ZMax << std::endl;
                Grainplot3.open(FName3);
                Grainplot4.open(FName4);
                for (int n = 0; n < GrainsThisLayer; n++) {
                    int ThisGrainID = GIDVals_ThisLayer[n];
                    int TempTopX = XMin;
                    int TempTopY = YMin;
                    int TempBottomX = XMax;
                    int TempBottomY = YMax;
                    for (int i = XMin; i <= XMax; i++) {
                        for (int j = YMin; j <= YMax; j++) {
                            if (GrainID(k, i, j) == ThisGrainID) {
                                if (i > TempTopX)
                                    TempTopX = i;
                                if (i < TempBottomX)
                                    TempBottomX = i;
                                if (j > TempTopY)
                                    TempTopY = j;
                                if (j < TempBottomY)
                                    TempBottomY = j;
                            }
                        }
                    }
                    float GrainExtentX = TempTopX - TempBottomX + 1;
                    float GrainExtentY = TempTopY - TempBottomY + 1;
                    Grainplot3 << GrainExtentX * deltax / pow(10, -6) << std::endl;
                    Grainplot4 << GrainExtentY * deltax / pow(10, -6) << std::endl;
                }
                Grainplot3.close();
                Grainplot4.close();
            }
        }
        if (((AnalysisTypes[3]) && (k % 5 == 0)) || (k == ZMax)) {
            long int AreaXArea = 0;
            for (int l = 0; l < GrainsThisLayer; l++) {
                long int MyGrainArea = 0;
                for (int ll = 0; ll < LayerArea; ll++) {
                    if (GIDVals_ThisLayer[l] == GIDAllVals_ThisLayer[ll])
                        MyGrainArea++;
                }
                AreaXArea += MyGrainArea * MyGrainArea;
            }
            double WeightedArea = ((double)(AreaXArea) / (double)(LayerArea));
            if (AnalysisTypes[3])
                Grainplot2 << WeightedArea * deltax * deltax / pow(10, -12) << std::endl;
            if (k == ZMax)
                std::cout << "-- The mean weighted grain area at the representative region top (Z = " << ZMax << ") is "
                          << WeightedArea * deltax * deltax / pow(10, -12) << " square microns" << std::endl;
        }
        if (AnalysisTypes[4])
            Grainplot1 << MeanGrainAreaThisLayer * deltax * deltax / pow(10, -12) << std::endl;
        if (k == ZMax)
            std::cout << "-- The mean grain area at the representative region top (Z = " << ZMax << ") is "
                      << MeanGrainAreaThisLayer * deltax * deltax / pow(10, -12) << " square microns" << std::endl;
    }
    if (AnalysisTypes[3])
        Grainplot1.close();
    if (AnalysisTypes[4])
        Grainplot2.close();
}

//*****************************************************************************/
void PrintPoleFigureData(bool *AnalysisTypes, std::string BaseFileName, int NumberOfOrientations, int XMin, int XMax,
                         int YMin, int YMax, int ZMin, int ZMax, ViewI3D_H GrainID, ViewI3D_H LayerID,
                         bool NewOrientationFormatYN, ViewF_H GrainEulerAngles) {

    if (AnalysisTypes[7]) {

        // Histogram of orientations for texture determination
        ViewI_H GOHistogram("GOHistogram", NumberOfOrientations);
        // frequency data on grain ids
        for (int k = ZMin; k <= ZMax; k++) {
            for (int j = YMin; j <= YMax; j++) {
                for (int i = XMin; i <= XMax; i++) {
                    if (LayerID(k, i, j) != -1) {
                        int GOVal = (abs(GrainID(k, i, j)) - 1) % NumberOfOrientations;
                        GOHistogram(GOVal)++;
                    }
                }
            }
        }
        if (NewOrientationFormatYN) {
            // Write pole figure data for this region using the new format
            std::string FNameM = BaseFileName + "_PFVolumeX" + std::to_string(XMin) + "-" + std::to_string(XMax) + "Y" +
                                 std::to_string(YMin) + "-" + std::to_string(YMax) + "Z" + std::to_string(ZMin) + "-" +
                                 std::to_string(ZMax) + ".txt";
            WritePoleFigureDataToFile(FNameM, NumberOfOrientations, GrainEulerAngles, GOHistogram);
        }
        else {
            // Write pole figure data for this region using the old format
            std::string FNameM = BaseFileName + "_MTEXOrientations.csv";
            WritePoleFigureDataToFile_OldFormat(FNameM, NumberOfOrientations, GOHistogram);
        }
    }
}

//*****************************************************************************/
void PrintCrossSectionOrientationData(int NumberOfCrossSections, std::string BaseFileName,
                                      std::vector<int> CrossSectionPlane, std::vector<int> CrossSectionLocation, int nx,
                                      int ny, int nz, int NumberOfOrientations, ViewI3D_H GrainID,
                                      std::vector<bool> PrintSectionPF, std::vector<bool> PrintSectionIPF,
                                      bool NewOrientationFormatYN, double deltax, ViewF_H, ViewF_H GrainEulerAngles) {

    // Loop over each cross-section specified in the file "AnalysisOutputs.txt"
    for (int n = 0; n < NumberOfCrossSections; n++) {
        // Print data for pyEBSD/MTEX
        std::string ThisCrossSectionPlane; // Which kind of cross-section is this?
        int Index1Low = 0;
        int Index2Low = 0;
        int Index1High, Index2High; // Values depend on the cross-section axes: nx, ny, or nz
        int ThisCrossSectionLocation = CrossSectionLocation[n]; // Cross-section location out-of-plane
        if (CrossSectionPlane[n] == 0) {
            ThisCrossSectionPlane = "XZ";
            Index1High = nx;
            Index2High = nz;
        }
        else if (CrossSectionPlane[n] == 1) {
            ThisCrossSectionPlane = "YZ";
            Index1High = ny;
            Index2High = nz;
        }
        else {
            ThisCrossSectionPlane = "XY";
            Index1High = nx;
            Index2High = ny;
        }

        // Should pole figure data be printed for this cross-section?
        // Should inverse pole figure-mapping data be printed for this cross-section?
        if ((PrintSectionPF[n]) || (PrintSectionIPF[n])) {
            // If the option was given in Section 2 analysis file, pole figure data here should be printed with the new
            // format
            std::cout << "Printing cross-section data for " << ThisCrossSectionPlane << " and with "
                      << ThisCrossSectionLocation << " as the out of plane location" << std::endl;
            std::string FNameIPF = BaseFileName + "-" + ThisCrossSectionPlane +
                                   std::to_string(ThisCrossSectionLocation) + "_IPFCrossSection.txt";
            std::ofstream GrainplotIPF;
            ViewI_H GOHistogram("GOHistogram", NumberOfOrientations);
            if (PrintSectionIPF[n]) {
                GrainplotIPF.open(FNameIPF);
                GrainplotIPF << std::fixed << std::setprecision(6);
            }
            int NucleatedGrainCells = 0;
            int CrossSectionSize = (Index1High - Index1Low) * (Index2High - Index2Low);
            for (int Index1 = Index1Low; Index1 < Index1High; Index1++) {
                for (int Index2 = Index2Low; Index2 < Index2High; Index2++) {
                    int Index3 = ThisCrossSectionLocation;
                    // How do Index1, Index2, Index3 correspond to GrainID(Z loc, X loc, Yloc)?
                    int ZLoc, XLoc, YLoc;
                    if (ThisCrossSectionPlane == "XY") {
                        XLoc = Index1;
                        YLoc = Index2;
                        ZLoc = Index3;
                    }
                    else if (ThisCrossSectionPlane == "YZ") {
                        XLoc = Index3;
                        YLoc = Index1;
                        ZLoc = Index2;
                    }
                    else {
                        XLoc = Index1;
                        YLoc = Index3;
                        ZLoc = Index2;
                    }
                    int GOVal = (abs(GrainID(ZLoc, XLoc, YLoc)) - 1) % NumberOfOrientations;
                    // Count number of cells in this cross-section have GrainID < 0 (grains formed via nucleation)
                    if (GrainID(ZLoc, XLoc, YLoc) < 0)
                        NucleatedGrainCells++;
                    // If constructing pole figure data from these orientations, add this value to the frequency data
                    if (PrintSectionPF[n])
                        GOHistogram(GOVal)++;
                    if (PrintSectionIPF[n]) {
                        if (NewOrientationFormatYN) {
                            // The grain structure is phase "1" - any unindexed points (which are possible from regions
                            // that didn't undergo melting) are assigned phase "0"
                            if (GOVal == -1)
                                GrainplotIPF << "0 0 0 0 " << Index1 * deltax * pow(10, 6) << " "
                                             << Index2 * deltax * pow(10, 6) << std::endl;
                            else
                                GrainplotIPF << GrainEulerAngles(3 * GOVal) << " " << GrainEulerAngles(3 * GOVal + 1)
                                             << " " << GrainEulerAngles(3 * GOVal + 2) << " 1 "
                                             << Index1 * deltax * pow(10, 6) << " " << Index2 * deltax * pow(10, 6)
                                             << std::endl;
                        }
                        else {
                            GrainplotIPF << Index1 * deltax * pow(10, 6) << "," << Index2 * deltax * pow(10, 6) << ","
                                         << GOVal << std::endl;
                        }
                    }
                }
            }
            std::cout << "The fraction of grains in this cross-section formed via nucleation events is "
                      << static_cast<double>(NucleatedGrainCells) / static_cast<double>(CrossSectionSize) << std::endl;
            if (PrintSectionIPF[n])
                GrainplotIPF.close();
            if (PrintSectionPF[n]) {
                std::string FNamePF = BaseFileName + "-" + ThisCrossSectionPlane +
                                      std::to_string(ThisCrossSectionLocation) + "_PFCrossSection.txt";
                WritePoleFigureDataToFile(FNamePF, NumberOfOrientations, GrainEulerAngles, GOHistogram);
            }
        }
    }
}

//*****************************************************************************/
void WritePoleFigureDataToFile(std::string Filename, int NumberOfOrientations, ViewF_H GrainEulerAngles,
                               ViewI_H GOHistogram) {

    // Using new format, write pole figure data to "Filename"
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
void WritePoleFigureDataToFile_OldFormat(std::string Filename, int NumberOfOrientations, ViewI_H GOHistogram) {

    // Using old format, write histogram data (used by a second post-processing script to construct pole figures) to
    // "Filename"
    std::ofstream MTEXPlot;
    MTEXPlot.open(Filename);
    for (int i = 0; i < NumberOfOrientations; i++) {
        MTEXPlot << GOHistogram[i] << std::endl;
    }
    MTEXPlot.close();
}

//*****************************************************************************/
void PrintExaConstitRVEData(int NumberOfRVEs, std::string BaseFileName, int, int, int, double deltax, ViewI3D_H GrainID,
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
            for (int i = XLow_RVE[n]; i <= XHigh_RVE[n]; i++) {
                for (int j = YLow_RVE[n]; j <= YHigh_RVE[n]; j++) {
                    GrainplotE << i << "," << j << "," << k << "," << GrainID(k, i, j) << std::endl;
                    if (GrainID(k, i, j) < 0)
                        NucleatedGrainCells++;
                }
            }
        }
        GrainplotE.close();
        int RVESize = (XHigh_RVE[n] - XLow_RVE[n]) * (YHigh_RVE[n] - YLow_RVE[n]) * (ZHigh_RVE[n] - ZLow_RVE[n]);
        std::cout << "The fraction of grains formed via nucleation events in this RVE is "
                  << static_cast<double>(NucleatedGrainCells) / static_cast<double>(RVESize) << std::endl;
    }
}
