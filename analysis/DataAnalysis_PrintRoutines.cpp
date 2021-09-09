// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT
#include "header.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>

//*****************************************************************************/
void PrintMisorientationData(bool *AnalysisTypes, std::string BaseFileName, int XMin, int XMax, int YMin, int YMax,
                             int ZMin, int ZMax, int ***Melted, float ***GrainUnitVector, int ***GrainID,
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
    for (int k = ZMin; k <= ZMax; k++) {
        for (int i = XMin; i <= XMax; i++) {
            for (int j = YMin; j <= YMax; j++) {
                // Only take data from cells in the representative area that underwent melting
                if (Melted[k][i][j] == 1) {
                    int MyOrientation = ((abs(GrainID[k][i][j]) - 1) % NumberOfOrientations);
                    double AngleZmin = 62.7;
                    for (int ll = 0; ll < 3; ll++) {
                        double AngleZ = abs((180 / M_PI) * acos(GrainUnitVector[MyOrientation][ll][2]));
                        if (AngleZ < AngleZmin) {
                            AngleZmin = AngleZ;
                        }
                    }
                    if (AnalysisTypes[0])
                        MisorientationPlot << AngleZmin << std::endl;
                    MisorientationSum += AngleZmin;
                    if (k == ZMax) {
                        MisorientationSumTop += AngleZmin;
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
                   int ZMax, int nx, int ny, int nz, int ***, int ***GrainID, double deltax) {

    // Get vector of unique GrainIDs
    int Counter = 0;
    int NucleatedGrainCells = 0;
    int DomainVol = (ZMax - ZMin + 1) * (YMax - YMin + 1) * (XMax - XMin + 1);
    std::vector<int> UniqueGrainIDs(DomainVol);
    for (int k = ZMin; k <= ZMax; k++) {
        for (int i = XMin; i <= XMax; i++) {
            for (int j = YMin; j <= YMax; j++) {
                UniqueGrainIDs[Counter] = GrainID[k][i][j];
                Counter++;
                if (GrainID[k][i][j] < 0)
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
                    if (GrainID[k][i][j] == ThisGrainID) {
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
                        int YMax, int ZMin, int ZMax, int ***GrainID) {

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
                GIDAllVals_ThisLayer[Counter] = GrainID[k][i][j];
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
                std::string FName3 = BaseFileName + "_GrainWidthDistribution.csv";
                std::ofstream Grainplot3;
                std::cout << "Printing file " << FName3 << " of grain width distribution (in microns) at Z = " << ZMax
                          << std::endl;
                Grainplot3.open(FName3);
                for (int n = 0; n < GrainsThisLayer; n++) {
                    int ThisGrainID = GIDVals_ThisLayer[n];
                    int TempTopX = XMin;
                    int TempTopY = YMin;
                    int TempBottomX = XMax;
                    int TempBottomY = YMax;
                    for (int i = XMin; i <= XMax; i++) {
                        for (int j = YMin; j <= YMax; j++) {
                            if (GrainID[k][i][j] == ThisGrainID) {
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
                }
                Grainplot3.close();
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
void PrintPoleFigureData(bool *AnalysisTypes, std::string BaseFileName, int NumberOfOrientations, float **, int XMin,
                         int XMax, int YMin, int YMax, int ZMin, int ZMax, int ***GrainID, int ***Melted) {

    if (AnalysisTypes[7]) {
        //        writing pyEBSD print routines in progress
        //        std::string FNamePYEBSD = BaseFileName + "_pyEBSDOrientations.csv";
        //        std::ofstream PyEBSDPlot;
        //        PyEBSDPlot.open(FNamePYEBSD);
        //        int RepX = XMax-XMin;
        //        for (int k=ZMin; k<=ZMax; k++) {
        //            for (int i=XMin; i<=XMax; i++) {
        //                int XMultiplier = (k-ZMin)*RepX; // used to offset each XY region in space
        //                for (int j=YMin; j<=YMax; j++) {
        //                    int GOVal = (abs(GrainID[k][i][j]) - 1) % NumberOfOrientations;
        //                    PyEBSDPlot << XMultiplier*(i-XMin) << "," << j-YMin << " " << GrainEulerAngles[GOVal][0]
        //                    << " " << GrainEulerAngles[GOVal][1] << " " << GrainEulerAngles[GOVal][2] << std::endl;
        //                }
        //            }
        //        }
        //        PyEBSDPlot.close();

        // Histogram of orientations for texture determination
        std::ofstream MTEXPlot;
        std::string FNameM = BaseFileName + "_MTEXOrientations.csv";
        MTEXPlot.open(FNameM);
        int *GOHistogram = new int[NumberOfOrientations];
        for (int i = 0; i < NumberOfOrientations; i++) {
            GOHistogram[i] = 0;
        }
        // frequency data on grain ids
        for (int k = ZMin; k <= ZMax; k++) {
            for (int j = YMin; j <= YMax; j++) {
                for (int i = XMin; i <= XMax; i++) {
                    if (Melted[k][i][j]) {
                        int GOVal = (abs(GrainID[k][i][j]) - 1) % NumberOfOrientations;
                        GOHistogram[GOVal]++;
                    }
                }
            }
        }
        for (int i = 0; i < NumberOfOrientations; i++) {
            MTEXPlot << GOHistogram[i] << std::endl;
        }
        MTEXPlot.close();
    }
}

//*****************************************************************************/
void PrintInversePoleFigureCrossSections(int NumberOfCrossSections, std::string BaseFileName,
                                         std::vector<int> CrossSectionPlane, std::vector<int> CrossSectionLocation,
                                         int nx, int ny, int nz, int NumberOfOrientations, int ***GrainID, float **) {

    // Loop over each cross-section specified in the file "AnalysisOutputs.txt"
    for (int n = 0; n < NumberOfCrossSections; n++) {
        // Print data for pyEBSD/MTEX
        std::string CSType;
        if (CrossSectionPlane[n] == 0)
            CSType = "XZ";
        else if (CrossSectionPlane[n] == 1)
            CSType = "YZ";
        else
            CSType = "XY";
        //        pyEBSD print routines in progress
        //        std::string FNameCS1 = BaseFileName + std::to_string(n) + " " + CSType + "pyEBSDCrossSection.ang";
        //        std::cout << "Cross-section number " << n+1 << " being printed to file " << FNameCS1 << " for pyEBSD"
        //        << std::endl; std::ofstream GrainplotCS1; GrainplotCS1.open(FNameCS1); GrainplotCS1.close();

        std::string FNameCS2 = BaseFileName + "-" + std::to_string(n) + "_" + CSType + "MTEXCrossSection.txt";
        std::cout << "Cross-section number " << n + 1 << " being printed to file " << FNameCS2 << " for MTEX"
                  << std::endl;
        std::ofstream GrainplotCS2;
        GrainplotCS2.open(FNameCS2);
        if (CSType == "XZ") {
            for (int i = 0; i < nx; i++) {
                for (int k = 1; k < nz - 1; k++) {
                    int GOVal = (abs(GrainID[k][i][CrossSectionLocation[n]]) - 1) % NumberOfOrientations;
                    GrainplotCS2 << i << "," << k << "," << GOVal << std::endl;
                }
            }
        }
        else if (CSType == "YZ") {
            for (int j = 0; j < ny; j++) {
                for (int k = 1; k < nz - 1; k++) {
                    int GOVal = (abs(GrainID[k][CrossSectionLocation[n]][j]) - 1) % NumberOfOrientations;
                    GrainplotCS2 << j << "," << k << "," << GOVal << std::endl;
                }
            }
        }
        else {
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    int GOVal = (abs(GrainID[CrossSectionLocation[n]][i][j]) - 1) % NumberOfOrientations;
                    GrainplotCS2 << i << "," << j << "," << GOVal << std::endl;
                }
            }
        }
        GrainplotCS2.close();
    }
}

//*****************************************************************************/
void PrintExaConstitRVEData(int NumberOfRVEs, std::string BaseFileName, int, int, int, double deltax,
                            int ***GrainID, std::vector<int> XLow_RVE, std::vector<int> XHigh_RVE,
                            std::vector<int> YLow_RVE, std::vector<int> YHigh_RVE, std::vector<int> ZLow_RVE,
                            std::vector<int> ZHigh_RVE) {

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
        for (int k = ZLow_RVE[n]; k <= ZHigh_RVE[n]; k++) {
            for (int i = XLow_RVE[n]; i <= XHigh_RVE[n]; i++) {
                for (int j = YLow_RVE[n]; j <= YHigh_RVE[n]; j++) {
                    GrainplotE << i << "," << j << "," << k << "," << GrainID[k][i][j] << std::endl;
                }
            }
        }
        GrainplotE.close();
    }
}
