// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT
#include "header.h"
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

//*****************************************************************************/
void PrintMisorientationData(bool* AnalysisTypes, std::string BaseFileName, int XMin, int XMax, int YMin, int YMax, int nz, int*** Melted, float*** GrainUnitVector, int*** GrainID, int NumberOfOrientations, int &NumberOfMeltedCells) {

    // Frequency of misorientations in the selected region
    std::ofstream MisrientationPlot;
    std::string FNameM = BaseFileName + "_MisorientationFrequency.txt";
    if (AnalysisTypes[0]) {
        MisrientationPlot.open(FNameM);
        std::cout << "Printing file of grain misorientations relative to the +Z direction, for all cells simulated by ExaCA (not including substrate)" << std::endl;
    }
    NumberOfMeltedCells = 0;
    long double MisorientationSum = 0.0;
    
    for (int k = 0; k < nz; k++) {
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
                    if (AnalysisTypes[0]) MisrientationPlot << AngleZmin << std::endl;
                    MisorientationSum += AngleZmin;
                    NumberOfMeltedCells++;
                }
            }
        }
    }
    if (AnalysisTypes[0]) MisrientationPlot.close();
    std::cout << "Within the representative region consisting of " << (XMax-XMin+1)*(YMax-YMin+1)*(nz) << " cells, " << NumberOfMeltedCells << " underwent melting and:" << std::endl;
    std::cout << "the mean misorientation relative to the +Z direction is " << MisorientationSum/((double)(NumberOfMeltedCells)) << " degrees" << std::endl;
        
}

//*****************************************************************************/
void PrintRVEData(int NumberOfRVEs, std::string BaseFileName, int nx, int ny, int nz, double deltax, int*** GrainID, std::vector <int> CenterX_RVE, std::vector <int> CenterY_RVE, std::vector <int> CenterZ_RVE, std::vector <int> Size_RVE, std::vector <bool> ExaConstitPrint, std::vector <bool> pyEBSDPrint) {

    // Loop over each RVE specified in the file "AnalysisOutputs.txt"
    for (int n=0; n<NumberOfRVEs; n++) {
        int xbound_low = CenterX_RVE[n] - Size_RVE[n]/2;
        int ybound_low = CenterY_RVE[n] - Size_RVE[n]/2;
        int zbound_low = CenterZ_RVE[n] - Size_RVE[n]/2;
        int xbound_high = Size_RVE[n] + xbound_low - 1;
        int ybound_high = Size_RVE[n] + ybound_low - 1;
        int zbound_high = Size_RVE[n] + zbound_low - 1;
        if ((xbound_low < 0)||(ybound_low < 0)||(zbound_low < 0)||(xbound_high > nx-1)||(ybound_high > ny-1)||(zbound_high > nz-1)) {
            std::cout << "WARNING: RVE number " << n << " cannot be printed, the selected center location and domain size go out of the bounds of the ExaCA microstructure data" << std::endl;
        }
        else {
            if (!(ExaConstitPrint[n])&&(!(pyEBSDPrint[n]))) {
                std::cout << "WARNING: No print options specified for RVE number " << n << std::endl;
            }
            if (ExaConstitPrint[n]) {
                // Print data for ExaConstit
                std::string FNameE = BaseFileName + std::to_string(n) + "_ExaConstit.csv";
                std::cout << "RVE number " << n+1 << " with center location (" << CenterX_RVE[n] << "," << CenterY_RVE[n] << "," << CenterZ_RVE[n] << "), and size " << Size_RVE[n] << " being printed to file for ExaConstit" << std::endl;
                std::ofstream GrainplotE;
                GrainplotE.open(FNameE);
                GrainplotE << "Coordinates are in CA units, 1 cell = " << deltax << " microns. Data is cell-centered. Origin at " << xbound_low << "," << ybound_low << "," << zbound_low << " , domain size is " << xbound_high-xbound_low << " by " << ybound_high-ybound_low << " by " << zbound_high-zbound_low << " cells" << std::endl;
                GrainplotE << "X coord, Y coord, Z coord, Grain ID" << std::endl;
                for (int k = zbound_low; k <=zbound_high; k++) {
                    for (int i = xbound_low; i <= xbound_high; i++) {
                        for (int j = ybound_low; j <= ybound_high; j++) {
                            GrainplotE << i << "," << j<< "," << k << "," << GrainID[k][i][j] << std::endl;
                        }
                    }
                }
                GrainplotE.close();
            }
            if (pyEBSDPrint[n]) {
                // work in progress - print data in a way to match .ang file format from pyEBSD example problems
                // but spread the volume out in a single cross-section
            }
        }
    }

}

//*****************************************************************************/

void PrintVolumeData(bool* AnalysisTypes, int XMin, int XMax, int YMin, int YMax, int nz, int NumberOfMeltedCells, int*** Melted, int*** GrainID) {

    int NucleatedGrainCells = 0;
    for (int k = 0; k < nz; k++) {
        for (int i = XMin; i <= XMax; i++) {
            for (int j = YMin; j <= YMax; j++) {
                // Only take data from cells in the representative area that underwent melting
                if (Melted[k][i][j] == 1) {
                    if (GrainID[k][i][j] < 0) NucleatedGrainCells++;
                }
            }
        }
    }
    std::cout << "The volume fraction consisting of nucleated grains is " << ((double)(NucleatedGrainCells))/((double)(NumberOfMeltedCells)) << std::endl;

}

//*****************************************************************************/
void PrintGrainAreaData(bool* AnalysisTypes, std::string BaseFileName, double deltax, int XMin, int XMax, int YMin, int YMax, int nz, int*** GrainID) {

    std::string FName1 = BaseFileName + "_GrainAreas.csv";
    std::string FName2 = BaseFileName + "_WeightedGrainAreas.csv";
    std::ofstream Grainplot1, Grainplot2;
    if (AnalysisTypes[2]) {
        std::cout << "Printing file of grain area values (in square microns) for all Z coordinates" << std::endl;
        Grainplot1.open(FName1);
    }
    if (AnalysisTypes[3]) {
        std::cout << "Printing file of weighted grain area values (in square microns) for every 5th Z coordinate" << std::endl;
        Grainplot2.open(FName2);
    }

    for (int k = 0; k < nz; k++) {
        std::vector<int> GIDVals_ThisLayer;
        for (int i = XMin; i <= XMax; i++) {
            for (int j = YMin; j <= YMax; j++) {
                GIDVals_ThisLayer.push_back(GrainID[k][i][j]);
            }
        }
        std::vector<int>::iterator it;
        sort(GIDVals_ThisLayer.begin(), GIDVals_ThisLayer.end());
        it = std::unique(GIDVals_ThisLayer.begin(), GIDVals_ThisLayer.end());
        int CellsThisLayer = GIDVals_ThisLayer.size();
        GIDVals_ThisLayer.resize(std::distance(GIDVals_ThisLayer.begin(), it));
        int GrainsThisLayer = GIDVals_ThisLayer.size();
        double MeanGrainAreaThisLayer = (double)(CellsThisLayer) / (double)(GrainsThisLayer);
        if ((AnalysisTypes[2])&&(k % 5 == 0)) {
            std::vector<int> GIDAllVals_ThisLayer;
            GIDAllVals_ThisLayer = GIDVals_ThisLayer;
            long int AreaXArea = 0;
            for (int l = 0; l < GrainsThisLayer; l++) {
                long int MyGrainArea = 0;
                for (int ll = 0; ll < CellsThisLayer; ll++) {
                    if (GIDVals_ThisLayer[l] == GIDAllVals_ThisLayer[ll])
                        MyGrainArea++;
                }
                AreaXArea += MyGrainArea * MyGrainArea;
            }
            double WeightedArea = ((double)(AreaXArea) / (double)((XMax - XMin + 1) * (YMax - YMin + 1)));
            Grainplot2 << WeightedArea * deltax * deltax / pow(10, -12) << std::endl;
        }
        Grainplot1 << MeanGrainAreaThisLayer * deltax * deltax / pow(10, -12) << std::endl;
    }
    if (AnalysisTypes[2])  Grainplot1.close();
    if (AnalysisTypes[3])  Grainplot2.close();
    
}
