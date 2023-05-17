// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "GArepresentativeregion.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
// Tests of representative region structure
//---------------------------------------------------------------------------//

void writeTestArea() {
    // Write out example data
    std::ofstream TestData;
    TestData.open("TestArea.json");
    TestData << "{" << std::endl;
    TestData << "   \"Regions\": {" << std::endl;
    TestData << "       \"RepresentativeArea\": {" << std::endl;
    TestData << "          \"units\": \"Meters\"," << std::endl;
    TestData << "          \"zBounds\": [0.000005,0.000005]," << std::endl;
    TestData << "          \"printPoleFigureData\": true," << std::endl;
    TestData << "          \"printInversePoleFigureData\": true," << std::endl;
    TestData << "          \"printStats\": [\"GrainTypeFractions\"]," << std::endl;
    TestData << "          \"printPerGrainStats\": [\"IPFZ-RGB\", \"Size\"]" << std::endl;
    TestData << "      }" << std::endl;
    TestData << "   }" << std::endl;
    TestData << "}" << std::endl;
    TestData.close();
}

void writeTestVolume() {
    // Write out example data
    std::ofstream TestData;
    TestData.open("TestVolume.json");
    TestData << "{" << std::endl;
    TestData << "   \"Regions\": {" << std::endl;
    TestData << "       \"RepresentativeVolume\": {" << std::endl;
    TestData << "          \"units\": \"Cells\"," << std::endl;
    TestData << "          \"xBounds\": [0, 4]," << std::endl;
    TestData << "          \"yBounds\": [1, 10]," << std::endl;
    TestData << "          \"zBounds\": [0, 9]," << std::endl;
    TestData << "          \"printExaConstit\": false," << std::endl;
    TestData << "          \"printPoleFigureData\": true," << std::endl;
    TestData << "          \"printStats\": [\"GrainTypeFractions\", \"Misorientation\", \"Size\", "
                "\"BuildTransAspectRatio\", \"XExtent\", \"YExtent\", \"ZExtent\"],"
             << std::endl;
    TestData << "          \"printPerGrainStats\": [\"IPFZ-RGB\", \"Size\"]," << std::endl;
    TestData << "          \"printLayerwiseData\": [\"MeanGrainArea\", \"MeanWeightedGrainArea\"]" << std::endl;
    TestData << "      }" << std::endl;
    TestData << "   }" << std::endl;
    TestData << "}" << std::endl;
    TestData.close();
}

void testConstructRepresentativeRegion_Volume() {
    // Write out example data
    writeTestVolume();

    // Read in test data
    std::ifstream AnalysisDataStream("TestVolume.json");
    nlohmann::json AnalysisData = nlohmann::json::parse(AnalysisDataStream);
    nlohmann::json RegionData = AnalysisData["Regions"]["RepresentativeVolume"];
    // Domain bounds
    const int nx = 20;
    const int ny = 30;
    const int nz = 40;
    double deltax = 1.25 * pow(10, -6);
    double XMin = -1 * pow(10, -6);
    double YMin = 0;
    double ZMin = 0;
    double XMax = XMin + (nx - 1) * deltax;
    double YMax = YMin + (ny - 1) * deltax;
    double ZMax = ZMin + (nz - 1) * deltax;
    std::vector<double> XYZBounds = {XMin, YMin, ZMin, XMax, YMax, ZMax};
    // View for storing grain ID data
    ViewI3D_H GrainID(Kokkos::ViewAllocateWithoutInitializing("GrainID"), nz, nx, ny);

    // Construct region
    RepresentativeRegion representativeRegion(RegionData, nx, ny, nz, deltax, XYZBounds, GrainID);

    // Check results
    EXPECT_TRUE(representativeRegion.regionType == "volume");
    EXPECT_TRUE(representativeRegion.regionOrientation == "XYZ");
    EXPECT_DOUBLE_EQ(representativeRegion.xBounds_Meters[0], XMin);
    EXPECT_DOUBLE_EQ(representativeRegion.yBounds_Meters[0], deltax);
    EXPECT_DOUBLE_EQ(representativeRegion.zBounds_Meters[0], ZMin);
    EXPECT_DOUBLE_EQ(representativeRegion.xBounds_Meters[1], XMin + deltax * 4);
    EXPECT_DOUBLE_EQ(representativeRegion.yBounds_Meters[1], YMin + deltax * 10);
    EXPECT_DOUBLE_EQ(representativeRegion.zBounds_Meters[1], ZMin + deltax * 9);
    double ExpectedRegionSize_Meters =
        (representativeRegion.xBounds_Meters[1] - representativeRegion.xBounds_Meters[0] + deltax) *
        (representativeRegion.yBounds_Meters[1] - representativeRegion.yBounds_Meters[0] + deltax) *
        (representativeRegion.zBounds_Meters[1] - representativeRegion.zBounds_Meters[0] + deltax);
    double ExpectedRegionSize_Microns = ExpectedRegionSize_Meters * pow(10, 18);
    EXPECT_DOUBLE_EQ(representativeRegion.regionSize_Microns, ExpectedRegionSize_Microns);
    EXPECT_DOUBLE_EQ(representativeRegion.regionSize_Meters, ExpectedRegionSize_Meters);
    EXPECT_EQ(representativeRegion.xBounds_Cells[0], 0);
    EXPECT_EQ(representativeRegion.yBounds_Cells[0], 1);
    EXPECT_EQ(representativeRegion.zBounds_Cells[0], 0);
    EXPECT_EQ(representativeRegion.xBounds_Cells[1], 4);
    EXPECT_EQ(representativeRegion.yBounds_Cells[1], 10);
    EXPECT_EQ(representativeRegion.zBounds_Cells[1], 9);
    int ExpectedRegionSize_Cells = 500;
    EXPECT_EQ(representativeRegion.regionSize_Cells, ExpectedRegionSize_Cells);
    int NumAnalysisOptions_Stats = representativeRegion.AnalysisOptions_StatsYN.size();
    for (int i = 0; i < NumAnalysisOptions_Stats; i++)
        EXPECT_TRUE(representativeRegion.AnalysisOptions_StatsYN[i]);
    int NumAnalysisOptions_PerGrainStats = representativeRegion.AnalysisOptions_PerGrainStatsYN.size();
    for (int i = 0; i < NumAnalysisOptions_PerGrainStats; i++) {
        if ((i == 1) || (i == 6))
            EXPECT_TRUE(representativeRegion.AnalysisOptions_PerGrainStatsYN[i]);
        else
            EXPECT_FALSE(representativeRegion.AnalysisOptions_PerGrainStatsYN[i]);
    }
    int NumAnalysisOptions_LayerwiseStats = representativeRegion.AnalysisOptions_LayerwiseStatsYN.size();
    for (int i = 0; i < NumAnalysisOptions_LayerwiseStats; i++)
        EXPECT_TRUE(representativeRegion.AnalysisOptions_LayerwiseStatsYN[i]);
    EXPECT_TRUE(representativeRegion.PrintPerGrainStatsYN);
    EXPECT_TRUE(representativeRegion.PrintPoleFigureYN);
    EXPECT_FALSE(representativeRegion.PrintExaConstitYN);
    EXPECT_FALSE(representativeRegion.PrintInversePoleFigureMapYN);
}

// Initialize a representative area from json format
void testConstructRepresentativeRegion_Area() {

    // Write out example data
    writeTestArea();

    // Read in test data
    std::ifstream AnalysisDataStream("TestArea.json");
    nlohmann::json AnalysisData = nlohmann::json::parse(AnalysisDataStream);
    nlohmann::json RegionData = AnalysisData["Regions"]["RepresentativeArea"];

    // Domain bounds
    const int nx = 10;
    const int ny = 20;
    const int nz = 20;
    double deltax = 1.25 * pow(10, -6);
    double XMin = 0;
    double YMin = 0;
    double ZMin = 0;
    double XMax = XMin + (nx - 1) * deltax;
    double YMax = YMin + (ny - 1) * deltax;
    double ZMax = ZMin + (nz - 1) * deltax;
    std::vector<double> XYZBounds = {XMin, YMin, ZMin, XMax, YMax, ZMax};
    // View for storing grain ID data
    ViewI3D_H GrainID(Kokkos::ViewAllocateWithoutInitializing("GrainID"), nz, nx, ny);

    // Construct region
    RepresentativeRegion representativeRegion(RegionData, nx, ny, nz, deltax, XYZBounds, GrainID);

    // Check results
    EXPECT_TRUE(representativeRegion.regionType == "area");
    EXPECT_TRUE(representativeRegion.regionOrientation == "XY");
    EXPECT_DOUBLE_EQ(representativeRegion.xBounds_Meters[0], XMin);
    EXPECT_DOUBLE_EQ(representativeRegion.yBounds_Meters[0], YMin);
    EXPECT_DOUBLE_EQ(representativeRegion.zBounds_Meters[0], 0.000005);
    EXPECT_DOUBLE_EQ(representativeRegion.xBounds_Meters[1], XMax);
    EXPECT_DOUBLE_EQ(representativeRegion.yBounds_Meters[1], YMax);
    EXPECT_DOUBLE_EQ(representativeRegion.zBounds_Meters[1], 0.000005);
    EXPECT_EQ(representativeRegion.xBounds_Cells[0], 0);
    EXPECT_EQ(representativeRegion.yBounds_Cells[0], 0);
    EXPECT_EQ(representativeRegion.zBounds_Cells[0], 4);
    EXPECT_EQ(representativeRegion.xBounds_Cells[1], 9);
    EXPECT_EQ(representativeRegion.yBounds_Cells[1], 19);
    EXPECT_EQ(representativeRegion.zBounds_Cells[1], 4);
    int ExpectedRegionSize_Cells = 200;
    EXPECT_EQ(representativeRegion.regionSize_Cells, ExpectedRegionSize_Cells);
    double ExpectedRegionSize_Meters = ExpectedRegionSize_Cells * pow(deltax, 2);
    double ExpectedRegionSize_Microns = ExpectedRegionSize_Meters * pow(10, 12);
    EXPECT_DOUBLE_EQ(representativeRegion.regionSize_Microns, ExpectedRegionSize_Microns);
    EXPECT_DOUBLE_EQ(representativeRegion.regionSize_Meters, ExpectedRegionSize_Meters);
    EXPECT_TRUE(representativeRegion.AnalysisOptions_StatsYN[0]);
    int NumAnalysisOptions_Stats = representativeRegion.AnalysisOptions_StatsYN.size();
    for (int i = 1; i < NumAnalysisOptions_Stats; i++)
        EXPECT_FALSE(representativeRegion.AnalysisOptions_StatsYN[i]);
    int NumAnalysisOptions_PerGrainStats = representativeRegion.AnalysisOptions_PerGrainStatsYN.size();
    for (int i = 0; i < NumAnalysisOptions_PerGrainStats; i++) {
        if ((i == 1) || (i == 6))
            EXPECT_TRUE(representativeRegion.AnalysisOptions_PerGrainStatsYN[i]);
        else
            EXPECT_FALSE(representativeRegion.AnalysisOptions_PerGrainStatsYN[i]);
    }
    int NumAnalysisOptions_LayerwiseStats = representativeRegion.AnalysisOptions_LayerwiseStatsYN.size();
    for (int i = 0; i < NumAnalysisOptions_LayerwiseStats; i++)
        EXPECT_FALSE(representativeRegion.AnalysisOptions_LayerwiseStatsYN[i]);
    EXPECT_TRUE(representativeRegion.PrintPerGrainStatsYN);
    EXPECT_TRUE(representativeRegion.PrintPoleFigureYN);
    EXPECT_FALSE(representativeRegion.PrintExaConstitYN);
    EXPECT_TRUE(representativeRegion.PrintInversePoleFigureMapYN);
}

// For a 3D domain, obtain a vector of the grain IDs that appear in a representative portion, and the associated grain
// statistics from the representative region
void testCollectGrainStats() {

    // Write out example data
    writeTestVolume();

    const int nx = 5;  // Representative region spans 0-4 (entire domain in X)
    const int ny = 11; // Representative region spans 1-10 (one less than domain size in Y)
    const int nz = 12; // Representative region spans 0-9 (two less than domain size in Z)
    const double deltax = 1.25 * pow(10, -6);
    std::vector<double> XYZBounds = {0.0, 0.0, 0.0, nx * deltax, ny * deltax, nz * deltax};

    // View for storing grain ID data
    ViewI3D_H GrainID(Kokkos::ViewAllocateWithoutInitializing("GrainID"), nz, nx, ny);
    // Assign grain ID using the Z coordinate
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                GrainID(k, i, j) = k;
            }
        }
    }

    // Representative region creation
    std::ifstream AnalysisDataStream("TestVolume.json");
    nlohmann::json AnalysisData = nlohmann::json::parse(AnalysisDataStream);
    nlohmann::json RegionsData = AnalysisData["Regions"]["RepresentativeVolume"];
    RepresentativeRegion representativeRegion(RegionsData, nx, ny, nz, deltax, XYZBounds, GrainID);

    // Check results of grain ID vector (50 of each grain ID should exist)
    EXPECT_EQ(representativeRegion.regionSize_Cells, nx * (ny - 1) * (nz - 2));
    for (int n = 0; n < representativeRegion.regionSize_Cells; n++) {
        EXPECT_EQ(representativeRegion.GrainIDVector[n], n / 50);
    }

    // Check number of grains and unique grain IDs
    EXPECT_EQ(representativeRegion.NumberOfGrains, 10);
    for (int n = 0; n < representativeRegion.NumberOfGrains; n++) {
        EXPECT_EQ(representativeRegion.UniqueGrainIDVector[n], n);
    }

    // Check against expected grain sizes, in cubic microns (each grain spans 50 cells)
    for (int n = 0; n < representativeRegion.NumberOfGrains; n++) {
        EXPECT_FLOAT_EQ(representativeRegion.GrainSizeVector_Microns[n], 50 * pow(deltax, 3) * pow(10, 18));
    }

    // Obtain the grain extents
    // Extent of each grain in X
    std::vector<float> GrainExtentX(representativeRegion.NumberOfGrains);
    representativeRegion.calcGrainExtent(GrainExtentX, GrainID, "X", deltax);
    // Extent of each grain in Y
    std::vector<float> GrainExtentY(representativeRegion.NumberOfGrains);
    representativeRegion.calcGrainExtent(GrainExtentY, GrainID, "Y", deltax);
    // Extent of each grain in Z
    std::vector<float> GrainExtentZ(representativeRegion.NumberOfGrains);
    representativeRegion.calcGrainExtent(GrainExtentZ, GrainID, "Z", deltax);
    // Check grain extents
    for (int n = 0; n < representativeRegion.NumberOfGrains; n++) {
        // Expected value should be in microns
        EXPECT_FLOAT_EQ(GrainExtentX[n], nx * deltax * pow(10, 6));
        EXPECT_FLOAT_EQ(GrainExtentY[n], (ny - 1) * deltax * pow(10, 6));
        EXPECT_FLOAT_EQ(GrainExtentZ[n], deltax * pow(10, 6));
    }
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, representative_region) {
    testConstructRepresentativeRegion_Volume();
    testConstructRepresentativeRegion_Area();
    testCollectGrainStats();
}
} // end namespace Test
