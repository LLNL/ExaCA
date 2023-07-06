// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAcelldata.hpp"
#include "CAfunctions.hpp"
#include "CAinitialize.hpp"
#include "CAnucleation.hpp"
#include "CAparsefiles.hpp"
#include "CAtypes.hpp"

#include <gtest/gtest.h>

#include "mpi.h"

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
// file_read_tests
//---------------------------------------------------------------------------//
void WriteTestData(std::string InputFilename, int PrintVersion) {

    std::ofstream TestDataFile;
    TestDataFile.open(InputFilename);
    // Write required inputs to all files
    TestDataFile << "{" << std::endl;
    TestDataFile << "   \"SimulationType\": \"R\"," << std::endl;
    TestDataFile << "   \"MaterialFileName\": \"Inconel625.json\"," << std::endl;
    TestDataFile << "   \"GrainOrientationFile\": \"GrainOrientationVectors.csv\"," << std::endl;
    TestDataFile << "   \"RandomSeed\": 2," << std::endl;
    TestDataFile << "   \"Domain\": {" << std::endl;
    TestDataFile << "      \"CellSize\": 1," << std::endl;
    TestDataFile << "      \"TimeStep\": 1.5," << std::endl;
    TestDataFile << "      \"NumberOfLayers\": 2," << std::endl;
    TestDataFile << "      \"LayerOffset\": 1" << std::endl;
    TestDataFile << "   }," << std::endl;
    TestDataFile << "   \"Nucleation\": {" << std::endl;
    TestDataFile << "      \"Density\": 10," << std::endl;
    TestDataFile << "      \"MeanUndercooling\": 5," << std::endl;
    TestDataFile << "      \"StDev\": 0.5" << std::endl;
    TestDataFile << "   }," << std::endl;
    TestDataFile << "   \"TemperatureData\": {" << std::endl;
    TestDataFile << "      \"TemperatureFiles\": [\".//1DummyTemperature.txt\",\".//2DummyTemperature.txt\"]"
                 << std::endl;
    TestDataFile << "   }," << std::endl;
    TestDataFile << "   \"Substrate\": {" << std::endl;
    TestDataFile << "      \"SubstrateFilename\": \"DummySubstrate.txt\"," << std::endl;
    TestDataFile << "      \"PowderDensity\": 1000," << std::endl;
    TestDataFile << "      \"PowderFirstLayer\": true" << std::endl;
    TestDataFile << "   }," << std::endl;
    TestDataFile << "   \"Printing\": {" << std::endl;
    TestDataFile << "      \"PathToOutput\": \"ExaCA\"," << std::endl;
    TestDataFile << "      \"OutputFile\": \"Test\"," << std::endl;
    TestDataFile << "      \"PrintBinary\": true," << std::endl;
    // two different permutations of print outputs
    if (PrintVersion == 0) {
        TestDataFile << "      \"PrintExaConstitSize\": 0," << std::endl;
        TestDataFile << "      \"PrintFieldsInit\": "
                        "[\"UndercoolingChange\",\"GrainID\",\"LayerID\",\"MeltTimeStep\",\"CritTimeStep\"],"
                     << std::endl;
        TestDataFile
            << "      \"PrintFieldsFinal\": [\"UndercoolingCurrent\",\"GrainMisorientation\",\"GrainID\",\"LayerID\"],"
            << std::endl;
        TestDataFile << "      \"PrintIntermediateOutput\": {" << std::endl;
        TestDataFile << "          \"Frequency\": 300," << std::endl;
        TestDataFile << "          \"PrintIdleFrames\": true" << std::endl;
        TestDataFile << "       }" << std::endl;
    }
    else if (PrintVersion == 1) {
        TestDataFile << "      \"PrintExaConstitSize\": 500," << std::endl;
        TestDataFile << "      \"PrintFieldsInit\": []," << std::endl;
        TestDataFile << "      \"PrintFieldsFinal\": [\"GrainMisorientation\",\"GrainID\",\"LayerID\"]" << std::endl;
    }
    TestDataFile << "   }" << std::endl;
    TestDataFile << "}" << std::endl;
    TestDataFile.close();
}

void testInputReadFromFile(int PrintVersion) {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // Three input files - one of each type
    // Inp_DirSolidification.txt and Inp_SpotMelt.txt were installed from the examples directory
    // Since no temperature files exist in the repo, and there is no ability to write temperature files to a different
    // directory ( would need examples/Temperatures) using the C++11 standard, dummy input files are written and parsed
    // to test an example problem that uses temperature data from a file.
    std::vector<std::string> InputFilenames = {"Inp_DirSolidification", "Inp_SpotMelt", "Inp_TemperatureTest"};
    for (int n = 0; n < 3; n++) {
        InputFilenames[n] += ".json";
    }
    std::vector<std::string> TemperatureFNames = {"1DummyTemperature.txt", "2DummyTemperature.txt"};

    // On rank 0, write dummy input files for using read temperature data (InputFilenames[2] and [3])
    if (id == 0) {
        WriteTestData(InputFilenames[2], PrintVersion);

        // Create test temperature files "1DummyTemperature.txt", "2DummyTemperature.txt"
        std::ofstream TestTemp1, TestTemp2;
        TestTemp1.open(TemperatureFNames[0]);
        TestTemp2.open(TemperatureFNames[1]);
        TestTemp1 << "X" << std::endl;
        TestTemp2 << "X" << std::endl;
        TestTemp1.close();
        TestTemp2.close();

        // Create test substrate file "DummySubstrate.txt"
        std::ofstream TestSub;
        TestSub.open("DummySubstrate.txt");
        TestSub << "X" << std::endl;
        TestSub.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Read and parse each input file
    for (auto FileName : InputFilenames) {
        int TempFilesInSeries, NumberOfLayers, LayerHeight, nx, ny, nz, NSpotsX, NSpotsY, SpotOffset, SpotRadius;
        float SubstrateGrainSpacing;
        double deltax, NMax, dTN, dTsigma, HT_deltax, deltat, G, R, FractSurfaceSitesActive, RNGSeed,
            PowderActiveFraction;
        bool BaseplateThroughPowder, LayerwiseTempRead, UseSubstrateFile, PowderFirstLayer;
        std::string SimulationType, GrainOrientationFile, temppath, tempfile, SubstrateFileName, MaterialFileName;
        std::vector<std::string> temp_paths;
        std::cout << "Reading " << FileName << std::endl;
        // Data printing structure - contains print options (false by default) and functions
        Print print(np);
        InputReadFromFile(id, FileName, SimulationType, deltax, NMax, dTN, dTsigma, GrainOrientationFile,
                          TempFilesInSeries, temp_paths, HT_deltax, deltat, NumberOfLayers, LayerHeight,
                          MaterialFileName, SubstrateFileName, SubstrateGrainSpacing, UseSubstrateFile, G, R, nx, ny,
                          nz, FractSurfaceSitesActive, NSpotsX, NSpotsY, SpotOffset, SpotRadius, RNGSeed,
                          BaseplateThroughPowder, PowderActiveFraction, LayerwiseTempRead, PowderFirstLayer, print);
        InterfacialResponseFunction irf(0, MaterialFileName, deltat, deltax);
        MPI_Barrier(MPI_COMM_WORLD);

        // Check the results
        // The existence of the specified orientation, substrate, and temperature filenames was already checked within
        // InputReadFromFile
        // These should be the same for all 3 test problems
        EXPECT_DOUBLE_EQ(deltax, 1.0 * pow(10, -6));
        EXPECT_DOUBLE_EQ(NMax, 1.0 * pow(10, 13));
        EXPECT_DOUBLE_EQ(dTN, 5.0);
        EXPECT_DOUBLE_EQ(dTsigma, 0.5);
        EXPECT_DOUBLE_EQ(irf.A, -0.00000010302 * deltat / deltax);
        EXPECT_DOUBLE_EQ(irf.B, 0.00010533 * deltat / deltax);
        EXPECT_DOUBLE_EQ(irf.C, 0.0022196 * deltat / deltax);
        EXPECT_DOUBLE_EQ(irf.D, 0);
        EXPECT_DOUBLE_EQ(irf.FreezingRange, 210);

        // These are different for all 3 test problems
        if (FileName == InputFilenames[0]) {
            EXPECT_TRUE(print.PrintTimeSeries);
            EXPECT_EQ(print.TimeSeriesInc, 5250);
            EXPECT_FALSE(print.PrintIdleTimeSeriesFrames);
            EXPECT_DOUBLE_EQ(G, 500000.0);
            EXPECT_DOUBLE_EQ(R, 300000.0);
            // compare with float to avoid floating point error with irrational number
            float deltat_comp = static_cast<float>(deltat);
            EXPECT_FLOAT_EQ(deltat_comp, pow(10, -6) / 15.0);
            EXPECT_EQ(nx, 200);
            EXPECT_EQ(ny, 200);
            EXPECT_EQ(nz, 200);
            EXPECT_DOUBLE_EQ(FractSurfaceSitesActive, 0.08);
            EXPECT_TRUE(print.BaseFileName == "TestProblemDirS");
            EXPECT_TRUE(print.PrintInitCritTimeStep);
            EXPECT_FALSE(print.PrintInitGrainID);
            EXPECT_FALSE(print.PrintInitLayerID);
            EXPECT_FALSE(print.PrintInitMeltTimeStep);
            EXPECT_FALSE(print.PrintInitUndercoolingChange);
            EXPECT_TRUE(print.PrintFinalGrainID);
            EXPECT_TRUE(print.PrintFinalLayerID);
            EXPECT_TRUE(print.PrintFinalMisorientation);
            EXPECT_FALSE(print.PrintFinalUndercoolingCurrent);
            EXPECT_FALSE(print.PrintFinalMeltTimeStep);
            EXPECT_FALSE(print.PrintFinalCritTimeStep);
            EXPECT_FALSE(print.PrintFinalUndercoolingChange);
            EXPECT_FALSE(print.PrintDefaultRVE);
            EXPECT_DOUBLE_EQ(RNGSeed, 0.0);
            EXPECT_FALSE(print.PrintBinary);
        }
        else if (FileName == InputFilenames[1]) {
            EXPECT_TRUE(print.PrintTimeSeries);
            EXPECT_EQ(print.TimeSeriesInc, 37500);
            EXPECT_TRUE(print.PrintIdleTimeSeriesFrames);
            EXPECT_DOUBLE_EQ(G, 500000.0);
            EXPECT_DOUBLE_EQ(R, 300000.0);
            // compare with float to avoid floating point error with irrational number
            float deltat_comp = static_cast<float>(deltat);
            EXPECT_FLOAT_EQ(deltat_comp, pow(10, -6) / 15.0);
            EXPECT_EQ(NSpotsX, 3);
            EXPECT_EQ(NSpotsY, 2);
            EXPECT_EQ(SpotOffset, 100);
            EXPECT_EQ(SpotRadius, 75);
            EXPECT_EQ(NumberOfLayers, 2);
            EXPECT_EQ(LayerHeight, 20);
            EXPECT_FALSE(UseSubstrateFile);
            EXPECT_FALSE(BaseplateThroughPowder);
            // Option defaults to false
            EXPECT_FALSE(PowderFirstLayer);
            EXPECT_FLOAT_EQ(SubstrateGrainSpacing, 25.0);
            EXPECT_TRUE(print.BaseFileName == "TestProblemSpot");
            EXPECT_TRUE(print.PrintInitCritTimeStep);
            EXPECT_TRUE(print.PrintInitMeltTimeStep);
            EXPECT_FALSE(print.PrintInitGrainID);
            EXPECT_FALSE(print.PrintInitLayerID);
            EXPECT_FALSE(print.PrintInitUndercoolingChange);
            EXPECT_TRUE(print.PrintFinalGrainID);
            EXPECT_TRUE(print.PrintFinalLayerID);
            EXPECT_TRUE(print.PrintFinalMisorientation);
            EXPECT_FALSE(print.PrintFinalUndercoolingCurrent);
            EXPECT_FALSE(print.PrintFinalMeltTimeStep);
            EXPECT_FALSE(print.PrintFinalCritTimeStep);
            EXPECT_FALSE(print.PrintFinalUndercoolingChange);
            EXPECT_FALSE(print.PrintDefaultRVE);
            EXPECT_DOUBLE_EQ(RNGSeed, 0.0);
        }
        else if (FileName == InputFilenames[2]) {
            EXPECT_DOUBLE_EQ(deltat, 1.5 * pow(10, -6));
            EXPECT_EQ(TempFilesInSeries, 2);
            EXPECT_EQ(NumberOfLayers, 2);
            EXPECT_EQ(LayerHeight, 1);
            EXPECT_TRUE(UseSubstrateFile);
            EXPECT_FALSE(LayerwiseTempRead);
            EXPECT_DOUBLE_EQ(PowderActiveFraction, 0.001);
            EXPECT_TRUE(PowderFirstLayer);
            EXPECT_DOUBLE_EQ(HT_deltax, deltax);
            EXPECT_TRUE(print.BaseFileName == "Test");
            EXPECT_TRUE(temp_paths[0] == ".//1DummyTemperature.txt");
            EXPECT_TRUE(temp_paths[1] == ".//2DummyTemperature.txt");
            if (PrintVersion == 0) {
                EXPECT_TRUE(print.PrintInitCritTimeStep);
                EXPECT_TRUE(print.PrintInitMeltTimeStep);
                EXPECT_TRUE(print.PrintInitGrainID);
                EXPECT_TRUE(print.PrintInitLayerID);
                EXPECT_TRUE(print.PrintInitUndercoolingChange);
                EXPECT_TRUE(print.PrintFinalMisorientation);
                EXPECT_TRUE(print.PrintFinalUndercoolingCurrent);
                EXPECT_TRUE(print.PrintFinalLayerID);
                EXPECT_TRUE(print.PrintFinalGrainID);
                EXPECT_TRUE(print.PrintTimeSeries);
                EXPECT_EQ(print.TimeSeriesInc, 200); // value from file divided by deltat
                EXPECT_TRUE(print.PrintIdleTimeSeriesFrames);
                EXPECT_FALSE(print.PrintDefaultRVE);
            }
            else if (PrintVersion == 1) {
                EXPECT_FALSE(print.PrintInitCritTimeStep);
                EXPECT_FALSE(print.PrintInitMeltTimeStep);
                EXPECT_FALSE(print.PrintInitGrainID);
                EXPECT_FALSE(print.PrintInitLayerID);
                EXPECT_FALSE(print.PrintInitUndercoolingChange);
                EXPECT_TRUE(print.PrintFinalMisorientation);
                EXPECT_FALSE(print.PrintFinalUndercoolingCurrent);
                EXPECT_TRUE(print.PrintFinalLayerID);
                EXPECT_TRUE(print.PrintFinalGrainID);
                EXPECT_TRUE(print.PrintDefaultRVE);
                EXPECT_FALSE(print.PrintTimeSeries);
            }
            EXPECT_FALSE(LayerwiseTempRead);
            EXPECT_DOUBLE_EQ(RNGSeed, 2.0);
            EXPECT_TRUE(print.PrintBinary);
        }
    }
}
//---------------------------------------------------------------------------//
// cell_init_tests
//---------------------------------------------------------------------------//
void testCellDataInit_ConstrainedGrowth() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // Create test data
    int nz = 2;
    int ZBound_Low = 0;
    int nzActive = 2;
    int nx = 1;
    double deltax = 1 * pow(10, -6);
    double *ZMaxLayer = new double[1];
    ZMaxLayer[0] = (nz - 1) * deltax;
    // Domain size in Y depends on the number of ranks - each rank has 4 cells in Y
    // Each rank is assigned a different portion of the domain in Y
    int ny = 4 * np;
    int MyYSlices = 4;
    int MyYOffset = 4 * id;
    int LocalActiveDomainSize = nx * MyYSlices * nzActive;
    int LocalDomainSize = nx * MyYSlices * nz;

    double FractSurfaceSitesActive = 0.5; // Each rank will have 2 active cells each, on average
    double RNGSeed = 0.0;
    std::string GrainOrientationFile = checkFileInstalled("GrainOrientationVectors.csv", id);
    int NGrainOrientations = 10000; // Number of grain orientations considered in the simulation
    ViewF GrainUnitVector(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), 9 * NGrainOrientations);
    OrientationInit(id, NGrainOrientations, GrainUnitVector, GrainOrientationFile);

    // Initialize neighbor lists
    NList NeighborX, NeighborY, NeighborZ;
    NeighborListInit(NeighborX, NeighborY, NeighborZ);

    // Initialize views - set initial GrainID values to 0, all CellType values to liquid
    ViewF DiagonalLength(Kokkos::ViewAllocateWithoutInitializing("DiagonalLength"), LocalActiveDomainSize);
    ViewI NumberOfSolidificationEvents(Kokkos::ViewAllocateWithoutInitializing("NumberOfSolidificationEvents"),
                                       LocalActiveDomainSize);
    ViewF DOCenter(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * LocalActiveDomainSize);
    ViewF CritDiagonalLength(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength"), 26 * LocalActiveDomainSize);
    CellData<TEST_MEMSPACE> cellData(LocalDomainSize, LocalActiveDomainSize, nx, MyYSlices, ZBound_Low);
    cellData.init_substrate(id, FractSurfaceSitesActive, MyYSlices, nx, ny, MyYOffset, NeighborX, NeighborY, NeighborZ,
                            GrainUnitVector, NGrainOrientations, DiagonalLength, DOCenter, CritDiagonalLength, RNGSeed);
    // Copy CellType, GrainID views to host to check values
    ViewI_H CellType_AllLayers_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cellData.CellType_AllLayers);
    ViewI_H GrainID_AllLayers_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cellData.GrainID_AllLayers);
    for (int i = 0; i < LocalDomainSize; i++) {
        if (i >= nx * MyYSlices) {
            // Not at bottom surface - should be liquid cells with GrainID still equal to 0
            EXPECT_EQ(GrainID_AllLayers_Host(i), 0);
            EXPECT_EQ(CellType_AllLayers_Host(i), Liquid);
        }
        else {
            // Check that active cells have GrainIDs > 0, and less than 2 * np + 1 (there are 2 * np different positive
            // GrainIDs used for epitaxial grain seeds)
            if (CellType_AllLayers_Host(i) == Active) {
                EXPECT_GT(GrainID_AllLayers_Host(i), 0);
                EXPECT_LT(GrainID_AllLayers_Host(i), 2 * np + 1);
            }
            else {
                // Liquid cells should still have GrainID = 0
                EXPECT_EQ(GrainID_AllLayers_Host(i), 0);
            }
        }
    }
}

void testCellDataInit(bool PowderFirstLayer) {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Create test data - 3 layers, starting on layer 0
    double deltax = 1 * pow(10, -6);
    int nz = 5;
    double ZMinLayer[3];
    double ZMaxLayer[3];

    // First layer: Z = 0 through 2
    int nzActive = 3;
    int ZBound_Low = 0;
    ZMinLayer[0] = 0;
    ZMaxLayer[0] = 2 * deltax;
    double ZMin = ZMinLayer[0];
    // Second layer: Z = 1 through 3
    ZMinLayer[1] = 1 * deltax;
    ZMaxLayer[1] = 3 * deltax;
    // Third layer: Z = 2 through 4
    ZMinLayer[2] = 2 * deltax;
    ZMaxLayer[2] = 4 * deltax;

    int nx = 3;
    // Each rank is assigned a different portion of the domain in Y
    int ny = 3 * np;
    int MyYSlices = 3;
    int MyYOffset = 3 * id;

    // If there is a powder layer, the baseplate should be Z = 0 through 1 w/ powder for the top row of cells, otherwise
    // it should be 0 through 2
    int BaseplateSize;
    int LayerHeight = 1;
    if (PowderFirstLayer)
        BaseplateSize = nx * MyYSlices * (round((ZMaxLayer[0] - ZMin) / deltax) + 1 - LayerHeight);
    else
        BaseplateSize = nx * MyYSlices * (round((ZMaxLayer[0] - ZMin) / deltax) + 1);
    int LocalActiveDomainSize = nx * MyYSlices * nzActive;
    int LocalDomainSize = nx * MyYSlices * nz;
    // There are 45 * np total cells in this domain (nx * ny * nz)
    // Each rank has 45 cells - the bottom 27 cells are assigned baseplate Grain ID values, unless the powder layer of
    // height 1 is used, in which case only the bottom 18 cells should be assigned baseplate Grain ID values. The top
    // cells (Z > 2) are outside the first layer of the domain and are not assigned Grain IDs with the rest of the
    // baseplate. This grain spacing ensures that there will be 1 grain per number of MPI ranks present (larger when
    // powder layer is present as the baseplate will only have a third as many cells)
    double SubstrateGrainSpacing;
    if (PowderFirstLayer)
        SubstrateGrainSpacing = 2.62;
    else
        SubstrateGrainSpacing = 3.0;
    double RNGSeed = 0.0;

    // Used if Powderfirstlayer = true
    double PowderActiveFraction = 1.0;

    // unused views in constructor
    NList NeighborX, NeighborY, NeighborZ;
    NeighborListInit(NeighborX, NeighborY, NeighborZ);
    ViewF GrainUnitVector(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), 0);
    ViewF DiagonalLength(Kokkos::ViewAllocateWithoutInitializing("DiagonalLength"), 0);
    ViewF DOCenter(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 0);
    ViewF CritDiagonalLength(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength"), 0);

    // Create dummy temperature data
    ViewI NumberOfSolidificationEvents(Kokkos::ViewAllocateWithoutInitializing("NumberOfSolidificationEvents"),
                                       LocalActiveDomainSize);
    Kokkos::parallel_for(
        "InitTestTemperatureData", LocalActiveDomainSize, KOKKOS_LAMBDA(const int &D3D1ConvPosition) {
            int Rem = D3D1ConvPosition % (nx * MyYSlices);
            int RankX = Rem / MyYSlices;
            int RankY = Rem % MyYSlices;
            // Assign some of these a value of 0 (these will be solid cells), and others a positive value
            // (these will be tempsolid cells)
            if (RankX + RankY % 2 == 0)
                NumberOfSolidificationEvents(D3D1ConvPosition) = 0;
            else
                NumberOfSolidificationEvents(D3D1ConvPosition) = 1;
        });

    // Call constructor
    CellData<TEST_MEMSPACE> cellData(LocalDomainSize, LocalActiveDomainSize, nx, MyYSlices, ZBound_Low);
    cellData.init_substrate("", false, false, PowderFirstLayer, nx, ny, nz, LayerHeight, LocalActiveDomainSize,
                            ZMaxLayer, ZMin, deltax, MyYSlices, MyYOffset, ZBound_Low, id, RNGSeed,
                            SubstrateGrainSpacing, PowderActiveFraction, NumberOfSolidificationEvents);

    // Copy GrainID results back to host to check first layer's initialization
    ViewI_H GrainID_AllLayers_H = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cellData.GrainID_AllLayers);

    // Baseplate grains - cells should have GrainIDs between 1 and np (inclusive)
    for (int i = 0; i < BaseplateSize; i++) {
        EXPECT_GT(GrainID_AllLayers_H(i), 0);
        EXPECT_LT(GrainID_AllLayers_H(i), np + 1);
    }

    // In the case where there is powder on the first layer, the baseplate size is "LayerHeight" fewer cells in the Z
    // direction, and there should be LayerHeight * nx * MyYSlices * np powder grain IDs
    int ExpectedNumPowderGrainsPerLayer = LayerHeight * nx * MyYSlices * np;
    if (PowderFirstLayer) {
        EXPECT_EQ(BaseplateSize, nx * MyYSlices * (nzActive - LayerHeight));
        // Powder grains should have unique Grain ID values larger than 0 and smaller than
        // NextLayer_FirstEpitaxialGrainID
        EXPECT_EQ(cellData.NextLayer_FirstEpitaxialGrainID, np + 1 + ExpectedNumPowderGrainsPerLayer);
        int BottomPowderLayer = nx * MyYSlices * (ZBound_Low + nzActive - LayerHeight);
        int TopPowderLayer = nx * MyYSlices * (ZBound_Low + nzActive);
        for (int i = BottomPowderLayer; i < TopPowderLayer; i++) {
            EXPECT_GT(GrainID_AllLayers_H(i), 0);
            EXPECT_LT(GrainID_AllLayers_H(i), cellData.NextLayer_FirstEpitaxialGrainID);
        }
    }
    else {
        EXPECT_EQ(BaseplateSize, nx * MyYSlices * nzActive);
        // Next unused GrainID should be the number of grains present in the baseplate plus 1 (since GrainID = 0 is not
        // used for any baseplate grains)
        EXPECT_EQ(cellData.NextLayer_FirstEpitaxialGrainID, np + 1);
    }

    // Copy cell types back to host to check
    ViewI_H CellType_AllLayers_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cellData.CellType_AllLayers);
    for (int D3D1ConvPosition = 0; D3D1ConvPosition < LocalActiveDomainSize; D3D1ConvPosition++) {
        int Rem = D3D1ConvPosition % (nx * MyYSlices);
        int RankX = Rem / MyYSlices;
        int RankY = Rem % MyYSlices;
        // Cells with no associated solidification events should be solid, others TempSolid
        if (RankX + RankY % 2 == 0)
            EXPECT_EQ(CellType_AllLayers_Host(D3D1ConvPosition), Solid);
        else
            EXPECT_EQ(CellType_AllLayers_Host(D3D1ConvPosition), TempSolid);
    }
    int PreviousLayer_FirstEpitaxialGrainID = cellData.NextLayer_FirstEpitaxialGrainID;
    ZBound_Low = 1;
    // Initialize the next layer using the same time-temperature history - powder should span cells at Z = 3
    cellData.init_next_layer(1, id, nx, ny, MyYSlices, MyYOffset, ZBound_Low, LayerHeight, LocalActiveDomainSize, false,
                             false, ZMin, ZMaxLayer, deltax, RNGSeed, PowderActiveFraction,
                             NumberOfSolidificationEvents);

    // Copy all grain IDs for all layers back to the host to check that they match
    // and that the powder layer was initialized correctly
    GrainID_AllLayers_H = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cellData.GrainID_AllLayers);

    // Check that the right number of powder grains were initialized
    EXPECT_EQ(cellData.NextLayer_FirstEpitaxialGrainID,
              PreviousLayer_FirstEpitaxialGrainID + ExpectedNumPowderGrainsPerLayer);

    // Powder grains should have unique Grain ID values between PreviousLayer_FirstEpitaxialGrainID and
    // NextLayer_FirstEpitaxialGrainID - 1
    int BottomPowderLayer = nx * MyYSlices * (ZBound_Low + nzActive - LayerHeight);
    int TopPowderLayer = nx * MyYSlices * (ZBound_Low + nzActive);
    for (int GlobalD3D1ConvPosition = BottomPowderLayer; GlobalD3D1ConvPosition < TopPowderLayer;
         GlobalD3D1ConvPosition++) {
        EXPECT_GT(GrainID_AllLayers_H(GlobalD3D1ConvPosition), PreviousLayer_FirstEpitaxialGrainID - 1);
        EXPECT_LT(GrainID_AllLayers_H(GlobalD3D1ConvPosition), cellData.NextLayer_FirstEpitaxialGrainID);
    }

    // Subview grain IDs should match the grain IDs overall
    ViewI GrainID = cellData.getGrainIDSubview();
    ViewI_H GrainID_H = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);
    for (int D3D1ConvPosition = 0; D3D1ConvPosition < LocalActiveDomainSize; D3D1ConvPosition++) {
        int RankZ = D3D1ConvPosition / (nx * MyYSlices);
        int GlobalZ = RankZ + ZBound_Low;
        int Rem = D3D1ConvPosition % (nx * MyYSlices);
        int RankX = Rem / MyYSlices;
        int RankY = Rem % MyYSlices;
        int GlobalD3D1ConvPosition = GlobalZ * nx * MyYSlices + RankX * MyYSlices + RankY;
        EXPECT_EQ(GrainID_H(D3D1ConvPosition), GrainID_AllLayers_H(GlobalD3D1ConvPosition));
    }

    // Copy cell types back to host to check - should be the same as the previous layer as the same time-temperature
    // history was used
    CellType_AllLayers_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cellData.CellType_AllLayers);
    for (int D3D1ConvPosition = 0; D3D1ConvPosition < LocalActiveDomainSize; D3D1ConvPosition++) {
        int Rem = D3D1ConvPosition % (nx * MyYSlices);
        int RankX = Rem / MyYSlices;
        int RankY = Rem % MyYSlices;
        // Cells with no associated solidification events should be solid, others TempSolid
        if (RankX + RankY % 2 == 0)
            EXPECT_EQ(CellType_AllLayers_Host(D3D1ConvPosition), Solid);
        else
            EXPECT_EQ(CellType_AllLayers_Host(D3D1ConvPosition), TempSolid);
    }
}
//---------------------------------------------------------------------------//
// nuclei_init_tests
//---------------------------------------------------------------------------//
void testNucleiInit() {

    using memory_space = TEST_MEMSPACE;
    using view_int_host = Kokkos::View<int *, Kokkos::HostSpace>;
    using view_float3d_host = Kokkos::View<float ***, Kokkos::HostSpace>;

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Create test data
    // Only Z = 1 through 4 is part of this layer
    int nzActive = 4;
    int ZBound_Low = 1;
    int nx = 4;
    // Each rank is assigned a different portion of the domain in Y
    int ny = 2 * np;
    int MyYSlices = 2;
    int MyYOffset = 2 * id;
    double deltax = 1;
    int LocalActiveDomainSize = nx * MyYSlices * nzActive;
    // MPI rank locations relative to the global grid
    bool AtNorthBoundary, AtSouthBoundary;
    if (id == 0)
        AtSouthBoundary = true;
    else
        AtSouthBoundary = false;
    if (id == np - 1)
        AtNorthBoundary = true;
    else
        AtNorthBoundary = false;

    // There are 40 * np total cells in this domain (nx * ny * nz)
    // Each rank has 40 cells - the top 32 cells are part of the active layer and are candidates for nucleation
    // assignment
    double NMax = 0.125; // This nucleation density ensures there will be 4 potential nuclei per MPI rank present
                         // without remelting (each cell solidifies once)
    int MaxPotentialNuclei_PerPass = 4 * np;
    // A cell can solidify 1-3 times
    int MaxSolidificationEvents_Count = 3;
    double dTN = 1;
    double dTsigma = 0.0001;
    double RNGSeed = 0.0;

    // Initialize MaxSolidificationEvents to 3 for each layer. LayerTimeTempHistory and NumberOfSolidificationEvents are
    // initialized for each cell on the host and copied to the device
    view_int_host MaxSolidificationEvents_Host(Kokkos::ViewAllocateWithoutInitializing("MaxSolidificationEvents_Host"),
                                               2);
    view_int_host NumberOfSolidificationEvents_Host(
        Kokkos::ViewAllocateWithoutInitializing("NumberOfSolidificationEvents_Host"), LocalActiveDomainSize);
    view_float3d_host LayerTimeTempHistory_Host(Kokkos::ViewAllocateWithoutInitializing("LayerTimeTempHistory_Host"),
                                                LocalActiveDomainSize, MaxSolidificationEvents_Count, 3);
    MaxSolidificationEvents_Host(0) = MaxSolidificationEvents_Count;
    MaxSolidificationEvents_Host(1) = MaxSolidificationEvents_Count;
    // Cells solidify 1, 2, or 3 times, depending on their X coordinate
    for (int k = 0; k < nzActive; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                int D3D1ConvPosition = k * nx * MyYSlices + i * MyYSlices + j;
                if (i < nx / 2 - 1)
                    NumberOfSolidificationEvents_Host(D3D1ConvPosition) = 3;
                else if (i < nx / 2)
                    NumberOfSolidificationEvents_Host(D3D1ConvPosition) = 2;
                else
                    NumberOfSolidificationEvents_Host(D3D1ConvPosition) = 1;
            }
        }
    }
    for (int n = 0; n < MaxSolidificationEvents_Count; n++) {
        for (int RankZ = 0; RankZ < nzActive; RankZ++) {
            for (int RankX = 0; RankX < nx; RankX++) {
                for (int RankY = 0; RankY < MyYSlices; RankY++) {
                    int D3D1ConvPosition = RankZ * nx * MyYSlices + RankX * MyYSlices + RankY;
                    int GlobalZ = RankZ + ZBound_Low;
                    if (n < NumberOfSolidificationEvents_Host(D3D1ConvPosition)) {
                        LayerTimeTempHistory_Host(D3D1ConvPosition, n, 0) =
                            GlobalZ + RankY + MyYOffset +
                            (LocalActiveDomainSize * n); // melting time step depends on solidification event number
                        LayerTimeTempHistory_Host(D3D1ConvPosition, n, 1) =
                            GlobalZ + RankY + MyYOffset + 1 +
                            (LocalActiveDomainSize * n); // liquidus time stemp depends on solidification event number
                        LayerTimeTempHistory_Host(D3D1ConvPosition, n, 2) =
                            1.2; // ensures that a cell's nucleation time will be 1 time step after its CritTimeStep
                                 // value
                    }
                }
            }
        }
    }

    auto MaxSolidificationEvents = Kokkos::create_mirror_view_and_copy(memory_space(), MaxSolidificationEvents_Host);
    auto NumberOfSolidificationEvents =
        Kokkos::create_mirror_view_and_copy(memory_space(), NumberOfSolidificationEvents_Host);
    auto LayerTimeTempHistory = Kokkos::create_mirror_view_and_copy(memory_space(), LayerTimeTempHistory_Host);

    // Nucleation data structure, containing views of nuclei locations, time steps, and ids, and nucleation event
    // counters - initialized with an estimate on the number of nuclei in the layer Without knowing
    // PossibleNuclei_ThisRankThisLayer yet, initialize nucleation data structures to estimated sizes, resize inside of
    // NucleiInit when the number of nuclei per rank is known
    int EstimatedNuclei_ThisRankThisLayer = NMax * pow(deltax, 3) * LocalActiveDomainSize;
    Nucleation<memory_space> nucleation(
        EstimatedNuclei_ThisRankThisLayer, NMax, deltax,
        100); // NucleiGrainID should start at -101 - supply optional input arg to constructor

    // Fill in nucleation data structures, and assign nucleation undercooling values to potential nucleation events
    // Potential nucleation grains are only associated with liquid cells in layer 1 - they will be initialized for each
    // successive layer when layer 1 in complete
    nucleation.placeNuclei(MaxSolidificationEvents, NumberOfSolidificationEvents, LayerTimeTempHistory, RNGSeed, 1, nx,
                           ny, nzActive, dTN, dTsigma, MyYSlices, MyYOffset, ZBound_Low, id, AtNorthBoundary,
                           AtSouthBoundary);

    // Copy results back to host to check
    auto NucleiLocation_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), nucleation.NucleiLocations);
    auto NucleiGrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), nucleation.NucleiGrainID);

    // Was the nucleation counter initialized to zero?
    EXPECT_EQ(nucleation.NucleationCounter, 0);

    // Is the total number of nuclei in the system correct, based on the number of remelting events? Equal probability
    // of creating a nucleus each time a cell resolidifies
    int ExpectedNucleiPerRank = 100 + MaxSolidificationEvents_Count * MaxPotentialNuclei_PerPass;
    EXPECT_EQ(nucleation.Nuclei_WholeDomain, ExpectedNucleiPerRank);
    for (int n = 0; n < nucleation.PossibleNuclei; n++) {
        // Are the nuclei grain IDs negative numbers in the expected range based on the inputs?
        EXPECT_GT(NucleiGrainID_Host(n), -(100 + ExpectedNucleiPerRank * np + 1));
        EXPECT_LT(NucleiGrainID_Host(n), -100);
        // Are the correct undercooling values associated with the correct cell locations?
        // Cell location is a local position (relative to the bottom of the layer)
        int CellLocation = NucleiLocation_Host(n);
        int RankZ = CellLocation / (nx * MyYSlices);
        int Rem = CellLocation % (nx * MyYSlices);
        int RankY = Rem % MyYSlices;
        int GlobalZ = RankZ + ZBound_Low;
        // Expected nucleation time with remelting can be one of 3 possibilities, depending on the associated
        // solidification event
        int Expected_NucleationTimeNoRM = GlobalZ + RankY + MyYOffset + 2;
        int AssociatedSEvent = nucleation.NucleationTimes_Host(n) / LocalActiveDomainSize;
        int Expected_NucleationTimeRM = Expected_NucleationTimeNoRM + AssociatedSEvent * LocalActiveDomainSize;
        EXPECT_EQ(nucleation.NucleationTimes_Host(n), Expected_NucleationTimeRM);

        // Are the nucleation events in order of the time steps at which they may occur?
        if (n < nucleation.PossibleNuclei - 2) {
            EXPECT_LE(nucleation.NucleationTimes_Host(n), nucleation.NucleationTimes_Host(n + 1));
        }
    }
}

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, fileread_test) {
    // input argument: test for two different print versions for last input file
    testInputReadFromFile(0);
    testInputReadFromFile(1);
}
TEST(TEST_CATEGORY, cell_init_tests) {
    testCellDataInit_ConstrainedGrowth();
    // For non-constrained solidification problems, test w/ and w/o space left for powder layer
    testCellDataInit(true);
    testCellDataInit(false);
}
TEST(TEST_CATEGORY, nuclei_init_test) { testNucleiInit(); }
} // end namespace Test
