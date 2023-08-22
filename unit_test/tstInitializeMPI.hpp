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
        int TempFilesInSeries, NumberOfLayers, LayerHeight, nx, ny, nz, NSpotsX, NSpotsY, SpotOffset, SpotRadius,
            singleGrainOrientation;
        float SubstrateGrainSpacing;
        double deltax, NMax, dTN, dTsigma, HT_deltax, deltat, G, R, FractSurfaceSitesActive, RNGSeed,
            PowderActiveFraction, initUndercooling;
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
                          BaseplateThroughPowder, PowderActiveFraction, LayerwiseTempRead, PowderFirstLayer, print,
                          initUndercooling, singleGrainOrientation);
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
            EXPECT_DOUBLE_EQ(initUndercooling, 0.0);
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
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, fileread_test) {
    // input argument: test for two different print versions for last input file
    testInputReadFromFile(0);
    testInputReadFromFile(1);
}
} // end namespace Test
