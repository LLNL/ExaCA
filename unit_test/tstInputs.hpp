// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAinputs.hpp"
#include "CAtypes.hpp"

#include <gtest/gtest.h>

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
    TestDataFile << "      \"BaseplateTopZ\": -0.00625" << std::endl;
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

void testInputs(int PrintVersion) {

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
    std::vector<std::string> InputFilenames = {"Inp_DirSolidification", "Inp_SpotMelt", "Inp_TemperatureTest",
                                               "Inp_TwoGrainDirSolidification"};
    for (int n = 0; n < 4; n++) {
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
        std::cout << "Reading " << FileName << std::endl;
        // Data printing structure - contains print options (false by default) and functions
        Inputs inputs(id, FileName);
        InterfacialResponseFunction irf(0, inputs.MaterialFileName, inputs.domain.deltat, inputs.domain.deltax);
        MPI_Barrier(MPI_COMM_WORLD);

        // Check the results
        // The existence of the specified orientation, substrate, and temperature filenames was already checked within
        // InputReadFromFile
        // These should be the same for all test problems (except the 4th one, which has 0 nucleation density)
        EXPECT_DOUBLE_EQ(inputs.domain.deltax, 1.0 * pow(10, -6));
        if (FileName == InputFilenames[3])
            EXPECT_DOUBLE_EQ(inputs.nucleation.NMax, 0.0);
        else
            EXPECT_DOUBLE_EQ(inputs.nucleation.NMax, 1.0 * pow(10, 13));
        EXPECT_DOUBLE_EQ(inputs.nucleation.dTN, 5.0);
        EXPECT_DOUBLE_EQ(inputs.nucleation.dTsigma, 0.5);
        EXPECT_DOUBLE_EQ(irf.A, -0.00000010302 * inputs.domain.deltat / inputs.domain.deltax);
        EXPECT_DOUBLE_EQ(irf.B, 0.00010533 * inputs.domain.deltat / inputs.domain.deltax);
        EXPECT_DOUBLE_EQ(irf.C, 0.0022196 * inputs.domain.deltat / inputs.domain.deltax);
        EXPECT_DOUBLE_EQ(irf.D, 0);
        EXPECT_DOUBLE_EQ(irf.FreezingRange, 210);

        // These are different for all 3 test problems
        if ((FileName == InputFilenames[0]) || (FileName == InputFilenames[3])) {
            EXPECT_TRUE(inputs.print.PrintTimeSeries);
            EXPECT_EQ(inputs.print.TimeSeriesInc, 5250);
            EXPECT_FALSE(inputs.print.PrintIdleTimeSeriesFrames);
            EXPECT_DOUBLE_EQ(inputs.temperature.G, 500000.0);
            EXPECT_DOUBLE_EQ(inputs.temperature.R, 300000.0);
            EXPECT_DOUBLE_EQ(inputs.temperature.initUndercooling, 0.0);
            // compare with float to avoid floating point error with irrational number
            float deltat_comp = static_cast<float>(inputs.domain.deltat);
            EXPECT_FLOAT_EQ(deltat_comp, pow(10, -6) / 15.0);
            EXPECT_EQ(inputs.domain.nx, 200);
            EXPECT_EQ(inputs.domain.ny, 200);
            EXPECT_EQ(inputs.domain.nz, 200);
            if (FileName == InputFilenames[0]) {
                EXPECT_DOUBLE_EQ(inputs.substrate.FractSurfaceSitesActive, 0.08);
                EXPECT_TRUE(inputs.print.BaseFileName == "TestProblemDirS");
                EXPECT_FALSE(inputs.substrate.CustomGrainLocationsIDs);
            }
            else {
                EXPECT_EQ(inputs.substrate.GrainLocationsX[0], 100);
                EXPECT_EQ(inputs.substrate.GrainLocationsY[0], 50);
                EXPECT_EQ(inputs.substrate.GrainIDs[0], 25);
                EXPECT_EQ(inputs.substrate.GrainLocationsX[1], 100);
                EXPECT_EQ(inputs.substrate.GrainLocationsY[1], 150);
                EXPECT_EQ(inputs.substrate.GrainIDs[1], 9936);
                EXPECT_TRUE(inputs.print.BaseFileName == "TestProblemTwoGrainDirS");
                EXPECT_TRUE(inputs.substrate.CustomGrainLocationsIDs);
            }
            EXPECT_TRUE(inputs.print.PrintInitCritTimeStep);
            EXPECT_FALSE(inputs.print.PrintInitGrainID);
            EXPECT_FALSE(inputs.print.PrintInitLayerID);
            EXPECT_FALSE(inputs.print.PrintInitMeltTimeStep);
            EXPECT_FALSE(inputs.print.PrintInitUndercoolingChange);
            EXPECT_TRUE(inputs.print.PrintFinalGrainID);
            EXPECT_TRUE(inputs.print.PrintFinalLayerID);
            EXPECT_TRUE(inputs.print.PrintFinalMisorientation);
            EXPECT_FALSE(inputs.print.PrintFinalUndercoolingCurrent);
            EXPECT_FALSE(inputs.print.PrintFinalMeltTimeStep);
            EXPECT_FALSE(inputs.print.PrintFinalCritTimeStep);
            EXPECT_FALSE(inputs.print.PrintFinalUndercoolingChange);
            EXPECT_FALSE(inputs.print.PrintDefaultRVE);
            EXPECT_DOUBLE_EQ(inputs.RNGSeed, 0.0);
            EXPECT_FALSE(inputs.print.PrintBinary);
        }
        else if (FileName == InputFilenames[1]) {
            EXPECT_TRUE(inputs.print.PrintTimeSeries);
            EXPECT_EQ(inputs.print.TimeSeriesInc, 37500);
            EXPECT_TRUE(inputs.print.PrintIdleTimeSeriesFrames);
            EXPECT_DOUBLE_EQ(inputs.temperature.G, 500000.0);
            EXPECT_DOUBLE_EQ(inputs.temperature.R, 300000.0);
            // compare with float to avoid floating point error with irrational number
            float deltat_comp = static_cast<float>(inputs.domain.deltat);
            EXPECT_FLOAT_EQ(deltat_comp, pow(10, -6) / 15.0);
            EXPECT_EQ(inputs.domain.NSpotsX, 3);
            EXPECT_EQ(inputs.domain.NSpotsY, 2);
            EXPECT_EQ(inputs.domain.SpotOffset, 100);
            EXPECT_EQ(inputs.domain.SpotRadius, 75);
            EXPECT_EQ(inputs.domain.NumberOfLayers, 2);
            EXPECT_EQ(inputs.domain.LayerHeight, 20);
            EXPECT_FALSE(inputs.substrate.UseSubstrateFile);
            EXPECT_FALSE(inputs.substrate.BaseplateThroughPowder);
            // Option defaults to 0.0
            EXPECT_DOUBLE_EQ(inputs.substrate.BaseplateTopZ, 0.0);
            EXPECT_FLOAT_EQ(inputs.substrate.SubstrateGrainSpacing, 25.0);
            EXPECT_TRUE(inputs.print.BaseFileName == "TestProblemSpot");
            EXPECT_TRUE(inputs.print.PrintInitCritTimeStep);
            EXPECT_TRUE(inputs.print.PrintInitMeltTimeStep);
            EXPECT_FALSE(inputs.print.PrintInitGrainID);
            EXPECT_FALSE(inputs.print.PrintInitLayerID);
            EXPECT_FALSE(inputs.print.PrintInitUndercoolingChange);
            EXPECT_TRUE(inputs.print.PrintFinalGrainID);
            EXPECT_TRUE(inputs.print.PrintFinalLayerID);
            EXPECT_TRUE(inputs.print.PrintFinalMisorientation);
            EXPECT_FALSE(inputs.print.PrintFinalUndercoolingCurrent);
            EXPECT_FALSE(inputs.print.PrintFinalMeltTimeStep);
            EXPECT_FALSE(inputs.print.PrintFinalCritTimeStep);
            EXPECT_FALSE(inputs.print.PrintFinalUndercoolingChange);
            EXPECT_FALSE(inputs.print.PrintDefaultRVE);
            EXPECT_DOUBLE_EQ(inputs.RNGSeed, 0.0);
        }
        else if (FileName == InputFilenames[2]) {
            EXPECT_DOUBLE_EQ(inputs.domain.deltat, 1.5 * pow(10, -6));
            EXPECT_EQ(inputs.temperature.TempFilesInSeries, 2);
            EXPECT_EQ(inputs.domain.NumberOfLayers, 2);
            EXPECT_EQ(inputs.domain.LayerHeight, 1);
            EXPECT_TRUE(inputs.substrate.UseSubstrateFile);
            EXPECT_FALSE(inputs.temperature.LayerwiseTempRead);
            EXPECT_DOUBLE_EQ(inputs.substrate.PowderActiveFraction, 0.001);
            // -0.00625 was input
            EXPECT_DOUBLE_EQ(inputs.substrate.BaseplateTopZ, -0.00625);
            EXPECT_DOUBLE_EQ(inputs.temperature.HT_deltax, inputs.domain.deltax);
            EXPECT_TRUE(inputs.print.BaseFileName == "Test");
            EXPECT_TRUE(inputs.temperature.temp_paths[0] == ".//1DummyTemperature.txt");
            EXPECT_TRUE(inputs.temperature.temp_paths[1] == ".//2DummyTemperature.txt");
            if (PrintVersion == 0) {
                EXPECT_TRUE(inputs.print.PrintInitCritTimeStep);
                EXPECT_TRUE(inputs.print.PrintInitMeltTimeStep);
                EXPECT_TRUE(inputs.print.PrintInitGrainID);
                EXPECT_TRUE(inputs.print.PrintInitLayerID);
                EXPECT_TRUE(inputs.print.PrintInitUndercoolingChange);
                EXPECT_TRUE(inputs.print.PrintFinalMisorientation);
                EXPECT_TRUE(inputs.print.PrintFinalUndercoolingCurrent);
                EXPECT_TRUE(inputs.print.PrintFinalLayerID);
                EXPECT_TRUE(inputs.print.PrintFinalGrainID);
                EXPECT_TRUE(inputs.print.PrintTimeSeries);
                EXPECT_EQ(inputs.print.TimeSeriesInc, 200); // value from file divided by deltat
                EXPECT_TRUE(inputs.print.PrintIdleTimeSeriesFrames);
                EXPECT_FALSE(inputs.print.PrintDefaultRVE);
            }
            else if (PrintVersion == 1) {
                EXPECT_FALSE(inputs.print.PrintInitCritTimeStep);
                EXPECT_FALSE(inputs.print.PrintInitMeltTimeStep);
                EXPECT_FALSE(inputs.print.PrintInitGrainID);
                EXPECT_FALSE(inputs.print.PrintInitLayerID);
                EXPECT_FALSE(inputs.print.PrintInitUndercoolingChange);
                EXPECT_TRUE(inputs.print.PrintFinalMisorientation);
                EXPECT_FALSE(inputs.print.PrintFinalUndercoolingCurrent);
                EXPECT_TRUE(inputs.print.PrintFinalLayerID);
                EXPECT_TRUE(inputs.print.PrintFinalGrainID);
                EXPECT_TRUE(inputs.print.PrintDefaultRVE);
                EXPECT_FALSE(inputs.print.PrintTimeSeries);
            }
            EXPECT_FALSE(inputs.temperature.LayerwiseTempRead);
            EXPECT_DOUBLE_EQ(inputs.RNGSeed, 2.0);
            EXPECT_TRUE(inputs.print.PrintBinary);
        }
    }
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, inputs) {
    // input argument: test for two different print versions for last input file
    testInputs(0);
    testInputs(1);
}
} // end namespace Test
