// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>
#include <nlohmann/json.hpp>

#include "CAinterfacialresponse.hpp"

#include <gtest/gtest.h>

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
// Tests the interfacial response function construction
//---------------------------------------------------------------------------//
void testInterfacialResponse() {

    // Test that the interfacial response can be read for the new file format
    std::vector<std::string> material_file_names = {"Inconel625.json", "Inconel625_Quadratic.json", "SS316.json"};
    double deltax = 0.5;
    double deltat = 1.0;
    for (auto file_name : material_file_names) {
        std::cout << "Reading " << file_name << std::endl;
        InterfacialResponseFunction irf(0, file_name, deltat, deltax);

        // Check that fitting parameters were correctly initialized and normalized
        // Fitting parameters should've been normalized by deltat / deltax, i.e. twice as large as the numbers in the
        // file
        double ATest, BTest, CTest, DTest, FreezingRangeTest, ExpectedV;
        double LocU = 11.0;
        if (file_name == "Inconel625.json") {
            ATest = -0.00000010302;
            BTest = 0.00010533;
            CTest = 0.0022196;
            DTest = 0;
            FreezingRangeTest = 210;
            ExpectedV = (deltat / deltax) * (ATest * pow(LocU, 3.0) + BTest * pow(LocU, 2.0) + CTest * LocU + DTest);
        }
        else if (file_name == "Inconel625_Quadratic.json") {
            ATest = 0.000072879;
            BTest = 0.004939;
            CTest = -0.047024;
            FreezingRangeTest = 210;
            ExpectedV = (deltat / deltax) * (ATest * pow(LocU, 2.0) + BTest * LocU + CTest);
        }
        else if (file_name == "SS316.json") {
            ATest = 0.000007325;
            BTest = 3.12;
            CTest = 0;
            FreezingRangeTest = 26.5;
            ExpectedV = (deltat / deltax) * (ATest * pow(LocU, BTest) + CTest);
        }
        else {
            throw std::runtime_error("File not set up for testing.");
        }
        // For all IRFs, A and C should be normalized by deltat/deltax (i.e., 2)
        // For the power law IRF (SS316), B is dimensionless and should not be normalized unlike the other IRFs where
        // all coefficients are normalized
        EXPECT_DOUBLE_EQ(irf.A, ATest * 2);
        if (file_name == "SS316.json")
            EXPECT_DOUBLE_EQ(irf.B, BTest);
        else
            EXPECT_DOUBLE_EQ(irf.B, BTest * 2);
        EXPECT_DOUBLE_EQ(irf.C, CTest * 2);
        if (file_name == "Inconel625.json") {
            EXPECT_DOUBLE_EQ(irf.D, DTest * 2);
        }
        EXPECT_DOUBLE_EQ(irf.FreezingRange, FreezingRangeTest);
        double ComputedV = irf.compute(LocU);
        EXPECT_DOUBLE_EQ(ComputedV, ExpectedV);
    }
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, interfacial_response) { testInterfacialResponse(); }
} // end namespace Test
