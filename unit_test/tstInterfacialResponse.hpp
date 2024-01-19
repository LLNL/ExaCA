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
        double a_test, b_test, c_test, d_test, freezing_range_test, expected_v;
        double loc_u = 11.0;
        if (file_name == "Inconel625.json") {
            a_test = -0.00000010302;
            b_test = 0.00010533;
            c_test = 0.0022196;
            d_test = 0;
            freezing_range_test = 210;
            expected_v =
                (deltat / deltax) * (a_test * pow(loc_u, 3.0) + b_test * pow(loc_u, 2.0) + c_test * loc_u + d_test);
        }
        else if (file_name == "Inconel625_Quadratic.json") {
            a_test = 0.000072879;
            b_test = 0.004939;
            c_test = -0.047024;
            freezing_range_test = 210;
            expected_v = (deltat / deltax) * (a_test * pow(loc_u, 2.0) + b_test * loc_u + c_test);
        }
        else if (file_name == "SS316.json") {
            a_test = 0.000007325;
            b_test = 3.12;
            c_test = 0;
            freezing_range_test = 26.5;
            expected_v = (deltat / deltax) * (a_test * pow(loc_u, b_test) + c_test);
        }
        else {
            throw std::runtime_error("File not set up for testing.");
        }
        // For all IRFs, A and C should be normalized by deltat/deltax (i.e., 2)
        // For the power law IRF (SS316), B is dimensionless and should not be normalized unlike the other IRFs where
        // all coefficients are normalized
        EXPECT_DOUBLE_EQ(irf.A, a_test * 2);
        if (file_name == "SS316.json")
            EXPECT_DOUBLE_EQ(irf.B, b_test);
        else
            EXPECT_DOUBLE_EQ(irf.B, b_test * 2);
        EXPECT_DOUBLE_EQ(irf.C, c_test * 2);
        if (file_name == "Inconel625.json") {
            EXPECT_DOUBLE_EQ(irf.D, d_test * 2);
        }
        EXPECT_DOUBLE_EQ(irf.freezing_range, freezing_range_test);
        double computer_v = irf.compute(loc_u);
        EXPECT_DOUBLE_EQ(computer_v, expected_v);
    }
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, interfacial_response) { testInterfacialResponse(); }
} // end namespace Test
