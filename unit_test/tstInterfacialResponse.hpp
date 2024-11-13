// Copyright Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>
#include <nlohmann/json.hpp>

#include "CAinputs.hpp"
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
    std::vector<std::string> material_file_names = {"Inconel625.json", "Inconel625_Quadratic.json", "SS316.json",
                                                    "SS316L_AusteniteFerrite.json"};
    double deltax = 0.5;
    double deltat = 1.0;
    for (auto file_name : material_file_names) {
        std::cout << "Reading " << file_name << std::endl;

        Inputs inputs;
        inputs.material_filename = file_name;
        inputs.parseIRF(0);
        InterfacialResponseFunction irf(deltat, deltax, inputs.irf);

        // Check that fitting parameters were correctly initialized and normalized for the single phase present or for
        // both phases (SS316L_FA material input file) Fitting parameters should've been normalized by deltat / deltax,
        // i.e. twice as large as the numbers in the file
        float a_test[2], b_test[2], c_test[2], d_test[2], freezing_range_test[2], expected_v[2];
        float loc_u = 11.0;
        if (file_name == material_file_names[0]) {
            a_test[0] = -0.00000010302;
            b_test[0] = 0.00010533;
            c_test[0] = 0.0022196;
            d_test[0] = 0;
            freezing_range_test[0] = 210;
            expected_v[0] = (deltat / deltax) *
                            (a_test[0] * pow(loc_u, 3.0) + b_test[0] * pow(loc_u, 2.0) + c_test[0] * loc_u + d_test[0]);
        }
        else if (file_name == material_file_names[1]) {
            a_test[0] = 0.000072879;
            b_test[0] = 0.004939;
            c_test[0] = -0.047024;
            freezing_range_test[0] = 210;
            expected_v[0] = (deltat / deltax) * (a_test[0] * pow(loc_u, 2.0) + b_test[0] * loc_u + c_test[0]);
        }
        else if (file_name == material_file_names[2]) {
            a_test[0] = 0.000007325;
            b_test[0] = 3.12;
            c_test[0] = 0;
            freezing_range_test[0] = 26.5;
            expected_v[0] = (deltat / deltax) * (a_test[0] * pow(loc_u, b_test[0]) + c_test[0]);
        }
        else if (file_name == material_file_names[3]) {
            // Austenite inputs
            a_test[0] = 3.1718E-18;
            b_test[0] = 12.024;
            c_test[0] = 0;
            freezing_range_test[0] = 24.9;
            expected_v[0] = (deltat / deltax) * (a_test[0] * pow(loc_u, b_test[0]) + c_test[0]);
            // Ferrite inputs
            a_test[1] = 5.6442E-5;
            b_test[1] = 0.37259;
            c_test[1] = 0;
            freezing_range_test[1] = 23.9;
            expected_v[1] = (deltat / deltax) * (a_test[1] * exp(loc_u * b_test[1]) + c_test[1]);
        }
        else {
            throw std::runtime_error("File not set up for testing.");
        }
        // For all IRFs, A and C should be normalized by deltat/deltax (i.e., 2)
        // For the power law IRF (SS316 and both SS316L phases), B is dimensionless and should not be normalized unlike
        // the other IRFs where all coefficients are normalized Check that the number of phases are correct
        if (file_name == material_file_names[3]) {
            EXPECT_EQ(irf.numPhases(), 2);
            // Check for correct phase names
            EXPECT_TRUE(irf._inputs.phase_names[0] == "austenite");
            EXPECT_TRUE(irf._inputs.phase_names[1] == "ferrite");
            // Check for correct transformation
            EXPECT_EQ(irf._inputs.transformation, 1);
            // Check for velocity cap for both phases
            EXPECT_TRUE(irf._inputs.velocity_cap[0]);
            EXPECT_TRUE(irf._inputs.velocity_cap[1]);
        }
        else
            EXPECT_EQ(irf.numPhases(), 1);
        for (int phase = 0; phase < irf.numPhases(); phase++) {
            EXPECT_FLOAT_EQ(irf.A(phase), a_test[phase] * 2);
            if ((irf._inputs.function[phase] == irf._inputs.power) ||
                (irf._inputs.function[phase] == irf._inputs.exponential))
                EXPECT_FLOAT_EQ(irf.B(phase), b_test[phase]);
            else
                EXPECT_FLOAT_EQ(irf.B(phase), b_test[phase] * 2);
            EXPECT_FLOAT_EQ(irf.C(phase), c_test[phase] * 2);
            if (file_name == material_file_names[0]) {
                EXPECT_FLOAT_EQ(irf.D(phase), d_test[phase] * 2);
            }
            EXPECT_FLOAT_EQ(irf.freezingRange(phase), freezing_range_test[phase]);
            float computed_v = irf.compute(loc_u, phase);
            EXPECT_FLOAT_EQ(computed_v, expected_v[phase]);
        }
    }
}

void testGetPreferredPhase() {

    // Single phase IRF
    Inputs inputs;
    inputs.irf.function[0] = inputs.irf.power;
    // IRF for phase 0: V = 1.0 * (loc_u)**0.5
    inputs.irf.A[0] = 1.0;
    inputs.irf.B[0] = 0.5;
    inputs.irf.C[0] = 0.0;
    inputs.domain.deltax = 1.0 * pow(10, -6);
    inputs.domain.deltat = 1.0 * pow(10, -6);
    InterfacialResponseFunction irf(inputs.domain.deltat, inputs.domain.deltax, inputs.irf);
    // Regardless of undercooling, the preferred phase should be 0 for single phase IRFs
    EXPECT_EQ(irf.getPreferredPhaseNucleation(1.0), 0);
    EXPECT_EQ(irf.getPreferredPhaseActivation(0), 0);

    // Two phase IRF: phase 1 IRF: V = 0.5 * (loc_u)**2.0
    irf._inputs.num_phases = 2;
    irf._inputs.transformation = 1;
    irf._inputs.function[1] = inputs.irf.power;
    irf._inputs.A[1] = 0.5;
    irf._inputs.B[1] = 2.0;
    irf._inputs.C[1] = 0.0;
    // Preferred phase is phase with larger V: should be phase 0 for loc_u of 1, phase 1 for loc_u of 2
    EXPECT_EQ(irf.getPreferredPhaseNucleation(1.0), 0);
    EXPECT_EQ(irf.getPreferredPhaseNucleation(2.0), 1);
    // With a solidification transformation, preferred phase for activation is always 0
    EXPECT_EQ(irf.getPreferredPhaseActivation(1), 0);
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, interfacial_response) {
    testInterfacialResponse();
    testGetPreferredPhase();
}
} // end namespace Test
