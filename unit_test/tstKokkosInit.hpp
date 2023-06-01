// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAfunctions.hpp"
#include "CAinitialize.hpp"
#include "CAparsefiles.hpp"
#include "CAtypes.hpp"

#include <gtest/gtest.h>

#include "mpi.h"

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
void testOrientationInit_Vectors() {

    int id;
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    int ValsPerLine = 9;
    int NGrainOrientations = 0;
    std::string GrainOrientationFile = checkFileInstalled("GrainOrientationVectors.csv", id);

    // View for storing orientation data
    ViewF GrainOrientationData(Kokkos::ViewAllocateWithoutInitializing("GrainOrientationData"), 0);

    // Call OrientationInit - without optional final argument
    OrientationInit(id, NGrainOrientations, GrainOrientationData, GrainOrientationFile);

    // Copy orientation data back to the host
    ViewF_H GrainOrientationData_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainOrientationData);

    // Check results
    EXPECT_EQ(NGrainOrientations, 10000);

    std::vector<float> ExpectedGrainOrientations = {0.848294,  0.493303,  0.19248,  -0.522525, 0.720911,  0.455253,
                                                    0.0858167, -0.486765, 0.869308, 0.685431,  0.188182,  0.7034,
                                                    -0.468504, 0.85348,   0.228203, -0.557394, -0.485963, 0.673166};
    for (int n = 0; n < 2 * ValsPerLine; n++) {
        EXPECT_FLOAT_EQ(GrainOrientationData_Host(n), ExpectedGrainOrientations[n]);
    }
}

void testOrientationInit_Angles() {

    int id;
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    int ValsPerLine = 3;
    int NGrainOrientations = 0;
    std::string GrainOrientationFile = checkFileInstalled("GrainOrientationEulerAnglesBungeZXZ.csv", id);

    // View for storing orientation data
    ViewF GrainOrientationData(Kokkos::ViewAllocateWithoutInitializing("GrainOrientationData"), 0);

    // Call OrientationInit - with optional final argument
    OrientationInit(id, NGrainOrientations, GrainOrientationData, GrainOrientationFile, ValsPerLine);

    // Copy orientation data back to the host
    ViewF_H GrainOrientationData_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainOrientationData);

    // Check results
    EXPECT_EQ(NGrainOrientations, 10000);

    // Check first two orientations
    std::vector<float> ExpectedGrainOrientations = {9.99854, 29.62172, 22.91854, 311.08350, 47.68814, 72.02547};
    for (int n = 0; n < 2 * ValsPerLine; n++) {
        EXPECT_FLOAT_EQ(GrainOrientationData_Host(n), ExpectedGrainOrientations[n]);
    }
}

// Test unidirectional solidification problem for either directional solidification or growth of a single grain seed,
// with thermal gradient G in the domain
void testTempInit_UnidirectionalGradient(std::string SimulationType, double G) {

    int nx = 2;
    int ny = 1;
    int nz = 5; // (Front is at Z = 0 for directional growth, single grain seed at Z = 2 for singlegrain problem)
    int kCenter = floorf(static_cast<float>(nz) / 2.0);
    int LocalDomainSize = nx * ny * nz;
    double InitUndercooling = 10; // Used only for single grain problem

    // For problems with non-zero thermal gradient, 1 K difference between each cell and its neighbor in Z
    double deltax;
    if (G == 0)
        deltax = 1 * pow(10, -6);
    else
        deltax = 1.0 / G;
    double GNorm = G * deltax;
    // Cells cool at rate of 1 K per time step
    double R = 1000000;
    double deltat = 1 * pow(10, -6);
    double RNorm = R * deltat;

    // Temperature views
    ViewI CritTimeStep(Kokkos::ViewAllocateWithoutInitializing("CritTimeStep"), LocalDomainSize);
    ViewF UndercoolingChange(Kokkos::ViewAllocateWithoutInitializing("UndercoolingChange"), LocalDomainSize);
    ViewI LayerID(Kokkos::ViewAllocateWithoutInitializing("LayerID"), LocalDomainSize);
    ViewI NumberOfSolidificationEvents(Kokkos::ViewAllocateWithoutInitializing("NumberOfSolidificationEvents"),
                                       LocalDomainSize);
    ViewI SolidificationEventCounter(Kokkos::ViewAllocateWithoutInitializing("SolidificationEventCounter"),
                                     LocalDomainSize);
    ViewI MeltTimeStep(Kokkos::ViewAllocateWithoutInitializing("MeltTimeStep"), LocalDomainSize);
    ViewI MaxSolidificationEvents(Kokkos::ViewAllocateWithoutInitializing("MaxSolidificationEvents"), 1);
    ViewF UndercoolingCurrent(Kokkos::ViewAllocateWithoutInitializing("UndercoolingCurrent"), LocalDomainSize);
    ViewF3D LayerTimeTempHistory(Kokkos::ViewAllocateWithoutInitializing("LayerTimeTempHistory"), LocalDomainSize, 1,
                                 3);

    TempInit_UnidirectionalGradient(G, R, 0, nx, ny, deltax, deltat, nz, LocalDomainSize, CritTimeStep,
                                    UndercoolingChange, LayerID, NumberOfSolidificationEvents,
                                    SolidificationEventCounter, MeltTimeStep, MaxSolidificationEvents,
                                    LayerTimeTempHistory, InitUndercooling, UndercoolingCurrent, SimulationType);

    // Copy temperature views back to host
    // Copy orientation data back to the host
    ViewI_H CritTimeStep_Host = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), CritTimeStep); // Copy orientation data back to the host
    ViewF_H UndercoolingChange_Host = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), UndercoolingChange); // Copy orientation data back to the host
    ViewI_H LayerID_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), LayerID); // Copy orientation data back to the host
    ViewI_H NumberOfSolidificationEvents_Host = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), NumberOfSolidificationEvents); // Copy orientation data back to the host
    ViewI_H SolidificationEventCounter_Host = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), SolidificationEventCounter); // Copy orientation data back to the host
    ViewI_H MeltTimeStep_Host = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), MeltTimeStep); // Copy orientation data back to the host
    ViewI_H MaxSolidificationEvents_Host = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), MaxSolidificationEvents); // Copy orientation data back to the host
    ViewF_H UndercoolingCurrent_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), UndercoolingCurrent);
    ViewF3D_H LayerTimeTempHistory_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), LayerTimeTempHistory);

    // Check results
    EXPECT_EQ(MaxSolidificationEvents_Host(0), 1);
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                int D3D1ConvPosition = k * nx * ny + i * ny + j;
                // Each cell solidifies once, and counter should start at 0, associated with the zeroth layer
                // MeltTimeStep should be -1 for all cells
                // Cells cool at 1 K per time step
                EXPECT_FLOAT_EQ(MeltTimeStep_Host(D3D1ConvPosition), -1.0);
                EXPECT_FLOAT_EQ(UndercoolingChange_Host(D3D1ConvPosition), RNorm);
                EXPECT_EQ(NumberOfSolidificationEvents_Host(D3D1ConvPosition), 1);
                EXPECT_EQ(SolidificationEventCounter(D3D1ConvPosition), 0);
                EXPECT_EQ(LayerID_Host(D3D1ConvPosition), 0);
                if (SimulationType == "C") {
                    // Directional solidification
                    // UndercoolingCurrent should be 0 for all cells
                    // CritTimeStep is equivalent to the Z coordinate of the cell
                    EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(D3D1ConvPosition), 0.0);
                    EXPECT_FLOAT_EQ(CritTimeStep_Host(D3D1ConvPosition), k);
                }
                else if (SimulationType == "SingleGrain") {
                    // Single grain
                    // UndercooliongCurrent depends on the undercooling set at the domain center and the thermal
                    // gradient CritTimeStep depends on the distance from the liquidus (all cells are already
                    // undercooled, so these values are negative)
                    EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(D3D1ConvPosition),
                                    InitUndercooling - GNorm * (k - kCenter));
                    EXPECT_FLOAT_EQ(CritTimeStep_Host(D3D1ConvPosition), -InitUndercooling + GNorm * (k - kCenter));
                }
                // Check that the data was copied into the time-temperature history view correctly
                EXPECT_FLOAT_EQ(static_cast<float>(MeltTimeStep_Host(D3D1ConvPosition)),
                                LayerTimeTempHistory_Host(D3D1ConvPosition, 0, 0));
                EXPECT_FLOAT_EQ(static_cast<float>(CritTimeStep_Host(D3D1ConvPosition)),
                                LayerTimeTempHistory_Host(D3D1ConvPosition, 0, 1));
                EXPECT_FLOAT_EQ(UndercoolingChange_Host(D3D1ConvPosition),
                                LayerTimeTempHistory_Host(D3D1ConvPosition, 0, 2));
            }
        }
    }
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, orientation_init_tests) {
    testOrientationInit_Vectors();
    testOrientationInit_Angles();
}
TEST(TEST_CATEGORY, temp_init_tests) {
    // Test for directional and single grain problems, and with and without a thermal gradient for the single grain
    // problem
    testTempInit_UnidirectionalGradient("C", 1000000);
    testTempInit_UnidirectionalGradient("SingleGrain", 0);
    testTempInit_UnidirectionalGradient("SingleGrain", 1000000);
}
} // end namespace Test
