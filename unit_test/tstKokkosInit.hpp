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

// Tests calcZBound_Low_NoRemelt and calcZBound_High_NoRemelt
void testcalcZBounds_NoRemelt() {

    // Domain setup
    int nx = 5;
    int ny = 4;
    int nz = 12;
    int DomainSize = nx * ny * nz;
    ViewI_H LayerID_Host(Kokkos::ViewAllocateWithoutInitializing("LayerID_Host"), DomainSize);

    int layernumber = 0;
    // Fill LayerID with -1s (not associated with a layer), assign a few cells the value "layernumber", assign a few
    // other cells the value "layernumber + 1"
    Kokkos::deep_copy(LayerID_Host, -1);
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                int Coordinate1D = k * nx * ny + i * ny + j;
                if ((i == 2) && (j == 2) && (k > 4)) {
                    if (k < 8)
                        LayerID_Host(Coordinate1D) = layernumber;
                    else
                        LayerID_Host(Coordinate1D) = layernumber + 1;
                }
            }
        }
    }

    // Copy view to device
    ViewI LayerID = Kokkos::create_mirror_view_and_copy(device_memory_space(), LayerID_Host);

    // Calculate ZBound_Low and ZBound_High
    int ZBound_Low = calcZBound_Low_NoRemelt(0, nx, ny, DomainSize, layernumber, LayerID);
    int ZBound_High = calcZBound_High_NoRemelt(0, nx, ny, DomainSize, layernumber, LayerID);

    // Check against expected values
    EXPECT_EQ(ZBound_Low, 5);
    EXPECT_EQ(ZBound_High, 7);
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, orientation_init_tests) {
    testOrientationInit_Vectors();
    testOrientationInit_Angles();
}
TEST(TEST_CATEGORY, activedomainsizecalc) { testcalcZBounds_NoRemelt(); }
} // end namespace Test
