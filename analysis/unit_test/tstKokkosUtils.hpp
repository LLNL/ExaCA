// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "GAutils.hpp"

#include <gtest/gtest.h>

#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
// GAutils tests with Kokkos
//---------------------------------------------------------------------------//
// For a 3D domain, obtain a vector of the grain IDs that appear in a representative portion
void testgetRepresentativeRegionGrainIDs() {

    const int nx = 5;
    const int ny = 4;
    const int nz = 3;

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

    // Representative region size
    const int XLow = 1;
    const int XHigh = 4;
    const int YLow = 2;
    const int YHigh = 3;
    const int ZLow = 0;
    const int ZHigh = 2;
    const int RepresentativeRegionSize = (ZHigh - ZLow + 1) * (XHigh - XLow + 1) * (YHigh - YLow + 1);
    // Number of cells from representative region at a given Z coordinate
    const int SliceSize = (XHigh - XLow + 1) * (YHigh - YLow + 1);

    // Fill vector
    std::vector<int> GrainIDVector =
        getRepresentativeRegionGrainIDs(GrainID, XLow, XHigh, YLow, YHigh, ZLow, ZHigh, RepresentativeRegionSize);

    // Check results
    for (int n = 0; n < SliceSize; n++)
        EXPECT_EQ(GrainIDVector[n], 0);
    for (int n = SliceSize; n < 2 * SliceSize; n++)
        EXPECT_EQ(GrainIDVector[n], 1);
    for (int n = 2 * SliceSize; n < 3 * SliceSize; n++)
        EXPECT_EQ(GrainIDVector[n], 2);
}

// For a 3D domain, obtain the extent in each grid direction of each unique grain ID
void testcalcGrainExtent() {

    const int nx = 5;
    const int ny = 4;
    const int nz = 3;
    const int XLow = 0;
    const int YLow = 0;
    const int ZLow = 0;
    const int XHigh = nx - 1;
    const int YHigh = ny - 1;
    const int ZHigh = nz - 1;
    double deltax = 2 * pow(10, -6);
    std::string RegionType = "volume";

    // View for storing grain ID data
    ViewI3D_H GrainID(Kokkos::ViewAllocateWithoutInitializing("GrainID"), nz, nx, ny);

    // 4 unique grains
    int NumberOfGrains = 4;
    std::vector<int> UniqueGrainIDVector = {1, 2, 3, 4};
    std::vector<float> GrainSizeVector(4, 0);
    for (int k = 0; k < nz - 1; k++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 2; j++) {
                GrainID(k, i, j) = 1;
                GrainSizeVector[0] += convertToMicrons(deltax, RegionType);
            }
        }
    }
    for (int k = 0; k < nz - 1; k++) {
        for (int i = 3; i < nx; i++) {
            for (int j = 0; j < 2; j++) {
                GrainID(k, i, j) = 2;
                GrainSizeVector[1] += convertToMicrons(deltax, RegionType);
            }
        }
    }
    for (int k = 0; k < nz - 1; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 2; j < ny; j++) {
                GrainID(k, i, j) = 3;
                GrainSizeVector[2] += convertToMicrons(deltax, RegionType);
            }
        }
    }
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            GrainID(nz - 1, i, j) = 4;
            GrainSizeVector[3] += convertToMicrons(deltax, RegionType);
        }
    }

    // Convert grain sizes to volumes
    // Extent of each grain in X
    std::vector<float> GrainExtentX(NumberOfGrains);
    calcGrainExtent(GrainExtentX, GrainID, UniqueGrainIDVector, GrainSizeVector, NumberOfGrains, XLow, XHigh, YLow,
                    YHigh, ZLow, ZHigh, "X", deltax, RegionType);
    // Extent of each grain in Y
    std::vector<float> GrainExtentY(NumberOfGrains);
    calcGrainExtent(GrainExtentY, GrainID, UniqueGrainIDVector, GrainSizeVector, NumberOfGrains, XLow, XHigh, YLow,
                    YHigh, ZLow, ZHigh, "Y", deltax, RegionType);
    // Extent of each grain in Z
    std::vector<float> GrainExtentZ(NumberOfGrains);
    calcGrainExtent(GrainExtentZ, GrainID, UniqueGrainIDVector, GrainSizeVector, NumberOfGrains, XLow, XHigh, YLow,
                    YHigh, ZLow, ZHigh, "Z", deltax, RegionType);

    // Check results
    std::vector<float> GrainExtentX_Expected = {3, 2, 5, 5}; // in cells
    std::vector<float> GrainExtentY_Expected = {2, 2, 2, 4};
    std::vector<float> GrainExtentZ_Expected = {2, 2, 2, 1};
    for (int n = 0; n < NumberOfGrains; n++) {
        // Expected value should be in microns
        float GrainExtentX_Expected_Microns = GrainExtentX_Expected[n] * convertToMicrons(deltax, "length");
        float GrainExtentY_Expected_Microns = GrainExtentY_Expected[n] * convertToMicrons(deltax, "length");
        float GrainExtentZ_Expected_Microns = GrainExtentZ_Expected[n] * convertToMicrons(deltax, "length");
        EXPECT_FLOAT_EQ(GrainExtentX[n], GrainExtentX_Expected_Microns);
        EXPECT_FLOAT_EQ(GrainExtentY[n], GrainExtentY_Expected_Microns);
        EXPECT_FLOAT_EQ(GrainExtentZ[n], GrainExtentZ_Expected_Microns);
    }
}

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, grain_util_tests_kokkos) {
    testgetRepresentativeRegionGrainIDs();
    testcalcGrainExtent();
}
} // end namespace Test
