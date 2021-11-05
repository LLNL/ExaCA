
#include <Kokkos_Core.hpp>

#include "CAinitialize.hpp"
#include "CAtypes.hpp"

#include <gtest/gtest.h>

#include "mpi.h"

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
void testSubstrateInit_FromGrainSpacing() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Create test data
    int nz = 5;
    int nzActive = 4;
    int nx = 3;
    int MyXSlices = 3;
    int MyXOffset = 0;
    // Each rank is assigned a different portion of the domain in Y
    int ny = 4;
    int MyYSlices = 2;
    int MyYOffset;
    if (id % 2 == 0)
        MyYOffset = 0;
    else
        MyYOffset = 2;
    double deltax = 1 * pow(10, -6);
    int LocalActiveDomainSize = MyXSlices * MyYSlices * nzActive;
    int LocalDomainSize = MyXSlices * MyYSlices * nz;
    // There are 48 total cells in this domain (nx * ny * nz)
    // Even ranks have 30, and odd ranks have the other 30
    double SubstrateGrainSpacing = 1.95; // This grain spacing ensures that there will be 3 substrate grains
    ViewI_H GrainID_H(Kokkos::ViewAllocateWithoutInitializing("GrainID"), LocalDomainSize);
    ViewI_H CritTimeStep_H(Kokkos::ViewAllocateWithoutInitializing("CritTimeStep"), LocalDomainSize);
    // Initialize GrainID to 0 everywhere
    // Initialize CritTimeStep to 0 everywhere to allow substrate generation for these cells
    for (int i = 0; i < MyXSlices * MyYSlices * nz; i++) {
        CritTimeStep_H(i) = 0;
        GrainID_H(i) = 0;
    }
    CritTimeStep_H(13) = 1; // Let 1 cell have a non-zero CritTimeStep - this cell should end up being assigned a
                            // GrainID of 0, while the others have GrainIDs of 1, 2, or 3
    SubstrateInit_FromGrainSpacing(SubstrateGrainSpacing, nx, ny, nz, nzActive, MyXSlices, MyYSlices, MyXOffset,
                                   MyYOffset, LocalActiveDomainSize, id, np, deltax, GrainID_H, CritTimeStep_H);

    // Check the results
    if (id % 2 == 0) {
        // Expected GrainID values: Cells 0-5 are walls, should have Grain ID 0
        // Cell 13 has CritTimeStep > 1, should also have GrainID 0
        // Other cells should have GrainID 1, 2, or 3, depending on the closest grain center
        ViewI_H ExpectedGID(Kokkos::ViewAllocateWithoutInitializing("ExGrainID"), LocalDomainSize);
        for (int i = 0; i < 6; i++) {
            ExpectedGID(i) = 0;
        }
        for (int i = 6; i < 19; i++) {
            ExpectedGID(i) = 1;
        }
        ExpectedGID(13) = 0;
        ExpectedGID(19) = 2;
        ExpectedGID(20) = 1;
        ExpectedGID(21) = 2;
        ExpectedGID(22) = 1;
        ExpectedGID(23) = 2;
        for (int i = 24; i < 30; i++) {
            ExpectedGID(i) = 3;
        }
        for (int i = 0; i < LocalDomainSize; i++) {
            EXPECT_EQ(GrainID_H(i), ExpectedGID(i));
        }
    }
    else {
        // Expected GrainID values: Cells 0-5 are walls, should have Grain ID 0
        // Cell 13 has CritTimeStep > 1, should also have GrainID 0
        // Other cells should have GrainID 1, 2, or 3, depending on the closest grain center
        ViewI_H ExpectedGID(Kokkos::ViewAllocateWithoutInitializing("ExGrainID"), LocalDomainSize);
        for (int i = 0; i < 6; i++) {
            ExpectedGID(i) = 0;
        }
        ExpectedGID(6) = 1;
        for (int i = 7; i < 24; i++) {
            ExpectedGID(i) = 2;
        }
        ExpectedGID(8) = 1;
        ExpectedGID(10) = 1;
        ExpectedGID(13) = 0;
        for (int i = 24; i < 30; i++) {
            ExpectedGID(i) = 3;
        }
        for (int i = 0; i < LocalDomainSize; i++) {
            EXPECT_EQ(GrainID_H(i), ExpectedGID(i));
        }
    }
}

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, substrate_init_test) { testSubstrateInit_FromGrainSpacing(); }

} // end namespace Test
