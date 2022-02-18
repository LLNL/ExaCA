
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
    // There are 60 total cells in this domain (nx * ny * nz)
    // Even ranks have 30, and odd ranks have the other 30
    // The bottom 24 cells on each rank are assigned baseplate Grain ID values
    // The cells at Z = 4 are outside the "active" portion of the domain and are not assigned Grain IDs with the rest of
    // the baseplate
    double SubstrateGrainSpacing = 2.2; // This grain spacing ensures that there will be > 3 substrate grains
    ViewI_H GrainID_H(Kokkos::ViewAllocateWithoutInitializing("GrainID"), LocalDomainSize);
    // Initialize GrainID to 0 everywhere
    for (int i = 0; i < MyXSlices * MyYSlices * nz; i++) {
        GrainID_H(i) = 0;
    }
    SubstrateInit_FromGrainSpacing(SubstrateGrainSpacing, nx, ny, nz, nzActive, MyXSlices, MyYSlices, MyXOffset,
                                   MyYOffset, LocalActiveDomainSize, id, np, deltax, GrainID_H);

    // Check the results for baseplate grains
    if (id % 2 == 0) {
        // All cells that are part of the baseplate should have GrainIDs between 1 and 6, depending on the closest
        // baseplate grain center
        ViewI_H ExpectedGID(Kokkos::ViewAllocateWithoutInitializing("ExGrainID"), LocalDomainSize);
        ExpectedGID(0) = 3;
        ExpectedGID(1) = 1;
        for (int i = 2; i < 6; i++) {
            ExpectedGID(i) = 3;
        }
        ExpectedGID(6) = 5;
        ExpectedGID(7) = 1;
        for (int i = 8; i < 12; i++) {
            ExpectedGID(i) = 3;
        }
        for (int i = 12; i < 16; i++) {
            ExpectedGID(i) = 5;
        }
        ExpectedGID(16) = 3;
        ExpectedGID(17) = 3;
        ExpectedGID(18) = 5;
        for (int i = 19; i < 24; i++) {
            ExpectedGID(i) = 5;
        }
        // Powder layer cells (at domain edge, Y = 0) - should have GrainIDs > 6 (exact value depends on MPI rank id)
        // Cells that the code "thinks" are part of the ghost nodes (Y = 1) should not have been assigned a Grain ID
        ExpectedGID(24) = 7 + (id / 2) * 6;
        ExpectedGID(25) = 0;
        ExpectedGID(26) = 8 + (id / 2) * 6;
        ExpectedGID(27) = 0;
        ExpectedGID(28) = 9 + (id / 2) * 6;
        ExpectedGID(29) = 0;
        for (int i = 0; i < LocalDomainSize; i++) {
            EXPECT_EQ(GrainID_H(i), ExpectedGID(i));
        }
    }
    else {
        // All cells that are part of the baseplate should have GrainIDs between 1 and 6, depending on the closest
        // baseplate grain center
        ViewI_H ExpectedGID(Kokkos::ViewAllocateWithoutInitializing("ExGrainID"), LocalDomainSize);
        ExpectedGID(0) = 1;
        ExpectedGID(1) = 1;
        ExpectedGID(2) = 2;
        ExpectedGID(3) = 2;
        ExpectedGID(4) = 3;
        ExpectedGID(5) = 2;
        ExpectedGID(6) = 1;
        ExpectedGID(7) = 1;
        ExpectedGID(8) = 2;
        ExpectedGID(9) = 2;
        ExpectedGID(10) = 3;
        ExpectedGID(11) = 2;
        for (int i = 12; i < 16; i++) {
            ExpectedGID(i) = 4;
        }
        ExpectedGID(16) = 6;
        ExpectedGID(17) = 6;
        ExpectedGID(18) = 5;
        ExpectedGID(19) = 4;
        ExpectedGID(20) = 5;
        for (int i = 21; i < 24; i++) {
            ExpectedGID(i) = 6;
        }
        // Powder layer cells (at domain edge, Y = 1) - should have GrainIDs > 6 (exact value depends on MPI rank id)
        // Cells that the code "thinks" are part of the ghost nodes (Y = 0) should not have been assigned a Grain ID
        ExpectedGID(24) = 0;
        ExpectedGID(25) = 10 + (id / 2) * 6;
        ExpectedGID(26) = 0;
        ExpectedGID(27) = 11 + (id / 2) * 6;
        ExpectedGID(28) = 0;
        ExpectedGID(29) = 12 + (id / 2) * 6;
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
