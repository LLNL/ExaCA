
#include <Kokkos_Core.hpp>

#include <gtest/gtest.h>

namespace Test {
//---------------------------------------------------------------------------//
void testGhost1D() {
    // Setup test data.

    // Send/recv ghost neighbor data.
    // GhostNodes1D(0, 0, MyLeft, MyRight, MyXSlices, MyYSlices, MyXOffset, MyYOffset, NeighborX_G, NeighborY_G,
    // NeighborZ_G, CellType_G, DOCenter_G,GrainID_G, GrainUnitVector_G, GrainOrientation_G, DiagonalLength_G,
    // CritDiagonalLength_G, NGrainOrientations, BufferA, BufferB, BufferAR, BufferBR, BufSizeX,  BufSizeY, BufSizeZ,
    // Locks, ZBound_Low);

    // Check the results.
    // EXPECT_EQ( );
}

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, ghost_1d_test) { testGhost1D(); }

} // end namespace Test
