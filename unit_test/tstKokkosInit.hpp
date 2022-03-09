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

void testGrainInit() {

    // Tests GrainInit subroutine for 1D decomposition - initializing grain structure one layer of a multilayer domain
    int DecompositionStrategy = 1;
    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // Create test grid
    double deltax = 1.0 * pow(10, -5);
    int layernumber = -id; // this lets the RNG seed for nuclei locations be the same on every rank, since the seed is
                           // set as id + layernumber
    int nx = 1;
    int ny = 3 * np;
    int nz = 2;
    int ZBound_Low = 0;
    int ZBound_High = 0; // the k = 1 data belongs to a different layer - only modifying the k = 0 data
    int MyXSlices, MyYSlices, MyXOffset, MyYOffset, NeighborRank_North, NeighborRank_South, NeighborRank_East,
        NeighborRank_West, NeighborRank_NorthEast, NeighborRank_NorthWest, NeighborRank_SouthEast,
        NeighborRank_SouthWest, ProcessorsInXDirection, ProcessorsInYDirection;
    long int LocalDomainSize;

    // Perform grid decomposition
    DomainDecomposition(DecompositionStrategy, id, np, MyXSlices, MyYSlices, MyXOffset, MyYOffset, NeighborRank_North,
                        NeighborRank_South, NeighborRank_East, NeighborRank_West, NeighborRank_NorthEast,
                        NeighborRank_NorthWest, NeighborRank_SouthEast, NeighborRank_SouthWest, nx, ny, nz,
                        ProcessorsInXDirection, ProcessorsInYDirection, LocalDomainSize);
    int LocalActiveDomainSize = (ZBound_High - ZBound_Low + 1) * MyXSlices * MyYSlices;
    // LocalActiveDomainSize = 3 for rank 0 if only 1 MPI rank, will equal 4 otherwise
    // LocalActiveDomainSize = 5 if > 1 MPI ranks and 0 < id < np-1
    // LocalActiveDomainSize = 4 if > 1 MPI ranks and id = np-1

    // Initialize neighbor lists
    ViewI_H NeighborX(Kokkos::ViewAllocateWithoutInitializing("NeighborX"), 26);
    ViewI_H NeighborY(Kokkos::ViewAllocateWithoutInitializing("NeighborY"), 26);
    ViewI_H NeighborZ(Kokkos::ViewAllocateWithoutInitializing("NeighborZ"), 26);
    ViewI2D_H ItList(Kokkos::ViewAllocateWithoutInitializing("ItList"), 9, 26);
    NeighborListInit(NeighborX, NeighborY, NeighborZ, ItList);

    // Initialize an orientation for the grain being tested
    int NGrainOrientations = 1;
    ViewI_H GrainOrientation(Kokkos::ViewAllocateWithoutInitializing("GrainOrientation"), NGrainOrientations);
    GrainOrientation(0) = 0;
    ViewF_H GrainUnitVector(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), 9 * NGrainOrientations);
    GrainUnitVector(0) = 0.848294;  // x1
    GrainUnitVector(1) = 0.493303;  // y1
    GrainUnitVector(2) = 0.19248;   // z1
    GrainUnitVector(3) = -0.522525; // x2
    GrainUnitVector(4) = 0.720911;  // y2
    GrainUnitVector(5) = 0.455253;  // z2
    GrainUnitVector(6) = 0.0858167; // x3
    GrainUnitVector(7) = -0.486765; // y3
    GrainUnitVector(8) = 0.869308;  // z3

    // Initialize all cells as liquid, except for one on each rank (which is assigned the grain ID)
    ViewI_H CellType(Kokkos::ViewAllocateWithoutInitializing("CellType"), LocalDomainSize);
    ViewI_H GrainID(Kokkos::ViewAllocateWithoutInitializing("GrainID"), LocalDomainSize);
    for (int D3D1ConvPositionGlobal = 0; D3D1ConvPositionGlobal < LocalDomainSize; D3D1ConvPositionGlobal++) {
        CellType(D3D1ConvPositionGlobal) = Liquid;
        GrainID(D3D1ConvPositionGlobal) = 0;
    }
    CellType(1) = Active;
    GrainID(1) = 1;

    // Set Layer ID for the k = 0 cells to be equal to this layernumber, and set the k = 1 cells to something else
    ViewI_H LayerID(Kokkos::ViewAllocateWithoutInitializing("LayerID"), LocalDomainSize);
    for (int D3D1ConvPositionGlobal = 0; D3D1ConvPositionGlobal < LocalDomainSize; D3D1ConvPositionGlobal++) {
        int k = D3D1ConvPositionGlobal / (MyXSlices * MyYSlices);
        if (k == 0)
            LayerID(D3D1ConvPositionGlobal) = layernumber;
        else
            LayerID(D3D1ConvPositionGlobal) = -1;
    }

    double NMax = 2.0 * pow(10, 14); // Nucleation density and cell size set so each rank has a nucleated grain
    int PossibleNuclei_ThisRank, NextLayer_FirstNucleatedGrainID;

    // Active cell data structures
    ViewF_H DiagonalLength(Kokkos::ViewAllocateWithoutInitializing("DiagonalLength"), LocalActiveDomainSize);
    ViewF_H CritDiagonalLength(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength"),
                               26 * LocalActiveDomainSize);
    ViewF_H DOCenter(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * LocalActiveDomainSize);

    // Call GrainInit subroutine
    GrainInit(layernumber, NGrainOrientations, DecompositionStrategy, nx, ny, LocalActiveDomainSize, MyXSlices,
              MyYSlices, MyXOffset, MyYOffset, id, np, NeighborRank_North, NeighborRank_South, NeighborRank_East,
              NeighborRank_West, NeighborRank_NorthEast, NeighborRank_NorthWest, NeighborRank_SouthEast,
              NeighborRank_SouthWest, NeighborX, NeighborY, NeighborZ, GrainOrientation, GrainUnitVector,
              DiagonalLength, CellType, GrainID, CritDiagonalLength, DOCenter, deltax, NMax,
              NextLayer_FirstNucleatedGrainID, PossibleNuclei_ThisRank, ZBound_High, ZBound_Low, LayerID);

    // Check the results
    // Does the active cell on each rank have the correct diagonal length, octahedron center, and critical diagonal
    // lengths?
    EXPECT_FLOAT_EQ(DiagonalLength(1), 0.01); // initial octahedron diagonal length
    // Octahedron center should be at the origin of the active cell on the grid in x, y, z (0, 1, 0) plus 0.5 to each
    // coordinate, since the center is coincident with the cell center
    EXPECT_FLOAT_EQ(DOCenter(3), 0.5);             // X position of octahedron center
    EXPECT_FLOAT_EQ(DOCenter(4), MyYOffset + 1.5); // Y position of octahedron center
    EXPECT_FLOAT_EQ(DOCenter(5), 0.5);             // Z position of octahedron center
    // Based on the orientation of the octahedron and the distance to each neighboring cell on the grid, these critical
    // diagonal lengths were calculated to be the expected values required for the octahedron to engulf the centers of
    // each neighboring cells
    EXPECT_FLOAT_EQ(CritDiagonalLength(26), 1.7009771);
    EXPECT_FLOAT_EQ(CritDiagonalLength(27), 2.1710062);
    EXPECT_FLOAT_EQ(CritDiagonalLength(28), 1.9409304);
    EXPECT_FLOAT_EQ(CritDiagonalLength(29), 1.9225504);
    EXPECT_FLOAT_EQ(CritDiagonalLength(30), 2.2444904);
    EXPECT_FLOAT_EQ(CritDiagonalLength(31), 2.6762383);
    EXPECT_FLOAT_EQ(CritDiagonalLength(32), 2.7775407);
    EXPECT_FLOAT_EQ(CritDiagonalLength(33), 2.6560771);
    EXPECT_FLOAT_EQ(CritDiagonalLength(34), 2.1579251);
    EXPECT_FLOAT_EQ(CritDiagonalLength(35), 1.5170407);
    EXPECT_FLOAT_EQ(CritDiagonalLength(36), 1.5170407);
    EXPECT_FLOAT_EQ(CritDiagonalLength(37), 2.0631697);
    EXPECT_FLOAT_EQ(CritDiagonalLength(38), 2.4170816);
    EXPECT_FLOAT_EQ(CritDiagonalLength(39), 2.4170816);
    EXPECT_FLOAT_EQ(CritDiagonalLength(40), 2.0631697);
    EXPECT_FLOAT_EQ(CritDiagonalLength(41), 1.4566354);
    EXPECT_FLOAT_EQ(CritDiagonalLength(42), 1.4566354);
    EXPECT_FLOAT_EQ(CritDiagonalLength(43), 1.7009771);
    EXPECT_FLOAT_EQ(CritDiagonalLength(44), 1.9409304);
    EXPECT_FLOAT_EQ(CritDiagonalLength(45), 2.1710062);
    EXPECT_FLOAT_EQ(CritDiagonalLength(46), 2.2444904);
    EXPECT_FLOAT_EQ(CritDiagonalLength(47), 1.9225504);
    EXPECT_FLOAT_EQ(CritDiagonalLength(48), 2.6560771);
    EXPECT_FLOAT_EQ(CritDiagonalLength(49), 2.1579251);
    EXPECT_FLOAT_EQ(CritDiagonalLength(50), 2.6762383);
    EXPECT_FLOAT_EQ(CritDiagonalLength(51), 2.7775407);

    // Each rank had one nucleated grain initialized, along with
    // Does the nucleated grain on each rank exist in the right spot, based on the RNG seed?
    EXPECT_EQ(PossibleNuclei_ThisRank, 1);
    // On rank 0, since no ghost nodes exist at the first position, this was the spot assigned for the nucleated grain
    // On ranks 1+, since the first position is the halo region, which should NOT be assigned a nuclei location, and the
    // second position is an active cell (which is not an allowable nucleation site), the nuclei should've been assigned
    // the third position
    if (id == 0)
        EXPECT_EQ(CellType(0), TemporaryInit);
    else
        EXPECT_EQ(CellType(2), TemporaryInit);
}

void testNucleiInit() {

    // Tests NucleiInit subroutine for 1D decomposition - initializing 1 nuclei per rank, with the exception of rank 1,
    // which gets 2 nuclei
    int DecompositionStrategy = 1;
    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // Create test grid
    int layernumber = -id; // this lets the RNG seed for nuclei locations be the same on every rank, since the seed is
                           // set as id + layernumber
    int nx = 1;
    int ny = 3 * np;
    int nz = 2;
    int ZBound_Low = 0;
    int ZBound_High = 0; // the k = 1 data belongs to a different layer - only modifying the k = 0 data
    int MyXSlices, MyYSlices, MyXOffset, MyYOffset, NeighborRank_North, NeighborRank_South, NeighborRank_East,
        NeighborRank_West, NeighborRank_NorthEast, NeighborRank_NorthWest, NeighborRank_SouthEast,
        NeighborRank_SouthWest, ProcessorsInXDirection, ProcessorsInYDirection;
    long int LocalDomainSize;

    // Perform grid decomposition
    DomainDecomposition(DecompositionStrategy, id, np, MyXSlices, MyYSlices, MyXOffset, MyYOffset, NeighborRank_North,
                        NeighborRank_South, NeighborRank_East, NeighborRank_West, NeighborRank_NorthEast,
                        NeighborRank_NorthWest, NeighborRank_SouthEast, NeighborRank_SouthWest, nx, ny, nz,
                        ProcessorsInXDirection, ProcessorsInYDirection, LocalDomainSize);
    int LocalActiveDomainSize = (ZBound_High - ZBound_Low + 1) * MyXSlices * MyYSlices;
    // LocalActiveDomainSize = 3 for rank 0 if only 1 MPI rank, will equal 4 otherwise
    // LocalActiveDomainSize = 5 if > 1 MPI ranks and 0 < id < np-1
    // LocalActiveDomainSize = 4 if > 1 MPI ranks and id = np-1

    // Cell data structures (CellType should be initialized to 0s)
    ViewI_H CellType("CellType", LocalDomainSize);
    ViewI_H GrainID(Kokkos::ViewAllocateWithoutInitializing("GrainID"), LocalDomainSize);
    ViewI_H CritTimeStep(Kokkos::ViewAllocateWithoutInitializing("CritTimeStep"), LocalDomainSize);
    ViewF_H UndercoolingChange(Kokkos::ViewAllocateWithoutInitializing("UndercoolingChange"), LocalDomainSize);

    // Let CritTimeStep be equal to the cell's Y position relative to the overall domain bounds
    // Let UndercoolingChange be 0.5 for all ranks
    for (int i = 0; i < LocalActiveDomainSize; i++) {
        CritTimeStep(i) = i + MyYOffset;
        UndercoolingChange(i) = 0.5;
    }

    // TemporaryInit-type cell(s), with assigned nuclei grain ID - third cell on each rank
    CellType(2) = TemporaryInit;
    GrainID(2) = -id;

    // Nuclei data structures - one nuclei per rank, with one extra nucleus on rank 1 to testing sending/receiving of
    // data
    int PossibleNuclei_ThisRank = 1;
    if (id == 1)
        PossibleNuclei_ThisRank = 2;
    ViewI_H NucleiLocation(Kokkos::ViewAllocateWithoutInitializing("NucleiLocation"), PossibleNuclei_ThisRank);
    ViewI_H NucleationTimes(Kokkos::ViewAllocateWithoutInitializing("NucleationTimes"), PossibleNuclei_ThisRank);

    // Nucleation parameters
    double dTN = 5;
    double dTsigma = 0.5;

    // Call NucleiInit subroutine
    NucleiInit(layernumber, DecompositionStrategy, MyXSlices, MyYSlices, ZBound_Low, LocalActiveDomainSize, id, dTN,
               dTsigma, NeighborRank_North, NeighborRank_South, NeighborRank_East, NeighborRank_West,
               NeighborRank_NorthEast, NeighborRank_NorthWest, NeighborRank_SouthEast, NeighborRank_SouthWest,
               NucleiLocation, NucleationTimes, CellType, GrainID, CritTimeStep, UndercoolingChange);

    // Check the results - were the NucleiLocations initialized correctly, with the correct NucleationTimes?
    EXPECT_EQ(NucleiLocation(0), 2);
    EXPECT_EQ(NucleationTimes(0), MyYOffset + 12);
    if (id == 1) {
        EXPECT_EQ(NucleiLocation(1), 0);
        EXPECT_EQ(NucleationTimes(1), 12);
    }
}

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, grain_init_tests) {
    testSubstrateInit_FromGrainSpacing();
    testGrainInit();
    testNucleiInit();
}

} // end namespace Test
