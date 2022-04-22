#include <Kokkos_Core.hpp>

#include "CAfunctions.hpp"
#include "CAinitialize.hpp"
#include "CAtypes.hpp"

#include <gtest/gtest.h>

#include "mpi.h"

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
void testSubstrateInit_ConstrainedGrowth() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // Create test data
    int nz = 2;
    int nzActive = 2;
    int nx = 1;
    int MyXSlices = 1;
    int MyXOffset = 0;
    // Domain size in Y depends on the number of ranks - each rank has 4 cells in Y
    // Each rank is assigned a different portion of the domain in Y
    int ny = 4 * np;
    int MyYSlices = 4;
    int MyYOffset = 4 * id;
    int LocalActiveDomainSize = MyXSlices * MyYSlices * nzActive;
    int LocalDomainSize = MyXSlices * MyYSlices * nz;

    double FractSurfaceSitesActive = 0.5; // Each rank will have 2 active cells each, on average
    double RNGSeed = 0.0;
    // Initialize grain orientations
    std::string GrainOrientationFile = "GrainOrientationVectors_Robert.csv";
    int NGrainOrientations = 10000; // Number of grain orientations considered in the simulation
    ViewF_H GrainUnitVector_Host(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector_Host"),
                                 9 * NGrainOrientations);
    OrientationInit(id, NGrainOrientations, GrainUnitVector_Host, GrainOrientationFile);
    using memory_space = Kokkos::DefaultExecutionSpace::memory_space;
    ViewF GrainUnitVector = Kokkos::create_mirror_view_and_copy(memory_space(), GrainUnitVector_Host);

    // Initialize neighbor lists
    ViewI_H NeighborX_Host(Kokkos::ViewAllocateWithoutInitializing("NeighborX_Host"), 26);
    ViewI_H NeighborY_Host(Kokkos::ViewAllocateWithoutInitializing("NeighborY_Host"), 26);
    ViewI_H NeighborZ_Host(Kokkos::ViewAllocateWithoutInitializing("NeighborZ_Host"), 26);
    NeighborListInit(NeighborX_Host, NeighborY_Host, NeighborZ_Host);
    ViewI NeighborX = Kokkos::create_mirror_view_and_copy(memory_space(), NeighborX_Host);
    ViewI NeighborY = Kokkos::create_mirror_view_and_copy(memory_space(), NeighborY_Host);
    ViewI NeighborZ = Kokkos::create_mirror_view_and_copy(memory_space(), NeighborZ_Host);

    // Initialize views - set initial GrainID values to 0, all CellType values to liquid
    ViewI CellType(Kokkos::ViewAllocateWithoutInitializing("CellType"), LocalDomainSize);
    Kokkos::deep_copy(CellType, Liquid);
    ViewI GrainID("GrainID", LocalDomainSize);
    ViewF DiagonalLength(Kokkos::ViewAllocateWithoutInitializing("DiagonalLength"), LocalActiveDomainSize);
    ViewF DOCenter(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * LocalActiveDomainSize);
    ViewF CritDiagonalLength(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength"), 26 * LocalActiveDomainSize);

    SubstrateInit_ConstrainedGrowth(id, FractSurfaceSitesActive, MyXSlices, MyYSlices, nx, ny, MyXOffset, MyYOffset,
                                    NeighborX, NeighborY, NeighborZ, GrainUnitVector, NGrainOrientations, CellType,
                                    GrainID, DiagonalLength, DOCenter, CritDiagonalLength, RNGSeed);

    // Copy CellType, GrainID views to host to check values
    ViewI CellType_Host(Kokkos::ViewAllocateWithoutInitializing("CellType_Host"), LocalDomainSize);
    ViewI GrainID_Host(Kokkos::ViewAllocateWithoutInitializing("GrainID_Host"), LocalDomainSize);
    Kokkos::deep_copy(GrainID_Host, GrainID);
    Kokkos::deep_copy(CellType_Host, CellType);
    for (int i = 0; i < LocalDomainSize; i++) {
        if (i >= MyXSlices * MyYSlices) {
            // Not at bottom surface - should be liquid cells with GrainID still equal to 0
            EXPECT_EQ(GrainID_Host(i), 0);
            EXPECT_EQ(CellType_Host(i), Liquid);
        }
        else {
            // Check that active cells have GrainIDs > 0, and less than 2 * np + 1 (there are 2 * np different positive
            // GrainIDs used for epitaxial grain seeds)
            if (CellType_Host(i) == Active) {
                EXPECT_GT(GrainID_Host(i), 0);
                EXPECT_LT(GrainID_Host(i), 2 * np + 1);
            }
            else {
                // Liquid cells should still have GrainID = 0
                EXPECT_EQ(GrainID_Host(i), 0);
            }
        }
    }
}

void testBaseplateInit_FromGrainSpacing() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Create test data
    int nz = 4;
    float ZMinLayer[1];
    float ZMaxLayer[1];
    ZMinLayer[0] = 0;
    ZMaxLayer[0] = 2 * pow(10, -6);
    int nx = 3;
    int MyXSlices = 3;
    int MyXOffset = 0;
    // Each rank is assigned a different portion of the domain in Y
    int ny = 3 * np;
    int MyYSlices = 3;
    int MyYOffset = 3 * id;
    double deltax = 1 * pow(10, -6);
    int BaseplateSize = MyXSlices * MyYSlices * (round((ZMaxLayer[0] - ZMinLayer[0]) / deltax) + 1);
    int LocalDomainSize = MyXSlices * MyYSlices * nz;
    // There are 36 * np total cells in this domain (nx * ny * nz)
    // Each rank has 36 cells - the bottom 27 cells are assigned baseplate Grain ID values
    // The top cells (Z = 4) are outside the "active" portion of the domain and are not assigned Grain IDs with the rest
    // of the baseplate
    double SubstrateGrainSpacing =
        3.0; // This grain spacing ensures that there will be 1 grain per number of MPI ranks present
    double RNGSeed = 0.0;
    int NextLayer_FirstEpitaxialGrainID;
    // Initialize GrainIDs to 0 on device
    ViewI GrainID("GrainID_Device", LocalDomainSize);

    BaseplateInit_FromGrainSpacing(SubstrateGrainSpacing, nx, ny, ZMinLayer, ZMaxLayer, MyXSlices, MyYSlices, MyXOffset,
                                   MyYOffset, id, deltax, GrainID, RNGSeed, NextLayer_FirstEpitaxialGrainID);

    // Copy results back to host to check
    ViewI_H GrainID_H("GrainID_Host", LocalDomainSize);
    Kokkos::deep_copy(GrainID_H, GrainID);

    // Check the results for baseplate grains - cells should have GrainIDs between 1 and np (inclusive) if they are part
    // of the active domain, or 0 (unassigned) if not part of the active domain
    for (int i = 0; i < BaseplateSize; i++) {
        EXPECT_GT(GrainID_H(i), 0);
        EXPECT_LT(GrainID_H(i), np + 1);
    }
    for (int i = BaseplateSize; i < LocalDomainSize; i++) {
        EXPECT_EQ(GrainID_H(i), 0);
    }
    // Next unused GrainID should be the number of grains present in the baseplate plus 2 (since GrainID = 0 is not used
    // for any baseplate grains)
    EXPECT_EQ(NextLayer_FirstEpitaxialGrainID, np + 2);
}

void testPowderInit() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // Create test data
    int nz = 6;
    int layernumber = 1;
    // Z = 4 and Z = 5 should be the region seeded with powder layer Grain IDs
    float ZMaxLayer[2];
    ZMaxLayer[1] = 5.0 * pow(10, -6);
    double deltax = 1.0 * pow(10, -6);
    float ZMin = 0;
    int LayerHeight = 2;
    int nx = 1;
    int MyXSlices = 1;
    int MyXOffset = 0;
    // Each rank is assigned a different portion of the domain in Y
    int ny = 2 * np;
    int MyYSlices = 2;
    int MyYOffset = 2 * id;
    int LocalDomainSize = MyXSlices * MyYSlices * nz;

    // Initialize GrainIDs to 0 on device
    ViewI GrainID("GrainID_Device", LocalDomainSize);

    // Starting Grain ID for powder grains
    int NextLayer_FirstEpitaxialGrainID = 10;
    int PreviousLayer_FirstEpitaxialGrainID = NextLayer_FirstEpitaxialGrainID;

    // Seed used to shuffle powder layer grain IDs
    double RNGSeed = 0.0;

    PowderInit(layernumber, nx, ny, LayerHeight, ZMaxLayer, ZMin, deltax, MyXSlices, MyYSlices, MyXOffset, MyYOffset,
               id, GrainID, RNGSeed, NextLayer_FirstEpitaxialGrainID);

    // Copy results back to host to check
    ViewI_H GrainID_H("GrainID_Host", LocalDomainSize);
    Kokkos::deep_copy(GrainID_H, GrainID);

    // Check the results - were the right number of powder grain IDs used?
    int Expected_NextLayer_FirstEpitaxialGrainID = PreviousLayer_FirstEpitaxialGrainID + nx * ny * LayerHeight;
    EXPECT_EQ(NextLayer_FirstEpitaxialGrainID, Expected_NextLayer_FirstEpitaxialGrainID);

    // Check the results - powder grains should have unique Grain ID values larger than
    // PreviousLayer_FirstEpitaxialGrainID-1 and smaller than NextLayer_FirstEpitaxialGrainID. Other cells should still
    // have Grain IDs of 0
    int TopPowderLayer = MyXSlices * MyYSlices * (nz - LayerHeight);
    for (int i = 0; i < TopPowderLayer; i++) {
        EXPECT_EQ(GrainID_H(i), 0);
    }
    for (int i = TopPowderLayer; i < LocalDomainSize; i++) {
        EXPECT_GT(GrainID_H(i), PreviousLayer_FirstEpitaxialGrainID - 1);
        EXPECT_LT(GrainID_H(i), NextLayer_FirstEpitaxialGrainID);
    }
}
void testCellTypeInit() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    int layernumber = 0;
    // Create test grid
    int nz = 7;
    // Each MPI rank has a different portion of the domain in the Z direction "active" for initialization, starting at Z
    // = 5
    int ZBound_High = 5;
    int ZBound_Low = ZBound_High - (5 - (id % nz));
    int nzActive = ZBound_High - ZBound_Low + 1;
    // Each rank is assigned the same X and Y regions
    int MyXSlices = 7;
    int MyXOffset = 0;
    // Each rank is assigned a different portion of the domain in Y
    int MyYSlices = 4;
    int MyYOffset = 0;
    int LocalActiveDomainSize = MyXSlices * MyYSlices * nzActive;
    int LocalDomainSize = MyXSlices * MyYSlices * nz;
    int DecompositionStrategy = 1;
    // MPI rank locations relative to the global grid
    bool AtNorthBoundary, AtSouthBoundary;
    bool AtEastBoundary = false;
    bool AtWestBoundary = false;
    if (id == 0)
        AtSouthBoundary = true;
    else
        AtSouthBoundary = false;
    if (id == np - 1)
        AtNorthBoundary = true;
    else
        AtNorthBoundary = false;

    // Initialize neighbor lists
    ViewI_H NeighborX_Host(Kokkos::ViewAllocateWithoutInitializing("NeighborX_Host"), 26);
    ViewI_H NeighborY_Host(Kokkos::ViewAllocateWithoutInitializing("NeighborY_Host"), 26);
    ViewI_H NeighborZ_Host(Kokkos::ViewAllocateWithoutInitializing("NeighborZ_Host"), 26);
    NeighborListInit(NeighborX_Host, NeighborY_Host, NeighborZ_Host);
    using memory_space = Kokkos::DefaultExecutionSpace::memory_space;
    ViewI NeighborX = Kokkos::create_mirror_view_and_copy(memory_space(), NeighborX_Host);
    ViewI NeighborY = Kokkos::create_mirror_view_and_copy(memory_space(), NeighborY_Host);
    ViewI NeighborZ = Kokkos::create_mirror_view_and_copy(memory_space(), NeighborZ_Host);

    // Initialize test CellType, CritTimeStep, GrainID, LayerID data
    ViewI_H CritTimeStep_Host(Kokkos::ViewAllocateWithoutInitializing("CritTimeStep_Host"), LocalDomainSize);
    ViewI_H GrainID_Host(Kokkos::ViewAllocateWithoutInitializing("GrainID_Host"), LocalDomainSize);
    ViewI_H LayerID_Host(Kokkos::ViewAllocateWithoutInitializing("LayerID_Host"), LocalDomainSize);
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < MyXSlices; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                int D3D1ConvPositionGlobal = k * MyXSlices * MyYSlices + i * MyYSlices + j;
                GrainID_Host(D3D1ConvPositionGlobal) = 1;
                // Let the top portion of the cells be part of a different layernumber
                if (k == nz - 1)
                    LayerID_Host(D3D1ConvPositionGlobal) = 1;
                else
                    LayerID_Host(D3D1ConvPositionGlobal) = 0;
                // Bottom left portion of the cells have no temperature data
                if (i + k <= 5)
                    CritTimeStep_Host(D3D1ConvPositionGlobal) = 0;
                else
                    CritTimeStep_Host(D3D1ConvPositionGlobal) = 1;
            }
        }
    }
    ViewI CritTimeStep = Kokkos::create_mirror_view_and_copy(memory_space(), CritTimeStep_Host);
    ViewI GrainID = Kokkos::create_mirror_view_and_copy(memory_space(), GrainID_Host);
    ViewI LayerID = Kokkos::create_mirror_view_and_copy(memory_space(), LayerID_Host);

    // Initialize an orientation for the grain being tested
    int NGrainOrientations = 1;
    ViewF_H GrainUnitVector_Host(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector_Host"),
                                 9 * NGrainOrientations);
    GrainUnitVector_Host(0) = 0.848294;  // x1
    GrainUnitVector_Host(1) = 0.493303;  // y1
    GrainUnitVector_Host(2) = 0.19248;   // z1
    GrainUnitVector_Host(3) = -0.522525; // x2
    GrainUnitVector_Host(4) = 0.720911;  // y2
    GrainUnitVector_Host(5) = 0.455253;  // z2
    GrainUnitVector_Host(6) = 0.0858167; // x3
    GrainUnitVector_Host(7) = -0.486765; // y3
    GrainUnitVector_Host(8) = 0.869308;  // z3
    ViewF GrainUnitVector = Kokkos::create_mirror_view_and_copy(memory_space(), GrainUnitVector_Host);

    // Active cell data structures
    ViewF DiagonalLength(Kokkos::ViewAllocateWithoutInitializing("DiagonalLength"), LocalActiveDomainSize);
    ViewF CritDiagonalLength(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength"), 26 * LocalActiveDomainSize);
    ViewF DOCenter(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * LocalActiveDomainSize);

    // Cell types to be initialized
    ViewI CellType(Kokkos::ViewAllocateWithoutInitializing("CellType"), LocalDomainSize);

    // Buffers for ghost node data (fixed size)
    int BufSizeX = MyXSlices;
    int BufSizeY = 0;

    // Send buffers for ghost node data should be initialized with zeros
    Buffer2D BufferSouthSend("BufferSouthSend", BufSizeX * nzActive, 5);
    Buffer2D BufferNorthSend("BufferNorthSend", BufSizeX * nzActive, 5);
    Buffer2D BufferEastSend("BufferEastSend", BufSizeY * nzActive, 5);
    Buffer2D BufferWestSend("BufferWestSend", BufSizeY * nzActive, 5);
    Buffer2D BufferNorthEastSend("BufferNorthEastSend", nzActive, 5);
    Buffer2D BufferNorthWestSend("BufferNorthWestSend", nzActive, 5);
    Buffer2D BufferSouthEastSend("BufferSouthEastSend", nzActive, 5);
    Buffer2D BufferSouthWestSend("BufferSouthWestSend", nzActive, 5);

    // Initialize cell types and active cell data structures
    CellTypeInit(layernumber, id, np, DecompositionStrategy, MyXSlices, MyYSlices, MyXOffset, MyYOffset, ZBound_Low, nz,
                 LocalActiveDomainSize, LocalDomainSize, CellType, CritTimeStep, NeighborX, NeighborY, NeighborZ,
                 NGrainOrientations, GrainUnitVector, DiagonalLength, GrainID, CritDiagonalLength, DOCenter, LayerID,
                 BufferWestSend, BufferEastSend, BufferNorthSend, BufferSouthSend, BufferNorthEastSend,
                 BufferNorthWestSend, BufferSouthEastSend, BufferSouthWestSend, BufSizeX, BufSizeY, AtNorthBoundary,
                 AtSouthBoundary, AtEastBoundary, AtWestBoundary);

    // Copy views back to host to check the results
    ViewF_H DiagonalLength_Host(Kokkos::ViewAllocateWithoutInitializing("DiagonalLength_Host"), LocalActiveDomainSize);
    ViewF_H CritDiagonalLength_Host(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength_Host"),
                                    26 * LocalActiveDomainSize);
    ViewF_H DOCenter_Host(Kokkos::ViewAllocateWithoutInitializing("DOCenter_Host"), 3 * LocalActiveDomainSize);
    ViewI_H CellType_Host(Kokkos::ViewAllocateWithoutInitializing("CellType_Host"), LocalDomainSize);
    Kokkos::deep_copy(DiagonalLength_Host, DiagonalLength);
    Kokkos::deep_copy(CritDiagonalLength_Host, CritDiagonalLength);
    Kokkos::deep_copy(DOCenter_Host, DOCenter);
    Kokkos::deep_copy(CellType_Host, CellType);

    // Solid cells where no temperature data existed
    // Active cells separate solid-liquid cells, as well as at the active region bottom (implicitly borders a solid
    // cell, since the layers below should be solid) Liquid cells are all other cells
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < MyXSlices; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                int D3D1ConvPositionGlobal = k * MyXSlices * MyYSlices + i * MyYSlices + j;
                if (i + k <= 5) {
                    EXPECT_EQ(CellType_Host(D3D1ConvPositionGlobal), Solid);
                }
                else if (i + k <= 7) {
                    EXPECT_EQ(CellType_Host(D3D1ConvPositionGlobal), Active);
                    // Check that active cell data structures were initialized properly for cells in the active portion
                    // of the domain
                    if ((k >= ZBound_Low) && (k <= ZBound_High)) {
                        int D3D1ConvPosition = (k - ZBound_Low) * MyXSlices * MyYSlices + i * MyYSlices + j;
                        EXPECT_FLOAT_EQ(DiagonalLength_Host(D3D1ConvPosition),
                                        0.01); // initial octahedron diagonal length
                        // Octahedron center should be at the origin of the active cell on the grid in x, y, z (0, 1, 0)
                        // plus 0.5 to each coordinate, since the center is coincident with the cell center
                        EXPECT_FLOAT_EQ(DOCenter_Host(3 * D3D1ConvPosition),
                                        i + MyXOffset + 0.5); // X position of octahedron center
                        EXPECT_FLOAT_EQ(DOCenter_Host(3 * D3D1ConvPosition + 1),
                                        j + MyYOffset + 0.5); // Y position of octahedron center
                        EXPECT_FLOAT_EQ(DOCenter_Host(3 * D3D1ConvPosition + 2),
                                        k + 0.5); // Z position of octahedron center
                        // Based on the orientation of the octahedron and the distance to each neighboring cell on the
                        // grid, these critical diagonal lengths were calculated to be the expected values required for
                        // the octahedron to engulf the centers of each neighboring cells
                        std::vector<float> CritDiagonalLength_Expected{
                            1.7009771, 2.1710062, 1.9409304, 1.9225504, 2.2444904, 2.6762383, 2.7775407,
                            2.6560771, 2.1579251, 1.5170407, 1.5170407, 2.0631697, 2.4170816, 2.4170816,
                            2.0631697, 1.4566354, 1.4566354, 1.7009771, 1.9409304, 2.1710062, 2.2444904,
                            1.9225504, 2.6560771, 2.1579251, 2.6762383, 2.7775407};
                        for (int n = 0; n < 26; n++) {
                            EXPECT_FLOAT_EQ(CritDiagonalLength_Host(26 * D3D1ConvPosition + n),
                                            CritDiagonalLength_Expected[n]);
                        }
                    }
                }
                else {
                    EXPECT_EQ(CellType_Host(D3D1ConvPositionGlobal), Liquid);
                }
            }
        }
    }

    // Further check that active cell data was properly loaded into send buffers, and that locations in the send buffers
    // not corresponding to active cells were left alone (should still be 0s)
    for (int k = ZBound_Low; k <= ZBound_High; k++) {
        for (int i = 0; i < MyXSlices; i++) {
            int GNPosition = (k - ZBound_Low) * BufSizeX + i; // Position of cell in buffer
            // Check the south buffer - Data being sent to the "south" (BufferSouthSend) is from active cells at Y = 1
            int D3D1ConvPositionGlobal_South =
                k * MyXSlices * MyYSlices + i * MyYSlices + 1; // Position of cell on grid
            if ((CellType_Host(D3D1ConvPositionGlobal_South) == Active) && (!(AtSouthBoundary))) {
                EXPECT_FLOAT_EQ(BufferSouthSend(GNPosition, 0), 1);
                EXPECT_FLOAT_EQ(BufferSouthSend(GNPosition, 1), i + MyXOffset + 0.5);
                EXPECT_FLOAT_EQ(BufferSouthSend(GNPosition, 2), 1 + MyYOffset + 0.5);
                EXPECT_FLOAT_EQ(BufferSouthSend(GNPosition, 3), k + 0.5);
                EXPECT_FLOAT_EQ(BufferSouthSend(GNPosition, 4), 0.01);
            }
            else {
                EXPECT_FLOAT_EQ(BufferSouthSend(GNPosition, 0), 0.0);
                EXPECT_FLOAT_EQ(BufferSouthSend(GNPosition, 1), 0.0);
                EXPECT_FLOAT_EQ(BufferSouthSend(GNPosition, 2), 0.0);
                EXPECT_FLOAT_EQ(BufferSouthSend(GNPosition, 3), 0.0);
                EXPECT_FLOAT_EQ(BufferSouthSend(GNPosition, 4), 0.0);
            }
            // Check the north buffer - Data being sent to the "north" (BufferNorthSend) is from active cells at Y = 2
            int D3D1ConvPositionGlobal_North = k * MyXSlices * MyYSlices + i * MyYSlices + 2;
            if ((CellType_Host(D3D1ConvPositionGlobal_North) == Active) && (!(AtNorthBoundary))) {
                EXPECT_FLOAT_EQ(BufferNorthSend(GNPosition, 1), i + MyXOffset + 0.5);
                EXPECT_FLOAT_EQ(BufferNorthSend(GNPosition, 2), 2 + MyYOffset + 0.5);
                EXPECT_FLOAT_EQ(BufferNorthSend(GNPosition, 3), k + 0.5);
                EXPECT_FLOAT_EQ(BufferNorthSend(GNPosition, 4), 0.01);
            }
            else {
                EXPECT_FLOAT_EQ(BufferNorthSend(GNPosition, 0), 0.0);
                EXPECT_FLOAT_EQ(BufferNorthSend(GNPosition, 1), 0.0);
                EXPECT_FLOAT_EQ(BufferNorthSend(GNPosition, 2), 0.0);
                EXPECT_FLOAT_EQ(BufferNorthSend(GNPosition, 3), 0.0);
                EXPECT_FLOAT_EQ(BufferNorthSend(GNPosition, 4), 0.0);
            }
        }
    }
}

void testNucleiInit() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Create test data
    int nz = 5;
    // Only Z = 1 through 4 is part of this layer
    int nzActive = 4;
    int ZBound_Low = 1;
    int layernumber = 1;
    int nx = 4;
    int MyXSlices = 4;
    int MyXOffset = 0;
    // Each rank is assigned a different portion of the domain in Y
    int ny = 2 * np;
    int MyYSlices = 2;
    int MyYOffset = 2 * id;
    double deltax = 1;
    int LocalActiveDomainSize = MyXSlices * MyYSlices * nzActive;
    int LocalDomainSize = MyXSlices * MyYSlices * nz;
    // MPI rank locations relative to the global grid
    bool AtNorthBoundary, AtSouthBoundary;
    bool AtEastBoundary = true;
    bool AtWestBoundary = true;
    if (id == 0)
        AtSouthBoundary = true;
    else
        AtSouthBoundary = false;
    if (id == np - 1)
        AtNorthBoundary = true;
    else
        AtNorthBoundary = false;

    // There are 40 * np total cells in this domain (nx * ny * nz)
    // Each rank has 40 cells - the top 32 cells are part of the active layer and are candidates for nucleation
    // assignment
    double NMax = 0.125; // This nucleation density ensures there will be 4 potential nuclei per MPI rank present
    double dTN = 1;
    double dTsigma = 0.0001;
    double RNGSeed = 0.0;

    // Initialize CellType to liquid, LayerID to 1 on device
    ViewI CellType("CellType_Device", LocalDomainSize);
    Kokkos::deep_copy(CellType, Liquid);
    ViewI LayerID("LayerID_Device", LocalDomainSize);
    Kokkos::deep_copy(LayerID, 1);

    // Initialize CritTimeStep values and UndercoolingChange values to values depending on cell coordinates relative to
    // global grid, then copy to the device
    ViewI_H CritTimeStep_Host(Kokkos::ViewAllocateWithoutInitializing("CritTimeStep_Host"), LocalDomainSize);
    ViewF_H UndercoolingChange_Host(Kokkos::ViewAllocateWithoutInitializing("UndercoolingChange_Host"),
                                    LocalDomainSize);
    for (int i = 0; i < LocalDomainSize; i++) {
        int RankZ = i / (MyXSlices * MyYSlices);
        int Rem = i % (MyXSlices * MyYSlices);
        int RankY = Rem % MyYSlices;
        CritTimeStep_Host(i) = RankZ + RankY + MyYOffset + 1;
        UndercoolingChange_Host(i) =
            1.2; // ensures that a cell's nucleation time will be 1 time step after its CritTimeStep value
    }
    using memory_space = Kokkos::DefaultExecutionSpace::memory_space;
    ViewI CritTimeStep = Kokkos::create_mirror_view_and_copy(memory_space(), CritTimeStep_Host);
    ViewF UndercoolingChange = Kokkos::create_mirror_view_and_copy(memory_space(), UndercoolingChange_Host);

    // Without knowing PossibleNuclei_ThisRankThisLayer yet, allocate nucleation data structures based on an estimated
    // size NucleationTimes only exists on the host - used to decide whether or not to call nucleation subroutine
    int EstimatedNuclei_ThisRankThisLayer = NMax * pow(deltax, 3) * LocalActiveDomainSize;
    ViewI_H NucleationTimes_Host(Kokkos::ViewAllocateWithoutInitializing("NucleationTimes_Host"),
                                 EstimatedNuclei_ThisRankThisLayer);
    ViewI NucleiLocation(Kokkos::ViewAllocateWithoutInitializing("NucleiLocation"), EstimatedNuclei_ThisRankThisLayer);
    ViewI NucleiGrainID(Kokkos::ViewAllocateWithoutInitializing("NucleiGrainID"), EstimatedNuclei_ThisRankThisLayer);

    // Counters for nucleation events
    int NucleationCounter;                // total nucleation events checked
    int PossibleNuclei_ThisRankThisLayer; // total successful nucleation events
    int Nuclei_WholeDomain = 100;         // NucleiGrainID should start at -101

    NucleiInit(layernumber, RNGSeed, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nx, ny, nzActive, ZBound_Low, id, NMax,
               dTN, dTsigma, deltax, NucleiLocation, NucleationTimes_Host, NucleiGrainID, CellType, CritTimeStep,
               UndercoolingChange, LayerID, PossibleNuclei_ThisRankThisLayer, Nuclei_WholeDomain, AtNorthBoundary,
               AtSouthBoundary, AtEastBoundary, AtWestBoundary, NucleationCounter);

    // Copy results back to host to check
    ViewI_H NucleiLocation_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), NucleiLocation);
    ViewI_H NucleiGrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), NucleiGrainID);

    // Was the nucleation counter initialized to zero?
    EXPECT_EQ(NucleationCounter, 0);

    // Is the total number of nuclei in the system correct?
    EXPECT_EQ(Nuclei_WholeDomain, 100 + 4 * np);
    for (int n = 0; n < PossibleNuclei_ThisRankThisLayer; n++) {
        // Are the nuclei grain IDs negative numbers in the expected range based on the inputs?
        EXPECT_GT(NucleiGrainID_Host(n), -(100 + 4 * np + 1));
        EXPECT_LT(NucleiGrainID_Host(n), -100);
        // Are the correct undercooling values associated with the correct cell locations?
        int CellLocation = NucleiLocation_Host(n);
        int RankZ = CellLocation / (MyXSlices * MyYSlices);
        int Rem = CellLocation % (MyXSlices * MyYSlices);
        int RankY = Rem % MyYSlices;
        int Expected_NucleationTime = RankZ + RankY + MyYOffset + 2;
        EXPECT_EQ(NucleationTimes_Host(n), Expected_NucleationTime);
        // Are the nucleation events in order of the time steps at which they may occur?
        if (n < PossibleNuclei_ThisRankThisLayer - 2) {
            EXPECT_LE(NucleationTimes_Host(n), NucleationTimes_Host(n + 1));
        }
    }
}

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, grain_init_tests) {
    testSubstrateInit_ConstrainedGrowth();
    testBaseplateInit_FromGrainSpacing();
    testPowderInit();
    testCellTypeInit();
    testNucleiInit();
}

} // end namespace Test
