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
    int DecompositionStrategy = 1;
    int ny = 4 * np;
    int MyYSlices = 4;
    int MyYOffset = 4 * id;
    int LocalActiveDomainSize = MyXSlices * MyYSlices * nzActive;
    int LocalDomainSize = MyXSlices * MyYSlices * nz;
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

    double FractSurfaceSitesActive = 0.5; // Each rank will have 2 active cells each, on average
    double RNGSeed = 0.0;
    std::string GrainOrientationFile = checkFileInstalled("GrainOrientationVectors.csv", id);
    int NGrainOrientations = 10000; // Number of grain orientations considered in the simulation
    ViewF GrainUnitVector(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), 9 * NGrainOrientations);
    OrientationInit(id, NGrainOrientations, GrainUnitVector, GrainOrientationFile);

    // Initialize neighbor lists
    NList NeighborX, NeighborY, NeighborZ;
    NeighborListInit(NeighborX, NeighborY, NeighborZ);

    // Initialize views - set initial GrainID values to 0, all CellType values to liquid
    ViewI CellType(Kokkos::ViewAllocateWithoutInitializing("CellType"), LocalDomainSize);
    Kokkos::deep_copy(CellType, Liquid);
    ViewI GrainID("GrainID", LocalDomainSize);
    ViewF DiagonalLength(Kokkos::ViewAllocateWithoutInitializing("DiagonalLength"), LocalActiveDomainSize);
    ViewF DOCenter(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * LocalActiveDomainSize);
    ViewF CritDiagonalLength(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength"), 26 * LocalActiveDomainSize);

    // Buffer sizes
    int BufSizeX = MyXSlices;
    int BufSizeY = 0;
    int BufSizeZ = nzActive;

    // Send/recv buffers for ghost node data should be initialized with zeros
    Buffer2D BufferSouthSend("BufferSouthSend", BufSizeX * BufSizeZ, 5);
    Buffer2D BufferNorthSend("BufferNorthSend", BufSizeX * BufSizeZ, 5);
    Buffer2D BufferEastSend("BufferEastSend", BufSizeY * BufSizeZ, 5);
    Buffer2D BufferWestSend("BufferWestSend", BufSizeY * BufSizeZ, 5);
    Buffer2D BufferNorthEastSend("BufferNorthEastSend", BufSizeZ, 5);
    Buffer2D BufferNorthWestSend("BufferNorthWestSend", BufSizeZ, 5);
    Buffer2D BufferSouthEastSend("BufferSouthEastSend", BufSizeZ, 5);
    Buffer2D BufferSouthWestSend("BufferSouthWestSend", BufSizeZ, 5);
    Buffer2D BufferSouthRecv("BufferSouthRecv", BufSizeX * BufSizeZ, 5);
    Buffer2D BufferNorthRecv("BufferNorthRecv", BufSizeX * BufSizeZ, 5);
    Buffer2D BufferEastRecv("BufferEastRecv", BufSizeY * BufSizeZ, 5);
    Buffer2D BufferWestRecv("BufferWestRecv", BufSizeY * BufSizeZ, 5);
    Buffer2D BufferNorthEastRecv("BufferNorthEastRecv", BufSizeZ, 5);
    Buffer2D BufferNorthWestRecv("BufferNorthWestRecv", BufSizeZ, 5);
    Buffer2D BufferSouthEastRecv("BufferSouthEastRecv", BufSizeZ, 5);
    Buffer2D BufferSouthWestRecv("BufferSouthWestRecv", BufSizeZ, 5);
    SubstrateInit_ConstrainedGrowth(
        id, FractSurfaceSitesActive, MyXSlices, MyYSlices, nx, ny, MyXOffset, MyYOffset, NeighborX, NeighborY,
        NeighborZ, GrainUnitVector, NGrainOrientations, CellType, GrainID, DiagonalLength, DOCenter, CritDiagonalLength,
        RNGSeed, np, DecompositionStrategy, BufferWestSend, BufferEastSend, BufferNorthSend, BufferSouthSend,
        BufferNorthEastSend, BufferNorthWestSend, BufferSouthEastSend, BufferSouthWestSend, BufSizeX, BufSizeY,
        AtNorthBoundary, AtSouthBoundary, AtEastBoundary, AtWestBoundary);

    // Copy CellType, GrainID views to host to check values
    ViewI_H CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
    ViewI_H GrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);
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
                                   MyYOffset, id, deltax, GrainID, RNGSeed, NextLayer_FirstEpitaxialGrainID, nz, false);

    // Copy results back to host to check
    ViewI_H GrainID_H = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);

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
               id, GrainID, RNGSeed, NextLayer_FirstEpitaxialGrainID, 1.0);

    // Copy results back to host to check
    ViewI_H GrainID_H = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);

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
void testCellTypeInit_NoRemelt() {

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
    using memory_space = Kokkos::DefaultExecutionSpace::memory_space;
    NList NeighborX, NeighborY, NeighborZ;
    NeighborListInit(NeighborX, NeighborY, NeighborZ);

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
    CellTypeInit_NoRemelt(layernumber, id, np, DecompositionStrategy, MyXSlices, MyYSlices, MyXOffset, MyYOffset,
                          ZBound_Low, nz, LocalActiveDomainSize, LocalDomainSize, CellType, CritTimeStep, NeighborX,
                          NeighborY, NeighborZ, NGrainOrientations, GrainUnitVector, DiagonalLength, GrainID,
                          CritDiagonalLength, DOCenter, LayerID, BufferWestSend, BufferEastSend, BufferNorthSend,
                          BufferSouthSend, BufferNorthEastSend, BufferNorthWestSend, BufferSouthEastSend,
                          BufferSouthWestSend, BufSizeX, BufSizeY, AtNorthBoundary, AtSouthBoundary, AtEastBoundary,
                          AtWestBoundary);

    // Copy views back to host to check the results
    ViewF_H DiagonalLength_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), DiagonalLength);
    ViewF_H CritDiagonalLength_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CritDiagonalLength);
    ViewF_H DOCenter_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), DOCenter);
    ViewI_H CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);

    // Solid cells where no temperature data existed
    // Active cells separate solid-liquid cells, as well as at the active region bottom (implicitly borders a solid
    // cell, since the layers below should be solid) Liquid cells are all other cells
    // Based on the orientation of the octahedron and the distance to each neighboring cell on the
    // grid, these critical diagonal lengths were calculated to be the expected values required for
    // the octahedron to engulf the centers of each neighboring cells
    std::vector<float> CritDiagonalLength_Expected{
        1.7009791, 2.1710088, 1.9409314, 1.9225541, 2.2444904, 2.6762404, 2.7775438, 2.6560757, 2.1579263,
        1.5170407, 1.5170407, 2.0631709, 2.4170833, 2.4170833, 2.0631709, 1.4566354, 1.4566354, 1.7009791,
        1.9409314, 2.1710088, 2.2444904, 1.9225541, 2.6560757, 2.1579263, 2.6762404, 2.7775438};
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
                        // Check critical diagonal length values against expected values
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

    auto BufferSouthSend_H = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), BufferSouthSend);
    auto BufferNorthSend_H = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), BufferNorthSend);
    // Further check that active cell data was properly loaded into send buffers, and that locations in the send buffers
    // not corresponding to active cells were left alone (should still be 0s)
    for (int k = ZBound_Low; k <= ZBound_High; k++) {
        for (int i = 0; i < MyXSlices; i++) {
            int GNPosition = (k - ZBound_Low) * BufSizeX + i; // Position of cell in buffer
            // Check the south buffer - Data being sent to the "south" (BufferSouthSend) is from active cells at Y = 1
            int D3D1ConvPositionGlobal_South =
                k * MyXSlices * MyYSlices + i * MyYSlices + 1; // Position of cell on grid
            if ((CellType_Host(D3D1ConvPositionGlobal_South) == Active) && (!(AtSouthBoundary))) {
                EXPECT_FLOAT_EQ(BufferSouthSend_H(GNPosition, 0), 1);
                EXPECT_FLOAT_EQ(BufferSouthSend_H(GNPosition, 1), i + MyXOffset + 0.5);
                EXPECT_FLOAT_EQ(BufferSouthSend_H(GNPosition, 2), 1 + MyYOffset + 0.5);
                EXPECT_FLOAT_EQ(BufferSouthSend_H(GNPosition, 3), k + 0.5);
                EXPECT_FLOAT_EQ(BufferSouthSend_H(GNPosition, 4), 0.01);
            }
            else {
                for (int l = 0; l < 5; l++) {
                    EXPECT_FLOAT_EQ(BufferSouthSend_H(GNPosition, l), 0.0);
                }
            }
            // Check the north buffer - Data being sent to the "north" (BufferNorthSend) is from active cells at Y = 2
            int D3D1ConvPositionGlobal_North = k * MyXSlices * MyYSlices + i * MyYSlices + 2;
            if ((CellType_Host(D3D1ConvPositionGlobal_North) == Active) && (!(AtNorthBoundary))) {
                EXPECT_FLOAT_EQ(BufferNorthSend_H(GNPosition, 0), 1);
                EXPECT_FLOAT_EQ(BufferNorthSend_H(GNPosition, 1), i + MyXOffset + 0.5);
                EXPECT_FLOAT_EQ(BufferNorthSend_H(GNPosition, 2), 2 + MyYOffset + 0.5);
                EXPECT_FLOAT_EQ(BufferNorthSend_H(GNPosition, 3), k + 0.5);
                EXPECT_FLOAT_EQ(BufferNorthSend_H(GNPosition, 4), 0.01);
            }
            else {
                for (int l = 0; l < 5; l++) {
                    EXPECT_FLOAT_EQ(BufferNorthSend_H(GNPosition, l), 0.0);
                }
            }
        }
    }
}

void testCellTypeInit_Remelt() {

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // Domain for each rank
    int MyXSlices = id + 5;
    int MyYSlices = id + 5;
    int nzActive = 5;
    int ZBound_Low = 2;
    int nz = nzActive + ZBound_Low;
    int LocalActiveDomainSize = MyXSlices * MyYSlices * nzActive;
    int LocalDomainSize = MyXSlices * MyYSlices * nz;

    // Temporary host view for initializing CritTimeStep
    ViewI_H CritTimeStep_Host(Kokkos::ViewAllocateWithoutInitializing("CritTimeStep_Host"), LocalDomainSize);
    for (int GlobalZ = 0; GlobalZ < nz; GlobalZ++) {
        for (int RankX = 0; RankX < MyXSlices; RankX++) {
            for (int RankY = 0; RankY < MyYSlices; RankY++) {
                int GlobalD3D1ConvPosition = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                if (GlobalZ < ZBound_Low) {
                    // Not in active domain, assign these a negative value
                    CritTimeStep_Host(GlobalD3D1ConvPosition) = -1;
                }
                else {
                    // Assign some of these a value of 0, and others a positive value
                    if (RankX + RankY % 2 == 0)
                        CritTimeStep_Host(GlobalD3D1ConvPosition) = 0;
                    else
                        CritTimeStep_Host(GlobalD3D1ConvPosition) = 1;
                }
            }
        }
    }

    ViewI CritTimeStep = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), CritTimeStep_Host);
    // Start with cell type of 0
    ViewI CellType("CellType", LocalDomainSize);

    // Initialize cell type values
    CellTypeInit_Remelt(MyXSlices, MyYSlices, LocalActiveDomainSize, CellType, CritTimeStep, id, ZBound_Low);

    // Copy cell types back to host to check
    ViewI_H CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);

    for (int GlobalZ = 0; GlobalZ < nz; GlobalZ++) {
        for (int RankX = 0; RankX < MyXSlices; RankX++) {
            for (int RankY = 0; RankY < MyYSlices; RankY++) {
                int GlobalD3D1ConvPosition = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                if (GlobalZ < ZBound_Low) {
                    // These cells should still be 0 (Wall) - untouched by CellTypeInit_Remelt
                    EXPECT_EQ(CellType_Host(GlobalD3D1ConvPosition), Wall);
                }
                else {
                    // Cells with CritTimeStep of 0 should be solid, others TempSolid
                    if (RankX + RankY % 2 == 0)
                        EXPECT_EQ(CellType_Host(GlobalD3D1ConvPosition), Solid);
                    else
                        EXPECT_EQ(CellType_Host(GlobalD3D1ConvPosition), TempSolid);
                }
            }
        }
    }
}

void testNucleiInit(bool RemeltingYN) {

    // Test is different depending on whether or not remelting is considered
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
                         // without remelting (each cell solidifies once)
    int MaxPotentialNuclei_PerPass = 4 * np;
    // Without remelting, each cell solidifies once, with remelting, a cell can solidify 1-3 times
    int MaxSolidificationEvents_Count;
    if (RemeltingYN)
        MaxSolidificationEvents_Count = 3;
    else
        MaxSolidificationEvents_Count = 1;
    double dTN = 1;
    double dTsigma = 0.0001;
    double RNGSeed = 0.0;

    // Initialize CellType to liquid, LayerID to 1 on device
    ViewI CellType("CellType_Device", LocalDomainSize);
    Kokkos::deep_copy(CellType, Liquid);
    ViewI LayerID("LayerID_Device", LocalDomainSize);
    Kokkos::deep_copy(LayerID, 1);

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

    // This part of the test is different depending on whether remelting is considered
    // Without remelting, initialize MaxSolidificationEvents to 1 for each layer, initialize CritTimeStep values and
    // UndercoolingChange values to values depending on cell coordinates relative to global grid, then copy to the
    // device. LayerTimeTempHistory and NumberOfSolidificationEvents are unused on the device With remelting, initialize
    // MaxSolidificationEvents to 3 for each layer. LayerTimeTempHistory and NumberOfSolidificationEvents are
    // initialized for each cell on the host and copied to the device, while CritTimeStep and UndercoolingChange go
    // unused
    ViewI_H CritTimeStep_Host(Kokkos::ViewAllocateWithoutInitializing("CritTimeStep_Host"), LocalDomainSize);
    ViewF_H UndercoolingChange_Host(Kokkos::ViewAllocateWithoutInitializing("UndercoolingChange_Host"),
                                    LocalDomainSize);
    ViewI_H MaxSolidificationEvents_Host(Kokkos::ViewAllocateWithoutInitializing("MaxSolidificationEvents_Host"), 2);
    ViewI_H NumberOfSolidificationEvents_Host(
        Kokkos::ViewAllocateWithoutInitializing("NumberOfSolidificationEvents_Host"), LocalActiveDomainSize);
    ViewF3D LayerTimeTempHistory_Host(Kokkos::ViewAllocateWithoutInitializing("LayerTimeTempHistory_Host"),
                                      LocalActiveDomainSize, MaxSolidificationEvents_Count, 3);
    MaxSolidificationEvents_Host(0) = MaxSolidificationEvents_Count;
    MaxSolidificationEvents_Host(1) = MaxSolidificationEvents_Count;
    if (RemeltingYN) {
        // Cells solidify 1, 2, or 3 times, depending on their X coordinate
        for (int k = 0; k < nzActive; k++) {
            for (int i = 0; i < MyXSlices; i++) {
                for (int j = 0; j < MyYSlices; j++) {
                    int D3D1ConvPosition = k * MyXSlices * MyYSlices + i * MyYSlices + j;
                    if (i < MyXSlices / 2 - 1)
                        NumberOfSolidificationEvents_Host(D3D1ConvPosition) = 3;
                    else if (i < MyXSlices / 2)
                        NumberOfSolidificationEvents_Host(D3D1ConvPosition) = 2;
                    else
                        NumberOfSolidificationEvents_Host(D3D1ConvPosition) = 1;
                }
            }
        }
        for (int n = 0; n < MaxSolidificationEvents_Count; n++) {
            for (int RankZ = 0; RankZ < nzActive; RankZ++) {
                for (int RankX = 0; RankX < MyXSlices; RankX++) {
                    for (int RankY = 0; RankY < MyYSlices; RankY++) {
                        int D3D1ConvPosition = RankZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                        int GlobalZ = RankZ + ZBound_Low;
                        if (n < NumberOfSolidificationEvents_Host(D3D1ConvPosition)) {
                            LayerTimeTempHistory_Host(D3D1ConvPosition, n, 0) =
                                GlobalZ + RankY + MyYOffset +
                                (LocalActiveDomainSize * n); // melting time step depends on solidification event number
                            LayerTimeTempHistory_Host(D3D1ConvPosition, n, 1) =
                                GlobalZ + RankY + MyYOffset + 1 +
                                (LocalActiveDomainSize *
                                 n); // liquidus time stemp depends on solidification event number
                            LayerTimeTempHistory_Host(D3D1ConvPosition, n, 2) =
                                1.2; // ensures that a cell's nucleation time will be 1 time step after its CritTimeStep
                                     // value
                        }
                    }
                }
            }
        }
    }
    else {
        for (int i = 0; i < LocalDomainSize; i++) {
            int GlobalZ = i / (MyXSlices * MyYSlices);
            int Rem = i % (MyXSlices * MyYSlices);
            int RankY = Rem % MyYSlices;
            CritTimeStep_Host(i) = GlobalZ + RankY + MyYOffset + 1;
            UndercoolingChange_Host(i) =
                1.2; // ensures that a cell's nucleation time will be 1 time step after its CritTimeStep value
        }
    }
    using memory_space = Kokkos::DefaultExecutionSpace::memory_space;
    ViewI CritTimeStep = Kokkos::create_mirror_view_and_copy(memory_space(), CritTimeStep_Host);
    ViewF UndercoolingChange = Kokkos::create_mirror_view_and_copy(memory_space(), UndercoolingChange_Host);
    ViewI MaxSolidificationEvents = Kokkos::create_mirror_view_and_copy(memory_space(), MaxSolidificationEvents_Host);
    ViewI NumberOfSolidificationEvents =
        Kokkos::create_mirror_view_and_copy(memory_space(), NumberOfSolidificationEvents_Host);
    ViewF3D LayerTimeTempHistory = Kokkos::create_mirror_view_and_copy(memory_space(), LayerTimeTempHistory_Host);

    NucleiInit(layernumber, RNGSeed, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nx, ny, nzActive, ZBound_Low, id, NMax,
               dTN, dTsigma, deltax, NucleiLocation, NucleationTimes_Host, NucleiGrainID, CellType, CritTimeStep,
               UndercoolingChange, LayerID, PossibleNuclei_ThisRankThisLayer, Nuclei_WholeDomain, AtNorthBoundary,
               AtSouthBoundary, AtEastBoundary, AtWestBoundary, RemeltingYN, NucleationCounter, MaxSolidificationEvents,
               NumberOfSolidificationEvents, LayerTimeTempHistory);

    // Copy results back to host to check
    ViewI_H NucleiLocation_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), NucleiLocation);
    ViewI_H NucleiGrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), NucleiGrainID);
    MaxSolidificationEvents_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), MaxSolidificationEvents);
    NumberOfSolidificationEvents_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), NumberOfSolidificationEvents);

    // Was the nucleation counter initialized to zero?
    EXPECT_EQ(NucleationCounter, 0);

    // Is the total number of nuclei in the system correct, based on the number of remelting events? Equal probability
    // of creating a nucleus each time a cell resolidifies
    int ExpectedNucleiPerRank = 100 + MaxSolidificationEvents_Count * MaxPotentialNuclei_PerPass;
    EXPECT_EQ(Nuclei_WholeDomain, ExpectedNucleiPerRank);
    for (int n = 0; n < PossibleNuclei_ThisRankThisLayer; n++) {
        // Are the nuclei grain IDs negative numbers in the expected range based on the inputs?
        EXPECT_GT(NucleiGrainID_Host(n), -(100 + ExpectedNucleiPerRank * np + 1));
        EXPECT_LT(NucleiGrainID_Host(n), -100);
        // Are the correct undercooling values associated with the correct cell locations?
        // Cell location is a global position (relative to the bottom of the whole domain, not the layer)
        int GlobalCellLocation = NucleiLocation_Host(n);
        int GlobalZ = GlobalCellLocation / (MyXSlices * MyYSlices);
        int Rem = GlobalCellLocation % (MyXSlices * MyYSlices);
        int RankY = Rem % MyYSlices;
        // Expected nucleation time is known exactly if no remelting
        int Expected_NucleationTimeNoRM = GlobalZ + RankY + MyYOffset + 2;
        if (RemeltingYN) {
            // Expected nucleation time with remelting can be one of 3 possibilities, depending on the associated
            // solidification event
            int AssociatedSEvent = NucleationTimes_Host(n) / LocalActiveDomainSize;
            int Expected_NucleationTimeRM = Expected_NucleationTimeNoRM + AssociatedSEvent * LocalActiveDomainSize;
            EXPECT_EQ(NucleationTimes_Host(n), Expected_NucleationTimeRM);
        }
        else {
            EXPECT_EQ(NucleationTimes_Host(n), Expected_NucleationTimeNoRM);
        }
        // Are the nucleation events in order of the time steps at which they may occur?
        if (n < PossibleNuclei_ThisRankThisLayer - 2) {
            EXPECT_LE(NucleationTimes_Host(n), NucleationTimes_Host(n + 1));
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
TEST(TEST_CATEGORY, grain_init_tests) {
    testSubstrateInit_ConstrainedGrowth();
    testBaseplateInit_FromGrainSpacing();
    testPowderInit();
}
TEST(TEST_CATEGORY, cell_init_test) {
    testCellTypeInit_NoRemelt();
    testCellTypeInit_Remelt();
}
TEST(TEST_CATEGORY, nuclei_init_test) {
    // w/ and w/o remelting
    testNucleiInit(true);
    testNucleiInit(false);
}
} // end namespace Test
