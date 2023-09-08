// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAcelldata.hpp"
#include "CAfunctions.hpp"
#include "CAinitialize.hpp"
#include "CAnucleation.hpp"
#include "CAparsefiles.hpp"
#include "CAtypes.hpp"

#include <gtest/gtest.h>

#include "mpi.h"

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
// cell_init_tests
//---------------------------------------------------------------------------//
void testCellDataInit_SingleGrain() {

    using memory_space = TEST_MEMSPACE;

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    int singleGrainOrientation = 0;

    // Let overall domain be 5 cells in X and Z, 50 cells in Y
    // This should in turn place the single grain in the cell at X = Y = 2 and Z = 24 (domain center)
    // Domain for each rank
    int nx = 5;
    int ny = 50;
    int nz = 5;
    int expectedGrainX = floorf(static_cast<float>(nx) / 2.0);
    int expectedGrainY = floorf(static_cast<float>(ny) / 2.0);
    int expectedGrainZ = floorf(static_cast<float>(nz) / 2.0);

    // Decompose domain
    int y_offset, ny_local, NeighborRank_North, NeighborRank_South;
    int DomainSize;
    bool AtNorthBoundary, AtSouthBoundary;
    DomainDecomposition(id, np, ny_local, y_offset, NeighborRank_North, NeighborRank_South, nx, ny, nz, DomainSize,
                        AtNorthBoundary, AtSouthBoundary);

    // Cell data struct
    CellData<memory_space> cellData(DomainSize, DomainSize, nx, ny_local, 0);

    // Init grain
    cellData.init_substrate(id, singleGrainOrientation, nx, ny, nz, ny_local, y_offset, DomainSize);

    // Copy cell type and grain ID back to host to check if the values match - only 1 cell should've been assigned type
    // active and GrainID = 1 (though it may be duplicated in the ghost nodes of other ranks)
    auto GrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cellData.GrainID_AllLayers);
    auto CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cellData.CellType_AllLayers);
    for (int coord_z = 0; coord_z < nz; coord_z++) {
        for (int coord_x = 0; coord_x < nx; coord_x++) {
            for (int coord_y = 0; coord_y < ny_local; coord_y++) {
                int coord_y_global = coord_y + y_offset;
                int coord_1d = get1Dindex(coord_x, coord_y, coord_z, nx, ny_local);
                if ((coord_z == expectedGrainZ) && (coord_x == expectedGrainX) && (coord_y_global == expectedGrainY)) {
                    EXPECT_EQ(GrainID_Host(coord_1d), singleGrainOrientation + 1);
                    EXPECT_EQ(CellType_Host(coord_1d), FutureActive);
                }
                else {
                    EXPECT_EQ(GrainID_Host(coord_1d), 0);
                    EXPECT_EQ(CellType_Host(coord_1d), Liquid);
                }
            }
        }
    }
}

void testCellDataInit_ConstrainedGrowth() {

    using memory_space = TEST_MEMSPACE;

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // Create test data
    int nz = 2;
    int z_layer_bottom = 0;
    int nz_layer = 2;
    int nx = 1;
    double deltax = 1 * pow(10, -6);
    double *ZMaxLayer = new double[1];
    ZMaxLayer[0] = (nz - 1) * deltax;
    // Domain size in Y depends on the number of ranks - each rank has 4 cells in Y
    // Each rank is assigned a different portion of the domain in Y
    int ny = 4 * np;
    int ny_local = 4;
    int y_offset = 4 * id;
    int DomainSize = nx * ny_local * nz_layer;
    int DomainSize_AllLayers = nx * ny_local * nz;

    double FractSurfaceSitesActive = 0.5; // Each rank will have 2 active cells each, on average
    double RNGSeed = 0.0;
    CellData<memory_space> cellData(DomainSize_AllLayers, DomainSize, nx, ny_local, z_layer_bottom);
    cellData.init_substrate(id, FractSurfaceSitesActive, ny_local, nx, ny, y_offset, RNGSeed);
    // Copy CellType, GrainID views to host to check values
    auto CellType_AllLayers_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cellData.CellType_AllLayers);
    auto GrainID_AllLayers_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cellData.GrainID_AllLayers);
    for (int index = 0; index < DomainSize; index++) {
        if (index >= nx * ny_local) {
            // Not at bottom surface - should be liquid cells with GrainID still equal to 0
            EXPECT_EQ(GrainID_AllLayers_Host(index), 0);
            EXPECT_EQ(CellType_AllLayers_Host(index), Liquid);
        }
        else {
            // Check that active cells have GrainIDs > 0, and less than 2 * np + 1 (there are 2 * np different positive
            // GrainIDs used for epitaxial grain seeds)
            if (CellType_AllLayers_Host(index) == FutureActive) {
                EXPECT_GT(GrainID_AllLayers_Host(index), 0);
                EXPECT_LT(GrainID_AllLayers_Host(index), 2 * np + 1);
            }
            else {
                // Liquid cells should still have GrainID = 0
                EXPECT_EQ(GrainID_AllLayers_Host(index), 0);
            }
        }
    }
}

void testCellDataInit(bool PowderFirstLayer) {

    using memory_space = TEST_MEMSPACE;
    using view_int = Kokkos::View<int *, memory_space>;
    using view_float = Kokkos::View<float *, memory_space>;

    int id, np;
    // Get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Get individual process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Create test data - 3 layers, starting on layer 0
    double deltax = 1 * pow(10, -6);
    int nz = 5;
    double ZMinLayer[3];
    double ZMaxLayer[3];

    // First layer: Z = 0 through 2
    int nz_layer = 3;
    int z_layer_bottom = 0;
    ZMinLayer[0] = 0;
    ZMaxLayer[0] = 2 * deltax;
    double ZMin = ZMinLayer[0];
    // Second layer: Z = 1 through 3
    ZMinLayer[1] = 1 * deltax;
    ZMaxLayer[1] = 3 * deltax;
    // Third layer: Z = 2 through 4
    ZMinLayer[2] = 2 * deltax;
    ZMaxLayer[2] = 4 * deltax;

    int nx = 3;
    // Each rank is assigned a different portion of the domain in Y
    int ny = 3 * np;
    int ny_local = 3;
    int y_offset = 3 * id;

    // If there is a powder layer, the baseplate should be Z = 0 through 1 w/ powder for the top row of cells, otherwise
    // it should be 0 through 2
    // If there is no powder layer, the baseplate should be Z = 0 through 2 with no powder
    double BaseplateTopZ;
    int BaseplateSize, ExpectedNumPowderGrainsPerLayer;
    if (PowderFirstLayer) {
        BaseplateTopZ = deltax;
        BaseplateSize = nx * ny_local * (round((ZMaxLayer[0] - ZMin) / deltax));
        ExpectedNumPowderGrainsPerLayer = nx * ny_local * np;
    }
    else {
        BaseplateTopZ = 2 * deltax;
        BaseplateSize = nx * ny_local * (round((ZMaxLayer[0] - ZMin) / deltax) + 1);
        ExpectedNumPowderGrainsPerLayer = 0;
    }
    int DomainSize = nx * ny_local * nz_layer;
    int DomainSize_AllLayers = nx * ny_local * nz;
    // There are 45 * np total cells in this domain (nx * ny * nz)
    // Each rank has 45 cells - the bottom 27 cells are assigned baseplate Grain ID values, unless the powder layer of
    // height 1 is used, in which case only the bottom 18 cells should be assigned baseplate Grain ID values. The top
    // cells (Z > 2) are outside the first layer of the domain and are not assigned Grain IDs with the rest of the
    // baseplate. This grain spacing ensures that there will be 1 grain per number of MPI ranks present (larger when
    // powder layer is present as the baseplate will only have a third as many cells)
    double SubstrateGrainSpacing;
    if (PowderFirstLayer)
        SubstrateGrainSpacing = 2.62;
    else
        SubstrateGrainSpacing = 3.0;
    double RNGSeed = 0.0;

    // Used if Powderfirstlayer = true
    double PowderActiveFraction = 1.0;

    // unused views in constructor
    NList NeighborX, NeighborY, NeighborZ;
    NeighborListInit(NeighborX, NeighborY, NeighborZ);
    view_float GrainUnitVector(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), 0);
    view_float DiagonalLength(Kokkos::ViewAllocateWithoutInitializing("DiagonalLength"), 0);
    view_float DOCenter(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 0);
    view_float CritDiagonalLength(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength"), 0);

    // Create dummy temperature data
    view_int NumberOfSolidificationEvents(Kokkos::ViewAllocateWithoutInitializing("NumberOfSolidificationEvents"),
                                          DomainSize);
    Kokkos::parallel_for(
        "InitTestTemperatureData", DomainSize, KOKKOS_LAMBDA(const int &index) {
            int coord_x = getCoordX(index, nx, ny_local);
            int coord_y = getCoordY(index, nx, ny_local);
            // Assign some of these a value of 0 (these will be solid cells), and others a positive value
            // (these will be tempsolid cells)
            if (coord_x + coord_y % 2 == 0)
                NumberOfSolidificationEvents(index) = 0;
            else
                NumberOfSolidificationEvents(index) = 1;
        });
    Kokkos::fence();

    // Call constructor
    CellData<memory_space> cellData(DomainSize_AllLayers, DomainSize, nx, ny_local, z_layer_bottom);
    cellData.init_substrate("", false, false, nx, ny, nz, DomainSize, ZMaxLayer, ZMin, deltax, ny_local, y_offset,
                            z_layer_bottom, id, RNGSeed, SubstrateGrainSpacing, PowderActiveFraction,
                            NumberOfSolidificationEvents, BaseplateTopZ);

    // Copy GrainID results back to host to check first layer's initialization
    auto GrainID_AllLayers_H = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cellData.GrainID_AllLayers);

    // Baseplate grains - cells should have GrainIDs between 1 and np (inclusive)
    for (int i = 0; i < BaseplateSize; i++) {
        EXPECT_GT(GrainID_AllLayers_H(i), 0);
        EXPECT_LT(GrainID_AllLayers_H(i), np + 1);
    }

    if (PowderFirstLayer) {
        // Powder grains should have unique Grain ID values larger than 0 and smaller than
        // NextLayer_FirstEpitaxialGrainID
        EXPECT_EQ(cellData.NextLayer_FirstEpitaxialGrainID, np + 1 + ExpectedNumPowderGrainsPerLayer);
        // Powder should only exist at cells corresponding to Z = 2
        int BottomPowderLayer = nx * ny_local * 2;
        int TopPowderLayer = nx * ny_local * 3;
        for (int i = BottomPowderLayer; i < TopPowderLayer; i++) {
            EXPECT_GT(GrainID_AllLayers_H(i), 0);
            EXPECT_LT(GrainID_AllLayers_H(i), cellData.NextLayer_FirstEpitaxialGrainID);
        }
    }
    else {
        // Next unused GrainID should be the number of grains present in the baseplate plus 1 (since GrainID = 0 is not
        // used for any baseplate grains)
        EXPECT_EQ(cellData.NextLayer_FirstEpitaxialGrainID, np + 1);
    }

    // Copy cell types back to host to check
    auto CellType_AllLayers_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cellData.CellType_AllLayers);
    for (int index = 0; index < DomainSize; index++) {
        int coord_x = getCoordX(index, nx, ny_local);
        int coord_y = getCoordY(index, nx, ny_local);
        // Cells with no associated solidification events should be solid, others TempSolid
        if (coord_x + coord_y % 2 == 0)
            EXPECT_EQ(CellType_AllLayers_Host(index), Solid);
        else
            EXPECT_EQ(CellType_AllLayers_Host(index), TempSolid);
    }
    int PreviousLayer_FirstEpitaxialGrainID = cellData.NextLayer_FirstEpitaxialGrainID;
    z_layer_bottom = 1;
    // Initialize the next layer using the same time-temperature history - powder should span cells at Z = 3
    ExpectedNumPowderGrainsPerLayer = nx * ny_local * np;
    cellData.init_next_layer(1, id, nx, ny, ny_local, y_offset, z_layer_bottom, DomainSize, false, ZMin, ZMaxLayer,
                             deltax, RNGSeed, PowderActiveFraction, NumberOfSolidificationEvents);

    // Copy all grain IDs for all layers back to the host to check that they match
    // and that the powder layer was initialized correctly
    GrainID_AllLayers_H = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cellData.GrainID_AllLayers);

    // Check that the right number of powder grains were initialized
    EXPECT_EQ(cellData.NextLayer_FirstEpitaxialGrainID,
              PreviousLayer_FirstEpitaxialGrainID + ExpectedNumPowderGrainsPerLayer);

    // Powder grains should have unique Grain ID values between PreviousLayer_FirstEpitaxialGrainID and
    // NextLayer_FirstEpitaxialGrainID - 1
    int BottomPowderLayer = nx * ny_local * 3;
    int TopPowderLayer = nx * ny_local * 4;
    for (int index_AllLayers = BottomPowderLayer; index_AllLayers < TopPowderLayer; index_AllLayers++) {
        EXPECT_GT(GrainID_AllLayers_H(index_AllLayers), PreviousLayer_FirstEpitaxialGrainID - 1);
        EXPECT_LT(GrainID_AllLayers_H(index_AllLayers), cellData.NextLayer_FirstEpitaxialGrainID);
    }

    // Subview grain IDs should match the grain IDs overall
    auto GrainID = cellData.getGrainIDSubview();
    auto GrainID_H = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);
    for (int index = 0; index < DomainSize; index++) {
        int index_AllLayers = index + z_layer_bottom * nx * ny_local;
        EXPECT_EQ(GrainID_H(index), GrainID_AllLayers_H(index_AllLayers));
    }

    // Copy cell types back to host to check - should be the same as the previous layer as the same time-temperature
    // history was used
    CellType_AllLayers_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cellData.CellType_AllLayers);
    for (int index = 0; index < DomainSize; index++) {
        int coord_x = getCoordX(index, nx, ny_local);
        int coord_y = getCoordY(index, nx, ny_local);
        // Cells with no associated solidification events should be solid, others TempSolid
        if (coord_x + coord_y % 2 == 0)
            EXPECT_EQ(CellType_AllLayers_Host(index), Solid);
        else
            EXPECT_EQ(CellType_AllLayers_Host(index), TempSolid);
    }
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, cell_init_tests) {
    testCellDataInit_SingleGrain();
    testCellDataInit_ConstrainedGrowth();
    // For non-constrained solidification problems, test w/ and w/o space left for powder layer
    testCellDataInit(true);
    testCellDataInit(false);
}
} // end namespace Test
