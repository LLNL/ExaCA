// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAcelldata.hpp"
#include "CAfunctions.hpp"
#include "CAinitialize.hpp"
#include "CAinputs.hpp"
#include "CAnucleation.hpp"
#include "CAtypes.hpp"
#include "CAupdate.hpp"

#include <gtest/gtest.h>

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
void testFillSteeringVector_Remelt() {

    // Create views - each rank has 125 cells, 75 of which are part of the active region of the domain
    int nx = 5;
    int ny_local = 5;
    int y_offset = 5;
    int nz = 5;
    int z_layer_bottom = 2;
    int nz_layer = 3;
    int DomainSize = nx * ny_local * nz_layer;
    int DomainSize_AllLayers = nx * ny_local * nz;

    // Initialize neighbor lists
    NList NeighborX, NeighborY, NeighborZ;
    NeighborListInit(NeighborX, NeighborY, NeighborZ);

    // default inputs struct
    Inputs inputs;
    // Initialize cell/temperature structures
    CellData<device_memory_space> cellData(DomainSize_AllLayers, DomainSize, nx, ny_local, z_layer_bottom,
                                           inputs.substrate);
    auto CellType = cellData.getCellTypeSubview();
    auto GrainID = cellData.getGrainIDSubview();
    Temperature<device_memory_space> temperature(DomainSize, 1, inputs.temperature);

    // Fill temperature structure
    Kokkos::parallel_for(
        "TempDataInit", DomainSize, KOKKOS_LAMBDA(const int &index) {
            GrainID(index) = 1;
            CellType(index) = TempSolid;
            // Cell coordinates on this rank in X, Y, and Z (GlobalZ = relative to domain bottom)
            int coord_z = getCoordZ(index, nx, ny_local);
            int coord_y = getCoordY(index, nx, ny_local);
            // Cells at Z = 0 through Z = 2 are Solid, Z = 3 and 4 are TempSolid
            if (coord_z <= 2) {
                // Solid cells should have -1 assigned as their melt/crit time steps
                temperature.LayerTimeTempHistory(index, 0, 0) = -1;
                temperature.LayerTimeTempHistory(index, 0, 1) = -1;
                temperature.LayerTimeTempHistory(index, 0, 2) = 0;
            }
            else {
                // Cells "melt" at a time step corresponding to their Y location in the overall domain (depends on
                // MyYOffset of the rank)
                temperature.LayerTimeTempHistory(index, 0, 0) = coord_y + y_offset + 1;
                // Cells reach liquidus during cooling 2 time steps after melting
                temperature.LayerTimeTempHistory(index, 0, 1) = temperature.LayerTimeTempHistory(index, 0, 0) + 2;
            }
            temperature.LayerTimeTempHistory(index, 0, 2) = 0.2;
            temperature.NumberOfSolidificationEvents(index) = 1;
        });

    // Steering Vector
    ViewI SteeringVector(Kokkos::ViewAllocateWithoutInitializing("SteeringVector"), DomainSize);
    ViewI_H numSteer_Host(Kokkos::ViewAllocateWithoutInitializing("SteeringVectorSize"), 1);
    numSteer_Host(0) = 0;
    // Copy views to device for test
    ViewI numSteer = Kokkos::create_mirror_view_and_copy(device_memory_space(), numSteer_Host);

    int numcycles = 15;
    for (int cycle = 1; cycle <= numcycles; cycle++) {
        // Update cell types, local undercooling each time step, and fill the steering vector
        std::cout << cycle << std::endl;
        FillSteeringVector_Remelt(cycle, DomainSize, nx, ny_local, NeighborX, NeighborY, NeighborZ, temperature,
                                  cellData, nz_layer, SteeringVector, numSteer, numSteer_Host);
    }

    // Copy CellType, SteeringVector, numSteer, UndercoolingCurrent, Buffers back to host to check steering vector
    // construction results
    CellType = cellData.getCellTypeSubview();
    ViewI_H CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
    ViewI_H SteeringVector_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), SteeringVector);
    numSteer_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), numSteer);
    ViewF_H UndercoolingCurrent_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), temperature.UndercoolingCurrent);
    ViewF3D_H LayerTimeTempHistory_Host =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), temperature.LayerTimeTempHistory);

    // Check the modified CellType and UndercoolingCurrent values on the host:
    // Check that the cells corresponding to outside of the "active" portion of the domain have unchanged values
    // Check that the cells corresponding to the "active" portion of the domain have potentially changed values
    // Z = 3: Rank 0, with all cells having GrainID = 0, should have Liquid cells everywhere, with undercooling
    // depending on the Y position Z = 3, Rank > 0, with GrainID > 0, should have TempSolid cells (if not enough time
    // steps to melt), Liquid cells (if enough time steps to melt but not reach the liquidus again), or FutureActive
    // cells (if below liquidus time step). If the cell is FutureActive type, it should also have a local undercooling
    // based on the Rank ID and the time step that it reached the liquidus relative to numcycles Z = 4, all ranks:
    // should have TempSolid cells (if not enough time steps to melt) or Liquid (if enough time steps to melt). Local
    // undercooling should be based on the Rank ID and the time step that it reached the liquidus relative to numcycles
    int FutureActiveCells = 0;
    for (int index = 0; index < DomainSize; index++) {
        int coord_z = getCoordZ(index, nx, ny_local);
        if (coord_z <= 2)
            EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(index), 0.0);
        else {
            if (numcycles < LayerTimeTempHistory_Host(index, 0, 0)) {
                EXPECT_EQ(CellType_Host(index), TempSolid);
                EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(index), 0.0);
            }
            else if ((numcycles >= LayerTimeTempHistory_Host(index, 0, 0)) &&
                     (numcycles <= LayerTimeTempHistory_Host(index, 0, 1))) {
                EXPECT_EQ(CellType_Host(index), Liquid);
                EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(index), 0.0);
            }
            else {
                EXPECT_FLOAT_EQ(UndercoolingCurrent_Host(index),
                                (numcycles - LayerTimeTempHistory_Host(index, 0, 1)) * 0.2);
                if (coord_z == 4)
                    EXPECT_EQ(CellType_Host(index), Liquid);
                else {
                    EXPECT_EQ(CellType_Host(index), FutureActive);
                    FutureActiveCells++;
                }
            }
        }
    }
    // Check the steering vector values on the host
    EXPECT_EQ(FutureActiveCells, numSteer_Host(0));
    for (int i = 0; i < FutureActiveCells; i++) {
        // This cell should correspond to a cell at GlobalZ = 3 (RankZ = 1), and some X and Y
        int LowerBoundCellLocation = nx * ny_local - 1;
        int UpperBoundCellLocation = 2 * nx * ny_local;
        EXPECT_GT(SteeringVector_Host(i), LowerBoundCellLocation);
        EXPECT_LT(SteeringVector_Host(i), UpperBoundCellLocation);
    }
}

void testcalcCritDiagonalLength() {
    using memory_space = TEST_MEMSPACE;
    using view_type = Kokkos::View<float *, memory_space>;

    // 2 "cells" in the domain
    int DomainSize = 2;
    // First orientation is orientation 12 (starting indexing at 0) of GrainOrientationVectors.csv
    // Second orientation is orientation 28 (starting indexing at 0)
    std::vector<float> GrainUnitVectorV{0.52877,  0.651038, -0.544565, -0.572875, 0.747163,  0.336989,
                                        0.626272, 0.133778, 0.768041,  0.736425,  -0.530983, 0.419208,
                                        0.664512, 0.683971, -0.301012, -0.126894, 0.500241,  0.856538};
    // Load unit vectors into view
    view_type GrainUnitVector(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), 9 * DomainSize);
    auto GrainUnitVector_Host = Kokkos::create_mirror_view(Kokkos::HostSpace(), GrainUnitVector);
    for (int i = 0; i < 9 * DomainSize; i++) {
        GrainUnitVector_Host(i) = GrainUnitVectorV[i];
    }
    // Copy octahedron center and grain unit vector data to the device
    GrainUnitVector = Kokkos::create_mirror_view_and_copy(memory_space(), GrainUnitVector_Host);

    // Initialize neighbor lists
    NList NeighborX, NeighborY, NeighborZ;
    NeighborListInit(NeighborX, NeighborY, NeighborZ);

    // Load octahedron centers into view
    view_type DOCenter(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * DomainSize);
    Kokkos::parallel_for(
        "TestInitCritDiagonalLength", 1, KOKKOS_LAMBDA(const int) {
            // Octahedron center for first cell (does not align with CA cell center, which is at 31.5, 3.5, 1.5)
            DOCenter(0) = 30.611142;
            DOCenter(1) = 3.523741;
            DOCenter(2) = 0.636301;
            // Octahedron center for second cell (aligns with CA cell center at 110.5, 60.5, 0.5))
            DOCenter(3) = 110.5;
            DOCenter(4) = 60.5;
            DOCenter(5) = 0.5;
        });

    // View for storing calculated critical diagonal lengths
    view_type CritDiagonalLength(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength"), 26 * DomainSize);

    // Expected results of calculations
    std::vector<float> CritDiagonalLength_Expected{
        2.732952,  3.403329,  2.062575, 3.708569,  1.7573346, 3.038192, 4.378946,  1.7094444, 3.2407162,
        2.646351,  1.5504522, 3.114523, 3.121724,  3.2783692, 0.177466, 3.1648562, 1.472129,  2.497145,
        3.2025092, 1.914978,  3.057610, 1.9366806, 3.639777,  3.351627, 3.3160222, 1.4418885, 1.715195,
        1.914002,  1.927272,  2.291471, 1.851513,  2.346452,  2.236490, 2.902006,  2.613426,  1.576758,
        1.576758,  2.248777,  2.266173, 2.266173,  2.248777,  1.527831, 1.527831,  1.715195,  1.927272,
        1.914002,  1.851513,  2.291471, 2.902006,  2.613426,  2.346452, 2.236490};

    // For first grain, calculate critical diagonal lenghts
    calcCritDiagonalLength(0, 31.5, 3.5, 1.5, DOCenter(0), DOCenter(1), DOCenter(2), NeighborX, NeighborY, NeighborZ, 0,
                           GrainUnitVector, CritDiagonalLength);
    // For second grain, calculate critical diagonal lenghts
    calcCritDiagonalLength(1, 110.5, 60.5, 0.5, DOCenter(3), DOCenter(4), DOCenter(5), NeighborX, NeighborY, NeighborZ,
                           1, GrainUnitVector, CritDiagonalLength);

    // Copy calculated critical diagonal lengths to the host to check against expected values
    auto CritDiagonalLength_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CritDiagonalLength);

    // Check results
    for (int i = 0; i < 26 * DomainSize; i++) {
        EXPECT_FLOAT_EQ(CritDiagonalLength_Host(i), CritDiagonalLength_Expected[i]);
    }
}

void testcreateNewOctahedron() {
    using memory_space = TEST_MEMSPACE;
    using view_type = Kokkos::View<float *, memory_space>;

    // Create a 18-cell domain with origin at (10, 10, 0)
    int nx = 3;
    int ny_local = 2;
    int nz_layer = 3;
    int DomainSize = nx * ny_local * nz_layer;
    int y_offset = 10;

    // Create diagonal length and octahedron center data structures on device
    view_type DiagonalLength(Kokkos::ViewAllocateWithoutInitializing("DiagonalLength"), DomainSize);
    view_type DOCenter(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * DomainSize);

    // Octahedra now use the layer coordinates, not the coordinates of the multilayer domain
    for (int coord_z = 0; coord_z < nz_layer; coord_z++) {
        for (int coord_x = 0; coord_x < nx; coord_x++) {
            for (int coord_y = 0; coord_y < ny_local; coord_y++) {
                int index = get1Dindex(coord_x, coord_y, coord_z, nx, ny_local);
                createNewOctahedron(index, DiagonalLength, DOCenter, coord_x, coord_y, y_offset, coord_z);
            }
        }
    }

    // Copy back to host and check values
    auto DiagonalLength_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), DiagonalLength);
    auto DOCenter_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), DOCenter);
    for (int coord_z = 0; coord_z < nz_layer; coord_z++) {
        for (int coord_x = 0; coord_x < nx; coord_x++) {
            for (int coord_y = 0; coord_y < ny_local; coord_y++) {
                int index = get1Dindex(coord_x, coord_y, coord_z, nx, ny_local);
                EXPECT_FLOAT_EQ(DiagonalLength_Host(index), 0.01);
                EXPECT_FLOAT_EQ(DOCenter_Host(3 * index), coord_x + 0.5);
                EXPECT_FLOAT_EQ(DOCenter_Host(3 * index + 1), coord_y + y_offset + 0.5);
                EXPECT_FLOAT_EQ(DOCenter_Host(3 * index + 2), coord_z + 0.5);
            }
        }
    }
}

void testConvertGrainIDForBuffer() {
    using memory_space = TEST_MEMSPACE;
    using view_type = Kokkos::View<int *, memory_space>;

    // Create a list of integer grain ID values
    // Test positive and negative values
    int NGrainIDValues = 7;
    int TotalNGrainIDValues = 2 * NGrainIDValues;
    std::vector<int> GrainIDV(TotalNGrainIDValues);
    GrainIDV[0] = 1;
    GrainIDV[1] = 2;
    GrainIDV[2] = 3;
    GrainIDV[3] = 10000;
    GrainIDV[4] = 10001;
    GrainIDV[5] = 19999;
    GrainIDV[6] = 30132129;
    for (int n = NGrainIDValues; n < TotalNGrainIDValues; n++) {
        GrainIDV[n] = -GrainIDV[n - NGrainIDValues];
    }
    view_type GrainID(Kokkos::ViewAllocateWithoutInitializing("GrainID"), TotalNGrainIDValues);
    auto GrainID_Host = Kokkos::create_mirror_view(Kokkos::HostSpace(), GrainID);
    for (int n = 0; n < TotalNGrainIDValues; n++) {
        GrainID_Host(n) = GrainIDV[n];
    }
    GrainID = Kokkos::create_mirror_view_and_copy(memory_space(), GrainID_Host);

    int NGrainOrientations = 10000;

    // Check that these were converted to their components and back correctly
    view_type GrainID_Converted(Kokkos::ViewAllocateWithoutInitializing("GrainID_Converted"), TotalNGrainIDValues);

    Kokkos::parallel_for(
        "TestInitGrainIDs", TotalNGrainIDValues, KOKKOS_LAMBDA(const int &n) {
            int MyGrainOrientation = getGrainOrientation(GrainID(n), NGrainOrientations, false);
            int MyGrainNumber = getGrainNumber(GrainID(n), NGrainOrientations);
            GrainID_Converted[n] = getGrainID(NGrainOrientations, MyGrainOrientation, MyGrainNumber);
        });
    auto GrainID_Converted_Host = Kokkos::create_mirror_view(Kokkos::HostSpace(), GrainID_Converted);
    for (int n = 0; n < TotalNGrainIDValues; n++)
        EXPECT_EQ(GrainIDV[n], GrainID_Converted_Host(n));
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, cell_update_tests) {
    testFillSteeringVector_Remelt();
    testcalcCritDiagonalLength();
    testcreateNewOctahedron();
    testConvertGrainIDForBuffer();
}

} // end namespace Test
