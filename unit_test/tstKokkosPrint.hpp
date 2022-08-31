
#include <Kokkos_Core.hpp>

#include "CAprint.hpp"
#include "CAtypes.hpp"

#include <gtest/gtest.h>

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
void testPrintExaConstitDefaultRVE() {

    // File name/path for test RVE output (each rank writes/reads different file)
    std::string BaseFileName = "TestRVERank_0";
    std::string PathToOutput = "";

    // Create test grid - set up so that the RVE is 5 cells in X, Y, and Z
    int nx = 10;
    int ny = 10;
    int nz = 10;
    int NumberOfLayers = 10;
    double deltax = 0.0001; // in meters
    int RVESize = 0.0005 / deltax;

    // Create test data
    ViewI3D_H GrainID_WholeDomain(Kokkos::ViewAllocateWithoutInitializing("GrainID_WholeDomain"), nz, nx, ny);
    ViewI3D_H LayerID_WholeDomain(Kokkos::ViewAllocateWithoutInitializing("LayerID_WholeDomain"), nz, nx, ny);
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                LayerID_WholeDomain(k, i, j) = k;
                GrainID_WholeDomain(k, i, j) = i + j;
            }
        }
    }

    // Print RVE
    PrintExaConstitDefaultRVE(BaseFileName, PathToOutput, nx, ny, nz, LayerID_WholeDomain, GrainID_WholeDomain, deltax,
                              NumberOfLayers, RVESize);

    // Check printed RVE
    std::ifstream GrainplotE;
    std::string ExpectedFilename = BaseFileName + "_ExaConstit.csv";
    GrainplotE.open(ExpectedFilename);
    std::string line;
    std::getline(GrainplotE, line);
    EXPECT_TRUE(line == "Coordinates are in CA units, 1 cell = 0.0001 m. Data is cell-centered. Origin at 3,3,4 , "
                        "domain size is 5 by 5 by 5 cells");
    std::getline(GrainplotE, line);
    EXPECT_TRUE(line == "X coord, Y coord, Z coord, Grain ID");
    for (int k = 4; k < 9; k++) {
        for (int i = 3; i < 8; i++) {
            for (int j = 3; j < 8; j++) {
                std::string ExpectedLine =
                    std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + "," + std::to_string(i + j);
                std::getline(GrainplotE, line);
                EXPECT_TRUE(line == ExpectedLine);
            }
        }
    }
    GrainplotE.close();
}

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, print_test) { testPrintExaConstitDefaultRVE(); }

} // end namespace Test
