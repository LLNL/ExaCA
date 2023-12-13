// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_CELLDATA_HPP
#define EXACA_CELLDATA_HPP

#include "CAconfig.hpp"
#include "CAfunctions.hpp"
#include "CAgrid.hpp"
#include "CAinputs.hpp"
#include "CAparsefiles.hpp"
#include "CAtypes.hpp"
#include "mpi.h"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

// Data that belongs to individual cells governing the current state and various ID values
template <typename MemorySpace>
struct CellData {

    using memory_space = MemorySpace;
    using view_type_int = Kokkos::View<int *, memory_space>;
    using view_type_int_host = typename view_type_int::HostMirror;
    using view_type_int_2d_host = Kokkos::View<int **, layout, Kokkos::HostSpace>;
    using view_type_int_unmanaged = Kokkos::View<int *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using view_type_short = Kokkos::View<short *, memory_space>;
    using view_type_float = Kokkos::View<float *, memory_space>;

    // Using the default exec space for this memory space.
    using execution_space = typename memory_space::execution_space;

    int NextLayer_FirstEpitaxialGrainID;
    view_type_int GrainID_AllLayers, CellType_AllLayers;
    view_type_short LayerID_AllLayers;
    // Substrate inputs from file
    SubstrateInputs _inputs;

    // Constructor for views and view bounds for current layer
    // TODO: CellType is only needed for the current layer, and LayerID is only needed for all layers/no subview is
    // needed. Leaving them as-is for now each with the "AllLayers" and "Current layer/subview slice" to minimize
    // changes to the rest of the code to accomodate. GrainID is initialized to zeros, while others are not initialized
    CellData(int domain_size_all_layers, SubstrateInputs inputs)
        : GrainID_AllLayers(view_type_int("GrainID", domain_size_all_layers))
        , CellType_AllLayers(view_type_int(Kokkos::ViewAllocateWithoutInitializing("CellType"), domain_size_all_layers))
        , LayerID_AllLayers(view_type_short(Kokkos::ViewAllocateWithoutInitializing("LayerID"), domain_size_all_layers))
        , _inputs(inputs) {}

    // Initializes the single active cell and associated active cell data structures for the single grain at the domain
    // center
    void init_substrate(int id, const Grid &grid) {

        // Location of the single grain
        int grainLocationX = floorf(static_cast<float>(grid.nx) / 2.0);
        int grainLocationY = floorf(static_cast<float>(grid.ny) / 2.0);
        int grainLocationZ = floorf(static_cast<float>(grid.nz) / 2.0);

        // Local copies for lambda capture.
        auto CellType_AllLayers_local = CellType_AllLayers;
        auto GrainID_AllLayers_local = GrainID_AllLayers;
        int singleGrainOrientation_local = _inputs.singleGrainOrientation;

        auto policy = Kokkos::RangePolicy<execution_space>(0, grid.domain_size);
        Kokkos::parallel_for(
            "SingleGrainInit", policy, KOKKOS_LAMBDA(const int &index) {
                int coord_x = grid.get_coord_X(index);
                int coord_y = grid.get_coord_Y(index);
                int coord_z = grid.get_coord_Z(index);
                int coord_y_global = coord_y + grid.y_offset;
                if ((coord_x == grainLocationX) && (coord_y_global == grainLocationY) && (coord_z == grainLocationZ)) {
                    CellType_AllLayers_local(index) = FutureActive;
                    GrainID_AllLayers_local(index) = singleGrainOrientation_local + 1;
                }
                else
                    CellType_AllLayers_local(index) = Liquid;
            });
        if ((grainLocationY >= grid.y_offset) && (grainLocationY < grid.y_offset + grid.ny_local))
            std::cout << "Rank " << id << " initialized a grain with orientation " << singleGrainOrientation_local
                      << " initialized at X = " << grainLocationX << ", Y = " << grainLocationY
                      << ", Z = " << grainLocationZ << std::endl;
    }

    // Get the X, Y coordinates and grain ID values for grains at the bottom surface for problem type C
    view_type_int_2d_host getSurfaceActiveCellData(int &SubstrateActCells, int nx, int ny, double RNGSeed) {
        // First get number of substrate grains
        if (_inputs.CustomGrainLocationsIDs)
            SubstrateActCells = _inputs.GrainLocationsX.size();
        else
            SubstrateActCells = std::round(_inputs.FractSurfaceSitesActive * nx * ny);
        // View for storing surface grain locations and IDs
        view_type_int_2d_host ActCellData_Host(Kokkos::ViewAllocateWithoutInitializing("ActCellData_Host"),
                                               SubstrateActCells, 3);
        if (_inputs.CustomGrainLocationsIDs) {
            // Values were already given in the inputs struct, copy these to the view
            for (int n = 0; n < SubstrateActCells; n++) {
                ActCellData_Host(n, 0) = _inputs.GrainLocationsX[n];
                ActCellData_Host(n, 1) = _inputs.GrainLocationsY[n];
                ActCellData_Host(n, 2) = _inputs.GrainIDs[n];
            }
        }
        else {
            // Calls to Xdist(gen) and Y dist(gen) return random locations for grain seeds
            // Since X = 0 and X = nx-1 are the cell centers of the last cells in X, locations are evenly scattered
            // between X = -0.49999 and X = nx - 0.5, as the cells have a half width of 0.5
            std::mt19937_64 gen(RNGSeed);
            std::uniform_real_distribution<double> Xdist(-0.49999, nx - 0.5);
            std::uniform_real_distribution<double> Ydist(-0.49999, ny - 0.5);
            // Randomly locate substrate grain seeds for cells in the interior of this subdomain (at the k = 0 bottom
            // surface)
            for (int n = 0; n < SubstrateActCells; n++) {
                double XLocation = Xdist(gen);
                double YLocation = Ydist(gen);
                // Randomly select integer coordinates between 0 and nx-1 or ny-1
                ActCellData_Host(n, 0) = round(XLocation);
                ActCellData_Host(n, 1) = round(YLocation);
                ActCellData_Host(n, 2) = n + 1; // grain ID for epitaxial seeds must be > 0
            }
        }
        return ActCellData_Host;
    }

    // Initializes cell types and epitaxial Grain ID values where substrate grains are future active cells on the bottom
    // surface of the constrained domain
    void init_substrate(int id, const Grid &grid, double RNGSeed) {

        // Fill the view of cell X, Y, and ID values, updating the number of substrate active cells appropriately
        // TODO: Could generate random numbers on GPU, instead of using host view and copying over - but would also need
        // inputs struct to store device data for grain locations in X, Y, and GrainIDs
        int SubstrateActCells;
        view_type_int_2d_host ActCellData_Host = getSurfaceActiveCellData(SubstrateActCells, grid.nx, grid.ny, RNGSeed);
        // Copy views of substrate grain locations and IDs back to the device
        auto ActCellData = Kokkos::create_mirror_view_and_copy(memory_space(), ActCellData_Host);

        // Start with all cells as liquid prior to locating substrate grain seeds
        // All cells have LayerID = 0 as this is not a multilayer problem
        Kokkos::deep_copy(CellType_AllLayers, Liquid);
        Kokkos::deep_copy(LayerID_AllLayers, 0);

        // Local copies for lambda capture.
        auto CellType_AllLayers_local = CellType_AllLayers;
        auto GrainID_AllLayers_local = GrainID_AllLayers;

        // Determine which grains/active cells belong to which MPI ranks
        auto policy = Kokkos::RangePolicy<execution_space>(0, SubstrateActCells);
        Kokkos::parallel_for(
            "ConstrainedGrainInit", policy, KOKKOS_LAMBDA(const int &n) {
                // What are the X and Y coordinates of this active cell relative to the X and Y bounds of this rank?
                if ((ActCellData(n, 1) >= grid.y_offset) && (ActCellData(n, 1) < grid.y_offset + grid.ny_local)) {
                    // Convert X and Y coordinates to values relative to this MPI rank's grid (Z = 0 for these active
                    // cells, at bottom surface) GrainIDs come from the position on the list of substrate active cells
                    // to avoid reusing the same value
                    int coord_x = ActCellData(n, 0);
                    int coord_y = ActCellData(n, 1) - grid.y_offset;
                    int coord_z = 0;
                    int index = grid.get_1D_index(coord_x, coord_y, coord_z);
                    CellType_AllLayers_local(index) = FutureActive;
                    GrainID_AllLayers_local(index) = ActCellData(n, 2); // assign GrainID > 0 to epitaxial seeds
                }
            });
        // Option to fill empty sites at bottom surface with the grain ID of the nearest grain
        if (_inputs.FillBottomSurface) {
            auto md_policy =
                Kokkos::MDRangePolicy<execution_space, Kokkos::Rank<2, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>(
                    {0, 0}, {grid.nx, grid.ny_local});
            // For cells that are not associated with grain centers, optionally assign them the GrainID of the nearest
            // grain center
            Kokkos::parallel_for(
                "BaseplateGen", md_policy, KOKKOS_LAMBDA(const int coord_x, const int coord_y) {
                    int index_AllLayers = grid.get_1D_index(coord_x, coord_y, 0);
                    if (GrainID_AllLayers_local(index_AllLayers) == 0) {
                        // This cell needs to be assigned a GrainID value
                        // Check each possible baseplate grain center to find the closest one
                        float MinDistanceToThisGrain = grid.nx * grid.ny;
                        int MinDistanceToThisGrain_GrainID = 0;
                        for (int n = 0; n < SubstrateActCells; n++) {
                            // Substrate grain center at coord_x_grain, coord_y_grain - how far is the cell at i,
                            // j+y_offset?
                            int coord_y_grain_global = ActCellData(n, 1);
                            int coord_x_grain = ActCellData(n, 0);
                            int coord_y_global = coord_y + grid.y_offset;
                            float DistanceToThisGrainX = coord_x - coord_x_grain;
                            float DistanceToThisGrainY = coord_y_global - coord_y_grain_global;
                            float DistanceToThisGrain = sqrtf(DistanceToThisGrainX * DistanceToThisGrainX +
                                                              DistanceToThisGrainY * DistanceToThisGrainY);
                            if (DistanceToThisGrain < MinDistanceToThisGrain) {
                                // This is the closest grain center to cell at "CAGridLocation" - update values
                                MinDistanceToThisGrain = DistanceToThisGrain;
                                MinDistanceToThisGrain_GrainID = ActCellData(n, 2);
                            }
                        }
                        // GrainID associated with the closest baseplate grain center
                        GrainID_AllLayers_local(index_AllLayers) = MinDistanceToThisGrain_GrainID;
                        // These cells are also future active cells
                        CellType_AllLayers_local(index_AllLayers) = FutureActive;
                    }
                });
        }
        if (id == 0)
            std::cout << "Number of substrate active cells across all ranks: " << SubstrateActCells << std::endl;
    }

    void init_substrate(int id, const Grid &grid, double RNGSeed, view_type_int NumberOfSolidificationEvents) {

        // Determine the number of cells in the Z direction that are part of the baseplate
        int BaseplateSizeZ = get_baseplate_size_z(id, grid);
        // Generate the baseplate microstructure, or read it from a file, to initialize the grain ID values from Z = 0
        // up to but not including Z = BaseplateTopZ
        if (_inputs.UseSubstrateFile)
            init_baseplate_grainid(id, grid, BaseplateSizeZ);
        else
            init_baseplate_grainid(id, grid, RNGSeed, BaseplateSizeZ);

        // Powder layer extends from Z = PowderBottomZ up to but not including Z = PowderTopZ
        // Bottom of layer is the next coordinate up from the baseplate
        int PowderBottomZ = round((_inputs.BaseplateTopZ - grid.z_min) / grid.deltax) + 1;
        int PowderTopZ = round((grid.z_max_layer[0] - grid.z_min) / grid.deltax) + 1;
        // Generate powder grain structure grain IDs for top of layer 0 if needed (i.e, if the powder layer height is
        // more than zero cells)
        if (PowderTopZ > PowderBottomZ)
            init_powder_grainid(0, id, RNGSeed, grid, PowderBottomZ, PowderTopZ);

        // LayerID starts at -1 for all cells
        Kokkos::deep_copy(LayerID_AllLayers, -1);

        // Initialize cell types and layer IDs based on whether cells will solidify in layer 0 or not
        init_celltype_layerid(0, id, grid, NumberOfSolidificationEvents);
    }

    // Determine the height of the baseplate, in CA cells
    // The baseplate always starts at the simulation bottom (Z coordinate corresponding to z_min, Z index = 0),
    // regardless of whether the first layer melts the cells at the bottom or not. If BaseplateThroughPowder is true,
    // the baseplate microstructure extends through the entire simulation domain in Z (size nz). If
    // BaseplateThroughPowder is false, the baseplate top from the input file is used, or it is assumed that the top of
    // the baseplate is at Z = 0 microns
    int get_baseplate_size_z(int id, const Grid &grid) {
        int BaseplateSizeZ;
        if (_inputs.BaseplateThroughPowder)
            BaseplateSizeZ = grid.nz;
        else {
            BaseplateSizeZ = round((_inputs.BaseplateTopZ - grid.z_min) / grid.deltax) + 1;
            int MaxBaseplateSizeZ = round((grid.z_max_layer[0] - grid.z_min) / grid.deltax) + 1;
            if (BaseplateSizeZ > MaxBaseplateSizeZ) {
                BaseplateSizeZ = MaxBaseplateSizeZ;
                if (id == 0)
                    std::cout << "Warning: The specified location of the baseplate top is located above the "
                                 "temperature data for the first layer of the problem; it will be truncated to stop at "
                                 "the top Z coordinate of the first layer's temperature data"
                              << std::endl;
            }
        }
        if ((id == 0) && (BaseplateSizeZ < 1))
            std::cout << "Warning: no baseplate microstructure will be used, as the temperature data does not overlap "
                         "with the input value for the baseplate top"
                      << std::endl;
        return BaseplateSizeZ;
    }

    // Initializes Grain ID values where the substrate comes from a file
    void init_baseplate_grainid(int id, const Grid &grid, int BaseplateSizeZ) {

        if (id == 0)
            std::cout << "Warning: Reading substrate data from a file will require a vtk file of GrainID values in a "
                         "future release"
                      << std::endl;
        // Assign GrainID values to cells that are part of the substrate - read values from file and initialize using
        // temporary host view
        view_type_int_host GrainID_AllLayers_Host(Kokkos::ViewAllocateWithoutInitializing("GrainID_Host"),
                                                  grid.domain_size_all_layers);
        std::ifstream Substrate;
        Substrate.open(_inputs.SubstrateFileName);
        int Substrate_LowY = grid.y_offset;
        int Substrate_HighY = grid.y_offset + grid.ny_local;
        int nxS, nyS, nzS;
        if ((id == 0) && (BaseplateSizeZ < 1))
            std::cout
                << "Warning: No substrate data from the file be used, as the powder layer extends further in the Z "
                   "direction than the first layer's temperature data"
                << std::endl;
        std::string s;
        getline(Substrate, s);
        std::size_t found = s.find("=");
        std::string str = s.substr(found + 1, s.length() - 1);
        nzS = stoi(str, nullptr, 10);
        getline(Substrate, s);
        found = s.find("=");
        str = s.substr(found + 1, s.length() - 1);
        nyS = stoi(str, nullptr, 10);
        getline(Substrate, s);
        found = s.find("=");
        str = s.substr(found + 1, s.length() - 1);
        nxS = stoi(str, nullptr, 10);
        if ((id == 0) && (nzS < BaseplateSizeZ)) {
            // Do not allow simulation if there is inssufficient substrate data in the specified file
            std::string error =
                "Error: only " + std::to_string(nzS) + " layers of substrate data are present in the file " +
                _inputs.SubstrateFileName + " ; at least " + std::to_string(BaseplateSizeZ) +
                " layers of substrate data are required to simulate the specified solidification problem";
            throw std::runtime_error(error);
        }

        // Assign GrainID values to all cells - cells that will be part of the melt pool footprint may still need their
        // initial GrainID
        for (int k = 0; k < nzS; k++) {
            if (k == BaseplateSizeZ)
                break;
            for (int j = 0; j < nyS; j++) {
                for (int i = 0; i < nxS; i++) {
                    std::string GIDVal;
                    getline(Substrate, GIDVal);
                    if ((j >= Substrate_LowY) && (j < Substrate_HighY)) {
                        int CAGridLocation;
                        CAGridLocation = k * grid.nx * grid.ny_local + i * grid.ny_local + (j - grid.y_offset);
                        GrainID_AllLayers_Host(CAGridLocation) = stoi(GIDVal, nullptr, 10);
                    }
                }
            }
        }
        Substrate.close();
        // Copy GrainIDs read from file to device
        GrainID_AllLayers = Kokkos::create_mirror_view_and_copy(memory_space(), GrainID_AllLayers_Host);
        if (id == 0)
            std::cout << "Substrate file read complete" << std::endl;
    }

    // Initializes Grain ID values where the baseplate is generated using an input grain spacing and a Voronoi
    // Tessellation
    void init_baseplate_grainid(int id, const Grid &grid, double RNGSeed, int BaseplateSizeZ) {

        std::mt19937_64 gen(RNGSeed);

        // Based on the baseplate volume (convert to cubic microns to match units) and the substrate grain spacing,
        // determine the number of baseplate grains
        int BaseplateVolume = grid.nx * grid.ny * BaseplateSizeZ;
        double BaseplateVolume_microns = BaseplateVolume * pow(grid.deltax, 3) * pow(10, 18);
        double SubstrateMeanGrainVolume_microns = pow(_inputs.SubstrateGrainSpacing, 3);
        int NumberOfBaseplateGrains = round(BaseplateVolume_microns / SubstrateMeanGrainVolume_microns);
        // Need at least 1 baseplate grain, cannot have more baseplate grains than cells in the baseplate
        NumberOfBaseplateGrains = std::max(NumberOfBaseplateGrains, 1);
        NumberOfBaseplateGrains = std::min(NumberOfBaseplateGrains, grid.nx * grid.ny * BaseplateSizeZ);
        // TODO: Use device RNG to generate baseplate grain locations, instead of host with copy
        // List of potential grain IDs (starting at 1) - index corresponds to the associated CA cell location
        // Assign positive values for indices 0 through NumberOfBaseplateGrains-1, assign zeros to the remaining indices
        std::vector<int> BaseplateGrainLocations(BaseplateVolume);
        std::vector<int> BaseplateGrainIDs(BaseplateVolume);
        for (int i = 0; i < BaseplateVolume; i++) {
            BaseplateGrainLocations[i] = i;
            if (i < NumberOfBaseplateGrains)
                BaseplateGrainIDs[i] = i + 1;
            else
                BaseplateGrainIDs[i] = 0;
        }
        // Shuffle list of grain IDs
        std::shuffle(BaseplateGrainIDs.begin(), BaseplateGrainIDs.end(), gen);

        // Create views of baseplate grain IDs and locations - copying BaseplateGrainIDs and BaseplateGrainLocations
        // values only for indices where BaseplateGrainIDs(i) != 0 (cells with BaseplateGrainIDs(i) = 0 were not
        // assigned baseplate center locations)
        view_type_int_host BaseplateGrainLocations_Host(
            Kokkos::ViewAllocateWithoutInitializing("BaseplateGrainLocations_Host"), NumberOfBaseplateGrains);
        view_type_int_host BaseplateGrainIDs_Host(Kokkos::ViewAllocateWithoutInitializing("BaseplateGrainIDs_Host"),
                                                  NumberOfBaseplateGrains);
        int count = 0;
        for (int i = 0; i < BaseplateVolume; i++) {
            if (BaseplateGrainIDs[i] != 0) {
                BaseplateGrainLocations_Host(count) = BaseplateGrainLocations[i];
                BaseplateGrainIDs_Host(count) = BaseplateGrainIDs[i];
                count++;
            }
        }

        // Copy baseplate views to the device
        auto BaseplateGrainIDs_Device = Kokkos::create_mirror_view_and_copy(memory_space(), BaseplateGrainIDs_Host);
        auto BaseplateGrainLocations_Device =
            Kokkos::create_mirror_view_and_copy(memory_space(), BaseplateGrainLocations_Host);
        if (id == 0) {
            std::cout << "Baseplate spanning domain coordinates Z = 0 through " << BaseplateSizeZ - 1 << std::endl;
            std::cout << "Number of baseplate grains: " << NumberOfBaseplateGrains << std::endl;
        }

        // Local copies for lambda capture.
        auto GrainID_AllLayers_local = GrainID_AllLayers;

        // First, assign cells that are associated with grain centers the appropriate non-zero GrainID values (assumes
        // GrainID values were initialized to zeros)
        auto policy = Kokkos::RangePolicy<execution_space>(0, NumberOfBaseplateGrains);
        Kokkos::parallel_for(
            "BaseplateInit", policy, KOKKOS_LAMBDA(const int &n) {
                int BaseplateGrainLoc = BaseplateGrainLocations_Device(n);
                // x, y, z associated with baseplate grain "n", at 1D coordinate "BaseplateGrainLoc"
                int coord_z_AllLayers = grid.get_coord_Z_global(BaseplateGrainLoc);
                int coord_y_global = grid.get_coord_Y_global(BaseplateGrainLoc);
                int coord_x = grid.get_coord_X_global(BaseplateGrainLoc);
                if ((coord_y_global >= grid.y_offset) && (coord_y_global < grid.y_offset + grid.ny_local)) {
                    // This grain is associated with a cell on this MPI rank
                    int coord_y = coord_y_global - grid.y_offset;
                    int index_AllLayers = grid.get_1D_index(coord_x, coord_y, coord_z_AllLayers);
                    GrainID_AllLayers_local(index_AllLayers) = BaseplateGrainIDs_Device(n);
                }
            });
        Kokkos::fence();

        auto md_policy =
            Kokkos::MDRangePolicy<execution_space, Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>(
                {0, 0, 0}, {BaseplateSizeZ, grid.nx, grid.ny_local});

        // For cells that are not associated with grain centers, assign them the GrainID of the nearest grain center
        Kokkos::parallel_for(
            "BaseplateGen", md_policy,
            KOKKOS_LAMBDA(const int coord_z_AllLayers, const int coord_x, const int coord_y) {
                int index_AllLayers = grid.get_1D_index(coord_x, coord_y, coord_z_AllLayers);
                if (GrainID_AllLayers_local(index_AllLayers) == 0) {
                    // This cell needs to be assigned a GrainID value
                    // Check each possible baseplate grain center to find the closest one
                    float MinDistanceToThisGrain = grid.nx * grid.ny * BaseplateSizeZ;
                    int MinDistanceToThisGrain_GrainID = 0;
                    for (int n = 0; n < NumberOfBaseplateGrains; n++) {
                        // Baseplate grain center at x_n, y_n, z_n - how far is the cell at i, j+y_offset, k?
                        int coord_z_grain_AllLayers = grid.get_coord_Z_global(BaseplateGrainLocations_Device(n));
                        int coord_y_grain_global = grid.get_coord_Y_global(BaseplateGrainLocations_Device(n));
                        int coord_x_grain = grid.get_coord_X_global(BaseplateGrainLocations_Device(n));
                        int coord_y_global = coord_y + grid.y_offset;
                        float DistanceToThisGrainX = coord_x - coord_x_grain;
                        float DistanceToThisGrainY = coord_y_global - coord_y_grain_global;
                        float DistanceToThisGrainZ = coord_z_grain_AllLayers - coord_z_AllLayers;
                        float DistanceToThisGrain = sqrtf(DistanceToThisGrainX * DistanceToThisGrainX +
                                                          DistanceToThisGrainY * DistanceToThisGrainY +
                                                          DistanceToThisGrainZ * DistanceToThisGrainZ);
                        if (DistanceToThisGrain < MinDistanceToThisGrain) {
                            // This is the closest grain center to cell at "CAGridLocation" - update values
                            MinDistanceToThisGrain = DistanceToThisGrain;
                            MinDistanceToThisGrain_GrainID = BaseplateGrainIDs_Device(n);
                        }
                    }
                    // GrainID associated with the closest baseplate grain center
                    GrainID_AllLayers_local(index_AllLayers) = MinDistanceToThisGrain_GrainID;
                }
            });

        NextLayer_FirstEpitaxialGrainID =
            NumberOfBaseplateGrains + 1; // avoid reusing GrainID in next layer's powder grain structure
        if (id == 0)
            std::cout << "Baseplate grain structure initialized" << std::endl;
    }

    // Each layer's top Z coordinates are seeded with CA-cell sized substrate grains (emulating bulk nucleation
    // alongside the edges of partially melted powder particles). These Z coordinates span PowderBottomZ up to but not
    // including PowderTopZ
    void init_powder_grainid(int layernumber, int id, double RNGSeed, const Grid &grid, int PowderBottomZ,
                             int PowderTopZ) {

        // On all ranks, generate list of powder grain IDs (starting with NextLayer_FirstEpitaxialGrainID, and shuffle
        // them so that their locations aren't sequential and depend on the RNGSeed (different for each layer)
        std::mt19937_64 gen(RNGSeed + layernumber);
        std::uniform_real_distribution<double> dis(0.0, 1.0);

        // TODO: This should be performed on the device, rather than the host
        int PowderLayerHeight = PowderTopZ - PowderBottomZ;
        int PowderLayerCells = grid.nx * grid.ny * PowderLayerHeight;
        int PowderLayerAssignedCells = round(static_cast<double>(PowderLayerCells) * _inputs.PowderActiveFraction);
        std::vector<int> PowderGrainIDs(PowderLayerCells, 0);
        for (int n = 0; n < PowderLayerAssignedCells; n++) {
            PowderGrainIDs[n] = n + NextLayer_FirstEpitaxialGrainID; // assigned a nonzero GrainID
        }
        std::shuffle(PowderGrainIDs.begin(), PowderGrainIDs.end(), gen);
        // Wrap powder layer GrainIDs into an unmanaged view, then copy to the device
        view_type_int_unmanaged PowderGrainIDs_Host(PowderGrainIDs.data(), PowderLayerCells);
        auto PowderGrainIDs_Device = Kokkos::create_mirror_view_and_copy(memory_space(), PowderGrainIDs_Host);
        // Associate powder grain IDs with CA cells in the powder layer
        MPI_Barrier(MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Initializing powder layer for Z = " << PowderBottomZ << " through " << PowderTopZ - 1 << " ("
                      << grid.nx * grid.ny * PowderLayerHeight << " cells)" << std::endl;

        int PowderStart = grid.nx * grid.ny * PowderBottomZ;
        if (id == 0)
            std::cout << "Powder layer has " << PowderLayerAssignedCells
                      << " cells assigned new grain ID values, ranging from " << NextLayer_FirstEpitaxialGrainID
                      << " through " << NextLayer_FirstEpitaxialGrainID + PowderLayerAssignedCells - 1 << std::endl;

        // Iterate over all cells in the powder layer, on each rank loading the powder grain ID data for local cell
        // locations
        auto GrainID_AllLayers_local = GrainID_AllLayers;
        auto powder_policy =
            Kokkos::MDRangePolicy<execution_space, Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>(
                {PowderBottomZ, 0, 0}, {PowderTopZ, grid.nx, grid.ny});
        Kokkos::parallel_for(
            "PowderGrainInit", powder_policy,
            KOKKOS_LAMBDA(const int coord_z_AllLayers, const int coord_x, const int coord_y_global) {
                // Is this powder coordinate in X and Y in bounds for this rank? Is the grain id of this site unassigned
                // (wasn't captured during solidification of the previous layer)?
                if ((coord_y_global >= grid.y_offset) && (coord_y_global < grid.y_offset + grid.ny_local)) {
                    int coord_y = coord_y_global - grid.y_offset;
                    int index_AllLayers = grid.get_1D_index(coord_x, coord_y, coord_z_AllLayers);
                    if (GrainID_AllLayers_local(index_AllLayers) == 0) {
                        int index_AllRanksAllLayers =
                            coord_z_AllLayers * grid.nx * grid.ny + coord_x * grid.ny + coord_y_global;
                        GrainID_AllLayers_local(index_AllLayers) =
                            PowderGrainIDs_Device(index_AllRanksAllLayers - PowderStart);
                    }
                }
            });
        Kokkos::fence();

        // Update NextLayer_FirstEpitaxialGrainID for next layer
        NextLayer_FirstEpitaxialGrainID += PowderLayerAssignedCells;
        MPI_Barrier(MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Initialized powder grain structure for layer " << layernumber << std::endl;
    }

    // Sets up views, powder layer (if necessary), and cell types for the next layer of a multilayer problem
    //*****************************************************************************/
    void init_next_layer(int nextlayernumber, int id, const Grid &grid, double RNGSeed,
                         view_type_int NumberOfSolidificationEvents) {

        // Subviews for the next layer's grain id, layer id, cell type are constructed based on updated layer bound
        // z_layer_bottom
        // Powder layer extends from Z = PowderBottomZ (1 cell above the top of the previous layer) up to but not
        // including Z = PowderTopZ
        int PowderBottomZ = round((grid.z_max_layer[nextlayernumber - 1] - grid.z_min) / grid.deltax) + 1;
        int PowderTopZ = round((grid.z_max_layer[nextlayernumber] - grid.z_min) / grid.deltax) + 1;
        if (!(_inputs.BaseplateThroughPowder))
            init_powder_grainid(nextlayernumber, id, RNGSeed, grid, PowderBottomZ, PowderTopZ);

        // Initialize active cell data structures and nuclei locations for the next layer "layernumber + 1"
        init_celltype_layerid(nextlayernumber, id, grid, NumberOfSolidificationEvents);
    }

    //*****************************************************************************/
    // Initializes cells for the current layer as either solid (don't resolidify) or tempsolid (will melt and
    // resolidify)
    void init_celltype_layerid(int layernumber, int id, const Grid &grid, view_type_int NumberOfSolidificationEvents) {

        int MeltPoolCellCount;
        // Local copies for lambda capture.
        auto CellType_AllLayers_local = CellType_AllLayers;
        auto LayerID_AllLayers_local = LayerID_AllLayers;

        auto policy = Kokkos::RangePolicy<execution_space>(0, grid.domain_size);

        Kokkos::parallel_reduce(
            "CellTypeInitSolidRM", policy,
            KOKKOS_LAMBDA(const int &index, int &local_count) {
                int index_AllLayers = index + grid.z_layer_bottom * grid.nx * grid.ny_local;
                if (NumberOfSolidificationEvents(index) > 0) {
                    CellType_AllLayers_local(index_AllLayers) = TempSolid;
                    LayerID_AllLayers_local(index_AllLayers) = layernumber;
                    local_count++;
                }
                else
                    CellType_AllLayers_local(index_AllLayers) = Solid;
            },
            MeltPoolCellCount);
        int TotalMeltPoolCellCount;
        MPI_Reduce(&MeltPoolCellCount, &TotalMeltPoolCellCount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Number of cells across all ranks to undergo solidification at least once: "
                      << TotalMeltPoolCellCount << std::endl;
    }

    // Stores/returns the volume fraction of nucleated grains to the console
    // Moved from CAfunctions.hpp
    float calcVolFractionNucleated(int id, const Grid &grid) {

        // For interior cells, add the number of cells that underwent melting/solidification and the number of cells
        // with sub-zero grain IDs
        int MeltedCells_Local = 0;
        int NucleatedGrainCells_Local = 0;
        // Local copies for lambda capture.
        auto GrainID_AllLayers_local = GrainID_AllLayers;
        auto LayerID_AllLayers_local = LayerID_AllLayers;
        Kokkos::parallel_reduce(
            "NumSolidifiedCells", grid.domain_size,
            KOKKOS_LAMBDA(const int index, int &update_meltcount, int &update_nucleatecount) {
                int coord_y = grid.get_coord_Y(index);
                // Is this Y coordinate in the halo region? If so, do not increment counter
                bool InHaloRegion = false;
                if (((coord_y == 0) && (!grid.at_south_boundary)) ||
                    ((coord_y == grid.ny_local - 1) && (!grid.at_north_boundary)))
                    InHaloRegion = true;
                if ((GrainID_AllLayers_local(index) < 0) && (!InHaloRegion))
                    update_nucleatecount++;
                if ((LayerID_AllLayers_local(index) != -1) && (!InHaloRegion))
                    update_meltcount++;
            },
            MeltedCells_Local, NucleatedGrainCells_Local);

        // Reduce the values by summing over all ranks
        int MeltedCells_Global, NucleatedGrainCells_Global;
        MPI_Allreduce(&MeltedCells_Local, &MeltedCells_Global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&NucleatedGrainCells_Local, &NucleatedGrainCells_Global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        // Calculate nucleated grain fraction
        float VolFractionNucleated =
            static_cast<float>(NucleatedGrainCells_Global) / static_cast<float>(MeltedCells_Global);
        if (id == 0)
            std::cout << "The fraction of the solidified material consisting of nucleated grains is "
                      << VolFractionNucleated << std::endl;
        return VolFractionNucleated;
    }

    // Take a view consisting of data for all layers, and return a subview of the same type consisting of just the cells
    // corresponding to the current layer of a multilayer problem
    auto getGrainIDSubview(const Grid &grid) const { return Kokkos::subview(GrainID_AllLayers, grid.layer_range); }
    auto getLayerIDSubview(const Grid &grid) const { return Kokkos::subview(LayerID_AllLayers, grid.layer_range); }
    auto getCellTypeSubview(const Grid &grid) const { return Kokkos::subview(CellType_AllLayers, grid.layer_range); }
};

#endif
