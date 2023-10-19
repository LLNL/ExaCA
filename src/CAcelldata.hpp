// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_CELLDATA_HPP
#define EXACA_CELLDATA_HPP

#include "CAconfig.hpp"
#include "CAfunctions.hpp"
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
    using view_type_int_unmanaged = Kokkos::View<int *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using view_type_short = Kokkos::View<short *, memory_space>;
    using view_type_float = Kokkos::View<float *, memory_space>;

    int NextLayer_FirstEpitaxialGrainID, BottomOfCurrentLayer, TopOfCurrentLayer;
    std::pair<int, int> LayerRange;
    view_type_int GrainID_AllLayers, CellType_AllLayers;
    view_type_short LayerID_AllLayers;
    // Substrate inputs from file
    SubstrateInputs _inputs;

    // Constructor for views and view bounds for current layer
    // TODO: CellType is only needed for the current layer, and LayerID is only needed for all layers/no subview is
    // needed. Leaving them as-is for now each with the "AllLayers" and "Current layer/subview slice" to minimize
    // changes to the rest of the code to accomodate. GrainID is initialized to zeros, while others are not initialized
    CellData(int DomainSize_AllLayers, int DomainSize, int nx, int ny_local, int z_layer_bottom, SubstrateInputs inputs)
        : GrainID_AllLayers(view_type_int("GrainID", DomainSize_AllLayers))
        , CellType_AllLayers(view_type_int(Kokkos::ViewAllocateWithoutInitializing("CellType"), DomainSize_AllLayers))
        , LayerID_AllLayers(view_type_short(Kokkos::ViewAllocateWithoutInitializing("LayerID"), DomainSize_AllLayers))
        , _inputs(inputs) {

        BottomOfCurrentLayer = z_layer_bottom * nx * ny_local;
        TopOfCurrentLayer = BottomOfCurrentLayer + DomainSize;
        LayerRange = std::make_pair(BottomOfCurrentLayer, TopOfCurrentLayer);
    }

    // Initializes the single active cell and associated active cell data structures for the single grain at the domain
    // center
    void init_substrate(int id, int nx, int ny, int nz, int ny_local, int y_offset, int DomainSize) {

        // Location of the single grain
        int grainLocationX = floorf(static_cast<float>(nx) / 2.0);
        int grainLocationY = floorf(static_cast<float>(ny) / 2.0);
        int grainLocationZ = floorf(static_cast<float>(nz) / 2.0);

        auto CellType_AllLayers_local = CellType_AllLayers;
        auto GrainID_AllLayers_local = GrainID_AllLayers;
        int singleGrainOrientation_local = _inputs.singleGrainOrientation;
        Kokkos::parallel_for(
            "SingleGrainInit", DomainSize, KOKKOS_LAMBDA(const int &index) {
                int coord_x = getCoordX(index, nx, ny_local);
                int coord_y = getCoordY(index, nx, ny_local);
                int coord_z = getCoordZ(index, nx, ny_local);
                int coord_y_global = coord_y + y_offset;
                if ((coord_x == grainLocationX) && (coord_y_global == grainLocationY) && (coord_z == grainLocationZ)) {
                    CellType_AllLayers_local(index) = FutureActive;
                    GrainID_AllLayers_local(index) = singleGrainOrientation_local + 1;
                }
                else
                    CellType_AllLayers_local(index) = Liquid;
            });
        if ((grainLocationY >= y_offset) && (grainLocationY < y_offset + ny_local))
            std::cout << "Rank " << id << " initialized a grain with orientation " << singleGrainOrientation_local
                      << " initialized at X = " << grainLocationX << ", Y = " << grainLocationY
                      << ", Z = " << grainLocationZ << std::endl;
    }

    // Initializes cell types and epitaxial Grain ID values where substrate grains are future active cells on the bottom
    // surface of the constrained domain
    void init_substrate(int id, int ny_local, int nx, int ny, int y_offset, double RNGSeed) {

        // Calls to Xdist(gen) and Y dist(gen) return random locations for grain seeds
        // Since X = 0 and X = nx-1 are the cell centers of the last cells in X, locations are evenly scattered between
        // X = -0.49999 and X = nx - 0.5, as the cells have a half width of 0.5
        std::mt19937_64 gen(RNGSeed);
        std::uniform_real_distribution<double> Xdist(-0.49999, nx - 0.5);
        std::uniform_real_distribution<double> Ydist(-0.49999, ny - 0.5);

        // Determine number of active cells from the fraction of sites active and the number of sites on the bottom
        // domain surface
        int SubstrateActCells = std::round(_inputs.FractSurfaceSitesActive * nx * ny);

        // On all ranks, generate active site locations - same list on every rank
        // TODO: Generate random numbers on GPU, instead of using host view and copying over - ensure that locations are
        // in the same order every time based on the RNGSeed
        view_type_int_host ActCellX_Host(Kokkos::ViewAllocateWithoutInitializing("ActCellX_Host"), SubstrateActCells);
        view_type_int_host ActCellY_Host(Kokkos::ViewAllocateWithoutInitializing("ActCellY_Host"), SubstrateActCells);

        // Randomly locate substrate grain seeds for cells in the interior of this subdomain (at the k = 0 bottom
        // surface)
        for (int n = 0; n < SubstrateActCells; n++) {
            double XLocation = Xdist(gen);
            double YLocation = Ydist(gen);
            // Randomly select integer coordinates between 0 and nx-1 or ny-1
            ActCellX_Host(n) = round(XLocation);
            ActCellY_Host(n) = round(YLocation);
        }

        // Copy views of substrate grain locations back to the device
        auto ActCellX_Device = Kokkos::create_mirror_view_and_copy(memory_space(), ActCellX_Host);
        auto ActCellY_Device = Kokkos::create_mirror_view_and_copy(memory_space(), ActCellY_Host);

        // Start with all cells as liquid prior to locating substrate grain seeds
        // All cells have LayerID = 0 as this is not a multilayer problem
        Kokkos::deep_copy(CellType_AllLayers, Liquid);
        Kokkos::deep_copy(LayerID_AllLayers, 0);
        auto CellType_AllLayers_local = CellType_AllLayers;
        auto GrainID_AllLayers_local = GrainID_AllLayers;

        // Determine which grains/active cells belong to which MPI ranks
        Kokkos::parallel_for(
            "ConstrainedGrainInit", SubstrateActCells, KOKKOS_LAMBDA(const int &n) {
                // What are the X and Y coordinates of this active cell relative to the X and Y bounds of this rank?
                if ((ActCellY_Device(n) >= y_offset) && (ActCellY_Device(n) < y_offset + ny_local)) {
                    // Convert X and Y coordinates to values relative to this MPI rank's grid (Z = 0 for these active
                    // cells, at bottom surface) GrainIDs come from the position on the list of substrate active cells
                    // to avoid reusing the same value
                    int coord_x = ActCellX_Device(n);
                    int coord_y = ActCellY_Device(n) - y_offset;
                    int coord_z = 0;
                    int index = get1Dindex(coord_x, coord_y, coord_z, nx, ny_local);
                    CellType_AllLayers_local(index) = FutureActive;
                    GrainID_AllLayers_local(index) = n + 1; // assign GrainID > 0 to epitaxial seeds
                }
            });
        if (id == 0)
            std::cout << "Number of substrate active cells across all ranks: " << SubstrateActCells << std::endl;
    }

    void init_substrate(int nx, int ny, int nz, int DomainSize, double *ZMaxLayer, double ZMin, double deltax,
                        int ny_local, int y_offset, int z_layer_bottom, int id, double RNGSeed,
                        view_type_int NumberOfSolidificationEvents) {

        // Determine the number of cells in the Z direction that are part of the baseplate
        int BaseplateSizeZ = get_baseplate_size_z(id, nz, ZMin, ZMaxLayer, deltax);

        // Generate the baseplate microstructure, or read it from a file, to initialize the grain ID values from Z = 0
        // up to but not including Z = BaseplateTopZ
        if (_inputs.UseSubstrateFile)
            init_baseplate_grainid(nz, nx, ny_local, y_offset, id, BaseplateSizeZ);
        else
            init_baseplate_grainid(nx, ny, ny_local, y_offset, id, deltax, RNGSeed, BaseplateSizeZ);

        // Powder layer extends from Z = PowderBottomZ up to but not including Z = PowderTopZ
        // Bottom of layer is the next coordinate up from the baseplate
        int PowderBottomZ = round((_inputs.BaseplateTopZ - ZMin) / deltax) + 1;
        int PowderTopZ = round((ZMaxLayer[0] - ZMin) / deltax) + 1;
        // Generate powder grain structure grain IDs for top of layer 0 if needed (i.e, if the powder layer height is
        // more than zero cells)
        if (PowderTopZ > PowderBottomZ)
            init_powder_grainid(0, nx, ny, ny_local, y_offset, id, RNGSeed, PowderBottomZ, PowderTopZ);

        // LayerID starts at -1 for all cells
        Kokkos::deep_copy(LayerID_AllLayers, -1);

        // Initialize cell types and layer IDs based on whether cells will solidify in layer 0 or not
        init_celltype_layerid(0, nx, ny_local, DomainSize, NumberOfSolidificationEvents, id, z_layer_bottom);
    }

    // Determine the height of the baseplate, in CA cells
    // The baseplate always starts at the simulation bottom (Z coordinate corresponding to ZMin, Z index = 0),
    // regardless of whether the first layer melts the cells at the bottom or not. If BaseplateThroughPowder is true,
    // the baseplate microstructure extends through the entire simulation domain in Z (size nz). If
    // BaseplateThroughPowder is false, the baseplate top from the input file is used, or it is assumed that the top of
    // the baseplate is at Z = 0 microns
    int get_baseplate_size_z(int id, int nz, double ZMin, double *ZMaxLayer, double deltax) {
        int BaseplateSizeZ;
        if (_inputs.BaseplateThroughPowder)
            BaseplateSizeZ = nz;
        else {
            BaseplateSizeZ = round((_inputs.BaseplateTopZ - ZMin) / deltax) + 1;
            int MaxBaseplateSizeZ = round((ZMaxLayer[0] - ZMin) / deltax) + 1;
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
    void init_baseplate_grainid(int nz, int nx, int ny_local, int y_offset, int id, int BaseplateSizeZ) {

        if (id == 0)
            std::cout << "Warning: Reading substrate data from a file will require a vtk file of GrainID values in a "
                         "future release"
                      << std::endl;
        // Assign GrainID values to cells that are part of the substrate - read values from file and initialize using
        // temporary host view
        view_type_int_host GrainID_AllLayers_Host(Kokkos::ViewAllocateWithoutInitializing("GrainID_Host"),
                                                  nx * ny_local * nz);
        std::ifstream Substrate;
        Substrate.open(_inputs.SubstrateFileName);
        int Substrate_LowY = y_offset;
        int Substrate_HighY = y_offset + ny_local;
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
                        CAGridLocation = k * nx * ny_local + i * ny_local + (j - y_offset);
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
    void init_baseplate_grainid(int nx, int ny, int ny_local, int y_offset, int id, double deltax, double RNGSeed,
                                int BaseplateSizeZ) {

        std::mt19937_64 gen(RNGSeed);

        // Based on the baseplate volume (convert to cubic microns to match units) and the substrate grain spacing,
        // determine the number of baseplate grains
        int BaseplateVolume = nx * ny * BaseplateSizeZ;
        double BaseplateVolume_microns = BaseplateVolume * pow(deltax, 3) * pow(10, 18);
        double SubstrateMeanGrainVolume_microns = pow(_inputs.SubstrateGrainSpacing, 3);
        int NumberOfBaseplateGrains = round(BaseplateVolume_microns / SubstrateMeanGrainVolume_microns);
        // Need at least 1 baseplate grain, cannot have more baseplate grains than cells in the baseplate
        NumberOfBaseplateGrains = std::max(NumberOfBaseplateGrains, 1);
        NumberOfBaseplateGrains = std::min(NumberOfBaseplateGrains, nx * ny * BaseplateSizeZ);
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

        // First, assign cells that are associated with grain centers the appropriate non-zero GrainID values (assumes
        // GrainID values were initialized to zeros)
        auto GrainID_AllLayers_local = GrainID_AllLayers;
        Kokkos::parallel_for(
            "BaseplateInit", NumberOfBaseplateGrains, KOKKOS_LAMBDA(const int &n) {
                int BaseplateGrainLoc = BaseplateGrainLocations_Device(n);
                // x, y, z associated with baseplate grain "n", at 1D coordinate "BaseplateGrainLoc"
                int coord_z_AllLayers = getCoordZ(BaseplateGrainLoc, nx, ny);
                int coord_y_global = getCoordY(BaseplateGrainLoc, nx, ny);
                int coord_x = getCoordX(BaseplateGrainLoc, nx, ny);
                if ((coord_y_global >= y_offset) && (coord_y_global < y_offset + ny_local)) {
                    // This grain is associated with a cell on this MPI rank
                    int coord_y = coord_y_global - y_offset;
                    int index_AllLayers = get1Dindex(coord_x, coord_y, coord_z_AllLayers, nx, ny_local);
                    GrainID_AllLayers_local(index_AllLayers) = BaseplateGrainIDs_Device(n);
                }
            });
        Kokkos::fence();
        // For cells that are not associated with grain centers, assign them the GrainID of the nearest grain center
        Kokkos::parallel_for(
            "BaseplateGen",
            Kokkos::MDRangePolicy<Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>(
                {0, 0, 0}, {BaseplateSizeZ, nx, ny_local}),
            KOKKOS_LAMBDA(const int coord_z_AllLayers, const int coord_x, const int coord_y) {
                int index_AllLayers = get1Dindex(coord_x, coord_y, coord_z_AllLayers, nx, ny_local);
                if (GrainID_AllLayers_local(index_AllLayers) == 0) {
                    // This cell needs to be assigned a GrainID value
                    // Check each possible baseplate grain center to find the closest one
                    float MinDistanceToThisGrain = nx * ny * BaseplateSizeZ;
                    int MinDistanceToThisGrain_GrainID = 0;
                    for (int n = 0; n < NumberOfBaseplateGrains; n++) {
                        // Baseplate grain center at x_n, y_n, z_n - how far is the cell at i, j+y_offset, k?
                        int coord_z_grain_AllLayers = getCoordZ(BaseplateGrainLocations_Device(n), nx, ny);
                        int coord_y_grain_global = getCoordY(BaseplateGrainLocations_Device(n), nx, ny);
                        int coord_x_grain = getCoordX(BaseplateGrainLocations_Device(n), nx, ny);
                        int coord_y_global = coord_y + y_offset;
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
    void init_powder_grainid(int layernumber, int nx, int ny, int ny_local, int y_offset, int id, double RNGSeed,
                             int PowderBottomZ, int PowderTopZ) {

        // On all ranks, generate list of powder grain IDs (starting with NextLayer_FirstEpitaxialGrainID, and shuffle
        // them so that their locations aren't sequential and depend on the RNGSeed (different for each layer)
        std::mt19937_64 gen(RNGSeed + layernumber);
        std::uniform_real_distribution<double> dis(0.0, 1.0);

        // TODO: This should be performed on the device, rather than the host
        int PowderLayerHeight = PowderTopZ - PowderBottomZ;
        int PowderLayerCells = nx * ny * PowderLayerHeight;
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
                      << nx * ny * PowderLayerHeight << " cells)" << std::endl;

        int PowderStart = nx * ny * PowderBottomZ;
        int PowderEnd = nx * ny * PowderTopZ;
        if (id == 0)
            std::cout << "Powder layer has " << PowderLayerAssignedCells
                      << " cells assigned new grain ID values, ranging from " << NextLayer_FirstEpitaxialGrainID
                      << " through " << NextLayer_FirstEpitaxialGrainID + PowderLayerAssignedCells - 1 << std::endl;
        auto GrainID_AllLayers_local = GrainID_AllLayers;
        Kokkos::parallel_for(
            "PowderGrainInit", Kokkos::RangePolicy<>(PowderStart, PowderEnd),
            KOKKOS_LAMBDA(const int &index_global_AllLayers) {
                int coord_z_AllLayers = getCoordZ(index_global_AllLayers, nx, ny);
                int coord_y_global = getCoordY(index_global_AllLayers, nx, ny);
                int coord_y = coord_y_global - y_offset;
                int coord_x = getCoordX(index_global_AllLayers, nx, ny);
                int index_AllLayers = get1Dindex(coord_x, coord_y, coord_z_AllLayers, nx, ny_local);
                // Is this powder coordinate in X and Y in bounds for this rank? Is the grain id of this site unassigned
                // (wasn't captured during solidification of the previous layer)?
                if ((coord_y_global >= y_offset) && (coord_y_global < y_offset + ny_local) &&
                    (GrainID_AllLayers_local(index_AllLayers) == 0))
                    GrainID_AllLayers_local(index_AllLayers) =
                        PowderGrainIDs_Device(index_global_AllLayers - PowderStart);
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
    void init_next_layer(int nextlayernumber, int id, int nx, int ny, int ny_local, int y_offset, int z_layer_bottom,
                         int DomainSize, double RNGSeed, double ZMin, double *ZMaxLayer, double deltax,
                         view_type_int NumberOfSolidificationEvents) {

        // Subviews for the next layer's grain id, layer id, cell type are constructed based on updated layer bound
        // z_layer_bottom
        BottomOfCurrentLayer = z_layer_bottom * nx * ny_local;
        TopOfCurrentLayer = BottomOfCurrentLayer + DomainSize;
        LayerRange = std::make_pair(BottomOfCurrentLayer, TopOfCurrentLayer);

        // Powder layer extends from Z = PowderBottomZ (1 cell above the top of the previous layer) up to but not
        // including Z = PowderTopZ
        int PowderBottomZ = round((ZMaxLayer[nextlayernumber - 1] - ZMin) / deltax) + 1;
        int PowderTopZ = round((ZMaxLayer[nextlayernumber] - ZMin) / deltax) + 1;
        if (!(_inputs.BaseplateThroughPowder))
            init_powder_grainid(nextlayernumber, nx, ny, ny_local, y_offset, id, RNGSeed, PowderBottomZ, PowderTopZ);

        // Initialize active cell data structures and nuclei locations for the next layer "layernumber + 1"
        init_celltype_layerid(nextlayernumber, nx, ny_local, DomainSize, NumberOfSolidificationEvents, id,
                              z_layer_bottom);
    }

    //*****************************************************************************/
    // Initializes cells for the current layer as either solid (don't resolidify) or tempsolid (will melt and
    // resolidify)
    void init_celltype_layerid(int layernumber, int nx, int ny_local, int DomainSize,
                               view_type_int NumberOfSolidificationEvents, int id, int z_layer_bottom) {

        int MeltPoolCellCount;
        auto CellType_AllLayers_local = CellType_AllLayers;
        auto LayerID_AllLayers_local = LayerID_AllLayers;
        Kokkos::parallel_reduce(
            "CellTypeInitSolidRM", DomainSize,
            KOKKOS_LAMBDA(const int &index, int &local_count) {
                int index_AllLayers = index + z_layer_bottom * nx * ny_local;
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

    // Take a view consisting of data for all layers, and return a subview of the same type consisting of just the cells
    // corresponding to the current layer of a multilayer problem
    auto getGrainIDSubview() { return Kokkos::subview(GrainID_AllLayers, LayerRange); }
    auto getLayerIDSubview() { return Kokkos::subview(LayerID_AllLayers, LayerRange); }
    auto getCellTypeSubview() { return Kokkos::subview(CellType_AllLayers, LayerRange); }
};

#endif
