// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_CELLDATA_HPP
#define EXACA_CELLDATA_HPP

#include "CAconfig.hpp"
#include "CAfunctions.hpp"
#include "CAtypes.hpp"
#include "mpi.h"

#include <Kokkos_Core.hpp>

#include <nlohmann/json.hpp>

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
    using view_type_float = Kokkos::View<float *, memory_space>;

    int NextLayer_FirstEpitaxialGrainID, BottomOfCurrentLayer, TopOfCurrentLayer;
    std::pair<int, int> LayerRange;
    view_type_int GrainID_AllLayers, CellType_AllLayers, LayerID_AllLayers;

    // Constructor for views and view bounds for current layer
    // TODO: CellType is only needed for the current layer, and LayerID is only needed for all layers/no subview is
    // needed. Leaving them as-is for now each with the "AllLayers" and "Current layer/subview slice" to minimize
    // changes to the rest of the code to accomodate. GrainID is initialized to zeros, while others are not initialized
    CellData(int LocalDomainSize, int LocalActiveDomainSize, int nx, int MyYSlices, int ZBound_Low)
        : GrainID_AllLayers(view_type_int("GrainID", LocalDomainSize))
        , CellType_AllLayers(view_type_int(Kokkos::ViewAllocateWithoutInitializing("CellType"), LocalDomainSize))
        , LayerID_AllLayers(view_type_int(Kokkos::ViewAllocateWithoutInitializing("LayerID"), LocalDomainSize)) {

        BottomOfCurrentLayer = ZBound_Low * nx * MyYSlices;
        TopOfCurrentLayer = BottomOfCurrentLayer + LocalActiveDomainSize;
        LayerRange = std::make_pair(BottomOfCurrentLayer, TopOfCurrentLayer);
    }

    // Initializes cell types and epitaxial Grain ID values where substrate grains are active cells on the bottom
    // surface of the constrained domain. Also initialize active cell data structures associated with the substrate
    // grains
    void init_substrate(int id, double FractSurfaceSitesActive, int MyYSlices, int nx, int ny, int MyYOffset,
                        NList NeighborX, NList NeighborY, NList NeighborZ, view_type_float GrainUnitVector,
                        int NGrainOrientations, view_type_float DiagonalLength, view_type_float DOCenter,
                        view_type_float CritDiagonalLength, double RNGSeed) {

        // Calls to Xdist(gen) and Y dist(gen) return random locations for grain seeds
        // Since X = 0 and X = nx-1 are the cell centers of the last cells in X, locations are evenly scattered between
        // X = -0.49999 and X = nx - 0.5, as the cells have a half width of 0.5
        std::mt19937_64 gen(RNGSeed);
        std::uniform_real_distribution<double> Xdist(-0.49999, nx - 0.5);
        std::uniform_real_distribution<double> Ydist(-0.49999, ny - 0.5);

        // Determine number of active cells from the fraction of sites active and the number of sites on the bottom
        // domain surface
        int SubstrateActCells = std::round(FractSurfaceSitesActive * nx * ny);

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
        auto ActCellX_Device = Kokkos::create_mirror_view_and_copy(device_memory_space(), ActCellX_Host);
        auto ActCellY_Device = Kokkos::create_mirror_view_and_copy(device_memory_space(), ActCellY_Host);

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
                if ((ActCellY_Device(n) >= MyYOffset) && (ActCellY_Device(n) < MyYOffset + MyYSlices)) {
                    // Convert X and Y coordinates to values relative to this MPI rank's grid (Z = 0 for these active
                    // cells, at bottom surface) GrainIDs come from the position on the list of substrate active cells
                    // to avoid reusing the same value
                    int GlobalX = ActCellX_Device(n);
                    int LocalY = ActCellY_Device(n) - MyYOffset;
                    int D3D1ConvPosition = GlobalX * MyYSlices + LocalY;
                    CellType_AllLayers_local(D3D1ConvPosition) = Active;
                    GrainID_AllLayers_local(D3D1ConvPosition) = n + 1; // assign GrainID > 0 to epitaxial seeds
                    // Initialize active cell data structures
                    int GlobalY = LocalY + MyYOffset;
                    int GlobalZ = 0;
                    // Initialize new octahedron
                    createNewOctahedron(D3D1ConvPosition, DiagonalLength, DOCenter, GlobalX, GlobalY, GlobalZ);

                    // The orientation for the new grain will depend on its Grain ID
                    int MyOrientation = getGrainOrientation(n + 1, NGrainOrientations);
                    float cx = GlobalX + 0.5;
                    float cy = GlobalY + 0.5;
                    float cz = GlobalZ + 0.5;
                    // Calculate critical values at which this active cell leads to the activation of a neighboring
                    // liquid cell. Octahedron center and cell center overlap for octahedra created as part of a new
                    // grain
                    calcCritDiagonalLength(D3D1ConvPosition, cx, cy, cz, cx, cy, cz, NeighborX, NeighborY, NeighborZ,
                                           MyOrientation, GrainUnitVector, CritDiagonalLength);
                }
            });
        if (id == 0)
            std::cout << "Number of substrate active cells across all ranks: " << SubstrateActCells << std::endl;
    }

    void init_substrate(std::string SubstrateFileName, bool UseSubstrateFile, bool BaseplateThroughPowder,
                        bool PowderFirstLayer, int nx, int ny, int nz, int LayerHeight, int LocalActiveDomainSize,
                        double *ZMaxLayer, double ZMin, double deltax, int MyYSlices, int MyYOffset, int ZBound_Low,
                        int id, double RNGSeed, double SubstrateGrainSpacing, double PowderActiveFraction,
                        view_type_int NumberOfSolidificationEvents) {

        // Generate the baseplate microstructure, or read it from a file, to initialize the grain ID values
        if (UseSubstrateFile)
            init_baseplate_grainid(SubstrateFileName, nz, nx, MyYSlices, MyYOffset, id, ZMin, ZMaxLayer,
                                   BaseplateThroughPowder, PowderFirstLayer, LayerHeight, deltax);
        else
            init_baseplate_grainid(SubstrateGrainSpacing, nx, ny, ZMin, ZMaxLayer, MyYSlices, MyYOffset, id, deltax,
                                   RNGSeed, nz, BaseplateThroughPowder, PowderFirstLayer, LayerHeight);

        // Optionally generate powder grain structure grain IDs for top of layer 0
        if (PowderFirstLayer)
            init_powder_grainid(0, nx, ny, LayerHeight, ZMaxLayer, ZMin, deltax, MyYSlices, MyYOffset, id, RNGSeed,
                                PowderActiveFraction);

        // LayerID starts at -1 for all cells
        Kokkos::deep_copy(LayerID_AllLayers, -1);

        // Initialize cell types and layer IDs based on whether cells will solidify in layer 0 or not
        init_celltype_layerid(0, nx, MyYSlices, LocalActiveDomainSize, NumberOfSolidificationEvents, id, ZBound_Low);
    }

    // Determine the height of the baseplate, in CA cells
    // The baseplate always starts at the simulation bottom (Z coordinate corresponding to ZMin), regardless of whether
    // the first layer melts the cells at the bottom or not If BaseplateThroughPowder is true, the baseplate
    // microstructure extends through the entire simulation domain in Z If BaseplateThroughPowder is false and
    // PowderFirstLayer is true, the baseplate extends to the Z coordinate that is LayerHeight cells short of the top of
    // the first layer (the final LayerHeight cells of the first layer will be initialized using the powder options
    // specified)) If BaseplateThroughPowder is false and PowderFirstLayer is false, the baseplate extends to the top of
    // the first layer
    int get_baseplate_size_z(int nz, double *ZMaxLayer, double ZMin, double deltax, int LayerHeight,
                             bool BaseplateThroughPowder, bool PowderFirstLayer) {
        int BaseplateSizeZ;
        if (BaseplateThroughPowder)
            BaseplateSizeZ = nz;
        else {
            if (PowderFirstLayer)
                BaseplateSizeZ = round((ZMaxLayer[0] - ZMin) / deltax) + 1 - LayerHeight;

            else
                BaseplateSizeZ = round((ZMaxLayer[0] - ZMin) / deltax) + 1;
        }
        return BaseplateSizeZ;
    }

    // Initializes Grain ID values where the substrate comes from a file
    void init_baseplate_grainid(std::string SubstrateFileName, int nz, int nx, int MyYSlices, int MyYOffset, int id,
                                double ZMin, double *ZMaxLayer, bool BaseplateThroughPowder, bool PowderFirstLayer,
                                int LayerHeight, double deltax) {

        // Assign GrainID values to cells that are part of the substrate - read values from file and initialize using
        // temporary host view
        view_type_int_host GrainID_AllLayers_Host(Kokkos::ViewAllocateWithoutInitializing("GrainID_Host"),
                                                  nx * MyYSlices * nz);
        std::ifstream Substrate;
        Substrate.open(SubstrateFileName);
        int Substrate_LowY = MyYOffset;
        int Substrate_HighY = MyYOffset + MyYSlices;
        int nxS, nyS, nzS;
        int BaseplateSizeZ = get_baseplate_size_z(nz, ZMaxLayer, ZMin, deltax, LayerHeight, BaseplateThroughPowder,
                                                  PowderFirstLayer); // in CA cells
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
                SubstrateFileName + " ; at least " + std::to_string(BaseplateSizeZ) +
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
                        CAGridLocation = k * nx * MyYSlices + i * MyYSlices + (j - MyYOffset);
                        GrainID_AllLayers_Host(CAGridLocation) = stoi(GIDVal, nullptr, 10);
                    }
                }
            }
        }
        Substrate.close();
        // Copy GrainIDs read from file to device
        GrainID_AllLayers = Kokkos::create_mirror_view_and_copy(device_memory_space(), GrainID_AllLayers_Host);
        if (id == 0)
            std::cout << "Substrate file read complete" << std::endl;
    }

    // Initializes Grain ID values where the baseplate is generated using an input grain spacing and a Voronoi
    // Tessellation
    void init_baseplate_grainid(float SubstrateGrainSpacing, int nx, int ny, double ZMin, double *ZMaxLayer,
                                int MyYSlices, int MyYOffset, int id, double deltax, double RNGSeed, int nz,
                                bool BaseplateThroughPowder, bool PowderFirstLayer, int LayerHeight) {

        // Number of cells to assign GrainID
        int BaseplateSizeZ = get_baseplate_size_z(nz, ZMaxLayer, ZMin, deltax, LayerHeight, BaseplateThroughPowder,
                                                  PowderFirstLayer); // in CA cells
        if ((id == 0) && (BaseplateSizeZ < 1))
            std::cout
                << "Warning: no baseplate microstructure will be used, as the powder layer extends further in the Z "
                   "direction than the first layer's temperature data"
                << std::endl;
        std::mt19937_64 gen(RNGSeed);

        // Based on the baseplate volume (convert to cubic microns to match units) and the substrate grain spacing,
        // determine the number of baseplate grains
        int BaseplateVolume = nx * ny * BaseplateSizeZ;
        double BaseplateVolume_microns = BaseplateVolume * pow(deltax, 3) * pow(10, 18);
        double SubstrateMeanGrainVolume_microns = pow(SubstrateGrainSpacing, 3);
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
        auto BaseplateGrainIDs_Device =
            Kokkos::create_mirror_view_and_copy(device_memory_space(), BaseplateGrainIDs_Host);
        auto BaseplateGrainLocations_Device =
            Kokkos::create_mirror_view_and_copy(device_memory_space(), BaseplateGrainLocations_Host);
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
                int z_n = BaseplateGrainLoc / (nx * ny);
                int Rem = BaseplateGrainLoc % (nx * ny);
                int x_n = Rem / ny;
                int y_n = Rem % ny;
                if ((y_n >= MyYOffset) && (y_n < MyYOffset + MyYSlices)) {
                    // This grain is associated with a cell on this MPI rank
                    int CAGridLocation = z_n * nx * MyYSlices + x_n * MyYSlices + (y_n - MyYOffset);
                    GrainID_AllLayers_local(CAGridLocation) = BaseplateGrainIDs_Device(n);
                }
            });
        Kokkos::fence();
        // For cells that are not associated with grain centers, assign them the GrainID of the nearest grain center
        Kokkos::parallel_for(
            "BaseplateGen",
            Kokkos::MDRangePolicy<Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>(
                {0, 0, 0}, {BaseplateSizeZ, nx, MyYSlices}),
            KOKKOS_LAMBDA(const int k, const int i, const int j) {
                int CAGridLocation = k * nx * MyYSlices + i * MyYSlices + j;
                if (GrainID_AllLayers_local(CAGridLocation) == 0) {
                    // This cell needs to be assigned a GrainID value
                    // Check each possible baseplate grain center to find the closest one
                    float MinDistanceToThisGrain = nx * ny * BaseplateSizeZ;
                    int MinDistanceToThisGrain_GrainID = 0;
                    for (int n = 0; n < NumberOfBaseplateGrains; n++) {
                        // Baseplate grain center at x_n, y_n, z_n - how far is the cell at i, j+MyYOffset, k?
                        int z_n = BaseplateGrainLocations_Device(n) / (nx * ny);
                        int Rem = BaseplateGrainLocations_Device(n) % (nx * ny);
                        int x_n = Rem / ny;
                        int y_n = Rem % ny;
                        float DistanceToThisGrainX = i - x_n;
                        float DistanceToThisGrainY = (j + MyYOffset) - y_n;
                        float DistanceToThisGrainZ = k - z_n;
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
                    GrainID_AllLayers_local(CAGridLocation) = MinDistanceToThisGrain_GrainID;
                }
            });

        NextLayer_FirstEpitaxialGrainID =
            NumberOfBaseplateGrains + 1; // avoid reusing GrainID in next layer's powder grain structure
        if (id == 0)
            std::cout << "Baseplate grain structure initialized" << std::endl;
    }

    // Each layer's top Z coordinates are seeded with CA-cell sized substrate grains (emulating bulk nucleation
    // alongside the edges of partially melted powder particles)
    void init_powder_grainid(int layernumber, int nx, int ny, int LayerHeight, double *ZMaxLayer, double ZMin,
                             double deltax, int MyYSlices, int MyYOffset, int id, double RNGSeed,
                             double PowderActiveFraction) {

        // On all ranks, generate list of powder grain IDs (starting with NextLayer_FirstEpitaxialGrainID, and shuffle
        // them so that their locations aren't sequential and depend on the RNGSeed (different for each layer)
        std::mt19937_64 gen(RNGSeed + layernumber);
        std::uniform_real_distribution<double> dis(0.0, 1.0);

        // TODO: This should be performed on the device, rather than the host
        int PowderLayerCells = nx * ny * LayerHeight;
        int PowderLayerAssignedCells = round(static_cast<double>(PowderLayerCells) * PowderActiveFraction);
        std::vector<int> PowderGrainIDs(PowderLayerCells, 0);
        for (int n = 0; n < PowderLayerAssignedCells; n++) {
            PowderGrainIDs[n] = n + NextLayer_FirstEpitaxialGrainID; // assigned a nonzero GrainID
        }
        std::shuffle(PowderGrainIDs.begin(), PowderGrainIDs.end(), gen);
        // Wrap powder layer GrainIDs into an unmanaged view, then copy to the device
        view_type_int_unmanaged PowderGrainIDs_Host(PowderGrainIDs.data(), PowderLayerCells);
        auto PowderGrainIDs_Device = Kokkos::create_mirror_view_and_copy(device_memory_space(), PowderGrainIDs_Host);
        // Associate powder grain IDs with CA cells in the powder layer
        // Use bounds from temperature field for this layer to determine which cells are part of the powder
        int PowderTopZ = round((ZMaxLayer[layernumber] - ZMin) / deltax) + 1;
        int PowderBottomZ = PowderTopZ - LayerHeight;
        MPI_Barrier(MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Initializing powder layer for Z = " << PowderBottomZ << " through " << PowderTopZ - 1 << " ("
                      << nx * ny * (PowderTopZ - PowderBottomZ) << " cells)" << std::endl;

        int PowderStart = nx * ny * PowderBottomZ;
        int PowderEnd = nx * ny * PowderTopZ;
        if (id == 0)
            std::cout << "Powder layer has " << PowderLayerAssignedCells
                      << " cells assigned new grain ID values, ranging from " << NextLayer_FirstEpitaxialGrainID
                      << " through " << NextLayer_FirstEpitaxialGrainID + PowderLayerAssignedCells - 1 << std::endl;
        auto GrainID_AllLayers_local = GrainID_AllLayers;
        Kokkos::parallel_for(
            "PowderGrainInit", Kokkos::RangePolicy<>(PowderStart, PowderEnd), KOKKOS_LAMBDA(const int &n) {
                int GlobalZ = n / (nx * ny);
                int Rem = n % (nx * ny);
                int GlobalX = Rem / ny;
                int GlobalY = Rem % ny;
                // Is this powder coordinate in X and Y in bounds for this rank? Is the grain id of this site unassigned
                // (wasn't captured during solidification of the previous layer)?
                int GlobalD3D1ConvPosition = GlobalZ * nx * MyYSlices + GlobalX * MyYSlices + (GlobalY - MyYOffset);
                if ((GlobalY >= MyYOffset) && (GlobalY < MyYOffset + MyYSlices) &&
                    (GrainID_AllLayers_local(GlobalD3D1ConvPosition) == 0))
                    GrainID_AllLayers_local(GlobalD3D1ConvPosition) = PowderGrainIDs_Device(n - PowderStart);
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
    void init_next_layer(int nextlayernumber, int id, int nx, int ny, int MyYSlices, int MyYOffset, int ZBound_Low,
                         int LayerHeight, int LocalActiveDomainSize, bool UseSubstrateFile, bool BaseplateThroughPowder,
                         double ZMin, double *ZMaxLayer, double deltax, double RNGSeed, double PowderActiveFraction,
                         view_type_int NumberOfSolidificationEvents) {

        // Subviews for the next layer's grain id, layer id, cell type are constructed based on updated layer bound
        // ZBound_Low
        BottomOfCurrentLayer = ZBound_Low * nx * MyYSlices;
        TopOfCurrentLayer = BottomOfCurrentLayer + LocalActiveDomainSize;
        LayerRange = std::make_pair(BottomOfCurrentLayer, TopOfCurrentLayer);

        // If the baseplate was initialized from a substrate grain spacing, initialize powder layer grain structure
        // for the next layer "layernumber + 1" Otherwise, the entire substrate (baseplate + powder) was read from a
        // file, and the powder layers have already been initialized
        if ((!(UseSubstrateFile)) && (!(BaseplateThroughPowder)))
            init_powder_grainid(nextlayernumber, nx, ny, LayerHeight, ZMaxLayer, ZMin, deltax, MyYSlices, MyYOffset, id,
                                RNGSeed, PowderActiveFraction);

        // Initialize active cell data structures and nuclei locations for the next layer "layernumber + 1"
        init_celltype_layerid(nextlayernumber, nx, MyYSlices, LocalActiveDomainSize, NumberOfSolidificationEvents, id,
                              ZBound_Low);
    }

    //*****************************************************************************/
    // Initializes cells for the current layer as either solid (don't resolidify) or tempsolid (will melt and
    // resolidify)
    void init_celltype_layerid(int layernumber, int nx, int MyYSlices, int LocalActiveDomainSize,
                               view_type_int NumberOfSolidificationEvents, int id, int ZBound_Low) {

        int MeltPoolCellCount;
        auto CellType_AllLayers_local = CellType_AllLayers;
        auto LayerID_AllLayers_local = LayerID_AllLayers;
        Kokkos::parallel_reduce(
            "CellTypeInitSolidRM", LocalActiveDomainSize,
            KOKKOS_LAMBDA(const int &D3D1ConvPosition, int &local_count) {
                int GlobalD3D1ConvPosition = D3D1ConvPosition + ZBound_Low * nx * MyYSlices;
                if (NumberOfSolidificationEvents(D3D1ConvPosition) > 0) {
                    CellType_AllLayers_local(GlobalD3D1ConvPosition) = TempSolid;
                    LayerID_AllLayers_local(GlobalD3D1ConvPosition) = layernumber;
                    local_count++;
                }
                else
                    CellType_AllLayers_local(GlobalD3D1ConvPosition) = Solid;
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
