// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_NUCLEATION_HPP
#define EXACA_NUCLEATION_HPP

#include "CAcelldata.hpp"
#include "CAconfig.hpp"
#include "CAtemperature.hpp"
#include "CAtypes.hpp"
#include "mpi.h"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cmath>
#include <random>
#include <string>
#include <vector>

// Data regarding nucleation events in the domain
template <typename MemorySpace>
struct Nucleation {

    using memory_space = MemorySpace;
    using view_type_int = Kokkos::View<int *, memory_space>;
    using view_type_int_host = typename view_type_int::HostMirror;

    // Four counters tracked here:
    // 1. Nuclei_WholeDomain - tracks all nuclei (regardless of whether an event would be possible based on the layer ID
    // and cell type), used for Grain ID assignment to ensure that no Grain ID get reused - same on all MPI ranks
    // 2. x - the subset of Nuclei_ThisLayer that are located within the bounds of a
    // given MPI rank that may possibly occur (nuclei locations are associated with a liquid cell with a layer ID that
    // matches this layer number). Starts at 0 each layer
    // 3. NucleationCounter - the number of nucleation events that have actually either failed or succeeded on a given
    // MPI rank. Starts at 0 each layer
    // 4. SuccessfulNucleationCounter - the number of nucleation events that have successfully occured at a given point
    // in the simulation of a layer. Starts at 0 each layer.
    int Nuclei_WholeDomain, PossibleNuclei, NucleationCounter, SuccessfulNucleationCounter;

    // The probability that a given liquid site will be a potential nucleus location
    double BulkProb;

    // The time steps at which nucleation events will occur in the given layer, on the host
    view_type_int_host NucleationTimes_Host;
    // The locations and grain IDs of the potential nuclei in the layer
    view_type_int NucleiLocations, NucleiGrainID;

    // Constructor - initialize CA views using input for initial guess at number of possible events
    // Default is that no nucleation has occurred to this point, optional input argument to start with the counter at a
    // specific value
    Nucleation(const int EstimatedNuclei_ThisRankThisLayer, const double NMax, const double deltax,
               int NumPriorNuclei = 0)
        : NucleationTimes_Host(view_type_int_host(Kokkos::ViewAllocateWithoutInitializing("NucleationTimes_Host"),
                                                  EstimatedNuclei_ThisRankThisLayer))
        , NucleiLocations(view_type_int(Kokkos::ViewAllocateWithoutInitializing("NucleiLocations"),
                                        EstimatedNuclei_ThisRankThisLayer))
        , NucleiGrainID(view_type_int(Kokkos::ViewAllocateWithoutInitializing("NucleiGrainID"),
                                      EstimatedNuclei_ThisRankThisLayer)) {
        Nuclei_WholeDomain = NumPriorNuclei;
        BulkProb = NMax * deltax * deltax * deltax;
        resetNucleiCounters(); // start counters at 0
    }

    // Reset the nuclei counters prior to nuclei initialization for a layer
    void resetNucleiCounters() {
        // Init appropriate counters for the layer to 0 - possible nuclei for this layer will be calculated in this
        // function
        PossibleNuclei = 0;
        NucleationCounter = 0;
        SuccessfulNucleationCounter = 0;
    }

    // Initialize nucleation site locations, GrainID values, and time at which nucleation events will potentially occur,
    // accounting for multiple possible nucleation events in cells that melt and solidify multiple times
    template <class... Params>
    void placeNuclei(Temperature<memory_space> &temperature, double RNGSeed, int layernumber, int nx, int ny,
                     int nz_layer, double dTN, double dTsigma, int ny_local, int y_offset, int, int id,
                     bool AtNorthBoundary, bool AtSouthBoundary) {

        // TODO: convert this subroutine into kokkos kernels, rather than copying data back to the host, and nucleation
        // data back to the device again. This is currently performed on the device due to heavy usage of standard
        // library algorithm functions Copy temperature data into temporary host views for this subroutine
        auto MaxSolidificationEvents_Host =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), temperature.MaxSolidificationEvents);
        auto NumberOfSolidificationEvents_Host =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), temperature.NumberOfSolidificationEvents);
        auto LayerTimeTempHistory_Host =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), temperature.LayerTimeTempHistory);

        // Use new RNG seed for each layer
        std::mt19937_64 generator(RNGSeed + layernumber);
        // Uniform distribution for nuclei location assignment
        std::uniform_real_distribution<double> Xdist(-0.49999, nx - 0.5);
        std::uniform_real_distribution<double> Ydist(-0.49999, ny - 0.5);
        std::uniform_real_distribution<double> Zdist(-0.49999, nz_layer - 0.5);
        // Gaussian distribution of nucleation undercooling
        std::normal_distribution<double> Gdistribution(dTN, dTsigma);

        // Max number of nucleated grains in this layer
        int Nuclei_ThisLayerSingle = BulkProb * (nx * ny * nz_layer); // equivalent to Nuclei_ThisLayer if no remelting
        // Multiplier for the number of nucleation events per layer, based on the number of solidification events
        int NucleiMultiplier = MaxSolidificationEvents_Host(layernumber);
        int Nuclei_ThisLayer = Nuclei_ThisLayerSingle * NucleiMultiplier;
        // Temporary vectors for storing nucleated grain IDs and undercooling values
        // Nuclei Grain ID are assigned to avoid reusing values from previous layers
        std::vector<int> NucleiGrainID_WholeDomain_V(Nuclei_ThisLayer);
        std::vector<double> NucleiUndercooling_WholeDomain_V(Nuclei_ThisLayer);
        // Views for storing potential nucleated grain coordinates
        view_type_int_host NucleiX(Kokkos::ViewAllocateWithoutInitializing("NucleiX"), Nuclei_ThisLayer);
        view_type_int_host NucleiY(Kokkos::ViewAllocateWithoutInitializing("NucleiY"), Nuclei_ThisLayer);
        view_type_int_host NucleiZ(Kokkos::ViewAllocateWithoutInitializing("NucleiZ"), Nuclei_ThisLayer);

        for (int meltevent = 0; meltevent < NucleiMultiplier; meltevent++) {
            for (int n = 0; n < Nuclei_ThisLayerSingle; n++) {
                int NEvent = meltevent * Nuclei_ThisLayerSingle + n;
                // Generate possible nuclei locations
                double NucleiX_unrounded = Xdist(generator);
                double NucleiY_unrounded = Ydist(generator);
                double NucleiZ_unrounded = Zdist(generator);
                // Round these coordinates so they're associated with a specific cell on the grid
                NucleiX(NEvent) = std::round(NucleiX_unrounded);
                NucleiY(NEvent) = std::round(NucleiY_unrounded);
                NucleiZ(NEvent) = std::round(NucleiZ_unrounded);
                // Assign each nuclei a Grain ID (negative values used for nucleated grains) and an undercooling
                NucleiGrainID_WholeDomain_V[NEvent] = -(Nuclei_WholeDomain + NEvent + 1); // avoid using grain ID 0
                NucleiUndercooling_WholeDomain_V[NEvent] = Gdistribution(generator);
            }
        }

        // Shuffle these vectors to make sure the same grain IDs and undercooling don't end up in the same spots each
        // layer
        std::shuffle(NucleiGrainID_WholeDomain_V.begin(), NucleiGrainID_WholeDomain_V.end(), generator);
        std::shuffle(NucleiUndercooling_WholeDomain_V.begin(), NucleiUndercooling_WholeDomain_V.end(), generator);
        if ((id == 0) && (Nuclei_ThisLayer > 0))
            std::cout << "Range of Grain IDs from which layer " << layernumber
                      << " nucleation events were selected: " << Nuclei_WholeDomain + 1 << " through "
                      << Nuclei_WholeDomain + Nuclei_ThisLayer << std::endl;
        // Update number of nuclei counter for whole domain based on the number of nuclei in this layer
        Nuclei_WholeDomain += Nuclei_ThisLayer;

        // Loop through nuclei for this layer - each MPI rank storing the nucleation events that are possible (i.e,
        // nucleation event is associated with a CA cell on that MPI rank's subdomain, the cell is liquid type, and the
        // cell is associated with the current layer of the multilayer problem) Don't put nuclei in "ghost" cells -
        // those nucleation events occur on other ranks and the existing halo exchange functionality will handle this
        std::vector<int> NucleiGrainID_MyRank_V(Nuclei_ThisLayer), NucleiLocation_MyRank_V(Nuclei_ThisLayer),
            NucleationTimes_MyRank_V(Nuclei_ThisLayer);
        for (int meltevent = 0; meltevent < NucleiMultiplier; meltevent++) {
            for (int n = 0; n < Nuclei_ThisLayerSingle; n++) {
                int NEvent = meltevent * Nuclei_ThisLayerSingle + n;
                if (((NucleiY(NEvent) > y_offset) || (AtSouthBoundary)) &&
                    ((NucleiY(NEvent) < y_offset + ny_local - 1) || (AtNorthBoundary))) {
                    // Convert 3D location (using global X and Y coordinates) into a 1D location (using local X and Y
                    // coordinates) for the possible nucleation event, both as relative to the bottom of this layer
                    int NucleiLocation_ThisLayer =
                        get1Dindex(NucleiX(NEvent), NucleiY(NEvent) - y_offset, NucleiZ(NEvent), nx, ny_local);
                    // Criteria for placing a nucleus - whether or not this nuclei is associated with a solidification
                    // event
                    if (meltevent < NumberOfSolidificationEvents_Host(NucleiLocation_ThisLayer)) {
                        // Nucleation event is possible - cell undergoes solidification at least once, this nucleation
                        // event is associated with one of the time periods during which the associated cell undergoes
                        // solidification
                        NucleiLocation_MyRank_V[PossibleNuclei] = NucleiLocation_ThisLayer;

                        float CritTimeStep_ThisEvent =
                            LayerTimeTempHistory_Host(NucleiLocation_ThisLayer, meltevent, 1);
                        float UndercoolingChange_ThisEvent =
                            LayerTimeTempHistory_Host(NucleiLocation_ThisLayer, meltevent, 2);
                        float TimeToNucUnd = round(CritTimeStep_ThisEvent + NucleiUndercooling_WholeDomain_V[NEvent] /
                                                                                UndercoolingChange_ThisEvent);
                        if (CritTimeStep_ThisEvent > TimeToNucUnd)
                            TimeToNucUnd = CritTimeStep_ThisEvent;
                        NucleationTimes_MyRank_V[PossibleNuclei] = round(TimeToNucUnd);
                        // Assign this cell the potential nucleated grain ID
                        NucleiGrainID_MyRank_V[PossibleNuclei] = NucleiGrainID_WholeDomain_V[NEvent];
                        // Increment counter on this MPI rank
                        PossibleNuclei++;
                    }
                }
            }
        }

        // How many nucleation events are actually possible (associated with a cell in this layer that will undergo
        // solidification)?
        int PossibleNuclei_AllRanksThisLayer;
        MPI_Reduce(&PossibleNuclei, &PossibleNuclei_AllRanksThisLayer, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Number of potential nucleation events in layer " << layernumber << " : "
                      << PossibleNuclei_AllRanksThisLayer << std::endl;

        // Now that the number of nucleation events on each rank is known, resize these vectors
        NucleiLocation_MyRank_V.resize(PossibleNuclei);
        NucleiGrainID_MyRank_V.resize(PossibleNuclei);
        NucleationTimes_MyRank_V.resize(PossibleNuclei);

        // Sort the list of time steps at which nucleation occurs, keeping the time steps paired with the corresponding
        // locations for nucleation events and grain IDs
        std::vector<std::tuple<int, int, int>> NucleationTimeLocID;
        NucleationTimeLocID.reserve(PossibleNuclei);
        for (int n = 0; n < PossibleNuclei; n++) {
            NucleationTimeLocID.push_back(
                std::make_tuple(NucleationTimes_MyRank_V[n], NucleiLocation_MyRank_V[n], NucleiGrainID_MyRank_V[n]));
        }
        // Sorting from low to high
        std::sort(NucleationTimeLocID.begin(), NucleationTimeLocID.end());

        // With PossibleNuclei_ThisRankThisLayer now known, resize views appropriately
        // Resize nucleation views now that PossibleNuclei_ThisRank is known for all MPI ranks
        Kokkos::resize(NucleationTimes_Host, PossibleNuclei);
        Kokkos::resize(NucleiLocations, PossibleNuclei);
        Kokkos::resize(NucleiGrainID, PossibleNuclei);

        // Create temporary view to store nucleation locations, grain ID data initialized on the host
        // NucleationTimes_H are stored using a host view that is passed to Nucleation subroutine and later used - don't
        // need a temporary host view
        view_type_int_host NucleiLocations_Host(Kokkos::ViewAllocateWithoutInitializing("NucleiLocations_Host"),
                                                PossibleNuclei);
        view_type_int_host NucleiGrainID_Host(Kokkos::ViewAllocateWithoutInitializing("NucleiGrainID_Host"),
                                              PossibleNuclei);
        for (int n = 0; n < PossibleNuclei; n++) {
            NucleationTimes_Host(n) = std::get<0>(NucleationTimeLocID[n]);
            NucleiLocations_Host(n) = std::get<1>(NucleationTimeLocID[n]);
            NucleiGrainID_Host(n) = std::get<2>(NucleationTimeLocID[n]);
        }
        // Copy nucleation data to the device
        NucleiLocations = Kokkos::create_mirror_view_and_copy(memory_space(), NucleiLocations_Host);
        NucleiGrainID = Kokkos::create_mirror_view_and_copy(memory_space(), NucleiGrainID_Host);
    }

    // Compute velocity from local undercooling.
    // functional form is assumed to be cubic if not explicitly given in input file
    void nucleate_grain(int cycle, CellData<memory_space> &cellData, int, int, int, view_type_int SteeringVector,
                        view_type_int numSteer_G) {

        auto CellType = cellData.getCellTypeSubview();
        auto GrainID = cellData.getGrainIDSubview();

        // Is there nucleation left in this layer to check?
        if (NucleationCounter < PossibleNuclei) {
            // Is there at least one potential nucleation event on this rank, at this time step?
            if (cycle == NucleationTimes_Host(NucleationCounter)) {
                bool NucleationCheck = true;
                int FirstEvent = NucleationCounter; // first potential nucleation event to check
                // Are there any other nucleation events this time step to check?
                while (NucleationCheck) {
                    NucleationCounter++;
                    // If the previous nucleation event was the last one for this layer of the simulation, exit loop
                    if (NucleationCounter == PossibleNuclei)
                        break;
                    // If the next nucleation event corresponds to a future time step, finish check
                    if (cycle != NucleationTimes_Host(NucleationCounter))
                        NucleationCheck = false;
                }
                int LastEvent = NucleationCounter;
                // parallel_reduce checks each potential nucleation event this time step (FirstEvent, up to but not
                // including LastEvent)
                int NucleationThisDT = 0; // return number of successful event from parallel_reduce
                auto NucleiLocations_local = NucleiLocations;
                auto NucleiGrainID_local = NucleiGrainID;
                // Launch kokkos kernel - check if the corresponding CA cell location is liquid
                Kokkos::parallel_reduce(
                    "NucleiUpdateLoop", Kokkos::RangePolicy<>(FirstEvent, LastEvent),
                    KOKKOS_LAMBDA(const int NucleationCounter_Device, int &update) {
                        int NucleationEventLocation = NucleiLocations_local(NucleationCounter_Device);
                        int update_val =
                            FutureActive; // added to steering vector to become a new active cell as part of cellcapture
                        int old_val = Liquid;
                        int OldCellTypeValue =
                            Kokkos::atomic_compare_exchange(&CellType(NucleationEventLocation), old_val, update_val);
                        if (OldCellTypeValue == Liquid) {
                            // Successful nucleation event - atomic update of cell type, proceeded if the atomic
                            // exchange is successful (cell was liquid) Add future active cell location to steering
                            // vector and change cell type, assign new Grain ID
                            GrainID(NucleationEventLocation) = NucleiGrainID_local(NucleationCounter_Device);
                            SteeringVector(Kokkos::atomic_fetch_add(&numSteer_G(0), 1)) = NucleationEventLocation;
                            // This undercooled liquid cell is now a nuclei (no nuclei are in the ghost nodes - halo
                            // exchange routine GhostNodes1D or GhostNodes2D is used to fill these)
                            update++;
                        }
                    },
                    NucleationThisDT);
                // Update the number of successful nuclei counter with the number of successful nucleation events from
                // this time step (NucleationThisDT)
                SuccessfulNucleationCounter += NucleationThisDT;
            }
        }
    }
};

#endif
