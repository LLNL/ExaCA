// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_TEMPS_HPP
#define EXACA_TEMPS_HPP

#include "CAconfig.hpp"
#include "CAparsefiles.hpp"
#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Interfacial repsonse function with various functional forms.
struct Temperatures {

    // Maximum number of times each CA cell in a given layer undergoes solidification
    ViewI MaxSolidificationEvents;
    // For each cell in the current layer (index 1), and each time solidification happens (index 2), hold the values
    // that will be used for MeltTimeStep, CritTimeStep, and UndercoolingChange (index 3)
    ViewF3D LayerTimeTempHistory;
    // The number of times that each CA cell will undergo solidification during this layer
    ViewI NumberOfSolidificationEvents;
    // A counter for the number of times each CA cell has undergone solidification so far this layer
    ViewI SolidificationEventCounter;
    // The current undercooling of a CA cell (if superheated liquid or hasn't undergone solidification yet, equals 0)
    ViewF UndercoolingCurrent;

    // Constructor creates views with size based on the inputs and calls init routines based on the simulation type
    Temperatures(std::string SimulationType, int id, const int NumberOfLayers, const int LocalActiveDomainSize) {

        Kokkos::resize(MaxSolidificationEvents, NumberOfLayers);
        int MaxSolidificationEvents_ThisLayer;
        if ((SimulationType == "C") || (SimulationType == "SingleGrain"))
            MaxSolidificationEvents_ThisLayer = 1;
        else if (SimulationType == "S")
            MaxSolidificationEvents_ThisLayer = calcMaxSolidificationEvents_Spots();
        else if (SimulationType == "R")
            MaxSolidificationEvents_ThisLayer = calcMaxSolidificationEvents_InputDataFromFile();
        Kokkos::resize(LayerTimeTempHistory, LocalActiveDomainSize, MaxSolidificationEvents_ThisLayer, 3);
        Kokkos::resize(NumberOfSolidificationEvents, LocalActiveDomainSize);
        Kokkos::resize(SolidificationEventCounter, LocalActiveDomainSize);
        Kokkos::resize(UndercoolingCurrent, LocalActiveDomainSize);
        if ((SimulationType == "C") || (SimulationType == "SingleGrain"))
            initUnidirectional();
        else if (SimulationType == "S")
            initSpots();
        else if (SimulationType == "R")
            initInputDataFromFile();
    }

    // Used for calculating the maximum number of times a spot will undergo solidification in each layer
    int calcMaxSolidificationEvents_Spots() {}

    // Reset local cell undercooling to 0
    KOKKOS_INLINE_FUNCTION
    void reset_undercooling(const int CellIndex) const { UndercoolingCurrent(CellIndex) = 0.0; }

    // Update local cell undercooling for the current melt-resolidification event
    KOKKOS_INLINE_FUNCTION
    void update_undercooling(const int CellIndex) const {
        UndercoolingCurrent(CellIndex) += LayerTimeTempHistory(CellIndex, SolidificationEventCounter(CellIndex), 2);
    }
};

#endif
