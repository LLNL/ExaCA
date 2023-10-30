// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAinitialize.hpp"

#include "mpi.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <regex>

void checkPowderOverflow(int nx, int ny, int LayerHeight, int NumberOfLayers, Inputs &inputs) {

    // Check to make sure powder grain density is compatible with the number of powder sites
    // If this problem type includes a powder layer of some grain density, ensure that integer overflow won't occur when
    // assigning powder layer GrainIDs
    if (!(inputs.substrate.BaseplateThroughPowder)) {
        long int NumCellsPowderLayers =
            (long int)(nx) * (long int)(ny) * (long int)(LayerHeight) * (long int)(NumberOfLayers - 1);
        long int NumAssignedCellsPowderLayers =
            std::lround(round(static_cast<double>(NumCellsPowderLayers) * inputs.substrate.PowderActiveFraction));
        if (NumAssignedCellsPowderLayers > INT_MAX)
            throw std::runtime_error("Error: A smaller value for powder density is required to avoid potential integer "
                                     "overflow when assigning powder layer GrainID");
    }
}
