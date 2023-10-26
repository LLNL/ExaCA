// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAfunctions.hpp"

#include "mpi.h"

#include <cmath>
#include <iostream>

/*************************** FUNCTIONS CALLED THROUGH MAIN SUBROUTINES ***************************/

//*****************************************************************************/

// Create a view of size "NumberOfOrientation" of the misorientation of each possible grain orientation with the X, Y,
// or Z directions (dir = 0, 1, or 2, respectively)
ViewF_H MisorientationCalc(int NumberOfOrientations, ViewF_H GrainUnitVector, int dir) {

    ViewF_H GrainMisorientation(Kokkos::ViewAllocateWithoutInitializing("GrainMisorientation"), NumberOfOrientations);
    // Find the smallest possible misorientation between the specified direction, and this grain orientations' 6
    // possible 001 directions (where 62.7 degrees is the largest possible misorientation between two 001 directions
    // for a cubic crystal system)
    for (int n = 0; n < NumberOfOrientations; n++) {
        float MisorientationAngleMin = 62.7;
        for (int ll = 0; ll < 3; ll++) {
            float Misorientation = std::abs((180 / M_PI) * acos(GrainUnitVector(9 * n + 3 * ll + dir)));
            if (Misorientation < MisorientationAngleMin) {
                MisorientationAngleMin = Misorientation;
            }
        }
        GrainMisorientation(n) = MisorientationAngleMin;
    }
    return GrainMisorientation;
}
