// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_FUNCTIONS_HPP
#define EXACA_FUNCTIONS_HPP

#include <Kokkos_Core.hpp>

//*****************************************************************************/
// Inline functions
// Get the orientation of a grain from a given grain ID and the number of possible orientations
// By default, start indexing at 0 (GrainID of 1 has Orientation number 0), optionally starting at 1
// GrainID of 0 is a special case - has orientation 0 no matter what (only used when reconstructing grain ID in
// getGrainID)
KOKKOS_INLINE_FUNCTION int getGrainOrientation(int MyGrainID, int NGrainOrientations, bool StartAtZero = true) {
    int MyOrientation;
    if (MyGrainID == 0)
        MyOrientation = 0;
    else {
        MyOrientation = (abs(MyGrainID) - 1) % NGrainOrientations;
        if (!(StartAtZero))
            MyOrientation++;
    }
    return MyOrientation;
}
// Get the repeat number for the orientation of a grain with a given grain ID
// 1, 2, 3... or -1, -2, -3...
KOKKOS_INLINE_FUNCTION int getGrainNumber(int MyGrainID, int NGrainOrientations) {
    int MyGrainNumber = (abs(MyGrainID) - 1) / NGrainOrientations + 1;
    if (MyGrainID < 0)
        MyGrainNumber = -MyGrainNumber;
    return MyGrainNumber;
}
// Get the grain ID from the repeat number of a grain from a given grain ID and the number of possible orientations
KOKKOS_INLINE_FUNCTION int getGrainID(int NGrainOrientations, int MyGrainOrientation, int MyGrainNumber) {
    int MyGrainID = NGrainOrientations * (abs(MyGrainNumber) - 1) + MyGrainOrientation;
    if (MyGrainNumber < 0)
        MyGrainID = -MyGrainID;
    return MyGrainID;
}
//*****************************************************************************/
// Create a view of size "NumberOfOrientation" of the misorientation of each possible grain orientation with the X, Y,
// or Z directions (dir = 0, 1, or 2, respectively)
template <typename ViewTypeFloatHost>
auto MisorientationCalc(int NumberOfOrientations, ViewTypeFloatHost GrainUnitVector, int dir) {

    ViewTypeFloatHost GrainMisorientation(Kokkos::ViewAllocateWithoutInitializing("GrainMisorientation"),
                                          NumberOfOrientations);
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

#endif
