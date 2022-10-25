// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_FUNCTIONS_HPP
#define EXACA_FUNCTIONS_HPP

#include <CAtypes.hpp>
#include <Kokkos_Core.hpp>
//*****************************************************************************/
// Inline functions

// Get the orientation of a grain from a given grain ID and the number of possible orientations
KOKKOS_INLINE_FUNCTION int getGrainOrientation(int MyGrainID, int NGrainOrientations) {
    int MyOrientation = (abs(MyGrainID) - 1) % NGrainOrientations;
    return MyOrientation;
}
//*****************************************************************************/
int YMPSlicesCalc(int p, int ny, int np);
int YOffsetCalc(int p, int ny, int np);
void AddGhostNodes(int NeighborRank_North, int NeighborRank_South, int &MyYSlices, int &MyYOffset);
double MaxVal(double TestVec3[6], int NVals);
void InitialDecomposition(int id, int np, int &NeighborRank_North, int &NeighborRank_South, bool &AtNorthBoundary,
                          bool &AtSouthBoundary);
ViewF_H MisorientationCalc(int NumberOfOrientations, ViewF_H GrainUnitVector, int dir);

#endif
