// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
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
int YMPSlicesCalc(int p, int ny, int np);
int YOffsetCalc(int p, int ny, int np);
void AddGhostNodes(int NeighborRank_North, int NeighborRank_South, int &MyYSlices, int &MyYOffset);
double MaxVal(double TestVec3[6], int NVals);
void InitialDecomposition(int id, int np, int &NeighborRank_North, int &NeighborRank_South, bool &AtNorthBoundary,
                          bool &AtSouthBoundary);
ViewF_H MisorientationCalc(int NumberOfOrientations, ViewF_H GrainUnitVector, int dir);
float calcVolFractionNucleated(int id, int nx, int MyYSlices, int LocalDomainSize, ViewI LayerID, ViewI GrainID,
                               bool AtNorthBoundary, bool AtSouthBoundary);

#endif
