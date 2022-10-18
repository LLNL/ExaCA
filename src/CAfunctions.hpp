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
int XMPSlicesCalc(int p, int nx, int ProcessorsInXDirection, int ProcessorsInYDirection, int DecompositionStrategy);
int XOffsetCalc(int p, int nx, int ProcessorsInXDirection, int ProcessorsInYDirection, int DecompositionStrategy);
int YMPSlicesCalc(int p, int ny, int ProcessorsInYDirection, int np, int DecompositionStrategy);
int YOffsetCalc(int p, int ny, int ProcessorsInYDirection, int np, int DecompositionStrategy);
void AddGhostNodes(int DecompositionStrategy, int NeighborRank_West, int NeighborRank_East, int NeighborRank_North,
                   int NeighborRank_South, int &XRemoteMPSlices, int &RemoteXOffset, int &YRemoteMPSlices,
                   int &RemoteYOffset);
void InitialDecomposition(int &DecompositionStrategy, int nx, int ny, int &ProcessorsInXDirection,
                          int &ProcessorsInYDirection, int id, int np, int &NeighborRank_North, int &NeighborRank_South,
                          int &NeighborRank_East, int &NeighborRank_West, int &NeighborRank_NorthEast,
                          int &NeighborRank_NorthWest, int &NeighborRank_SouthEast, int &NeighborRank_SouthWest,
                          bool &AtNorthBoundary, bool &AtSouthBoundary, bool &AtEastBoundary, bool &AtWestBoundary);

#endif
