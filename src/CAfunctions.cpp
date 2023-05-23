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
int YMPSlicesCalc(int p, int ny, int np) {
    int YRemoteMPSlices = 0;
    int YSlicesPerP = ny / np;
    int YRemainder = ny % np;
    if (YRemainder == 0) {
        YRemoteMPSlices = YSlicesPerP;
    }
    else {
        if (YRemainder > p) {
            YRemoteMPSlices = YSlicesPerP + 1;
        }
        else {
            YRemoteMPSlices = YSlicesPerP;
        }
    }
    return YRemoteMPSlices;
}

//*****************************************************************************/
int YOffsetCalc(int p, int ny, int np) {
    int RemoteYOffset = 0;
    int YSlicesPerP = ny / np;
    int YRemainder = ny % np;
    if (YRemainder == 0) {
        RemoteYOffset = p * YSlicesPerP;
    }
    else {
        if (YRemainder > p) {
            RemoteYOffset = p * (YSlicesPerP + 1);
        }
        else {
            RemoteYOffset = (YSlicesPerP + 1) * (YRemainder - 1) + YSlicesPerP + 1 + (p - YRemainder) * YSlicesPerP;
        }
    }
    return RemoteYOffset;
}

// Add ghost nodes to the appropriate subdomains (added where the subdomains overlap, but not at edges of physical
// domain)
void AddGhostNodes(int NeighborRank_North, int NeighborRank_South, int &MyYSlices, int &MyYOffset) {

    // Add halo regions in Y direction if this subdomain borders subdomains on other processors
    // If only 1 rank in the y direction, no halo regions - subdomain is coincident with overall simulation domain
    // If multiple ranks in the y direction, either 1 halo region (borders another rank's subdomain in either the +y or
    // -y direction) or 2 halo regions (if it borders other rank's subdomains in both the +y and -y directions)
    if (NeighborRank_North != MPI_PROC_NULL)
        MyYSlices++;
    if (NeighborRank_South != MPI_PROC_NULL) {
        MyYSlices++;
        // Also adjust subdomain offset, as these ghost nodes were added on the -y side of the subdomain
        MyYOffset--;
    }
}

//*****************************************************************************/
int FindItBounds(int RankX, int RankY, int MyXSlices, int MyYSlices) {
    int ItBounds;
    // If X and Y coordinates are not on edges, Case 0: iteratation over neighbors 0-25 possible
    // If Y coordinate is on lower edge, Case 1: iteration over only neighbors 9-25 possible
    // If Y coordinate is on upper edge, Case 2: iteration over only neighbors 0-16 possible
    // If X coordinate is on lower edge, Case 3: iteration over only neighbors
    // 0,1,3,4,6,8,9,10,11,13,15,17,18,20,21,22,24 If X coordinate is on upper edge, Case 4: iteration over only
    // neighbors 0,2,3,4,5,7,9,10,12,14,16,17,19,20,21,23,25 If X/Y coordinates are on lower edge, Case 5: iteration
    // over only neighbors 9,10,11,13,15,17,18,20,21,22,24 If X coordinate is on upper edge/Y on lower edge, Case 6: If
    // X coordinate is on lower edge/Y on upper edge, Case 7: If X/Y coordinates are on upper edge, Case 8:
    if (RankY == 0) {
        if (RankX == 0) {
            ItBounds = 5;
        }
        else if (RankX == MyXSlices - 1) {
            ItBounds = 6;
        }
        else {
            ItBounds = 1;
        }
    }
    else if (RankY == MyYSlices - 1) {
        if (RankX == 0) {
            ItBounds = 7;
        }
        else if (RankX == MyXSlices - 1) {
            ItBounds = 8;
        }
        else {
            ItBounds = 2;
        }
    }
    else {
        if (RankX == 0) {
            ItBounds = 3;
        }
        else if (RankX == MyXSlices - 1) {
            ItBounds = 4;
        }
        else {
            ItBounds = 0;
        }
    }
    return ItBounds;
}

//*****************************************************************************/
// Determine the mapping of processors to grid data
void InitialDecomposition(int id, int np, int &NeighborRank_North, int &NeighborRank_South, bool &AtNorthBoundary,
                          bool &AtSouthBoundary) {

    if (np > 1) {
        NeighborRank_South = id - 1;
        NeighborRank_North = id + 1;
        if (id == 0)
            NeighborRank_South = MPI_PROC_NULL;
        if (id == np - 1)
            NeighborRank_North = MPI_PROC_NULL;
    }
    else {
        // No MPI communication
        NeighborRank_North = MPI_PROC_NULL;
        NeighborRank_South = MPI_PROC_NULL;
    }

    // Based on the decomposition, store whether each MPI rank is each boundary or not
    if (NeighborRank_North == MPI_PROC_NULL)
        AtNorthBoundary = true;
    else
        AtNorthBoundary = false;
    if (NeighborRank_South == MPI_PROC_NULL)
        AtSouthBoundary = true;
    else
        AtSouthBoundary = false;
}

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
