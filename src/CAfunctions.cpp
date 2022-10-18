// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAfunctions.hpp"

#include "mpi.h"

#include <cmath>
#include <iostream>

/*************************** FUNCTIONS CALLED THROUGH MAIN SUBROUTINES ***************************/

//*****************************************************************************/
int XMPSlicesCalc(int p, int nx, int ProcessorsInXDirection, int ProcessorsInYDirection, int DecompositionStrategy) {
    int XRemoteMPSlices = 0;
    if (DecompositionStrategy == 1) {
        XRemoteMPSlices = nx; // no ghost nodes, domain not divided up in the x direction
    }
    else if (DecompositionStrategy >= 2) {
        int XSlicesPerP = nx / ProcessorsInXDirection;
        int XRemainder = nx % ProcessorsInXDirection;
        if (XRemainder == 0) {
            XRemoteMPSlices = XSlicesPerP;
        }
        else {
            int XPosition = p / ProcessorsInYDirection;
            if (XPosition < XRemainder) {
                XRemoteMPSlices = XSlicesPerP + 1;
            }
            else {
                XRemoteMPSlices = XSlicesPerP;
            }
        }
    }
    return XRemoteMPSlices;
}

//*****************************************************************************/
int XOffsetCalc(int p, int nx, int ProcessorsInXDirection, int ProcessorsInYDirection, int DecompositionStrategy) {
    int RemoteXOffset = 0;
    if (DecompositionStrategy == 1) {
        // No decomposition in the x direction, so no ghost nodes to affect the local domain offset
        RemoteXOffset = 0;
    }
    else if (DecompositionStrategy >= 2) {
        int XSlicesPerP = nx / ProcessorsInXDirection;
        int XRemainder = nx % ProcessorsInXDirection;
        int XPosition = p / ProcessorsInYDirection;
        if (XRemainder == 0) {
            RemoteXOffset = XPosition * XSlicesPerP;
        }
        else {
            if (XRemainder > XPosition) {
                RemoteXOffset = XPosition * (XSlicesPerP + 1);
            }
            else {
                RemoteXOffset =
                    (XSlicesPerP + 1) * (XRemainder - 1) + XSlicesPerP + 1 + (XPosition - XRemainder) * XSlicesPerP;
            }
        }
    }

    return RemoteXOffset;
}

//*****************************************************************************/
int YMPSlicesCalc(int p, int ny, int ProcessorsInYDirection, int np, int DecompositionStrategy) {
    int YRemoteMPSlices = 0;
    if (DecompositionStrategy == 1) {
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
    }
    else if (DecompositionStrategy >= 2) {
        int YSlicesPerP = ny / ProcessorsInYDirection;
        int YRemainder = ny % ProcessorsInYDirection;
        if (YRemainder == 0) {
            YRemoteMPSlices = YSlicesPerP;
        }
        else {
            int YPosition = p % ProcessorsInYDirection;
            if (YPosition < YRemainder) {
                YRemoteMPSlices = YSlicesPerP + 1;
            }
            else {
                YRemoteMPSlices = YSlicesPerP;
            }
        }
    }

    return YRemoteMPSlices;
}

//*****************************************************************************/
int YOffsetCalc(int p, int ny, int ProcessorsInYDirection, int np, int DecompositionStrategy) {
    int RemoteYOffset = 0;
    if (DecompositionStrategy == 1) {
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
    }
    else if (DecompositionStrategy >= 2) {
        int YSlicesPerP = ny / ProcessorsInYDirection;
        int YRemainder = ny % ProcessorsInYDirection;
        int YPosition = p % ProcessorsInYDirection;
        if (YRemainder == 0) {
            RemoteYOffset = YPosition * YSlicesPerP;
        }
        else {
            if (YRemainder > YPosition) {
                RemoteYOffset = YPosition * (YSlicesPerP + 1);
            }
            else {
                RemoteYOffset =
                    (YSlicesPerP + 1) * (YRemainder - 1) + YSlicesPerP + 1 + (YPosition - YRemainder) * YSlicesPerP;
            }
        }
    }
    return RemoteYOffset;
}

// Add ghost nodes to the appropriate subdomains (added where the subdomains overlap, but not at edges of physical
// domain)
void AddGhostNodes(int DecompositionStrategy, int NeighborRank_West, int NeighborRank_East, int NeighborRank_North,
                   int NeighborRank_South, int &XRemoteMPSlices, int &RemoteXOffset, int &YRemoteMPSlices,
                   int &RemoteYOffset) {

    if (DecompositionStrategy > 1) {
        // Decomposing domain in X and Y directions
        // Add halo regions in X if this subdomain borders subdomains on other processors
        // If only 1 rank in the x direction, no halo regions - subdomain is coincident with overall simulation domain
        // If multiple ranks in the x direction, either 1 halo region (borders another rank's subdomain in either the +x
        // or -x direction) or 2 halo regions (if it borders other rank's subdomains in both the +x and -x directions)
        if (NeighborRank_West != MPI_PROC_NULL) {
            XRemoteMPSlices++;
            // Also adjust subdomain offset, as these ghost nodes were added on the -x side of the subdomain
            RemoteXOffset--;
        }
        if (NeighborRank_East != MPI_PROC_NULL)
            XRemoteMPSlices++;
    }
    // Add halo regions in Y direction if this subdomain borders subdomains on other processors
    // If only 1 rank in the y direction, no halo regions - subdomain is coincident with overall simulation domain
    // If multiple ranks in the y direction, either 1 halo region (borders another rank's subdomain in either the +y or
    // -y direction) or 2 halo regions (if it borders other rank's subdomains in both the +y and -y directions)
    if (NeighborRank_North != MPI_PROC_NULL)
        YRemoteMPSlices++;
    if (NeighborRank_South != MPI_PROC_NULL) {
        YRemoteMPSlices++;
        // Also adjust subdomain offset, as these ghost nodes were added on the -y side of the subdomain
        RemoteYOffset--;
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
void InitialDecomposition(int &DecompositionStrategy, int nx, int ny, int &ProcessorsInXDirection,
                          int &ProcessorsInYDirection, int id, int np, int &NeighborRank_North, int &NeighborRank_South,
                          int &NeighborRank_East, int &NeighborRank_West, int &NeighborRank_NorthEast,
                          int &NeighborRank_NorthWest, int &NeighborRank_SouthEast, int &NeighborRank_SouthWest,
                          bool &AtNorthBoundary, bool &AtSouthBoundary, bool &AtEastBoundary, bool &AtWestBoundary) {

    if (DecompositionStrategy == 1) {
    OneDim:
        ProcessorsInXDirection = 1;
        ProcessorsInYDirection = np;
        if (np > 1) {
            NeighborRank_South = id - 1;
            NeighborRank_North = id + 1;
            if (id == 0)
                NeighborRank_South = MPI_PROC_NULL;
            if (id == np - 1)
                NeighborRank_North = MPI_PROC_NULL;
            NeighborRank_West = MPI_PROC_NULL;
            NeighborRank_NorthWest = MPI_PROC_NULL;
            NeighborRank_NorthEast = MPI_PROC_NULL;
            NeighborRank_East = MPI_PROC_NULL;
            NeighborRank_SouthWest = MPI_PROC_NULL;
            NeighborRank_SouthEast = MPI_PROC_NULL;
        }
        else {
            // No MPI communication
            NeighborRank_West = MPI_PROC_NULL;
            NeighborRank_East = MPI_PROC_NULL;
            NeighborRank_North = MPI_PROC_NULL;
            NeighborRank_NorthWest = MPI_PROC_NULL;
            NeighborRank_NorthEast = MPI_PROC_NULL;
            NeighborRank_South = MPI_PROC_NULL;
            NeighborRank_SouthWest = MPI_PROC_NULL;
            NeighborRank_SouthEast = MPI_PROC_NULL;
        }
    }
    else if (DecompositionStrategy == 2) {
        if (np % 2 != 0) {
            if (id == 0)
                std::cout << "Warning: Number of processors not divisible by two- defaulting to 1D decomposition"
                          << std::endl;
            DecompositionStrategy = 1;
            goto OneDim;
        }
        ProcessorsInXDirection = 2;
        ProcessorsInYDirection = np / 2;
        NeighborRank_South = id - 1;
        NeighborRank_North = id + 1;
        NeighborRank_East = id + ProcessorsInYDirection;
        NeighborRank_West = id - ProcessorsInYDirection;
        NeighborRank_SouthEast = NeighborRank_East - 1;
        NeighborRank_NorthEast = NeighborRank_East + 1;
        NeighborRank_SouthWest = NeighborRank_West - 1;
        NeighborRank_NorthWest = NeighborRank_West + 1;
        if (id % ProcessorsInYDirection == 0)
            NeighborRank_South = MPI_PROC_NULL;
        if (id % ProcessorsInYDirection == ProcessorsInYDirection - 1)
            NeighborRank_North = MPI_PROC_NULL;
        if (id < ProcessorsInYDirection)
            NeighborRank_West = MPI_PROC_NULL;
        if (id >= (np - ProcessorsInYDirection))
            NeighborRank_East = MPI_PROC_NULL;
        if ((NeighborRank_South == MPI_PROC_NULL) || (NeighborRank_East == MPI_PROC_NULL))
            NeighborRank_SouthEast = MPI_PROC_NULL;
        if ((NeighborRank_North == MPI_PROC_NULL) || (NeighborRank_East == MPI_PROC_NULL))
            NeighborRank_NorthEast = MPI_PROC_NULL;
        if ((NeighborRank_South == MPI_PROC_NULL) || (NeighborRank_West == MPI_PROC_NULL))
            NeighborRank_SouthWest = MPI_PROC_NULL;
        if ((NeighborRank_North == MPI_PROC_NULL) || (NeighborRank_West == MPI_PROC_NULL))
            NeighborRank_NorthWest = MPI_PROC_NULL;
    }
    else if (DecompositionStrategy == 3) {

        // Is the number of ranks a perfect square or divisible by 2?
        bool DivisibleBy2, PerfectSquare;
        if (np % 2 != 0) {
            DivisibleBy2 = false;
        }
        else {
            DivisibleBy2 = true;
        }
        double Square = sqrt(np);
        if (Square == round(Square)) {
            PerfectSquare = true;
        }
        else {
            PerfectSquare = false;
        }

        if (((!PerfectSquare) && (!DivisibleBy2)) || (np == 2)) {
            if (id == 0)
                std::cout
                    << "Note: Number of processors is either 2, not divisible by two, or not divisible by itself- "
                       "defaulting to 1D decomposition"
                    << std::endl;
            DecompositionStrategy = 1;
            goto OneDim;
        }
        int PX, PY;

        if (PerfectSquare) {
            PX = Square;
            PY = Square;
        }
        else {
            int RoundSquare = floor(Square);
            int NotSquare = 1;
            while (NotSquare == 1) {
                if (np % RoundSquare == 0) {
                    PY = np / RoundSquare;
                    PX = RoundSquare;
                    NotSquare = 0;
                }
                else {
                    RoundSquare--;
                }
            }
        }
        if (((nx > ny) && (PY > PX)) || ((ny > nx) && (PX > PY))) {
            ProcessorsInXDirection = PY;
            ProcessorsInYDirection = PX;
        }
        else {
            ProcessorsInXDirection = PX;
            ProcessorsInYDirection = PY;
        }
        if ((ProcessorsInXDirection == 2) && (id == 0))
            std::cout << "Note: Decomposition pattern 3 is equivalent to 2" << std::endl;
        if (id == 0)
            std::cout << " Processors in X: " << ProcessorsInXDirection
                      << " Processors in Y: " << ProcessorsInYDirection << std::endl;
        NeighborRank_South = id - 1;
        NeighborRank_North = id + 1;
        NeighborRank_East = id + ProcessorsInYDirection;
        NeighborRank_West = id - ProcessorsInYDirection;
        if (id % ProcessorsInYDirection == 0)
            NeighborRank_South = MPI_PROC_NULL;
        if (id % ProcessorsInYDirection == ProcessorsInYDirection - 1)
            NeighborRank_North = MPI_PROC_NULL;
        if (id < ProcessorsInYDirection)
            NeighborRank_West = MPI_PROC_NULL;
        if (id >= (np - ProcessorsInYDirection))
            NeighborRank_East = MPI_PROC_NULL;
        if ((NeighborRank_South == MPI_PROC_NULL) || (NeighborRank_East == MPI_PROC_NULL))
            NeighborRank_SouthEast = MPI_PROC_NULL;
        else
            NeighborRank_SouthEast = NeighborRank_East - 1;

        if ((NeighborRank_North == MPI_PROC_NULL) || (NeighborRank_East == MPI_PROC_NULL))
            NeighborRank_NorthEast = MPI_PROC_NULL;
        else
            NeighborRank_NorthEast = NeighborRank_East + 1;

        if ((NeighborRank_South == MPI_PROC_NULL) || (NeighborRank_West == MPI_PROC_NULL))
            NeighborRank_SouthWest = MPI_PROC_NULL;
        else
            NeighborRank_SouthWest = NeighborRank_West - 1;

        if ((NeighborRank_North == MPI_PROC_NULL) || (NeighborRank_West == MPI_PROC_NULL))
            NeighborRank_NorthWest = MPI_PROC_NULL;
        else
            NeighborRank_NorthWest = NeighborRank_West + 1;
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
    if (NeighborRank_West == MPI_PROC_NULL)
        AtWestBoundary = true;
    else
        AtWestBoundary = false;
    if (NeighborRank_East == MPI_PROC_NULL)
        AtEastBoundary = true;
    else
        AtEastBoundary = false;
}

//*****************************************************************************/
// Determine the scan limits for each rank in X and Y directions
void XYLimitCalc(int &LLX, int &LLY, int &ULX, int &ULY, int MyXSlices, int MyYSlices, int NeighborRank_South,
                 int NeighborRank_North, int NeighborRank_East, int NeighborRank_West) {
    LLX = 0;
    LLY = 0;
    ULX = MyXSlices - 1;
    ULY = MyYSlices - 1;
    if (NeighborRank_South == MPI_PROC_NULL)
        LLY = 1;
    if (NeighborRank_North == MPI_PROC_NULL)
        ULY = MyYSlices - 2;
    if (NeighborRank_East == MPI_PROC_NULL)
        ULX = MyXSlices - 2;
    if (NeighborRank_West == MPI_PROC_NULL)
        LLX = 1;
}
