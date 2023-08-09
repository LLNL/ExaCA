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
// Subdivide ny into ny_local across np ranks as evenly as possible, return ny_local on rank p
int get_nylocal(int p, int ny, int np) {
    int ny_local = 0;
    int ny_local_est = ny / np;
    int YRemainder = ny % np;
    if (YRemainder == 0) {
        ny_local = ny_local_est;
    }
    else {
        if (YRemainder > p) {
            ny_local = ny_local_est + 1;
        }
        else {
            ny_local = ny_local_est;
        }
    }
    return ny_local;
}

//*****************************************************************************/
// Get the offset in the Y direction of the local grid on rank p from the global grid of np ranks, which has ny total
// cells in Y
int get_yoffset(int p, int ny, int np) {
    int y_offset = 0;
    int y_offset_est = ny / np;
    int YRemainder = ny % np;
    if (YRemainder == 0) {
        y_offset = p * y_offset_est;
    }
    else {
        if (YRemainder > p) {
            y_offset = p * (y_offset_est + 1);
        }
        else {
            y_offset = (y_offset_est + 1) * (YRemainder - 1) + y_offset_est + 1 + (p - YRemainder) * y_offset_est;
        }
    }
    return y_offset;
}

// Add ghost nodes to the appropriate subdomains (added where the subdomains overlap, but not at edges of physical
// domain)
void AddGhostNodes(int NeighborRank_North, int NeighborRank_South, int &ny_local, int &y_offset) {

    // Add halo regions in Y direction if this subdomain borders subdomains on other processors
    // If only 1 rank in the y direction, no halo regions - subdomain is coincident with overall simulation domain
    // If multiple ranks in the y direction, either 1 halo region (borders another rank's subdomain in either the +y or
    // -y direction) or 2 halo regions (if it borders other rank's subdomains in both the +y and -y directions)
    if (NeighborRank_North != MPI_PROC_NULL)
        ny_local++;
    if (NeighborRank_South != MPI_PROC_NULL) {
        ny_local++;
        // Also adjust subdomain offset, as these ghost nodes were added on the -y side of the subdomain
        y_offset--;
    }
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

// Stores/returns the volume fraction of nucleated grains to the console
float calcVolFractionNucleated(int id, int nx, int ny_local, int DomainSize, ViewS LayerID, ViewI GrainID,
                               bool AtNorthBoundary, bool AtSouthBoundary) {

    // For interior cells, add the number of cells that underwent melting/solidification and the number of cells with
    // sub-zero grain IDs
    int MeltedCells_Local = 0;
    int NucleatedGrainCells_Local = 0;
    Kokkos::parallel_reduce(
        "NumSolidifiedCells", DomainSize,
        KOKKOS_LAMBDA(const int index, int &update_meltcount, int &update_nucleatecount) {
            int coord_y = getCoordY(index, nx, ny_local);
            // Is this Y coordinate in the halo region? If so, do not increment counter
            bool InHaloRegion = false;
            if (((coord_y == 0) && (!AtSouthBoundary)) || ((coord_y == ny_local - 1) && (!AtNorthBoundary)))
                InHaloRegion = true;
            if ((GrainID(index) < 0) && (!InHaloRegion))
                update_nucleatecount++;
            if ((LayerID(index) != -1) && (!InHaloRegion))
                update_meltcount++;
        },
        MeltedCells_Local, NucleatedGrainCells_Local);

    // Reduce the values by summing over all ranks
    int MeltedCells_Global, NucleatedGrainCells_Global;
    MPI_Allreduce(&MeltedCells_Local, &MeltedCells_Global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&NucleatedGrainCells_Local, &NucleatedGrainCells_Global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // Calculate nucleated grain fraction
    float VolFractionNucleated =
        static_cast<float>(NucleatedGrainCells_Global) / static_cast<float>(MeltedCells_Global);
    if (id == 0)
        std::cout << "The fraction of the solidified material consisting of nucleated grains is "
                  << VolFractionNucleated << std::endl;
    return VolFractionNucleated;
}
