// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAupdate.hpp"
#include "CAfunctions.hpp"
#include "CAghostnodes.hpp"
#include "CAprint.hpp"
#include "mpi.h"

#include <cmath>

// Using for compatibility with device math functions.
using std::max;
using std::min;

//*****************************************************************************/
void Nucleation(int cycle, int &SuccessfulNucEvents_ThisRank, int &NucleationCounter, int PossibleNuclei_ThisRank,
                ViewI_H NucleationTimes_H, ViewI NucleiLocations, ViewI NucleiGrainID, ViewI CellType, ViewI GrainID,
                int ZBound_Low, int MyXSlices, int MyYSlices, ViewI SteeringVector, ViewI numSteer_G) {

    // Is there nucleation left in this layer to check?
    if (NucleationCounter < PossibleNuclei_ThisRank) {
        // Is there at least one potential nucleation event on this rank, at this time step?
        if (cycle == NucleationTimes_H(NucleationCounter)) {
            bool NucleationCheck = true;
            int FirstEvent = NucleationCounter; // first potential nucleation event to check
            // Are there any other nucleation events this time step to check?
            while (NucleationCheck) {
                NucleationCounter++;
                // If the previous nucleation event was the last one for this layer of the simulation, exit loop
                if (NucleationCounter == PossibleNuclei_ThisRank)
                    break;
                // If the next nucleation event corresponds to a future time step, finish check
                if (cycle != NucleationTimes_H(NucleationCounter))
                    NucleationCheck = false;
            }
            int LastEvent = NucleationCounter;
            // parallel_reduce checks each potential nucleation event this time step (FirstEvent, up to but not
            // including LastEvent)
            int NucleationThisDT = 0; // return number of successful event from parallel_reduce
            // Launch kokkos kernel - check if the corresponding CA cell location is liquid
            Kokkos::parallel_reduce(
                "NucleiUpdateLoop", Kokkos::RangePolicy<>(FirstEvent, LastEvent),
                KOKKOS_LAMBDA(const int NucleationCounter_Device, int &update) {
                    int NucleationEventLocation_GlobalGrid = NucleiLocations(NucleationCounter_Device);
                    int update_val =
                        FutureActive; // added to steering vector to become a new active cell as part of cellcapture
                    int old_val = Liquid;
                    int OldCellTypeValue = Kokkos::atomic_compare_exchange(
                        &CellType(NucleationEventLocation_GlobalGrid), old_val, update_val);
                    if (OldCellTypeValue == Liquid) {
                        // Successful nucleation event - atomic update of cell type, proceeded if the atomic
                        // exchange is successful (cell was liquid) Add future active cell location to steering
                        // vector and change cell type, assign new Grain ID
                        GrainID(NucleationEventLocation_GlobalGrid) = NucleiGrainID(NucleationCounter_Device);
                        int GlobalZ = NucleationEventLocation_GlobalGrid / (MyXSlices * MyYSlices);
                        int Rem = NucleationEventLocation_GlobalGrid % (MyXSlices * MyYSlices);
                        int RankX = Rem / MyYSlices;
                        int RankY = Rem % MyYSlices;
                        int RankZ = GlobalZ - ZBound_Low;
                        int NucleationEventLocation_LocalGrid =
                            RankZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                        SteeringVector(Kokkos::atomic_fetch_add(&numSteer_G(0), 1)) = NucleationEventLocation_LocalGrid;
                        // This undercooled liquid cell is now a nuclei (no nuclei are in the ghost nodes - halo
                        // exchange routine GhostNodes1D or GhostNodes2D is used to fill these)
                        update++;
                    }
                },
                NucleationThisDT);
            // Update the number of successful nuclei counter with the number of successful nucleation events from this
            // time step (NucleationThisDT)
            SuccessfulNucEvents_ThisRank += NucleationThisDT;
        }
    }
}

//*****************************************************************************/
// Determine which cells are associated with the "steering vector" of cells that are either active, or becoming active
// this time step
void FillSteeringVector_NoRemelt(int cycle, int LocalActiveDomainSize, int MyXSlices, int MyYSlices, ViewI CritTimeStep,
                                 ViewF UndercoolingCurrent, ViewF UndercoolingChange, ViewI CellType, int ZBound_Low,
                                 int layernumber, ViewI LayerID, ViewI SteeringVector, ViewI numSteer,
                                 ViewI_H numSteer_Host) {

    // Cells associated with this layer that are not solid type but have passed the liquidus (crit time step) have their
    // undercooling values updated Cells that meet the aforementioned criteria and are active type should be added to
    // the steering vector
    Kokkos::parallel_for(
        "FillSV", LocalActiveDomainSize, KOKKOS_LAMBDA(const int &D3D1ConvPosition) {
            // Cells of interest for the CA
            int RankZ = D3D1ConvPosition / (MyXSlices * MyYSlices);
            int Rem = D3D1ConvPosition % (MyXSlices * MyYSlices);
            int RankX = Rem / MyYSlices;
            int RankY = Rem % MyYSlices;
            int GlobalZ = RankZ + ZBound_Low;
            int GlobalD3D1ConvPosition = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
            int cellType = CellType(GlobalD3D1ConvPosition);

            int layerCheck = (LayerID(GlobalD3D1ConvPosition) <= layernumber);
            int isNotSolid = (cellType != Solid);
            int pastCritTime = (cycle > CritTimeStep(GlobalD3D1ConvPosition));

            int cell_Liquid = (cellType == Liquid);
            int cell_Active = (cellType == Active);

            if (layerCheck && isNotSolid && pastCritTime) {
                UndercoolingCurrent(GlobalD3D1ConvPosition) +=
                    UndercoolingChange(GlobalD3D1ConvPosition) * (cell_Liquid + cell_Active);
                if (cell_Active) {
                    SteeringVector(Kokkos::atomic_fetch_add(&numSteer(0), 1)) = D3D1ConvPosition;
                }
            }
        });
    Kokkos::deep_copy(numSteer_Host, numSteer);
}

//*****************************************************************************/
// Determine which cells are associated with the "steering vector" of cells that are either active, or becoming active
// this time step - version with remelting
void FillSteeringVector_Remelt(int cycle, int LocalActiveDomainSize, int MyXSlices, int MyYSlices, NList NeighborX,
                               NList NeighborY, NList NeighborZ, ViewI CritTimeStep, ViewF UndercoolingCurrent,
                               ViewF UndercoolingChange, ViewI CellType, ViewI GrainID, int ZBound_Low, int nzActive,
                               ViewI SteeringVector, ViewI numSteer, ViewI_H numSteer_Host, ViewI MeltTimeStep,
                               int BufSizeX, int BufSizeY, bool AtNorthBoundary, bool AtSouthBoundary,
                               bool AtEastBoundary, bool AtWestBoundary, Buffer2D BufferWestSend,
                               Buffer2D BufferEastSend, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend,
                               Buffer2D BufferNorthEastSend, Buffer2D BufferNorthWestSend, Buffer2D BufferSouthEastSend,
                               Buffer2D BufferSouthWestSend, int DecompositionStrategy) {

    Kokkos::parallel_for(
        "FillSV_RM", LocalActiveDomainSize, KOKKOS_LAMBDA(const int &D3D1ConvPosition) {
            // Coordinate of this cell on the "global" (all cells in the Z direction) grid
            int RankZ = D3D1ConvPosition / (MyXSlices * MyYSlices);
            int Rem = D3D1ConvPosition % (MyXSlices * MyYSlices);
            int RankX = Rem / MyYSlices;
            int RankY = Rem % MyYSlices;
            int GlobalZ = RankZ + ZBound_Low;
            int GlobalD3D1ConvPosition = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;

            int cellType = CellType(GlobalD3D1ConvPosition);
            bool isNotSolid = ((cellType != TempSolid) && (cellType != Solid));
            bool atMeltTime = (cycle == MeltTimeStep(GlobalD3D1ConvPosition));
            bool pastCritTime = (cycle > CritTimeStep(GlobalD3D1ConvPosition));
            if ((atMeltTime) && ((cellType == TempSolid) || (cellType == Active))) {
                // This cell should be a liquid cell
                CellType(GlobalD3D1ConvPosition) = Liquid;
                // Reset current undercooling to zero
                UndercoolingCurrent(GlobalD3D1ConvPosition) = 0.0;
                // Remove solid cell data from the buffer
                if (DecompositionStrategy == 1)
                    loadghostnodes(0, 0, 0, 0, 0, BufSizeX, MyYSlices, RankX, RankY, RankZ, AtNorthBoundary,
                                   AtSouthBoundary, BufferSouthSend, BufferNorthSend);
                else
                    loadghostnodes(0, 0, 0, 0, 0, BufSizeX, BufSizeY, MyXSlices, MyYSlices, RankX, RankY, RankZ,
                                   AtNorthBoundary, AtSouthBoundary, AtWestBoundary, AtEastBoundary, BufferSouthSend,
                                   BufferNorthSend, BufferWestSend, BufferEastSend, BufferNorthEastSend,
                                   BufferSouthEastSend, BufferSouthWestSend, BufferNorthWestSend);
            }
            else if ((isNotSolid) && (pastCritTime)) {
                // Update cell undercooling
                UndercoolingCurrent(GlobalD3D1ConvPosition) += UndercoolingChange(GlobalD3D1ConvPosition);
                if (cellType == Active) {
                    // Add active cells below liquidus to steering vector
                    SteeringVector(Kokkos::atomic_fetch_add(&numSteer(0), 1)) = D3D1ConvPosition;
                }
                else if ((cellType == Liquid) && (GrainID(GlobalD3D1ConvPosition) != 0)) {
                    // If this cell borders at least one solid/tempsolid cell and is part of a grain, it should become
                    // active
                    for (int l = 0; l < 26; l++) {
                        // "l" correpsponds to the specific neighboring cell
                        // Local coordinates of adjacent cell center
                        int MyNeighborX = RankX + NeighborX[l];
                        int MyNeighborY = RankY + NeighborY[l];
                        int MyNeighborZ = RankZ + NeighborZ[l];
                        if ((MyNeighborX >= 0) && (MyNeighborX < MyXSlices) && (MyNeighborY >= 0) &&
                            (MyNeighborY < MyYSlices) && (MyNeighborZ < nzActive) && (MyNeighborZ >= 0)) {
                            int GlobalNeighborD3D1ConvPosition = (MyNeighborZ + ZBound_Low) * MyXSlices * MyYSlices +
                                                                 MyNeighborX * MyYSlices + MyNeighborY;
                            if ((CellType(GlobalNeighborD3D1ConvPosition) == TempSolid) ||
                                (CellType(GlobalNeighborD3D1ConvPosition) == Solid) || (RankZ == 0)) {
                                // Cell activation to be performed as part of steering vector
                                l = 26;
                                SteeringVector(Kokkos::atomic_fetch_add(&numSteer(0), 1)) = D3D1ConvPosition;
                                CellType(GlobalD3D1ConvPosition) =
                                    FutureActive; // this cell cannot be captured - is being activated
                            }
                        }
                    }
                }
            }
        });
    Kokkos::fence();

    // Copy size of steering vector (containing positions of undercooled liquid/active cells) to the host
    Kokkos::deep_copy(numSteer_Host, numSteer);
}

// Decentered octahedron algorithm for the capture of new interface cells by grains
void CellCapture(int, int np, int, int DecompositionStrategy, int, int, int MyXSlices, int MyYSlices,
                 InterfacialResponseFunction irf, int MyXOffset, int MyYOffset, NList NeighborX, NList NeighborY,
                 NList NeighborZ, ViewI CritTimeStep, ViewF UndercoolingCurrent, ViewF UndercoolingChange,
                 ViewF GrainUnitVector, ViewF CritDiagonalLength, ViewF DiagonalLength, ViewI CellType, ViewF DOCenter,
                 ViewI GrainID, int NGrainOrientations, Buffer2D BufferWestSend, Buffer2D BufferEastSend,
                 Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, Buffer2D BufferNorthEastSend,
                 Buffer2D BufferNorthWestSend, Buffer2D BufferSouthEastSend, Buffer2D BufferSouthWestSend, int BufSizeX,
                 int BufSizeY, int ZBound_Low, int nzActive, int, ViewI SteeringVector, ViewI numSteer,
                 ViewI_H numSteer_Host, bool AtNorthBoundary, bool AtSouthBoundary, bool AtEastBoundary,
                 bool AtWestBoundary, ViewI SolidificationEventCounter, ViewI MeltTimeStep,
                 ViewF3D LayerTimeTempHistory, ViewI NumberOfSolidificationEvents, bool RemeltingYN) {

    // Loop over list of active and soon-to-be active cells, potentially performing cell capture events and updating
    // cell types
    Kokkos::parallel_for(
        "CellCapture", numSteer_Host(0), KOKKOS_LAMBDA(const int &num) {
            numSteer(0) = 0;
            int D3D1ConvPosition = SteeringVector(num);
            // Cells of interest for the CA - active cells and future active cells
            int RankZ = D3D1ConvPosition / (MyXSlices * MyYSlices);
            int Rem = D3D1ConvPosition % (MyXSlices * MyYSlices);
            int RankX = Rem / MyYSlices;
            int RankY = Rem % MyYSlices;
            int GlobalZ = RankZ + ZBound_Low;
            int GlobalD3D1ConvPosition = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
            if (CellType(GlobalD3D1ConvPosition) == Active) {
                // Update local diagonal length of active cell
                double LocU = UndercoolingCurrent(GlobalD3D1ConvPosition);
                LocU = min(210.0, LocU);
                double V = irf.compute(LocU);
                DiagonalLength(D3D1ConvPosition) += min(0.045, V); // Max amount the diagonal can grow per time step
                // Cycle through all neigboring cells on this processor to see if they have been captured
                // Cells in ghost nodes cannot capture cells on other processors
                bool DeactivateCell = true; // switch that becomes false if the cell has at least 1 liquid type neighbor
                // Which neighbors should be iterated over?
                for (int l = 0; l < 26; l++) {
                    // Local coordinates of adjacent cell center
                    int MyNeighborX = RankX + NeighborX[l];
                    int MyNeighborY = RankY + NeighborY[l];
                    int MyNeighborZ = RankZ + NeighborZ[l];
                    // Check if neighbor is in bounds
                    if ((MyNeighborX >= 0) && (MyNeighborX < MyXSlices) && (MyNeighborY >= 0) &&
                        (MyNeighborY < MyYSlices) && (MyNeighborZ < nzActive) && (MyNeighborZ >= 0)) {
                        long int NeighborD3D1ConvPosition =
                            MyNeighborZ * MyXSlices * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
                        long int GlobalNeighborD3D1ConvPosition =
                            (MyNeighborZ + ZBound_Low) * MyXSlices * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
                        if (CellType(GlobalNeighborD3D1ConvPosition) == Liquid)
                            DeactivateCell = false;
                        // Capture of cell located at "NeighborD3D1ConvPosition" if this condition is satisfied
                        if ((DiagonalLength(D3D1ConvPosition) >= CritDiagonalLength(26 * D3D1ConvPosition + l)) &&
                            (CellType(GlobalNeighborD3D1ConvPosition) == Liquid)) {
                            // Use of atomic_compare_exchange
                            // (https://github.com/kokkos/kokkos/wiki/Kokkos%3A%3Aatomic_compare_exchange) old_val =
                            // atomic_compare_exchange(ptr_to_value,comparison_value, new_value); Atomicly sets the
                            // value at the address given by ptr_to_value to new_value if the current value at
                            // ptr_to_value is equal to comparison_value Returns the previously stored value at the
                            // address independent on whether the exchange has happened. If this cell's is a liquid
                            // cell, change it to "TemporaryUpdate" type and return a value of "liquid" If this cell has
                            // already been changed to "TemporaryUpdate" type, return a value of "0"
                            int update_val = TemporaryUpdate;
                            int old_val = Liquid;
                            int OldCellTypeValue = Kokkos::atomic_compare_exchange(
                                &CellType(GlobalNeighborD3D1ConvPosition), old_val, update_val);
                            // Only proceed if CellType was previously liquid (this current thread changed the value to
                            // TemporaryUpdate)
                            if (OldCellTypeValue == Liquid) {
                                int GlobalX = RankX + MyXOffset;
                                int GlobalY = RankY + MyYOffset;
                                int h = GrainID(GlobalD3D1ConvPosition);
                                int MyOrientation = getGrainOrientation(h, NGrainOrientations);

                                // The new cell is captured by this cell's growing octahedron (Grain "h")
                                GrainID(GlobalNeighborD3D1ConvPosition) = h;

                                // (cxold, cyold, czold) are the coordiantes of this decentered octahedron
                                float cxold = DOCenter((long int)(3) * D3D1ConvPosition);
                                float cyold = DOCenter((long int)(3) * D3D1ConvPosition + (long int)(1));
                                float czold = DOCenter((long int)(3) * D3D1ConvPosition + (long int)(2));

                                // (xp,yp,zp) are the global coordinates of the new cell's center
                                float xp = GlobalX + NeighborX[l] + 0.5;
                                float yp = GlobalY + NeighborY[l] + 0.5;
                                float zp = GlobalZ + NeighborZ[l] + 0.5;

                                // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the image
                                // of the center of the new cell
                                float x0 = xp - cxold;
                                float y0 = yp - cyold;
                                float z0 = zp - czold;

                                // mag0 is the magnitude of (x0,y0,z0)
                                float mag0 = sqrtf(x0 * x0 + y0 * y0 + z0 * z0);

                                // Calculate unit vectors for the octahedron that intersect the new cell center
                                float Angle1 = (GrainUnitVector(9 * MyOrientation) * x0 +
                                                GrainUnitVector(9 * MyOrientation + 1) * y0 +
                                                GrainUnitVector(9 * MyOrientation + 2) * z0) /
                                               mag0;
                                float Angle2 = (GrainUnitVector(9 * MyOrientation + 3) * x0 +
                                                GrainUnitVector(9 * MyOrientation + 4) * y0 +
                                                GrainUnitVector(9 * MyOrientation + 5) * z0) /
                                               mag0;
                                float Angle3 = (GrainUnitVector(9 * MyOrientation + 6) * x0 +
                                                GrainUnitVector(9 * MyOrientation + 7) * y0 +
                                                GrainUnitVector(9 * MyOrientation + 8) * z0) /
                                               mag0;
                                float Diag1X = GrainUnitVector(9 * MyOrientation) * (2 * (Angle1 < 0) - 1);
                                float Diag1Y = GrainUnitVector(9 * MyOrientation + 1) * (2 * (Angle1 < 0) - 1);
                                float Diag1Z = GrainUnitVector(9 * MyOrientation + 2) * (2 * (Angle1 < 0) - 1);

                                float Diag2X = GrainUnitVector(9 * MyOrientation + 3) * (2 * (Angle2 < 0) - 1);
                                float Diag2Y = GrainUnitVector(9 * MyOrientation + 4) * (2 * (Angle2 < 0) - 1);
                                float Diag2Z = GrainUnitVector(9 * MyOrientation + 5) * (2 * (Angle2 < 0) - 1);

                                float Diag3X = GrainUnitVector(9 * MyOrientation + 6) * (2 * (Angle3 < 0) - 1);
                                float Diag3Y = GrainUnitVector(9 * MyOrientation + 7) * (2 * (Angle3 < 0) - 1);
                                float Diag3Z = GrainUnitVector(9 * MyOrientation + 8) * (2 * (Angle3 < 0) - 1);

                                float U1[3], U2[3];
                                U1[0] = Diag2X - Diag1X;
                                U1[1] = Diag2Y - Diag1Y;
                                U1[2] = Diag2Z - Diag1Z;
                                U2[0] = Diag3X - Diag1X;
                                U2[1] = Diag3Y - Diag1Y;
                                U2[2] = Diag3Z - Diag1Z;
                                float UU[3];
                                UU[0] = U1[1] * U2[2] - U1[2] * U2[1];
                                UU[1] = U1[2] * U2[0] - U1[0] * U2[2];
                                UU[2] = U1[0] * U2[1] - U1[1] * U2[0];
                                float NDem = sqrtf(UU[0] * UU[0] + UU[1] * UU[1] + UU[2] * UU[2]);
                                float Norm[3];
                                Norm[0] = UU[0] / NDem;
                                Norm[1] = UU[1] / NDem;
                                Norm[2] = UU[2] / NDem;
                                // normal to capturing plane
                                double norm[3], TriangleX[3], TriangleY[3], TriangleZ[3], ParaT;
                                norm[0] = Norm[0];
                                norm[1] = Norm[1];
                                norm[2] = Norm[2];
                                ParaT = (norm[0] * x0 + norm[1] * y0 + norm[2] * z0) /
                                        (norm[0] * Diag1X + norm[1] * Diag1Y + norm[2] * Diag1Z);

                                TriangleX[0] = cxold + ParaT * Diag1X;
                                TriangleY[0] = cyold + ParaT * Diag1Y;
                                TriangleZ[0] = czold + ParaT * Diag1Z;

                                TriangleX[1] = cxold + ParaT * Diag2X;
                                TriangleY[1] = cyold + ParaT * Diag2Y;
                                TriangleZ[1] = czold + ParaT * Diag2Z;

                                TriangleX[2] = cxold + ParaT * Diag3X;
                                TriangleY[2] = cyold + ParaT * Diag3Y;
                                TriangleZ[2] = czold + ParaT * Diag3Z;

                                // Determine which of the 3 corners of the capturing face is closest to the captured
                                // cell center
                                float DistToCorner[3];
                                DistToCorner[0] = sqrtf(((TriangleX[0] - xp) * (TriangleX[0] - xp)) +
                                                        ((TriangleY[0] - yp) * (TriangleY[0] - yp)) +
                                                        ((TriangleZ[0] - zp) * (TriangleZ[0] - zp)));
                                DistToCorner[1] = sqrtf(((TriangleX[1] - xp) * (TriangleX[1] - xp)) +
                                                        ((TriangleY[1] - yp) * (TriangleY[1] - yp)) +
                                                        ((TriangleZ[1] - zp) * (TriangleZ[1] - zp)));
                                DistToCorner[2] = sqrtf(((TriangleX[2] - xp) * (TriangleX[2] - xp)) +
                                                        ((TriangleY[2] - yp) * (TriangleY[2] - yp)) +
                                                        ((TriangleZ[2] - zp) * (TriangleZ[2] - zp)));

                                int x, y, z;
                                x = (DistToCorner[0] < DistToCorner[1]);
                                y = (DistToCorner[1] < DistToCorner[2]);
                                z = (DistToCorner[2] < DistToCorner[0]);

                                int idx = 2 * (z - y) * z + (y - x) * y;
                                float mindisttocorner = DistToCorner[idx];
                                float xc = TriangleX[idx];
                                float yc = TriangleY[idx];
                                float zc = TriangleZ[idx];

                                float x1 = TriangleX[(idx + 1) % 3];
                                float y1 = TriangleY[(idx + 1) % 3];
                                float z1 = TriangleZ[(idx + 1) % 3];
                                float x2 = TriangleX[(idx + 2) % 3];
                                float y2 = TriangleY[(idx + 2) % 3];
                                float z2 = TriangleZ[(idx + 2) % 3];

                                float D1 =
                                    sqrtf(((xp - x2) * (xp - x2)) + ((yp - y2) * (yp - y2)) + ((zp - z2) * (zp - z2)));
                                float D2 =
                                    sqrtf(((xc - x2) * (xc - x2)) + ((yc - y2) * (yc - y2)) + ((zc - z2) * (zc - z2)));
                                float D3 =
                                    sqrtf(((xp - x1) * (xp - x1)) + ((yp - y1) * (yp - y1)) + ((zp - z1) * (zp - z1)));
                                float D4 =
                                    sqrtf(((xc - x1) * (xc - x1)) + ((yc - y1) * (yc - y1)) + ((zc - z1) * (zc - z1)));

                                float I1 = 0;
                                float I2 = D2;
                                float J1 = 0;
                                float J2 = D4;
                                // If minimum distance to corner = 0, the octahedron corner captured the new cell center
                                if (mindisttocorner != 0) {
                                    I1 = D1 * ((xp - x2) * (xc - x2) + (yp - y2) * (yc - y2) + (zp - z2) * (zc - z2)) /
                                         (D1 * D2);
                                    I2 = D2 - I1;
                                    J1 = D3 * ((xp - x1) * (xc - x1) + (yp - y1) * (yc - y1) + (zp - z1) * (zc - z1)) /
                                         (D3 * D4);
                                    J2 = D4 - J1;
                                }
                                float L12 = 0.5 * (fmin(I1, sqrtf(3.0)) + fmin(I2, sqrtf(3.0)));
                                float L13 = 0.5 * (fmin(J1, sqrtf(3.0)) + fmin(J2, sqrtf(3.0)));
                                float NewODiagL = sqrtf(2.0) * fmax(L12, L13); // half diagonal length of new octahedron

                                DiagonalLength(NeighborD3D1ConvPosition) = NewODiagL;
                                // Calculate coordinates of new decentered octahedron center
                                float CaptDiag[3];
                                CaptDiag[0] = xc - cxold;
                                CaptDiag[1] = yc - cyold;
                                CaptDiag[2] = zc - czold;

                                float CaptDiagMagnitude = sqrt(CaptDiag[0] * CaptDiag[0] + CaptDiag[1] * CaptDiag[1] +
                                                               CaptDiag[2] * CaptDiag[2]);
                                float CaptDiagUV[3];
                                CaptDiagUV[0] = CaptDiag[0] / CaptDiagMagnitude;
                                CaptDiagUV[1] = CaptDiag[1] / CaptDiagMagnitude;
                                CaptDiagUV[2] = CaptDiag[2] / CaptDiagMagnitude;
                                // (cx, cy, cz) are the coordiantes of the new active cell's decentered octahedron
                                float cx = xc - NewODiagL * CaptDiagUV[0];
                                float cy = yc - NewODiagL * CaptDiagUV[1];
                                float cz = zc - NewODiagL * CaptDiagUV[2];

                                DOCenter((long int)(3) * NeighborD3D1ConvPosition) = cx;
                                DOCenter((long int)(3) * NeighborD3D1ConvPosition + (long int)(1)) = cy;
                                DOCenter((long int)(3) * NeighborD3D1ConvPosition + (long int)(2)) = cz;

                                // Get new critical diagonal length values for the newly activated cell (at array
                                // position "NeighborD3D1ConvPosition")
                                calcCritDiagonalLength(NeighborD3D1ConvPosition, xp, yp, zp, cx, cy, cz, NeighborX,
                                                       NeighborY, NeighborZ, MyOrientation, GrainUnitVector,
                                                       CritDiagonalLength);

                                if (np > 1) {

                                    double GhostGID = static_cast<double>(h);
                                    double GhostDOCX = cx;
                                    double GhostDOCY = cy;
                                    double GhostDOCZ = cz;
                                    double GhostDL = NewODiagL;
                                    // Collect data for the ghost nodes, if necessary
                                    // Data loaded into the ghost nodes is for the cell that was just captured
                                    if (DecompositionStrategy == 1)
                                        loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, BufSizeX,
                                                       MyYSlices, MyNeighborX, MyNeighborY, MyNeighborZ,
                                                       AtNorthBoundary, AtSouthBoundary, BufferSouthSend,
                                                       BufferNorthSend);
                                    else
                                        loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, BufSizeX,
                                                       BufSizeY, MyXSlices, MyYSlices, MyNeighborX, MyNeighborY,
                                                       MyNeighborZ, AtNorthBoundary, AtSouthBoundary, AtWestBoundary,
                                                       AtEastBoundary, BufferSouthSend, BufferNorthSend, BufferWestSend,
                                                       BufferEastSend, BufferNorthEastSend, BufferSouthEastSend,
                                                       BufferSouthWestSend, BufferNorthWestSend);
                                } // End if statement for serial/parallel code
                                // Only update the new cell's type once Critical Diagonal Length, Triangle Index, and
                                // Diagonal Length values have been assigned to it Avoids the race condition in which
                                // the new cell is activated, and another thread acts on the new active cell before the
                                // cell's new critical diagonal length/triangle index/diagonal length values are
                                // assigned
                                CellType(GlobalNeighborD3D1ConvPosition) = Active;
                            } // End if statement within locked capture loop
                        }     // End if statement for outer capture loop
                    }         // End if statement over neighbors on the active grid
                }             // End loop over all neighbors of this active cell
                if (DeactivateCell) {
                    // This active cell has no more neighboring cells to be captured
                    if (RemeltingYN) {
                        // Update the counter for the number of times this cell went from liquid to active to solid
                        SolidificationEventCounter(D3D1ConvPosition)++;
                        // Did the cell solidify for the last time in the layer?
                        // If so, this cell is solid - ignore until next layer (if needed)
                        // If not, update MeltTimeStep, CritTimeStep, and UndercoolingChange with values for the next
                        // solidification event, and change cell type to TempSolid
                        if (SolidificationEventCounter(D3D1ConvPosition) ==
                            NumberOfSolidificationEvents(D3D1ConvPosition)) {
                            CellType(GlobalD3D1ConvPosition) = Solid;
                        }
                        else {
                            CellType(GlobalD3D1ConvPosition) = TempSolid;
                            MeltTimeStep(GlobalD3D1ConvPosition) = (int)(LayerTimeTempHistory(
                                D3D1ConvPosition, SolidificationEventCounter(D3D1ConvPosition), 0));
                            CritTimeStep(GlobalD3D1ConvPosition) = (int)(LayerTimeTempHistory(
                                D3D1ConvPosition, SolidificationEventCounter(D3D1ConvPosition), 1));
                            UndercoolingChange(GlobalD3D1ConvPosition) =
                                LayerTimeTempHistory(D3D1ConvPosition, SolidificationEventCounter(D3D1ConvPosition), 2);
                        }
                    }
                    else {
                        // If no remelting, this cell becomes solid type - it will not change type again
                        CellType(GlobalD3D1ConvPosition) = Solid;
                    }
                }
            }
            else if (CellType(GlobalD3D1ConvPosition) == FutureActive) {
                // Successful nucleation event - this cell is becoming a new active cell
                CellType(GlobalD3D1ConvPosition) = TemporaryUpdate; // avoid operating on the new active cell before its
                                                                    // associated octahedron data is initialized

                // Location of this cell on the global grid
                int GlobalX = RankX + MyXOffset;
                int GlobalY = RankY + MyYOffset;
                int MyGrainID = GrainID(GlobalD3D1ConvPosition); // GrainID was assigned as part of Nucleation

                // Initialize new octahedron
                createNewOctahedron(D3D1ConvPosition, DiagonalLength, DOCenter, GlobalX, GlobalY, GlobalZ);
                // The orientation for the new grain will depend on its Grain ID (nucleated grains have negative GrainID
                // values)
                int MyOrientation = getGrainOrientation(MyGrainID, NGrainOrientations);
                float cx = GlobalX + 0.5;
                float cy = GlobalY + 0.5;
                float cz = GlobalZ + 0.5;
                // Calculate critical values at which this active cell leads to the activation of a neighboring liquid
                // cell. Octahedron center and cell center overlap for octahedra created as part of a new grain
                calcCritDiagonalLength(D3D1ConvPosition, cx, cy, cz, cx, cy, cz, NeighborX, NeighborY, NeighborZ,
                                       MyOrientation, GrainUnitVector, CritDiagonalLength);
                if (np > 1) {

                    double GhostGID = static_cast<double>(MyGrainID);
                    double GhostDOCX = static_cast<double>(GlobalX + 0.5);
                    double GhostDOCY = static_cast<double>(GlobalY + 0.5);
                    double GhostDOCZ = static_cast<double>(GlobalZ + 0.5);
                    double GhostDL = 0.01;
                    // Collect data for the ghost nodes, if necessary
                    if (DecompositionStrategy == 1)
                        loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, BufSizeX, MyYSlices, RankX,
                                       RankY, RankZ, AtNorthBoundary, AtSouthBoundary, BufferSouthSend,
                                       BufferNorthSend);
                    else
                        loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, BufSizeX, BufSizeY,
                                       MyXSlices, MyYSlices, RankX, RankY, RankZ, AtNorthBoundary, AtSouthBoundary,
                                       AtWestBoundary, AtEastBoundary, BufferSouthSend, BufferNorthSend, BufferWestSend,
                                       BufferEastSend, BufferNorthEastSend, BufferSouthEastSend, BufferSouthWestSend,
                                       BufferNorthWestSend);
                } // End if statement for serial/parallel code
                // Cell activation is now finished - cell type can be changed from TemporaryUpdate to Active
                CellType(GlobalD3D1ConvPosition) = Active;
            }
        });
    Kokkos::fence();
}

//*****************************************************************************/
// Jump to the next time step with work to be done, if nothing left to do in the near future
// Without remelting, the cells of interest are undercooled liquid cells, and the view checked for future work is
// CritTimeStep With remelting, the cells of interest are active cells, and the view checked for future work is
// MeltTimeStep Print intermediate output during this jump if PrintIdleMovieFrames = true
void JumpTimeStep(int &cycle, unsigned long int RemainingCellsOfInterest, unsigned long int LocalIncompleteCells,
                  ViewI FutureWorkView, int LocalActiveDomainSize, int MyXSlices, int MyYSlices, int ZBound_Low,
                  bool RemeltingYN, ViewI CellType, ViewI LayerID, int id, int layernumber, int np, int nx, int ny,
                  int nz, int MyXOffset, int MyYOffset, int ProcessorsInXDirection, int ProcessorsInYDirection,
                  ViewI GrainID, ViewI CritTimeStep, ViewF GrainUnitVector, ViewF UndercoolingChange,
                  ViewF UndercoolingCurrent, std::string OutputFile, int DecompositionStrategy, int NGrainOrientations,
                  std::string PathToOutput, int &IntermediateFileCounter, int nzActive, double deltax, float XMin,
                  float YMin, float ZMin, int NumberOfLayers, int &XSwitch, std::string TemperatureDataType,
                  bool PrintIdleMovieFrames, int MovieFrameInc, int FinishTimeStep = 0) {

    MPI_Bcast(&RemainingCellsOfInterest, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    if (RemainingCellsOfInterest == 0) {
        // If this rank still has cells that will later undergo transformation (LocalIncompleteCells > 0), check when
        // the next superheated cells go below the liquidus (no remelting) or check when the next solid cells go above
        // the liquidus (remelting) Otherwise, assign the largest possible time step as the next time work needs to be
        // done on the rank
        unsigned long int NextWorkTimeStep;
        if (LocalIncompleteCells > 0) {
            Kokkos::parallel_reduce(
                "CheckNextTSForWork", LocalActiveDomainSize,
                KOKKOS_LAMBDA(const int &D3D1ConvPosition, unsigned long int &tempv) {
                    int RankZ = D3D1ConvPosition / (MyXSlices * MyYSlices);
                    int Rem = D3D1ConvPosition % (MyXSlices * MyYSlices);
                    int RankX = Rem / MyYSlices;
                    int RankY = Rem % MyYSlices;
                    int GlobalZ = RankZ + ZBound_Low;
                    int GlobalD3D1ConvPosition = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                    unsigned long int NextWorkTimeStep_ThisCell =
                        (unsigned long int)(FutureWorkView(GlobalD3D1ConvPosition));
                    // remelting/no remelting criteria for a cell to be associated with future work
                    if (((!(RemeltingYN)) && (CellType(GlobalD3D1ConvPosition) == Liquid) &&
                         (LayerID(GlobalD3D1ConvPosition) == layernumber)) ||
                        ((RemeltingYN) && (CellType(GlobalD3D1ConvPosition) == TempSolid))) {
                        if (NextWorkTimeStep_ThisCell < tempv)
                            tempv = NextWorkTimeStep_ThisCell;
                    }
                },
                Kokkos::Min<unsigned long int>(NextWorkTimeStep));
        }
        else
            NextWorkTimeStep = INT_MAX;

        unsigned long int GlobalNextWorkTimeStep;
        MPI_Allreduce(&NextWorkTimeStep, &GlobalNextWorkTimeStep, 1, MPI_UNSIGNED_LONG, MPI_MIN, MPI_COMM_WORLD);
        if ((GlobalNextWorkTimeStep - cycle) > 5000) {
            if (PrintIdleMovieFrames) {
                // Print any movie frames that occur during the skipped time steps
                for (unsigned long int cycle_jump = cycle + 1; cycle_jump < GlobalNextWorkTimeStep; cycle_jump++) {
                    if (cycle_jump % MovieFrameInc == 0) {
                        // Print current state of ExaCA simulation (up to and including the current layer's data)
                        // Host mirrors of CellType and GrainID are not maintained - pass device views and perform
                        // copy inside of subroutine
                        PrintExaCAData(id, layernumber, np, nx, ny, nz, MyXSlices, MyYSlices, MyXOffset, MyYOffset,
                                       ProcessorsInXDirection, ProcessorsInYDirection, GrainID, CritTimeStep,
                                       GrainUnitVector, LayerID, CellType, UndercoolingChange, UndercoolingCurrent,
                                       OutputFile, DecompositionStrategy, NGrainOrientations, PathToOutput, 0, false,
                                       false, false, true, false, IntermediateFileCounter, ZBound_Low, nzActive, deltax,
                                       XMin, YMin, ZMin, NumberOfLayers);
                        IntermediateFileCounter++;
                    }
                }
            }
            // Jump to next time step when solidification starts again
            cycle = GlobalNextWorkTimeStep - 1;
            if (id == 0)
                std::cout << "Jumping to cycle " << cycle + 1 << std::endl;
        }
        // If all cells have cooled below solidus, jump to next layer (no remelting/using input temp data only)
        if ((!(RemeltingYN)) && (TemperatureDataType == "R") && (cycle >= FinishTimeStep))
            XSwitch = 1;
    }
}

//*****************************************************************************/
// Prints intermediate code output to stdout, checks to see if solidification is complete
void IntermediateOutputAndCheck(int id, int np, int &cycle, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset,
                                int LocalDomainSize, int LocalActiveDomainSize, int nx, int ny, int nz, int nzActive,
                                double deltax, float XMin, float YMin, float ZMin, int DecompositionStrategy,
                                int ProcessorsInXDirection, int ProcessorsInYDirection,
                                int SuccessfulNucEvents_ThisRank, int &XSwitch, ViewI CellType, ViewI CritTimeStep,
                                ViewI GrainID, std::string TemperatureDataType, int *FinishTimeStep, int layernumber,
                                int, int ZBound_Low, int NGrainOrientations, ViewI LayerID, ViewF GrainUnitVector,
                                ViewF UndercoolingChange, ViewF UndercoolingCurrent, std::string PathToOutput,
                                std::string OutputFile, bool PrintIdleMovieFrames, int MovieFrameInc,
                                int &IntermediateFileCounter, int NumberOfLayers) {

    unsigned long int LocalSuperheatedCells;
    unsigned long int LocalUndercooledCells;
    unsigned long int LocalActiveCells;
    unsigned long int LocalSolidCells;
    unsigned long int LocalRemainingLiquidCells;
    Kokkos::parallel_reduce(
        LocalDomainSize,
        KOKKOS_LAMBDA(const int &D3D1ConvPosition, unsigned long int &sum_superheated,
                      unsigned long int &sum_undercooled, unsigned long int &sum_active, unsigned long int &sum_solid,
                      unsigned long int &sum_remaining_liquid) {
            if (LayerID(D3D1ConvPosition) == layernumber) {
                if (CellType(D3D1ConvPosition) == Liquid) {
                    if (CritTimeStep(D3D1ConvPosition) > cycle)
                        sum_superheated += 1;
                    else
                        sum_undercooled += 1;
                }
                else if (CellType(D3D1ConvPosition) == Active)
                    sum_active += 1;
                else if (CellType(D3D1ConvPosition) == Solid)
                    sum_solid += 1;
            }
            else {
                if (CellType(D3D1ConvPosition) == Liquid)
                    sum_remaining_liquid += 1;
            }
        },
        LocalSuperheatedCells, LocalUndercooledCells, LocalActiveCells, LocalSolidCells, LocalRemainingLiquidCells);

    unsigned long int Global_SuccessfulNucEvents_ThisRank = 0;
    unsigned long int GlobalSuperheatedCells, GlobalUndercooledCells, GlobalActiveCells, GlobalSolidCells,
        GlobalRemainingLiquidCells;
    MPI_Reduce(&LocalSuperheatedCells, &GlobalSuperheatedCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalUndercooledCells, &GlobalUndercooledCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalActiveCells, &GlobalActiveCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalSolidCells, &GlobalSolidCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalRemainingLiquidCells, &GlobalRemainingLiquidCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&SuccessfulNucEvents_ThisRank, &Global_SuccessfulNucEvents_ThisRank, 1, MPI_INT, MPI_SUM, 0,
               MPI_COMM_WORLD);

    if (id == 0) {
        std::cout << "cycle = " << cycle << " Superheated liquid cells = " << GlobalSuperheatedCells
                  << " Undercooled liquid cells = " << GlobalUndercooledCells
                  << " Number of nucleation events this layer " << Global_SuccessfulNucEvents_ThisRank
                  << " Remaining liquid cells in future layers of this simulation = " << GlobalRemainingLiquidCells
                  << std::endl;
        if (GlobalSuperheatedCells + GlobalUndercooledCells == 0)
            XSwitch = 1;
    }
    MPI_Bcast(&XSwitch, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // If an appropraite problem type/solidification is not finished, jump to the next time step with work to be done,
    // if nothing left to do in the near future
    if ((XSwitch == 0) && ((TemperatureDataType == "R") || (TemperatureDataType == "S")))
        JumpTimeStep(cycle, GlobalUndercooledCells, LocalSuperheatedCells, CritTimeStep, LocalActiveDomainSize,
                     MyXSlices, MyYSlices, ZBound_Low, false, CellType, LayerID, id, layernumber, np, nx, ny, nz,
                     MyXOffset, MyYOffset, ProcessorsInXDirection, ProcessorsInYDirection, GrainID, CritTimeStep,
                     GrainUnitVector, UndercoolingChange, UndercoolingCurrent, OutputFile, DecompositionStrategy,
                     NGrainOrientations, PathToOutput, IntermediateFileCounter, nzActive, deltax, XMin, YMin, ZMin,
                     NumberOfLayers, XSwitch, TemperatureDataType, PrintIdleMovieFrames, MovieFrameInc,
                     FinishTimeStep[layernumber]);
}

//*****************************************************************************/
// Prints intermediate code output to stdout (intermediate output collected and printed is different than without
// remelting) and checks to see if solidification is complete in the case where cells can solidify multiple times
void IntermediateOutputAndCheck_Remelt(
    int id, int np, int &cycle, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int LocalActiveDomainSize,
    int nx, int ny, int nz, int nzActive, double deltax, float XMin, float YMin, float ZMin, int DecompositionStrategy,
    int ProcessorsInXDirection, int ProcessorsInYDirection, int SuccessfulNucEvents_ThisRank, int &XSwitch,
    ViewI CellType, ViewI CritTimeStep, ViewI GrainID, std::string TemperatureDataType, int layernumber, int,
    int ZBound_Low, int NGrainOrientations, ViewI LayerID, ViewF GrainUnitVector, ViewF UndercoolingChange,
    ViewF UndercoolingCurrent, std::string PathToOutput, std::string OutputFile, bool PrintIdleMovieFrames,
    int MovieFrameInc, int &IntermediateFileCounter, int NumberOfLayers, ViewI MeltTimeStep) {

    unsigned long int LocalSuperheatedCells;
    unsigned long int LocalUndercooledCells;
    unsigned long int LocalActiveCells;
    unsigned long int LocalTempSolidCells;
    unsigned long int LocalFinishedSolidCells;
    Kokkos::parallel_reduce(
        LocalActiveDomainSize,
        KOKKOS_LAMBDA(const int &D3D1ConvPosition, unsigned long int &sum_superheated,
                      unsigned long int &sum_undercooled, unsigned long int &sum_active,
                      unsigned long int &sum_temp_solid, unsigned long int &sum_finished_solid) {
            int GlobalD3D1ConvPosition = D3D1ConvPosition + ZBound_Low * MyXSlices * MyYSlices;
            if (CellType(GlobalD3D1ConvPosition) == Liquid) {
                if (CritTimeStep(GlobalD3D1ConvPosition) > cycle)
                    sum_superheated += 1;
                else
                    sum_undercooled += 1;
            }
            else if (CellType(GlobalD3D1ConvPosition) == Active)
                sum_active += 1;
            else if (CellType(GlobalD3D1ConvPosition) == TempSolid)
                sum_temp_solid += 1;
            else if (CellType(GlobalD3D1ConvPosition) == Solid)
                sum_finished_solid += 1;
        },
        LocalSuperheatedCells, LocalUndercooledCells, LocalActiveCells, LocalTempSolidCells, LocalFinishedSolidCells);

    unsigned long int Global_SuccessfulNucEvents_ThisRank = 0;
    unsigned long int GlobalSuperheatedCells, GlobalUndercooledCells, GlobalActiveCells, GlobalTempSolidCells,
        GlobalFinishedSolidCells;
    MPI_Reduce(&LocalSuperheatedCells, &GlobalSuperheatedCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalUndercooledCells, &GlobalUndercooledCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalActiveCells, &GlobalActiveCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalTempSolidCells, &GlobalTempSolidCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalFinishedSolidCells, &GlobalFinishedSolidCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&SuccessfulNucEvents_ThisRank, &Global_SuccessfulNucEvents_ThisRank, 1, MPI_INT, MPI_SUM, 0,
               MPI_COMM_WORLD);

    if (id == 0) {
        std::cout << "cycle = " << cycle << " in layer " << layernumber
                  << ": Superheated liquid cells = " << GlobalSuperheatedCells
                  << " Undercooled liquid cells = " << GlobalUndercooledCells
                  << " Number of nucleation events this layer " << Global_SuccessfulNucEvents_ThisRank
                  << " cells to undergo at least one more solidification event = " << GlobalTempSolidCells
                  << " cells that are finished with solidification = " << GlobalFinishedSolidCells << std::endl;
        if (GlobalSuperheatedCells + GlobalUndercooledCells + GlobalActiveCells + GlobalTempSolidCells == 0)
            XSwitch = 1;
    }
    MPI_Bcast(&XSwitch, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if ((XSwitch == 0) && ((TemperatureDataType == "R") || (TemperatureDataType == "S")))
        JumpTimeStep(cycle, GlobalActiveCells, LocalTempSolidCells, MeltTimeStep, LocalActiveDomainSize, MyXSlices,
                     MyYSlices, ZBound_Low, true, CellType, LayerID, id, layernumber, np, nx, ny, nz, MyXOffset,
                     MyYOffset, ProcessorsInXDirection, ProcessorsInYDirection, GrainID, CritTimeStep, GrainUnitVector,
                     UndercoolingChange, UndercoolingCurrent, OutputFile, DecompositionStrategy, NGrainOrientations,
                     PathToOutput, IntermediateFileCounter, nzActive, deltax, XMin, YMin, ZMin, NumberOfLayers, XSwitch,
                     TemperatureDataType, PrintIdleMovieFrames, MovieFrameInc);
}
