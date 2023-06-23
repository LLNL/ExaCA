// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
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
// Determine which cells are associated with the "steering vector" of cells that are either active, or becoming active
// this time step
void FillSteeringVector_NoRemelt(int cycle, int LocalActiveDomainSize, int nx, int MyYSlices, ViewI CritTimeStep,
                                 ViewF UndercoolingCurrent, ViewF UndercoolingChange, CellData &cellData,
                                 int ZBound_Low, int layernumber, ViewI SteeringVector, ViewI numSteer,
                                 ViewI_H numSteer_Host) {

    // Cells associated with this layer that are not solid type but have passed the liquidus (crit time step) have their
    // undercooling values updated Cells that meet the aforementioned criteria and are active type should be added to
    // the steering vector
    ViewI CellType = cellData.getCellTypeSubview();
    ViewI LayerID = cellData.getLayerIDSubview();
    Kokkos::parallel_for(
        "FillSV", LocalActiveDomainSize, KOKKOS_LAMBDA(const int &D3D1ConvPosition) {
            // Cells of interest for the CA
            int RankZ = D3D1ConvPosition / (nx * MyYSlices);
            int Rem = D3D1ConvPosition % (nx * MyYSlices);
            int RankX = Rem / MyYSlices;
            int RankY = Rem % MyYSlices;
            int GlobalZ = RankZ + ZBound_Low;
            int GlobalD3D1ConvPosition = GlobalZ * nx * MyYSlices + RankX * MyYSlices + RankY;
            int cellType = CellType(D3D1ConvPosition);
            // TODO: layer check no longer needed here
            int layerCheck = (LayerID(D3D1ConvPosition) <= layernumber);
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
void FillSteeringVector_Remelt(int cycle, int LocalActiveDomainSize, int nx, int MyYSlices, NList NeighborX,
                               NList NeighborY, NList NeighborZ, ViewI CritTimeStep, ViewF UndercoolingCurrent,
                               ViewF UndercoolingChange, CellData &cellData, int ZBound_Low, int nzActive,
                               ViewI SteeringVector, ViewI numSteer, ViewI_H numSteer_Host, ViewI MeltTimeStep,
                               ViewI SolidificationEventCounter, ViewI NumberOfSolidificationEvents,
                               ViewF3D LayerTimeTempHistory) {

    ViewI CellType = cellData.getCellTypeSubview();
    ViewI GrainID = cellData.getGrainIDSubview();
    Kokkos::parallel_for(
        "FillSV_RM", LocalActiveDomainSize, KOKKOS_LAMBDA(const int &D3D1ConvPosition) {
            // Coordinate of this cell on the "global" (all cells in the Z direction) grid
            int RankZ = D3D1ConvPosition / (nx * MyYSlices);
            int Rem = D3D1ConvPosition % (nx * MyYSlices);
            int RankX = Rem / MyYSlices;
            int RankY = Rem % MyYSlices;
            int GlobalZ = RankZ + ZBound_Low;
            int GlobalD3D1ConvPosition = GlobalZ * nx * MyYSlices + RankX * MyYSlices + RankY;

            int cellType = CellType(D3D1ConvPosition);
            bool isNotSolid = ((cellType != TempSolid) && (cellType != Solid));
            bool atMeltTime = (cycle == MeltTimeStep(GlobalD3D1ConvPosition));
            bool atCritTime = (cycle == CritTimeStep(GlobalD3D1ConvPosition));
            bool pastCritTime = (cycle > CritTimeStep(GlobalD3D1ConvPosition));
            if (atMeltTime) {
                CellType(D3D1ConvPosition) = Liquid;
                UndercoolingCurrent(GlobalD3D1ConvPosition) = 0.0;
                // If this cell melts at least one more time after the melting event that just took place, replace the
                // value for melt time step with the next time step this cell goes above the liquidus
                if (cellType != TempSolid) {
                    // This cell either hasn't started or hasn't finished the previous solidification event, but has
                    // undergone melting - increment the solidification counter to move on to the next
                    // melt-solidification event and replace CritTimeStep and UndercoolingChange with the values
                    // corresponding to this next event
                    SolidificationEventCounter(D3D1ConvPosition)++;
                    CritTimeStep(GlobalD3D1ConvPosition) = static_cast<int>(
                        LayerTimeTempHistory(D3D1ConvPosition, SolidificationEventCounter(D3D1ConvPosition), 1));
                    UndercoolingChange(GlobalD3D1ConvPosition) =
                        LayerTimeTempHistory(D3D1ConvPosition, SolidificationEventCounter(D3D1ConvPosition), 2);
                    // If this cell melts at least one more time after the melting event that just took place, replace
                    // the value for melt time step with the next time step this cell goes above the liquidus
                    if ((SolidificationEventCounter(D3D1ConvPosition) + 1) !=
                        NumberOfSolidificationEvents(D3D1ConvPosition)) {
                        MeltTimeStep(GlobalD3D1ConvPosition) = static_cast<int>(LayerTimeTempHistory(
                            D3D1ConvPosition, SolidificationEventCounter(D3D1ConvPosition) + 1, 0));
                    }
                }
                else if ((SolidificationEventCounter(D3D1ConvPosition) + 1) !=
                         NumberOfSolidificationEvents(D3D1ConvPosition)) {
                    // If this cell melts at least one more time after the melting event that just took place, replace
                    // the value for melt time step with the next time step this cell goes above the liquidus
                    MeltTimeStep(GlobalD3D1ConvPosition) = static_cast<int>(
                        LayerTimeTempHistory(D3D1ConvPosition, SolidificationEventCounter(D3D1ConvPosition) + 1, 0));
                }
                // Any adjacent active cells should also be remelted, as these cells are more likely heating up than
                // cooling down These are converted to the temporary FutureLiquid state, to be later iterated over and
                // loaded into the steering vector as necessary
                for (int l = 0; l < 26; l++) {
                    // "l" correpsponds to the specific neighboring cell
                    // Local coordinates of adjacent cell center
                    int MyNeighborX = RankX + NeighborX[l];
                    int MyNeighborY = RankY + NeighborY[l];
                    int MyNeighborZ = RankZ + NeighborZ[l];
                    if ((MyNeighborX >= 0) && (MyNeighborX < nx) && (MyNeighborY >= 0) && (MyNeighborY < MyYSlices) &&
                        (MyNeighborZ < nzActive) && (MyNeighborZ >= 0)) {
                        int NeighborD3D1ConvPosition =
                            MyNeighborZ * nx * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
                        if (CellType(NeighborD3D1ConvPosition) == Active) {
                            CellType(NeighborD3D1ConvPosition) = FutureLiquid;
                            SteeringVector(Kokkos::atomic_fetch_add(&numSteer(0), 1)) = NeighborD3D1ConvPosition;
                        }
                    }
                }
            }
            else if ((isNotSolid) && (pastCritTime)) {
                // Update cell undercooling
                UndercoolingCurrent(GlobalD3D1ConvPosition) += UndercoolingChange(GlobalD3D1ConvPosition);
                if (cellType == Active) {
                    // Add active cells below liquidus to steering vector
                    SteeringVector(Kokkos::atomic_fetch_add(&numSteer(0), 1)) = D3D1ConvPosition;
                }
            }
            else if ((atCritTime) && (cellType == Liquid) && (GrainID(D3D1ConvPosition) != 0)) {
                // If this cell has cooled to the liquidus temperature, borders at least one solid/tempsolid cell, and
                // is part of a grain, it should become active. This only needs to be checked on the time step where the
                // cell reaches the liquidus, not every time step beyond this
                for (int l = 0; l < 26; l++) {
                    // "l" correpsponds to the specific neighboring cell
                    // Local coordinates of adjacent cell center
                    int MyNeighborX = RankX + NeighborX[l];
                    int MyNeighborY = RankY + NeighborY[l];
                    int MyNeighborZ = RankZ + NeighborZ[l];
                    if ((MyNeighborX >= 0) && (MyNeighborX < nx) && (MyNeighborY >= 0) && (MyNeighborY < MyYSlices) &&
                        (MyNeighborZ < nzActive) && (MyNeighborZ >= 0)) {
                        int NeighborD3D1ConvPosition =
                            MyNeighborZ * nx * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
                        if ((CellType(NeighborD3D1ConvPosition) == TempSolid) ||
                            (CellType(NeighborD3D1ConvPosition) == Solid) || (RankZ == 0)) {
                            // Cell activation to be performed as part of steering vector
                            l = 26;
                            SteeringVector(Kokkos::atomic_fetch_add(&numSteer(0), 1)) = D3D1ConvPosition;
                            CellType(D3D1ConvPosition) =
                                FutureActive; // this cell cannot be captured - is being activated
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
void CellCapture(int, int np, int, int, int, int nx, int MyYSlices, InterfacialResponseFunction irf, int MyYOffset,
                 NList NeighborX, NList NeighborY, NList NeighborZ, ViewI CritTimeStep, ViewF UndercoolingCurrent,
                 ViewF UndercoolingChange, ViewF GrainUnitVector, ViewF CritDiagonalLength, ViewF DiagonalLength,
                 CellData &cellData, ViewF DOCenter, int NGrainOrientations, Buffer2D BufferNorthSend,
                 Buffer2D BufferSouthSend, ViewI SendSizeNorth, ViewI SendSizeSouth, int ZBound_Low, int nzActive, int,
                 ViewI SteeringVector, ViewI numSteer, ViewI_H numSteer_Host, bool AtNorthBoundary,
                 bool AtSouthBoundary, ViewI SolidificationEventCounter, ViewF3D LayerTimeTempHistory,
                 ViewI NumberOfSolidificationEvents, int &BufSize) {

    // Loop over list of active and soon-to-be active cells, potentially performing cell capture events and updating
    // cell types
    ViewI CellType = cellData.getCellTypeSubview();
    ViewI GrainID = cellData.getGrainIDSubview();
    Kokkos::parallel_for(
        "CellCapture", numSteer_Host(0), KOKKOS_LAMBDA(const int &num) {
            numSteer(0) = 0;
            int D3D1ConvPosition = SteeringVector(num);
            // Cells of interest for the CA - active cells and future active cells
            int RankZ = D3D1ConvPosition / (nx * MyYSlices);
            int Rem = D3D1ConvPosition % (nx * MyYSlices);
            int GlobalX = Rem / MyYSlices;
            int RankY = Rem % MyYSlices;
            int GlobalZ = RankZ + ZBound_Low;
            int GlobalD3D1ConvPosition = GlobalZ * nx * MyYSlices + GlobalX * MyYSlices + RankY;
            if (CellType(D3D1ConvPosition) == Active) {
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
                    int MyNeighborX = GlobalX + NeighborX[l];
                    int MyNeighborY = RankY + NeighborY[l];
                    int MyNeighborZ = RankZ + NeighborZ[l];
                    // Check if neighbor is in bounds
                    if ((MyNeighborX >= 0) && (MyNeighborX < nx) && (MyNeighborY >= 0) && (MyNeighborY < MyYSlices) &&
                        (MyNeighborZ < nzActive) && (MyNeighborZ >= 0)) {
                        long int NeighborD3D1ConvPosition =
                            MyNeighborZ * nx * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
                        if (CellType(NeighborD3D1ConvPosition) == Liquid)
                            DeactivateCell = false;
                        // Capture of cell located at "NeighborD3D1ConvPosition" if this condition is satisfied
                        if ((DiagonalLength(D3D1ConvPosition) >= CritDiagonalLength(26 * D3D1ConvPosition + l)) &&
                            (CellType(NeighborD3D1ConvPosition) == Liquid)) {
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
                            int OldCellTypeValue = Kokkos::atomic_compare_exchange(&CellType(NeighborD3D1ConvPosition),
                                                                                   old_val, update_val);
                            // Only proceed if CellType was previously liquid (this current thread changed the value to
                            // TemporaryUpdate)
                            if (OldCellTypeValue == Liquid) {
                                int GlobalY = RankY + MyYOffset;
                                int h = GrainID(D3D1ConvPosition);
                                int MyOrientation = getGrainOrientation(h, NGrainOrientations);

                                // The new cell is captured by this cell's growing octahedron (Grain "h")
                                GrainID(NeighborD3D1ConvPosition) = h;

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
                                    // TODO: Test loading ghost nodes in a separate kernel, potentially adopting this
                                    // change if the slowdown is minor
                                    int GhostGID = h;
                                    float GhostDOCX = cx;
                                    float GhostDOCY = cy;
                                    float GhostDOCZ = cz;
                                    float GhostDL = NewODiagL;
                                    // Collect data for the ghost nodes, if necessary
                                    // Data loaded into the ghost nodes is for the cell that was just captured
                                    bool DataFitsInBuffer =
                                        loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL,
                                                       SendSizeNorth, SendSizeSouth, MyYSlices, MyNeighborX,
                                                       MyNeighborY, MyNeighborZ, AtNorthBoundary, AtSouthBoundary,
                                                       BufferSouthSend, BufferNorthSend, NGrainOrientations, BufSize);
                                    if (!(DataFitsInBuffer)) {
                                        // This cell's data did not fit in the buffer with current size BufSize - mark
                                        // with temporary type
                                        CellType(NeighborD3D1ConvPosition) = ActiveFailedBufferLoad;
                                    }
                                    else {
                                        // Cell activation is now finished - cell type can be changed from
                                        // TemporaryUpdate to Active
                                        CellType(NeighborD3D1ConvPosition) = Active;
                                    }
                                } // End if statement for serial/parallel code
                                else {
                                    // Only update the new cell's type once Critical Diagonal Length, Triangle Index,
                                    // and Diagonal Length values have been assigned to it Avoids the race condition in
                                    // which the new cell is activated, and another thread acts on the new active cell
                                    // before the cell's new critical diagonal length/triangle index/diagonal length
                                    // values are assigned
                                    CellType(NeighborD3D1ConvPosition) = Active;
                                }
                            } // End if statement within locked capture loop
                        }     // End if statement for outer capture loop
                    }         // End if statement over neighbors on the active grid
                }             // End loop over all neighbors of this active cell
                if (DeactivateCell) {
                    // This active cell has no more neighboring cells to be captured
                    // Update the counter for the number of times this cell went from liquid to active to solid
                    SolidificationEventCounter(D3D1ConvPosition)++;
                    // Did the cell solidify for the last time in the layer?
                    // If so, this cell is solid - ignore until next layer (if needed)
                    // If not, update CritTimeStep, and UndercoolingChange with values for the next
                    // solidification event, and change cell type to TempSolid
                    if (SolidificationEventCounter(D3D1ConvPosition) ==
                        NumberOfSolidificationEvents(D3D1ConvPosition)) {
                        CellType(D3D1ConvPosition) = Solid;
                    }
                    else {
                        CellType(D3D1ConvPosition) = TempSolid;
                        CritTimeStep(GlobalD3D1ConvPosition) = (int)(LayerTimeTempHistory(
                            D3D1ConvPosition, SolidificationEventCounter(D3D1ConvPosition), 1));
                        UndercoolingChange(GlobalD3D1ConvPosition) =
                            LayerTimeTempHistory(D3D1ConvPosition, SolidificationEventCounter(D3D1ConvPosition), 2);
                    }
                }
            }
            else if (CellType(D3D1ConvPosition) == FutureActive) {
                // Successful nucleation event - this cell is becoming a new active cell
                CellType(D3D1ConvPosition) = TemporaryUpdate; // avoid operating on the new active cell before its
                                                              // associated octahedron data is initialized

                // Location of this cell on the global grid
                int GlobalY = RankY + MyYOffset;
                int MyGrainID = GrainID(D3D1ConvPosition); // GrainID was assigned as part of Nucleation

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
                    // TODO: Test loading ghost nodes in a separate kernel, potentially adopting this change if the
                    // slowdown is minor
                    int GhostGID = MyGrainID;
                    float GhostDOCX = GlobalX + 0.5;
                    float GhostDOCY = GlobalY + 0.5;
                    float GhostDOCZ = GlobalZ + 0.5;
                    float GhostDL = 0.01;
                    // Collect data for the ghost nodes, if necessary
                    bool DataFitsInBuffer =
                        loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, SendSizeNorth, SendSizeSouth,
                                       MyYSlices, GlobalX, RankY, RankZ, AtNorthBoundary, AtSouthBoundary,
                                       BufferSouthSend, BufferNorthSend, NGrainOrientations, BufSize);
                    if (!(DataFitsInBuffer)) {
                        // This cell's data did not fit in the buffer with current size BufSize - mark with temporary
                        // type
                        CellType(D3D1ConvPosition) = ActiveFailedBufferLoad;
                    }
                    else {
                        // Cell activation is now finished - cell type can be changed from TemporaryUpdate to Active
                        CellType(D3D1ConvPosition) = Active;
                    }
                } // End if statement for serial/parallel code
                else {
                    // Cell activation is now finished - cell type can be changed from TemporaryUpdate to Active
                    CellType(D3D1ConvPosition) = Active;
                } // End if statement for serial/parallel code
            }
            else if (CellType(D3D1ConvPosition) == FutureLiquid) {
                // This type was assigned to a cell that was recently transformed from active to liquid, due to its
                // bordering of a cell above the liquidus. This information may need to be sent to other MPI ranks
                // Dummy values for first 4 arguments (Grain ID and octahedron center coordinates), 0 for diagonal
                // length
                bool DataFitsInBuffer = loadghostnodes(
                    -1, -1.0, -1.0, -1.0, 0.0, SendSizeNorth, SendSizeSouth, MyYSlices, GlobalX, RankY, RankZ,
                    AtNorthBoundary, AtSouthBoundary, BufferSouthSend, BufferNorthSend, NGrainOrientations, BufSize);
                if (!(DataFitsInBuffer)) {
                    // This cell's data did not fit in the buffer with current size BufSize - mark with temporary
                    // type
                    CellType(D3D1ConvPosition) = LiquidFailedBufferLoad;
                }
                else {
                    // Cell activation is now finished - cell type can be changed from FutureLiquid to Active
                    CellType(D3D1ConvPosition) = Liquid;
                }
            }
        });
    Kokkos::fence();
}

//*****************************************************************************/
// Jump to the next time step with work to be done, if nothing left to do in the near future
// The cells of interest are active cells, and the view checked for future work is
// MeltTimeStep Print intermediate output during this jump if PrintIdleMovieFrames = true
void JumpTimeStep(int &cycle, unsigned long int RemainingCellsOfInterest, unsigned long int LocalTempSolidCells,
                  ViewI MeltTimeStep, int LocalActiveDomainSize, int MyYSlices, int ZBound_Low, ViewI CellType,
                  ViewI LayerID, int id, int layernumber, int np, int nx, int ny, ViewI GrainID, ViewF GrainUnitVector,
                  Print print, int NGrainOrientations, int nzActive, double deltax, double XMin, double YMin,
                  double ZMin) {

    MPI_Bcast(&RemainingCellsOfInterest, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    if (RemainingCellsOfInterest == 0) {
        // If this rank still has cells that will later undergo transformation (LocalIncompleteCells > 0), check when
        // the next solid cells go above the liquidus (remelting) Otherwise, assign the largest possible time step as
        // the next time work needs to be done on the rank
        unsigned long int NextMeltTimeStep;
        if (LocalTempSolidCells > 0) {
            Kokkos::parallel_reduce(
                "CheckNextTSForWork", LocalActiveDomainSize,
                KOKKOS_LAMBDA(const int &D3D1ConvPosition, unsigned long int &tempv) {
                    int RankZ = D3D1ConvPosition / (nx * MyYSlices);
                    int Rem = D3D1ConvPosition % (nx * MyYSlices);
                    int RankX = Rem / MyYSlices;
                    int RankY = Rem % MyYSlices;
                    int GlobalZ = RankZ + ZBound_Low;
                    int GlobalD3D1ConvPosition = GlobalZ * nx * MyYSlices + RankX * MyYSlices + RankY;
                    unsigned long int NextMeltTimeStep_ThisCell =
                        static_cast<unsigned long int>(MeltTimeStep(GlobalD3D1ConvPosition));
                    // criteria for a cell to be associated with future work
                    if (CellType(D3D1ConvPosition) == TempSolid) {
                        if (NextMeltTimeStep_ThisCell < tempv)
                            tempv = NextMeltTimeStep_ThisCell;
                    }
                },
                Kokkos::Min<unsigned long int>(NextMeltTimeStep));
        }
        else
            NextMeltTimeStep = INT_MAX;

        unsigned long int GlobalNextMeltTimeStep;
        MPI_Allreduce(&NextMeltTimeStep, &GlobalNextMeltTimeStep, 1, MPI_UNSIGNED_LONG, MPI_MIN, MPI_COMM_WORLD);
        if ((GlobalNextMeltTimeStep - cycle) > 5000) {
            if (print.PrintIdleTimeSeriesFrames) {
                // Print any movie frames that occur during the skipped time steps
                for (unsigned long int cycle_jump = cycle + 1; cycle_jump < GlobalNextMeltTimeStep; cycle_jump++) {
                    if (cycle_jump % print.TimeSeriesInc == 0) {
                        // Print current state of ExaCA simulation (up to and including the current layer's data)
                        print.printIntermediateGrainMisorientation(
                            id, np, cycle, nx, ny, MyYSlices, nzActive, deltax, XMin, YMin, ZMin, GrainID, LayerID,
                            CellType, GrainUnitVector, NGrainOrientations, layernumber, ZBound_Low);
                    }
                }
            }
            // Jump to next time step when solidification starts again
            cycle = GlobalNextMeltTimeStep - 1;
            if (id == 0)
                std::cout << "Jumping to cycle " << cycle + 1 << std::endl;
        }
    }
}

//*****************************************************************************/
// Prints intermediate code output to stdout (intermediate output collected and printed is different than without
// remelting) and checks to see if solidification is complete in the case where cells can solidify multiple times
void IntermediateOutputAndCheck(int id, int np, int &cycle, int MyYSlices, int LocalActiveDomainSize, int nx, int ny,
                                int nzActive, double deltax, double XMin, double YMin, double ZMin,
                                int SuccessfulNucEvents_ThisRank, int &XSwitch, CellData &cellData, ViewI CritTimeStep,
                                std::string TemperatureDataType, int layernumber, int, int ZBound_Low,
                                int NGrainOrientations, ViewF GrainUnitVector, Print print, ViewI MeltTimeStep) {

    ViewI CellType = cellData.getCellTypeSubview();
    ViewI GrainID = cellData.getGrainIDSubview();
    ViewI LayerID = cellData.getLayerIDSubview();
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
            int GlobalD3D1ConvPosition = D3D1ConvPosition + ZBound_Low * nx * MyYSlices;
            if (CellType(D3D1ConvPosition) == Liquid) {
                if (CritTimeStep(GlobalD3D1ConvPosition) > cycle)
                    sum_superheated += 1;
                else
                    sum_undercooled += 1;
            }
            else if (CellType(D3D1ConvPosition) == Active)
                sum_active += 1;
            else if (CellType(D3D1ConvPosition) == TempSolid)
                sum_temp_solid += 1;
            else if (CellType(D3D1ConvPosition) == Solid)
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
    // Cells of interest are those currently undergoing a melting-solidification cycle
    unsigned long int RemainingCellsOfInterest = GlobalActiveCells + GlobalSuperheatedCells + GlobalUndercooledCells;
    if ((XSwitch == 0) && ((TemperatureDataType == "R") || (TemperatureDataType == "S")))
        JumpTimeStep(cycle, RemainingCellsOfInterest, LocalTempSolidCells, MeltTimeStep, LocalActiveDomainSize,
                     MyYSlices, ZBound_Low, CellType, LayerID, id, layernumber, np, nx, ny, GrainID, GrainUnitVector,
                     print, NGrainOrientations, nzActive, deltax, XMin, YMin, ZMin);
}
