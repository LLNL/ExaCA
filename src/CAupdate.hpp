// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_UPDATE_HPP
#define EXACA_UPDATE_HPP

#include "CAfunctions.hpp"
#include "CAghostnodes.hpp"
#include "CAinterfacialresponse.hpp"
#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

#include <string>

// Using for compatibility with device math functions.
using std::min;

// Assign octahedron a small initial size, and a center location
template <typename ViewType>
KOKKOS_INLINE_FUNCTION void createNewOctahedron(int D3D1ConvPosition, ViewType DiagonalLength, ViewType DOCenter,
                                                int GlobalX, int GlobalY, int GlobalZ) {
    DiagonalLength(D3D1ConvPosition) = 0.01;
    DOCenter(3 * D3D1ConvPosition) = GlobalX + 0.5;
    DOCenter(3 * D3D1ConvPosition + 1) = GlobalY + 0.5;
    DOCenter(3 * D3D1ConvPosition + 2) = GlobalZ + 0.5;
}

// For the newly active cell located at 1D array position D3D1ConvPosition (3D center coordinate of xp, yp, zp), update
// CritDiagonalLength values for cell capture of neighboring cells. The octahedron has a center located at (cx, cy, cz)
template <typename ViewType>
KOKKOS_INLINE_FUNCTION void calcCritDiagonalLength(int D3D1ConvPosition, float xp, float yp, float zp, float cx,
                                                   float cy, float cz, NList NeighborX, NList NeighborY,
                                                   NList NeighborZ, int MyOrientation, ViewType GrainUnitVector,
                                                   ViewType CritDiagonalLength) {
    // Calculate critical octahedron diagonal length to activate nearest neighbor.
    // First, calculate the unique planes (4) associated with all octahedron faces (8)
    // Then just look at distance between face and the point of interest (cell center of
    // neighbor). The critical diagonal length will be the maximum of these (since all other
    // planes will have passed over the point by then
    // ... meaning it must be in the octahedron)
    double Fx[4], Fy[4], Fz[4];

    Fx[0] = GrainUnitVector(9 * MyOrientation) + GrainUnitVector(9 * MyOrientation + 3) +
            GrainUnitVector(9 * MyOrientation + 6);
    Fx[1] = GrainUnitVector(9 * MyOrientation) - GrainUnitVector(9 * MyOrientation + 3) +
            GrainUnitVector(9 * MyOrientation + 6);
    Fx[2] = GrainUnitVector(9 * MyOrientation) + GrainUnitVector(9 * MyOrientation + 3) -
            GrainUnitVector(9 * MyOrientation + 6);
    Fx[3] = GrainUnitVector(9 * MyOrientation) - GrainUnitVector(9 * MyOrientation + 3) -
            GrainUnitVector(9 * MyOrientation + 6);

    Fy[0] = GrainUnitVector(9 * MyOrientation + 1) + GrainUnitVector(9 * MyOrientation + 4) +
            GrainUnitVector(9 * MyOrientation + 7);
    Fy[1] = GrainUnitVector(9 * MyOrientation + 1) - GrainUnitVector(9 * MyOrientation + 4) +
            GrainUnitVector(9 * MyOrientation + 7);
    Fy[2] = GrainUnitVector(9 * MyOrientation + 1) + GrainUnitVector(9 * MyOrientation + 4) -
            GrainUnitVector(9 * MyOrientation + 7);
    Fy[3] = GrainUnitVector(9 * MyOrientation + 1) - GrainUnitVector(9 * MyOrientation + 4) -
            GrainUnitVector(9 * MyOrientation + 7);

    Fz[0] = GrainUnitVector(9 * MyOrientation + 2) + GrainUnitVector(9 * MyOrientation + 5) +
            GrainUnitVector(9 * MyOrientation + 8);
    Fz[1] = GrainUnitVector(9 * MyOrientation + 2) - GrainUnitVector(9 * MyOrientation + 5) +
            GrainUnitVector(9 * MyOrientation + 8);
    Fz[2] = GrainUnitVector(9 * MyOrientation + 2) + GrainUnitVector(9 * MyOrientation + 5) -
            GrainUnitVector(9 * MyOrientation + 8);
    Fz[3] = GrainUnitVector(9 * MyOrientation + 2) - GrainUnitVector(9 * MyOrientation + 5) -
            GrainUnitVector(9 * MyOrientation + 8);

    for (int n = 0; n < 26; n++) {
        float x0 = xp + NeighborX[n] - cx;
        float y0 = yp + NeighborY[n] - cy;
        float z0 = zp + NeighborZ[n] - cz;
        float D0 = x0 * Fx[0] + y0 * Fy[0] + z0 * Fz[0];
        float D1 = x0 * Fx[1] + y0 * Fy[1] + z0 * Fz[1];
        float D2 = x0 * Fx[2] + y0 * Fy[2] + z0 * Fz[2];
        float D3 = x0 * Fx[3] + y0 * Fy[3] + z0 * Fz[3];
        float Dfabs = fmax(fmax(fabs(D0), fabs(D1)), fmax(fabs(D2), fabs(D3)));
        CritDiagonalLength((long int)(26) * D3D1ConvPosition + (long int)(n)) = Dfabs;
    }
}

// Decentered octahedron algorithm for the capture of new interface cells by grains
template <typename IRFtype>
void CellCapture(int, int np, int, int, int, int nx, int MyYSlices, IRFtype irf, int MyYOffset, NList NeighborX,
                 NList NeighborY, NList NeighborZ, ViewI CritTimeStep, ViewF UndercoolingCurrent,
                 ViewF UndercoolingChange, ViewF GrainUnitVector, ViewF CritDiagonalLength, ViewF DiagonalLength,
                 ViewI CellType, ViewF DOCenter, ViewI GrainID, int NGrainOrientations, Buffer2D BufferNorthSend,
                 Buffer2D BufferSouthSend, int BufSizeX, int ZBound_Low, int nzActive, int, ViewI SteeringVector,
                 ViewI numSteer, ViewI_H numSteer_Host, bool AtNorthBoundary, bool AtSouthBoundary,
                 ViewI SolidificationEventCounter, ViewI MeltTimeStep, ViewF3D LayerTimeTempHistory,
                 ViewI NumberOfSolidificationEvents, bool RemeltingYN) {

    // Loop over list of active and soon-to-be active cells, potentially performing cell capture events and updating
    // cell types
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
            if (CellType(GlobalD3D1ConvPosition) == Active) {
                // Update local diagonal length of active cell
                double LocU = UndercoolingCurrent(GlobalD3D1ConvPosition);
                LocU = min(210.0, LocU);
                double V = irf->compute(LocU);
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
                        long int GlobalNeighborD3D1ConvPosition =
                            (MyNeighborZ + ZBound_Low) * nx * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
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
                                    loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, BufSizeX,
                                                   MyYSlices, MyNeighborX, MyNeighborY, MyNeighborZ, AtNorthBoundary,
                                                   AtSouthBoundary, BufferSouthSend, BufferNorthSend);
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
                    loadghostnodes(GhostGID, GhostDOCX, GhostDOCY, GhostDOCZ, GhostDL, BufSizeX, MyYSlices, GlobalX,
                                   RankY, RankZ, AtNorthBoundary, AtSouthBoundary, BufferSouthSend, BufferNorthSend);
                } // End if statement for serial/parallel code
                // Cell activation is now finished - cell type can be changed from TemporaryUpdate to Active
                CellType(GlobalD3D1ConvPosition) = Active;
            }
        });
    Kokkos::fence();
}

void Nucleation(int cycle, int &SuccessfulNucEvents_ThisRank, int &NucleationCounter, int PossibleNuclei_ThisRank,
                ViewI_H NucleationTimes_H, ViewI NucleiLocations, ViewI NucleiGrainID, ViewI CellType, ViewI GrainID,
                int ZBound_Low, int nx, int MyYSlices, ViewI SteeringVector, ViewI numSteer_G);
void FillSteeringVector_NoRemelt(int cycle, int LocalActiveDomainSize, int nx, int MyYSlices, ViewI CritTimeStep,
                                 ViewF UndercoolingCurrent, ViewF UndercoolingChange, ViewI CellType, int ZBound_Low,
                                 int layernumber, ViewI LayerID, ViewI SteeringVector, ViewI numSteer_G,
                                 ViewI_H numSteer_H);
void FillSteeringVector_Remelt(int cycle, int LocalActiveDomainSize, int nx, int MyYSlices, NList NeighborX,
                               NList NeighborY, NList NeighborZ, ViewI CritTimeStep, ViewF UndercoolingCurrent,
                               ViewF UndercoolingChange, ViewI CellType, ViewI GrainID, int ZBound_Low, int nzActive,
                               ViewI SteeringVector, ViewI numSteer, ViewI_H numSteer_Host, ViewI MeltTimeStep,
                               int BufSizeX, bool AtNorthBoundary, bool AtSouthBoundary, Buffer2D BufferNorthSend,
                               Buffer2D BufferSouthSend);
void JumpTimeStep(int &cycle, unsigned long int RemainingCellsOfInterest, ViewI FutureWorkView,
                  unsigned long int LocalIncompleteCells, int LocalActiveDomainSize, int MyYSlices, int ZBound_Low,
                  bool RemeltingYN, ViewI CellType, ViewI LayerID, int id, int layernumber, int np, int nx, int ny,
                  int nz, int MyYOffset, ViewI GrainID, ViewI CritTimeStep, ViewF GrainUnitVector,
                  ViewF UndercoolingChange, ViewF UndercoolingCurrent, std::string OutputFile,
                  int DecompositionStrategy, int NGrainOrientations, std::string PathToOutput,
                  int &IntermediateFileCounter, int nzActive, double deltax, float XMin, float YMin, float ZMin,
                  int NumberOfLayers, int &XSwitch, std::string TemperatureDataType, bool PrintIdleMovieFrames,
                  int MovieFrameInc, bool PrintBinary, int FinishTimeStep);
void IntermediateOutputAndCheck(int id, int np, int &cycle, int MyYSlices, int MyYOffset, int LocalDomainSize,
                                int LocalActiveDomainSize, int nx, int ny, int nz, int nzActive, double deltax,
                                float XMin, float YMin, float ZMin, int SuccessfulNucEvents_ThisRank, int &XSwitch,
                                ViewI CellType, ViewI CritTimeStep, ViewI GrainID, std::string TemperatureDataType,
                                int *FinishTimeStep, int layernumber, int, int ZBound_Low, int NGrainOrientations,
                                ViewI LayerID, ViewF GrainUnitVector, ViewF UndercoolingChange,
                                ViewF UndercoolingCurrent, std::string PathToOutput, std::string OutputFile,
                                bool PrintIdleMovieFrames, int MovieFrameInc, int &IntermediateFileCounter,
                                int NumberOfLayers, bool PrintBinary);
void IntermediateOutputAndCheck_Remelt(
    int id, int np, int &cycle, int MyYSlices, int MyYOffset, int LocalActiveDomainSize, int nx, int ny, int nz,
    int nzActive, double deltax, float XMin, float YMin, float ZMin, int SuccessfulNucEvents_ThisRank, int &XSwitch,
    ViewI CellType, ViewI CritTimeStep, ViewI GrainID, std::string TemperatureDataType, int layernumber, int,
    int ZBound_Low, int NGrainOrientations, ViewI LayerID, ViewF GrainUnitVector, ViewF UndercoolingChange,
    ViewF UndercoolingCurrent, std::string PathToOutput, std::string OutputFile, bool PrintIdleMovieFrames,
    int MovieFrameInc, int &IntermediateFileCounter, int NumberOfLayers, ViewI MeltTimeStep, bool PrintBinary);

#endif
