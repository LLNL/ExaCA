// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAupdate.hpp"

#include "mpi.h"

#include <cmath>
#include <iostream>

// Using for compatibility with device math functions.
using std::max;
using std::min;

namespace sample { // namespace helps with name resolution in reduction identity
template <class ScalarType, int N>
struct array_type {
    ScalarType the_array[N];

    KOKKOS_INLINE_FUNCTION // Default constructor - Initialize to 0's
    array_type() {
        for (int i = 0; i < N; i++) {
            the_array[i] = 0;
        }
    }
    KOKKOS_INLINE_FUNCTION // Copy Constructor
    array_type(const array_type &rhs) {
        for (int i = 0; i < N; i++) {
            the_array[i] = rhs.the_array[i];
        }
    }
    KOKKOS_INLINE_FUNCTION // add operator
        array_type &
        operator+=(const array_type &src) {
        for (int i = 0; i < N; i++) {
            the_array[i] += src.the_array[i];
        }
        return *this;
    }
    KOKKOS_INLINE_FUNCTION // volatile add operator
        void
        operator+=(const volatile array_type &src) volatile {
        for (int i = 0; i < N; i++) {
            the_array[i] += src.the_array[i];
        }
    }
};
typedef array_type<unsigned long int, 5> ValueType; // used to simplify code below
} // namespace sample
namespace Kokkos { // reduction identity must be defined in Kokkos namespace
template <>
struct reduction_identity<sample::ValueType> {
    KOKKOS_FORCEINLINE_FUNCTION static sample::ValueType sum() { return sample::ValueType(); }
};
} // namespace Kokkos

//*****************************************************************************/
void Nucleation(int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int cycle, int &nn, ViewI CellType,
                ViewI NucleiLocations, ViewI NucleationTimes, ViewI GrainID, ViewI GrainOrientation, ViewF DOCenter,
                ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewF GrainUnitVector, ViewF CritDiagonalLength,
                ViewF DiagonalLength, int NGrainOrientations, int PossibleNuclei_ThisRank, int ZBound_Low,
                int layernumber, ViewI LayerID) {

    // Loop through local list of nucleation events - has the time step exceeded the time step for nucleation at the
    // sites?
    int NucleationThisDT = 0;
    Kokkos::parallel_reduce(
        "NucleiUpdateLoop", PossibleNuclei_ThisRank,
        KOKKOS_LAMBDA(const int &NucCounter, int &update) {
            if ((cycle >= NucleationTimes(NucCounter)) && (CellType(NucleiLocations(NucCounter)) == Liquid) &&
                (LayerID(NucleiLocations(NucCounter)) <= layernumber)) {
                // (X,Y,Z) coordinates of nucleation event, on active cell grid (RankX,RankY,RankZ) and global grid
                // (RankX,RankY,GlobalZ)
                long int GlobalD3D1ConvPosition = NucleiLocations(NucCounter);
                CellType(GlobalD3D1ConvPosition) = TemporaryUpdate;
                int GlobalZ = GlobalD3D1ConvPosition / (MyXSlices * MyYSlices);
                int Rem = GlobalD3D1ConvPosition % (MyXSlices * MyYSlices);
                int RankX = Rem / MyYSlices;
                int RankY = Rem % MyYSlices;
                int RankZ = GlobalZ - ZBound_Low;
                int D3D1ConvPosition = RankZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                // This undercooled liquid cell is now a nuclei (add to count if it isn't in the ghost nodes, to avoid
                // double counting)

                if ((RankX > 0) && (RankX < MyXSlices - 1) && (RankY > 0) && (RankY < MyYSlices - 1))
                    update++;
                int GlobalX = RankX + MyXOffset;
                int GlobalY = RankY + MyYOffset;
                int MyGrainID = GrainID(GlobalD3D1ConvPosition);

                DiagonalLength(D3D1ConvPosition) = 0.01;
                long int DOX = (long int)(3) * D3D1ConvPosition;
                long int DOY = (long int)(3) * D3D1ConvPosition + (long int)(1);
                long int DOZ = (long int)(3) * D3D1ConvPosition + (long int)(2);
                DOCenter(DOX) = GlobalX + 0.5;
                DOCenter(DOY) = GlobalY + 0.5;
                DOCenter(DOZ) = GlobalZ + 0.5;
                // The orientation for the new grain will depend on its Grain ID (nucleated grains have negative GrainID
                // values)
                int MyOrientation = GrainOrientation(((abs(MyGrainID) - 1) % NGrainOrientations));
                // Calculate critical values at which this active cell leads to the activation of a neighboring liquid
                // cell
                for (int n = 0; n < 26; n++) {

                    // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the image of the center
                    // of a neighbor cell
                    double x0 = NeighborX(n);
                    double y0 = NeighborY(n);
                    double z0 = NeighborZ(n);

                    // mag0 is the magnitude of (x0,y0,z0)
                    double mag0 = pow(pow(x0, 2.0) + pow(y0, 2.0) + pow(z0, 2.0), 0.5);

                    // Calculate unit vectors for the octahedron that intersect the new cell center
                    double Diag1X, Diag1Y, Diag1Z, Diag2X, Diag2Y, Diag2Z, Diag3X, Diag3Y, Diag3Z;
                    double Angle1 =
                        (GrainUnitVector(9 * MyOrientation) * x0 + GrainUnitVector(9 * MyOrientation + 1) * y0 +
                         GrainUnitVector(9 * MyOrientation + 2) * z0) /
                        mag0;
                    if (Angle1 < 0) {
                        Diag1X = GrainUnitVector(9 * MyOrientation);
                        Diag1Y = GrainUnitVector(9 * MyOrientation + 1);
                        Diag1Z = GrainUnitVector(9 * MyOrientation + 2);
                    }
                    else {
                        Diag1X = -GrainUnitVector(9 * MyOrientation);
                        Diag1Y = -GrainUnitVector(9 * MyOrientation + 1);
                        Diag1Z = -GrainUnitVector(9 * MyOrientation + 2);
                    }

                    double Angle2 =
                        (GrainUnitVector(9 * MyOrientation + 3) * x0 + GrainUnitVector(9 * MyOrientation + 4) * y0 +
                         GrainUnitVector(9 * MyOrientation + 5) * z0) /
                        mag0;
                    if (Angle2 < 0) {
                        Diag2X = GrainUnitVector(9 * MyOrientation + 3);
                        Diag2Y = GrainUnitVector(9 * MyOrientation + 4);
                        Diag2Z = GrainUnitVector(9 * MyOrientation + 5);
                    }
                    else {
                        Diag2X = -GrainUnitVector(9 * MyOrientation + 3);
                        Diag2Y = -GrainUnitVector(9 * MyOrientation + 4);
                        Diag2Z = -GrainUnitVector(9 * MyOrientation + 5);
                    }

                    double Angle3 =
                        (GrainUnitVector(9 * MyOrientation + 6) * x0 + GrainUnitVector(9 * MyOrientation + 7) * y0 +
                         GrainUnitVector(9 * MyOrientation + 8) * z0) /
                        mag0;
                    if (Angle3 < 0) {
                        Diag3X = GrainUnitVector(9 * MyOrientation + 6);
                        Diag3Y = GrainUnitVector(9 * MyOrientation + 7);
                        Diag3Z = GrainUnitVector(9 * MyOrientation + 8);
                    }
                    else {
                        Diag3X = -GrainUnitVector(9 * MyOrientation + 6);
                        Diag3Y = -GrainUnitVector(9 * MyOrientation + 7);
                        Diag3Z = -GrainUnitVector(9 * MyOrientation + 8);
                    }

                    double U1[3], U2[3], UU[3], Norm[3];
                    U1[0] = Diag2X - Diag1X;
                    U1[1] = Diag2Y - Diag1Y;
                    U1[2] = Diag2Z - Diag1Z;
                    U2[0] = Diag3X - Diag1X;
                    U2[1] = Diag3Y - Diag1Y;
                    U2[2] = Diag3Z - Diag1Z;
                    UU[0] = U1[1] * U2[2] - U1[2] * U2[1];
                    UU[1] = U1[2] * U2[0] - U1[0] * U2[2];
                    UU[2] = U1[0] * U2[1] - U1[1] * U2[0];
                    double NDem = sqrt(UU[0] * UU[0] + UU[1] * UU[1] + UU[2] * UU[2]);
                    Norm[0] = UU[0] / NDem;
                    Norm[1] = UU[1] / NDem;
                    Norm[2] = UU[2] / NDem;
                    // normal to capturing plane
                    double normx = Norm[0];
                    double normy = Norm[1];
                    double normz = Norm[2];
                    double ParaT =
                        (normx * x0 + normy * y0 + normz * z0) / (normx * Diag1X + normy * Diag1Y + normz * Diag1Z);
                    float CDLVal =
                        pow(pow(ParaT * Diag1X, 2.0) + pow(ParaT * Diag1Y, 2.0) + pow(ParaT * Diag1Z, 2.0), 0.5);
                    long int CDLIndex = (long int)(26) * D3D1ConvPosition + (long int)(n);
                    CritDiagonalLength(CDLIndex) = CDLVal;
                }
                CellType(GlobalD3D1ConvPosition) = Active;
            }
        },
        NucleationThisDT);
    nn += NucleationThisDT;
}

//*****************************************************************************/
// Decentered octahedron algorithm for the capture of new interface cells by grains
void CellCapture(int np, int cycle, int DecompositionStrategy, int LocalActiveDomainSize, int, int MyXSlices,
                 int MyYSlices, double AConst, double BConst, double CConst, double DConst, int MyXOffset,
                 int MyYOffset, ViewI2D ItList, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewI CritTimeStep,
                 ViewF UndercoolingCurrent, ViewF UndercoolingChange, ViewF GrainUnitVector, ViewF CritDiagonalLength,
                 ViewF DiagonalLength, ViewI GrainOrientation, ViewI CellType, ViewF DOCenter, ViewI GrainID,
                 int NGrainOrientations, Buffer2D BufferWestSend, Buffer2D BufferEastSend, Buffer2D BufferNorthSend,
                 Buffer2D BufferSouthSend, Buffer2D BufferNorthEastSend, Buffer2D BufferNorthWestSend,
                 Buffer2D BufferSouthEastSend, Buffer2D BufferSouthWestSend, int BufSizeX, int BufSizeY, int ZBound_Low,
                 int nzActive, int, int layernumber, ViewI LayerID, ViewI SteeringVector, ViewI numSteer_G,
                 ViewI_H numSteer_H) {

    // Cell capture - parallel reduce loop over all type Active cells, counting number of ghost node cells that need to
    // be accounted for
    Kokkos::parallel_for(
        "CellCapture", LocalActiveDomainSize, KOKKOS_LAMBDA(const long int &D3D1ConvPosition) {
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
                    SteeringVector(Kokkos::atomic_fetch_add(&numSteer_G(0), 1)) = D3D1ConvPosition;
                }
            }
        });
    Kokkos::deep_copy(numSteer_H, numSteer_G);

    Kokkos::parallel_for(
        "CellCapture", numSteer_H(0), KOKKOS_LAMBDA(const int &num) {
            numSteer_G(0) = 0;
            int D3D1ConvPosition = SteeringVector(num);
            // Cells of interest for the CA
            int RankZ = D3D1ConvPosition / (MyXSlices * MyYSlices);
            int Rem = D3D1ConvPosition % (MyXSlices * MyYSlices);
            int RankX = Rem / MyYSlices;
            int RankY = Rem % MyYSlices;
            int GlobalZ = RankZ + ZBound_Low;
            int GlobalD3D1ConvPosition = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
            // Update local diagonal length of active cell
            double LocU = UndercoolingCurrent(GlobalD3D1ConvPosition);
            LocU = min(210.0, LocU);
            double V = AConst * pow(LocU, 3.0) + BConst * pow(LocU, 2.0) + CConst * LocU + DConst;
            V = max(0.0, V);
            DiagonalLength(D3D1ConvPosition) += min(0.045, V); // Max amount the diagonal can grow per time step
            // Cycle through all neigboring cells on this processor to see if they have been captured
            // Cells in ghost nodes cannot capture cells on other processors
            int LCount = 0;
            // Which neighbors should be iterated over?
            int ItBounds;
            // If X and Y coordinates are not on edges, Case 0: iteratation over neighbors 0-25 possible
            // If Y coordinate is on lower edge, Case 1: iteration over only neighbors 9-25 possible
            // If Y coordinate is on upper edge, Case 2: iteration over only neighbors 0-16 possible
            // If X coordinate is on lower edge, Case 3: iteration over only neighbors
            // 0,1,3,4,6,8,9,10,11,13,15,17,18,20,21,22,24 If X coordinate is on upper edge, Case 4: iteration over only
            // neighbors 0,2,3,4,5,7,9,10,12,14,16,17,19,20,21,23,25 If X/Y coordinates are on lower edge, Case 5:
            // iteration over only neighbors 9,10,11,13,15,17,18,20,21,22,24 If X coordinate is on upper edge/Y on lower
            // edge, Case 6: If X coordinate is on lower edge/Y on upper edge, Case 7: If X/Y coordinates are on upper
            // edge, Case 8:
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
            int NListLength;
            if (ItBounds == 0) {
                NListLength = 26;
            }
            else if (ItBounds > 4) {
                NListLength = 11;
            }
            else {
                NListLength = 17;
            }
            // "ll" corresponds to the specific position on the list of neighboring cells
            for (int ll = 0; ll < NListLength; ll++) {
                // "l" correpsponds to the specific neighboring cell
                int l = ItList(ItBounds, ll);
                // Local coordinates of adjacent cell center
                int MyNeighborX = RankX + NeighborX(l);
                int MyNeighborY = RankY + NeighborY(l);
                int MyNeighborZ = RankZ + NeighborZ(l);

                if (MyNeighborZ < nzActive) {
                    long int NeighborD3D1ConvPosition =
                        MyNeighborZ * MyXSlices * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
                    long int GlobalNeighborD3D1ConvPosition =
                        (MyNeighborZ + ZBound_Low) * MyXSlices * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
                    if (CellType(GlobalNeighborD3D1ConvPosition) == Liquid)
                        LCount = 1;
                    // Capture of cell located at "NeighborD3D1ConvPosition" if this condition is satisfied
                    if ((DiagonalLength(D3D1ConvPosition) >= CritDiagonalLength(26 * D3D1ConvPosition + l)) &&
                        (CellType(GlobalNeighborD3D1ConvPosition) == Liquid)) {
                        // Use of atomic_compare_exchange
                        // (https://github.com/kokkos/kokkos/wiki/Kokkos%3A%3Aatomic_compare_exchange) old_val =
                        // atomic_compare_exchange(ptr_to_value,comparison_value, new_value); Atomicly sets the value at
                        // the address given by ptr_to_value to new_value if the current value at ptr_to_value is equal
                        // to comparison_value Returns the previously stored value at the address independent on whether
                        // the exchange has happened. If this cell's is a liquid cell, change it to "TemporaryUpdate"
                        // type and return a value of "liquid" If this cell has already been changed to
                        // "TemporaryUpdate" type, return a value of "0"
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
                            int MyOrientation = GrainOrientation(((abs(h) - 1) % NGrainOrientations));

                            // The new cell is captured by this cell's growing octahedron (Grain "h")
                            GrainID(GlobalNeighborD3D1ConvPosition) = h;
                            // (cxold, cyold, czold) are the coordiantes of this decentered octahedron
                            double cxold = DOCenter((long int)(3) * D3D1ConvPosition);
                            double cyold = DOCenter((long int)(3) * D3D1ConvPosition + (long int)(1));
                            double czold = DOCenter((long int)(3) * D3D1ConvPosition + (long int)(2));

                            // (xp,yp,zp) are the global coordinates of the new cell's center
                            double xp = GlobalX + NeighborX(l) + 0.5;
                            double yp = GlobalY + NeighborY(l) + 0.5;
                            double zp = GlobalZ + NeighborZ(l) + 0.5;

                            // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the image of
                            // the center of the new cell
                            double x0 = xp - cxold;
                            double y0 = yp - cyold;
                            double z0 = zp - czold;

                            // mag0 is the magnitude of (x0,y0,z0)
                            double mag0 = sqrtf(x0 * x0 + y0 * y0 + z0 * z0);

                            // Calculate unit vectors for the octahedron that intersect the new cell center
                            double Diag1X, Diag1Y, Diag1Z, Diag2X, Diag2Y, Diag2Z, Diag3X, Diag3Y, Diag3Z;

                            double Angle1 =
                                (GrainUnitVector(9 * MyOrientation) * x0 + GrainUnitVector(9 * MyOrientation + 1) * y0 +
                                 GrainUnitVector(9 * MyOrientation + 2) * z0) /
                                mag0;
                            double Angle2 = (GrainUnitVector(9 * MyOrientation + 3) * x0 +
                                             GrainUnitVector(9 * MyOrientation + 4) * y0 +
                                             GrainUnitVector(9 * MyOrientation + 5) * z0) /
                                            mag0;
                            double Angle3 = (GrainUnitVector(9 * MyOrientation + 6) * x0 +
                                             GrainUnitVector(9 * MyOrientation + 7) * y0 +
                                             GrainUnitVector(9 * MyOrientation + 8) * z0) /
                                            mag0;

                            Diag1X = GrainUnitVector(9 * MyOrientation) * (2 * (Angle1 < 0) - 1);
                            Diag1Y = GrainUnitVector(9 * MyOrientation + 1) * (2 * (Angle1 < 0) - 1);
                            Diag1Z = GrainUnitVector(9 * MyOrientation + 2) * (2 * (Angle1 < 0) - 1);

                            Diag2X = GrainUnitVector(9 * MyOrientation + 3) * (2 * (Angle2 < 0) - 1);
                            Diag2Y = GrainUnitVector(9 * MyOrientation + 4) * (2 * (Angle2 < 0) - 1);
                            Diag2Z = GrainUnitVector(9 * MyOrientation + 5) * (2 * (Angle2 < 0) - 1);

                            Diag3X = GrainUnitVector(9 * MyOrientation + 6) * (2 * (Angle3 < 0) - 1);
                            Diag3Y = GrainUnitVector(9 * MyOrientation + 7) * (2 * (Angle3 < 0) - 1);
                            Diag3Z = GrainUnitVector(9 * MyOrientation + 8) * (2 * (Angle3 < 0) - 1);

                            double U1[3], U2[3], UU[3], Norm[3];
                            U1[0] = Diag2X - Diag1X;
                            U1[1] = Diag2Y - Diag1Y;
                            U1[2] = Diag2Z - Diag1Z;
                            U2[0] = Diag3X - Diag1X;
                            U2[1] = Diag3Y - Diag1Y;
                            U2[2] = Diag3Z - Diag1Z;
                            UU[0] = U1[1] * U2[2] - U1[2] * U2[1];
                            UU[1] = U1[2] * U2[0] - U1[0] * U2[2];
                            UU[2] = U1[0] * U2[1] - U1[1] * U2[0];
                            double NDem = sqrtf(UU[0] * UU[0] + UU[1] * UU[1] + UU[2] * UU[2]);
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

                            // Determine which of the 3 corners of the capturing face is closest to the captured cell
                            // center
                            double DistToCorner[3];
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
                            double mindisttocorner, xc, yc, zc;

                            mindisttocorner = DistToCorner[idx];
                            xc = TriangleX[idx], yc = TriangleY[idx], zc = TriangleZ[idx];

                            double x1, y1, z1, x2, y2, z2;
                            x1 = TriangleX[(idx + 1) % 3], y1 = TriangleY[(idx + 1) % 3], z1 = TriangleZ[(idx + 1) % 3];
                            x2 = TriangleX[(idx + 2) % 3], y2 = TriangleY[(idx + 2) % 3], z2 = TriangleZ[(idx + 2) % 3];

                            double D1 =
                                sqrtf(((xp - x2) * (xp - x2)) + ((yp - y2) * (yp - y2)) + ((zp - z2) * (zp - z2)));
                            double D2 =
                                sqrtf(((xc - x2) * (xc - x2)) + ((yc - y2) * (yc - y2)) + ((zc - z2) * (zc - z2)));
                            double D3 =
                                sqrtf(((xp - x1) * (xp - x1)) + ((yp - y1) * (yp - y1)) + ((zp - z1) * (zp - z1)));
                            double D4 =
                                sqrtf(((xc - x1) * (xc - x1)) + ((yc - y1) * (yc - y1)) + ((zc - z1) * (zc - z1)));

                            double I1, I2, J1, J2;
                            I1 = 0;
                            I2 = D2;
                            J1 = 0;
                            J2 = D4;
                            // If minimum distance to corner = 0, the octahedron corner captured the new cell center
                            if (mindisttocorner != 0) {
                                I1 = D1 * ((xp - x2) * (xc - x2) + (yp - y2) * (yc - y2) + (zp - z2) * (zc - z2)) /
                                     (D1 * D2);
                                I2 = D2 - I1;
                                J1 = D3 * ((xp - x1) * (xc - x1) + (yp - y1) * (yc - y1) + (zp - z1) * (zc - z1)) /
                                     (D3 * D4);
                                J2 = D4 - J1;
                            }
                            double L12 = 0.5 * (min(I1, sqrt(3.0)) + min(I2, sqrt(3.0)));
                            double L13 = 0.5 * (min(J1, sqrt(3.0)) + min(J2, sqrt(3.0)));
                            double NewODiagL = sqrt(2.0) * max(L12, L13); // half diagonal length of new octahedron

                            DiagonalLength(NeighborD3D1ConvPosition) = NewODiagL;
                            // Calculate coordinates of new decentered octahedron center
                            double CaptDiag[3], CaptDiagUV[3];
                            CaptDiag[0] = xc - cxold;
                            CaptDiag[1] = yc - cyold;
                            CaptDiag[2] = zc - czold;

                            NDem =
                                sqrt(CaptDiag[0] * CaptDiag[0] + CaptDiag[1] * CaptDiag[1] + CaptDiag[2] * CaptDiag[2]);
                            CaptDiagUV[0] = CaptDiag[0] / NDem;
                            CaptDiagUV[1] = CaptDiag[1] / NDem;
                            CaptDiagUV[2] = CaptDiag[2] / NDem;
                            // (cx, cy, cz) are the coordiantes of the new active cell's decentered octahedron
                            double cx = xc - NewODiagL * CaptDiagUV[0];
                            double cy = yc - NewODiagL * CaptDiagUV[1];
                            double cz = zc - NewODiagL * CaptDiagUV[2];

                            DOCenter((long int)(3) * NeighborD3D1ConvPosition) = cx;
                            DOCenter((long int)(3) * NeighborD3D1ConvPosition + (long int)(1)) = cy;
                            DOCenter((long int)(3) * NeighborD3D1ConvPosition + (long int)(2)) = cz;

                            // Calculate critical octahedron diagonal length to activate nearest neighbor.
                            // First, calculate the unique planes (4) associated with all octahedron faces (8)
                            // Then just look at distance between face and the point of interest (cell center of
                            // neighbor). The critical diagonal length will be the maximum of these (since all other
                            // planes will have passed over the point by then
                            // ... meaning it must be in the octahedron)

                            double Fx[4], Fy[4], Fz[4], D[4], Dfabs;

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
                                double x0 = xp + NeighborX[n] - cx;
                                double y0 = yp + NeighborY[n] - cy;
                                double z0 = zp + NeighborZ[n] - cz;
                                D[0] = x0 * Fx[0] + y0 * Fy[0] + z0 * Fz[0];
                                D[1] = x0 * Fx[1] + y0 * Fy[1] + z0 * Fz[1];
                                D[2] = x0 * Fx[2] + y0 * Fy[2] + z0 * Fz[2];
                                D[3] = x0 * Fx[3] + y0 * Fy[3] + z0 * Fz[3];
                                Dfabs = max(max(fabs(D[0]), fabs(D[1])), max(fabs(D[2]), fabs(D[3])));
                                CritDiagonalLength((long int)(26) * NeighborD3D1ConvPosition + (long int)(n)) = Dfabs;
                            }

                            if (np > 1) {

                                float GhostGID = (float)(GrainID(GlobalNeighborD3D1ConvPosition));
                                float GhostDOCX = DOCenter((long int)(3) * NeighborD3D1ConvPosition);
                                float GhostDOCY = DOCenter((long int)(3) * NeighborD3D1ConvPosition + (long int)(1));
                                float GhostDOCZ = DOCenter((long int)(3) * NeighborD3D1ConvPosition + (long int)(2));
                                float GhostDL = DiagonalLength(NeighborD3D1ConvPosition);
                                // Collect data for the ghost nodes:
                                if (DecompositionStrategy == 1) {
                                    int GNPosition = MyNeighborZ * BufSizeX + MyNeighborX;
                                    if (MyNeighborY == 1) {
                                        BufferSouthSend(GNPosition, 0) = GhostGID;
                                        BufferSouthSend(GNPosition, 1) = GhostDOCX;
                                        BufferSouthSend(GNPosition, 2) = GhostDOCY;
                                        BufferSouthSend(GNPosition, 3) = GhostDOCZ;
                                        BufferSouthSend(GNPosition, 4) = GhostDL;
                                    }
                                    else if (MyNeighborY == MyYSlices - 2) {
                                        int GNPosition = MyNeighborZ * BufSizeX + MyNeighborX;
                                        BufferNorthSend(GNPosition, 0) = GhostGID;
                                        BufferNorthSend(GNPosition, 1) = GhostDOCX;
                                        BufferNorthSend(GNPosition, 2) = GhostDOCY;
                                        BufferNorthSend(GNPosition, 3) = GhostDOCZ;
                                        BufferNorthSend(GNPosition, 4) = GhostDL;
                                    }
                                }
                                else {
                                    if (MyNeighborY == 1) {
                                        if (MyNeighborX == MyXSlices - 2) {
                                            int GNPosition = MyNeighborZ * BufSizeX + MyNeighborX - 1;
                                            BufferSouthSend(GNPosition, 0) = GhostGID;
                                            BufferSouthSend(GNPosition, 1) = GhostDOCX;
                                            BufferSouthSend(GNPosition, 2) = GhostDOCY;
                                            BufferSouthSend(GNPosition, 3) = GhostDOCZ;
                                            BufferSouthSend(GNPosition, 4) = GhostDL;
                                            GNPosition = MyNeighborZ * BufSizeY + MyNeighborY - 1;
                                            BufferEastSend(GNPosition, 0) = GhostGID;
                                            BufferEastSend(GNPosition, 1) = GhostDOCX;
                                            BufferEastSend(GNPosition, 2) = GhostDOCY;
                                            BufferEastSend(GNPosition, 3) = GhostDOCZ;
                                            BufferEastSend(GNPosition, 4) = GhostDL;
                                            GNPosition = MyNeighborZ;
                                            BufferSouthEastSend(GNPosition, 0) = GhostGID;
                                            BufferSouthEastSend(GNPosition, 1) = GhostDOCX;
                                            BufferSouthEastSend(GNPosition, 2) = GhostDOCY;
                                            BufferSouthEastSend(GNPosition, 3) = GhostDOCZ;
                                            BufferSouthEastSend(GNPosition, 4) = GhostDL;
                                        }
                                        else if (MyNeighborX == 1) {
                                            int GNPosition = MyNeighborZ * BufSizeX + MyNeighborX - 1;
                                            BufferSouthSend(GNPosition, 0) = GhostGID;
                                            BufferSouthSend(GNPosition, 1) = GhostDOCX;
                                            BufferSouthSend(GNPosition, 2) = GhostDOCY;
                                            BufferSouthSend(GNPosition, 3) = GhostDOCZ;
                                            BufferSouthSend(GNPosition, 4) = GhostDL;
                                            GNPosition = MyNeighborZ * BufSizeY + MyNeighborY - 1;
                                            BufferWestSend(GNPosition, 0) = GhostGID;
                                            BufferWestSend(GNPosition, 1) = GhostDOCX;
                                            BufferWestSend(GNPosition, 2) = GhostDOCY;
                                            BufferWestSend(GNPosition, 3) = GhostDOCZ;
                                            BufferWestSend(GNPosition, 4) = GhostDL;
                                            GNPosition = MyNeighborZ;
                                            BufferSouthWestSend(GNPosition, 0) = GhostGID;
                                            BufferSouthWestSend(GNPosition, 1) = GhostDOCX;
                                            BufferSouthWestSend(GNPosition, 2) = GhostDOCY;
                                            BufferSouthWestSend(GNPosition, 3) = GhostDOCZ;
                                            BufferSouthWestSend(GNPosition, 4) = GhostDL;
                                        }
                                        else if ((MyNeighborX > 1) && (MyNeighborX < MyXSlices - 2)) {
                                            // This is being sent to MyLeft
                                            int GNPosition = MyNeighborZ * BufSizeX + MyNeighborX - 1;
                                            BufferSouthSend(GNPosition, 0) = GhostGID;
                                            BufferSouthSend(GNPosition, 1) = GhostDOCX;
                                            BufferSouthSend(GNPosition, 2) = GhostDOCY;
                                            BufferSouthSend(GNPosition, 3) = GhostDOCZ;
                                            BufferSouthSend(GNPosition, 4) = GhostDL;
                                        }
                                    }
                                    else if (MyNeighborY == MyYSlices - 2) {
                                        // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                                        if (MyNeighborX == MyXSlices - 2) {
                                            int GNPosition = MyNeighborZ * BufSizeX + MyNeighborX - 1;
                                            BufferNorthSend(GNPosition, 0) = GhostGID;
                                            BufferNorthSend(GNPosition, 1) = GhostDOCX;
                                            BufferNorthSend(GNPosition, 2) = GhostDOCY;
                                            BufferNorthSend(GNPosition, 3) = GhostDOCZ;
                                            BufferNorthSend(GNPosition, 4) = GhostDL;
                                            GNPosition = MyNeighborZ * BufSizeY + MyNeighborY - 1;
                                            BufferEastSend(GNPosition, 0) = GhostGID;
                                            BufferEastSend(GNPosition, 1) = GhostDOCX;
                                            BufferEastSend(GNPosition, 2) = GhostDOCY;
                                            BufferEastSend(GNPosition, 3) = GhostDOCZ;
                                            BufferEastSend(GNPosition, 4) = GhostDL;
                                            GNPosition = MyNeighborZ;
                                            BufferNorthEastSend(GNPosition, 0) = GhostGID;
                                            BufferNorthEastSend(GNPosition, 1) = GhostDOCX;
                                            BufferNorthEastSend(GNPosition, 2) = GhostDOCY;
                                            BufferNorthEastSend(GNPosition, 3) = GhostDOCZ;
                                            BufferNorthEastSend(GNPosition, 4) = GhostDL;
                                        }
                                        else if (MyNeighborX == 1) {
                                            int GNPosition = MyNeighborZ * BufSizeX + MyNeighborX - 1;
                                            BufferNorthSend(GNPosition, 0) = GhostGID;
                                            BufferNorthSend(GNPosition, 1) = GhostDOCX;
                                            BufferNorthSend(GNPosition, 2) = GhostDOCY;
                                            BufferNorthSend(GNPosition, 3) = GhostDOCZ;
                                            BufferNorthSend(GNPosition, 4) = GhostDL;
                                            GNPosition = MyNeighborZ * BufSizeY + MyNeighborY - 1;
                                            BufferWestSend(GNPosition, 0) = GhostGID;
                                            BufferWestSend(GNPosition, 1) = GhostDOCX;
                                            BufferWestSend(GNPosition, 2) = GhostDOCY;
                                            BufferWestSend(GNPosition, 3) = GhostDOCZ;
                                            BufferWestSend(GNPosition, 4) = GhostDL;
                                            GNPosition = MyNeighborZ;
                                            BufferNorthWestSend(GNPosition, 0) = GhostGID;
                                            BufferNorthWestSend(GNPosition, 1) = GhostDOCX;
                                            BufferNorthWestSend(GNPosition, 2) = GhostDOCY;
                                            BufferNorthWestSend(GNPosition, 3) = GhostDOCZ;
                                            BufferNorthWestSend(GNPosition, 4) = GhostDL;
                                        }
                                        else if ((MyNeighborX > 1) && (MyNeighborX < MyXSlices - 2)) {
                                            int GNPosition = MyNeighborZ * BufSizeX + MyNeighborX - 1;
                                            BufferNorthSend(GNPosition, 0) = GhostGID;
                                            BufferNorthSend(GNPosition, 1) = GhostDOCX;
                                            BufferNorthSend(GNPosition, 2) = GhostDOCY;
                                            BufferNorthSend(GNPosition, 3) = GhostDOCZ;
                                            BufferNorthSend(GNPosition, 4) = GhostDL;
                                        }
                                    }
                                    else if ((MyNeighborX == 1) && (MyNeighborY > 1) && (MyNeighborY < MyYSlices - 2)) {
                                        int GNPosition = MyNeighborZ * BufSizeY + MyNeighborY - 1;
                                        BufferWestSend(GNPosition, 0) = GhostGID;
                                        BufferWestSend(GNPosition, 1) = GhostDOCX;
                                        BufferWestSend(GNPosition, 2) = GhostDOCY;
                                        BufferWestSend(GNPosition, 3) = GhostDOCZ;
                                        BufferWestSend(GNPosition, 4) = GhostDL;
                                    }
                                    else if ((MyNeighborX == MyXSlices - 2) && (MyNeighborY > 1) &&
                                             (MyNeighborY < MyYSlices - 2)) {
                                        int GNPosition = MyNeighborZ * BufSizeY + MyNeighborY - 1;
                                        BufferEastSend(GNPosition, 0) = GhostGID;
                                        BufferEastSend(GNPosition, 1) = GhostDOCX;
                                        BufferEastSend(GNPosition, 2) = GhostDOCY;
                                        BufferEastSend(GNPosition, 3) = GhostDOCZ;
                                        BufferEastSend(GNPosition, 4) = GhostDL;
                                    }
                                } // End if statement for ghost node marking
                            }     // End if statement for serial/parallel code
                            // Only update the new cell's type once Critical Diagonal Length, Triangle Index, and
                            // Diagonal Length values have been assigned to it Avoids the race condition in which the
                            // new cell is activated, and another thread acts on the new active cell before the cell's
                            // new critical diagonal length/triangle index/diagonal length values are assigned
                            CellType(GlobalNeighborD3D1ConvPosition) = Active;
                        } // End if statement within locked capture loop
                    }     // End if statement for outer capture loop
                }         // End if statement over neighbors on the active grid
            }             // End loop over all neighbors of this active cell
            if (LCount == 0) {
                // This active cell has no more neighboring cells to be captured, becomes solid
                CellType(GlobalD3D1ConvPosition) = Solid;
            }
        });
    Kokkos::fence();
}

//*****************************************************************************/
// Prints intermediate code output to stdout, checks to see if solidification is complete
void IntermediateOutputAndCheck(int id, int &cycle, int MyXSlices, int MyYSlices, int LocalDomainSize,
                                int LocalActiveDomainSize, int nn, int &XSwitch, ViewI CellType, ViewI CritTimeStep,
                                std::string TemperatureDataType, int *FinishTimeStep, int layernumber, int,
                                int ZBound_Low, ViewI LayerID) {

    sample::ValueType CellTypeStorage;
    Kokkos::parallel_reduce(
        LocalDomainSize,
        KOKKOS_LAMBDA(const int &D3D1ConvPosition, sample::ValueType &upd) {
            if (LayerID(D3D1ConvPosition) == layernumber) {
                if (CellType(D3D1ConvPosition) == Liquid) {
                    if (CritTimeStep(D3D1ConvPosition) > cycle)
                        upd.the_array[0] += 1;
                    else
                        upd.the_array[1] += 1;
                }
                else if (CellType(D3D1ConvPosition) == Active)
                    upd.the_array[2] += 1;
                else if (CellType(D3D1ConvPosition) == Solid)
                    upd.the_array[3] += 1;
            }
            else {
                if (CellType(D3D1ConvPosition) == Liquid)
                    upd.the_array[4] += 1;
            }
        },
        Kokkos::Sum<sample::ValueType>(CellTypeStorage));
    unsigned long int Global_nn = 0;
    unsigned long int LocalSuperheatedCells = CellTypeStorage.the_array[0];
    unsigned long int LocalUndercooledCells = CellTypeStorage.the_array[1];
    unsigned long int LocalActiveCells = CellTypeStorage.the_array[2];
    unsigned long int LocalSolidCells = CellTypeStorage.the_array[3];
    unsigned long int LocalRemainingLiquidCells = CellTypeStorage.the_array[4];

    unsigned long int GlobalSuperheatedCells, GlobalUndercooledCells, GlobalActiveCells, GlobalSolidCells,
        GlobalRemainingLiquidCells;
    MPI_Reduce(&LocalSuperheatedCells, &GlobalSuperheatedCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalUndercooledCells, &GlobalUndercooledCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalActiveCells, &GlobalActiveCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalSolidCells, &GlobalSolidCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalRemainingLiquidCells, &GlobalRemainingLiquidCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&nn, &Global_nn, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (id == 0) {
        std::cout << "cycle = " << cycle << " Superheated liquid cells = " << GlobalSuperheatedCells
                  << " Undercooled liquid cells = " << GlobalUndercooledCells
                  << " Number of nucleation events this layer " << Global_nn
                  << " Remaining liquid cells in future layers of this simulation = " << GlobalRemainingLiquidCells
                  << std::endl;
        if (GlobalSuperheatedCells + GlobalUndercooledCells == 0)
            XSwitch = 1;
    }
    MPI_Bcast(&XSwitch, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if ((XSwitch == 0) && (TemperatureDataType == "R")) {
        MPI_Bcast(&GlobalUndercooledCells, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        if (GlobalUndercooledCells == 0) {
            // Check when the next superheated cells go below the liquidus
            unsigned long int NextCTS;
            Kokkos::parallel_reduce(
                "CellCapture", LocalActiveDomainSize,
                KOKKOS_LAMBDA(const int &D3D1ConvPosition, unsigned long int &tempv) {
                    int RankZ = D3D1ConvPosition / (MyXSlices * MyYSlices);
                    int Rem = D3D1ConvPosition % (MyXSlices * MyYSlices);
                    int RankX = Rem / MyYSlices;
                    int RankY = Rem % MyYSlices;
                    int GlobalZ = RankZ + ZBound_Low;
                    int GlobalD3D1ConvPosition = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                    unsigned long int CritTimeStep_ThisCell = (unsigned long int)(CritTimeStep(GlobalD3D1ConvPosition));
                    if ((CellType(GlobalD3D1ConvPosition) == Liquid) &&
                        (LayerID(GlobalD3D1ConvPosition) == layernumber)) {
                        if (CritTimeStep_ThisCell < tempv)
                            tempv = CritTimeStep_ThisCell;
                    }
                },
                Kokkos::Min<unsigned long int>(NextCTS));

            unsigned long int GlobalNextCTS;
            MPI_Allreduce(&NextCTS, &GlobalNextCTS, 1, MPI_UNSIGNED_LONG, MPI_MIN, MPI_COMM_WORLD);
            if ((GlobalNextCTS - cycle) > 5000) {
                // Jump to next time step when solidification starts again
                cycle = GlobalNextCTS - 1;
                if (id == 0)
                    std::cout << "Jumping to cycle " << cycle + 1 << std::endl;
            }
            if (cycle >= FinishTimeStep[layernumber])
                XSwitch = 1;
        }
    }
}
