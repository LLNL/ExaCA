// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAupdate.hpp"
#include "CAprint.hpp"

#include "mpi.h"

#include <cmath>

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
    KOKKOS_INLINE_FUNCTION // Assignment operator
        array_type &
        operator=(const array_type &src) {
        for (int i = 0; i < N; i++) {
            the_array[i] = src.the_array[i];
        }
        return *this;
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
//void Nucleation(int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int cycle, int &nn,
//                ViewI CellType, ViewI NucleiLocations, ViewI NucleationTimes, ViewI NucleiGrainID, ViewI GrainID,
//                ViewI GrainOrientation, ViewF DOCenter, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ,
//                ViewF GrainUnitVector, int NGrainOrientations, double AConst, double BConst, double CConst, double DConst,
//                int PossibleNuclei_ThisRank, int ZBound_Low, int layernumber, ViewI LayerID, ViewI CritTimeStep, ViewF UndercoolingChange, ViewF DiagonalLength, ViewF CritDiagonalLength) {
void Nucleation(int MyXSlices, int MyYSlices, int cycle, int &nn,
                    ViewI CellType, ViewI NucleiLocations, ViewI NucleationTimes, ViewI NucleiGrainID, ViewI GrainID,
                    int PossibleNuclei_ThisRank) {
    // Loop through local list of nucleation events - has the time step exceeded the time step for nucleation at the
    // sites?
    int NucleationThisDT = 0;
    Kokkos::parallel_reduce(
        "NucleiUpdateLoop", PossibleNuclei_ThisRank,
        KOKKOS_LAMBDA(const int &NucCounter, int &update) {
            int TimeStepCondition = cycle == NucleationTimes(NucCounter);
            int CellTypeCondition = CellType(NucleiLocations(NucCounter)) == Liquid;
            if (TimeStepCondition) {
                // This nucleation event is either about to happen - or the nucleation site is extinct
                // If remelting is considered - set "NucleationTimes" to this nucleation event to a negative number
                // if the site is extinct - this is to avoid a race condition in case this site solidifies
                // multiple times and is the site of multiple nucleation events
                if (!(CellTypeCondition))
                    NucleationTimes(NucCounter) = -1;
                else {
                    // (X,Y,Z) coordinates of nucleation event, on active cell grid (RankX,RankY,RankZ) and global grid
                    // (RankX,RankY,GlobalZ)
                    // Global grid used for cell type/nuclei location if no remelting, local active (this layer only)
                    // grid used if remelting Grain ID is stored inside the NucleiGrainID view in simulations with
                    // remelting (in case multiple nucleation events are possible in the same cell For simulations
                    // without remelting, the GrainID view was initialized with the new Grain ID
                    //int RankX, RankY, RankZ, Rem, GlobalZ, D3D1ConvPosition, GlobalD3D1ConvPosition, MyGrainID;
                    int GlobalD3D1ConvPosition = NucleiLocations(NucCounter);
                    CellType(GlobalD3D1ConvPosition) = TemporaryUpdate;
                    //int GlobalZ = GlobalD3D1ConvPosition / (MyXSlices * MyYSlices);
                    int Rem = GlobalD3D1ConvPosition % (MyXSlices * MyYSlices);
                    int RankX = Rem / MyYSlices;
                    int RankY = Rem % MyYSlices;
//                    RankZ = GlobalZ - ZBound_Low;
//                    D3D1ConvPosition = RankZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
// //                   if (Remelting) {
//                        MyGrainID = NucleiGrainID(NucCounter);
                        GrainID(GlobalD3D1ConvPosition) = NucleiGrainID(NucCounter);
////                    }
////                    else {
////                        MyGrainID = GrainID(GlobalD3D1ConvPosition);
////                    }
                    // This undercooled liquid cell is now a nuclei (add to count if it isn't in the ghost nodes, to
                    // avoid double counting)
                    if ((RankX > 0) && (RankX < MyXSlices - 1) && (RankY > 0) && (RankY < MyYSlices - 1))
                        update++;
//                    int GlobalX = RankX + MyXOffset;
//                    int GlobalY = RankY + MyYOffset;
//                    // DiagonalLength(D3D1ConvPosition) = 0.01;
//                    long int DOX = (long int)(3) * D3D1ConvPosition;
//                    long int DOY = (long int)(3) * D3D1ConvPosition + (long int)(1);
//                    long int DOZ = (long int)(3) * D3D1ConvPosition + (long int)(2);
//                    DOCenter(DOX) = GlobalX + 0.5;
//                    DOCenter(DOY) = GlobalY + 0.5;
//                    DOCenter(DOZ) = GlobalZ + 0.5;
//                    // The orientation for the new grain will depend on its Grain ID (nucleated grains have negative
//                    // GrainID values)
//                    int MyOrientation = GrainOrientation(((abs(MyGrainID) - 1) % NGrainOrientations));
//                    // Calculate critical values at which this active cell leads to the activation of a neighboring
//                    // liquid cell
//                    for (int n = 0; n < 26; n++) {
//
//                        // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the image of the
//                        // center of a neighbor cell
//                        double x0 = NeighborX(n);
//                        double y0 = NeighborY(n);
//                        double z0 = NeighborZ(n);
//
//                        // mag0 is the magnitude of (x0,y0,z0)
//                        double mag0 = pow(pow(x0, 2.0) + pow(y0, 2.0) + pow(z0, 2.0), 0.5);
//
//                        // Calculate unit vectors for the octahedron that intersect the new cell center
//                        double Diag1X, Diag1Y, Diag1Z, Diag2X, Diag2Y, Diag2Z, Diag3X, Diag3Y, Diag3Z;
//                        double Angle1 =
//                            (GrainUnitVector(9 * MyOrientation) * x0 + GrainUnitVector(9 * MyOrientation + 1) * y0 +
//                             GrainUnitVector(9 * MyOrientation + 2) * z0) /
//                            mag0;
//                        if (Angle1 < 0) {
//                            Diag1X = GrainUnitVector(9 * MyOrientation);
//                            Diag1Y = GrainUnitVector(9 * MyOrientation + 1);
//                            Diag1Z = GrainUnitVector(9 * MyOrientation + 2);
//                        }
//                        else {
//                            Diag1X = -GrainUnitVector(9 * MyOrientation);
//                            Diag1Y = -GrainUnitVector(9 * MyOrientation + 1);
//                            Diag1Z = -GrainUnitVector(9 * MyOrientation + 2);
//                        }
//
//                        double Angle2 =
//                            (GrainUnitVector(9 * MyOrientation + 3) * x0 + GrainUnitVector(9 * MyOrientation + 4) * y0 +
//                             GrainUnitVector(9 * MyOrientation + 5) * z0) /
//                            mag0;
//                        if (Angle2 < 0) {
//                            Diag2X = GrainUnitVector(9 * MyOrientation + 3);
//                            Diag2Y = GrainUnitVector(9 * MyOrientation + 4);
//                            Diag2Z = GrainUnitVector(9 * MyOrientation + 5);
//                        }
//                        else {
//                            Diag2X = -GrainUnitVector(9 * MyOrientation + 3);
//                            Diag2Y = -GrainUnitVector(9 * MyOrientation + 4);
//                            Diag2Z = -GrainUnitVector(9 * MyOrientation + 5);
//                        }
//
//                        double Angle3 =
//                            (GrainUnitVector(9 * MyOrientation + 6) * x0 + GrainUnitVector(9 * MyOrientation + 7) * y0 +
//                             GrainUnitVector(9 * MyOrientation + 8) * z0) /
//                            mag0;
//                        if (Angle3 < 0) {
//                            Diag3X = GrainUnitVector(9 * MyOrientation + 6);
//                            Diag3Y = GrainUnitVector(9 * MyOrientation + 7);
//                            Diag3Z = GrainUnitVector(9 * MyOrientation + 8);
//                        }
//                        else {
//                            Diag3X = -GrainUnitVector(9 * MyOrientation + 6);
//                            Diag3Y = -GrainUnitVector(9 * MyOrientation + 7);
//                            Diag3Z = -GrainUnitVector(9 * MyOrientation + 8);
//                        }
//
//                        double U1[3], U2[3], UU[3], Norm[3];
//                        U1[0] = Diag2X - Diag1X;
//                        U1[1] = Diag2Y - Diag1Y;
//                        U1[2] = Diag2Z - Diag1Z;
//                        U2[0] = Diag3X - Diag1X;
//                        U2[1] = Diag3Y - Diag1Y;
//                        U2[2] = Diag3Z - Diag1Z;
//                        UU[0] = U1[1] * U2[2] - U1[2] * U2[1];
//                        UU[1] = U1[2] * U2[0] - U1[0] * U2[2];
//                        UU[2] = U1[0] * U2[1] - U1[1] * U2[0];
//                        double NDem = sqrt(UU[0] * UU[0] + UU[1] * UU[1] + UU[2] * UU[2]);
//                        Norm[0] = UU[0] / NDem;
//                        Norm[1] = UU[1] / NDem;
//                        Norm[2] = UU[2] / NDem;
//                        // normal to capturing plane
//                        double normx = Norm[0];
//                        double normy = Norm[1];
//                        double normz = Norm[2];
//                        double ParaT =
//                            (normx * x0 + normy * y0 + normz * z0) / (normx * Diag1X + normy * Diag1Y + normz * Diag1Z);
//                        float CDL =
//                            pow(pow(ParaT * Diag1X, 2.0) + pow(ParaT * Diag1Y, 2.0) + pow(ParaT * Diag1Z, 2.0), 0.5);
//                        // Use critical diagonal length value, current undercooling, rate of undercooling, and distance to each neighbor to determine the time step at which each neighboring cell would be captured by this one
//                        int future_cycle = max(cycle,CritTimeStep(GlobalD3D1ConvPosition));
//                        float DL = 0.01; // initial diagonal length for this cell
//                        // Undercooling was 0 at the critical time step, and was incremented by UndercoolingChange every time step since
//                        double Undercooling = (cycle - CritTimeStep(GlobalD3D1ConvPosition)) * UndercoolingChange(GlobalD3D1ConvPosition);
//                        bool Inc_DL = true;
//                        while (Inc_DL) {
//                            if (DL >= CDL) {
//                                CaptureTimeStep(26 * D3D1ConvPosition + n) = future_cycle;
//                                Inc_DL = false;
//                            }
//                            else {
//                                future_cycle++;
//                                // Update local undercooling based on local cooling rate
//                                Undercooling += UndercoolingChange(GlobalD3D1ConvPosition);
//                                Undercooling = min(210.0, Undercooling);
//                                double V = AConst * pow(Undercooling, 3.0) + BConst * pow(Undercooling, 2.0) + CConst * Undercooling + DConst;
//                                V = max(0.0, V);
//                                DL += min(0.045, V); // Max amount the diagonal can grow per time step
//                            }
//                        }
//                    }
//                    CellType(GlobalD3D1ConvPosition) = Active;
                }
            }
        },
        NucleationThisDT);
    nn += NucleationThisDT;
    Kokkos::fence();
}

void CellCapture_RM(int cycle, int LocalActiveDomainSize, int,
                 int MyXSlices, int MyYSlices, int nx, int ny, double AConst, double BConst, double CConst, double DConst,
                 int MyXOffset, int MyYOffset, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ,
                 ViewI CritTimeStep, ViewF UndercoolingChange, ViewF GrainUnitVector, ViewI GrainOrientation, ViewI CellType,
                 ViewF DOCenter, ViewI GrainID, int NGrainOrientations, int ZBound_Low, int nzActive, int, ViewI SteeringVector,
                 ViewI numSteer_G, ViewI_H numSteer_H, ViewI MeltTimeStep, ViewI SolidificationEventCounter,
                 ViewI NumberOfSolidificationEvents, ViewF3D LayerTimeTempHistory, ViewF DiagonalLength, ViewF CritDiagonalLength) {

// Perform melting events, maintaining active cells at the solid-liquid interface
// For all active and liquid cells that have undercooled, record their positions in the steering vector
Kokkos::parallel_for(
    "preCellCapture_RM", LocalActiveDomainSize, KOKKOS_LAMBDA(const long int &D3D1ConvPosition) {
        
        // Coordinate of this cell on the "global" (all cells in the Z direction) grid
        int RankZ = D3D1ConvPosition / (MyXSlices * MyYSlices);
        int Rem = D3D1ConvPosition % (MyXSlices * MyYSlices);
        int RankX = Rem / MyYSlices;
        int RankY = Rem % MyYSlices;
        int GlobalZ = RankZ + ZBound_Low;
        int GlobalD3D1ConvPosition = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;

        int cellType = CellType(GlobalD3D1ConvPosition);
        bool atMeltTime = (cycle == MeltTimeStep(GlobalD3D1ConvPosition));
        bool atpastCritTime = (cycle >= CritTimeStep(GlobalD3D1ConvPosition));
        
        if ((atMeltTime) && (cellType != Liquid)) {
            // This cell should be a liquid cell
            CellType(GlobalD3D1ConvPosition) = Liquid;
        }
        else if (atpastCritTime) {
            if (cellType == Active)  {
                // Add to steering vector
                SteeringVector(Kokkos::atomic_fetch_add(&numSteer_G(0), 1)) = D3D1ConvPosition;
            }
            else if (cellType == Liquid) {
                // If this cell borders at least one solid/tempsolid cell, it should become active
                for (int l = 0; l < 26; l++) {
                    // "l" correpsponds to the specific neighboring cell
                    // Local coordinates of adjacent cell center
                    int MyNeighborX = RankX + NeighborX(l);
                    int MyNeighborY = RankY + NeighborY(l);
                    int MyNeighborZ = RankZ + NeighborZ(l);
                    bool InBoundsThisRank = (MyNeighborX >= 0) * (MyNeighborX < MyXSlices) * (MyNeighborY >= 0) * (MyNeighborY < MyYSlices) * (MyNeighborZ >= 0) * (MyNeighborZ < nzActive);
                    if (InBoundsThisRank) {
                        //int NeighborD3D1ConvPosition = MyNeighborZ * MyXSlices * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
                        int GlobalNeighborD3D1ConvPosition = (MyNeighborZ + ZBound_Low) * MyXSlices * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
                        if ((CellType(GlobalNeighborD3D1ConvPosition) == TempSolid) || (CellType(GlobalNeighborD3D1ConvPosition) == Solid) || (RankZ == 0)) {
                            // Cell activation to be performed as part of steering vector
                            l = 26;
                            SteeringVector(Kokkos::atomic_fetch_add(&numSteer_G(0), 1)) = D3D1ConvPosition;
                            CellType(GlobalD3D1ConvPosition) = TemporaryUpdate; // this cell cannot be captured - is being activated
                        }
                    }
                }
            }
            else if (cellType == TemporaryUpdate) {
                // This was a nucleation event - add to steering vector for active cell calculations
                SteeringVector(Kokkos::atomic_fetch_add(&numSteer_G(0), 1)) = D3D1ConvPosition;
            }
        }
    });
    Kokkos::fence();

    // Copy size of steering vector (containing positions of undercooled liquid/active cells) to the host
    Kokkos::deep_copy(numSteer_H, numSteer_G);

    // Iterate over all cells in the steering vector
    // Two kinds of cells: active type cells that are checking to potentially capture any liquid neighbors, and cells marked with TemporaryUpdate type to become new active cells
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

            if (CellType(GlobalD3D1ConvPosition) == TemporaryUpdate) {
                // This liquid cell is becoming an active cell via cell activation
                int MyGrainID = GrainID(GlobalD3D1ConvPosition);
                int GlobalX = RankX + MyXOffset;
                int GlobalY = RankY + MyYOffset;
                DiagonalLength(D3D1ConvPosition) = 0.01;
                long int DOX = (long int)(3) * D3D1ConvPosition;
                long int DOY = (long int)(3) * D3D1ConvPosition + (long int)(1);
                long int DOZ = (long int)(3) * D3D1ConvPosition + (long int)(2);
                DOCenter(DOX) = GlobalX + 0.5;
                DOCenter(DOY) = GlobalY + 0.5;
                DOCenter(DOZ) = GlobalZ + 0.5;
                // The orientation for the new grain will depend on its Grain ID (nucleated grains have negative
                // GrainID values)
                int MyOrientation = GrainOrientation(((abs(MyGrainID) - 1) % NGrainOrientations));
                // Calculate critical values at which this active cell leads to the activation of a neighboring
                // liquid cell
                for (int n = 0; n < 26; n++) {

                    // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the image of the
                    // center of a neighbor cell
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
                    float CDL =
                        pow(pow(ParaT * Diag1X, 2.0) + pow(ParaT * Diag1Y, 2.0) + pow(ParaT * Diag1Z, 2.0), 0.5);
                    CritDiagonalLength(26 * D3D1ConvPosition + n) = CDL;
                }
                CellType(GlobalD3D1ConvPosition) = Active;
            }
            else  {
                // This cell is active - check neighbors to see if they can be captured
                // Update local diagonal length of active cell
                double LocU = (cycle - CritTimeStep(GlobalD3D1ConvPosition)) * UndercoolingChange(GlobalD3D1ConvPosition);
                LocU = min(210.0, LocU);
                double V = AConst * pow(LocU, 3.0) + BConst * pow(LocU, 2.0) + CConst * LocU + DConst;
                V = max(0.0, V);
                DiagonalLength(D3D1ConvPosition) += min(0.045, V); // Max amount the diagonal can grow per time step
                // Cycle through all neigboring cells on this processor to see if they have been captured
                bool DeactivateCell = true;
                for (int l = 0; l < 26; l++) {
                    // Local coordinates of adjacent cell center
                    int MyNeighborX = RankX + NeighborX(l);
                    int MyNeighborY = RankY + NeighborY(l);
                    int MyNeighborZ = RankZ + NeighborZ(l);

                    if ((MyNeighborZ < nzActive) && (MyNeighborZ >= 0) && (MyNeighborX >= 0) && (MyNeighborX < MyXSlices) && (MyNeighborY >= 0) && (MyNeighborY < MyYSlices)) {
                        int NeighborD3D1ConvPosition =
                            MyNeighborZ * MyXSlices * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
                        int GlobalNeighborD3D1ConvPosition =
                            (MyNeighborZ + ZBound_Low) * MyXSlices * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
                        if (CellType(GlobalNeighborD3D1ConvPosition) == Liquid)
                            DeactivateCell = false;
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
                                CellType(GlobalNeighborD3D1ConvPosition) = Active;
                            } // end if statement for locked capture loop
                        } // end if statement for capture loop
                    } // end if statement for checking if neighboring cells are in bounds
                } // end loop over neighboring cells
                // Deactivate cell only if its not in the ghost nodes (RankX/RankY = 0 or MyXSlices-1/MyYSlices-1... with the exception of being at the domain boundary (global X = 0 or nx and global Y = 0 or ny)
                if ((DeactivateCell) && ((RankX != 0) || (MyXOffset == 0)) && ((RankX != MyXSlices-1) || (MyXSlices+MyXOffset == nx)) && ((RankY != 0) || (MyYOffset == 0)) && ((RankY != MyYSlices-1) || (MyYSlices+MyYOffset == ny))) {
                    // This cell is no longer active - all liquid neighbors have been captured
                    SolidificationEventCounter(D3D1ConvPosition)++;
                    // Did the cell do this for the last time in the layer?
                    if (SolidificationEventCounter(D3D1ConvPosition) ==
                       NumberOfSolidificationEvents(D3D1ConvPosition)) {
                       // If so, this cell is done solidifying for now - change to wall type to ignore until next layer
                        // (if needed)
                        CellType(GlobalD3D1ConvPosition) = Solid;
                    }
                    else {
                       // Update MeltTimeStep, CritTimeStep, and UndercoolingChange with values for the next
                       // solidification event
                       CellType(GlobalD3D1ConvPosition) = TempSolid;
                       MeltTimeStep(GlobalD3D1ConvPosition) = (int)(LayerTimeTempHistory(
                            D3D1ConvPosition, SolidificationEventCounter(D3D1ConvPosition), 0));
                       CritTimeStep(GlobalD3D1ConvPosition) = (int)(LayerTimeTempHistory(
                           D3D1ConvPosition, SolidificationEventCounter(D3D1ConvPosition), 1));
                       UndercoolingChange(GlobalD3D1ConvPosition) =
                            LayerTimeTempHistory(D3D1ConvPosition, SolidificationEventCounter(D3D1ConvPosition), 2);
                    }
                }
            }   // End if statement cell capture
        }); // End parallel for loop over all undercooled active/aliquid cells
    Kokkos::fence();
    
}

//*****************************************************************************/
// Decentered octahedron algorithm for the capture of new interface cells by grains
void CellCapture_NoRM(int np, int cycle, int DecompositionStrategy, int LocalActiveDomainSize, int,
                 int MyXSlices, int MyYSlices, double AConst, double BConst, double CConst, double DConst,
                 int MyXOffset, int MyYOffset, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewI OppositeNeighbor,
                 ViewI CritTimeStep, ViewF UndercoolingCurrent, ViewF UndercoolingChange, ViewF GrainUnitVector,
                 ViewF CritDiagonalLength, ViewF DiagonalLength, ViewI GrainOrientation, ViewI CellType, ViewF DOCenter,
                 ViewI GrainID, int NGrainOrientations, Buffer2D BufferWestSend, Buffer2D BufferEastSend,
                 Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, Buffer2D BufferNorthEastSend,
                 Buffer2D BufferNorthWestSend, Buffer2D BufferSouthEastSend, Buffer2D BufferSouthWestSend, int BufSizeX,
                 int BufSizeY, int ZBound_Low, int nzActive, int, ViewI SteeringVector,
                 ViewI numSteer_G, ViewI_H numSteer_H, ViewI CaptureTimeStep) {

    // Cell capture - determing coordinates of all active cells, update undercooling in undercooled liquid cells/active
    // cells be accounted for
    Kokkos::parallel_for(
        "preCellCapture_NoRM", LocalActiveDomainSize, KOKKOS_LAMBDA(const long int &D3D1ConvPosition) {
            // Coordinate of this cell on the "global" (all cells in the Z direction) grid
            int RankZ = D3D1ConvPosition / (MyXSlices * MyYSlices);
            int Rem = D3D1ConvPosition % (MyXSlices * MyYSlices);
            int RankX = Rem / MyYSlices;
            int RankY = Rem % MyYSlices;
            int GlobalZ = RankZ + ZBound_Low;
            int GlobalD3D1ConvPosition = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;

            int cellType = CellType(GlobalD3D1ConvPosition);
            bool atCritTime = (cycle == CritTimeStep(GlobalD3D1ConvPosition));
            bool pastCritTime = (cycle >= CritTimeStep(GlobalD3D1ConvPosition));

            if ((cellType == Active) || (cellType == Liquid)) {
                if (atCritTime) {
                    // Cell just reached the liquidus temperature during cooling - set undercooling to 0
                    UndercoolingCurrent(GlobalD3D1ConvPosition) = 0.0;
                    SteeringVector(Kokkos::atomic_fetch_add(&numSteer_G(), 1)) = D3D1ConvPosition;
//                    if ((RankX == 7) && (RankY == 16) && (RankZ == 16))
//                        printf("Start solidification DL = %f cycle = %d \n",DiagonalLength(D3D1ConvPosition),cycle);
                }
                else if (pastCritTime) {
                    // Update local undercooling based on local cooling rate
                    UndercoolingCurrent(GlobalD3D1ConvPosition) += UndercoolingChange(GlobalD3D1ConvPosition);
                    if (cellType == Active) {
                        // Update local diagonal length of active cell
                        double LocU = UndercoolingCurrent(GlobalD3D1ConvPosition);
                        LocU = min(210.0, LocU);
                        double V = AConst * pow(LocU, 3.0) + BConst * pow(LocU, 2.0) + CConst * LocU + DConst;
                        V = max(0.0, V);
                        DiagonalLength(D3D1ConvPosition) += min(0.045, V); // Max amount the diagonal can grow per time step
//                        if ((RankX == 7) && (RankY == 16) && (RankZ == 16))
//                            printf("%f\n",DiagonalLength(D3D1ConvPosition));
                    }
                    SteeringVector(Kokkos::atomic_fetch_add(&numSteer_G(), 1)) = D3D1ConvPosition;
                }
            }
        });

    Kokkos::deep_copy(numSteer_H, numSteer_G);
 //   std::cout << "SV size: " << numSteer_H(0) << std::endl;
    Kokkos::parallel_for(
        "CellCapture", numSteer_H(), KOKKOS_LAMBDA(const int &num) {
            numSteer_G(0) = 0;
            int D3D1ConvPosition = SteeringVector(num);
            
            // Cells of interest for the CA
            int RankZ = D3D1ConvPosition / (MyXSlices * MyYSlices);
            int Rem = D3D1ConvPosition % (MyXSlices * MyYSlices);
            int RankX = Rem / MyYSlices;
            int RankY = Rem % MyYSlices;
            int GlobalZ = RankZ + ZBound_Low;
            int GlobalD3D1ConvPosition = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
            if (CellType(GlobalD3D1ConvPosition) == Active) {
                bool DeactivateCell = true;
               // if ((cycle == 10000) && (RankX <= 119)) printf("Cell at %d %d %d has diagonal length %f \n",RankX,RankY,RankZ,DiagonalLength(D3D1ConvPosition));
                for (int l = 0; l < 26; l++) {
                    // "l" correpsponds to the specific neighboring cell
                    // Local coordinates of adjacent cell center
                    int MyNeighborX = RankX + NeighborX(l);
                    int MyNeighborY = RankY + NeighborY(l);
                    int MyNeighborZ = RankZ + NeighborZ(l);
                    bool InBoundsThisRank = (MyNeighborX >= 0) * (MyNeighborX < MyXSlices) * (MyNeighborY >= 0) * (MyNeighborY < MyYSlices) * (MyNeighborZ >= 0) * (MyNeighborZ < nzActive);
                    if (InBoundsThisRank) {
//                        if ((cycle == 10000) && (RankX <= 119))
//                            printf("Checking neighbor %d %d %d \n",MyNeighborX,MyNeighborY,MyNeighborZ);
                        int GlobalNeighborD3D1ConvPosition = (MyNeighborZ + ZBound_Low) * MyXSlices * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
                        if (CellType(GlobalNeighborD3D1ConvPosition) == Liquid) {
                            // This cell will be remaining active - not all neighboring cells have been captured
                            DeactivateCell = false;
                            // Also, exit this loop over all the cell's neighbors and set check value to true
                            //l = 27;
//                            if ((cycle == 10000) && (RankX <= 119))  {
//                                printf("Diagonal length to capture cell at %d %d %d is %f \n",MyNeighborX,MyNeighborY,MyNeighborZ,CritDiagonalLength(26 * D3D1ConvPosition + l));
//                            }
                        }
                    }
                }
                if (DeactivateCell) {
                    CellType(GlobalD3D1ConvPosition) = Solid;
                }
            }
            else {
                // Loop over neighboring active cells to see if one of them, if any, captured this cell
                // Neighbor that will capture this cell (l = -1 if no capture event)
               // if ((cycle == 10000) && (RankX <= 119)) printf("Liquid cell at %d %d %d \n",RankX,RankY,RankZ);
                int CapturingNeighbor = -1;
                // Excess diagonal length (by how much did the potential capturing neighbor's diagonal length exceed the critical diagonal length)
                float CapturingNeighbor_ExcessDL = 0.0;
                int Alt_CapturingNeighbor = -1;
                int CaptureTimeExcess = 0;
                for (int l = 0; l < 26; l++) {
                    // "l" correpsponds to the specific neighboring cell
                    // Local coordinates of adjacent cell center
                    int MyNeighborX = RankX + NeighborX(l);
                    int MyNeighborY = RankY + NeighborY(l);
                    int MyNeighborZ = RankZ + NeighborZ(l);
                    bool InBoundsThisRank = (MyNeighborX >= 0) * (MyNeighborX < MyXSlices) * (MyNeighborY >= 0) * (MyNeighborY < MyYSlices) * (MyNeighborZ >= 0) * (MyNeighborZ < nzActive);
                    if (InBoundsThisRank) {
                        int NeighborD3D1ConvPosition = MyNeighborZ * MyXSlices * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
                        int GlobalNeighborD3D1ConvPosition = (MyNeighborZ + ZBound_Low) * MyXSlices * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
                        // If a neighboring cell is active, and that neighbor's diagonal length was large enough to capture this cell (the opposite of the vector from this cell to the neighbor):
                        // This neighbor may have captured this cell - only if the "overshoot" of the critical diagonal length relative to the diagonal length is larger than for any potential others
                        float ExcessDL = DiagonalLength(NeighborD3D1ConvPosition) - CritDiagonalLength(26 * NeighborD3D1ConvPosition + OppositeNeighbor(l));
                        int CaptureTimeExcess_ThisNeighbor = cycle - CaptureTimeStep(26 * NeighborD3D1ConvPosition + OppositeNeighbor(l));
                        if ((CellType(GlobalNeighborD3D1ConvPosition) == Active) && (ExcessDL >= 0.0) && (ExcessDL >= CapturingNeighbor_ExcessDL)) {
                            CapturingNeighbor = l;
                            CapturingNeighbor_ExcessDL = ExcessDL;
                          //  printf("cycle %d\n",cycle);//
                          //  printf("Potential cell capture of %d %d %d via DL criteria from cell %d %d %d (neighbor %d): ExcessDL = %f DL = %f CDL = %f (CaptureTimeExcess = %d )\n", RankX,RankY,RankZ,MyNeighborX,MyNeighborY,MyNeighborZ,l,ExcessDL,DiagonalLength(NeighborD3D1ConvPosition),CritDiagonalLength(26 * NeighborD3D1ConvPosition + OppositeNeighbor(l)),CaptureTimeExcess_ThisNeighbor);
                        }
                        if ((CellType(GlobalNeighborD3D1ConvPosition) == Active) && (CaptureTimeExcess_ThisNeighbor >= 0) && (CaptureTimeExcess_ThisNeighbor >= CaptureTimeExcess)) {
                            Alt_CapturingNeighbor = l;
                            CaptureTimeExcess = CaptureTimeExcess_ThisNeighbor;
                          //  printf("cycle %d\n",cycle);
                          //  printf("Potential cell capture of %d %d %d via time criteria from cell %d %d %d (neighbor %d): CaptureTimeExcess = %d (ExcessDL = %f DL = %f CDL = %f \n",RankX,RankY,RankZ,MyNeighborX,MyNeighborY,MyNeighborZ,l,CaptureTimeExcess,ExcessDL,DiagonalLength(NeighborD3D1ConvPosition), CritDiagonalLength(26 * NeighborD3D1ConvPosition + OppositeNeighbor(l)));
                        }
                    }
                }
//                if (CapturingNeighbor != Alt_CapturingNeighbor) {
//                    printf("Captured neighbor %d when should've captured neighbor %d \n",CapturingNeighbor,Alt_CapturingNeighbor);
//                }
                if (Alt_CapturingNeighbor != -1) {
                   // Cell capture event occurred - this liquid cell is now active

                   CellType(GlobalD3D1ConvPosition) = TemporaryUpdate; // don't turn to active until routine is complete to avoid race condition
                   // Perform cell capture of this cell from "CapturingCell"
                   // Local/Global coordinates of adjacent capturing cell center
                   int CapturingCellX = RankX + NeighborX(Alt_CapturingNeighbor);
                   int CapturingCellY = RankY + NeighborY(Alt_CapturingNeighbor);
                   int CapturingCellZ = RankZ + NeighborZ(Alt_CapturingNeighbor);
                   int LocalPositionCapturingCell = CapturingCellZ * MyXSlices * MyYSlices + CapturingCellX * MyYSlices + CapturingCellY;
                   int GlobalPositionCapturingCell = (CapturingCellZ + ZBound_Low) * MyXSlices * MyYSlices + CapturingCellX * MyYSlices + CapturingCellY;
                   int h = GrainID(GlobalPositionCapturingCell);
                   int MyOrientation = GrainOrientation(((abs(h) - 1) % NGrainOrientations));

                   // Cell is captured by the capturing cell's growing octahedron (Grain "h")
                   GrainID(GlobalD3D1ConvPosition) = h;

                   //printf("%d %d %d captured by GrainID %d (cycle %d) Capturing Neighbor %d AltCN was %d \n",RankX,RankY,RankZ,h,cycle,CapturingNeighbor,Alt_CapturingNeighbor);

                   // (cxold, cyold, czold) are the coordinates of the capturing cell's decentered octahedron
                   double cxold = DOCenter((long int)(3) * LocalPositionCapturingCell);
                   double cyold = DOCenter((long int)(3) * LocalPositionCapturingCell + (long int)(1));
                   double czold = DOCenter((long int)(3) * LocalPositionCapturingCell + (long int)(2));
                   // (xp,yp,zp) are the global coordinates of the cell's center
                   double xp = RankX + MyXOffset + 0.5;
                   double yp = RankY + MyYOffset + 0.5;
                   double zp = GlobalZ + 0.5;

                   // (x0,y0,z0) is a vector pointing from the capturing cell's decentered octahedron center to the image of
                   // the center of this cell
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

                    DiagonalLength(D3D1ConvPosition) = NewODiagL;
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

                    DOCenter((long int)(3) * D3D1ConvPosition) = cx;
                    DOCenter((long int)(3) * D3D1ConvPosition + (long int)(1)) = cy;
                    DOCenter((long int)(3) * D3D1ConvPosition + (long int)(2)) = cz;

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
                        double x0 = xp + NeighborX(n) - cx;
                        double y0 = yp + NeighborY(n) - cy;
                        double z0 = zp + NeighborZ(n) - cz;
                        D[0] = x0 * Fx[0] + y0 * Fy[0] + z0 * Fz[0];
                        D[1] = x0 * Fx[1] + y0 * Fy[1] + z0 * Fz[1];
                        D[2] = x0 * Fx[2] + y0 * Fy[2] + z0 * Fz[2];
                        D[3] = x0 * Fx[3] + y0 * Fy[3] + z0 * Fz[3];
                        Dfabs = max(max(fabs(D[0]), fabs(D[1])), max(fabs(D[2]), fabs(D[3])));
                        CritDiagonalLength((long int)(26) * D3D1ConvPosition + (long int)(n)) = Dfabs;
                       // printf("Ways to go: %f \n",Dfabs - NewODiagL);
//                        if (CritDiagonalLength((long int)(26) * D3D1ConvPosition + (long int)(n)) == 0.000000) {
//                            printf("Zero CDL %d %d \n",Alt_CapturingNeighbor,n);
//                        }
                        // Figure out which time step the diagonal length will reach the critical diagonal length
                        int future_cycle = cycle;
                        float CDL = Dfabs;
                        float DL = DiagonalLength(D3D1ConvPosition);
                        double Undercooling = UndercoolingCurrent(GlobalD3D1ConvPosition);
                        bool Inc_DL = true;
                        while (Inc_DL) {
                            if (DL >= CDL) {
                                CaptureTimeStep((long int)(26) * D3D1ConvPosition + (long int)(n)) = future_cycle;
//                                if (future_cycle <= cycle) {
//                                    CaptureTimeStep((long int)(26) * D3D1ConvPosition + (long int)(n)) = 100000000;
//                                }
                                Inc_DL = false;
                            }
                            else {
                                future_cycle++;
                                // Update local undercooling based on local cooling rate
                                Undercooling += UndercoolingChange(GlobalD3D1ConvPosition);
                                Undercooling = min(210.0, Undercooling);
                                double V = AConst * pow(Undercooling, 3.0) + BConst * pow(Undercooling, 2.0) + CConst * Undercooling + DConst;
                                V = max(0.0, V);
                                DL += min(0.045, V); // Max amount the diagonal can grow per time step
                            }
                        }
                    }
                    // Update ghost nodes due to cell activation
                    if (np > 1) {

                        float GhostGID = (float)(GrainID(GlobalD3D1ConvPosition));
                        float GhostDOCX = DOCenter((long int)(3) * D3D1ConvPosition);
                        float GhostDOCY = DOCenter((long int)(3) * D3D1ConvPosition + (long int)(1));
                        float GhostDOCZ = DOCenter((long int)(3) * D3D1ConvPosition + (long int)(2));
                        float GhostDL = DiagonalLength(D3D1ConvPosition);
                        
                        // Collect data for the ghost nodes:
                        if (DecompositionStrategy == 1) {
                            int GNPosition = RankZ * BufSizeX + RankX;
                            if (RankY == 1) {
                                BufferSouthSend(GNPosition, 0) = GhostGID;
                                BufferSouthSend(GNPosition, 1) = GhostDOCX;
                                BufferSouthSend(GNPosition, 2) = GhostDOCY;
                                BufferSouthSend(GNPosition, 3) = GhostDOCZ;
                                BufferSouthSend(GNPosition, 4) = GhostDL;
                            }
                            else if (RankY == MyYSlices - 2) {
                                int GNPosition = RankZ * BufSizeX + RankX;
                                BufferNorthSend(GNPosition, 0) = GhostGID;
                                BufferNorthSend(GNPosition, 1) = GhostDOCX;
                                BufferNorthSend(GNPosition, 2) = GhostDOCY;
                                BufferNorthSend(GNPosition, 3) = GhostDOCZ;
                                BufferNorthSend(GNPosition, 4) = GhostDL;
                            }
                        }
                        else {
                            if (RankY == 1) {
                                if (RankX == MyXSlices - 2) {
                                    int GNPosition = RankZ * BufSizeX + RankX - 1;
                                    BufferSouthSend(GNPosition, 0) = GhostGID;
                                    BufferSouthSend(GNPosition, 1) = GhostDOCX;
                                    BufferSouthSend(GNPosition, 2) = GhostDOCY;
                                    BufferSouthSend(GNPosition, 3) = GhostDOCZ;
                                    BufferSouthSend(GNPosition, 4) = GhostDL;
                                    GNPosition = RankZ * BufSizeY + RankY - 1;
                                    BufferEastSend(GNPosition, 0) = GhostGID;
                                    BufferEastSend(GNPosition, 1) = GhostDOCX;
                                    BufferEastSend(GNPosition, 2) = GhostDOCY;
                                    BufferEastSend(GNPosition, 3) = GhostDOCZ;
                                    BufferEastSend(GNPosition, 4) = GhostDL;
                                    GNPosition = RankZ;
                                    BufferSouthEastSend(GNPosition, 0) = GhostGID;
                                    BufferSouthEastSend(GNPosition, 1) = GhostDOCX;
                                    BufferSouthEastSend(GNPosition, 2) = GhostDOCY;
                                    BufferSouthEastSend(GNPosition, 3) = GhostDOCZ;
                                    BufferSouthEastSend(GNPosition, 4) = GhostDL;
                                }
                                else if (RankX == 1) {
                                    int GNPosition = RankZ * BufSizeX + RankX - 1;
                                    BufferSouthSend(GNPosition, 0) = GhostGID;
                                    BufferSouthSend(GNPosition, 1) = GhostDOCX;
                                    BufferSouthSend(GNPosition, 2) = GhostDOCY;
                                    BufferSouthSend(GNPosition, 3) = GhostDOCZ;
                                    BufferSouthSend(GNPosition, 4) = GhostDL;
                                    GNPosition = RankZ * BufSizeY + RankY - 1;
                                    BufferWestSend(GNPosition, 0) = GhostGID;
                                    BufferWestSend(GNPosition, 1) = GhostDOCX;
                                    BufferWestSend(GNPosition, 2) = GhostDOCY;
                                    BufferWestSend(GNPosition, 3) = GhostDOCZ;
                                    BufferWestSend(GNPosition, 4) = GhostDL;
                                    GNPosition = RankZ;
                                    BufferSouthWestSend(GNPosition, 0) = GhostGID;
                                    BufferSouthWestSend(GNPosition, 1) = GhostDOCX;
                                    BufferSouthWestSend(GNPosition, 2) = GhostDOCY;
                                    BufferSouthWestSend(GNPosition, 3) = GhostDOCZ;
                                    BufferSouthWestSend(GNPosition, 4) = GhostDL;
                                }
                                else if ((RankX > 1) && (RankX < MyXSlices - 2)) {
                                    // This is being sent to MyLeft
                                    int GNPosition = RankZ * BufSizeX + RankX - 1;
                                    BufferSouthSend(GNPosition, 0) = GhostGID;
                                    BufferSouthSend(GNPosition, 1) = GhostDOCX;
                                    BufferSouthSend(GNPosition, 2) = GhostDOCY;
                                    BufferSouthSend(GNPosition, 3) = GhostDOCZ;
                                    BufferSouthSend(GNPosition, 4) = GhostDL;
                                }
                            }
                            else if (RankY == MyYSlices - 2) {
                                // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                                if (RankX == MyXSlices - 2) {
                                    int GNPosition = RankZ * BufSizeX + RankX - 1;
                                    BufferNorthSend(GNPosition, 0) = GhostGID;
                                    BufferNorthSend(GNPosition, 1) = GhostDOCX;
                                    BufferNorthSend(GNPosition, 2) = GhostDOCY;
                                    BufferNorthSend(GNPosition, 3) = GhostDOCZ;
                                    BufferNorthSend(GNPosition, 4) = GhostDL;
                                    GNPosition = RankZ * BufSizeY + RankY - 1;
                                    BufferEastSend(GNPosition, 0) = GhostGID;
                                    BufferEastSend(GNPosition, 1) = GhostDOCX;
                                    BufferEastSend(GNPosition, 2) = GhostDOCY;
                                    BufferEastSend(GNPosition, 3) = GhostDOCZ;
                                    BufferEastSend(GNPosition, 4) = GhostDL;
                                    GNPosition = RankZ;
                                    BufferNorthEastSend(GNPosition, 0) = GhostGID;
                                    BufferNorthEastSend(GNPosition, 1) = GhostDOCX;
                                    BufferNorthEastSend(GNPosition, 2) = GhostDOCY;
                                    BufferNorthEastSend(GNPosition, 3) = GhostDOCZ;
                                    BufferNorthEastSend(GNPosition, 4) = GhostDL;
                                }
                                else if (RankX == 1) {
                                    int GNPosition = RankZ * BufSizeX + RankX - 1;
                                    BufferNorthSend(GNPosition, 0) = GhostGID;
                                    BufferNorthSend(GNPosition, 1) = GhostDOCX;
                                    BufferNorthSend(GNPosition, 2) = GhostDOCY;
                                    BufferNorthSend(GNPosition, 3) = GhostDOCZ;
                                    BufferNorthSend(GNPosition, 4) = GhostDL;
                                    GNPosition = RankZ * BufSizeY + RankY - 1;
                                    BufferWestSend(GNPosition, 0) = GhostGID;
                                    BufferWestSend(GNPosition, 1) = GhostDOCX;
                                    BufferWestSend(GNPosition, 2) = GhostDOCY;
                                    BufferWestSend(GNPosition, 3) = GhostDOCZ;
                                    BufferWestSend(GNPosition, 4) = GhostDL;
                                    GNPosition = RankZ;
                                    BufferNorthWestSend(GNPosition, 0) = GhostGID;
                                    BufferNorthWestSend(GNPosition, 1) = GhostDOCX;
                                    BufferNorthWestSend(GNPosition, 2) = GhostDOCY;
                                    BufferNorthWestSend(GNPosition, 3) = GhostDOCZ;
                                    BufferNorthWestSend(GNPosition, 4) = GhostDL;
                                }
                                else if ((RankX > 1) && (RankX < MyXSlices - 2)) {
                                    int GNPosition = RankZ * BufSizeX + RankX - 1;
                                    BufferNorthSend(GNPosition, 0) = GhostGID;
                                    BufferNorthSend(GNPosition, 1) = GhostDOCX;
                                    BufferNorthSend(GNPosition, 2) = GhostDOCY;
                                    BufferNorthSend(GNPosition, 3) = GhostDOCZ;
                                    BufferNorthSend(GNPosition, 4) = GhostDL;
                                }
                            }
                            else if ((RankX == 1) && (RankY > 1) && (RankY < MyYSlices - 2)) {
                                int GNPosition = RankZ * BufSizeY + RankY - 1;
                                BufferWestSend(GNPosition, 0) = GhostGID;
                                BufferWestSend(GNPosition, 1) = GhostDOCX;
                                BufferWestSend(GNPosition, 2) = GhostDOCY;
                                BufferWestSend(GNPosition, 3) = GhostDOCZ;
                                BufferWestSend(GNPosition, 4) = GhostDL;
                            }
                            else if ((RankX == MyXSlices - 2) && (RankY > 1) &&
                                     (RankY < MyYSlices - 2)) {
                                int GNPosition = RankZ * BufSizeY + RankY - 1;
                                BufferEastSend(GNPosition, 0) = GhostGID;
                                BufferEastSend(GNPosition, 1) = GhostDOCX;
                                BufferEastSend(GNPosition, 2) = GhostDOCY;
                                BufferEastSend(GNPosition, 3) = GhostDOCZ;
                                BufferEastSend(GNPosition, 4) = GhostDL;
                            }
                        } // End if statement for ghost node marking
                    }   // End if statement for serial/parallel code
                    // Only update the new cell's type once Critical Diagonal Length and Diagonal Length values have been assigned to it
                    // Avoids the race condition in which the new cell is activated, and another thread acts on the new active cell before the cell's new critical diagonal length/diagonal length values are assigned
                    CellType(GlobalD3D1ConvPosition) = Active;
                }   // End if statement cell capture
            } // End if statement for undercooled liquid cells
        }); // End parallel for loop over all undercooled active/aliquid cells
    Kokkos::fence();
}

//*****************************************************************************/
// Prints intermediate code output to stdout, checks to see if solidification is complete
void IntermediateOutputAndCheck(int id, int np, int &cycle, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int LocalDomainSize,
                                int LocalActiveDomainSize, int nx, int ny, int nz, int nzActive, double deltax,
                                float XMin, float YMin, float ZMin, int DecompositionStrategy,
                                int ProcessorsInXDirection, int ProcessorsInYDirection, int nn, int &XSwitch,
                                ViewI CellType, ViewI_H CellType_H, ViewI CritTimeStep, ViewI_H CritTimeStep_H,
                                ViewI GrainID, ViewI_H GrainID_H, std::string TemperatureDataType, int *FinishTimeStep,
                                int layernumber, int, int ZBound_Low, int NGrainOrientations, bool *Melted,
                                ViewI LayerID, ViewI_H LayerID_H, ViewI_H GrainOrientation_H, ViewF_H GrainUnitVector_H,
                                ViewF_H UndercoolingChange_H, std::string PathToOutput,
                                std::string OutputFile, bool PrintIdleMovieFrames, int MovieFrameInc,
                                int &IntermediateFileCounter) {

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

    if ((XSwitch == 0) && ((TemperatureDataType == "R") || (TemperatureDataType == "S"))) {
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
                if (PrintIdleMovieFrames) {
                    // Print any movie frames that occur during the skipped time steps
                    for (unsigned long int cycle_jump = cycle + 1; cycle_jump < GlobalNextCTS; cycle_jump++) {
                        if (cycle_jump % MovieFrameInc == 0) {
                            // Print current state of ExaCA simulation (up to and including the current layer's data)
                            Kokkos::deep_copy(GrainID_H, GrainID);
                            Kokkos::deep_copy(CellType_H, CellType);
                            PrintExaCAData(
                                id, layernumber, np, nx, ny, nz, MyXSlices, MyYSlices, MyXOffset, MyYOffset, ProcessorsInXDirection,
                                ProcessorsInYDirection, GrainID_H, GrainOrientation_H, CritTimeStep_H,
                                GrainUnitVector_H, LayerID_H, CellType_H, UndercoolingChange_H,
                                OutputFile, DecompositionStrategy, NGrainOrientations, Melted, PathToOutput, 0, false,
                                false, true, IntermediateFileCounter, ZBound_Low, nzActive, deltax, XMin, YMin, ZMin);
                            IntermediateFileCounter++;
                        }
                    }
                }
                // Jump to next time step when solidification starts again
                cycle = GlobalNextCTS - 1;
                if (id == 0)
                    std::cout << "Jumping to cycle " << cycle + 1 << std::endl;
            }
            if ((TemperatureDataType == "R") && (cycle >= FinishTimeStep[layernumber]))
                XSwitch = 1;
        }
    }
}

// Prints intermediate code output to stdout, checks if solidification is complete in the case of an extended data
// format simulation
void IntermediateOutputAndCheck_Remelt(
    int id, int np, int &cycle, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int LocalActiveDomainSize, int nx, int ny, int nz,
    int nzActive, double deltax, float XMin, float YMin, float ZMin, int DecompositionStrategy,
    int ProcessorsInXDirection, int ProcessorsInYDirection, int nn, int &XSwitch, ViewI CellType, ViewI_H CellType_H,
    ViewI_H CritTimeStep_H, ViewI MeltTimeStep, ViewI GrainID, ViewI_H GrainID_H, int layernumber, int, int ZBound_Low,
    int NGrainOrientations, bool *Melted, ViewI LayerID, ViewI_H LayerID_H, ViewI_H GrainOrientation_H,
    ViewF_H GrainUnitVector_H, ViewF_H UndercoolingChange_H, std::string PathToOutput,
    std::string OutputFile, bool PrintIdleMovieFrames, int MovieFrameInc, int &IntermediateFileCounter) {

    sample::ValueType CellTypeStorage;
    Kokkos::parallel_reduce(
        LocalActiveDomainSize,
        KOKKOS_LAMBDA(const int &D3D1ConvPosition, sample::ValueType &upd) {
            int RankZ = D3D1ConvPosition / (MyXSlices * MyYSlices);
            int Rem = D3D1ConvPosition % (MyXSlices * MyYSlices);
            int RankX = Rem / MyYSlices;
            int RankY = Rem % MyYSlices;
            int GlobalZ = RankZ + ZBound_Low;
            int GlobalD3D1ConvPosition = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
            if (CellType(GlobalD3D1ConvPosition) == Liquid)
                upd.the_array[0] += 1;
            else if (CellType(GlobalD3D1ConvPosition) == Active)
                upd.the_array[1] += 1;
            else if (CellType(GlobalD3D1ConvPosition) == TempSolid)
                upd.the_array[2] += 1;
            else if (CellType(GlobalD3D1ConvPosition) == Solid)
                upd.the_array[3] += 1;
        },
        Kokkos::Sum<sample::ValueType>(CellTypeStorage));
    unsigned long int Global_nn = 0;
    unsigned long int LocalLiquidCells = CellTypeStorage.the_array[0];
    unsigned long int LocalActiveCells = CellTypeStorage.the_array[1];
    unsigned long int LocalSolidCells = CellTypeStorage.the_array[2];
    unsigned long int LocalWallCells = CellTypeStorage.the_array[3];

    unsigned long int GlobalLiquidCells, GlobalActiveCells, GlobalSolidCells, GlobalWallCells;
    MPI_Reduce(&LocalLiquidCells, &GlobalLiquidCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalActiveCells, &GlobalActiveCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalSolidCells, &GlobalSolidCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&LocalWallCells, &GlobalWallCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&nn, &Global_nn, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    int SteeringVectorCells = GlobalActiveCells + GlobalLiquidCells;
    if (id == 0) {
        std::cout << "In layer number " << layernumber << " cycle = " << cycle
                  << " Liquid cells = " << GlobalLiquidCells
                  << " Active cells = " << GlobalActiveCells
                  << " Solid cells = " << GlobalSolidCells
                  << " Cells finished with solidification = " << GlobalWallCells
                  << " Number of nucleation events this layer " << Global_nn << std::endl;
        // Layer is done when all solid cells have become wall cells (they've finished melting and resolidifying)
        // and when all liquid cells have solidified
        if (GlobalSolidCells + GlobalLiquidCells == 0)
            XSwitch = 1;
    }
    MPI_Bcast(&XSwitch, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (XSwitch == 0) {
        MPI_Bcast(&SteeringVectorCells, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (SteeringVectorCells == 0) {
            // Check when the next cells undergo melting
            unsigned long int NextMTS;
            Kokkos::parallel_reduce(
                "NextMeltCheck", LocalActiveDomainSize,
                KOKKOS_LAMBDA(const int &D3D1ConvPosition, unsigned long int &tempv) {
                    int RankZ = D3D1ConvPosition / (MyXSlices * MyYSlices);
                    int Rem = D3D1ConvPosition % (MyXSlices * MyYSlices);
                    int RankX = Rem / MyYSlices;
                    int RankY = Rem % MyYSlices;
                    int GlobalZ = RankZ + ZBound_Low;
                    int GlobalD3D1ConvPosition = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                    unsigned long int MeltTimeStep_ThisCell = (unsigned long int)(MeltTimeStep(GlobalD3D1ConvPosition));
                    if ((LayerID(GlobalD3D1ConvPosition) == layernumber) &&
                        (CellType(GlobalD3D1ConvPosition) != Solid)) {
                        if (MeltTimeStep_ThisCell < tempv) {
                            tempv = MeltTimeStep_ThisCell;
                        }
                    }
                },
                Kokkos::Min<unsigned long int>(NextMTS));

            unsigned long int GlobalNextMTS;
            MPI_Allreduce(&NextMTS, &GlobalNextMTS, 1, MPI_UNSIGNED_LONG, MPI_MIN, MPI_COMM_WORLD);
            if ((GlobalNextMTS - cycle) > 5000) {
                if (PrintIdleMovieFrames) {
                    // Print any movie frames that occur during the skipped time steps
                    for (unsigned long int cycle_jump = cycle + 1; cycle_jump < GlobalNextMTS; cycle_jump++) {
                        if (cycle_jump % MovieFrameInc == 0) {
                            // Print current state of ExaCA simulation (up to and including the current layer's data)
                            Kokkos::deep_copy(GrainID_H, GrainID);
                            Kokkos::deep_copy(CellType_H, CellType);
                            PrintExaCAData(
                                id, layernumber, np, nx, ny, nz, MyXSlices, MyYSlices, MyXOffset, MyYOffset, ProcessorsInXDirection,
                                ProcessorsInYDirection, GrainID_H, GrainOrientation_H, CritTimeStep_H,
                                GrainUnitVector_H, LayerID_H, CellType_H, UndercoolingChange_H,
                                OutputFile, DecompositionStrategy, NGrainOrientations, Melted, PathToOutput, 0, false,
                                false, true, IntermediateFileCounter, ZBound_Low, nzActive, deltax, XMin, YMin, ZMin);
                            IntermediateFileCounter++;
                        }
                    }
                }
                // Jump to next time step when melting/solidification starts again
                cycle = GlobalNextMTS - 1;
                if (id == 0)
                    std::cout << "Jumping to cycle " << cycle + 1 << std::endl;
            }
        }
    }
}
