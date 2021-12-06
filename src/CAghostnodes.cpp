// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAghostnodes.hpp"

#include "mpi.h"

#include <cmath>
#include <vector>

//*****************************************************************************/
// Initial placement of data in ghost nodes
void GhostNodesInit(int, int, int DecompositionStrategy, int NeighborRank_North, int NeighborRank_South,
                    int NeighborRank_East, int NeighborRank_West, int NeighborRank_NorthEast,
                    int NeighborRank_NorthWest, int NeighborRank_SouthEast, int NeighborRank_SouthWest, int MyXSlices,
                    int MyYSlices, int MyXOffset, int MyYOffset, int ZBound_Low, int nzActive,
                    int LocalActiveDomainSize, int NGrainOrientations, ViewI NeighborX, ViewI NeighborY,
                    ViewI NeighborZ, ViewF GrainUnitVector, ViewI GrainOrientation, ViewI GrainID, ViewI CellType,
                    ViewF DOCenter, ViewF DiagonalLength, ViewF CritDiagonalLength) {

    // Fill buffers with ghost node data following initialization of data on GPUs
    // Similar to the calls to GhostNodes1D/GhostNodes2D, but the information sent/received in the halo regions is
    // different Need to send and receive cell type data, grain id data from other ranks
    Buffer3D BufferSouthSend(Kokkos::ViewAllocateWithoutInitializing("BufferSouthSend"), MyXSlices, nzActive, 2);
    Buffer3D BufferNorthSend(Kokkos::ViewAllocateWithoutInitializing("BufferNorthSend"), MyXSlices, nzActive, 2);
    Buffer3D BufferSouthRecv(Kokkos::ViewAllocateWithoutInitializing("BufferSouthRecv"), MyXSlices, nzActive, 2);
    Buffer3D BufferNorthRecv(Kokkos::ViewAllocateWithoutInitializing("BufferNorthRecv"), MyXSlices, nzActive, 2);

    Buffer3D BufferWestSend, BufferEastSend, BufferWestRecv, BufferEastRecv;
    Buffer2D BufferNorthEastSend, BufferNorthWestSend, BufferSouthWestSend, BufferSouthEastSend;
    Buffer2D BufferNorthEastRecv, BufferNorthWestRecv, BufferSouthWestRecv, BufferSouthEastRecv;
    if (DecompositionStrategy > 1) {
        BufferWestSend = Buffer3D(Kokkos::ViewAllocateWithoutInitializing("BufferWestSend"), MyYSlices, nzActive, 2);
        BufferEastSend = Buffer3D(Kokkos::ViewAllocateWithoutInitializing("BufferEastSend"), MyYSlices, nzActive, 2);
        BufferNorthEastSend = Buffer2D(Kokkos::ViewAllocateWithoutInitializing("BufferNorthEastSend"), nzActive, 2);
        BufferNorthWestSend = Buffer2D(Kokkos::ViewAllocateWithoutInitializing("BufferNorthWestSend"), nzActive, 2);
        BufferSouthWestSend = Buffer2D(Kokkos::ViewAllocateWithoutInitializing("BufferSouthWestSend"), nzActive, 2);
        BufferSouthEastSend = Buffer2D(Kokkos::ViewAllocateWithoutInitializing("BufferSouthEastSend"), nzActive, 2);
        BufferWestRecv = Buffer3D(Kokkos::ViewAllocateWithoutInitializing("BufferWestRecv"), MyYSlices, nzActive, 2);
        BufferEastRecv = Buffer3D(Kokkos::ViewAllocateWithoutInitializing("BufferEastRecv"), MyYSlices, nzActive, 2);
        BufferNorthEastRecv = Buffer2D(Kokkos::ViewAllocateWithoutInitializing("BufferNorthEastRecv"), nzActive, 2);
        BufferNorthWestRecv = Buffer2D(Kokkos::ViewAllocateWithoutInitializing("BufferNorthWestRecv"), nzActive, 2);
        BufferSouthWestRecv = Buffer2D(Kokkos::ViewAllocateWithoutInitializing("BufferSouthWestRecv"), nzActive, 2);
        BufferSouthEastRecv = Buffer2D(Kokkos::ViewAllocateWithoutInitializing("BufferSouthEastRecv"), nzActive, 2);
    }

    // Load send buffers
    Kokkos::parallel_for(
        "GNInit", LocalActiveDomainSize, KOKKOS_LAMBDA(const long int &D3D1ConvPosition) {
            int RankZ = D3D1ConvPosition / (MyXSlices * MyYSlices);
            int Rem = D3D1ConvPosition % (MyXSlices * MyYSlices);
            int RankX = Rem / MyYSlices;
            int RankY = Rem % MyYSlices;
            int GlobalZ = RankZ + ZBound_Low;
            int D3D1ConvPositionGlobal = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
            int GhostGID = GrainID(D3D1ConvPositionGlobal);
            int GhostCT = CellType(D3D1ConvPositionGlobal);
            // Collect data for the ghost nodes:
            if (RankY == 1) {
                BufferSouthSend(RankX, RankZ, 0) = GhostGID;
                BufferSouthSend(RankX, RankZ, 1) = GhostCT;
            }
            else if (RankY == MyYSlices - 2) {
                BufferNorthSend(RankX, RankZ, 0) = GhostGID;
                BufferNorthSend(RankX, RankZ, 1) = GhostCT;
            }
            if (DecompositionStrategy > 1) {
                if (RankX == MyXSlices - 2) {
                    BufferEastSend(RankY, RankZ, 0) = GhostGID;
                    BufferEastSend(RankY, RankZ, 1) = GhostCT;
                }
                else if (RankX == 1) {
                    BufferWestSend(RankY, RankZ, 0) = GhostGID;
                    BufferWestSend(RankY, RankZ, 1) = GhostCT;
                }
                if ((RankY == 1) && (RankX == 1)) {
                    BufferSouthWestSend(RankZ, 0) = GhostGID;
                    BufferSouthWestSend(RankZ, 1) = GhostCT;
                }
                else if ((RankY == 1) && (RankX == MyXSlices - 2)) {
                    BufferSouthEastSend(RankZ, 0) = GhostGID;
                    BufferSouthEastSend(RankZ, 1) = GhostCT;
                }
                else if ((RankY == MyYSlices - 2) && (RankX == MyXSlices - 2)) {
                    BufferNorthEastSend(RankZ, 0) = GhostGID;
                    BufferNorthEastSend(RankZ, 1) = GhostCT;
                }
                else if ((RankY == MyYSlices - 2) && (RankX == 1)) {
                    BufferNorthWestSend(RankZ, 0) = GhostGID;
                    BufferNorthWestSend(RankZ, 1) = GhostCT;
                }
            }
        });
    Kokkos::fence();

    int NumberOfSends;
    if (DecompositionStrategy > 1)
        NumberOfSends = 8;
    else
        NumberOfSends = 2;
    std::vector<MPI_Request> SendRequests(NumberOfSends, MPI_REQUEST_NULL);
    std::vector<MPI_Request> RecvRequests(NumberOfSends, MPI_REQUEST_NULL);

    // Send data to each other rank (MPI_Isend)
    MPI_Isend(BufferSouthSend.data(), 2 * MyXSlices * nzActive, MPI_INT, NeighborRank_South, 0, MPI_COMM_WORLD,
              &SendRequests[0]);
    MPI_Isend(BufferNorthSend.data(), 2 * MyXSlices * nzActive, MPI_INT, NeighborRank_North, 0, MPI_COMM_WORLD,
              &SendRequests[1]);
    if (DecompositionStrategy > 1) {
        MPI_Isend(BufferEastSend.data(), 2 * MyYSlices * nzActive, MPI_INT, NeighborRank_East, 0, MPI_COMM_WORLD,
                  &SendRequests[2]);
        MPI_Isend(BufferWestSend.data(), 2 * MyYSlices * nzActive, MPI_INT, NeighborRank_West, 0, MPI_COMM_WORLD,
                  &SendRequests[3]);
        MPI_Isend(BufferNorthWestSend.data(), 2 * nzActive, MPI_INT, NeighborRank_NorthWest, 0, MPI_COMM_WORLD,
                  &SendRequests[4]);
        MPI_Isend(BufferNorthEastSend.data(), 2 * nzActive, MPI_INT, NeighborRank_NorthEast, 0, MPI_COMM_WORLD,
                  &SendRequests[5]);
        MPI_Isend(BufferSouthWestSend.data(), 2 * nzActive, MPI_INT, NeighborRank_SouthWest, 0, MPI_COMM_WORLD,
                  &SendRequests[6]);
        MPI_Isend(BufferSouthEastSend.data(), 2 * nzActive, MPI_INT, NeighborRank_SouthEast, 0, MPI_COMM_WORLD,
                  &SendRequests[7]);
    }

    // Receive buffers for all neighbors (MPI_Irecv)
    MPI_Irecv(BufferSouthRecv.data(), 2 * MyXSlices * nzActive, MPI_INT, NeighborRank_South, 0, MPI_COMM_WORLD,
              &RecvRequests[0]);
    MPI_Irecv(BufferNorthRecv.data(), 2 * MyXSlices * nzActive, MPI_INT, NeighborRank_North, 0, MPI_COMM_WORLD,
              &RecvRequests[1]);
    if (DecompositionStrategy > 1) {
        MPI_Irecv(BufferEastRecv.data(), 2 * MyYSlices * nzActive, MPI_INT, NeighborRank_East, 0, MPI_COMM_WORLD,
                  &RecvRequests[2]);
        MPI_Irecv(BufferWestRecv.data(), 2 * MyYSlices * nzActive, MPI_INT, NeighborRank_West, 0, MPI_COMM_WORLD,
                  &RecvRequests[3]);
        MPI_Irecv(BufferNorthWestRecv.data(), 2 * nzActive, MPI_INT, NeighborRank_NorthWest, 0, MPI_COMM_WORLD,
                  &RecvRequests[4]);
        MPI_Irecv(BufferNorthEastRecv.data(), 2 * nzActive, MPI_INT, NeighborRank_NorthEast, 0, MPI_COMM_WORLD,
                  &RecvRequests[5]);
        MPI_Irecv(BufferSouthWestRecv.data(), 2 * nzActive, MPI_INT, NeighborRank_SouthWest, 0, MPI_COMM_WORLD,
                  &RecvRequests[6]);
        MPI_Irecv(BufferSouthEastRecv.data(), 2 * nzActive, MPI_INT, NeighborRank_SouthEast, 0, MPI_COMM_WORLD,
                  &RecvRequests[7]);
    }

    // unpack in any order
    bool unpack_complete = false;
    while (!unpack_complete) {
        // Get the next buffer to unpack from rank "unpack_index"
        int unpack_index = MPI_UNDEFINED;
        MPI_Waitany(NumberOfSends, RecvRequests.data(), &unpack_index, MPI_STATUS_IGNORE);

        // If there are no more buffers to unpack, leave the while loop
        if (MPI_UNDEFINED == unpack_index) {
            unpack_complete = true;
        }
        // Otherwise unpack the next buffer.
        else {
            Kokkos::parallel_for(
                "BufferUnpack", nzActive, KOKKOS_LAMBDA(const int &RankZ) {
                    int RankX, RankY;
                    int GlobalZ = RankZ + ZBound_Low;

                    // Which rank was the data received from?
                    if (((unpack_index == 0) && (NeighborRank_South != MPI_PROC_NULL)) ||
                        ((unpack_index == 1) && (NeighborRank_North != MPI_PROC_NULL))) {
                        // Data receieved from North or South
                        if (unpack_index == 0)
                            RankY = 0;
                        else
                            RankY = MyYSlices - 1;
                        // If 1D decomposition, unpack all X coordinates from buffer
                        int XIterationRangeStart = 0;
                        int XIterationRangeEnd = MyXSlices;
                        // If 2D decomposition, unpack all X coordinates except the two on the ends from the buffer
                        // (those are filled using corner ghost nodes)
                        if (DecompositionStrategy > 1) {
                            XIterationRangeStart = 1;
                            XIterationRangeEnd = MyXSlices - 1;
                        }
                        for (int RankX = XIterationRangeStart; RankX < XIterationRangeEnd; RankX++) {
                            int GlobalCellLocation = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                            if (unpack_index == 0) {
                                GrainID(GlobalCellLocation) = BufferSouthRecv(RankX, RankZ, 0);
                                CellType(GlobalCellLocation) = BufferSouthRecv(RankX, RankZ, 1);
                            }
                            else {
                                GrainID(GlobalCellLocation) = BufferNorthRecv(RankX, RankZ, 0);
                                CellType(GlobalCellLocation) = BufferNorthRecv(RankX, RankZ, 1);
                            }
                            if (CellType(GlobalCellLocation) == Active) {
                                // Mark cell for later placement of octahedra, critical diagonal length calculations
                                CellType(GlobalCellLocation) = Ghost;
                            }
                        }
                    }

                    else if (((unpack_index == 2) && (NeighborRank_East != MPI_PROC_NULL)) ||
                             ((unpack_index == 3) && (NeighborRank_West != MPI_PROC_NULL))) {
                        // Data received from East or West
                        if (unpack_index == 2)
                            RankX = MyXSlices - 1;
                        else
                            RankX = 0;
                        // This is a 2D decomposition, only unpack Y coordinates other than the two on the ends from the
                        // buffer (those are filled using corner ghost nodes)
                        for (int RankY = 1; RankY < MyYSlices - 1; RankY++) {
                            int GlobalCellLocation = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                            if (unpack_index == 2) {
                                GrainID(GlobalCellLocation) = BufferEastRecv(RankY, RankZ, 0);
                                CellType(GlobalCellLocation) = BufferEastRecv(RankY, RankZ, 1);
                            }
                            else {
                                GrainID(GlobalCellLocation) = BufferWestRecv(RankY, RankZ, 0);
                                CellType(GlobalCellLocation) = BufferWestRecv(RankY, RankZ, 1);
                            }
                            if (CellType(GlobalCellLocation) == Active) {
                                // Mark cell for later placement of octahedra, critical diagonal length calculations
                                CellType(GlobalCellLocation) = Ghost;
                            }
                        }
                    }
                    else {
                        if ((unpack_index == 4) && (NeighborRank_NorthWest != MPI_PROC_NULL)) {
                            // Data received from NorthWest
                            RankX = 0;
                            RankY = MyYSlices - 1;
                            int GlobalCellLocation = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                            GrainID(GlobalCellLocation) = BufferNorthWestRecv(RankZ, 0);
                            CellType(GlobalCellLocation) = BufferNorthWestRecv(RankZ, 1);
                            // Mark cell for later placement of octahedra, critical diagonal length calculations
                            if (CellType(GlobalCellLocation) == Active)
                                CellType(GlobalCellLocation) = Ghost;
                        }
                        else if ((unpack_index == 5) && (NeighborRank_NorthEast != MPI_PROC_NULL)) {
                            // Data received from NorthEast
                            RankX = MyXSlices - 1;
                            RankY = MyYSlices - 1;
                            int GlobalCellLocation = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                            GrainID(GlobalCellLocation) = BufferNorthEastRecv(RankZ, 0);
                            CellType(GlobalCellLocation) = BufferNorthEastRecv(RankZ, 1);
                            // Mark cell for later placement of octahedra, critical diagonal length calculations
                            if (CellType(GlobalCellLocation) == Active)
                                CellType(GlobalCellLocation) = Ghost;
                        }
                        else if ((unpack_index == 6) && (NeighborRank_SouthWest != MPI_PROC_NULL)) {
                            // Data received from SouthWest
                            RankX = 0;
                            RankY = 0;
                            int GlobalCellLocation = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                            GrainID(GlobalCellLocation) = BufferSouthWestRecv(RankZ, 0);
                            CellType(GlobalCellLocation) = BufferSouthWestRecv(RankZ, 1);
                            // Mark cell for later placement of octahedra, critical diagonal length calculations
                            if (CellType(GlobalCellLocation) == Active)
                                CellType(GlobalCellLocation) = Ghost;
                        }
                        else if ((unpack_index == 7) && (NeighborRank_SouthEast != MPI_PROC_NULL)) {
                            // Data received from SouthEast
                            RankX = MyXSlices - 1;
                            RankY = 0;
                            int GlobalCellLocation = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                            GrainID(GlobalCellLocation) = BufferSouthEastRecv(RankZ, 0);
                            CellType(GlobalCellLocation) = BufferSouthEastRecv(RankZ, 1);
                            // Mark cell for later placement of octahedra, critical diagonal length calculations
                            if (CellType(GlobalCellLocation) == Active)
                                CellType(GlobalCellLocation) = Ghost;
                        }
                    }
                });
            Kokkos::fence();
        }
    }
    // Wait on send requests
    MPI_Waitall(NumberOfSends, SendRequests.data(), MPI_STATUSES_IGNORE);

    // Octahedra for active cells that were marked in the previous loop, initially have centers aligned with cell
    // centers, Diagonal lengths of 0.01
    Kokkos::parallel_for(
        "GNInit", LocalActiveDomainSize, KOKKOS_LAMBDA(const int &CellLocation) {
            int RankZ = CellLocation / (MyXSlices * MyYSlices);
            int Rem = CellLocation % (MyXSlices * MyYSlices);
            int RankX = Rem / MyYSlices;
            int RankY = Rem % MyYSlices;
            int GlobalZ = RankZ + ZBound_Low;
            int GlobalCellLocation = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
            if (CellType(GlobalCellLocation) == Ghost) {
                CellType(GlobalCellLocation) = Active;
                DOCenter((long int)(3) * CellLocation) = RankX + MyXOffset + 0.5;
                DOCenter((long int)(3) * CellLocation + (long int)(1)) = RankY + MyYOffset + 0.5;
                DOCenter((long int)(3) * CellLocation + (long int)(2)) = GlobalZ + 0.5;
                int MyOrientation = GrainOrientation(((abs(GrainID(GlobalCellLocation)) - 1) % NGrainOrientations));
                DiagonalLength(CellLocation) = 0.01;
                // Global coordinates of cell center
                double xp = RankX + MyXOffset + 0.5;
                double yp = RankY + MyYOffset + 0.5;
                double zp = GlobalZ + 0.5;
                // Calculate critical values at which this active cell leads to the activation of a neighboring liquid
                // cell
                for (int n = 0; n < 26; n++) {

                    int MyNeighborX = RankX + NeighborX(n);
                    int MyNeighborY = RankY + NeighborY(n);
                    int MyNeighborZ = RankZ + NeighborZ(n);
                    long int NeighborPosition =
                        MyNeighborZ * MyXSlices * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;

                    if (NeighborPosition == CellLocation) {
                        // Do not calculate critical diagonal length req'd for the newly captured cell to capture the
                        // original
                        CritDiagonalLength((long int)(26) * NeighborPosition + (long int)(n)) = 10000000.0;
                    }

                    // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the global coordinates
                    // of the center of a neighbor cell
                    double x0 = xp + NeighborX(n) - (RankX + MyXOffset + 0.5);
                    double y0 = yp + NeighborY(n) - (RankY + MyYOffset + 0.5);
                    double z0 = zp + NeighborZ(n) - (GlobalZ + 0.5);
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
                    CritDiagonalLength((long int)(26) * CellLocation + (long int)(n)) = CDLVal;
                }
            }
        });
    Kokkos::fence();
}

//*****************************************************************************/
// 2D domain decomposition: update ghost nodes with new cell data from Nucleation and CellCapture routines
void GhostNodes2D(int, int, int NeighborRank_North, int NeighborRank_South, int NeighborRank_East,
                  int NeighborRank_West, int NeighborRank_NorthEast, int NeighborRank_NorthWest,
                  int NeighborRank_SouthEast, int NeighborRank_SouthWest, int MyXSlices, int MyYSlices, int MyXOffset,
                  int MyYOffset, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewI CellType, ViewF DOCenter,
                  ViewI GrainID, ViewF GrainUnitVector, ViewI GrainOrientation, ViewF DiagonalLength,
                  ViewF CritDiagonalLength, int NGrainOrientations, Buffer2D BufferWestSend, Buffer2D BufferEastSend,
                  Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, Buffer2D BufferNorthEastSend,
                  Buffer2D BufferNorthWestSend, Buffer2D BufferSouthEastSend, Buffer2D BufferSouthWestSend,
                  Buffer2D BufferWestRecv, Buffer2D BufferEastRecv, Buffer2D BufferNorthRecv, Buffer2D BufferSouthRecv,
                  Buffer2D BufferNorthEastRecv, Buffer2D BufferNorthWestRecv, Buffer2D BufferSouthEastRecv,
                  Buffer2D BufferSouthWestRecv, int BufSizeX, int BufSizeY, int nzActive, int ZBound_Low) {

    Kokkos::fence();
    MPI_Barrier(MPI_COMM_WORLD);

    // Allocate requests.
    std::vector<MPI_Request> SendRequests(8, MPI_REQUEST_NULL);
    std::vector<MPI_Request> RecvRequests(8, MPI_REQUEST_NULL);

    // Send data to each other rank (MPI_Isend)
    MPI_Isend(BufferSouthSend.data(), 5 * BufSizeX * nzActive, MPI_FLOAT, NeighborRank_South, 0, MPI_COMM_WORLD,
              &SendRequests[0]);
    MPI_Isend(BufferNorthSend.data(), 5 * BufSizeX * nzActive, MPI_FLOAT, NeighborRank_North, 0, MPI_COMM_WORLD,
              &SendRequests[1]);
    MPI_Isend(BufferEastSend.data(), 5 * BufSizeY * nzActive, MPI_FLOAT, NeighborRank_East, 0, MPI_COMM_WORLD,
              &SendRequests[2]);
    MPI_Isend(BufferWestSend.data(), 5 * BufSizeY * nzActive, MPI_FLOAT, NeighborRank_West, 0, MPI_COMM_WORLD,
              &SendRequests[3]);
    MPI_Isend(BufferNorthWestSend.data(), 5 * nzActive, MPI_FLOAT, NeighborRank_NorthWest, 0, MPI_COMM_WORLD,
              &SendRequests[4]);
    MPI_Isend(BufferNorthEastSend.data(), 5 * nzActive, MPI_FLOAT, NeighborRank_NorthEast, 0, MPI_COMM_WORLD,
              &SendRequests[5]);
    MPI_Isend(BufferSouthWestSend.data(), 5 * nzActive, MPI_FLOAT, NeighborRank_SouthWest, 0, MPI_COMM_WORLD,
              &SendRequests[6]);
    MPI_Isend(BufferSouthEastSend.data(), 5 * nzActive, MPI_FLOAT, NeighborRank_SouthEast, 0, MPI_COMM_WORLD,
              &SendRequests[7]);

    // Receive buffers for all neighbors (MPI_Irecv)
    MPI_Irecv(BufferSouthRecv.data(), 5 * BufSizeX * nzActive, MPI_FLOAT, NeighborRank_South, 0, MPI_COMM_WORLD,
              &RecvRequests[0]);
    MPI_Irecv(BufferNorthRecv.data(), 5 * BufSizeX * nzActive, MPI_FLOAT, NeighborRank_North, 0, MPI_COMM_WORLD,
              &RecvRequests[1]);
    MPI_Irecv(BufferEastRecv.data(), 5 * BufSizeY * nzActive, MPI_FLOAT, NeighborRank_East, 0, MPI_COMM_WORLD,
              &RecvRequests[2]);
    MPI_Irecv(BufferWestRecv.data(), 5 * BufSizeY * nzActive, MPI_FLOAT, NeighborRank_West, 0, MPI_COMM_WORLD,
              &RecvRequests[3]);
    MPI_Irecv(BufferNorthWestRecv.data(), 5 * nzActive, MPI_FLOAT, NeighborRank_NorthWest, 0, MPI_COMM_WORLD,
              &RecvRequests[4]);
    MPI_Irecv(BufferNorthEastRecv.data(), 5 * nzActive, MPI_FLOAT, NeighborRank_NorthEast, 0, MPI_COMM_WORLD,
              &RecvRequests[5]);
    MPI_Irecv(BufferSouthWestRecv.data(), 5 * nzActive, MPI_FLOAT, NeighborRank_SouthWest, 0, MPI_COMM_WORLD,
              &RecvRequests[6]);
    MPI_Irecv(BufferSouthEastRecv.data(), 5 * nzActive, MPI_FLOAT, NeighborRank_SouthEast, 0, MPI_COMM_WORLD,
              &RecvRequests[7]);

    // unpack in any order
    bool unpack_complete = false;
    while (!unpack_complete) {
        // Get the next buffer to unpack from rank "unpack_index"
        int unpack_index = MPI_UNDEFINED;
        MPI_Waitany(8, RecvRequests.data(), &unpack_index, MPI_STATUS_IGNORE);

        // If there are no more buffers to unpack, leave the while loop
        if (MPI_UNDEFINED == unpack_index) {
            unpack_complete = true;
        }
        // Otherwise unpack the next buffer.
        else {
            int RecvBufSize;
            if (unpack_index <= 1) {
                RecvBufSize = BufSizeX * nzActive;
            }
            else if (unpack_index >= 4) {
                RecvBufSize = nzActive;
            }
            else {
                RecvBufSize = BufSizeY * nzActive;
            }
            Kokkos::parallel_for(
                "BufferUnpack", RecvBufSize, KOKKOS_LAMBDA(const int &BufPosition) {
                    int RankX, RankY, RankZ, NewGrainID, GlobalZ, GlobalCellLocation;
                    long int CellLocation;
                    bool Place = false;
                    float DOCenterX, DOCenterY, DOCenterZ, NewDiagonalLength;
                    // Which rank was the data received from?
                    if ((unpack_index == 0) && (NeighborRank_South != MPI_PROC_NULL)) {
                        // Data receieved from South
                        // Adjust X Position by +1 since X = 0 is not included in this buffer
                        RankX = BufPosition % BufSizeX + 1;
                        RankY = 0;
                        RankZ = BufPosition / BufSizeX;
                        CellLocation = RankZ * MyXSlices * MyYSlices + MyYSlices * RankX + RankY;
                        GlobalZ = RankZ + ZBound_Low;
                        GlobalCellLocation = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                        if ((BufferSouthRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferSouthRecv(BufPosition, 0));
                            DOCenterX = BufferSouthRecv(BufPosition, 1);
                            DOCenterY = BufferSouthRecv(BufPosition, 2);
                            DOCenterZ = BufferSouthRecv(BufPosition, 3);
                            NewDiagonalLength = BufferSouthRecv(BufPosition, 4);
                        }
                    }
                    else if ((unpack_index == 1) && (NeighborRank_North != MPI_PROC_NULL)) {
                        // Data received from North
                        // Adjust X Position by +1 since X = 0 is not included in this buffer
                        RankX = BufPosition % BufSizeX + 1;
                        RankY = MyYSlices - 1;
                        RankZ = BufPosition / BufSizeX;
                        CellLocation = RankZ * MyXSlices * MyYSlices + MyYSlices * RankX + RankY;
                        GlobalZ = RankZ + ZBound_Low;
                        GlobalCellLocation = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                        if ((BufferNorthRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferNorthRecv(BufPosition, 0));
                            DOCenterX = BufferNorthRecv(BufPosition, 1);
                            DOCenterY = BufferNorthRecv(BufPosition, 2);
                            DOCenterZ = BufferNorthRecv(BufPosition, 3);
                            NewDiagonalLength = BufferNorthRecv(BufPosition, 4);
                        }
                    }
                    else if ((unpack_index == 2) && (NeighborRank_East != MPI_PROC_NULL)) {
                        // Data received from East
                        // Adjust Y Position by +1 since Y = 0 is not included in this buffer
                        RankX = MyXSlices - 1;
                        RankY = BufPosition % BufSizeY + 1;
                        RankZ = BufPosition / BufSizeY;
                        CellLocation = RankZ * MyXSlices * MyYSlices + MyYSlices * RankX + RankY;
                        GlobalZ = RankZ + ZBound_Low;
                        GlobalCellLocation = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                        if ((BufferEastRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferEastRecv(BufPosition, 0));
                            DOCenterX = BufferEastRecv(BufPosition, 1);
                            DOCenterY = BufferEastRecv(BufPosition, 2);
                            DOCenterZ = BufferEastRecv(BufPosition, 3);
                            NewDiagonalLength = BufferEastRecv(BufPosition, 4);
                        }
                    }
                    else if ((unpack_index == 3) && (NeighborRank_West != MPI_PROC_NULL)) {
                        // Data received from West
                        // Adjust Y Position by +1 since Y = 0 is not included in this buffer
                        RankX = 0;
                        RankY = BufPosition % BufSizeY + 1;
                        RankZ = BufPosition / BufSizeY;
                        CellLocation = RankZ * MyXSlices * MyYSlices + MyYSlices * RankX + RankY;
                        GlobalZ = RankZ + ZBound_Low;
                        GlobalCellLocation = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                        if ((BufferWestRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferWestRecv(BufPosition, 0));
                            DOCenterX = BufferWestRecv(BufPosition, 1);
                            DOCenterY = BufferWestRecv(BufPosition, 2);
                            DOCenterZ = BufferWestRecv(BufPosition, 3);
                            NewDiagonalLength = BufferWestRecv(BufPosition, 4);
                        }
                    }
                    else if ((unpack_index == 4) && (NeighborRank_NorthWest != MPI_PROC_NULL)) {
                        // Data received from NorthWest
                        RankX = 0;
                        RankY = MyYSlices - 1;
                        RankZ = BufPosition;
                        CellLocation = RankZ * MyXSlices * MyYSlices + MyYSlices * RankX + RankY;
                        GlobalZ = RankZ + ZBound_Low;
                        GlobalCellLocation = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                        if ((BufferNorthWestRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferNorthWestRecv(BufPosition, 0));
                            DOCenterX = BufferNorthWestRecv(BufPosition, 1);
                            DOCenterY = BufferNorthWestRecv(BufPosition, 2);
                            DOCenterZ = BufferNorthWestRecv(BufPosition, 3);
                            NewDiagonalLength = BufferNorthWestRecv(BufPosition, 4);
                        }
                    }
                    else if ((unpack_index == 5) && (NeighborRank_NorthEast != MPI_PROC_NULL)) {
                        // Data received from NorthEast
                        RankX = MyXSlices - 1;
                        RankY = MyYSlices - 1;
                        RankZ = BufPosition;
                        CellLocation = RankZ * MyXSlices * MyYSlices + MyYSlices * RankX + RankY;
                        GlobalZ = RankZ + ZBound_Low;
                        GlobalCellLocation = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                        if ((BufferNorthEastRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferNorthEastRecv(BufPosition, 0));
                            DOCenterX = BufferNorthEastRecv(BufPosition, 1);
                            DOCenterY = BufferNorthEastRecv(BufPosition, 2);
                            DOCenterZ = BufferNorthEastRecv(BufPosition, 3);
                            NewDiagonalLength = BufferNorthEastRecv(BufPosition, 4);
                        }
                    }
                    else if ((unpack_index == 6) && (NeighborRank_SouthWest != MPI_PROC_NULL)) {
                        // Data received from SouthWest
                        RankX = 0;
                        RankY = 0;
                        RankZ = BufPosition;
                        CellLocation = RankZ * MyXSlices * MyYSlices + MyYSlices * RankX + RankY;
                        GlobalZ = RankZ + ZBound_Low;
                        GlobalCellLocation = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                        if ((BufferSouthWestRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferSouthWestRecv(BufPosition, 0));
                            DOCenterX = BufferSouthWestRecv(BufPosition, 1);
                            DOCenterY = BufferSouthWestRecv(BufPosition, 2);
                            DOCenterZ = BufferSouthWestRecv(BufPosition, 3);
                            NewDiagonalLength = BufferSouthWestRecv(BufPosition, 4);
                        }
                    }
                    else if ((unpack_index == 7) && (NeighborRank_SouthEast != MPI_PROC_NULL)) {
                        // Data received from SouthEast
                        RankX = MyXSlices - 1;
                        RankY = 0;
                        RankZ = BufPosition;
                        CellLocation = RankZ * MyXSlices * MyYSlices + MyYSlices * RankX + RankY;
                        GlobalZ = RankZ + ZBound_Low;
                        GlobalCellLocation = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                        if ((BufferSouthEastRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferSouthEastRecv(BufPosition, 0));
                            DOCenterX = BufferSouthEastRecv(BufPosition, 1);
                            DOCenterY = BufferSouthEastRecv(BufPosition, 2);
                            DOCenterZ = BufferSouthEastRecv(BufPosition, 3);
                            NewDiagonalLength = BufferSouthEastRecv(BufPosition, 4);
                        }
                    }

                    if (Place) {
                        // Update this ghost node cell's information with data from other rank
                        CellType(GlobalCellLocation) = Active;
                        GrainID(GlobalCellLocation) = NewGrainID;
                        DOCenter((long int)(3) * CellLocation) = DOCenterX;
                        DOCenter((long int)(3) * CellLocation + (long int)(1)) = DOCenterY;
                        DOCenter((long int)(3) * CellLocation + (long int)(2)) = DOCenterZ;
                        int MyOrientation =
                            GrainOrientation(((abs(GrainID(GlobalCellLocation)) - 1) % NGrainOrientations));
                        DiagonalLength(CellLocation) = NewDiagonalLength;
                        // Global coordinates of cell center
                        double xp = RankX + MyXOffset + 0.5;
                        double yp = RankY + MyYOffset + 0.5;
                        double zp = GlobalZ + 0.5;
                        // Calculate critical values at which this active cell leads to the activation of a neighboring
                        // liquid cell
                        for (int n = 0; n < 26; n++) {

                            int MyNeighborX = RankX + NeighborX(n);
                            int MyNeighborY = RankY + NeighborY(n);
                            int MyNeighborZ = RankZ + NeighborZ(n);
                            long int NeighborPosition =
                                MyNeighborZ * MyXSlices * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;

                            if (NeighborPosition == CellLocation) {
                                // Do not calculate critical diagonal length req'd for the newly captured cell to
                                // capture the original
                                CritDiagonalLength((long int)(26) * NeighborPosition + (long int)(n)) = 10000000.0;
                            }

                            // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the global
                            // coordinates of the center of a neighbor cell
                            double x0 = xp + NeighborX(n) - DOCenterX;
                            double y0 = yp + NeighborY(n) - DOCenterY;
                            double z0 = zp + NeighborZ(n) - DOCenterZ;
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

                            double Angle2 = (GrainUnitVector(9 * MyOrientation + 3) * x0 +
                                             GrainUnitVector(9 * MyOrientation + 4) * y0 +
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

                            double Angle3 = (GrainUnitVector(9 * MyOrientation + 6) * x0 +
                                             GrainUnitVector(9 * MyOrientation + 7) * y0 +
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
                            double ParaT = (normx * x0 + normy * y0 + normz * z0) /
                                           (normx * Diag1X + normy * Diag1Y + normz * Diag1Z);
                            float CDLVal = pow(
                                pow(ParaT * Diag1X, 2.0) + pow(ParaT * Diag1Y, 2.0) + pow(ParaT * Diag1Z, 2.0), 0.5);
                            CritDiagonalLength((long int)(26) * CellLocation + (long int)(n)) = CDLVal;
                        }
                    }
                });
        }
    }
    // Wait on send requests
    MPI_Waitall(8, SendRequests.data(), MPI_STATUSES_IGNORE);
}

//*****************************************************************************/
// 1D domain decomposition: update ghost nodes with new cell data from Nucleation and CellCapture routines
void GhostNodes1D(int, int, int NeighborRank_North, int NeighborRank_South, int MyXSlices, int MyYSlices, int MyXOffset,
                  int MyYOffset, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewI CellType, ViewF DOCenter,
                  ViewI GrainID, ViewF GrainUnitVector, ViewI GrainOrientation, ViewF DiagonalLength,
                  ViewF CritDiagonalLength, int NGrainOrientations, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend,
                  Buffer2D BufferNorthRecv, Buffer2D BufferSouthRecv, int BufSizeX, int, int nzActive, int ZBound_Low) {

    std::vector<MPI_Request> SendRequests(2, MPI_REQUEST_NULL);
    std::vector<MPI_Request> RecvRequests(2, MPI_REQUEST_NULL);

    // Send data to each other rank (MPI_Isend)
    MPI_Isend(BufferSouthSend.data(), 5 * BufSizeX * nzActive, MPI_FLOAT, NeighborRank_South, 0, MPI_COMM_WORLD,
              &SendRequests[0]);
    MPI_Isend(BufferNorthSend.data(), 5 * BufSizeX * nzActive, MPI_FLOAT, NeighborRank_North, 0, MPI_COMM_WORLD,
              &SendRequests[1]);

    // Receive buffers for all neighbors (MPI_Irecv)
    MPI_Irecv(BufferSouthRecv.data(), 5 * BufSizeX * nzActive, MPI_FLOAT, NeighborRank_South, 0, MPI_COMM_WORLD,
              &RecvRequests[0]);
    MPI_Irecv(BufferNorthRecv.data(), 5 * BufSizeX * nzActive, MPI_FLOAT, NeighborRank_North, 0, MPI_COMM_WORLD,
              &RecvRequests[1]);

    // unpack in any order
    bool unpack_complete = false;
    while (!unpack_complete) {
        // Get the next buffer to unpack from rank "unpack_index"
        int unpack_index = MPI_UNDEFINED;
        MPI_Waitany(2, RecvRequests.data(), &unpack_index, MPI_STATUS_IGNORE);
        // If there are no more buffers to unpack, leave the while loop
        if (MPI_UNDEFINED == unpack_index) {
            unpack_complete = true;
        }
        // Otherwise unpack the next buffer.
        else {
            int RecvBufSize = BufSizeX * nzActive;
            Kokkos::parallel_for(
                "BufferUnpack", RecvBufSize, KOKKOS_LAMBDA(const int &BufPosition) {
                    int RankX, RankY, RankZ, NewGrainID, GlobalZ, GlobalCellLocation;
                    long int CellLocation;
                    float DOCenterX, DOCenterY, DOCenterZ, NewDiagonalLength;
                    bool Place = false;
                    RankZ = BufPosition / BufSizeX;
                    RankX = BufPosition % BufSizeX;
                    // Which rank was the data received from?
                    if ((unpack_index == 0) && (NeighborRank_South != MPI_PROC_NULL)) {
                        // Data receieved from South
                        RankY = 0;
                        CellLocation = RankZ * MyXSlices * MyYSlices + MyYSlices * RankX + RankY;
                        GlobalZ = RankZ + ZBound_Low;
                        GlobalCellLocation = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                        if ((BufferSouthRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferSouthRecv(BufPosition, 0));
                            DOCenterX = BufferSouthRecv(BufPosition, 1);
                            DOCenterY = BufferSouthRecv(BufPosition, 2);
                            DOCenterZ = BufferSouthRecv(BufPosition, 3);
                            NewDiagonalLength = BufferSouthRecv(BufPosition, 4);
                        }
                    }
                    else if ((unpack_index == 1) && (NeighborRank_North != MPI_PROC_NULL)) {
                        // Data received from North
                        RankY = MyYSlices - 1;
                        CellLocation = RankZ * MyXSlices * MyYSlices + MyYSlices * RankX + RankY;
                        GlobalZ = RankZ + ZBound_Low;
                        GlobalCellLocation = GlobalZ * MyXSlices * MyYSlices + RankX * MyYSlices + RankY;
                        if ((BufferNorthRecv(BufPosition, 4) > 0) && (CellType(GlobalCellLocation) == Liquid)) {
                            Place = true;
                            NewGrainID = (int)(BufferNorthRecv(BufPosition, 0));
                            DOCenterX = BufferNorthRecv(BufPosition, 1);
                            DOCenterY = BufferNorthRecv(BufPosition, 2);
                            DOCenterZ = BufferNorthRecv(BufPosition, 3);
                            NewDiagonalLength = BufferNorthRecv(BufPosition, 4);
                        }
                    }
                    if (Place) {
                        // Update this ghost node cell's information with data from other rank
                        GrainID(GlobalCellLocation) = NewGrainID;
                        DOCenter((long int)(3) * CellLocation) = DOCenterX;
                        DOCenter((long int)(3) * CellLocation + (long int)(1)) = DOCenterY;
                        DOCenter((long int)(3) * CellLocation + (long int)(2)) = DOCenterZ;
                        int MyOrientation =
                            GrainOrientation(((abs(GrainID(GlobalCellLocation)) - 1) % NGrainOrientations));
                        DiagonalLength(CellLocation) = NewDiagonalLength;
                        // Global coordinates of cell center
                        double xp = RankX + MyXOffset + 0.5;
                        double yp = RankY + MyYOffset + 0.5;
                        double zp = GlobalZ + 0.5;
                        // Calculate critical values at which this active cell leads to the activation of a neighboring
                        // liquid cell
                        for (int n = 0; n < 26; n++) {

                            int MyNeighborX = RankX + NeighborX(n);
                            int MyNeighborY = RankY + NeighborY(n);
                            int MyNeighborZ = RankZ + NeighborZ(n);
                            long int NeighborPosition =
                                MyNeighborZ * MyXSlices * MyYSlices + MyNeighborX * MyYSlices + MyNeighborY;
                            if (NeighborPosition == CellLocation) {
                                // Do not calculate critical diagonal length req'd for the newly captured cell to
                                // capture the original
                                CritDiagonalLength((long int)(26) * NeighborPosition + (long int)(n)) = 10000000.0;
                            }

                            // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the global
                            // coordinates of the center of a neighbor cell
                            double x0 = xp + NeighborX(n) - DOCenterX;
                            double y0 = yp + NeighborY(n) - DOCenterY;
                            double z0 = zp + NeighborZ(n) - DOCenterZ;
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

                            double Angle2 = (GrainUnitVector(9 * MyOrientation + 3) * x0 +
                                             GrainUnitVector(9 * MyOrientation + 4) * y0 +
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

                            double Angle3 = (GrainUnitVector(9 * MyOrientation + 6) * x0 +
                                             GrainUnitVector(9 * MyOrientation + 7) * y0 +
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
                            double ParaT = (normx * x0 + normy * y0 + normz * z0) /
                                           (normx * Diag1X + normy * Diag1Y + normz * Diag1Z);
                            float CDLVal = pow(
                                pow(ParaT * Diag1X, 2.0) + pow(ParaT * Diag1Y, 2.0) + pow(ParaT * Diag1Z, 2.0), 0.5);
                            CritDiagonalLength((long int)(26) * CellLocation + (long int)(n)) = CDLVal;
                        }
                        CellType(GlobalCellLocation) = Active;
                    }
                });
        }
    }

    // Wait on send requests
    MPI_Waitall(2, SendRequests.data(), MPI_STATUSES_IGNORE);
}
