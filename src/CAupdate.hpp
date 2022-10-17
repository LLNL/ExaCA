// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_UPDATE_HPP
#define EXACA_UPDATE_HPP

#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

#include <string>

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

KOKKOS_INLINE_FUNCTION void loadghostnodes_1D(const int GhostGID, const float GhostDOCX, const float GhostDOCY,
                                              const float GhostDOCZ, const float GhostDL, const int BufSizeX,
                                              const int MyYSlices, const int RankX, const int RankY, const int RankZ,
                                              const bool AtNorthBoundary, const bool AtSouthBoundary,
                                              Buffer2D BufferSouthSend, Buffer2D BufferNorthSend);
KOKKOS_INLINE_FUNCTION void
loadghostnodes_2D(const int GhostGID, const float GhostDOCX, const float GhostDOCY, const float GhostDOCZ,
                  const float GhostDL, const int BufSizeX, const int BufSizeY, const int MyXSlices, const int MyYSlices,
                  const int RankX, const int RankY, const int RankZ, const bool AtNorthBoundary,
                  const bool AtSouthBoundary, const bool AtWestBoundary, const bool AtEastBoundary,
                  Buffer2D BufferSouthSend, Buffer2D BufferNorthSend, Buffer2D BufferWestSend, Buffer2D BufferEastSend,
                  Buffer2D BufferNorthEastSend, Buffer2D BufferSouthEastSend, Buffer2D BufferSouthWestSend,
                  Buffer2D BufferNorthWestSend);
void Nucleation(int cycle, int &SuccessfulNucEvents_ThisRank, int &NucleationCounter, int PossibleNuclei_ThisRank,
                ViewI_H NucleationTimes_H, ViewI NucleiLocations, ViewI NucleiGrainID, ViewI CellType, ViewI GrainID,
                int ZBound_Low, int MyXSlices, int MyYSlices, ViewI SteeringVector, ViewI numSteer_G);
void FillSteeringVector_NoRemelt(int cycle, int LocalActiveDomainSize, int MyXSlices, int MyYSlices, ViewI CritTimeStep,
                                 ViewF UndercoolingCurrent, ViewF UndercoolingChange, ViewI CellType, int ZBound_Low,
                                 int layernumber, ViewI LayerID, ViewI SteeringVector, ViewI numSteer_G,
                                 ViewI_H numSteer_H);
void FillSteeringVector_Remelt(int cycle, int LocalActiveDomainSize, int MyXSlices, int MyYSlices, NList NeighborX,
                               NList NeighborY, NList NeighborZ, ViewI CritTimeStep, ViewF UndercoolingCurrent,
                               ViewF UndercoolingChange, ViewI CellType, ViewI GrainID, int ZBound_Low, int nzActive,
                               ViewI SteeringVector, ViewI numSteer, ViewI_H numSteer_Host, ViewI MeltTimeStep,
                               int BufSizeX, int BufSizeY, bool AtNorthBoundary, bool AtSouthBoundary,
                               bool AtEastBoundary, bool AtWestBoundary, Buffer2D BufferWestSend,
                               Buffer2D BufferEastSend, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend,
                               Buffer2D BufferNorthEastSend, Buffer2D BufferNorthWestSend, Buffer2D BufferSouthEastSend,
                               Buffer2D BufferSouthWestSend, int DecompositionStrategy);
void CellCapture(int id, int np, int cycle, int DecompositionStrategy, int LocalActiveDomainSize, int LocalDomainSize,
                 int MyXSlices, int MyYSlices, double AConst, double BConst, double CConst, double DConst,
                 int MyXOffset, int MyYOffset, NList NeighborX, NList NeighborY, NList NeighborZ, ViewI CritTimeStep,
                 ViewF UndercoolingCurrent, ViewF UndercoolingChange, ViewF GrainUnitVector, ViewF CritDiagonalLength,
                 ViewF DiagonalLength, ViewI CellType, ViewF DOCenter, ViewI GrainID, int NGrainOrientations,
                 Buffer2D BufferWestSend, Buffer2D BufferEastSend, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend,
                 Buffer2D BufferNorthEastSend, Buffer2D BufferNorthWestSend, Buffer2D BufferSouthEastSend,
                 Buffer2D BufferSouthWestSend, int BufSizeX, int BufSizeY, int ZBound_Low, int nzActive, int nz,
                 ViewI SteeringVector, ViewI numSteer_G, ViewI_H numSteer_H, bool AtNorthBoundary, bool AtSouthBoundary,
                 bool AtEastBoundary, bool AtWestBoundary, ViewI SolidificationEventCounter, ViewI MeltTimeStep,
                 ViewF3D LayerTimeTempHistory, ViewI NumberOfSolidificationEvents, bool RemeltingYN);
void JumpTimeStep(int &cycle, unsigned long int RemainingCellsOfInterest, unsigned long int LocalIncompleteCells,
                  ViewI FutureWorkView, int LocalActiveDomainSize, int MyXSlices, int MyYSlices, int ZBound_Low,
                  bool RemeltingYN, ViewI CellType, ViewI LayerID, int id, int layernumber, int np, int nx, int ny,
                  int nz, int MyXOffset, int MyYOffset, int ProcessorsInXDirection, int ProcessorsInYDirection,
                  ViewI GrainID, ViewI CritTimeStep, ViewF GrainUnitVector, ViewF UndercoolingChange,
                  ViewF UndercoolingCurrent, std::string OutputFile, int DecompositionStrategy, int NGrainOrientations,
                  std::string PathToOutput, int &IntermediateFileCounter, int nzActive, double deltax, float XMin,
                  float YMin, float ZMin, int NumberOfLayers, int &XSwitch, std::string TemperatureDataType,
                  bool PrintIdleMovieFrames, int MovieFrameInc, int FinishTimeStep);
void IntermediateOutputAndCheck(int id, int np, int &cycle, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset,
                                int LocalDomainSize, int LocalActiveDomainSize, int nx, int ny, int nz, int nzActive,
                                double deltax, float XMin, float YMin, float ZMin, int DecompositionStrategy,
                                int ProcessorsInXDirection, int ProcessorsInYDirection,
                                int SuccessfulNucEvents_ThisRank, int &XSwitch, ViewI CellType, ViewI CritTimeStep,
                                ViewI GrainID, std::string TemperatureDataType, int *FinishTimeStep, int layernumber,
                                int, int ZBound_Low, int NGrainOrientations, ViewI LayerID, ViewF GrainUnitVector,
                                ViewF UndercoolingChange, ViewF UndercoolingCurrent, std::string PathToOutput,
                                std::string OutputFile, bool PrintIdleMovieFrames, int MovieFrameInc,
                                int &IntermediateFileCounter, int NumberOfLayers);
void IntermediateOutputAndCheck_Remelt(
    int id, int np, int &cycle, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int LocalActiveDomainSize,
    int nx, int ny, int nz, int nzActive, double deltax, float XMin, float YMin, float ZMin, int DecompositionStrategy,
    int ProcessorsInXDirection, int ProcessorsInYDirection, int SuccessfulNucEvents_ThisRank, int &XSwitch,
    ViewI CellType, ViewI CritTimeStep, ViewI GrainID, std::string TemperatureDataType, int layernumber, int,
    int ZBound_Low, int NGrainOrientations, ViewI LayerID, ViewF GrainUnitVector, ViewF UndercoolingChange,
    ViewF UndercoolingCurrent, std::string PathToOutput, std::string OutputFile, bool PrintIdleMovieFrames,
    int MovieFrameInc, int &IntermediateFileCounter, int NumberOfLayers, ViewI MeltTimeStep);

#endif
