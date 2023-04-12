// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_UPDATE_HPP
#define EXACA_UPDATE_HPP

#include "CAinterfacialresponse.hpp"
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
                               Buffer2D BufferSouthSend, int NGrainOrientations);
void CellCapture(int id, int np, int cycle, int LocalActiveDomainSize, int LocalDomainSize, int nx, int MyYSlices,
                 InterfacialResponseFunction irf, int MyYOffset, NList NeighborX, NList NeighborY, NList NeighborZ,
                 ViewI CritTimeStep, ViewF UndercoolingCurrent, ViewF UndercoolingChange, ViewF GrainUnitVector,
                 ViewF CritDiagonalLength, ViewF DiagonalLength, ViewI CellType, ViewF DOCenter, ViewI GrainID,
                 int NGrainOrientations, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, int BufSizeX,
                 int ZBound_Low, int nzActive, int nz, ViewI SteeringVector, ViewI numSteer_G, ViewI_H numSteer_H,
                 bool AtNorthBoundary, bool AtSouthBoundary, ViewI SolidificationEventCounter, ViewI MeltTimeStep,
                 ViewF3D LayerTimeTempHistory, ViewI NumberOfSolidificationEvents, bool RemeltingYN);
void JumpTimeStep(int &cycle, unsigned long int RemainingCellsOfInterest, ViewI FutureWorkView,
                  unsigned long int LocalIncompleteCells, int LocalActiveDomainSize, int MyYSlices, int ZBound_Low,
                  bool RemeltingYN, ViewI CellType, ViewI LayerID, int id, int layernumber, int np, int nx, int ny,
                  int nz, int MyYOffset, ViewI GrainID, ViewI CritTimeStep, ViewF GrainUnitVector,
                  ViewF UndercoolingChange, ViewF UndercoolingCurrent, std::string OutputFile,
                  int DecompositionStrategy, int NGrainOrientations, std::string PathToOutput,
                  int &IntermediateFileCounter, int nzActive, double deltax, double XMin, double YMin, double ZMin,
                  int NumberOfLayers, int &XSwitch, std::string TemperatureDataType, bool PrintIdleMovieFrames,
                  int MovieFrameInc, bool PrintBinary, int FinishTimeStep);
void IntermediateOutputAndCheck(int id, int np, int &cycle, int MyYSlices, int MyYOffset, int LocalDomainSize,
                                int LocalActiveDomainSize, int nx, int ny, int nz, int nzActive, double deltax,
                                double XMin, double YMin, double ZMin, int SuccessfulNucEvents_ThisRank, int &XSwitch,
                                ViewI CellType, ViewI CritTimeStep, ViewI GrainID, std::string TemperatureDataType,
                                int *FinishTimeStep, int layernumber, int, int ZBound_Low, int NGrainOrientations,
                                ViewI LayerID, ViewF GrainUnitVector, ViewF UndercoolingChange,
                                ViewF UndercoolingCurrent, std::string PathToOutput, std::string OutputFile,
                                bool PrintIdleMovieFrames, int MovieFrameInc, int &IntermediateFileCounter,
                                int NumberOfLayers, bool PrintBinary);
void IntermediateOutputAndCheck_Remelt(
    int id, int np, int &cycle, int MyYSlices, int MyYOffset, int LocalActiveDomainSize, int nx, int ny, int nz,
    int nzActive, double deltax, double XMin, double YMin, double ZMin, int SuccessfulNucEvents_ThisRank, int &XSwitch,
    ViewI CellType, ViewI CritTimeStep, ViewI GrainID, std::string TemperatureDataType, int layernumber, int,
    int ZBound_Low, int NGrainOrientations, ViewI LayerID, ViewF GrainUnitVector, ViewF UndercoolingChange,
    ViewF UndercoolingCurrent, std::string PathToOutput, std::string OutputFile, bool PrintIdleMovieFrames,
    int MovieFrameInc, int &IntermediateFileCounter, int NumberOfLayers, ViewI MeltTimeStep, bool PrintBinary);

#endif
