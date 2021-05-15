// Copyright 2021 Lawrence Livermore National Security, LLC and other
// ExaCA Project Developers. See the top-level LICENSE file for details.
// SPDX-License-Identifier: MIT

#ifndef EXACA_UPDATE_HPP
#define EXACA_UPDATE_HPP

#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

#include <string>

void Nucleation(int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int cycle, int &nn, ViewI CellType,
                ViewI NucleiLocations, ViewI NucleationTimes, ViewI GrainID, ViewI GrainOrientation, ViewF DOCenter,
                ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewF GrainUnitVector, ViewF CritDiagonalLength,
                ViewF DiagonalLength, int NGrainOrientations, int PossibleNuclei_ThisRank, int ZBound_Low,
                int layernumber, ViewI LayerID);
void CellCapture(int np, int cycle, int DecompositionStrategy, int LocalActiveDomainSize, int LocalDomainSize,
                 int MyXSlices, int MyYSlices, double AConst, double BConst, double CConst, double DConst,
                 int MyXOffset, int MyYOffset, ViewI2D ItList, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ,
                 ViewI CritTimeStep, ViewF UndercoolingCurrent, ViewF UndercoolingChange, ViewF GrainUnitVector,
                 ViewF CritDiagonalLength, ViewF DiagonalLength, ViewI GrainOrientation, ViewI CellType, ViewF DOCenter,
                 ViewI GrainID, int NGrainOrientations, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferC,
                 Buffer2D BufferD, Buffer2D BufferE, Buffer2D BufferF, Buffer2D BufferG, Buffer2D BufferH, int BufSizeX,
                 int BufSizeY, int ZBound_Low, int nzActive, int nz, int layernumber, ViewI LayerID,
                 ViewI SteeringVector, ViewI numSteer_G, ViewI_H numSteer_H);
void IntermediateOutputAndCheck(int pid, int &cycle, int MyXSlices, int MyYSlices, int LocalDomainSize,
                                int LocalActiveDomainSize, int nn, int &XSwitch, ViewI CellType, ViewI CritTimeStep,
                                std::string SimulationType, int *FinishTimeStep, int layernumber, int NumberOfLayers,
                                int ZBound_Low, ViewI LayerID);

#endif
