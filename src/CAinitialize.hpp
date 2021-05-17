// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_INIT_HPP
#define EXACA_INIT_HPP

#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

#include <string>
#include <vector>

void InputReadFromFile(int id, std::string InputFile, std::string &SimulationType, int &DecompositionStrategy,
                       double &AConst, double &BConst, double &CConst, double &DConst, double &FreezingRange,
                       double &deltax, double &NMax, double &dTN, double &dTsigma, std::string &OutputFile,
                       std::string &GrainOrientationFile, std::string &tempfile, int &TempFilesInSeries,
                       bool &ExtraWalls, double &HT_deltax, bool &RemeltingYN, double &deltat, int &NumberOfLayers,
                       int &LayerHeight, std::string &SubstrateFileName, float &SubstrateGrainSpacing,
                       bool &UseSubstrateFile, double &G, double &R, int &nx, int &ny, int &nz,
                       double &FractSurfaceSitesActive, std::string &PathToOutput, bool (&FilesToPrint)[6],
                       bool &PrintFilesYN);
void ParallelMeshInit(int DecompositionStrategy, ViewI_H NeighborX, ViewI_H NeighborY, ViewI_H NeighborZ,
                      ViewI2D_H ItList, std::string SimulationType, int id, int np, int &MyXSlices, int &MyYSlices,
                      int &MyXOffset, int &MyYOffset, int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn,
                      int &MyLeftOut, int &MyRightIn, int &MyRightOut, double &deltax, double HT_deltax, int &nx,
                      int &ny, int &nz, int &ProcessorsInXDirection, int &ProcessorsInYDirection, std::string tempfile,
                      float &XMin, float &XMax, float &YMin, float &YMax, float &ZMin, float &ZMax, float FreezingRange,
                      int &LayerHeight, int NumberOfLayers, int TempFilesInSeries,
                      unsigned int &NumberOfTemperatureDataPoints, float *ZMinLayer, float *ZMaxLayer, int *FirstValue,
                      int *LastValue, std::vector<float> &RawData);
void TempInit(int layernumber, double G, double R, std::string SimulationType, int id, int &MyXSlices, int &MyYSlices,
              int &MyXOffset, int &MyYOffset, double &deltax, double HT_deltax, double deltat, int &nx, int &ny,
              int &nz, ViewI_H CritTimeStep, ViewF_H UndercoolingChange, ViewF_H UndercoolingCurrent, float XMin,
              float YMin, float ZMin, bool *Melted, float *ZMinLayer, float *ZMaxLayer, int LayerHeight,
              int NumberOfLayers, int &nzActive, int &ZBound_Low, int &ZBound_High, int *FinishTimeStep,
              double FreezingRange, ViewI_H LayerID, int *FirstValue, int *LastValue, std::vector<float> RawData);
void OrientationInit(int id, int NGrainOrientations, ViewI_H GrainOrientation, ViewF_H GrainUnitVector,
                     std::string GrainOrientationFile);
void SubstrateInit_ConstrainedGrowth(double FractSurfaceSitesActive, int MyXSlices, int MyYSlices, int nx, int ny,
                                     int nz, int MyXOffset, int MyYOffset, int pid, int np, ViewI_H CellType,
                                     ViewI_H GrainID);
void SubstrateInit_FromFile(std::string SubstrateFileName, int nz, int MyXSlices, int MyYSlices, int MyXOffset,
                            int MyYOffset, int pid, ViewI_H CritTimeStep, ViewI_H GrainID);
void SubstrateInit_FromGrainSpacing(float SubstrateGrainSpacing, int nx, int ny, int nz, int nzActive, int MyXSlices,
                                    int MyYSlices, int MyXOffset, int MyYOffset, int LocalActiveDomainSize, int pid,
                                    int np, double deltax, ViewI_H GrainID, ViewI_H CritTimeStep);
void ActiveCellWallInit(int id, int MyXSlices, int MyYSlices, int nx, int ny, int nz, int MyXOffset, int MyYOffset,
                        ViewI_H CellType_H, ViewI_H GrainID_H, ViewI_H CritTimeStep_H, ViewI2D_H ItList_H,
                        ViewI_H NeighborX_H, ViewI_H NeighborY_H, ViewI_H NeighborZ_H, ViewF_H UndercoolingChange_H,
                        bool ExtraWalls);
void GrainInit(int layernumber, int NGrainOrientations, int DecompositionStrategy, int nz, int LocalActiveDomainSize,
               int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int id, int np, int MyLeft, int MyRight,
               int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, ViewI2D_H ItList,
               ViewI_H NeighborX, ViewI_H NeighborY, ViewI_H NeighborZ, ViewI_H GrainOrientation,
               ViewF_H GrainUnitVector, ViewF_H DiagonalLength, ViewI_H CellType, ViewI_H GrainID,
               ViewF_H CritDiagonalLength, ViewF_H DOCenter, ViewI_H CritTimeStep, double deltax, double NMax,
               int &NextLayer_FirstNucleatedGrainID, int &PossibleNuclei_ThisRank, int ZBound_High, int ZBound_Low);
void NucleiInit(int DecompositionStrategy, int MyXSlices, int MyYSlices, int nz, int id, double dTN, double dTsigma,
                int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut,
                int MyRightOut, ViewI_H NucleiLocation, ViewI_H NucleationTimes, ViewI_H CellType, ViewI_H GrainID,
                ViewI_H CritTimeStep, ViewF_H UndercoolingChange);
void DomainShiftAndResize(int id, int MyXSlices, int MyYSlices, int &ZShift, int &ZBound_Low, int &ZBound_High,
                          int &nzActive, int LocalDomainSize, int &LocalActiveDomainSize, int &BufSizeZ,
                          int LayerHeight, ViewI CellType, int layernumber, ViewI LayerID);
void LayerSetup(int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int LocalActiveDomainSize,
                ViewI GrainOrientation, int NGrainOrientations, ViewF GrainUnitVector, ViewI NeighborX, ViewI NeighborY,
                ViewI NeighborZ, ViewF DiagonalLength, ViewI CellType, ViewI GrainID, ViewF CritDiagonalLength,
                ViewF DOCenter, int DecompositionStrategy, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferC,
                Buffer2D BufferD, Buffer2D BufferE, Buffer2D BufferF, Buffer2D BufferG, Buffer2D BufferH,
                Buffer2D BufferAR, Buffer2D BufferBR, Buffer2D BufferCR, Buffer2D BufferDR, Buffer2D BufferER,
                Buffer2D BufferFR, Buffer2D BufferGR, Buffer2D BufferHR, int &ZBound_Low);

#endif
