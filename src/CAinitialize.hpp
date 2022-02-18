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
                       std::string &GrainOrientationFile, std::string &temppath, std::string &tempfile,
                       int &TempFilesInSeries, std::vector<std::string> &temp_paths, double &HT_deltax,
                       bool &RemeltingYN, double &deltat, int &NumberOfLayers, int &LayerHeight,
                       std::string &SubstrateFileName, float &SubstrateGrainSpacing, bool &UseSubstrateFile, double &G,
                       double &R, int &nx, int &ny, int &nz, double &FractSurfaceSitesActive, std::string &PathToOutput,
                       int &PrintDebug, bool &PrintMisorientation, bool &PrintFullOutput, int &NSpotsX, int &NSpotsY,
                       int &SpotOffset, int &SpotRadius, bool &PrintTimeSeries, int &TimeSeriesInc,
                       bool &PrintIdleTimeSeriesFrames);
void NeighborListInit(ViewI_H NeighborX, ViewI_H NeighborY, ViewI_H NeighborZ, ViewI2D_H ItList);
void FindXYZBounds(std::string SimulationType, int id, double &deltax, int &nx, int &ny, int &nz,
                   std::vector<std::string> &temp_paths, float &XMin, float &XMax, float &YMin, float &YMax,
                   float &ZMin, float &ZMax, int &LayerHeight, int NumberOfLayers, int TempFilesInSeries,
                   float *ZMinLayer, float *ZMaxLayer);
void DomainDecomposition(int DecompositionStrategy, int id, int np, int &MyXSlices, int &MyYSlices, int &MyXOffset,
                         int &MyYOffset, int &NeighborRank_North, int &NeighborRank_South, int &NeighborRank_East,
                         int &NeighborRank_West, int &NeighborRank_NorthEast, int &NeighborRank_NorthWest,
                         int &NeighborRank_SouthEast, int &NeighborRank_SouthWest, int &nx, int &ny, int &nz,
                         int &ProcessorsInXDirection, int &ProcessorsInYDirection, long int &LocalDomainSize);
void ReadTemperatureData(int id, double deltax, double HT_deltax, int MyXSlices, int MyYSlices, int MyXOffset,
                         int MyYOffset, float XMin, float YMin, std::vector<std::string> &temp_paths,
                         int NumberOfLayers, int TempFilesInSeries, unsigned int &NumberOfTemperatureDataPoints,
                         std::vector<double> &RawData, int *FirstValue, int *LastValue);
void TempInit_DirSolidification(double G, double R, int id, int &MyXSlices, int &MyYSlices, double deltax,
                                double deltat, int nz, ViewI_H CritTimeStep, ViewF_H UndercoolingChange,
                                ViewF_H UndercoolingCurrent, bool *Melted, int &nzActive, int &ZBound_Low,
                                int &ZBound_High, ViewI_H LayerID);
void TempInit_SpotMelt(double G, double R, std::string SimulationType, int id, int &MyXSlices, int &MyYSlices,
                       int &MyXOffset, int &MyYOffset, double deltax, double deltat, int nz, ViewI_H CritTimeStep,
                       ViewF_H UndercoolingChange, ViewF_H UndercoolingCurrent, bool *Melted, int LayerHeight,
                       int NumberOfLayers, int &nzActive, int &ZBound_Low, int &ZBound_High, double FreezingRange,
                       ViewI_H LayerID, int NSpotsX, int NSpotsY, int SpotRadius, int SpotOffset);
void TempInit_Reduced(int id, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset, double &deltax,
                      double HT_deltax, double deltat, int &nz, ViewI_H CritTimeStep, ViewF_H UndercoolingChange,
                      ViewF_H UndercoolingCurrent, float XMin, float YMin, float ZMin, bool *Melted, float *ZMinLayer,
                      float *ZMaxLayer, int LayerHeight, int NumberOfLayers, int &nzActive, int &ZBound_Low,
                      int &ZBound_High, int *FinishTimeStep, double FreezingRange, ViewI_H LayerID, int *FirstValue,
                      int *LastValue, std::vector<double> RawData);
void OrientationInit(int id, int NGrainOrientations, ViewI_H GrainOrientation, ViewF_H GrainUnitVector,
                     std::string GrainOrientationFile);
void SubstrateInit_ConstrainedGrowth(double FractSurfaceSitesActive, int MyXSlices, int MyYSlices, int nx, int ny,
                                     int nz, int MyXOffset, int MyYOffset, int pid, int np, ViewI_H CellType,
                                     ViewI_H GrainID);
void SubstrateInit_FromFile(std::string SubstrateFileName, int nz, int MyXSlices, int MyYSlices, int MyXOffset,
                            int MyYOffset, int pid, ViewI_H GrainID);
void SubstrateInit_FromGrainSpacing(float SubstrateGrainSpacing, int nx, int ny, int nz, int nzActive, int MyXSlices,
                                    int MyYSlices, int MyXOffset, int MyYOffset, int LocalActiveDomainSize, int pid,
                                    int np, double deltax, ViewI_H GrainID);
void ActiveCellInit(int id, int MyXSlices, int MyYSlices, int nz, ViewI_H CellType, ViewI_H CritTimeStep,
                    ViewI_H NeighborX, ViewI_H NeighborY, ViewI_H NeighborZ);
void GrainInit(int layernumber, int NGrainOrientations, int DecompositionStrategy, int nx, int ny, int nz,
               int LocalActiveDomainSize, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int id, int np,
               int NeighborRank_North, int NeighborRank_South, int NeighborRank_East, int NeighborRank_West,
               int NeighborRank_NorthEast, int NeighborRank_NorthWest, int NeighborRank_SouthEast,
               int NeighborRank_SouthWest, ViewI2D_H ItList, ViewI_H NeighborX, ViewI_H NeighborY, ViewI_H NeighborZ,
               ViewI_H GrainOrientation, ViewF_H GrainUnitVector, ViewF_H DiagonalLength, ViewI_H CellType,
               ViewI_H GrainID, ViewF_H CritDiagonalLength, ViewF_H DOCenter, double deltax, double NMax,
               int &NextLayer_FirstNucleatedGrainID, int &PossibleNuclei_ThisRank, int ZBound_High, int ZBound_Low);
void NucleiInit(int DecompositionStrategy, int MyXSlices, int MyYSlices, int nz, int id, double dTN, double dTsigma,
                int NeighborRank_North, int NeighborRank_South, int NeighborRank_East, int NeighborRank_West,
                int NeighborRank_NorthEast, int NeighborRank_NorthWest, int NeighborRank_SouthEast,
                int NeighborRank_SouthWest, ViewI_H NucleiLocation, ViewI_H NucleationTimes, ViewI_H CellType,
                ViewI_H GrainID, ViewI_H CritTimeStep, ViewF_H UndercoolingChange);
void DomainShiftAndResize(int id, int MyXSlices, int MyYSlices, int &ZShift, int &ZBound_Low, int &ZBound_High,
                          int &nzActive, int LocalDomainSize, int &LocalActiveDomainSize, int &BufSizeZ,
                          int LayerHeight, ViewI CellType, int layernumber, ViewI LayerID);
void LayerSetup(int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int LocalActiveDomainSize,
                ViewI GrainOrientation, int NGrainOrientations, ViewF GrainUnitVector, ViewI NeighborX, ViewI NeighborY,
                ViewI NeighborZ, ViewF DiagonalLength, ViewI CellType, ViewI GrainID, ViewF CritDiagonalLength,
                ViewF DOCenter, int DecompositionStrategy, Buffer2D BufferWestSend, Buffer2D BufferEastSend,
                Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, Buffer2D BufferNorthEastSend,
                Buffer2D BufferNorthWestSend, Buffer2D BufferSouthEastSend, Buffer2D BufferSouthWestSend,
                Buffer2D BufferWestRecv, Buffer2D BufferEastRecv, Buffer2D BufferNorthRecv, Buffer2D BufferSouthRecv,
                Buffer2D BufferNorthEastRecv, Buffer2D BufferNorthWestRecv, Buffer2D BufferSouthEastRecv,
                Buffer2D BufferSouthWestRecv, int &ZBound_Low);

void skipLines(std::ifstream &stream, std::string seperator = "*****");
std::string parseInput(std::ifstream &stream, std::string key);
bool parseInputFromList(std::string line, std::vector<std::string> Inputs, std::vector<std::string> &InputsRead,
                        int NumInputs);
bool getInputBool(std::string val);
std::string checkFileInstalled(const std::string name, const std::string type, const int id);
void checkFileNotEmpty(std::string testfilename);

#endif
