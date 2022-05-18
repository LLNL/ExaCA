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
                       std::string &GrainOrientationFile, int &TempFilesInSeries, std::vector<std::string> &temp_paths,
                       double &HT_deltax, bool &RemeltingYN, double &deltat, int &NumberOfLayers, int &LayerHeight,
                       std::string &SubstrateFileName, float &SubstrateGrainSpacing, bool &UseSubstrateFile, double &G,
                       double &R, int &nx, int &ny, int &nz, double &FractSurfaceSitesActive, std::string &PathToOutput,
                       int &PrintDebug, bool &PrintMisorientation, bool &PrintFinalUndercoolingVals,
                       bool &PrintFullOutput, int &NSpotsX, int &NSpotsY, int &SpotOffset, int &SpotRadius,
                       bool &PrintTimeSeries, int &TimeSeriesInc, bool &PrintIdleTimeSeriesFrames,
                       bool &PrintDefaultRVE, double &RNGSeed);
void NeighborListInit(NList &NeighborX, NList &NeighborY, NList &NeighborZ);
void FindXYZBounds(std::string SimulationType, int id, double &deltax, int &nx, int &ny, int &nz,
                   std::vector<std::string> &temp_paths, float &XMin, float &XMax, float &YMin, float &YMax,
                   float &ZMin, float &ZMax, int &LayerHeight, int NumberOfLayers, int TempFilesInSeries,
                   float *ZMinLayer, float *ZMaxLayer, int SpotRadius);
void DomainDecomposition(int DecompositionStrategy, int id, int np, int &MyXSlices, int &MyYSlices, int &MyXOffset,
                         int &MyYOffset, int &NeighborRank_North, int &NeighborRank_South, int &NeighborRank_East,
                         int &NeighborRank_West, int &NeighborRank_NorthEast, int &NeighborRank_NorthWest,
                         int &NeighborRank_SouthEast, int &NeighborRank_SouthWest, int &nx, int &ny, int &nz,
                         int &ProcessorsInXDirection, int &ProcessorsInYDirection, long int &LocalDomainSize,
                         bool &AtNorthBoundary, bool &AtSouthBoundary, bool &AtEastBoundary, bool &AtWestBoundary);
void ReadTemperatureData(int id, double &deltax, double HT_deltax, int &HTtoCAratio, int MyXSlices, int MyYSlices,
                         int MyXOffset, int MyYOffset, float XMin, float YMin, std::vector<std::string> &temp_paths,
                         int NumberOfLayers, int TempFilesInSeries, unsigned int &NumberOfTemperatureDataPoints,
                         std::vector<double> &RawData, int *FirstValue, int *LastValue, bool RemeltingYN,
                         ViewI &MaxSolidificationEvents, float *ZMinLayer, float *ZMaxLayer, int LayerHeight);
int calcZBound_Low(bool RemeltingYN, std::string SimulationType, int LayerHeight, int layernumber, float *ZMinLayer,
                   float ZMin, double deltax);
int calcZBound_High(std::string SimulationType, int SpotRadius, int LayerHeight, int layernumber, float ZMin,
                    double deltax, int nz, float *ZMaxLayer);
int calcnzActive(int ZBound_Low, int ZBound_High, int id, int layernumber);
int calcLocalActiveDomainSize(int MyXSlices, int MyYSlices, int nzActive);
void TempInit_DirSolidification(double G, double R, int id, int &MyXSlices, int &MyYSlices, double deltax,
                                double deltat, int nz, int LocalDomainSize, ViewI &CritTimeStep,
                                ViewF &UndercoolingChange, bool *Melted, ViewI &LayerID,
                                ViewI &MaxSolidificationEvents);
int calcMaxSolidificationEventsSpot(int MyXSlices, int MyYSlices, int NumberOfSpots, int NSpotsX, int SpotRadius,
                                    int SpotOffset, int MyXOffset, int MyYOffset);
void TempInit_SpotMeltRemelting(int layernumber, double G, double R, std::string, int id, int &MyXSlices,
                                int &MyYSlices, int &MyXOffset, int &MyYOffset, double deltax, double deltat,
                                int ZBound_Low, int nz, int LocalActiveDomainSize, int LocalDomainSize,
                                ViewI &CritTimeStep, ViewF &UndercoolingChange, ViewF &UndercoolingCurrent,
                                bool *Melted, int LayerHeight, double FreezingRange, ViewI &LayerID, int NSpotsX,
                                int NSpotsY, int SpotRadius, int SpotOffset, ViewF3D &LayerTimeTempHistory,
                                ViewI &NumberOfSolidificationEvents, ViewI &MeltTimeStep,
                                ViewI &MaxSolidificationEvents, ViewI &SolidificationEventCounter);
void TempInit_SpotMelt(double G, double R, std::string SimulationType, int id, int &MyXSlices, int &MyYSlices,
                       int &MyXOffset, int &MyYOffset, double deltax, double deltat, int nz, int LocalDomainSize,
                       ViewI &CritTimeStep, ViewF &UndercoolingChange, bool *Melted, int LayerHeight,
                       int NumberOfLayers, double FreezingRange, ViewI &LayerID, int NSpotsX, int NSpotsY,
                       int SpotRadius, int SpotOffset, ViewI &MaxSolidificationEvents);
void TempInit_Reduced(int id, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset, double deltax,
                      int HTtoCAratio, double deltat, int nz, int LocalDomainSize, ViewI &CritTimeStep,
                      ViewF &UndercoolingChange, float XMin, float YMin, float ZMin, bool *Melted, float *ZMinLayer,
                      float *ZMaxLayer, int LayerHeight, int NumberOfLayers, int *FinishTimeStep, double FreezingRange,
                      ViewI &LayerID, int *FirstValue, int *LastValue, std::vector<double> RawData);
void OrientationInit(int id, int &NGrainOrientations, ViewF &ReadOrientationData, std::string GrainOrientationFile);
void TempInit_Remelt(int layernumber, int id, int MyXSlices, int MyYSlices, int nz, int LocalActiveDomainSize,
                     int LocalDomainSize, int MyXOffset, int MyYOffset, double &deltax, double deltat,
                     double FreezingRange, ViewF3D &LayerTimeTempHistory, ViewI &NumberOfSolidificationEvents,
                     ViewI MaxSolidificationEvents, ViewI &MeltTimeStep, ViewI &CritTimeStep, ViewF &UndercoolingChange,
                     ViewF &UndercoolingCurrent, float XMin, float YMin, bool *Melted, float *ZMinLayer,
                     int LayerHeight, int nzActive, int ZBound_Low, int *FinishTimeStep, ViewI &LayerID,
                     int *FirstValue, int *LastValue, std::vector<double> RawData, ViewI &SolidificationEventCounter);
void SubstrateInit_ConstrainedGrowth(int id, double FractSurfaceSitesActive, int MyXSlices, int MyYSlices, int nx,
                                     int ny, int MyXOffset, int MyYOffset, NList NeighborX, NList NeighborY,
                                     NList NeighborZ, ViewF GrainUnitVector, int NGrainOrientations, ViewI CellType,
                                     ViewI GrainID, ViewF DiagonalLength, ViewF DOCenter, ViewF CritDiagonalLength,
                                     double RNGSeed, int np, int DecompositionStrategy, Buffer2D BufferWestSend,
                                     Buffer2D BufferEastSend, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend,
                                     Buffer2D BufferNorthEastSend, Buffer2D BufferNorthWestSend,
                                     Buffer2D BufferSouthEastSend, Buffer2D BufferSouthWestSend, int BufSizeX,
                                     int BufSizeY, bool AtNorthBoundary, bool AtSouthBoundary, bool AtEastBoundary,
                                     bool AtWestBoundary);
void SubstrateInit_FromFile(std::string SubstrateFileName, int nz, int MyXSlices, int MyYSlices, int MyXOffset,
                            int MyYOffset, int pid, ViewI &GrainID);
void BaseplateInit_FromGrainSpacing(float SubstrateGrainSpacing, int nx, int ny, float *ZMinLayer, float *ZMaxLayer,
                                    int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int id, double deltax,
                                    ViewI GrainID, double RNGSeed, int &NextLayer_FirstEpitaxialGrainID);
void PowderInit(int layernumber, int nx, int ny, int LayerHeight, float *ZMaxLayer, float ZMin, double deltax,
                int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int id, ViewI GrainID, double RNGSeed,
                int &NextLayer_FirstEpitaxialGrainID);
void CellTypeInit_Remelt(int MyXSlices, int MyYSlices, int LocalActiveDomainSize, ViewI CellType, ViewI CritTimeStep,
                         int id, int ZBound_Low);
void CellTypeInit(int layernumber, int id, int np, int DecompositionStrategy, int MyXSlices, int MyYSlices,
                  int MyXOffset, int MyYOffset, int ZBound_Low, int nz, int LocalActiveDomainSize, int LocalDomainSize,
                  ViewI CellType, ViewI CritTimeStep, NList NeighborX, NList NeighborY, NList NeighborZ,
                  int NGrainOrientations, ViewF GrainUnitVector, ViewF DiagonalLength, ViewI GrainID,
                  ViewF CritDiagonalLength, ViewF DOCenter, ViewI LayerID, Buffer2D BufferWestSend,
                  Buffer2D BufferEastSend, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend,
                  Buffer2D BufferNorthEastSend, Buffer2D BufferNorthWestSend, Buffer2D BufferSouthEastSend,
                  Buffer2D BufferSouthWestSend, int BufSizeX, int BufSizeY, bool AtNorthBoundary, bool AtSouthBoundary,
                  bool AtEastBoundary, bool AtWestBoundary);
void NucleiInit(int layernumber, double RNGSeed, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nx,
                int ny, int nzActive, int ZBound_Low, int id, double NMax, double dTN, double dTsigma, double deltax,
                ViewI &NucleiLocation, ViewI_H &NucleationTimes_Host, ViewI &NucleiGrainID, ViewI CellType,
                ViewI CritTimeStep, ViewF UndercoolingChange, ViewI LayerID, int &PossibleNuclei_ThisRankThisLayer,
                int &Nuclei_WholeDomain, bool AtNorthBoundary, bool AtSouthBoundary, bool AtEastBoundary,
                bool AtWestBoundary, bool RemeltingYN, int &NucleationCounter, ViewI &MaxSolidificationEvents,
                ViewI NumberOfSolidificationEvents, ViewF3D LayerTimeTempHistory);
void DomainShiftAndResize(int id, int MyXSlices, int MyYSlices, int &ZShift, int &ZBound_Low, int &ZBound_High,
                          int &nzActive, int LocalDomainSize, int &LocalActiveDomainSize, int &BufSizeZ,
                          int LayerHeight, ViewI CellType, int layernumber, ViewI LayerID);
void ZeroResetViews(int LocalActiveDomainSize, int BufSizeX, int BufSizeY, int BufSizeZ, ViewF &DiagonalLength,
                    ViewF &CritDiagonalLength, ViewF &DOCenter, int DecompositionStrategy, Buffer2D &BufferWestSend,
                    Buffer2D &BufferEastSend, Buffer2D &BufferNorthSend, Buffer2D &BufferSouthSend,
                    Buffer2D &BufferNorthEastSend, Buffer2D &BufferNorthWestSend, Buffer2D &BufferSouthEastSend,
                    Buffer2D &BufferSouthWestSend, Buffer2D &BufferWestRecv, Buffer2D &BufferEastRecv,
                    Buffer2D &BufferNorthRecv, Buffer2D &BufferSouthRecv, Buffer2D &BufferNorthEastRecv,
                    Buffer2D &BufferNorthWestRecv, Buffer2D &BufferSouthEastRecv, Buffer2D &BufferSouthWestRecv,
                    ViewI &SteeringVector);

void skipLines(std::ifstream &stream, std::string seperator = "*****");
std::string parseInput(std::ifstream &stream, std::string key);
bool parseInputFromList(std::string line, std::vector<std::string> Inputs, std::vector<std::string> &InputsRead,
                        int NumInputs);
bool getInputBool(std::string val);
std::string checkFileInstalled(const std::string name, const std::string type, const int id);
void checkFileNotEmpty(std::string testfilename);

#endif
