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
                       bool &PrintDefaultRVE, double &RNGSeed, bool &BaseplateThroughPowder, double &PowderDensity,
                       int &RVESize);
void checkPowderOverflow(int nx, int ny, int LayerHeight, int NumberOfLayers, bool BaseplateThroughPowder,
                         double PowderDensity);
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
                         std::vector<double> &RawData, int *FirstValue, int *LastValue);
int calcZBound_Low(bool RemeltingYN, std::string SimulationType, int LayerHeight, int layernumber, float *ZMinLayer,
                   float ZMin, double deltax);
int calcZBound_High(std::string SimulationType, int SpotRadius, int LayerHeight, int layernumber, float ZMin,
                    double deltax, int nz, float *ZMaxLayer);
int calcnzActive(int ZBound_Low, int ZBound_High, int id, int layernumber);
int calcLocalActiveDomainSize(int MyXSlices, int MyYSlices, int nzActive);
void TempInit_DirSolidification(double G, double R, int id, int &MyXSlices, int &MyYSlices, double deltax,
                                double deltat, int nz, int LocalDomainSize, ViewI &CritTimeStep,
                                ViewF &UndercoolingChange, ViewI &LayerID);
int calcMaxSolidificationEventsSpot(int MyXSlices, int MyYSlices, int NumberOfSpots, int NSpotsX, int SpotRadius,
                                    int SpotOffset, int MyXOffset, int MyYOffset);
void OrientationInit(int id, int &NGrainOrientations, ViewF &ReadOrientationData, std::string GrainOrientationFile,
                     int ValsPerLine = 9);
void TempInit_SpotRemelt(int layernumber, double G, double R, std::string, int id, int &MyXSlices, int &MyYSlices,
                         int &MyXOffset, int &MyYOffset, double deltax, double deltat, int ZBound_Low, int nz,
                         int LocalActiveDomainSize, int LocalDomainSize, ViewI &CritTimeStep, ViewF &UndercoolingChange,
                         ViewF &UndercoolingCurrent, int LayerHeight, double FreezingRange, ViewI &LayerID, int NSpotsX,
                         int NSpotsY, int SpotRadius, int SpotOffset, ViewF3D &LayerTimeTempHistory,
                         ViewI &NumberOfSolidificationEvents, ViewI &MeltTimeStep, ViewI &MaxSolidificationEvents,
                         ViewI &SolidificationEventCounter);
void TempInit_SpotNoRemelt(double G, double R, std::string SimulationType, int id, int &MyXSlices, int &MyYSlices,
                           int &MyXOffset, int &MyYOffset, double deltax, double deltat, int nz, int LocalDomainSize,
                           ViewI &CritTimeStep, ViewF &UndercoolingChange, int LayerHeight, int NumberOfLayers,
                           double FreezingRange, ViewI &LayerID, int NSpotsX, int NSpotsY, int SpotRadius,
                           int SpotOffset);
int getTempCoordX(int i, float XMin, double deltax, const std::vector<double> &RawData);
int getTempCoordY(int i, float YMin, double deltax, const std::vector<double> &RawData);
int getTempCoordZ(int i, double deltax, const std::vector<double> &RawData, int LayerHeight, int LayerCounter,
                  float *ZMinLayer);
double getTempCoordTM(int i, const std::vector<double> &RawData);
double getTempCoordTL(int i, const std::vector<double> &RawData);
double getTempCoordCR(int i, const std::vector<double> &RawData);
void TempInit_ReadDataNoRemelt(int id, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset, double deltax,
                               int HTtoCAratio, double deltat, int nz, int LocalDomainSize, ViewI &CritTimeStep,
                               ViewF &UndercoolingChange, float XMin, float YMin, float ZMin, float *ZMinLayer,
                               float *ZMaxLayer, int LayerHeight, int NumberOfLayers, int *FinishTimeStep,
                               double FreezingRange, ViewI &LayerID, int *FirstValue, int *LastValue,
                               std::vector<double> RawData);
void calcMaxSolidificationEventsR(int id, int layernumber, int TempFilesInSeries, ViewI_H MaxSolidificationEvents_Host,
                                  int StartRange, int EndRange, std::vector<double> RawData, float XMin, float YMin,
                                  double deltax, float *ZMinLayer, int LayerHeight, int MyXSlices, int MyYSlices,
                                  int MyXOffset, int MyYOffset, int LocalActiveDomainSize);
void TempInit_ReadDataRemelt(int layernumber, int id, int MyXSlices, int MyYSlices, int nz, int LocalActiveDomainSize,
                             int LocalDomainSize, int MyXOffset, int MyYOffset, double &deltax, double deltat,
                             double FreezingRange, ViewF3D &LayerTimeTempHistory, ViewI &NumberOfSolidificationEvents,
                             ViewI &MaxSolidificationEvents, ViewI &MeltTimeStep, ViewI &CritTimeStep,
                             ViewF &UndercoolingChange, ViewF &UndercoolingCurrent, float XMin, float YMin,
                             float *ZMinLayer, int LayerHeight, int nzActive, int ZBound_Low, int *FinishTimeStep,
                             ViewI &LayerID, int *FirstValue, int *LastValue, std::vector<double> RawData,
                             ViewI &SolidificationEventCounter, int TempFilesInSeries);
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
                            int MyYOffset, int pid, ViewI &GrainID, int nzActive, bool BaseplateThroughPowder);
void BaseplateInit_FromGrainSpacing(float SubstrateGrainSpacing, int nx, int ny, float *ZMinLayer, float *ZMaxLayer,
                                    int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int id, double deltax,
                                    ViewI GrainID, double RNGSeed, int &NextLayer_FirstEpitaxialGrainID, int nz,
                                    double BaseplateThroughPowder);
void PowderInit(int layernumber, int nx, int ny, int LayerHeight, float *ZMaxLayer, float ZMin, double deltax,
                int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int id, ViewI GrainID, double RNGSeed,
                int &NextLayer_FirstEpitaxialGrainID, double PowderDensity);
void CellTypeInit_Remelt(int MyXSlices, int MyYSlices, int LocalActiveDomainSize, ViewI CellType, ViewI CritTimeStep,
                         int id, int ZBound_Low);
void CellTypeInit_NoRemelt(int layernumber, int id, int np, int DecompositionStrategy, int MyXSlices, int MyYSlices,
                           int MyXOffset, int MyYOffset, int ZBound_Low, int nz, int LocalActiveDomainSize,
                           int LocalDomainSize, ViewI CellType, ViewI CritTimeStep, NList NeighborX, NList NeighborY,
                           NList NeighborZ, int NGrainOrientations, ViewF GrainUnitVector, ViewF DiagonalLength,
                           ViewI GrainID, ViewF CritDiagonalLength, ViewF DOCenter, ViewI LayerID,
                           Buffer2D BufferWestSend, Buffer2D BufferEastSend, Buffer2D BufferNorthSend,
                           Buffer2D BufferSouthSend, Buffer2D BufferNorthEastSend, Buffer2D BufferNorthWestSend,
                           Buffer2D BufferSouthEastSend, Buffer2D BufferSouthWestSend, int BufSizeX, int BufSizeY,
                           bool AtNorthBoundary, bool AtSouthBoundary, bool AtEastBoundary, bool AtWestBoundary);
void NucleiInit(int layernumber, double RNGSeed, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nx,
                int ny, int nzActive, int ZBound_Low, int id, double NMax, double dTN, double dTsigma, double deltax,
                ViewI &NucleiLocation, ViewI_H &NucleationTimes_Host, ViewI &NucleiGrainID, ViewI CellType,
                ViewI CritTimeStep, ViewF UndercoolingChange, ViewI LayerID, int &PossibleNuclei_ThisRankThisLayer,
                int &Nuclei_WholeDomain, bool AtNorthBoundary, bool AtSouthBoundary, bool AtEastBoundary,
                bool AtWestBoundary, bool RemeltingYN, int &NucleationCounter, ViewI &MaxSolidificationEvents,
                ViewI NumberOfSolidificationEvents, ViewF3D LayerTimeTempHistory);
void placeNucleiData_NoRemelt(int Nuclei_ThisLayerSingle, ViewI_H NucleiX, ViewI_H NucleiY, ViewI_H NucleiZ,
                              int MyXOffset, int MyYOffset, int MyXSlices, int MyYSlices, bool AtNorthBoundary,
                              bool AtSouthBoundary, bool AtWestBoundary, bool AtEastBoundary, int ZBound_Low,
                              ViewI_H CellType_Host, ViewI_H LayerID_Host, ViewI_H CritTimeStep_Host,
                              ViewF_H UndercoolingChange_Host, int layernumber,
                              std::vector<int> NucleiGrainID_WholeDomain_V,
                              std::vector<double> NucleiUndercooling_WholeDomain_V,
                              std::vector<int> &NucleiGrainID_MyRank_V, std::vector<int> &NucleiLocation_MyRank_V,
                              std::vector<int> &NucleationTimes_MyRank_V, int &PossibleNuclei_ThisRankThisLayer);
void placeNucleiData_Remelt(int NucleiMultiplier, int Nuclei_ThisLayerSingle, ViewI_H NucleiX, ViewI_H NucleiY,
                            ViewI_H NucleiZ, int MyXOffset, int MyYOffset, int MyXSlices, int MyYSlices,
                            bool AtNorthBoundary, bool AtSouthBoundary, bool AtWestBoundary, bool AtEastBoundary,
                            int ZBound_Low, ViewI_H NumberOfSolidificationEvents_Host,
                            ViewF3D_H LayerTimeTempHistory_Host, std::vector<int> NucleiGrainID_WholeDomain_V,
                            std::vector<double> NucleiUndercooling_WholeDomain_V,
                            std::vector<int> &NucleiGrainID_MyRank_V, std::vector<int> &NucleiLocation_MyRank_V,
                            std::vector<int> &NucleationTimes_MyRank_V, int &PossibleNuclei_ThisRankThisLayer);
void DomainShiftAndResize_NoRemelt(int id, int MyXSlices, int MyYSlices, int &ZShift, int &ZBound_Low, int &ZBound_High,
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

#endif
