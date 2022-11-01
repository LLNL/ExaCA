// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_INIT_HPP
#define EXACA_INIT_HPP

#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

#include <string>
#include <vector>

void InputReadFromFile(int id, std::string InputFile, std::string &SimulationType, double &deltax, double &NMax,
                       double &dTN, double &dTsigma, std::string &OutputFile, std::string &GrainOrientationFile,
                       int &TempFilesInSeries, std::vector<std::string> &temp_paths, double &HT_deltax,
                       bool &RemeltingYN, double &deltat, int &NumberOfLayers, int &LayerHeight,
                       std::string &MaterialFileName, std::string &SubstrateFileName, float &SubstrateGrainSpacing,
                       bool &UseSubstrateFile, double &G, double &R, int &nx, int &ny, int &nz,
                       double &FractSurfaceSitesActive, std::string &PathToOutput, int &PrintDebug,
                       bool &PrintMisorientation, bool &PrintFinalUndercoolingVals, bool &PrintFullOutput, int &NSpotsX,
                       int &NSpotsY, int &SpotOffset, int &SpotRadius, bool &PrintTimeSeries, int &TimeSeriesInc,
                       bool &PrintIdleTimeSeriesFrames, bool &PrintDefaultRVE, double &RNGSeed,
                       bool &BaseplateThroughPowder, double &PowderDensity, int &RVESize, bool &LayerwiseTempRead,
                       bool &PrintBinary);
void checkPowderOverflow(int nx, int ny, int LayerHeight, int NumberOfLayers, bool BaseplateThroughPowder,
                         double PowderDensity);
void NeighborListInit(NList &NeighborX, NList &NeighborY, NList &NeighborZ);
void FindXYZBounds(std::string SimulationType, int id, double &deltax, int &nx, int &ny, int &nz,
                   std::vector<std::string> &temp_paths, float &XMin, float &XMax, float &YMin, float &YMax,
                   float &ZMin, float &ZMax, int &LayerHeight, int NumberOfLayers, int TempFilesInSeries,
                   float *ZMinLayer, float *ZMaxLayer, int SpotRadius);
void DomainDecomposition(int id, int np, int &MyYSlices, int &MyYOffset, int &NeighborRank_North,
                         int &NeighborRank_South, int &nx, int &ny, int &nz, long int &LocalDomainSize,
                         bool &AtNorthBoundary, bool &AtSouthBoundary);
void ReadTemperatureData(int id, double &deltax, double HT_deltax, int &HTtoCAratio, int MyYSlices, int MyYOffset,
                         float YMin, std::vector<std::string> &temp_paths, int NumberOfLayers, int TempFilesInSeries,
                         unsigned int &NumberOfTemperatureDataPoints, std::vector<double> &RawData, int *FirstValue,
                         int *LastValue, bool LayerwiseTempRead, int layernumber);
int calcZBound_Low_Remelt(std::string SimulationType, int LayerHeight, int layernumber, float *ZMinLayer, float ZMin,
                          double deltax);
int calcZBound_High_Remelt(std::string SimulationType, int SpotRadius, int LayerHeight, int layernumber, float ZMin,
                           double deltax, int nz, float *ZMaxLayer);
int calcZBound_Low_NoRemelt(int id, int MyXSlices, int MyYSlices, int LocalDomainSize, int layernumber, ViewI LayerID);
int calcZBound_High_NoRemelt(int id, int MyXSlices, int MyYSlices, int LocalDomainSize, int layernumber, ViewI LayerID);
int calcnzActive(int ZBound_Low, int ZBound_High, int id, int layernumber);
int calcLocalActiveDomainSize(int nx, int MyYSlices, int nzActive);
void TempInit_DirSolidification(double G, double R, int id, int &nx, int &MyYSlices, double deltax, double deltat,
                                int nz, int LocalDomainSize, ViewI &CritTimeStep, ViewF &UndercoolingChange,
                                ViewI &LayerID);
int calcMaxSolidificationEventsSpot(int nx, int MyYSlices, int NumberOfSpots, int NSpotsX, int SpotRadius,
                                    int SpotOffset, int MyYOffset);
void OrientationInit(int id, int &NGrainOrientations, ViewF &ReadOrientationData, std::string GrainOrientationFile,
                     int ValsPerLine = 9);
void TempInit_SpotRemelt(int layernumber, double G, double R, std::string, int id, int &nx, int &MyYSlices,
                         int &MyYOffset, double deltax, double deltat, int ZBound_Low, int nz,
                         int LocalActiveDomainSize, int LocalDomainSize, ViewI &CritTimeStep, ViewF &UndercoolingChange,
                         ViewF &UndercoolingCurrent, int LayerHeight, double FreezingRange, ViewI &LayerID, int NSpotsX,
                         int NSpotsY, int SpotRadius, int SpotOffset, ViewF3D &LayerTimeTempHistory,
                         ViewI &NumberOfSolidificationEvents, ViewI &MeltTimeStep, ViewI &MaxSolidificationEvents,
                         ViewI &SolidificationEventCounter);
void TempInit_SpotNoRemelt(double G, double R, std::string SimulationType, int id, int &nx, int &MyYSlices,
                           int &MyYOffset, double deltax, double deltat, int nz, int LocalDomainSize,
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
void TempInit_ReadDataNoRemelt(int id, int &nx, int &MyYSlices, int &MyYOffset, double deltax, int HTtoCAratio,
                               double deltat, int nz, int LocalDomainSize, ViewI &CritTimeStep,
                               ViewF &UndercoolingChange, float XMin, float YMin, float ZMin, float *ZMinLayer,
                               float *ZMaxLayer, int LayerHeight, int NumberOfLayers, int *FinishTimeStep,
                               double FreezingRange, ViewI &LayerID, int *FirstValue, int *LastValue,
                               std::vector<double> RawData);
void calcMaxSolidificationEventsR(int id, int layernumber, int TempFilesInSeries, ViewI_H MaxSolidificationEvents_Host,
                                  int StartRange, int EndRange, std::vector<double> RawData, float XMin, float YMin,
                                  double deltax, float *ZMinLayer, int LayerHeight, int nx, int MyYSlices,
                                  int MyYOffset, int LocalActiveDomainSize);
void TempInit_ReadDataRemelt(int layernumber, int id, int nx, int MyYSlices, int nz, int LocalActiveDomainSize,
                             int LocalDomainSize, int MyYOffset, double &deltax, double deltat, double FreezingRange,
                             ViewF3D &LayerTimeTempHistory, ViewI &NumberOfSolidificationEvents,
                             ViewI &MaxSolidificationEvents, ViewI &MeltTimeStep, ViewI &CritTimeStep,
                             ViewF &UndercoolingChange, ViewF &UndercoolingCurrent, float XMin, float YMin,
                             float *ZMinLayer, int LayerHeight, int nzActive, int ZBound_Low, int *FinishTimeStep,
                             ViewI &LayerID, int *FirstValue, int *LastValue, std::vector<double> RawData,
                             ViewI &SolidificationEventCounter, int TempFilesInSeries);
void SubstrateInit_ConstrainedGrowth(int id, double FractSurfaceSitesActive, int MyYSlices, int nx, int ny,
                                     int MyYOffset, NList NeighborX, NList NeighborY, NList NeighborZ,
                                     ViewF GrainUnitVector, int NGrainOrientations, ViewI CellType, ViewI GrainID,
                                     ViewF DiagonalLength, ViewF DOCenter, ViewF CritDiagonalLength, double RNGSeed,
                                     int np, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, int BufSizeX,
                                     bool AtNorthBoundary, bool AtSouthBoundary);
void SubstrateInit_FromFile(std::string SubstrateFileName, int nz, int nx, int MyYSlices, int MyYOffset, int pid,
                            ViewI &GrainID, int nzActive, bool BaseplateThroughPowder);
void BaseplateInit_FromGrainSpacing(float SubstrateGrainSpacing, int nx, int ny, float *ZMinLayer, float *ZMaxLayer,
                                    int MyYSlices, int MyYOffset, int id, double deltax, ViewI GrainID, double RNGSeed,
                                    int &NextLayer_FirstEpitaxialGrainID, int nz, double BaseplateThroughPowder);
void PowderInit(int layernumber, int nx, int ny, int LayerHeight, float *ZMaxLayer, float ZMin, double deltax,
                int MyYSlices, int MyYOffset, int id, ViewI GrainID, double RNGSeed,
                int &NextLayer_FirstEpitaxialGrainID, double PowderDensity);
void CellTypeInit_Remelt(int nx, int MyYSlices, int LocalActiveDomainSize, ViewI CellType, ViewI CritTimeStep, int id,
                         int ZBound_Low);
void CellTypeInit_NoRemelt(int layernumber, int id, int np, int nx, int MyYSlices, int MyYOffset, int ZBound_Low,
                           int nz, int LocalActiveDomainSize, int LocalDomainSize, ViewI CellType, ViewI CritTimeStep,
                           NList NeighborX, NList NeighborY, NList NeighborZ, int NGrainOrientations,
                           ViewF GrainUnitVector, ViewF DiagonalLength, ViewI GrainID, ViewF CritDiagonalLength,
                           ViewF DOCenter, ViewI LayerID, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend,
                           int BufSizeX, bool AtNorthBoundary, bool AtSouthBoundary);
void NucleiInit(int layernumber, double RNGSeed, int MyYSlices, int MyYOffset, int nx, int ny, int nzActive,
                int ZBound_Low, int id, double NMax, double dTN, double dTsigma, double deltax, ViewI &NucleiLocation,
                ViewI_H &NucleationTimes_Host, ViewI &NucleiGrainID, ViewI CellType, ViewI CritTimeStep,
                ViewF UndercoolingChange, ViewI LayerID, int &PossibleNuclei_ThisRankThisLayer, int &Nuclei_WholeDomain,
                bool AtNorthBoundary, bool AtSouthBoundary, bool RemeltingYN, int &NucleationCounter,
                ViewI &MaxSolidificationEvents, ViewI NumberOfSolidificationEvents, ViewF3D LayerTimeTempHistory);
void placeNucleiData_NoRemelt(int Nuclei_ThisLayerSingle, ViewI_H NucleiX, ViewI_H NucleiY, ViewI_H NucleiZ,
                              int MyYOffset, int nx, int MyYSlices, bool AtNorthBoundary, bool AtSouthBoundary,
                              int ZBound_Low, ViewI_H CellType_Host, ViewI_H LayerID_Host, ViewI_H CritTimeStep_Host,
                              ViewF_H UndercoolingChange_Host, int layernumber,
                              std::vector<int> NucleiGrainID_WholeDomain_V,
                              std::vector<double> NucleiUndercooling_WholeDomain_V,
                              std::vector<int> &NucleiGrainID_MyRank_V, std::vector<int> &NucleiLocation_MyRank_V,
                              std::vector<int> &NucleationTimes_MyRank_V, int &PossibleNuclei_ThisRankThisLayer);
void placeNucleiData_Remelt(int NucleiMultiplier, int Nuclei_ThisLayerSingle, ViewI_H NucleiX, ViewI_H NucleiY,
                            ViewI_H NucleiZ, int MyYOffset, int nx, int MyYSlices, bool AtNorthBoundary,
                            bool AtSouthBoundary, int ZBound_Low, ViewI_H NumberOfSolidificationEvents_Host,
                            ViewF3D_H LayerTimeTempHistory_Host, std::vector<int> NucleiGrainID_WholeDomain_V,
                            std::vector<double> NucleiUndercooling_WholeDomain_V,
                            std::vector<int> &NucleiGrainID_MyRank_V, std::vector<int> &NucleiLocation_MyRank_V,
                            std::vector<int> &NucleationTimes_MyRank_V, int &PossibleNuclei_ThisRankThisLayer);
void ZeroResetViews(int LocalActiveDomainSize, int BufSizeX, int BufSizeZ, ViewF &DiagonalLength,
                    ViewF &CritDiagonalLength, ViewF &DOCenter, Buffer2D &BufferNorthSend, Buffer2D &BufferSouthSend,
                    Buffer2D &BufferNorthRecv, Buffer2D &BufferSouthRecv, ViewI &SteeringVector);

#endif
