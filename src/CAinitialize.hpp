// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_INIT_HPP
#define EXACA_INIT_HPP

#include "CAprint.hpp"
#include "CAtypes.hpp"
#include <Kokkos_Core.hpp>

#include <nlohmann/json.hpp>

#include <string>
#include <vector>

void InputReadFromFile(int id, std::string InputFile, std::string &SimulationType, double &deltax, double &NMax,
                       double &dTN, double &dTsigma, std::string &GrainOrientationFile, int &TempFilesInSeries,
                       std::vector<std::string> &temp_paths, double &HT_deltax, double &deltat, int &NumberOfLayers,
                       int &LayerHeight, std::string &MaterialFileName, std::string &SubstrateFileName,
                       float &SubstrateGrainSpacing, bool &UseSubstrateFile, double &G, double &R, int &nx, int &ny,
                       int &nz, double &FractSurfaceSitesActive, int &NSpotsX, int &NSpotsY, int &SpotOffset,
                       int &SpotRadius, double &RNGSeed, bool &BaseplateThroughPowder, double &PowderActiveFraction,
                       bool &LayerwiseTempRead, bool &PowderFirstLayer, Print &print);
void checkPowderOverflow(int nx, int ny, int LayerHeight, int NumberOfLayers, bool BaseplateThroughPowder,
                         double PowderDensity);
void NeighborListInit(NList &NeighborX, NList &NeighborY, NList &NeighborZ);
void FindXYZBounds(std::string SimulationType, int id, double &deltax, int &nx, int &ny, int &nz,
                   std::vector<std::string> &temp_paths, double &XMin, double &XMax, double &YMin, double &YMax,
                   double &ZMin, double &ZMax, int &LayerHeight, int NumberOfLayers, int TempFilesInSeries,
                   double *ZMinLayer, double *ZMaxLayer, int SpotRadius);
void DomainDecomposition(int id, int np, int &MyYSlices, int &MyYOffset, int &NeighborRank_North,
                         int &NeighborRank_South, int &nx, int &ny, int &nz, long int &LocalDomainSize,
                         bool &AtNorthBoundary, bool &AtSouthBoundary);
void ReadTemperatureData(int id, double &deltax, double HT_deltax, int &HTtoCAratio, int MyYSlices, int MyYOffset,
                         double YMin, std::vector<std::string> &temp_paths, int NumberOfLayers, int TempFilesInSeries,
                         int *FirstValue, int *LastValue, bool LayerwiseTempRead, int layernumber,
                         ViewD_H &RawTemperatureData);
int calcZBound_Low(std::string SimulationType, int LayerHeight, int layernumber, double *ZMinLayer, double ZMin,
                   double deltax);
int calcZBound_High(std::string SimulationType, int SpotRadius, int LayerHeight, int layernumber, double ZMin,
                    double deltax, int nz, double *ZMaxLayer);
int calcnzActive(int ZBound_Low, int ZBound_High, int id, int layernumber);
int calcLocalActiveDomainSize(int nx, int MyYSlices, int nzActive);
void TempInit_DirSolidification(double G, double R, int id, int &nx, int &MyYSlices, double deltax, double deltat,
                                int nz, int LocalDomainSize, ViewI &CritTimeStep, ViewF &UndercoolingChange,
                                ViewI &NumberOfSolidificationEvents, ViewI &SolidificationEventCounter,
                                ViewI &MeltTimeStep, ViewI MaxSolidificationEvents, ViewF3D &LayerTimeTempHistory);
int calcMaxSolidificationEventsSpot(int nx, int MyYSlices, int NumberOfSpots, int NSpotsX, int SpotRadius,
                                    int SpotOffset, int MyYOffset);
void TempInit_Spot(int layernumber, double G, double R, std::string, int id, int &nx, int &MyYSlices, int &MyYOffset,
                   double deltax, double deltat, int ZBound_Low, int nz, int LocalActiveDomainSize, int LocalDomainSize,
                   ViewI &CritTimeStep, ViewF &UndercoolingChange, ViewF &UndercoolingCurrent, int LayerHeight,
                   double FreezingRange, int NSpotsX, int NSpotsY, int SpotRadius, int SpotOffset,
                   ViewF3D &LayerTimeTempHistory, ViewI &NumberOfSolidificationEvents, ViewI &MeltTimeStep,
                   ViewI &MaxSolidificationEvents, ViewI &SolidificationEventCounter);
int getTempCoordX(int i, double XMin, double deltax, const ViewD_H RawTemperatureData);
int getTempCoordY(int i, double YMin, double deltax, const ViewD_H RawTemperatureData);
int getTempCoordZ(int i, double deltax, const ViewD_H RawTemperatureData, int LayerHeight, int LayerCounter,
                  double *ZMinLayer);
double getTempCoordTM(int i, const ViewD_H RawTemperatureData);
double getTempCoordTL(int i, const ViewD_H RawTemperatureData);
double getTempCoordCR(int i, const ViewD_H RawTemperatureData);
void calcMaxSolidificationEventsR(int id, int layernumber, int TempFilesInSeries, ViewI_H MaxSolidificationEvents_Host,
                                  int StartRange, int EndRange, ViewD_H RawTemperatureData, double XMin, double YMin,
                                  double deltax, double *ZMinLayer, int LayerHeight, int nx, int MyYSlices,
                                  int MyYOffset, int LocalActiveDomainSize);
void TempInit_ReadData(int layernumber, int id, int nx, int MyYSlices, int nz, int LocalActiveDomainSize,
                       int LocalDomainSize, int MyYOffset, double &deltax, double deltat, double FreezingRange,
                       ViewF3D &LayerTimeTempHistory, ViewI &NumberOfSolidificationEvents,
                       ViewI &MaxSolidificationEvents, ViewI &MeltTimeStep, ViewI &CritTimeStep,
                       ViewF &UndercoolingChange, ViewF &UndercoolingCurrent, double XMin, double YMin,
                       double *ZMinLayer, int LayerHeight, int nzActive, int ZBound_Low, int *FinishTimeStep,
                       int *FirstValue, int *LastValue, ViewD_H RawTemperatureData, ViewI &SolidificationEventCounter,
                       int TempFilesInSeries);
void ZeroResetViews(int LocalActiveDomainSize, ViewF &DiagonalLength, ViewF &CritDiagonalLength, ViewF &DOCenter,
                    ViewI &SteeringVector);

//*****************************************************************************/
// Initialize grain orientations and unit vectors
template <typename ViewTypeFloat>
void OrientationInit(int, int &NGrainOrientations, ViewTypeFloat &GrainOrientationData,
                     std::string GrainOrientationFile, int ValsPerLine = 9) {

    // Read file of grain orientations
    std::ifstream O;
    O.open(GrainOrientationFile);

    // Line 1 is the number of orientation values to read (if not specified already)
    std::string ValueRead;
    getline(O, ValueRead);
    NGrainOrientations = getInputInt(ValueRead);

    // Temporary host view for storing grain orientations read from file
    using view_type_host = typename ViewTypeFloat::HostMirror;
    view_type_host GrainOrientationData_Host(Kokkos::ViewAllocateWithoutInitializing("GrainOrientationData_H"),
                                             ValsPerLine * NGrainOrientations);
    // Populate data structure for grain orientation data
    for (int i = 0; i < NGrainOrientations; i++) {
        std::vector<std::string> ParsedLine(ValsPerLine);
        std::string ReadLine;
        if (!getline(O, ReadLine))
            break;
        splitString(ReadLine, ParsedLine, ValsPerLine);
        // Place the 3 grain orientation angles or 9 rotation matrix components into the orientation data view
        for (int Comp = 0; Comp < ValsPerLine; Comp++) {
            GrainOrientationData_Host(ValsPerLine * i + Comp) = getInputFloat(ParsedLine[Comp]);
        }
    }
    O.close();

    // Resize device view and orientation data to device
    Kokkos::realloc(GrainOrientationData, ValsPerLine * NGrainOrientations);
    using memory_space = typename ViewTypeFloat::memory_space;
    GrainOrientationData = Kokkos::create_mirror_view_and_copy(memory_space(), GrainOrientationData_Host);
}

#endif
