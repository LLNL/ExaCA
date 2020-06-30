#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <array>
#include <cstdarg>
#include <string>
#include <cmath>
#include <math.h>
#include <vector>
#include <random>
#include <Kokkos_Core.hpp>
#include "mpi.h"

using namespace std;

enum TypeNames { Solid = 0, LiqSol = 1, Liquid = 2,
                 Delayed = 3, Active = 4, Wall = 5,
                 Ghost1 = 6, Ghost2 = 7, Ghost3 = 8 };

typedef Kokkos::View<float*> ViewF;
typedef Kokkos::View<int*> ViewI;
typedef Kokkos::View<int*, Kokkos::MemoryTraits<Kokkos::Atomic>> View_a;
typedef Kokkos::View<float**,Kokkos::CudaSpace> Buffer2D;
typedef Kokkos::View<float*> TestView;


// Called from "main_main.cpp"
void RunProgram_Reduced(int id, int np, int ierr, string InputFile);
int main ( int argc, char *argv[] );

// Contained in "CAfunctions.cpp"
double CrossP1(double TestVec1[3], double TestVec2[3]);
double CrossP2(double TestVec1[3], double TestVec2[3]);
double CrossP3(double TestVec1[3], double TestVec2[3]);
int FindItBounds(int RankX, int RankY, int MyXSlices, int MyYSlices);
int MaxIndex(double TestVec3[6]);
int XMPSlicesCalc(int p, int nx, int ProcessorsInXDirection, int ProcessorsInYDirection, int np, int DecompositionStrategy);
int XOffsetCalc(int p, int nx, int ProcessorsInXDirection, int ProcessorsInYDirection, int np, int DecompositionStrategy);
int YMPSlicesCalc(int p, int ny, int ProcessorsInYDirection, int np, int DecompositionStrategy);
int YOffsetCalc(int p, int ny, int ProcessorsInYDirection, int np, int DecompositionStrategy);
double MaxVal(double TestVec3[6],int NVals);
void NewGrain(int MyXSlices, int MyYSlices, int nz, int RankX, int RankY, int RankZ, int MyGrainID, int GlobalX, int GlobalY, ViewI::HostMirror CellType, ViewI::HostMirror GrainID, ViewF::HostMirror DiagonalLength, ViewF::HostMirror DOCenter);
void CritDiagLengthCalc(double xp, double yp, double zp, int MyOrientation, int NewCellX, int NewCellY, int NewCellZ, long int CellLocation, double cx, double cy, double cz, int NeighborX[26], int NeighborY[26], int NeighborZ[26], float* GrainUnitVector, ViewF::HostMirror CritDiagonalLength);
void InitialDecomposition(int &DecompositionStrategy, int nx, int ny, int &ProcessorsInXDirection, int &ProcessorsInYDirection, int id, int np, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset,int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn, int &MyLeftOut, int &MyRightIn, int &MyRightOut);
void XYLimitCalc(int &LLX, int &LLY, int &ULX, int &ULY, int MyXSlices, int MyYSlices, int MyLeft, int MyRight, int MyIn, int MyOut);

// Contained in "CAInitialize.cpp"
void InputReadFromFile(int id, string InputFile, string &SimulationType, int &DecompositionStrategy, double &AConst, double &BConst, double &CConst, double &DConst, double &FreezingRange, double &deltax, double &NMax, double &dTN, double &dTsigma, string &OutputFile, string &GrainOrientationFile, string &tempfile, int &TempFilesInSeries, string& BurstBuffer, string& ExtraWalls, double &HT_deltax, string &TemperatureDataSource, double &deltat, int &NumberOfLayers, int &LayerHeight, string &SubstrateFileName, double &G, double &R, int &nx, int &ny, int &nz, double &FractSurfaceSitesActive);
void ParallelMeshInit(int DecompositionStrategy, int (&NeighborX)[26], int (&NeighborY)[26], int (&NeighborZ)[26], int (&ItList)[9][26], string SimulationType, int ierr, int id, int np, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset,int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn, int &MyLeftOut, int &MyRightIn, int &MyRightOut, double &deltax, double HT_deltax, double &deltat, int &nx, int &ny, int &nz, int &ProcessorsInXDirection, int &ProcessorsInYDirection, string tempfile, float &XMin, float &XMax, float &YMin, float &YMax, float &ZMin, float &ZMax, string TemperatureDataSource, int LayerHeight, int NumberOfLayers, int TempFilesInSeries, float* ZMinLayer, float* ZMaxLayer);
void TempInit(int layernumber, int TempFilesInSeries, double G, double R, int DecompositionStrategy, int (&NeighborX)[26], int (&NeighborY)[26], int (&NeighborZ)[26], int (&ItList)[9][26], string SimulationType, int ierr, int id, int np, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset,int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn, int &MyLeftOut, int &MyRightIn, int &MyRightOut, double deltax, double HT_deltax, double deltat, int &nx, int &ny, int &nz, int &ProcessorsInXDirection, int &ProcessorsInYDirection, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange, ViewF::HostMirror UndercoolingCurrent, string tempfile, float XMin, float XMax, float YMin, float YMax, float ZMin, float ZMax, bool* Melted, string TemperatureDataSource, float* ZMinLayer, float* ZMaxLayer, int LayerHeight, int NumberOfLayers, int &nzActive, int &ZBound_Low, int &ZBound_High, int* FinishTimeStep, double FreezingRange, ViewI::HostMirror LayerID);
void OrientationInit(int id, int NGrainOrientations, int* GrainOrientation, float* GrainUnitVector, string GrainOrientationFile);
void GrainInit(int layernumber, int LayerHeight, string SimulationType, string SubstrateFileName, double FractSurfaceSitesActive, int NGrainOrientations, int DecompositionStrategy, int ProcessorsInXDirection, int ProcessorsInYDirection, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int id, int np, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, int ItList[9][26], int NeighborX[26], int NeighborY[26], int NeighborZ[26], int* GrainOrientation, float* GrainUnitVector, ViewF::HostMirror DiagonalLength, ViewI::HostMirror CellType, ViewI::HostMirror GrainID, ViewF::HostMirror CritDiagonalLength, ViewF::HostMirror DOCenter, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange, bool* Melted, double deltax, double NMax, int &NextLayer_FirstNucleatedGrainID, int &PossibleNuclei_ThisRank, int ZBound_High, int NumberOfLayers, int TempFilesInSeries, double HT_deltax, double deltat, double XMin, double XMax, double YMin, double YMax, double ZMin, double ZMax, string tempfile, string TemperatureDataSource, int ZBound_Low, string ExtraWalls);
void NucleiInit(int DecompositionStrategy, int MyXSlices, int MyYSlices, int nz, int id, double dTN, double dTsigma, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, int &PossibleNuclei_ThisRank, ViewI::HostMirror NucleiLocation, ViewI::HostMirror NucleationTimes, int* GrainOrientation, ViewI::HostMirror CellType, ViewI::HostMirror GrainID, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange);
void DomainShiftAndResize(int id, int MyXSlices, int MyYSlices, int &ZShift, int &ZBound_Low, int &ZBound_High, int &nzActive, int LocalDomainSize, int &LocalActiveDomainSize, int &BufSizeZ, int LayerHeight, ViewI CellType, int layernumber, ViewI LayerID);
void LayerSetup(string SubstrateFileName, int layernumber, int LayerHeight, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nz, int LocalDomainSize, int LocalActiveDomainSize, int* GrainOrientation, int NGrainOrientations, float* GrainUnitVector, int NeighborX[26], int NeighborY[26], int NeighborZ[26], int id, int np, ViewF DiagonalLength, ViewI CellType, ViewI GrainID, ViewF CritDiagonalLength, ViewF DOCenter, ViewI CritTimeStep, ViewF UndercoolingChange, ViewF UndercoolingCurrent, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferC, Buffer2D BufferD, Buffer2D BufferE, Buffer2D BufferF, Buffer2D BufferG, Buffer2D BufferH, Buffer2D BufferAR, Buffer2D BufferBR, Buffer2D BufferCR, Buffer2D BufferDR, Buffer2D BufferER, Buffer2D BufferFR, Buffer2D BufferGR, Buffer2D BufferHR, int BufSizeX, int BufSizeY, int BufSizeZ, int TempFilesInSeries, int &ZBound_Low, int &ZBound_High, float ZMin, double deltax, int nzActive, ViewI Locks);
//void EraseTop(int NextLayerNumber, int TempFilesInSeries, int DecompositionStrategy, string SimulationType, int id, int np, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset, double deltax, double HT_deltax, double deltat, int &nx, int &ny, int &nz, int &ProcessorsInXDirection, int &ProcessorsInYDirection, ViewI::HostMirror GrainID, string tempfile, float XMin, float XMax, float YMin, float YMax, float ZMin, float ZMax, bool* Melted, string TemperatureDataSource, int LayersSimulatedAtOnce, int LayerHeight, int NumberOfLayers, ViewI::HostMirror CellType, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange);

// Contained in "CAupdate.cpp"
void Nucleation(int id, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, const int nz, int cycle, int &nn, ViewI CritTimeStep, ViewI CellType, ViewF UndercoolingCurrent, ViewF UndercoolingChange, ViewI NucleiLocations, ViewI NucleationTimes, ViewI GrainID, int* GrainOrientation, ViewF DOCenter, int NeighborX[26], int NeighborY[26], int NeighborZ[26], float* GrainUnitVector, ViewF CritDiagonalLength, ViewF DiagonalLength, int NGrainOrientations, int PossibleNuclei_ThisRank, ViewI Locks, int ZBound_Low, int layernumber, ViewI LayerID);
void CellCapture(int id, int np, int cycle, int DecompositionStrategy, int LocalActiveDomainSize, int MyXSlices, int MyYSlices, const int nz, double AConst, double BConst, double CConst, double DConst, int MyXOffset, int MyYOffset, int ItList[9][26], int NeighborX[26], int NeighborY[26], int NeighborZ[26], ViewI CritTimeStep, ViewF UndercoolingCurrent, ViewF UndercoolingChange, float* GrainUnitVector, ViewF CritDiagonalLength, ViewF DiagonalLength, int* GrainOrientation, ViewI CellType, ViewF DOCenter, ViewI GrainID, int NGrainOrientations, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferC, Buffer2D BufferD, Buffer2D BufferE, Buffer2D BufferF, Buffer2D BufferG, Buffer2D BufferH, int BufSizeX, int BufSizeY, ViewI Locks, int ZBound_Low, int nzActive, int layernumber, ViewI LayerID);
void IntermediateOutputAndCheck(int pid, int &cycle, int MyXSlices, int MyYSlices, int LocalDomainSize, int LocalActiveDomainSize, int nn, int &XSwitch, ViewI CellType, ViewI CritTimeStep, string SimulationType, int* FinishTimeStep, int layernumber, int NumberOfLayers, int ZBound_Low, ViewI Locks, ViewI LayerID);

// Contained in "CAghostnodes.cpp"
void GhostNodesInit_GPU(int DecompositionStrategy,int MyXSlices, int MyYSlices, ViewI GrainID, ViewI CellType, ViewF DOCenter, ViewF DiagonalLength, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferC, Buffer2D BufferD, Buffer2D BufferE, Buffer2D BufferF, Buffer2D BufferG, Buffer2D BufferH, int BufSizeX, int BufSizeY, int LocalActiveDomainSize, int ZBound_Low);
void GhostNodes1D_GPU(int cycle, int id, int MyLeft, int MyRight, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nz, int NeighborX[26], int NeighborY[26], int NeighborZ[26], ViewI CellType, ViewF DOCenter, ViewI GrainID, float* GrainUnitVector, int* GrainOrientation, ViewF DiagonalLength, ViewF CritDiagonalLength, int NGrainOrientations, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferAR, Buffer2D BufferBR, int BufSizeX, int BufSizeY, int BufSizeZ, ViewI Locks, int ZBound_Low);
void GhostNodes2D_GPU(int cycle, int id, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nz, int NeighborX[26], int NeighborY[26], int NeighborZ[26], ViewI CellType, ViewF DOCenter, ViewI GrainID, float* GrainUnitVector, int* GrainOrientation, ViewF DiagonalLength, ViewF CritDiagonalLength, int NGrainOrientations, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferC, Buffer2D BufferD, Buffer2D BufferE, Buffer2D BufferF, Buffer2D BufferG, Buffer2D BufferH, Buffer2D BufferAR, Buffer2D BufferBR, Buffer2D BufferCR, Buffer2D BufferDR, Buffer2D BufferER, Buffer2D BufferFR, Buffer2D BufferGR, Buffer2D BufferHR, int BufSizeX, int BufSizeY, int BufSizeZ, ViewI Locks, int ZBound_Low);

// Contained in "CAPrint.cpp"
void PrintValues(int pid, int np,int nx, int ny,int nz, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI::HostMirror GrainID, int* GrainOrientation,float* GrainUnitVector,string OutputFile, int DecompositionStrategy, int NGrainOrientations, bool* Melted);
void PrintTempValues(int pid, int np, int nx, int ny, int nz, int MyXSlices,int MyYSlices,int ProcessorsInXDirection,int ProcessorsInYDirection, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange, int DecompositionStrategy);
void PrintCT(int pid, int np,int nx, int ny,int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI::HostMirror CellType, string OutputFile, int DecompositionStrategy);
