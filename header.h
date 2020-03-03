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
void RunProgram_Reduced(int pid, int np, int ierr, int DecompositionStrategy, double deltax, double AConst, double BConst, double CConst, double DConst, double NMax, double dTN, double dTsigma, string BaseFileName, string GrainOrientationFile, string TemperatureDataType);
void RunProgram_Bulky(int pid, int np, int ierr, int DecompositionStrategy, double deltax, double AConst, double BConst, double CConst, double NMax, double dTN, double dTsigma, string BaseFileName, string GrainOrientationFile);
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
void CritDiagLengthCalc(double xp, double yp, double zp, int MyOrientation, int NewCellX, int NewCellY, int NewCellZ, int CellLocation, double cx, double cy, double cz, int NeighborX[26], int NeighborY[26], int NeighborZ[26], float* GrainUnitVector, ViewI::HostMirror TriangleIndex, ViewF::HostMirror CritDiagonalLength);
void InitialDecomposition(int &DecompositionStrategy, int nx, int ny, int &ProcessorsInXDirection, int &ProcessorsInYDirection, int id, int np, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset,int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn, int &MyLeftOut, int &MyRightIn, int &MyRightOut);
void XYLimitCalc(int &LLX, int &LLY, int &ULX, int &ULY, int MyXSlices, int MyYSlices, int MyLeft, int MyRight, int MyIn, int MyOut);

// Contained in "CAInitialize.cpp"
void MasterInputRead(int &DecompositionStrategy, double &deltax, double &AConst, double &BConst, double &CConst, double &DConst, double &NMax, double &dTN, double &dTsigma, string &BaseFileName, string &GrainOrientationFile, string &TemperatureDataType);
void CInputRead(double &G, double &R, int &nx, int &ny, int &nz, double deltax, double &deltat, double &FractSurfaceSitesActive);
void RInputRead(string &tempfile, double &HT_deltax, double &deltat, int &NumberOfLayers, int &LayerOffset, string &BaseFileName, string &TemperatureDataSource, string &SubstrateFileName, bool &LayerwiseTemperature, int &TempFilesInSeries);
void ParallelMeshInit(double &G, double &R, int &DecompositionStrategy, int (&NeighborX)[26], int (&NeighborY)[26], int (&NeighborZ)[26], int (&ItList)[9][26], string TemperatureDataType, int ierr, int id, int np, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset,int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn, int &MyLeftOut, int &MyRightIn, int &MyRightOut, double &deltax, double HT_deltax, double &deltat, int &nx, int &ny, int &nz, int &ProcessorsInXDirection, int &ProcessorsInYDirection, string tempfile, float &XMin, float &XMax, float &YMin, float &YMax, float &ZMin, float &ZMax, string TemperatureDataSource, bool LayerwiseTemperature);
void TempInit(int layernumber, int TempFilesInSeries, double G, double R, int DecompositionStrategy, int (&NeighborX)[26], int (&NeighborY)[26], int (&NeighborZ)[26], int (&ItList)[9][26], string TemperatureDataType, int ierr, int id, int np, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset,int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn, int &MyLeftOut, int &MyRightIn, int &MyRightOut, double deltax, double HT_deltax, double deltat, int &nx, int &ny, int &nz, int &ProcessorsInXDirection, int &ProcessorsInYDirection, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange, ViewF::HostMirror UndercoolingCurrent, string tempfile, float XMin, float XMax, float YMin, float YMax, float ZMin, float ZMax, bool* Melted, string TemperatureDataSource, bool LayerwiseTemperature);
void OrientationInit(int id, int NGrainOrientations, int* GrainOrientation, float* GrainUnitVector, string GrainOrientationFile);
void GrainInit(int layernumber, int LayerHeight, string TemperatureDataType, string SubstrateFileName, double FractSurfaceSitesActive, int NGrainOrientations, int DecompositionStrategy, int ProcessorsInXDirection, int ProcessorsInYDirection, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int id, int np, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, int ItList[9][26], int NeighborX[26], int NeighborY[26], int NeighborZ[26], int* GrainOrientation, float* GrainUnitVector, ViewF::HostMirror DiagonalLength, ViewI::HostMirror CellType, ViewI::HostMirror TriangleIndex, ViewI::HostMirror GrainID, ViewF::HostMirror CritDiagonalLength, ViewF::HostMirror DOCenter, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange, bool* Melted, double deltax, double NMax, int &NextLayer_FirstNucleatedGrainID, int &PossibleNuclei_ThisRank);
void NucleiInit(int DecompositionStrategy, int MyXSlices, int MyYSlices, int nz, int id, double dTN, double dTsigma, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, int &PossibleNuclei_ThisRank, ViewI::HostMirror NucleiLocation, ViewI::HostMirror NucleationTimes, int* GrainOrientation, ViewI::HostMirror CellType, ViewI::HostMirror GrainID, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange);
void LayerSetup(string SubstrateFileName, int layernumber, int LayerHeight, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nz, int LocalDomainSize, int id, int np, ViewF::HostMirror DiagonalLength, ViewI::HostMirror CellType, ViewI::HostMirror TriangleIndex, ViewI::HostMirror GrainID, int* GrainID_Stored, ViewF::HostMirror CritDiagonalLength, ViewF::HostMirror DOCenter, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange, ViewF::HostMirror UndercoolingCurrent, bool* Melted, bool* Melted_Stored, bool LayerwiseTemeperature, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferC, Buffer2D BufferD, Buffer2D BufferE, Buffer2D BufferF, Buffer2D BufferG, Buffer2D BufferH, Buffer2D BufferAR, Buffer2D BufferBR, Buffer2D BufferCR, Buffer2D BufferDR, Buffer2D BufferER, Buffer2D BufferFR, Buffer2D BufferGR, Buffer2D BufferHR, int BufSizeX, int BufSizeY, int BufSizeZ);

// Contained in "CAupdate.cpp"
void Nucleation(int id, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, const int nz, int cycle, int &nn, ViewI CritTimeStep, ViewI CellType, ViewF UndercoolingCurrent, ViewF UndercoolingChange, ViewI NucleiLocations, ViewI NucleationTimes, ViewI GrainID, int* GrainOrientation, ViewF DOCenter, int NeighborX[26], int NeighborY[26], int NeighborZ[26], float* GrainUnitVector, ViewI TriangleIndex, ViewF CritDiagonalLength, ViewF DiagonalLength, int NGrainOrientations, int PossibleNuclei_ThisRank);
void CellCapture(int id, int np, int cycle, int DecompositionStrategy, int LocalDomainSize, int MyXSlices, int MyYSlices, const int nz, double AConst, double BConst, double CConst, double DConst, int MyXOffset, int MyYOffset, int ItList[9][26], int NeighborX[26], int NeighborY[26], int NeighborZ[26], ViewI CritTimeStep, ViewF UndercoolingCurrent, ViewF UndercoolingChange, float* GrainUnitVector, ViewI TriangleIndex, ViewF CritDiagonalLength, ViewF DiagonalLength, int* GrainOrientation, ViewI CellType, ViewF DOCenter, ViewI GrainID, int NGrainOrientations, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferC, Buffer2D BufferD, Buffer2D BufferE, Buffer2D BufferF, Buffer2D BufferG, Buffer2D BufferH, int BufSizeX, int BufSizeY, ViewI Locks);
void IntermediateOutputAndCheck(int pid, int &cycle, int LocalDomainSize, int nn, int &XSwitch, ViewI CellType, ViewI CritTimeStep, string TemperatureDataType);

// Contained in "CAghostnodes.cpp"
void GhostNodes1D(int cycle, int id, int MyLeft, int MyRight, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nz, int NeighborX[26], int NeighborY[26], int NeighborZ[26], ViewI::HostMirror CellType, ViewF::HostMirror DOCenter, ViewI::HostMirror GrainID, float* GrainUnitVector, ViewI::HostMirror TriangleIndex, int* GrainOrientation, ViewF::HostMirror DiagonalLength, ViewF::HostMirror CritDiagonalLength, int NGrainOrientations);
void GhostNodes2D(int cycle, int id, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nz, int NeighborX[26], int NeighborY[26], int NeighborZ[26], ViewI::HostMirror CellType, ViewF::HostMirror DOCenter, ViewI::HostMirror GrainID, float* GrainUnitVector, ViewI::HostMirror TriangleIndex, int* GrainOrientation, ViewF::HostMirror DiagonalLength, ViewF::HostMirror CritDiagonalLength, int NGrainOrientations);
void GhostNodes1D_GPU(int cycle, int id, int MyLeft, int MyRight, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nz, int NeighborX[26], int NeighborY[26], int NeighborZ[26], ViewI CellType, ViewF DOCenter, ViewI GrainID, float* GrainUnitVector, ViewI TriangleIndex, int* GrainOrientation, ViewF DiagonalLength, ViewF CritDiagonalLength, int NGrainOrientations, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferAR, Buffer2D BufferBR, int BufSizeX, int BufSizeY, int BufSizeZ, ViewI Locks);
void GhostNodes2D_GPU(int cycle, int id, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nz, int NeighborX[26], int NeighborY[26], int NeighborZ[26], ViewI CellType, ViewF DOCenter, ViewI GrainID, float* GrainUnitVector, ViewI TriangleIndex, int* GrainOrientation, ViewF DiagonalLength, ViewF CritDiagonalLength, int NGrainOrientations, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferC, Buffer2D BufferD, Buffer2D BufferE, Buffer2D BufferF, Buffer2D BufferG, Buffer2D BufferH, Buffer2D BufferAR, Buffer2D BufferBR, Buffer2D BufferCR, Buffer2D BufferDR, Buffer2D BufferER, Buffer2D BufferFR, Buffer2D BufferGR, Buffer2D BufferHR, int BufSizeX, int BufSizeY, int BufSizeZ, ViewI Locks);

// Contained in "CAPrint.cpp"
void PrintValues(int pid, int np,int nx, int ny,int nz, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI::HostMirror GrainID, int* GrainOrientation,float* GrainUnitVector,string BaseFileName, int DecompositionStrategy, int NGrainOrientations, bool* Melted);
void PrintTempValues(int pid, int np, int nx, int ny, int nz, int MyXSlices,int MyYSlices,int ProcessorsInXDirection,int ProcessorsInYDirection, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange, int DecompositionStrategy);
void PrintValuesMultilayer(int NumberOfLayers, int LayerHeight, int id, int np, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI::HostMirror GrainID, int* GrainID_Stored, int* GrainOrientation, float* GrainUnitVector, string BaseFileName, int DecompositionStrategy, int NGrainOrientations, bool* Melted,bool* MeltedStored);
void PrintCT(int pid, int np,int nx, int ny,int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI::HostMirror CellType, string BaseFileName, int DecompositionStrategy);
