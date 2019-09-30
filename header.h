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
void InitialDecomposition(int DecompositionStrategy, int nx, int ny, int &ProcessorsInXDirection, int &ProcessorsInYDirection, int id, int np, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset,int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn, int &MyLeftOut, int &MyRightIn, int &MyRightOut);
void XYLimitCalc(int &LLX, int &LLY, int &ULX, int &ULY, int MyXSlices, int MyYSlices, int MyLeft, int MyRight, int MyIn, int MyOut);

// Contained in "CAInitialize.cpp"
void MasterInputRead(int &DecompositionStrategy, double &deltax, double &AConst, double &BConst, double &CConst, double &NMax, double &dTN, double &dTsigma, string &BaseFileName, string &GrainOrientationFile, string &TemperatureDataType, double &InitialGrainWidth);
void AInputRead(string &AFile1, string &AFile2, string &AFile3, int &NRatio, int &NumberOfLayers, int &LayerHeight);
void CInputRead(double &G, double &R, int &nx, int &ny, int &nz, int &NRatio, int &NumberOfLayers);
void UInputRead(double &G, double &R, double &DomUndercooling, int &nx, int &ny, int &nz, int &NumberOfLayers);
void RInputRead(string &tempfile, double &HT_deltax, double &deltat, int &NumberOfLayers, int &LayerHeight, string &BaseFileName);
void ParallelMeshInit(double &G, double &R, int DecompositionStrategy, int (&NeighborX)[26], int (&NeighborY)[26], int (&NeighborZ)[26], int (&ItList)[9][26], string TemperatureDataType, int ierr, int id, int np, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset,int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn, int &MyLeftOut, int &MyRightIn, int &MyRightOut, double &deltax, double HT_deltax, double &deltat, int &nx, int &ny, int &nz, int &ProcessorsInXDirection, int &ProcessorsInYDirection, string AFile3, string tempfile, float &XMin, float &XMax, float &YMin, float &YMax, float &ZMin, float &ZMax, double DomUndercooling, double NRatio);
void TempInit(double G, double R, int DecompositionStrategy, int (&NeighborX)[26], int (&NeighborY)[26], int (&NeighborZ)[26], int (&ItList)[9][26], string TemperatureDataType, int ierr, int id, int np, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset,int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn, int &MyLeftOut, int &MyRightIn, int &MyRightOut, double deltax, double HT_deltax, double deltat, int &nx, int &ny, int &nz, int &ProcessorsInXDirection, int &ProcessorsInYDirection, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange, ViewF::HostMirror UndercoolingCurrent, double DomUndercooling, string tempfile, string FileA1, string FileA2, float &XMin, float &XMax, float &YMin, float &YMax, float &ZMin, float &ZMax);
void OrientationInit(int pid, int NGrainOrientations, int* GrainOrientation, float* GrainUnitVector, string GrainOrientationFile);
void ConstrainedGrains(string TemperatureDataType, int InitialGrainWidth, int NGrainOrientations, int DecompositionStrategy, int ProcessorsInXDirection, int ProcessorsInYDirection, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int id, int np, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, int ItList[9][26], int NeighborX[26], int NeighborY[26], int NeighborZ[26], vector <int> &NucLocI, vector <int> &NucLocJ, vector <int> &NucLocK, vector <int> &NucleationTimes, vector <float> &NucleationUndercooling, int* GrainOrientation, float* GrainUnitVector, ViewF::HostMirror DiagonalLength, ViewI::HostMirror CellType, ViewI::HostMirror TriangleIndex, ViewI::HostMirror GrainID, ViewF::HostMirror CritDiagonalLength, ViewF::HostMirror DOCenter, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange, double deltax, double NMax, double dTN, double dTsigma, int &NextLayer_FirstSubstrateGrainID, int &NextLayer_FirstNucleatedGrainID, int &ACount, int &BCount, int &CCount, int &DCount, int &ECount, int &FCount, int &GCount, int &HCount);
void ConstrainedGrainsUpdate(int DecompositionStrategy, int MyXSlices, int MyYSlices, int id, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, vector <int> &NucLocI, vector <int> &NucLocJ, vector <int> &NucLocK, vector <int> &NucleationTimes, vector <float> &NucleationUndercooling, int* GrainOrientation, ViewI::HostMirror CellType, ViewI::HostMirror GrainID);
void LayerSetup(int layernumber, int LayerHeight, int InitialGrainWidth, int NGrainOrientations, int DecompositionStrategy, int ProcessorsInXDirection, int ProcessorsInYDirection, int nx, int ny, int nz, string TemperatureDataType, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int pid, int np, int ierr, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyLeftOut, int MyRightIn, int MyRightOut, int ItList[9][26], int NeighborX[26], int NeighborY[26], int NeighborZ[26], vector <int> &NucLocI, vector <int> &NucLocJ, vector <int> &NucLocK, vector <int> &NucleationTimes, vector <float> &NucleationUndercooling, int* GrainOrientation, float* GrainUnitVector, ViewF::HostMirror DiagonalLength, ViewI::HostMirror CellType, ViewI::HostMirror TriangleIndex, ViewI::HostMirror GrainID, int* GrainID_Stored, ViewF::HostMirror CritDiagonalLength, ViewF::HostMirror DOCenter, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange, ViewF::HostMirror UndercoolingCurrent, string tempfile, double deltax, double deltat, double NMax, double dTN, double dTsigma, string FileA1, string FileA2, int XMin, int XMax, int YMin, int YMax, int ZMin, int ZMax, int &NextLayer_FirstSubstrateGrainID, int &NextLayer_FirstNucleatedGrainID, int &ACount, int &BCount, int &CCount, int &DCount, int &ECount, int &FCount, int &GCount, int &HCount);

// Contained in "CAupdate.cpp"
void TemperatureUpdate(int pid, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nz, int cycle, int &nn, double AConst, double BConst, double CConst, ViewI CritTimeStep, ViewI CellType, ViewF UndercoolingCurrent, ViewF UndercoolingChange, vector <int> NucLocI, vector <int> NucLocJ, vector <int> NucLocK, vector <int> NucleationTimes, vector <float> NucleationUndercooling, ViewI GrainID, int* GrainOrientation, ViewF DOCenter, int NeighborX[26], int NeighborY[26], int NeighborZ[26], float* GrainUnitVector, ViewI TriangleIndex, ViewF CritDiagonalLength, ViewF DiagonalLength, int NGrainOrientations);
void CellCapture_ALT(int pid, int cycle, int DecompositionStrategy, int &ACount, int &BCount, int &CCount, int &DCount, int &ECount, int &FCount, int &GCount, int &HCount, int MyXSlices, int MyYSlices, int nz, int MyXOffset, int MyYOffset, int ItList[9][26], int NeighborX[26], int NeighborY[26], int NeighborZ[26], float* GrainUnitVector, ViewI::HostMirror  TriangleIndex, ViewF::HostMirror  CritDiagonalLength, ViewF::HostMirror  DiagonalLength, int* GrainOrientation, ViewI::HostMirror  CellType, ViewF::HostMirror  DOCenter, ViewI::HostMirror  GrainID, int NGrainOrientations);
void CellCapture(int pid, int cycle, int DecompositionStrategy, int &ACount, int &BCount, int &CCount, int &DCount, int &ECount, int &FCount, int &GCount, int &HCount, int MyXSlices, int MyYSlices, int nz, int MyXOffset, int MyYOffset, int ItList[9][26], int NeighborX[26], int NeighborY[26], int NeighborZ[26], float* GrainUnitVector, ViewI TriangleIndex, ViewF CritDiagonalLength, ViewF DiagonalLength, int* GrainOrientation, ViewI CellType, ViewF DOCenter, ViewI GrainID, int NGrainOrientations);
void IntermediateOutputAndCheck(int pid, int cycle, int MyXSlices, int MyYSlices, int nz, int nn, int &XSwitch, ViewI::HostMirror CellType);

// Contained in "CAghostnodes.cpp"
void GhostNodes1D(int cycle, int pid, int ACount, int BCount, int MyLeft, int MyRight, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nz, int NeighborX[26], int NeighborY[26], int NeighborZ[26], ViewI::HostMirror CellType, ViewF::HostMirror DOCenter, ViewI::HostMirror GrainID, float* GrainUnitVector, ViewI::HostMirror TriangleIndex, int* GrainOrientation, ViewF::HostMirror DiagonalLength, ViewF::HostMirror CritDiagonalLength, int NGrainOrientations);
void GhostNodes2D(int cycle, int pid, int ACount, int BCount, int CCount, int DCount, int ECount, int FCount, int GCount, int HCount, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nz, int NeighborX[26], int NeighborY[26], int NeighborZ[26], ViewI::HostMirror CellType, ViewF::HostMirror DOCenter, ViewI::HostMirror GrainID, float* GrainUnitVector, ViewI::HostMirror TriangleIndex, int* GrainOrientation, ViewF::HostMirror DiagonalLength, ViewF::HostMirror CritDiagonalLength, int NGrainOrientations);

// Contained in "CAmain.cpp"
void RunProgram();
int main ( int argc, char *argv[] );

// Contained in "CAPrint.cpp"
void PrintValues(int pid, int np,int nx, int ny,int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI::HostMirror GrainID, int* GrainOrientation,float* GrainUnitVector,string BaseFileName, int DecompositionStrategy, int NGrainOrientations);
void PrintTempValues(int pid, int np, int nx, int ny, int nz, int MyXSlices,int MyYSlices,int ProcessorsInXDirection,int ProcessorsInYDirection, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange, int DecompositionStrategy);
void PrintValuesMultilayer(int NumberOfLayers, int LayerHeight, int id, int np, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI::HostMirror GrainID, int* GrainID_Stored, int* GrainOrientation, float* GrainUnitVector, string BaseFileName, int DecompositionStrategy, int NGrainOrientations);
void PrintCT(int pid, int np,int nx, int ny,int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI::HostMirror CellType, string BaseFileName, int DecompositionStrategy);
