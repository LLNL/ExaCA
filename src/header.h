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
#include <vector>
#include <random>
#include <Kokkos_Core.hpp>
#include "mpi.h"

enum TypeNames { Wall = 0, Solid = 1, Active = 2,
                TemporaryUpdate = 3, TemporaryInit = 4, Liquid = 5,
                Ghost = 6};

// Use Kokkos::DefaultExecutionSpace
typedef Kokkos::View<float*> ViewF;
typedef Kokkos::View<int*> ViewI;
typedef Kokkos::View<int**> ViewI2D;
typedef Kokkos::View<int*, Kokkos::MemoryTraits<Kokkos::Atomic>> View_a;
typedef Kokkos::View<float**> Buffer2D;
typedef Kokkos::View<int***> Buffer3D; // Used in ghost node initialization of integer structures CellType and GrainID
typedef Kokkos::View<float*> TestView;

using exe_space = Kokkos::DefaultExecutionSpace::execution_space;
typedef typename exe_space::array_layout layout;
typedef Kokkos::View<float*, layout, Kokkos::HostSpace> ViewF_H;
typedef Kokkos::View<int*, layout, Kokkos::HostSpace> ViewI_H;
typedef Kokkos::View<int**, layout, Kokkos::HostSpace> ViewI2D_H;

// Called from "main_main.cpp"
void RunProgram_Reduced(int id, int np, std::string InputFile);
int main ( int argc, char *argv[] );

// Contained in "CAfunctions.cpp"
double CrossP1(double TestVec1[3], double TestVec2[3]);
double CrossP2(double TestVec1[3], double TestVec2[3]);
double CrossP3(double TestVec1[3], double TestVec2[3]);
int FindItBounds(int RankX, int RankY, int MyXSlices, int MyYSlices);
int MaxIndex(double TestVec3[6]);
int XMPSlicesCalc(int p, int nx, int ProcessorsInXDirection, int ProcessorsInYDirection, int DecompositionStrategy);
int XOffsetCalc(int p, int nx, int ProcessorsInXDirection, int ProcessorsInYDirection, int DecompositionStrategy);
int YMPSlicesCalc(int p, int ny, int ProcessorsInYDirection, int np, int DecompositionStrategy);
int YOffsetCalc(int p, int ny, int ProcessorsInYDirection, int np, int DecompositionStrategy);
double MaxVal(double TestVec3[6],int NVals);
void InitialDecomposition(int &DecompositionStrategy, int nx, int ny, int &ProcessorsInXDirection, int &ProcessorsInYDirection, int id, int np, int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn, int &MyLeftOut, int &MyRightIn, int &MyRightOut);
void XYLimitCalc(int &LLX, int &LLY, int &ULX, int &ULY, int MyXSlices, int MyYSlices, int MyLeft, int MyRight, int MyIn, int MyOut);

// Contained in "CAInitialize.cpp"
void InputReadFromFile(int id, std::string InputFile, std::string &SimulationType, int &DecompositionStrategy, double &AConst, double &BConst, double &CConst, double &DConst, double& FreezingRange, double &deltax, double &NMax, double &dTN, double &dTsigma, std::string &OutputFile, std::string &GrainOrientationFile, std::string &tempfile, int &TempFilesInSeries, bool &ExtraWalls, double &HT_deltax, bool &RemeltingYN, double &deltat, int &NumberOfLayers, int &LayerHeight, std::string &SubstrateFileName, float &SubstrateGrainSpacing, bool &UseSubstrateFile, double &G, double &R, int &nx, int &ny, int &nz, double &FractSurfaceSitesActive, std::string &PathToOutput, bool (&FilesToPrint)[6], bool &PrintFilesYN);
void ParallelMeshInit(int DecompositionStrategy, ViewI_H NeighborX, ViewI_H NeighborY, ViewI_H NeighborZ, ViewI2D_H ItList, std::string SimulationType, int id, int np, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset,int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn, int &MyLeftOut, int &MyRightIn, int &MyRightOut, double &deltax, double HT_deltax, int &nx, int &ny, int &nz, int &ProcessorsInXDirection, int &ProcessorsInYDirection, std::string tempfile, float &XMin, float &XMax, float &YMin, float &YMax, float &ZMin, float &ZMax, float FreezingRange, int &LayerHeight, int NumberOfLayers, int TempFilesInSeries, unsigned int &NumberOfTemperatureDataPoints, float* ZMinLayer, float* ZMaxLayer, int* FirstValue, int* LastValue, std::vector<float> &RawData);
void TempInit(int layernumber, double G, double R, std::string SimulationType, int id, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset, double deltax, double HT_deltax, double deltat, int &nx, int &ny, int &nz, ViewI_H CritTimeStep, ViewF_H UndercoolingChange, ViewF_H UndercoolingCurrent, float XMin, float YMin, float ZMin, bool* Melted, float* ZMinLayer, float* ZMaxLayer, int LayerHeight, int NumberOfLayers, int &nzActive, int &ZBound_Low, int &ZBound_High, int* FinishTimeStep, double FreezingRange, ViewI_H LayerID, int* FirstValue, int* LastValue, std::vector<float> RawData);
void OrientationInit(int id, int NGrainOrientations, ViewI_H GrainOrientation, ViewF_H GrainUnitVector, std::string GrainOrientationFile);
void SubstrateInit_ConstrainedGrowth(double FractSurfaceSitesActive, int MyXSlices, int MyYSlices, int nx, int ny, int nz, int MyXOffset, int MyYOffset, int pid, int np, ViewI_H CellType, ViewI_H GrainID);
void SubstrateInit_FromFile(std::string SubstrateFileName, int nz, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int pid, ViewI_H CritTimeStep, ViewI_H GrainID);
void SubstrateInit_FromGrainSpacing(float SubstrateGrainSpacing, int nx, int ny, int nz, int nzActive, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int LocalActiveDomainSize, int pid, int np, double deltax, ViewI_H GrainID, ViewI_H CritTimeStep);
void ActiveCellWallInit(int id, int MyXSlices, int MyYSlices, int nx, int ny, int nz, int MyXOffset, int MyYOffset, ViewI_H CellType_H, ViewI_H GrainID_H, ViewI_H CritTimeStep_H, ViewI2D_H ItList_H, ViewI_H NeighborX_H, ViewI_H NeighborY_H, ViewI_H NeighborZ_H, ViewF_H UndercoolingChange_H, bool ExtraWalls);
void GrainInit(int layernumber, int NGrainOrientations, int DecompositionStrategy, int nz, int LocalActiveDomainSize, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int id, int np, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, ViewI2D_H ItList, ViewI_H NeighborX, ViewI_H NeighborY, ViewI_H NeighborZ, ViewI_H GrainOrientation, ViewF_H GrainUnitVector, ViewF_H DiagonalLength, ViewI_H CellType, ViewI_H GrainID, ViewF_H CritDiagonalLength, ViewF_H DOCenter, ViewI_H CritTimeStep, double deltax, double NMax, int &NextLayer_FirstNucleatedGrainID, int &PossibleNuclei_ThisRank, int ZBound_High, int ZBound_Low);
void NucleiInit(int DecompositionStrategy, int MyXSlices, int MyYSlices, int nz, int id, double dTN, double dTsigma, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, ViewI_H NucleiLocation, ViewI_H NucleationTimes, ViewI_H CellType, ViewI_H GrainID, ViewI_H CritTimeStep, ViewF_H UndercoolingChange);
void DomainShiftAndResize(int id, int MyXSlices, int MyYSlices, int &ZShift, int &ZBound_Low, int &ZBound_High, int &nzActive, int LocalDomainSize, int &LocalActiveDomainSize, int &BufSizeZ, int LayerHeight, ViewI CellType, int layernumber, ViewI LayerID);
void LayerSetup(int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int LocalActiveDomainSize, ViewI GrainOrientation, int NGrainOrientations, ViewF GrainUnitVector, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewF DiagonalLength, ViewI CellType, ViewI GrainID, ViewF CritDiagonalLength, ViewF DOCenter, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferC, Buffer2D BufferD, Buffer2D BufferE, Buffer2D BufferF, Buffer2D BufferG, Buffer2D BufferH, Buffer2D BufferAR, Buffer2D BufferBR, Buffer2D BufferCR, Buffer2D BufferDR, Buffer2D BufferER, Buffer2D BufferFR, Buffer2D BufferGR, Buffer2D BufferHR, int BufSizeX, int BufSizeY, int BufSizeZ, int &ZBound_Low);

// Contained in "CAupdate.cpp"
void Nucleation(int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int cycle, int &nn, ViewI CellType, ViewI NucleiLocations, ViewI NucleationTimes, ViewI GrainID, ViewI GrainOrientation, ViewF DOCenter, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewF GrainUnitVector, ViewF CritDiagonalLength, ViewF DiagonalLength, int NGrainOrientations, int PossibleNuclei_ThisRank, int ZBound_Low, int layernumber, ViewI LayerID);
void CellCapture(int np, int cycle, int DecompositionStrategy, int LocalActiveDomainSize, int LocalDomainSize, int MyXSlices, int MyYSlices, double AConst, double BConst, double CConst, double DConst, int MyXOffset, int MyYOffset, ViewI2D ItList, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewI CritTimeStep, ViewF UndercoolingCurrent, ViewF UndercoolingChange, ViewF GrainUnitVector, ViewF CritDiagonalLength, ViewF DiagonalLength, ViewI GrainOrientation, ViewI CellType, ViewF DOCenter, ViewI GrainID, int NGrainOrientations, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferC, Buffer2D BufferD, Buffer2D BufferE, Buffer2D BufferF, Buffer2D BufferG, Buffer2D BufferH, int BufSizeX, int BufSizeY, int ZBound_Low, int nzActive, int nz, int layernumber, ViewI LayerID, ViewI SteeringVector, ViewI numSteer_G, ViewI_H numSteer_H);
void IntermediateOutputAndCheck(int pid, int &cycle, int MyXSlices, int MyYSlices, int LocalDomainSize, int LocalActiveDomainSize, int nn, int &XSwitch, ViewI CellType, ViewI CritTimeStep, std::string SimulationType, int* FinishTimeStep, int layernumber, int NumberOfLayers, int ZBound_Low, ViewI LayerID);

// Contained in "CAghostnodes.cpp"
void GhostNodesInit_GPU(int id, int np, int DecompositionStrategy, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int ZBound_Low, int nzActive, int LocalActiveDomainSize, int NGrainOrientations, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewF GrainUnitVector, ViewI GrainOrientation, ViewI GrainID, ViewI CellType, ViewF DOCenter, ViewF DiagonalLength, ViewF CritDiagonalLength);
void GhostNodes1D_GPU(int cycle, int id, int MyLeft, int MyRight, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewI CellType, ViewF DOCenter, ViewI GrainID, ViewF GrainUnitVector, ViewI GrainOrientation, ViewF DiagonalLength, ViewF CritDiagonalLength, int NGrainOrientations, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferAR, Buffer2D BufferBR, int BufSizeX, int BufSizeY, int BufSizeZ, int ZBound_Low);
void GhostNodes2D_GPU(int cycle, int id, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewI CellType, ViewF DOCenter, ViewI GrainID, ViewF GrainUnitVector, ViewI GrainOrientation, ViewF DiagonalLength, ViewF CritDiagonalLength, int NGrainOrientations, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferC, Buffer2D BufferD, Buffer2D BufferE, Buffer2D BufferF, Buffer2D BufferG, Buffer2D BufferH, Buffer2D BufferAR, Buffer2D BufferBR, Buffer2D BufferCR, Buffer2D BufferDR, Buffer2D BufferER, Buffer2D BufferFR, Buffer2D BufferGR, Buffer2D BufferHR, int BufSizeX, int BufSizeY, int BufSizeZ, int ZBound_Low);

// Contained in "CAPrint.cpp"
void CollectGrainData(int pid, int np,int nx, int ny,int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI_H GrainID, ViewI_H GrainOrientation, ViewF_H GrainUnitVector,std::string OutputFile, int DecompositionStrategy, int NGrainOrientations, bool* Melted, std::string PathToOutput, bool FilesToPrint[4], double deltax);
void PrintTempValues(int pid, int np, int nx, int ny, int nz, int MyXSlices,int MyYSlices,int ProcessorsInXDirection,int ProcessorsInYDirection, ViewI_H CritTimeStep, int DecompositionStrategy, std::string PathToOutput);
void PrintCT(int pid, int np,int nx, int ny,int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI_H CellType, std::string OutputFile, int DecompositionStrategy);
void PrintOrientations(std::string FName, int nx, int ny, int nz, std::vector <std::vector <std::vector <bool> > > Melted_WholeDomain, std::vector <std::vector <std::vector <int> > > GrainID_WholeDomain, int NGrainOrientations);
void PrintGrainIDs(std::string FName, int nx, int ny, int nz, std::vector <std::vector <std::vector <bool> > > Melted_WholeDomain, std::vector <std::vector <std::vector <int> > > GrainID_WholeDomain);
void PrintGrainIDsForExaConstit(std::string FName, int nx, int ny, int nz, std::vector <std::vector <std::vector <int> > > GrainID_WholeDomain, double deltax);
void PrintParaview(std::string FName, int nx, int ny, int nz, std::vector <std::vector <std::vector <bool> > > Melted_WholeDomain, std::vector <std::vector <std::vector <int> > > GrainID_WholeDomain, ViewI_H GrainOrientation, ViewF_H GrainUnitVector, int NGrainOrientations);
void PrintGrainAreas(std::string FName, double deltax, int nx, int ny, int nz, std::vector <std::vector <std::vector <int> > > GrainID_WholeDomain);
void PrintWeightedGrainAreas(std::string FName, double deltax, int nx, int ny, int nz, std::vector <std::vector <std::vector <int> > > GrainID_WholeDomain);
