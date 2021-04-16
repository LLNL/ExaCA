#ifndef EXACA_INIT_HPP
#define EXACA_INIT_HPP

#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

#include <string>
#include <vector>

// Initializes input parameters, mesh, temperature field, and grain structures for CA simulations
class InputParameters {

  public:
    InputParameters(int id, std::string InputFile);

    void parallelMeshInit(int id, int np, ViewI_H NeighborX, ViewI_H NeighborY, ViewI_H NeighborZ, ViewI2D_H ItList,
                          int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset, int &NeighborRank_North,
                          int &NeighborRank_South, int &NeighborRank_East, int &NeighborRank_West,
                          int &NeighborRank_NorthEast, int &NeighborRank_NorthWest, int &NeighborRank_SouthEast,
                          int &NeighborRank_SouthWest, int &ProcessorsInXDirection, int &ProcessorsInYDirection,
                          float &XMin, float &XMax, float &YMin, float &YMax, float &ZMin, float &ZMax,
                          unsigned int &NumberOfTemperatureDataPoints, float *ZMinLayer, float *ZMaxLayer,
                          int *FirstValue, int *LastValue, std::vector<float> &RawData);
    void tempInit(int id, int layernumber, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset,
                  ViewI_H CritTimeStep, ViewF_H UndercoolingChange, ViewF_H UndercoolingCurrent, float XMin, float YMin,
                  float ZMin, bool *Melted, float *ZMinLayer, float *ZMaxLayer, int &nzActive, int &ZBound_Low,
                  int &ZBound_High, int *FinishTimeStep, ViewI_H LayerID, int *FirstValue, int *LastValue,
                  std::vector<float> RawData);

    // Overall
    int DecompositionStrategy;
    double deltax;
    double deltat;
    bool ExtraWalls;
    int nx;
    int ny;
    int nz;
    bool RemeltingYN;
    std::string SimulationType;
    // CellCapture
    double AConst;
    double BConst;
    double CConst;
    double DConst;
    double FreezingRange;
    double NMax;
    // Nucleation
    double dTN;
    double dTsigma;
    double HT_deltax;
    // Substrate
    double FractSurfaceSitesActive;
    std::string SubstrateFileName;
    float SubstrateGrainSpacing;
    bool UseSubstrateFile;
    // Grain
    std::string GrainOrientationFile;
    // Outputs
    bool(FilesToPrint)[6];
    std::string OutputFile;
    std::string PathToOutput;
    bool PrintFilesYN;
    std::string tempfile;
    int TempFilesInSeries;
    // AM specific
    int LayerHeight;
    int NumberOfLayers;

  private:
    void readFromFile(int id, std::string InputFile);
    // Helper functions for input file read.
    void skipLines(std::ifstream &stream);
    std::string getKey(std::ifstream &stream, std::string &line, std::size_t &colon);
    std::string removeWhitespace(std::string line, std::size_t colon);
    std::string parseInputMultiple(std::ifstream &stream, std::string key1, std::string key2, int &WhichKey);
    bool parseInputBool(std::ifstream &stream, std::string key);
    std::string parseInput(std::ifstream &stream, std::string key);
    // Helper functions for temperature file read.
    void checkTemperatureCoordinateBound(std::string Label, float LowerBound, float UpperBound, float InputValue,
                                         int LineNumber, std::string TemperatureFilename);
    void checkTemperatureDataPoint(std::string Label, float InputValue, int LineNumber,
                                   std::string TemperatureFilename);

    double G;
    double R;
};

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
               int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int id, int np, int NeighborRank_North,
               int NeighborRank_South, int NeighborRank_East, int NeighborRank_West, int NeighborRank_NorthEast,
               int NeighborRank_NorthWest, int NeighborRank_SouthEast, int NeighborRank_SouthWest, ViewI2D_H ItList,
               ViewI_H NeighborX, ViewI_H NeighborY, ViewI_H NeighborZ, ViewI_H GrainOrientation,
               ViewF_H GrainUnitVector, ViewF_H DiagonalLength, ViewI_H CellType, ViewI_H GrainID,
               ViewF_H CritDiagonalLength, ViewF_H DOCenter, ViewI_H CritTimeStep, double deltax, double NMax,
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

#endif
