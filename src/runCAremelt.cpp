// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "runCAremelt.hpp"

#include "CAghostnodes.hpp"
#include "CAinitialize.hpp"
#include "CAprint.hpp"
#include "CAtypes.hpp"
#include "CAupdate.hpp"

#include "mpi.h"

#include <string>
#include <vector>

void RunProgram_Remelt(int id, int np, std::string SimulationType, std::string InputFile) {
    double NuclTime = 0.0, CaptureTime = 0.0, GhostTime = 0.0;
    double StartNuclTime, StartCaptureTime, StartGhostTime;
    double StartTime = MPI_Wtime();

    int nx, ny, nz, DecompositionStrategy, NumberOfLayers, LayerHeight, TempFilesInSeries;
    unsigned int NumberOfTemperatureDataPoints = 0; // Initialized to 0 - updated if/when temperature files are read
    bool ExtraWalls = false; // If simulating a spot melt problem where the side walls are not part of the substrate,
                             // this is changed to true in the input file
    bool PrintFilesYN, UseSubstrateFile;
    bool RemeltingYN = true;
    bool FilesToPrint[6] = {0}; // Which specific files to print are specified in the input file
    float SubstrateGrainSpacing;
    double HT_deltax, deltax, deltat, FractSurfaceSitesActive, G, R, AConst, BConst, CConst, DConst, FreezingRange,
        NMax, dTN, dTsigma;
    std::string SubstrateFileName, tempfile, OutputFile, GrainOrientationFile, PathToOutput;

    // Read input data
    InputReadFromFile(id, InputFile, SimulationType, DecompositionStrategy, AConst, BConst, CConst, DConst,
                      FreezingRange, deltax, NMax, dTN, dTsigma, OutputFile, GrainOrientationFile, tempfile,
                      TempFilesInSeries, ExtraWalls, HT_deltax, deltat, NumberOfLayers, LayerHeight, SubstrateFileName,
                      SubstrateGrainSpacing, UseSubstrateFile, G, R, nx, ny, nz, FractSurfaceSitesActive, PathToOutput,
                      FilesToPrint, PrintFilesYN);

    // Grid decomposition
    int ProcessorsInXDirection, ProcessorsInYDirection;
    // Variables characterizing local processor grids relative to global domain
    int MyXSlices, MyXOffset, MyYSlices, MyYOffset;
    // Variables characterizing process IDs of neighboring MPI ranks on the grid
    // Positive X/Negative X directions are West/East, Positive Y/NegativeY directions are North/South
    int NeighborRank_North, NeighborRank_South, NeighborRank_East, NeighborRank_West, NeighborRank_NorthEast,
        NeighborRank_NorthWest, NeighborRank_SouthEast, NeighborRank_SouthWest;
    // Neighbor lists for cells
    ViewI_H NeighborX_H(Kokkos::ViewAllocateWithoutInitializing("NeighborX"), 26);
    ViewI_H NeighborY_H(Kokkos::ViewAllocateWithoutInitializing("NeighborY"), 26);
    ViewI_H NeighborZ_H(Kokkos::ViewAllocateWithoutInitializing("NeighborZ"), 26);
    ViewI2D_H ItList_H(Kokkos::ViewAllocateWithoutInitializing("ItList"), 9, 26);
    float XMin, YMin, ZMin, XMax, YMax, ZMax; // OpenFOAM simulation bounds (if using OpenFOAM data)
    float *ZMinLayer = new float[NumberOfLayers];
    float *ZMaxLayer = new float[NumberOfLayers];
    int *FinishTimeStep = new int[NumberOfLayers];

    // Data structure for storing raw temperature data from file(s)
    // With no remelting, each data point has 6 values (X, Y, Z coordinates, melting time, liquidus time, and either
    // solidus time OR cooling rate) Initial estimate for size
    std::vector<float> RawData(1000000);
    ViewI DummyView_G(Kokkos::ViewAllocateWithoutInitializing("DG"), 0);
    // Max number of times a cell in each layer will undergo solidification
    ViewI_H MaxSolidificationEvents_H(Kokkos::ViewAllocateWithoutInitializing("NumberOfRemeltEvents"), NumberOfLayers);

    // Contains "NumberOfLayers" values corresponding to the location within "RawData" of the first data element in each
    // temperature file
    int *FirstValue = new int[NumberOfLayers];
    int *LastValue = new int[NumberOfLayers];
    // Initialization of the grid and decomposition, along with deltax and deltat
    // Read in temperature data
    ParallelMeshInit(DecompositionStrategy, MaxSolidificationEvents_H, NeighborX_H, NeighborY_H, NeighborZ_H, ItList_H,
                     SimulationType, id, np, MyXSlices, MyYSlices, MyXOffset, MyYOffset, NeighborRank_North,
                     NeighborRank_South, NeighborRank_East, NeighborRank_West, NeighborRank_NorthEast,
                     NeighborRank_NorthWest, NeighborRank_SouthEast, NeighborRank_SouthWest, deltax, HT_deltax, nx, ny,
                     nz, ProcessorsInXDirection, ProcessorsInYDirection, tempfile, XMin, XMax, YMin, YMax, ZMin, ZMax,
                     FreezingRange, LayerHeight, NumberOfLayers, TempFilesInSeries, NumberOfTemperatureDataPoints,
                     ZMinLayer, ZMaxLayer, FirstValue, LastValue, RawData);
    // By default, the active domain bounds are the same as the global domain bounds
    // For multilayer problems, this is not the case and ZBound_High and nzActive will be adjusted in TempInit to
    // account for only the first layer of solidification
    long int LocalDomainSize = MyXSlices * MyYSlices * nz; // Number of cells on this MPI rank

    // "ZMin" is the global Z coordinate that corresponds to cells at Z = 2 (Z = 0 is the domain's
    // bottom wall, Z = 1 are the active cells just outside of the melt pool) "ZMax" is the global Z coordinate
    // that corresponds to cells at Z = nz-2 (Z = nz-1 is the domain's top wall)
    std::cout << "ZMinLayer = " << ZMinLayer[0] << " ZMaxLayer = " << ZMaxLayer[0] << std::endl;
    int ZBound_Low = round((ZMinLayer[0] - ZMin) / deltax) + 2;
    int nzActive = round((ZMaxLayer[0] - ZMinLayer[0]) / deltax) +
                   1; // (note this doesn't include the 2 rows of wall/active cells at the bottom surface)
    int ZBound_High = ZBound_Low + nzActive;
    int LocalActiveDomainSize = MyXSlices * MyYSlices * nzActive;

    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Mesh initialized, max solidification events in layer 0 is " << MaxSolidificationEvents_H(0)
                  << std::endl;

    // Temperature fields characterized by these variables:
    // A view that holds melting time step, liquidus time step, and cooling rate/time step data for all cells this layer
    ViewF3D_H LayerTimeTempHistory_H(Kokkos::ViewAllocateWithoutInitializing("MaxSEvents"), LocalActiveDomainSize,
                                     MaxSolidificationEvents_H(0), 3);
    // The number of times that each CA cell will undergo solidification during this layer
    ViewI_H NumberOfSolidificationEvents_H(Kokkos::ViewAllocateWithoutInitializing("NumSEvents"),
                                           LocalActiveDomainSize);
    // A counter for the number of times each CA cell has undergone solidification so far this layer
    ViewI_H SolidificationEventCounter_H(Kokkos::ViewAllocateWithoutInitializing("SEventCounter"),
                                         LocalActiveDomainSize);
    // The next time that each cell will melt during this layer
    ViewI_H MeltTimeStep_H(Kokkos::ViewAllocateWithoutInitializing("MeltTimeStep"), LocalDomainSize);
    // The next time each cell will cool below the liquidus during this layer
    ViewI_H CritTimeStep_H(Kokkos::ViewAllocateWithoutInitializing("CritTimeStep"), LocalDomainSize);
    // Marks which layer a cell's final solidification is associated with
    ViewI_H LayerID_H(Kokkos::ViewAllocateWithoutInitializing("LayerID"), LocalDomainSize);
    // The rate of cooling from the liquidus temperature for each cell for this solidification event
    ViewF_H UndercoolingChange_H(Kokkos::ViewAllocateWithoutInitializing("UndercoolingChange"), LocalDomainSize);
    // The undercooling of each cell (> 0 if actively cooling from the liquidus, = 0 otherwise)
    ViewF_H UndercoolingCurrent_H(Kokkos::ViewAllocateWithoutInitializing("UndercoolingCurrent"), LocalDomainSize);
    // Whether or not a given cell underwent melting at any point during the simulation
    bool *Melted = new bool[LocalDomainSize];

    // Initialize the temperature fields -
    TempInitRemelt(0, id, MyXSlices, MyYSlices, nz, MyXOffset, MyYOffset, deltax, deltat, FreezingRange,
                   LayerTimeTempHistory_H, MaxSolidificationEvents_H, NumberOfSolidificationEvents_H,
                   SolidificationEventCounter_H, MeltTimeStep_H, CritTimeStep_H, UndercoolingChange_H,
                   UndercoolingCurrent_H, XMin, YMin, Melted, ZMinLayer, LayerHeight, nzActive, ZBound_Low, ZBound_High,
                   FinishTimeStep, LayerID_H, FirstValue, LastValue, RawData);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Done with temperature field initialization, active domain size is " << nzActive << " out of "
                  << nz << " cells in the Z direction" << std::endl;

    // PrintTempValues(id,np,nx,ny,nz, MyXSlices, MyYSlices, ProcessorsInXDirection, ProcessorsInYDirection,
    // CritTimeStep_H, DecompositionStrategy,PathToOutput);

    int NGrainOrientations = 10000; // Number of grain orientations considered in the simulation
    ViewF_H GrainUnitVector_H(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), 9 * NGrainOrientations);
    ViewI_H GrainOrientation_H(Kokkos::ViewAllocateWithoutInitializing("GrainOrientation"), NGrainOrientations);

    // Initialize grain orientations
    OrientationInit(id, NGrainOrientations, GrainOrientation_H, GrainUnitVector_H, GrainOrientationFile);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Done with orientation initialization " << std::endl;

    // CA cell variables
    ViewI_H GrainID_H(Kokkos::ViewAllocateWithoutInitializing("GrainID"), LocalDomainSize);
    ViewI_H CellType_H(Kokkos::ViewAllocateWithoutInitializing("CellType"), LocalDomainSize);

    // Variables characterizing the active cell region within each rank's grid
    ViewF_H DiagonalLength_H(Kokkos::ViewAllocateWithoutInitializing("DiagonalLength"), LocalActiveDomainSize);
    ViewF_H CritDiagonalLength_H(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength"),
                                 26 * LocalActiveDomainSize);
    ViewF_H DOCenter_H(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * LocalActiveDomainSize);

    // Initialize the grain structure - for either a constrained solidification problem, using a substrate from a file,
    // or generating a substrate using the existing CA algorithm
    int PossibleNuclei_ThisRank;
    int NextLayer_FirstNucleatedGrainID = -1;
    if (UseSubstrateFile)
        SubstrateInit_FromFile(SubstrateFileName, RemeltingYN, nz, MyXSlices, MyYSlices, MyXOffset, MyYOffset, id,
                               CritTimeStep_H, GrainID_H);
    else
        SubstrateInit_FromGrainSpacing(SubstrateGrainSpacing, RemeltingYN, nx, ny, nz, nzActive, MyXSlices, MyYSlices,
                                       MyXOffset, MyYOffset, LocalActiveDomainSize, id, np, deltax, GrainID_H,
                                       CritTimeStep_H);

    // Estimate number of nuclei on each rank (resize later when "PossibleNuclei_ThisRank" is known)
    ViewI_H NucleationTimes_H(Kokkos::ViewAllocateWithoutInitializing("NucleationTimes"), LocalActiveDomainSize);
    ViewI_H NucleiLocation_H(Kokkos::ViewAllocateWithoutInitializing("NucleiLocation"), LocalActiveDomainSize);
    ViewI_H NucleiGrainID_H(Kokkos::ViewAllocateWithoutInitializing("NucleiGrainID"), LocalActiveDomainSize);
    GrainNucleiInitRemelt(0, LocalActiveDomainSize, MyXSlices, MyYSlices, nzActive, id, np, DiagonalLength_H,
                          CellType_H, CritTimeStep_H, NumberOfSolidificationEvents_H, LayerTimeTempHistory_H, deltax,
                          NMax, dTN, dTsigma, NextLayer_FirstNucleatedGrainID, PossibleNuclei_ThisRank,
                          NucleationTimes_H, NucleiLocation_H, NucleiGrainID_H, ZBound_Low);
    Kokkos::resize(NucleationTimes_H, PossibleNuclei_ThisRank);
    Kokkos::resize(NucleiLocation_H, PossibleNuclei_ThisRank);
    Kokkos::resize(NucleiGrainID_H, PossibleNuclei_ThisRank);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Grain struct initialized" << std::endl;

    // Normalize solidification parameters
    AConst = AConst * deltat / deltax;
    BConst = BConst * deltat / deltax;
    CConst = CConst * deltat / deltax;
    int cycle;

    // Buffers for ghost node data (fixed size)
    int BufSizeX, BufSizeY, BufSizeZ;
    if (DecompositionStrategy == 1) {
        BufSizeX = MyXSlices;
        BufSizeY = 0;
        BufSizeZ = nzActive;
    }
    else {
        BufSizeX = MyXSlices - 2;
        BufSizeY = MyYSlices - 2;
        BufSizeZ = nzActive;
    }
    Buffer2D BufferSouthSend(Kokkos::ViewAllocateWithoutInitializing("BufferSouthSend"), BufSizeX * BufSizeZ, 5);
    Buffer2D BufferNorthSend(Kokkos::ViewAllocateWithoutInitializing("BufferNorthSend"), BufSizeX * BufSizeZ, 5);
    Buffer2D BufferEastSend(Kokkos::ViewAllocateWithoutInitializing("BufferEastSend"), BufSizeY * BufSizeZ, 5);
    Buffer2D BufferWestSend(Kokkos::ViewAllocateWithoutInitializing("BufferWestSend"), BufSizeY * BufSizeZ, 5);
    Buffer2D BufferNorthEastSend(Kokkos::ViewAllocateWithoutInitializing("BufferNorthEastSend"), BufSizeZ, 5);
    Buffer2D BufferNorthWestSend(Kokkos::ViewAllocateWithoutInitializing("BufferNorthWestSend"), BufSizeZ, 5);
    Buffer2D BufferSouthEastSend(Kokkos::ViewAllocateWithoutInitializing("BufferSouthEastSend"), BufSizeZ, 5);
    Buffer2D BufferSouthWestSend(Kokkos::ViewAllocateWithoutInitializing("BufferSouthWestSend"), BufSizeZ, 5);
    Buffer2D BufferSouthRecv(Kokkos::ViewAllocateWithoutInitializing("BufferSouthRecv"), BufSizeX * BufSizeZ, 5);
    Buffer2D BufferNorthRecv(Kokkos::ViewAllocateWithoutInitializing("BufferNorthRecv"), BufSizeX * BufSizeZ, 5);
    Buffer2D BufferEastRecv(Kokkos::ViewAllocateWithoutInitializing("BufferEastRecv"), BufSizeY * BufSizeZ, 5);
    Buffer2D BufferWestRecv(Kokkos::ViewAllocateWithoutInitializing("BufferWestRecv"), BufSizeY * BufSizeZ, 5);
    Buffer2D BufferNorthEastRecv(Kokkos::ViewAllocateWithoutInitializing("BufferNorthEastRecv"), BufSizeZ, 5);
    Buffer2D BufferNorthWestRecv(Kokkos::ViewAllocateWithoutInitializing("BufferNorthWestRecv"), BufSizeZ, 5);
    Buffer2D BufferSouthEastRecv(Kokkos::ViewAllocateWithoutInitializing("BufferSouthEastRecv"), BufSizeZ, 5);
    Buffer2D BufferSouthWestRecv(Kokkos::ViewAllocateWithoutInitializing("BufferSouthWestRecv"), BufSizeZ, 5);

    // Copy view data to GPU
    using memory_space = Kokkos::DefaultExecutionSpace::memory_space;
    ViewI GrainID_G = Kokkos::create_mirror_view_and_copy(memory_space(), GrainID_H);
    ViewI CellType_G = Kokkos::create_mirror_view_and_copy(memory_space(), CellType_H);
    ViewF DiagonalLength_G = Kokkos::create_mirror_view_and_copy(memory_space(), DiagonalLength_H);
    ViewF CritDiagonalLength_G = Kokkos::create_mirror_view_and_copy(memory_space(), CritDiagonalLength_H);
    ViewF DOCenter_G = Kokkos::create_mirror_view_and_copy(memory_space(), DOCenter_H);
    ViewI CritTimeStep_G = Kokkos::create_mirror_view_and_copy(memory_space(), CritTimeStep_H);
    ViewI LayerID_G = Kokkos::create_mirror_view_and_copy(memory_space(), LayerID_H);
    ViewF UndercoolingChange_G = Kokkos::create_mirror_view_and_copy(memory_space(), UndercoolingChange_H);
    ViewF UndercoolingCurrent_G = Kokkos::create_mirror_view_and_copy(memory_space(), UndercoolingCurrent_H);
    ViewI NucleiGrainID_G = Kokkos::create_mirror_view_and_copy(memory_space(), NucleiGrainID_H);
    ViewI NucleiLocation_G = Kokkos::create_mirror_view_and_copy(memory_space(), NucleiLocation_H);
    ViewI NucleationTimes_G = Kokkos::create_mirror_view_and_copy(memory_space(), NucleationTimes_H);
    ViewI NeighborX_G = Kokkos::create_mirror_view_and_copy(memory_space(), NeighborX_H);
    ViewI NeighborY_G = Kokkos::create_mirror_view_and_copy(memory_space(), NeighborY_H);
    ViewI NeighborZ_G = Kokkos::create_mirror_view_and_copy(memory_space(), NeighborZ_H);
    ViewI2D ItList_G = Kokkos::create_mirror_view_and_copy(memory_space(), ItList_H);
    ViewI GrainOrientation_G = Kokkos::create_mirror_view_and_copy(memory_space(), GrainOrientation_H);
    ViewF GrainUnitVector_G = Kokkos::create_mirror_view_and_copy(memory_space(), GrainUnitVector_H);
    ViewI MaxSolidificationEvents_G = Kokkos::create_mirror_view_and_copy(memory_space(), MaxSolidificationEvents_H);
    ViewI NumberOfSolidificationEvents_G =
        Kokkos::create_mirror_view_and_copy(memory_space(), NumberOfSolidificationEvents_H);
    ViewF3D LayerTimeTempHistory_G = Kokkos::create_mirror_view_and_copy(memory_space(), LayerTimeTempHistory_H);
    ViewI SolidificationEventCounter_G =
        Kokkos::create_mirror_view_and_copy(memory_space(), SolidificationEventCounter_H);
    ViewI MeltTimeStep_G = Kokkos::create_mirror_view_and_copy(memory_space(), MeltTimeStep_H);

    // The next time each cell will cool below the liquidus during this layer
    // Steering Vector
    ViewI SteeringVector(Kokkos::ViewAllocateWithoutInitializing("SteeringVector"), LocalActiveDomainSize);
    ViewI_H numSteer_H(Kokkos::ViewAllocateWithoutInitializing("SteeringVectorSize"), 1);
    numSteer_H(0) = 0;
    ViewI numSteer_G = Kokkos::create_mirror_view_and_copy(memory_space(), numSteer_H);

    if (np > 1) {
        // Ghost nodes for initial microstructure state and cell types
        GhostNodesInit(id, np, DecompositionStrategy, NeighborRank_North, NeighborRank_South, NeighborRank_East,
                       NeighborRank_West, NeighborRank_NorthEast, NeighborRank_NorthWest, NeighborRank_SouthEast,
                       NeighborRank_SouthWest, MyXSlices, MyYSlices, MyXOffset, MyYOffset, ZBound_Low, nzActive,
                       LocalActiveDomainSize, NGrainOrientations, NeighborX_G, NeighborY_G, NeighborZ_G,
                       GrainUnitVector_G, GrainOrientation_G, GrainID_G, CellType_G, DOCenter_G, DiagonalLength_G,
                       CritDiagonalLength_G);
    }

    double InitTime = MPI_Wtime() - StartTime;
    if (id == 0)
        std::cout << "\nData initialized: Time spent: " << InitTime << " s" << std::endl;
    cycle = 0;

    for (int layernumber = 0; layernumber < NumberOfLayers; layernumber++) {

        int nn = 0; // Counter for the number of nucleation events
        int XSwitch = 0;
        double LayerTime1 = MPI_Wtime();

        // Loop continues until all liquid cells claimed by solid grains
        do {
            cycle++;

            // Update cells on GPU - undercooling and diagonal length updates, nucleation
            StartNuclTime = MPI_Wtime();
            Nucleation(RemeltingYN, MyXSlices, MyYSlices, MyXOffset, MyYOffset, cycle, nn, CellType_G, NucleiLocation_G,
                       NucleationTimes_G, NucleiGrainID_G, GrainID_G, GrainOrientation_G, DOCenter_G, NeighborX_G,
                       NeighborY_G, NeighborZ_G, GrainUnitVector_G, CritDiagonalLength_G, DiagonalLength_G,
                       NGrainOrientations, PossibleNuclei_ThisRank, ZBound_Low, layernumber, LayerID_G);
            NuclTime += MPI_Wtime() - StartNuclTime;

            // Update cells on GPU - new active cells, solidification of old active cells
            StartCaptureTime = MPI_Wtime();
            CellCapture(np, cycle, RemeltingYN, DecompositionStrategy, LocalActiveDomainSize, id, MyXSlices, MyYSlices,
                        AConst, BConst, CConst, DConst, MyXOffset, MyYOffset, ItList_G, NeighborX_G, NeighborY_G,
                        NeighborZ_G, CritTimeStep_G, UndercoolingCurrent_G, UndercoolingChange_G, GrainUnitVector_G,
                        CritDiagonalLength_G, DiagonalLength_G, GrainOrientation_G, CellType_G, DOCenter_G, GrainID_G,
                        NGrainOrientations, BufferWestSend, BufferEastSend, BufferNorthSend, BufferSouthSend,
                        BufferNorthEastSend, BufferNorthWestSend, BufferSouthEastSend, BufferSouthWestSend, BufSizeX,
                        BufSizeY, ZBound_Low, nzActive, nz, layernumber, LayerID_G, SteeringVector, numSteer_G,
                        numSteer_H, MeltTimeStep_G, SolidificationEventCounter_G, NumberOfSolidificationEvents_G,
                        LayerTimeTempHistory_G);
            CaptureTime += MPI_Wtime() - StartCaptureTime;

            if (np > 1) {
                // Update ghost nodes
                StartGhostTime = MPI_Wtime();
                if (DecompositionStrategy == 1)
                    GhostNodes1D(cycle, id, NeighborRank_North, NeighborRank_South, MyXSlices, MyYSlices, MyXOffset,
                                 MyYOffset, NeighborX_G, NeighborY_G, NeighborZ_G, CellType_G, DOCenter_G, GrainID_G,
                                 GrainUnitVector_G, GrainOrientation_G, DiagonalLength_G, CritDiagonalLength_G,
                                 NGrainOrientations, BufferNorthSend, BufferSouthSend, BufferNorthRecv, BufferSouthRecv,
                                 BufSizeX, BufSizeY, BufSizeZ, ZBound_Low);
                else
                    GhostNodes2D(cycle, id, NeighborRank_North, NeighborRank_South, NeighborRank_East,
                                 NeighborRank_West, NeighborRank_NorthEast, NeighborRank_NorthWest,
                                 NeighborRank_SouthEast, NeighborRank_SouthWest, MyXSlices, MyYSlices, MyXOffset,
                                 MyYOffset, NeighborX_G, NeighborY_G, NeighborZ_G, CellType_G, DOCenter_G, GrainID_G,
                                 GrainUnitVector_G, GrainOrientation_G, DiagonalLength_G, CritDiagonalLength_G,
                                 NGrainOrientations, BufferWestSend, BufferEastSend, BufferNorthSend, BufferSouthSend,
                                 BufferNorthEastSend, BufferNorthWestSend, BufferSouthEastSend, BufferSouthWestSend,
                                 BufferWestRecv, BufferEastRecv, BufferNorthRecv, BufferSouthRecv, BufferNorthEastRecv,
                                 BufferNorthWestRecv, BufferSouthEastRecv, BufferSouthWestRecv, BufSizeX, BufSizeY,
                                 BufSizeZ, ZBound_Low);
                GhostTime += MPI_Wtime() - StartGhostTime;
            }

            if (cycle % 1000 == 0) {
                 IntermediateOutputAndCheck_Remelt(id, cycle, MyXSlices, MyYSlices, LocalActiveDomainSize,
                                                    nn, XSwitch, CellType_G, MeltTimeStep_G,
                                                    FinishTimeStep, layernumber, NumberOfLayers,
                                                   ZBound_Low, LayerID_G);
            }

        } while (XSwitch == 0);

        if (layernumber != NumberOfLayers - 1) {

            ZBound_Low = round((ZMinLayer[layernumber + 1] - ZMin) / deltax) + 2;
            nzActive = round((ZMaxLayer[layernumber + 1] - ZMinLayer[layernumber + 1]) / deltax) +
                       1; // (note this doesn't include the 2 rows of wall/active cells at the bottom surface)
            ZBound_High = ZBound_Low + nzActive;
            LocalActiveDomainSize = MyXSlices * MyYSlices * nzActive;

            // Initialize temperature data for next layer
            TempInitRemelt(layernumber + 1, id, MyXSlices, MyYSlices, nz, MyXOffset, MyYOffset, deltax, deltat,
                           FreezingRange, LayerTimeTempHistory_H, MaxSolidificationEvents_H,
                           NumberOfSolidificationEvents_H, SolidificationEventCounter_H, MeltTimeStep_H, CritTimeStep_H,
                           UndercoolingChange_H, UndercoolingCurrent_H, XMin, YMin, Melted, ZMinLayer, LayerHeight,
                           nzActive, ZBound_Low, ZBound_High, FinishTimeStep, LayerID_H, FirstValue, LastValue,
                           RawData);

            // Re-initialize solid cells (part of this layer that will melt/solidify) and wall cells (cells that are
            // ignored in this layer) Estimate number of nuclei on each rank (resize later when
            // "PossibleNuclei_ThisRank" is known)
            Kokkos::resize(NucleationTimes_H, LocalActiveDomainSize);
            Kokkos::resize(NucleiLocation_H, LocalActiveDomainSize);
            Kokkos::resize(NucleiGrainID_H, LocalActiveDomainSize);

            GrainNucleiInitRemelt(layernumber + 1, LocalActiveDomainSize, MyXSlices, MyYSlices, nzActive, id, np,
                                  DiagonalLength_H, CellType_H, CritTimeStep_H, NumberOfSolidificationEvents_H,
                                  LayerTimeTempHistory_H, deltax, NMax, dTN, dTsigma, NextLayer_FirstNucleatedGrainID,
                                  PossibleNuclei_ThisRank, NucleationTimes_H, NucleiLocation_H, NucleiGrainID_H,
                                  ZBound_Low);

            Kokkos::resize(NucleationTimes_H, PossibleNuclei_ThisRank);
            Kokkos::resize(NucleiLocation_H, PossibleNuclei_ThisRank);
            Kokkos::resize(NucleiGrainID_H, PossibleNuclei_ThisRank);
            Kokkos::resize(NucleationTimes_G, PossibleNuclei_ThisRank);
            Kokkos::resize(NucleiLocation_G, PossibleNuclei_ThisRank);
            Kokkos::resize(NucleiGrainID_G, PossibleNuclei_ThisRank);
            // Resize steering vector as LocalActiveDomainSize may have changed
            Kokkos::resize(SteeringVector, LocalActiveDomainSize);

            // Resize active cell data structures
            Kokkos::resize(DiagonalLength_G, LocalActiveDomainSize);
            Kokkos::resize(DOCenter_G, 3 * LocalActiveDomainSize);
            Kokkos::resize(CritDiagonalLength_G, 26 * LocalActiveDomainSize);

            Kokkos::resize(BufferNorthSend, BufSizeX * BufSizeZ, 5);
            Kokkos::resize(BufferSouthSend, BufSizeX * BufSizeZ, 5);
            Kokkos::resize(BufferEastSend, BufSizeY * BufSizeZ, 5);
            Kokkos::resize(BufferWestSend, BufSizeY * BufSizeZ, 5);
            Kokkos::resize(BufferNorthEastSend, BufSizeZ, 5);
            Kokkos::resize(BufferNorthWestSend, BufSizeZ, 5);
            Kokkos::resize(BufferSouthEastSend, BufSizeZ, 5);
            Kokkos::resize(BufferSouthWestSend, BufSizeZ, 5);

            Kokkos::resize(BufferNorthRecv, BufSizeX * BufSizeZ, 5);
            Kokkos::resize(BufferSouthRecv, BufSizeX * BufSizeZ, 5);
            Kokkos::resize(BufferEastRecv, BufSizeY * BufSizeZ, 5);
            Kokkos::resize(BufferWestRecv, BufSizeY * BufSizeZ, 5);
            Kokkos::resize(BufferNorthEastRecv, BufSizeZ, 5);
            Kokkos::resize(BufferNorthWestRecv, BufSizeZ, 5);
            Kokkos::resize(BufferSouthEastRecv, BufSizeZ, 5);
            Kokkos::resize(BufferSouthWestRecv, BufSizeZ, 5);

            MPI_Barrier(MPI_COMM_WORLD);
            if (id == 0)
                std::cout << "Resize executed" << std::endl;

            // Update active cell data structures for simulation of next layer
            LayerSetup(MyXSlices, MyYSlices, MyXOffset, MyYOffset, LocalActiveDomainSize, GrainOrientation_G,
                       NGrainOrientations, GrainUnitVector_G, NeighborX_G, NeighborY_G, NeighborZ_G, DiagonalLength_G,
                       CellType_G, GrainID_G, CritDiagonalLength_G, DOCenter_G, DecompositionStrategy, BufferWestRecv,
                       BufferEastSend, BufferNorthSend, BufferSouthSend, BufferNorthEastSend, BufferNorthWestSend,
                       BufferSouthEastSend, BufferSouthWestSend, BufferWestRecv, BufferEastRecv, BufferNorthRecv,
                       BufferSouthRecv, BufferNorthEastRecv, BufferNorthWestRecv, BufferSouthEastRecv,
                       BufferSouthWestRecv, ZBound_Low);

            if (id == 0)
                std::cout << "New layer setup, GN dimensions are " << BufSizeX << " " << BufSizeY << " " << BufSizeZ
                          << std::endl;
            std::cout << "New layer setup, GN dimensions are " << BufSizeX << " " << BufSizeY << " " << BufSizeZ
                      << std::endl;

            // Deep copy updated views from host
            Kokkos::deep_copy(CellType_G, CellType_H);
            Kokkos::deep_copy(CritTimeStep_G, CritTimeStep_H);
            Kokkos::deep_copy(MeltTimeStep_G, MeltTimeStep_H);
            Kokkos::deep_copy(UndercoolingChange_G, UndercoolingChange_H);
            Kokkos::deep_copy(UndercoolingCurrent_G, UndercoolingCurrent_H);
            Kokkos::deep_copy(LayerTimeTempHistory_G, LayerTimeTempHistory_H);
            Kokkos::deep_copy(LayerID_G, LayerID_H);
            Kokkos::deep_copy(NumberOfSolidificationEvents_G, NumberOfSolidificationEvents_H);
            Kokkos::deep_copy(SolidificationEventCounter_G, SolidificationEventCounter_H);
            Kokkos::deep_copy(NucleationTimes_G, NucleationTimes_H);
            Kokkos::deep_copy(NucleiLocation_G, NucleiLocation_H);
            Kokkos::deep_copy(NucleiGrainID_G, NucleiGrainID_H);

            // Update ghost nodes for grain locations and attributes
            MPI_Barrier(MPI_COMM_WORLD);
            if (id == 0)
                std::cout << "New layer ghost nodes initialized" << std::endl;
            if (np > 1) {
                GhostNodesInit(id, np, DecompositionStrategy, NeighborRank_North, NeighborRank_South, NeighborRank_East,
                               NeighborRank_West, NeighborRank_NorthEast, NeighborRank_NorthWest,
                               NeighborRank_SouthEast, NeighborRank_SouthWest, MyXSlices, MyYSlices, MyXOffset,
                               MyYOffset, ZBound_Low, nzActive, LocalActiveDomainSize, NGrainOrientations, NeighborX_G,
                               NeighborY_G, NeighborZ_G, GrainUnitVector_G, GrainOrientation_G, GrainID_G, CellType_G,
                               DOCenter_G, DiagonalLength_G, CritDiagonalLength_G);
            }
            XSwitch = 0;
            MPI_Barrier(MPI_COMM_WORLD);
            double LayerTime2 = MPI_Wtime();
            cycle = 0;
            if (id == 0)
                std::cout << "Time for layer number " << layernumber << " was " << LayerTime2 - LayerTime1
                          << " s, starting layer " << layernumber + 1 << std::endl;
        }
        else {
            MPI_Barrier(MPI_COMM_WORLD);
            double LayerTime2 = MPI_Wtime();
            if (id == 0)
                std::cout << "Time for final layer was " << LayerTime2 - LayerTime1 << " s" << std::endl;
        }
    }

    double RunTime = MPI_Wtime() - InitTime;

    // Copy GPU results for GrainID back to CPU for printing to file(s)
    Kokkos::deep_copy(GrainID_H, GrainID_G);
    Kokkos::deep_copy(CellType_H, CellType_G);

    MPI_Barrier(MPI_COMM_WORLD);
    if (PrintFilesYN) {
        if (id == 0)
            std::cout << "Collecting data on rank 0 and printing to files" << std::endl;
        CollectGrainData(id, np, nx, ny, nz, MyXSlices, MyYSlices, ProcessorsInXDirection, ProcessorsInYDirection,
                         GrainID_H, GrainOrientation_H, GrainUnitVector_H, OutputFile, DecompositionStrategy,
                         NGrainOrientations, Melted, PathToOutput, FilesToPrint, deltax);
    }
    else {
        if (id == 0)
            std::cout << "No output files to be printed, exiting program" << std::endl;
    }

    double OutTime = MPI_Wtime() - RunTime - InitTime;
    double InitMaxTime, InitMinTime, OutMaxTime, OutMinTime = 0.0;
    double NuclMaxTime, NuclMinTime, CaptureMaxTime, CaptureMinTime, GhostMaxTime, GhostMinTime = 0.0;
    MPI_Allreduce(&InitTime, &InitMaxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&InitTime, &InitMinTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&NuclTime, &NuclMaxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&NuclTime, &NuclMinTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&CaptureTime, &CaptureMaxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&CaptureTime, &CaptureMinTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&GhostTime, &GhostMaxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&GhostTime, &GhostMinTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&OutTime, &OutMaxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&OutTime, &OutMinTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    if (id == 0) {
        std::cout << "===================================================================================" << std::endl;
        std::cout << "Having run with = " << np << " processors" << std::endl;
        std::cout << "Output written at cycle = " << cycle << std::endl;
        std::cout << "Total time = " << InitTime + RunTime + OutTime << std::endl;
        std::cout << "Time spent initializing data = " << InitTime << " s" << std::endl;
        std::cout << "Time spent performing CA calculations = " << RunTime << " s" << std::endl;
        std::cout << "Time spent collecting and printing output data = " << OutTime << " s\n" << std::endl;

        std::cout << "Max/min rank time initializing data  = " << InitMaxTime << " / " << InitMinTime << " s"
                  << std::endl;
        std::cout << "Max/min rank time in CA nucleation   = " << NuclMaxTime << " / " << NuclMinTime << " s"
                  << std::endl;
        std::cout << "Max/min rank time in CA cell capture = " << CaptureMaxTime << " / " << CaptureMinTime << " s"
                  << std::endl;
        std::cout << "Max/min rank time in CA ghosting     = " << GhostMaxTime << " / " << GhostMinTime << " s"
                  << std::endl;
        std::cout << "Max/min rank time exporting data     = " << OutMaxTime << " / " << OutMinTime << " s\n"
                  << std::endl;

        std::cout << "===================================================================================" << std::endl;
    }
}
