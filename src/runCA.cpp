// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "runCA.hpp"

#include "CAghostnodes.hpp"
#include "CAinitialize.hpp"
#include "CAprint.hpp"
#include "CAtypes.hpp"
#include "CAupdate.hpp"

#include "mpi.h"

#include <string>
#include <vector>

void RunProgram_Reduced(int id, int np, std::string InputFile) {
    double NuclTime = 0.0, CaptureTime = 0.0, GhostTime = 0.0;
    double StartNuclTime, StartCaptureTime, StartGhostTime;
    double StartTime = MPI_Wtime();

    int nx, ny, nz, DecompositionStrategy, NumberOfLayers, LayerHeight, TempFilesInSeries;
    int NSpotsX, NSpotsY, SpotOffset, SpotRadius;
    unsigned int NumberOfTemperatureDataPoints = 0; // Initialized to 0 - updated if/when temperature files are read
    bool ExtraWalls = false; // If simulating a spot melt problem where the side walls are not part of the substrate,
                             // this is changed to true in the input file
    int PrintDebug;
    bool PrintMisorientation, PrintFullOutput, RemeltingYN, UseSubstrateFile;
    float SubstrateGrainSpacing;
    double HT_deltax, deltax, deltat, FractSurfaceSitesActive, G, R, AConst, BConst, CConst, DConst, FreezingRange,
        NMax, dTN, dTsigma;
    std::string SubstrateFileName, temppath, tempfile, SimulationType, OutputFile, GrainOrientationFile, PathToOutput;
    std::vector<std::string> temp_paths;

    // Read input data
    InputReadFromFile(id, InputFile, SimulationType, DecompositionStrategy, AConst, BConst, CConst, DConst,
                      FreezingRange, deltax, NMax, dTN, dTsigma, OutputFile, GrainOrientationFile, temppath, tempfile,
                      TempFilesInSeries, temp_paths, ExtraWalls, HT_deltax, RemeltingYN, deltat, NumberOfLayers,
                      LayerHeight, SubstrateFileName, SubstrateGrainSpacing, UseSubstrateFile, G, R, nx, ny, nz,
                      FractSurfaceSitesActive, PathToOutput, PrintDebug, PrintMisorientation, PrintFullOutput, NSpotsX,
                      NSpotsY, SpotOffset, SpotRadius);

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
    // Store data as double - needed for small time steps to resolve local differences in solidification conditions
    // With no remelting, each data point has 5 values (X, Y, Z coordinates, liquidus time, and either solidus time OR
    // cooling rate) Initial estimate for size
    std::vector<double> RawData(1000000);

    // Contains "NumberOfLayers" values corresponding to the location within "RawData" of the first data element in each
    // temperature file
    int *FirstValue = new int[NumberOfLayers];
    int *LastValue = new int[NumberOfLayers];
    // Initialization of the grid and decomposition, along with deltax and deltat
    // Read in temperature data
    ParallelMeshInit(DecompositionStrategy, NeighborX_H, NeighborY_H, NeighborZ_H, ItList_H, SimulationType, id, np,
                     MyXSlices, MyYSlices, MyXOffset, MyYOffset, NeighborRank_North, NeighborRank_South,
                     NeighborRank_East, NeighborRank_West, NeighborRank_NorthEast, NeighborRank_NorthWest,
                     NeighborRank_SouthEast, NeighborRank_SouthWest, deltax, HT_deltax, nx, ny, nz,
                     ProcessorsInXDirection, ProcessorsInYDirection, temp_paths, XMin, XMax, YMin, YMax, ZMin, ZMax,
                     FreezingRange, LayerHeight, NumberOfLayers, TempFilesInSeries, NumberOfTemperatureDataPoints,
                     ZMinLayer, ZMaxLayer, FirstValue, LastValue, RawData);

    long int LocalDomainSize = MyXSlices * MyYSlices * nz; // Number of cells on this MPI rank

    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Mesh initialized" << std::endl;

    // Temperature fields characterized by these variables:
    ViewI_H CritTimeStep_H(Kokkos::ViewAllocateWithoutInitializing("CritTimeStep"), LocalDomainSize);
    ViewI_H LayerID_H(Kokkos::ViewAllocateWithoutInitializing("LayerID"), LocalDomainSize);
    ViewF_H UndercoolingChange_H(Kokkos::ViewAllocateWithoutInitializing("UndercoolingChange"), LocalDomainSize);
    ViewF_H UndercoolingCurrent_H(Kokkos::ViewAllocateWithoutInitializing("UndercoolingCurrent"), LocalDomainSize);
    bool *Melted = new bool[LocalDomainSize];

    // By default, the active domain bounds are the same as the global domain bounds
    // For multilayer problems, this is not the case and ZBound_High and nzActive will be adjusted in TempInit to
    // account for only the first layer of solidification
    int ZBound_Low = -1;
    int ZBound_High = -1;
    int nzActive = -1;

    // Initialize the temperature fields
    if (SimulationType == "R") {
        // input temperature data from files using reduced/sparse data format
        TempInit_Reduced(id, MyXSlices, MyYSlices, MyXOffset, MyYOffset, deltax, HT_deltax, deltat, nx, ny, nz,
                         CritTimeStep_H, UndercoolingChange_H, UndercoolingCurrent_H, XMin, YMin, ZMin, Melted,
                         ZMinLayer, ZMaxLayer, LayerHeight, NumberOfLayers, nzActive, ZBound_Low, ZBound_High,
                         FinishTimeStep, FreezingRange, LayerID_H, FirstValue, LastValue, RawData);
    }
    else if (SimulationType == "S") {
        // spot melt array test problem
        TempInit_SpotMelt(G, R, SimulationType, id, MyXSlices, MyYSlices, MyXOffset, MyYOffset, deltax, deltat, nz,
                          CritTimeStep_H, UndercoolingChange_H, UndercoolingCurrent_H, Melted, LayerHeight,
                          NumberOfLayers, nzActive, ZBound_Low, ZBound_High, FreezingRange, LayerID_H, NSpotsX, NSpotsY,
                          SpotRadius, SpotOffset);
    }
    else if (SimulationType == "C") {
        // directional/constrained solidification test problem
        TempInit_DirSolidification(G, R, id, MyXSlices, MyYSlices, MyXOffset, MyYOffset, deltax, deltat, nx, ny, nz,
                                   CritTimeStep_H, UndercoolingChange_H, UndercoolingCurrent_H, Melted, nzActive,
                                   ZBound_Low, ZBound_High, LayerID_H);
    }
    // Delete temporary data structure for temperature data read
    RawData.clear();
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Done with temperature field initialization, active domain size is " << nzActive << " out of "
                  << nz << " cells in the Z direction" << std::endl;

    int LocalActiveDomainSize = MyXSlices * MyYSlices * nzActive; // Number of active cells on this MPI rank

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
    int PossibleNuclei_ThisRank, NextLayer_FirstNucleatedGrainID;
    if (SimulationType == "C") {
        SubstrateInit_ConstrainedGrowth(FractSurfaceSitesActive, MyXSlices, MyYSlices, nx, ny, nz, MyXOffset, MyYOffset,
                                        id, np, CellType_H, GrainID_H);
    }
    else {
        if (UseSubstrateFile)
            SubstrateInit_FromFile(SubstrateFileName, nz, MyXSlices, MyYSlices, MyXOffset, MyYOffset, id,
                                   CritTimeStep_H, GrainID_H);
        else
            SubstrateInit_FromGrainSpacing(SubstrateGrainSpacing, nx, ny, nz, nzActive, MyXSlices, MyYSlices, MyXOffset,
                                           MyYOffset, LocalActiveDomainSize, id, np, deltax, GrainID_H, CritTimeStep_H);
        ActiveCellWallInit(id, MyXSlices, MyYSlices, nx, ny, nz, MyXOffset, MyYOffset, CellType_H, GrainID_H,
                           CritTimeStep_H, ItList_H, NeighborX_H, NeighborY_H, NeighborZ_H, UndercoolingChange_H,
                           ExtraWalls);
    }
    GrainInit(-1, NGrainOrientations, DecompositionStrategy, nz, LocalActiveDomainSize, MyXSlices, MyYSlices, MyXOffset,
              MyYOffset, id, np, NeighborRank_North, NeighborRank_South, NeighborRank_East, NeighborRank_West,
              NeighborRank_NorthEast, NeighborRank_NorthWest, NeighborRank_SouthEast, NeighborRank_SouthWest, ItList_H,
              NeighborX_H, NeighborY_H, NeighborZ_H, GrainOrientation_H, GrainUnitVector_H, DiagonalLength_H,
              CellType_H, GrainID_H, CritDiagonalLength_H, DOCenter_H, CritTimeStep_H, deltax, NMax,
              NextLayer_FirstNucleatedGrainID, PossibleNuclei_ThisRank, ZBound_High, ZBound_Low);

    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Grain struct initialized" << std::endl;

    ViewI_H NucleationTimes_H(Kokkos::ViewAllocateWithoutInitializing("NucleationTimes"), PossibleNuclei_ThisRank);
    ViewI_H NucleiLocation_H(Kokkos::ViewAllocateWithoutInitializing("NucleiLocation"), PossibleNuclei_ThisRank);

    // Update nuclei on ghost nodes, fill in nucleation data structures, and assign nucleation undercooling values to
    // potential nucleation events
    if (id == 0)
        std::cout << " Possible nucleation events (rank: # events): " << std::endl;
    NucleiInit(DecompositionStrategy, MyXSlices, MyYSlices, nz, id, dTN, dTsigma, NeighborRank_North,
               NeighborRank_South, NeighborRank_East, NeighborRank_West, NeighborRank_NorthEast, NeighborRank_NorthWest,
               NeighborRank_SouthEast, NeighborRank_SouthWest, NucleiLocation_H, NucleationTimes_H, CellType_H,
               GrainID_H, CritTimeStep_H, UndercoolingChange_H);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Nucleation initialized" << std::endl;

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
    ViewI NucleiLocation_G = Kokkos::create_mirror_view_and_copy(memory_space(), NucleiLocation_H);
    ViewI NucleationTimes_G = Kokkos::create_mirror_view_and_copy(memory_space(), NucleationTimes_H);
    ViewI NeighborX_G = Kokkos::create_mirror_view_and_copy(memory_space(), NeighborX_H);
    ViewI NeighborY_G = Kokkos::create_mirror_view_and_copy(memory_space(), NeighborY_H);
    ViewI NeighborZ_G = Kokkos::create_mirror_view_and_copy(memory_space(), NeighborZ_H);
    ViewI2D ItList_G = Kokkos::create_mirror_view_and_copy(memory_space(), ItList_H);
    ViewI GrainOrientation_G = Kokkos::create_mirror_view_and_copy(memory_space(), GrainOrientation_H);
    ViewF GrainUnitVector_G = Kokkos::create_mirror_view_and_copy(memory_space(), GrainUnitVector_H);

    // Steering Vector
    ViewI SteeringVector(Kokkos::ViewAllocateWithoutInitializing("SteeringVector"), LocalActiveDomainSize);
    ViewI_H numSteer_H(Kokkos::ViewAllocateWithoutInitializing("SteeringVectorSize"), 1);
    numSteer_H(0) = 0;
    ViewI numSteer_G = Kokkos::create_mirror_view_and_copy(memory_space(), numSteer_H);

    if (np > 1) {
        // Ghost nodes for initial microstructure state
        GhostNodesInit(id, np, DecompositionStrategy, NeighborRank_North, NeighborRank_South, NeighborRank_East,
                       NeighborRank_West, NeighborRank_NorthEast, NeighborRank_NorthWest, NeighborRank_SouthEast,
                       NeighborRank_SouthWest, MyXSlices, MyYSlices, MyXOffset, MyYOffset, ZBound_Low, nzActive,
                       LocalActiveDomainSize, NGrainOrientations, NeighborX_G, NeighborY_G, NeighborZ_G,
                       GrainUnitVector_G, GrainOrientation_G, GrainID_G, CellType_G, DOCenter_G, DiagonalLength_G,
                       CritDiagonalLength_G);
    }

    // If specified, print initial values in some views for debugging purposes
    double InitTime = MPI_Wtime() - StartTime;
    if (id == 0)
        std::cout << "Data initialized: Time spent: " << InitTime << " s" << std::endl;
    if (PrintDebug) {
        PrintExaCAData(id, np, nx, ny, nz, MyXSlices, MyYSlices, ProcessorsInXDirection, ProcessorsInYDirection,
                       GrainID_H, GrainOrientation_H, CritTimeStep_H, GrainUnitVector_H, LayerID_H, CellType_H,
                       UndercoolingChange_H, UndercoolingCurrent_H, OutputFile, DecompositionStrategy,
                       NGrainOrientations, Melted, PathToOutput, PrintDebug, false, false);
        MPI_Barrier(MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Initialization data file(s) printed" << std::endl;
    }
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
            Nucleation(MyXSlices, MyYSlices, MyXOffset, MyYOffset, cycle, nn, CellType_G, NucleiLocation_G,
                       NucleationTimes_G, GrainID_G, GrainOrientation_G, DOCenter_G, NeighborX_G, NeighborY_G,
                       NeighborZ_G, GrainUnitVector_G, CritDiagonalLength_G, DiagonalLength_G, NGrainOrientations,
                       PossibleNuclei_ThisRank, ZBound_Low, layernumber, LayerID_G);
            NuclTime += MPI_Wtime() - StartNuclTime;

            // Update cells on GPU - new active cells, solidification of old active cells
            StartCaptureTime = MPI_Wtime();
            CellCapture(np, cycle, DecompositionStrategy, LocalActiveDomainSize, LocalDomainSize, MyXSlices, MyYSlices,
                        AConst, BConst, CConst, DConst, MyXOffset, MyYOffset, ItList_G, NeighborX_G, NeighborY_G,
                        NeighborZ_G, CritTimeStep_G, UndercoolingCurrent_G, UndercoolingChange_G, GrainUnitVector_G,
                        CritDiagonalLength_G, DiagonalLength_G, GrainOrientation_G, CellType_G, DOCenter_G, GrainID_G,
                        NGrainOrientations, BufferWestSend, BufferEastSend, BufferNorthSend, BufferSouthSend,
                        BufferNorthEastSend, BufferNorthWestSend, BufferSouthEastSend, BufferSouthWestSend, BufSizeX,
                        BufSizeY, ZBound_Low, nzActive, nz, layernumber, LayerID_G, SteeringVector, numSteer_G,
                        numSteer_H);
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
                IntermediateOutputAndCheck(id, cycle, MyXSlices, MyYSlices, LocalDomainSize, LocalActiveDomainSize, nn,
                                           XSwitch, CellType_G, CritTimeStep_G, SimulationType, FinishTimeStep,
                                           layernumber, NumberOfLayers, ZBound_Low, LayerID_G);
            }

        } while (XSwitch == 0);

        if (layernumber != NumberOfLayers - 1) {
            // Determine new active cell domain size and offset from bottom of global domain
            int ZShift;
            DomainShiftAndResize(id, MyXSlices, MyYSlices, ZShift, ZBound_Low, ZBound_High, nzActive, LocalDomainSize,
                                 LocalActiveDomainSize, BufSizeZ, LayerHeight, CellType_G, layernumber, LayerID_G);

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
    if ((PrintMisorientation) || (PrintFullOutput)) {
        if (id == 0)
            std::cout << "Collecting data on rank 0 and printing to files" << std::endl;
        PrintExaCAData(id, np, nx, ny, nz, MyXSlices, MyYSlices, ProcessorsInXDirection, ProcessorsInYDirection,
                       GrainID_H, GrainOrientation_H, CritTimeStep_H, GrainUnitVector_H, LayerID_H, CellType_H,
                       UndercoolingChange_H, UndercoolingCurrent_H, OutputFile, DecompositionStrategy,
                       NGrainOrientations, Melted, PathToOutput, 0, PrintMisorientation, PrintFullOutput);
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

    PrintExaCALog(id, np, InputFile, SimulationType, DecompositionStrategy, MyXSlices, MyYSlices, MyXOffset, MyYOffset,
                  AConst, BConst, CConst, DConst, FreezingRange, deltax, NMax, dTN, dTsigma, tempfile,
                  TempFilesInSeries, HT_deltax, RemeltingYN, deltat, NumberOfLayers, LayerHeight, SubstrateFileName,
                  SubstrateGrainSpacing, UseSubstrateFile, G, R, nx, ny, nz, FractSurfaceSitesActive, PathToOutput,
                  NSpotsX, NSpotsY, SpotOffset, SpotRadius, OutputFile, InitTime, RunTime, OutTime, cycle, InitMaxTime,
                  InitMinTime, NuclMaxTime, NuclMinTime, CaptureMaxTime, CaptureMinTime, GhostMaxTime, GhostMinTime,
                  OutMaxTime, OutMinTime);
}
