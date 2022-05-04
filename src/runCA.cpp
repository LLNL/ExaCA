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
    double StartInitTime = MPI_Wtime();

    int nx, ny, nz, DecompositionStrategy, NumberOfLayers, LayerHeight, TempFilesInSeries;
    int NSpotsX, NSpotsY, SpotOffset, SpotRadius, HTtoCAratio;
    unsigned int NumberOfTemperatureDataPoints = 0; // Initialized to 0 - updated if/when temperature files are read
    int PrintDebug, TimeSeriesInc;
    bool PrintMisorientation, PrintFinalUndercoolingVals, PrintFullOutput, RemeltingYN, UseSubstrateFile,
        PrintTimeSeries, PrintIdleTimeSeriesFrames, PrintDefaultRVE;
    float SubstrateGrainSpacing;
    double HT_deltax, deltax, deltat, FractSurfaceSitesActive, G, R, AConst, BConst, CConst, DConst, FreezingRange,
        NMax, dTN, dTsigma, RNGSeed;
    std::string SubstrateFileName, temppath, tempfile, SimulationType, OutputFile, GrainOrientationFile, PathToOutput;
    std::vector<std::string> temp_paths;

    // Read input data
    InputReadFromFile(id, InputFile, SimulationType, DecompositionStrategy, AConst, BConst, CConst, DConst,
                      FreezingRange, deltax, NMax, dTN, dTsigma, OutputFile, GrainOrientationFile, temppath, tempfile,
                      TempFilesInSeries, temp_paths, HT_deltax, RemeltingYN, deltat, NumberOfLayers, LayerHeight,
                      SubstrateFileName, SubstrateGrainSpacing, UseSubstrateFile, G, R, nx, ny, nz,
                      FractSurfaceSitesActive, PathToOutput, PrintDebug, PrintMisorientation,
                      PrintFinalUndercoolingVals, PrintFullOutput, NSpotsX, NSpotsY, SpotOffset, SpotRadius,
                      PrintTimeSeries, TimeSeriesInc, PrintIdleTimeSeriesFrames, PrintDefaultRVE, RNGSeed);

    // Grid decomposition
    int ProcessorsInXDirection, ProcessorsInYDirection;
    // Variables characterizing local processor grids relative to global domain
    int MyXSlices, MyXOffset, MyYSlices, MyYOffset;
    long int LocalDomainSize;
    // Variables characterizing process IDs of neighboring MPI ranks on the grid
    // Positive X/Negative X directions are West/East, Positive Y/NegativeY directions are North/South
    int NeighborRank_North, NeighborRank_South, NeighborRank_East, NeighborRank_West, NeighborRank_NorthEast,
        NeighborRank_NorthWest, NeighborRank_SouthEast, NeighborRank_SouthWest;
    // Variables denoting whether or not each MPI rank's grid is at a global domain boundary
    bool AtNorthBoundary, AtSouthBoundary, AtEastBoundary, AtWestBoundary;
    // Neighbor lists for cells
    ViewI_H NeighborX_H(Kokkos::ViewAllocateWithoutInitializing("NeighborX"), 26);
    ViewI_H NeighborY_H(Kokkos::ViewAllocateWithoutInitializing("NeighborY"), 26);
    ViewI_H NeighborZ_H(Kokkos::ViewAllocateWithoutInitializing("NeighborZ"), 26);
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

    // Intialize neighbor list structures (NeighborX, NeighborY, NeighborZ)
    NeighborListInit(NeighborX_H, NeighborY_H, NeighborZ_H);

    // Obtain the physical XYZ bounds of the domain, using either domain size from the input file, or reading
    // temperature data files and parsing the coordinates
    FindXYZBounds(SimulationType, id, deltax, nx, ny, nz, temp_paths, XMin, XMax, YMin, YMax, ZMin, ZMax, LayerHeight,
                  NumberOfLayers, TempFilesInSeries, ZMinLayer, ZMaxLayer);

    // Decompose the domain into subdomains on each MPI rank: Each subdomain contains "MyXSlices" cells in X, and
    // "MyYSlices" in Y. Each subdomain is offset from the full domain origin by "MyXOffset" cells in X, and "MyYOffset"
    // cells in Y
    DomainDecomposition(DecompositionStrategy, id, np, MyXSlices, MyYSlices, MyXOffset, MyYOffset, NeighborRank_North,
                        NeighborRank_South, NeighborRank_East, NeighborRank_West, NeighborRank_NorthEast,
                        NeighborRank_NorthWest, NeighborRank_SouthEast, NeighborRank_SouthWest, nx, ny, nz,
                        ProcessorsInXDirection, ProcessorsInYDirection, LocalDomainSize, AtNorthBoundary,
                        AtSouthBoundary, AtEastBoundary, AtWestBoundary);

    // Read in temperature data from files, stored in "RawData", with the appropriate MPI ranks storing the appropriate
    // data
    if (SimulationType == "R")
        ReadTemperatureData(id, deltax, HT_deltax, HTtoCAratio, MyXSlices, MyYSlices, MyXOffset, MyYOffset, XMin, YMin,
                            temp_paths, NumberOfLayers, TempFilesInSeries, NumberOfTemperatureDataPoints, RawData,
                            FirstValue, LastValue);

    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Mesh initialized and (if being used), temperature data read" << std::endl;

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
        TempInit_Reduced(id, MyXSlices, MyYSlices, MyXOffset, MyYOffset, deltax, HTtoCAratio, deltat, nz,
                         CritTimeStep_H, UndercoolingChange_H, UndercoolingCurrent_H, XMin, YMin, ZMin, Melted,
                         ZMinLayer, ZMaxLayer, LayerHeight, NumberOfLayers, nzActive, ZBound_Low, ZBound_High,
                         FinishTimeStep, FreezingRange, LayerID_H, FirstValue, LastValue, RawData);
    }
    else if (SimulationType == "S") {
        // spot melt array test problem
        TempInit_SpotMelt(G, R, SimulationType, id, MyXSlices, MyYSlices, MyXOffset, MyYOffset, deltax, ZMin, ZMinLayer,
                          ZMaxLayer, deltat, nz, CritTimeStep_H, UndercoolingChange_H, UndercoolingCurrent_H, Melted,
                          LayerHeight, NumberOfLayers, nzActive, ZBound_Low, ZBound_High, FreezingRange, LayerID_H,
                          NSpotsX, NSpotsY, SpotRadius, SpotOffset);
    }
    else if (SimulationType == "C") {
        // directional/constrained solidification test problem
        TempInit_DirSolidification(G, R, id, MyXSlices, MyYSlices, deltax, deltat, nz, CritTimeStep_H,
                                   UndercoolingChange_H, UndercoolingCurrent_H, Melted, nzActive, ZBound_Low,
                                   ZBound_High, LayerID_H);
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

    // Initialize grain orientations
    OrientationInit(id, NGrainOrientations, GrainUnitVector_H, GrainOrientationFile);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Done with orientation initialization " << std::endl;

    // The remainder of the initialization is performed on the device
    // Copy previously initialized host view data to the device
    using memory_space = Kokkos::DefaultExecutionSpace::memory_space;
    ViewI CritTimeStep_G = Kokkos::create_mirror_view_and_copy(memory_space(), CritTimeStep_H);
    ViewI LayerID_G = Kokkos::create_mirror_view_and_copy(memory_space(), LayerID_H);
    ViewF UndercoolingChange_G = Kokkos::create_mirror_view_and_copy(memory_space(), UndercoolingChange_H);
    ViewF UndercoolingCurrent_G = Kokkos::create_mirror_view_and_copy(memory_space(), UndercoolingCurrent_H);
    ViewI NeighborX_G = Kokkos::create_mirror_view_and_copy(memory_space(), NeighborX_H);
    ViewI NeighborY_G = Kokkos::create_mirror_view_and_copy(memory_space(), NeighborY_H);
    ViewI NeighborZ_G = Kokkos::create_mirror_view_and_copy(memory_space(), NeighborZ_H);
    ViewF GrainUnitVector_G = Kokkos::create_mirror_view_and_copy(memory_space(), GrainUnitVector_H);

    // Allocate device views: initialize GrainID to 0 for all cells (unassigned), assign CellType values later
    ViewI GrainID_G("GrainID_G", LocalDomainSize);
    ViewI CellType_G(Kokkos::ViewAllocateWithoutInitializing("CellType_G"), LocalDomainSize);
    // Variables characterizing the active cell region within each rank's grid
    // DiagonalLengths should be initialized to 0, DOCenter/CritDiagonalLength do not need initialization
    ViewF DiagonalLength_G("DiagonalLength_G", LocalActiveDomainSize);
    ViewF DOCenter_G(Kokkos::ViewAllocateWithoutInitializing("DOCenter_G"), 3 * LocalActiveDomainSize);
    ViewF CritDiagonalLength_G(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength_G"),
                               26 * LocalActiveDomainSize);

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
    // Send/recv buffers for ghost node data should be initialized with zeros
    Buffer2D BufferSouthSend("BufferSouthSend", BufSizeX * BufSizeZ, 5);
    Buffer2D BufferNorthSend("BufferNorthSend", BufSizeX * BufSizeZ, 5);
    Buffer2D BufferEastSend("BufferEastSend", BufSizeY * BufSizeZ, 5);
    Buffer2D BufferWestSend("BufferWestSend", BufSizeY * BufSizeZ, 5);
    Buffer2D BufferNorthEastSend("BufferNorthEastSend", BufSizeZ, 5);
    Buffer2D BufferNorthWestSend("BufferNorthWestSend", BufSizeZ, 5);
    Buffer2D BufferSouthEastSend("BufferSouthEastSend", BufSizeZ, 5);
    Buffer2D BufferSouthWestSend("BufferSouthWestSend", BufSizeZ, 5);
    Buffer2D BufferSouthRecv("BufferSouthRecv", BufSizeX * BufSizeZ, 5);
    Buffer2D BufferNorthRecv("BufferNorthRecv", BufSizeX * BufSizeZ, 5);
    Buffer2D BufferEastRecv("BufferEastRecv", BufSizeY * BufSizeZ, 5);
    Buffer2D BufferWestRecv("BufferWestRecv", BufSizeY * BufSizeZ, 5);
    Buffer2D BufferNorthEastRecv("BufferNorthEastRecv", BufSizeZ, 5);
    Buffer2D BufferNorthWestRecv("BufferNorthWestRecv", BufSizeZ, 5);
    Buffer2D BufferSouthEastRecv("BufferSouthEastRecv", BufSizeZ, 5);
    Buffer2D BufferSouthWestRecv("BufferSouthWestRecv", BufSizeZ, 5);

    // Initialize the grain structure and cell types - for either a constrained solidification problem, using a
    // substrate from a file, or generating a substrate using the existing CA algorithm
    int NextLayer_FirstEpitaxialGrainID;
    if (SimulationType == "C") {
        SubstrateInit_ConstrainedGrowth(
            id, FractSurfaceSitesActive, MyXSlices, MyYSlices, nx, ny, MyXOffset, MyYOffset, NeighborX_G, NeighborY_G,
            NeighborZ_G, GrainUnitVector_G, NGrainOrientations, CellType_G, GrainID_G, DiagonalLength_G, DOCenter_G,
            CritDiagonalLength_G, RNGSeed, np, DecompositionStrategy, BufferWestSend, BufferEastSend, BufferNorthSend,
            BufferSouthSend, BufferNorthEastSend, BufferNorthWestSend, BufferSouthEastSend, BufferSouthWestSend,
            BufSizeX, BufSizeY, AtNorthBoundary, AtSouthBoundary, AtEastBoundary, AtWestBoundary);
    }
    else {
        if (UseSubstrateFile)
            SubstrateInit_FromFile(SubstrateFileName, nz, MyXSlices, MyYSlices, MyXOffset, MyYOffset, id, GrainID_G);
        else
            BaseplateInit_FromGrainSpacing(SubstrateGrainSpacing, nx, ny, ZMinLayer, ZMaxLayer, MyXSlices, MyYSlices,
                                           MyXOffset, MyYOffset, id, deltax, GrainID_G, RNGSeed,
                                           NextLayer_FirstEpitaxialGrainID);
        // Separate routine for active cell data structure init for problems other than constrained solidification
        CellTypeInit(0, id, np, DecompositionStrategy, MyXSlices, MyYSlices, MyXOffset, MyYOffset, ZBound_Low, nz,
                     LocalActiveDomainSize, LocalDomainSize, CellType_G, CritTimeStep_G, NeighborX_G, NeighborY_G,
                     NeighborZ_G, NGrainOrientations, GrainUnitVector_G, DiagonalLength_G, GrainID_G,
                     CritDiagonalLength_G, DOCenter_G, LayerID_G, BufferWestSend, BufferEastSend, BufferNorthSend,
                     BufferSouthSend, BufferNorthEastSend, BufferNorthWestSend, BufferSouthEastSend,
                     BufferSouthWestSend, BufSizeX, BufSizeY, AtNorthBoundary, AtSouthBoundary, AtEastBoundary,
                     AtWestBoundary);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Grain struct initialized" << std::endl;

    int PossibleNuclei_ThisRankThisLayer, Nuclei_WholeDomain, NucleationCounter;
    // Without knowing PossibleNuclei_ThisRankThisLayer yet, initialize nucleation data structures to estimated sizes,
    // resize inside of NucleiInit when the number of nuclei per rank is known
    int EstimatedNuclei_ThisRankThisLayer = NMax * pow(deltax, 3) * LocalActiveDomainSize;
    ViewI_H NucleationTimes_Host(Kokkos::ViewAllocateWithoutInitializing("NucleationTimes_Host"),
                                 EstimatedNuclei_ThisRankThisLayer);
    ViewI NucleiLocation(Kokkos::ViewAllocateWithoutInitializing("NucleiLocation"), EstimatedNuclei_ThisRankThisLayer);
    ViewI NucleiGrainID(Kokkos::ViewAllocateWithoutInitializing("NucleiGrainID"), EstimatedNuclei_ThisRankThisLayer);

    // Fill in nucleation data structures, and assign nucleation undercooling values to potential nucleation events
    // Potential nucleation grains are only associated with liquid cells in layer 0 - they will be initialized for each
    // successive layer when layer 0 in complete
    NucleiInit(0, RNGSeed, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nx, ny, nzActive, ZBound_Low, id, NMax, dTN,
               dTsigma, deltax, NucleiLocation, NucleationTimes_Host, NucleiGrainID, CellType_G, CritTimeStep_G,
               UndercoolingChange_G, LayerID_G, PossibleNuclei_ThisRankThisLayer, Nuclei_WholeDomain, AtNorthBoundary,
               AtSouthBoundary, AtEastBoundary, AtWestBoundary, NucleationCounter);

    // Normalize solidification parameters
    AConst = AConst * deltat / deltax;
    BConst = BConst * deltat / deltax;
    CConst = CConst * deltat / deltax;
    int cycle;

    // Steering Vector
    ViewI SteeringVector(Kokkos::ViewAllocateWithoutInitializing("SteeringVector"), LocalActiveDomainSize);
    ViewI_H numSteer_H(Kokkos::ViewAllocateWithoutInitializing("SteeringVectorSize"), 1);
    numSteer_H(0) = 0;
    ViewI numSteer_G = Kokkos::create_mirror_view_and_copy(memory_space(), numSteer_H);

    // Update ghost node data for initial state of simulation
    if (np > 1) {
        if (DecompositionStrategy == 1)
            GhostNodes1D(-1, id, NeighborRank_North, NeighborRank_South, MyXSlices, MyYSlices, MyXOffset, MyYOffset,
                         NeighborX_G, NeighborY_G, NeighborZ_G, CellType_G, DOCenter_G, GrainID_G, GrainUnitVector_G,
                         DiagonalLength_G, CritDiagonalLength_G, NGrainOrientations, BufferNorthSend, BufferSouthSend,
                         BufferNorthRecv, BufferSouthRecv, BufSizeX, BufSizeY, BufSizeZ, ZBound_Low);
        else
            GhostNodes2D(-1, id, NeighborRank_North, NeighborRank_South, NeighborRank_East, NeighborRank_West,
                         NeighborRank_NorthEast, NeighborRank_NorthWest, NeighborRank_SouthEast, NeighborRank_SouthWest,
                         MyXSlices, MyYSlices, MyXOffset, MyYOffset, NeighborX_G, NeighborY_G, NeighborZ_G, CellType_G,
                         DOCenter_G, GrainID_G, GrainUnitVector_G, DiagonalLength_G, CritDiagonalLength_G,
                         NGrainOrientations, BufferWestSend, BufferEastSend, BufferNorthSend, BufferSouthSend,
                         BufferNorthEastSend, BufferNorthWestSend, BufferSouthEastSend, BufferSouthWestSend,
                         BufferWestRecv, BufferEastRecv, BufferNorthRecv, BufferSouthRecv, BufferNorthEastRecv,
                         BufferNorthWestRecv, BufferSouthEastRecv, BufferSouthWestRecv, BufSizeX, BufSizeY, BufSizeZ,
                         ZBound_Low);
    }

    // If specified, print initial values in some views for debugging purposes
    double InitTime = MPI_Wtime() - StartInitTime;
    if (id == 0)
        std::cout << "Data initialized: Time spent: " << InitTime << " s" << std::endl;
    if (PrintDebug) {
        // Host mirrors of CellType and GrainID are not maintained - pass device views and perform copy inside of
        // subroutine
        PrintExaCAData(id, -1, np, nx, ny, nz, MyXSlices, MyYSlices, MyXOffset, MyYOffset, ProcessorsInXDirection,
                       ProcessorsInYDirection, GrainID_G, CritTimeStep_H, GrainUnitVector_H, LayerID_H, CellType_G,
                       UndercoolingChange_H, UndercoolingCurrent_H, OutputFile, DecompositionStrategy,
                       NGrainOrientations, Melted, PathToOutput, PrintDebug, false, false, false, false, false, 0,
                       ZBound_Low, nzActive, deltax, XMin, YMin, ZMin, NumberOfLayers);
        MPI_Barrier(MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Initialization data file(s) printed" << std::endl;
    }
    cycle = 0;
    int IntermediateFileCounter = 0;
    double StartRunTime = MPI_Wtime();
    for (int layernumber = 0; layernumber < NumberOfLayers; layernumber++) {

        int SuccessfulNucEvents_ThisRank = 0;
        int XSwitch = 0;
        double LayerTime1 = MPI_Wtime();

        // Loop continues until all liquid cells claimed by solid grains
        do {
            // Start of time step - check and see if intermediate system output is to be printed to files
            if ((PrintTimeSeries) && (cycle % TimeSeriesInc == 0)) {
                // Print current state of ExaCA simulation (up to and including the current layer's data)
                // Host mirrors of CellType and GrainID are not maintained - pass device views and perform copy inside
                // of subroutine
                PrintExaCAData(id, layernumber, np, nx, ny, nz, MyXSlices, MyYSlices, MyXOffset, MyYOffset,
                               ProcessorsInXDirection, ProcessorsInYDirection, GrainID_G, CritTimeStep_H,
                               GrainUnitVector_H, LayerID_H, CellType_G, UndercoolingChange_H, UndercoolingCurrent_H,
                               OutputFile, DecompositionStrategy, NGrainOrientations, Melted, PathToOutput, 0, false,
                               false, false, true, false, IntermediateFileCounter, ZBound_Low, nzActive, deltax, XMin,
                               YMin, ZMin, NumberOfLayers);
                IntermediateFileCounter++;
            }
            cycle++;

            // Update cells on GPU - undercooling and diagonal length updates, nucleation
            StartNuclTime = MPI_Wtime();
            Nucleation(cycle, SuccessfulNucEvents_ThisRank, NucleationCounter, PossibleNuclei_ThisRankThisLayer,
                       NucleationTimes_Host, NucleiLocation, NucleiGrainID, CellType_G, GrainID_G, ZBound_Low,
                       MyXSlices, MyYSlices, SteeringVector, numSteer_G);
            NuclTime += MPI_Wtime() - StartNuclTime;

            // Update cells on GPU - new active cells, solidification of old active cells
            StartCaptureTime = MPI_Wtime();
            CellCapture(np, cycle, DecompositionStrategy, LocalActiveDomainSize, LocalDomainSize, MyXSlices, MyYSlices,
                        AConst, BConst, CConst, DConst, MyXOffset, MyYOffset, NeighborX_G, NeighborY_G, NeighborZ_G,
                        CritTimeStep_G, UndercoolingCurrent_G, UndercoolingChange_G, GrainUnitVector_G,
                        CritDiagonalLength_G, DiagonalLength_G, CellType_G, DOCenter_G, GrainID_G, NGrainOrientations,
                        BufferWestSend, BufferEastSend, BufferNorthSend, BufferSouthSend, BufferNorthEastSend,
                        BufferNorthWestSend, BufferSouthEastSend, BufferSouthWestSend, BufSizeX, BufSizeY, ZBound_Low,
                        nzActive, nz, layernumber, LayerID_G, SteeringVector, numSteer_G, numSteer_H, AtNorthBoundary,
                        AtSouthBoundary, AtEastBoundary, AtWestBoundary);
            CaptureTime += MPI_Wtime() - StartCaptureTime;

            if (np > 1) {
                // Update ghost nodes
                StartGhostTime = MPI_Wtime();
                if (DecompositionStrategy == 1)
                    GhostNodes1D(cycle, id, NeighborRank_North, NeighborRank_South, MyXSlices, MyYSlices, MyXOffset,
                                 MyYOffset, NeighborX_G, NeighborY_G, NeighborZ_G, CellType_G, DOCenter_G, GrainID_G,
                                 GrainUnitVector_G, DiagonalLength_G, CritDiagonalLength_G, NGrainOrientations,
                                 BufferNorthSend, BufferSouthSend, BufferNorthRecv, BufferSouthRecv, BufSizeX, BufSizeY,
                                 BufSizeZ, ZBound_Low);
                else
                    GhostNodes2D(cycle, id, NeighborRank_North, NeighborRank_South, NeighborRank_East,
                                 NeighborRank_West, NeighborRank_NorthEast, NeighborRank_NorthWest,
                                 NeighborRank_SouthEast, NeighborRank_SouthWest, MyXSlices, MyYSlices, MyXOffset,
                                 MyYOffset, NeighborX_G, NeighborY_G, NeighborZ_G, CellType_G, DOCenter_G, GrainID_G,
                                 GrainUnitVector_G, DiagonalLength_G, CritDiagonalLength_G, NGrainOrientations,
                                 BufferWestSend, BufferEastSend, BufferNorthSend, BufferSouthSend, BufferNorthEastSend,
                                 BufferNorthWestSend, BufferSouthEastSend, BufferSouthWestSend, BufferWestRecv,
                                 BufferEastRecv, BufferNorthRecv, BufferSouthRecv, BufferNorthEastRecv,
                                 BufferNorthWestRecv, BufferSouthEastRecv, BufferSouthWestRecv, BufSizeX, BufSizeY,
                                 BufSizeZ, ZBound_Low);
                GhostTime += MPI_Wtime() - StartGhostTime;
            }

            if (cycle % 1000 == 0) {
                IntermediateOutputAndCheck(
                    id, np, cycle, MyXSlices, MyYSlices, MyXOffset, MyYOffset, LocalDomainSize, LocalActiveDomainSize,
                    nx, ny, nz, nzActive, deltax, XMin, YMin, ZMin, DecompositionStrategy, ProcessorsInXDirection,
                    ProcessorsInYDirection, SuccessfulNucEvents_ThisRank, XSwitch, CellType_G, CritTimeStep_G,
                    CritTimeStep_H, GrainID_G, SimulationType, FinishTimeStep, layernumber, NumberOfLayers, ZBound_Low,
                    NGrainOrientations, Melted, LayerID_G, LayerID_H, GrainUnitVector_H, UndercoolingChange_H,
                    UndercoolingCurrent_H, PathToOutput, OutputFile, PrintIdleTimeSeriesFrames, TimeSeriesInc,
                    IntermediateFileCounter, NumberOfLayers);
            }

        } while (XSwitch == 0);
        if (layernumber != NumberOfLayers - 1) {
            // Reset intermediate file counter to zero if printing video files
            if (PrintTimeSeries)
                IntermediateFileCounter = 0;

            // Determine new active cell domain size and offset from bottom of global domain
            int ZShift;
            DomainShiftAndResize(id, MyXSlices, MyYSlices, ZShift, ZBound_Low, ZBound_High, nzActive, LocalDomainSize,
                                 LocalActiveDomainSize, BufSizeZ, LayerHeight, CellType_G, layernumber, LayerID_G);

            // Resize steering vector as LocalActiveDomainSize may have changed
            Kokkos::resize(SteeringVector, LocalActiveDomainSize);

            // Resize active cell data structures - host and device
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

            // Reset halo region views on device to 0 for the next layer
            // Also reset active region views on host to 0 for the next layer (setting up active cell data structures
            // and nuclei locations will be performed on the host, then copied back to the device)
            ZeroResetViews(DiagonalLength_G, CritDiagonalLength_G, DOCenter_G, DecompositionStrategy, BufferWestSend,
                           BufferEastSend, BufferNorthSend, BufferSouthSend, BufferNorthEastSend, BufferNorthWestSend,
                           BufferSouthEastSend, BufferSouthWestSend, BufferWestRecv, BufferEastRecv, BufferNorthRecv,
                           BufferSouthRecv, BufferNorthEastRecv, BufferNorthWestRecv, BufferSouthEastRecv,
                           BufferSouthWestRecv);

            if (id == 0)
                std::cout << "New layer setup, GN dimensions are " << BufSizeX << " " << BufSizeY << " " << BufSizeZ
                          << std::endl;

            // If the baseplate was initialized from a substrate grain spacing, initialize powder layer grain structure
            // for the next layer "layernumber + 1" Otherwise, the entire substrate (baseplate + powder) was read from a
            // file, and the powder layers have already been initialized
            if (!(UseSubstrateFile))
                PowderInit(layernumber + 1, nx, ny, LayerHeight, ZMaxLayer, ZMin, deltax, MyXSlices, MyYSlices,
                           MyXOffset, MyYOffset, id, GrainID_G, RNGSeed, NextLayer_FirstEpitaxialGrainID);

            // Initialize active cell data structures and nuclei locations for the next layer "layernumber + 1"
            CellTypeInit(layernumber + 1, id, np, DecompositionStrategy, MyXSlices, MyYSlices, MyXOffset, MyYOffset,
                         ZBound_Low, nz, LocalActiveDomainSize, LocalDomainSize, CellType_G, CritTimeStep_G,
                         NeighborX_G, NeighborY_G, NeighborZ_G, NGrainOrientations, GrainUnitVector_G, DiagonalLength_G,
                         GrainID_G, CritDiagonalLength_G, DOCenter_G, LayerID_G, BufferWestSend, BufferEastSend,
                         BufferNorthSend, BufferSouthSend, BufferNorthEastSend, BufferNorthWestSend,
                         BufferSouthEastSend, BufferSouthWestSend, BufSizeX, BufSizeY, AtNorthBoundary, AtSouthBoundary,
                         AtEastBoundary, AtWestBoundary);

            // Initialize potential nucleation event data for next layer "layernumber + 1"
            // Views containing nucleation data will be resized to the possible number of nuclei on a given MPI rank for
            // the next layer
            NucleiInit(layernumber + 1, RNGSeed, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nx, ny, nzActive,
                       ZBound_Low, id, NMax, dTN, dTsigma, deltax, NucleiLocation, NucleationTimes_Host, NucleiGrainID,
                       CellType_G, CritTimeStep_G, UndercoolingChange_G, LayerID_G, PossibleNuclei_ThisRankThisLayer,
                       Nuclei_WholeDomain, AtNorthBoundary, AtSouthBoundary, AtEastBoundary, AtWestBoundary,
                       NucleationCounter);

            // Update ghost nodes for grain locations and attributes
            MPI_Barrier(MPI_COMM_WORLD);
            if (np > 1) {
                if (DecompositionStrategy == 1)
                    GhostNodes1D(-1, id, NeighborRank_North, NeighborRank_South, MyXSlices, MyYSlices, MyXOffset,
                                 MyYOffset, NeighborX_G, NeighborY_G, NeighborZ_G, CellType_G, DOCenter_G, GrainID_G,
                                 GrainUnitVector_G, DiagonalLength_G, CritDiagonalLength_G, NGrainOrientations,
                                 BufferNorthSend, BufferSouthSend, BufferNorthRecv, BufferSouthRecv, BufSizeX, BufSizeY,
                                 BufSizeZ, ZBound_Low);
                else
                    GhostNodes2D(-1, id, NeighborRank_North, NeighborRank_South, NeighborRank_East, NeighborRank_West,
                                 NeighborRank_NorthEast, NeighborRank_NorthWest, NeighborRank_SouthEast,
                                 NeighborRank_SouthWest, MyXSlices, MyYSlices, MyXOffset, MyYOffset, NeighborX_G,
                                 NeighborY_G, NeighborZ_G, CellType_G, DOCenter_G, GrainID_G, GrainUnitVector_G,
                                 DiagonalLength_G, CritDiagonalLength_G, NGrainOrientations, BufferWestSend,
                                 BufferEastSend, BufferNorthSend, BufferSouthSend, BufferNorthEastSend,
                                 BufferNorthWestSend, BufferSouthEastSend, BufferSouthWestSend, BufferWestRecv,
                                 BufferEastRecv, BufferNorthRecv, BufferSouthRecv, BufferNorthEastRecv,
                                 BufferNorthWestRecv, BufferSouthEastRecv, BufferSouthWestRecv, BufSizeX, BufSizeY,
                                 BufSizeZ, ZBound_Low);
            }
            if (id == 0)
                std::cout << "New layer ghost nodes initialized" << std::endl;
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

    double RunTime = MPI_Wtime() - StartRunTime;
    double StartOutTime = MPI_Wtime();

    if (PrintFinalUndercoolingVals) {
        // Copy final undercooling values back to the CPU for potential printing to file
        Kokkos::deep_copy(UndercoolingCurrent_H, UndercoolingCurrent_G);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (((PrintMisorientation) || (PrintFinalUndercoolingVals) || (PrintFullOutput)) || (PrintDefaultRVE)) {
        if (id == 0)
            std::cout << "Collecting data on rank 0 and printing to files" << std::endl;
        // Host mirrors of CellType and GrainID are not maintained - pass device views and perform copy inside of
        // subroutine
        PrintExaCAData(id, NumberOfLayers - 1, np, nx, ny, nz, MyXSlices, MyYSlices, MyXOffset, MyYOffset,
                       ProcessorsInXDirection, ProcessorsInYDirection, GrainID_G, CritTimeStep_H, GrainUnitVector_H,
                       LayerID_H, CellType_G, UndercoolingChange_H, UndercoolingCurrent_H, OutputFile,
                       DecompositionStrategy, NGrainOrientations, Melted, PathToOutput, 0, PrintMisorientation,
                       PrintFinalUndercoolingVals, PrintFullOutput, false, PrintDefaultRVE, 0, ZBound_Low, nzActive,
                       deltax, XMin, YMin, ZMin, NumberOfLayers);
    }
    else {
        if (id == 0)
            std::cout << "No output files to be printed, exiting program" << std::endl;
    }

    double OutTime = MPI_Wtime() - StartOutTime;
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
