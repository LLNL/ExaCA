// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
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
    double NuclTime = 0.0, CreateSVTime = 0.0, CaptureTime = 0.0, GhostTime = 0.0;
    double StartNuclTime, StartCreateSVTime, StartCaptureTime, StartGhostTime;
    double StartInitTime = MPI_Wtime();

    int nx, ny, nz, DecompositionStrategy, NumberOfLayers, LayerHeight, TempFilesInSeries;
    int NSpotsX, NSpotsY, SpotOffset, SpotRadius, HTtoCAratio, RVESize;
    unsigned int NumberOfTemperatureDataPoints = 0; // Initialized to 0 - updated if/when temperature files are read
    int PrintDebug, TimeSeriesInc;
    bool PrintMisorientation, PrintFinalUndercoolingVals, PrintFullOutput, RemeltingYN, UseSubstrateFile,
        PrintTimeSeries, PrintIdleTimeSeriesFrames, PrintDefaultRVE, BaseplateThroughPowder;
    float SubstrateGrainSpacing;
    double HT_deltax, deltax, deltat, FractSurfaceSitesActive, G, R, AConst, BConst, CConst, DConst, FreezingRange,
        NMax, dTN, dTsigma, RNGSeed, PowderDensity;
    std::string SubstrateFileName, SimulationType, OutputFile, GrainOrientationFile, PathToOutput;
    std::vector<std::string> temp_paths;

    // Read input data
    InputReadFromFile(id, InputFile, SimulationType, DecompositionStrategy, AConst, BConst, CConst, DConst,
                      FreezingRange, deltax, NMax, dTN, dTsigma, OutputFile, GrainOrientationFile, TempFilesInSeries,
                      temp_paths, HT_deltax, RemeltingYN, deltat, NumberOfLayers, LayerHeight, SubstrateFileName,
                      SubstrateGrainSpacing, UseSubstrateFile, G, R, nx, ny, nz, FractSurfaceSitesActive, PathToOutput,
                      PrintDebug, PrintMisorientation, PrintFinalUndercoolingVals, PrintFullOutput, NSpotsX, NSpotsY,
                      SpotOffset, SpotRadius, PrintTimeSeries, TimeSeriesInc, PrintIdleTimeSeriesFrames,
                      PrintDefaultRVE, RNGSeed, BaseplateThroughPowder, PowderDensity, RVESize);

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
    // View initialization is performed on host views, but copied to the device views within subroutines
    // Neighbor lists for cells
    NList NeighborX, NeighborY, NeighborZ;
    float XMin, YMin, ZMin, XMax, YMax, ZMax; // OpenFOAM simulation bounds (if using OpenFOAM data)
    float *ZMinLayer = new float[NumberOfLayers];
    float *ZMaxLayer = new float[NumberOfLayers];
    int *FinishTimeStep = new int[NumberOfLayers];

    // Data structure for storing raw temperature data from file(s)
    // Store data as double - needed for small time steps to resolve local differences in solidification conditions
    // With no remelting, each data point has 5 values (X, Y, Z coordinates, liquidus time, and cooling rate)
    // With remelting, each data point has 6 values (X, Y, Z coordinates, melting time, liquidus time, and cooling rate)
    // Initial estimate for size
    std::vector<double> RawData(1000000);

    // Contains "NumberOfLayers" values corresponding to the location within "RawData" of the first data element in each
    // temperature file
    int *FirstValue = new int[NumberOfLayers];
    int *LastValue = new int[NumberOfLayers];

    // Intialize neighbor list structures (NeighborX, NeighborY, NeighborZ)
    NeighborListInit(NeighborX, NeighborY, NeighborZ);

    // Obtain the physical XYZ bounds of the domain, using either domain size from the input file, or reading
    // temperature data files and parsing the coordinates
    FindXYZBounds(SimulationType, id, deltax, nx, ny, nz, temp_paths, XMin, XMax, YMin, YMax, ZMin, ZMax, LayerHeight,
                  NumberOfLayers, TempFilesInSeries, ZMinLayer, ZMaxLayer, SpotRadius);

    // Ensure that input powder layer init options are compatible with this domain size, if needed for this problem type
    if ((SimulationType == "R") || (SimulationType == "S"))
        checkPowderOverflow(nx, ny, LayerHeight, NumberOfLayers, BaseplateThroughPowder, PowderDensity);

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
    // Current undercooling is initialized to 0
    ViewI CritTimeStep(Kokkos::ViewAllocateWithoutInitializing("CritTimeStep"), LocalDomainSize);
    ViewI LayerID(Kokkos::ViewAllocateWithoutInitializing("LayerID"), LocalDomainSize);
    ViewF UndercoolingChange(Kokkos::ViewAllocateWithoutInitializing("UndercoolingChange"), LocalDomainSize);
    ViewF UndercoolingCurrent("UndercoolingCurrent", LocalDomainSize);

    // With remelting, temperature fields are also characterized by these variables:
    // Maximum number of times a cell in a given layer undergoes solidification
    ViewI MaxSolidificationEvents(Kokkos::ViewAllocateWithoutInitializing("NumberOfRemeltEvents"), NumberOfLayers);
    // For each cell in the current layer (index 1), and each time solidification happens (index 2), hold the values
    // that will be used for MeltTimeStep, CritTimeStep, and UndercoolingChange (index 3)
    ViewF3D LayerTimeTempHistory(Kokkos::ViewAllocateWithoutInitializing("TimeTempHistory"), 0, 0, 0);
    // The number of times that each CA cell will undergo solidification during this layer
    ViewI NumberOfSolidificationEvents(Kokkos::ViewAllocateWithoutInitializing("NumSEvents"), 0);
    // A counter for the number of times each CA cell has undergone solidification so far this layer
    ViewI SolidificationEventCounter(Kokkos::ViewAllocateWithoutInitializing("SEventCounter"), 0);
    // The next time that each cell will melt during this layer
    ViewI MeltTimeStep(Kokkos::ViewAllocateWithoutInitializing("MeltTimeStep"), 0);

    // Bounds of the current layer: Z coordinates span ZBound_Low-ZBound_High, inclusive
    int ZBound_Low, ZBound_High, nzActive, LocalActiveDomainSize;
    if (RemeltingYN) {
        // For simulations with remelting, ZBound_Low/ZBound_High/nzActive/LocalActiveDomainSize are calculated in
        // advance of the first layer
        ZBound_Low = calcZBound_Low_Remelt(SimulationType, LayerHeight, 0, ZMinLayer, ZMin, deltax);
        ZBound_High = calcZBound_High_Remelt(SimulationType, SpotRadius, LayerHeight, 0, ZMin, deltax, nz, ZMaxLayer);
        nzActive = calcnzActive(ZBound_Low, ZBound_High, id, 0);
        LocalActiveDomainSize =
            calcLocalActiveDomainSize(MyXSlices, MyYSlices, nzActive); // Number of active cells on this MPI rank
        // Initialize the temperature fields:
        // R: input temperature data from files using reduced/sparse data format (with remelting)
        // S: spot melt array test problem (with remelting)
        if (SimulationType == "R")
            TempInit_ReadDataRemelt(0, id, MyXSlices, MyYSlices, nz, LocalActiveDomainSize, LocalDomainSize, MyXOffset,
                                    MyYOffset, deltax, deltat, FreezingRange, LayerTimeTempHistory,
                                    NumberOfSolidificationEvents, MaxSolidificationEvents, MeltTimeStep, CritTimeStep,
                                    UndercoolingChange, UndercoolingCurrent, XMin, YMin, ZMinLayer, LayerHeight,
                                    nzActive, ZBound_Low, FinishTimeStep, LayerID, FirstValue, LastValue, RawData,
                                    SolidificationEventCounter, TempFilesInSeries);
        else if (SimulationType == "S")
            TempInit_SpotRemelt(0, G, R, SimulationType, id, MyXSlices, MyYSlices, MyXOffset, MyYOffset, deltax, deltat,
                                ZBound_Low, nz, LocalActiveDomainSize, LocalDomainSize, CritTimeStep,
                                UndercoolingChange, UndercoolingCurrent, LayerHeight, FreezingRange, LayerID, NSpotsX,
                                NSpotsY, SpotRadius, SpotOffset, LayerTimeTempHistory, NumberOfSolidificationEvents,
                                MeltTimeStep, MaxSolidificationEvents, SolidificationEventCounter);
    }
    else {
        // For simulations without remelting, ZBound_Low/ZBound_High/nzActive/LocalActiveDomainSize are calculated after
        // the multilayer temperature field is initialized, to avoid unneeded inclusion of cells that solidify during
        // the first layer, but remelt during the second layer Initialize the temperature fields R: input temperature
        // data from files using reduced/sparse data format (without remelting) S: spot melt array test problem (without
        // remelting) C: directional/constrained solidification test problem
        if (SimulationType == "R")
            TempInit_ReadDataNoRemelt(id, MyXSlices, MyYSlices, MyXOffset, MyYOffset, deltax, HTtoCAratio, deltat, nz,
                                      LocalDomainSize, CritTimeStep, UndercoolingChange, XMin, YMin, ZMin, ZMinLayer,
                                      ZMaxLayer, LayerHeight, NumberOfLayers, FinishTimeStep, FreezingRange, LayerID,
                                      FirstValue, LastValue, RawData);
        else if (SimulationType == "S")
            TempInit_SpotNoRemelt(G, R, SimulationType, id, MyXSlices, MyYSlices, MyXOffset, MyYOffset, deltax, deltat,
                                  nz, LocalDomainSize, CritTimeStep, UndercoolingChange, LayerHeight, NumberOfLayers,
                                  FreezingRange, LayerID, NSpotsX, NSpotsY, SpotRadius, SpotOffset);
        else if (SimulationType == "C")
            TempInit_DirSolidification(G, R, id, MyXSlices, MyYSlices, deltax, deltat, nz, LocalDomainSize,
                                       CritTimeStep, UndercoolingChange, LayerID);
        // Assign values to the bounds of the current layer now that the multilayer temperature field is initialized
        ZBound_Low = calcZBound_Low_NoRemelt(id, MyXSlices, MyYSlices, LocalDomainSize, 0, LayerID);
        ZBound_High = calcZBound_High_NoRemelt(id, MyXSlices, MyYSlices, LocalDomainSize, 0, LayerID);
        nzActive = calcnzActive(ZBound_Low, ZBound_High, id, 0);
        LocalActiveDomainSize =
            calcLocalActiveDomainSize(MyXSlices, MyYSlices, nzActive); // Number of active cells on this MPI rank
    }
    // Delete temporary data structure for temperature data read if remelting is not performed (otherwise keep it to
    // avoid having to reread temperature files)
    if (!(RemeltingYN))
        RawData.clear();
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Done with temperature field initialization, active domain size is " << nzActive << " out of "
                  << nz << " cells in the Z direction" << std::endl;

    int NGrainOrientations = 0; // Number of grain orientations considered in the simulation
    // No initialize size yet, will be resized in OrientationInit
    ViewF GrainUnitVector(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), 0);

    // Initialize grain orientations
    OrientationInit(id, NGrainOrientations, GrainUnitVector, GrainOrientationFile);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Done with orientation initialization " << std::endl;

    // Allocate device views: initialize GrainID to 0 for all cells (unassigned), assign CellType values later
    ViewI GrainID("GrainID", LocalDomainSize);
    ViewI CellType(Kokkos::ViewAllocateWithoutInitializing("CellType"), LocalDomainSize);
    // Variables characterizing the active cell region within each rank's grid
    // DiagonalLengths should be initialized to 0, DOCenter/CritDiagonalLength do not need initialization
    ViewF DiagonalLength("DiagonalLength", LocalActiveDomainSize);
    ViewF DOCenter(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * LocalActiveDomainSize);
    ViewF CritDiagonalLength(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength"), 26 * LocalActiveDomainSize);

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
            id, FractSurfaceSitesActive, MyXSlices, MyYSlices, nx, ny, MyXOffset, MyYOffset, NeighborX, NeighborY,
            NeighborZ, GrainUnitVector, NGrainOrientations, CellType, GrainID, DiagonalLength, DOCenter,
            CritDiagonalLength, RNGSeed, np, DecompositionStrategy, BufferWestSend, BufferEastSend, BufferNorthSend,
            BufferSouthSend, BufferNorthEastSend, BufferNorthWestSend, BufferSouthEastSend, BufferSouthWestSend,
            BufSizeX, BufSizeY, AtNorthBoundary, AtSouthBoundary, AtEastBoundary, AtWestBoundary);
    }
    else {
        if (UseSubstrateFile)
            SubstrateInit_FromFile(SubstrateFileName, nz, MyXSlices, MyYSlices, MyXOffset, MyYOffset, id, GrainID,
                                   nzActive, BaseplateThroughPowder);
        else
            BaseplateInit_FromGrainSpacing(SubstrateGrainSpacing, nx, ny, ZMinLayer, ZMaxLayer, MyXSlices, MyYSlices,
                                           MyXOffset, MyYOffset, id, deltax, GrainID, RNGSeed,
                                           NextLayer_FirstEpitaxialGrainID, nz, BaseplateThroughPowder);
        // Separate routine for active cell data structure init for problems other than constrained solidification
        if (RemeltingYN)
            CellTypeInit_Remelt(MyXSlices, MyYSlices, LocalActiveDomainSize, CellType, CritTimeStep, id, ZBound_Low);
        else {
            CellTypeInit_NoRemelt(0, id, np, DecompositionStrategy, MyXSlices, MyYSlices, MyXOffset, MyYOffset,
                                  ZBound_Low, nz, LocalActiveDomainSize, LocalDomainSize, CellType, CritTimeStep,
                                  NeighborX, NeighborY, NeighborZ, NGrainOrientations, GrainUnitVector, DiagonalLength,
                                  GrainID, CritDiagonalLength, DOCenter, LayerID, BufferWestSend, BufferEastSend,
                                  BufferNorthSend, BufferSouthSend, BufferNorthEastSend, BufferNorthWestSend,
                                  BufferSouthEastSend, BufferSouthWestSend, BufSizeX, BufSizeY, AtNorthBoundary,
                                  AtSouthBoundary, AtEastBoundary, AtWestBoundary);
        }
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
               dTsigma, deltax, NucleiLocation, NucleationTimes_Host, NucleiGrainID, CellType, CritTimeStep,
               UndercoolingChange, LayerID, PossibleNuclei_ThisRankThisLayer, Nuclei_WholeDomain, AtNorthBoundary,
               AtSouthBoundary, AtEastBoundary, AtWestBoundary, RemeltingYN, NucleationCounter, MaxSolidificationEvents,
               NumberOfSolidificationEvents, LayerTimeTempHistory);

    // Normalize solidification parameters
    AConst = AConst * deltat / deltax;
    BConst = BConst * deltat / deltax;
    CConst = CConst * deltat / deltax;
    DConst = DConst * deltat / deltax;
    int cycle;

    // Steering Vector
    ViewI SteeringVector(Kokkos::ViewAllocateWithoutInitializing("SteeringVector"), LocalActiveDomainSize);
    ViewI_H numSteer_Host(Kokkos::ViewAllocateWithoutInitializing("SteeringVectorSize"), 1);
    numSteer_Host(0) = 0;
    ViewI numSteer = Kokkos::create_mirror_view_and_copy(device_memory_space(), numSteer_Host);

    // Update ghost node data for initial state of simulation - only needed if no remelting, as there are no active
    // cells inititially in simulations that directly model the melting process
    if ((np > 1) && (!(RemeltingYN))) {
        if (DecompositionStrategy == 1)
            GhostNodes1D(-1, id, NeighborRank_North, NeighborRank_South, MyXSlices, MyYSlices, MyXOffset, MyYOffset,
                         NeighborX, NeighborY, NeighborZ, CellType, DOCenter, GrainID, GrainUnitVector, DiagonalLength,
                         CritDiagonalLength, NGrainOrientations, BufferNorthSend, BufferSouthSend, BufferNorthRecv,
                         BufferSouthRecv, BufSizeX, BufSizeY, BufSizeZ, ZBound_Low);
        else
            GhostNodes2D(-1, id, NeighborRank_North, NeighborRank_South, NeighborRank_East, NeighborRank_West,
                         NeighborRank_NorthEast, NeighborRank_NorthWest, NeighborRank_SouthEast, NeighborRank_SouthWest,
                         MyXSlices, MyYSlices, MyXOffset, MyYOffset, NeighborX, NeighborY, NeighborZ, CellType,
                         DOCenter, GrainID, GrainUnitVector, DiagonalLength, CritDiagonalLength, NGrainOrientations,
                         BufferWestSend, BufferEastSend, BufferNorthSend, BufferSouthSend, BufferNorthEastSend,
                         BufferNorthWestSend, BufferSouthEastSend, BufferSouthWestSend, BufferWestRecv, BufferEastRecv,
                         BufferNorthRecv, BufferSouthRecv, BufferNorthEastRecv, BufferNorthWestRecv,
                         BufferSouthEastRecv, BufferSouthWestRecv, BufSizeX, BufSizeY, BufSizeZ, ZBound_Low);
    }

    // If specified, print initial values in some views for debugging purposes
    double InitTime = MPI_Wtime() - StartInitTime;
    if (id == 0)
        std::cout << "Data initialized: Time spent: " << InitTime << " s" << std::endl;
    if (PrintDebug) {
        // Host mirrors of CellType and GrainID are not maintained - pass device views and perform copy inside of
        // subroutine
        PrintExaCAData(id, -1, np, nx, ny, nz, MyXSlices, MyYSlices, MyXOffset, MyYOffset, ProcessorsInXDirection,
                       ProcessorsInYDirection, GrainID, CritTimeStep, GrainUnitVector, LayerID, CellType,
                       UndercoolingChange, UndercoolingCurrent, OutputFile, DecompositionStrategy, NGrainOrientations,
                       PathToOutput, PrintDebug, false, false, false, false, false, 0, ZBound_Low, nzActive, deltax,
                       XMin, YMin, ZMin, NumberOfLayers);
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
                               ProcessorsInXDirection, ProcessorsInYDirection, GrainID, CritTimeStep, GrainUnitVector,
                               LayerID, CellType, UndercoolingChange, UndercoolingCurrent, OutputFile,
                               DecompositionStrategy, NGrainOrientations, PathToOutput, 0, false, false, false, true,
                               false, IntermediateFileCounter, ZBound_Low, nzActive, deltax, XMin, YMin, ZMin,
                               NumberOfLayers);
                IntermediateFileCounter++;
            }
            cycle++;

            // Update cells on GPU - undercooling and diagonal length updates, nucleation
            // Cells with a successful nucleation event are marked and added to a steering vector, later dealt with in
            // CellCapture
            StartNuclTime = MPI_Wtime();
            Nucleation(cycle, SuccessfulNucEvents_ThisRank, NucleationCounter, PossibleNuclei_ThisRankThisLayer,
                       NucleationTimes_Host, NucleiLocation, NucleiGrainID, CellType, GrainID, ZBound_Low, MyXSlices,
                       MyYSlices, SteeringVector, numSteer);
            NuclTime += MPI_Wtime() - StartNuclTime;

            // Update cells on GPU - new active cells, solidification of old active cells
            // Cell capture performed in two steps - first, adding cells of interest to a steering vector (different
            // subroutine called with versus without remelting), and second, iterating over the steering vector to
            // perform active cell creation and capture operations
            StartCreateSVTime = MPI_Wtime();
            if (RemeltingYN)
                FillSteeringVector_Remelt(cycle, LocalActiveDomainSize, MyXSlices, MyYSlices, NeighborX, NeighborY,
                                          NeighborZ, CritTimeStep, UndercoolingCurrent, UndercoolingChange, CellType,
                                          GrainID, ZBound_Low, nzActive, SteeringVector, numSteer, numSteer_Host,
                                          MeltTimeStep, BufSizeX, BufSizeY, AtNorthBoundary, AtSouthBoundary,
                                          AtEastBoundary, AtWestBoundary, BufferWestSend, BufferEastSend,
                                          BufferNorthSend, BufferSouthSend, BufferNorthEastSend, BufferNorthWestSend,
                                          BufferSouthEastSend, BufferSouthWestSend, DecompositionStrategy);
            else
                FillSteeringVector_NoRemelt(cycle, LocalActiveDomainSize, MyXSlices, MyYSlices, CritTimeStep,
                                            UndercoolingCurrent, UndercoolingChange, CellType, ZBound_Low, layernumber,
                                            LayerID, SteeringVector, numSteer, numSteer_Host);
            CreateSVTime += MPI_Wtime() - StartCreateSVTime;

            StartCaptureTime = MPI_Wtime();
            CellCapture(id, np, cycle, DecompositionStrategy, LocalActiveDomainSize, LocalDomainSize, MyXSlices,
                        MyYSlices, AConst, BConst, CConst, DConst, MyXOffset, MyYOffset, NeighborX, NeighborY,
                        NeighborZ, CritTimeStep, UndercoolingCurrent, UndercoolingChange, GrainUnitVector,
                        CritDiagonalLength, DiagonalLength, CellType, DOCenter, GrainID, NGrainOrientations,
                        BufferWestSend, BufferEastSend, BufferNorthSend, BufferSouthSend, BufferNorthEastSend,
                        BufferNorthWestSend, BufferSouthEastSend, BufferSouthWestSend, BufSizeX, BufSizeY, ZBound_Low,
                        nzActive, nz, SteeringVector, numSteer, numSteer_Host, AtNorthBoundary, AtSouthBoundary,
                        AtEastBoundary, AtWestBoundary, SolidificationEventCounter, MeltTimeStep, LayerTimeTempHistory,
                        NumberOfSolidificationEvents, RemeltingYN);
            CaptureTime += MPI_Wtime() - StartCaptureTime;

            if (np > 1) {
                // Update ghost nodes
                StartGhostTime = MPI_Wtime();
                if (DecompositionStrategy == 1)
                    GhostNodes1D(cycle, id, NeighborRank_North, NeighborRank_South, MyXSlices, MyYSlices, MyXOffset,
                                 MyYOffset, NeighborX, NeighborY, NeighborZ, CellType, DOCenter, GrainID,
                                 GrainUnitVector, DiagonalLength, CritDiagonalLength, NGrainOrientations,
                                 BufferNorthSend, BufferSouthSend, BufferNorthRecv, BufferSouthRecv, BufSizeX, BufSizeY,
                                 BufSizeZ, ZBound_Low);
                else
                    GhostNodes2D(
                        cycle, id, NeighborRank_North, NeighborRank_South, NeighborRank_East, NeighborRank_West,
                        NeighborRank_NorthEast, NeighborRank_NorthWest, NeighborRank_SouthEast, NeighborRank_SouthWest,
                        MyXSlices, MyYSlices, MyXOffset, MyYOffset, NeighborX, NeighborY, NeighborZ, CellType, DOCenter,
                        GrainID, GrainUnitVector, DiagonalLength, CritDiagonalLength, NGrainOrientations,
                        BufferWestSend, BufferEastSend, BufferNorthSend, BufferSouthSend, BufferNorthEastSend,
                        BufferNorthWestSend, BufferSouthEastSend, BufferSouthWestSend, BufferWestRecv, BufferEastRecv,
                        BufferNorthRecv, BufferSouthRecv, BufferNorthEastRecv, BufferNorthWestRecv, BufferSouthEastRecv,
                        BufferSouthWestRecv, BufSizeX, BufSizeY, BufSizeZ, ZBound_Low);
                GhostTime += MPI_Wtime() - StartGhostTime;
            }

            if (cycle % 1000 == 0) {

                if (RemeltingYN)
                    IntermediateOutputAndCheck_Remelt(
                        id, np, cycle, MyXSlices, MyYSlices, MyXOffset, MyYOffset, LocalActiveDomainSize, nx, ny, nz,
                        nzActive, deltax, XMin, YMin, ZMin, DecompositionStrategy, ProcessorsInXDirection,
                        ProcessorsInYDirection, SuccessfulNucEvents_ThisRank, XSwitch, CellType, CritTimeStep, GrainID,
                        SimulationType, layernumber, NumberOfLayers, ZBound_Low, NGrainOrientations, LayerID,
                        GrainUnitVector, UndercoolingChange, UndercoolingCurrent, PathToOutput, OutputFile,
                        PrintIdleTimeSeriesFrames, TimeSeriesInc, IntermediateFileCounter, NumberOfLayers,
                        MeltTimeStep);
                else
                    IntermediateOutputAndCheck(
                        id, np, cycle, MyXSlices, MyYSlices, MyXOffset, MyYOffset, LocalDomainSize,
                        LocalActiveDomainSize, nx, ny, nz, nzActive, deltax, XMin, YMin, ZMin, DecompositionStrategy,
                        ProcessorsInXDirection, ProcessorsInYDirection, SuccessfulNucEvents_ThisRank, XSwitch, CellType,
                        CritTimeStep, GrainID, SimulationType, FinishTimeStep, layernumber, NumberOfLayers, ZBound_Low,
                        NGrainOrientations, LayerID, GrainUnitVector, UndercoolingChange, UndercoolingCurrent,
                        PathToOutput, OutputFile, PrintIdleTimeSeriesFrames, TimeSeriesInc, IntermediateFileCounter,
                        NumberOfLayers);
            }

        } while (XSwitch == 0);
        if (layernumber != NumberOfLayers - 1) {
            // Reset intermediate file counter to zero if printing video files
            if (PrintTimeSeries)
                IntermediateFileCounter = 0;

            // Determine new active cell domain size and offset from bottom of global domain
            if (RemeltingYN) {
                // Determine the bounds of the next layer: Z coordinates span ZBound_Low-ZBound_High, inclusive
                ZBound_Low =
                    calcZBound_Low_Remelt(SimulationType, LayerHeight, layernumber + 1, ZMinLayer, ZMin, deltax);
                ZBound_High = calcZBound_High_Remelt(SimulationType, SpotRadius, LayerHeight, layernumber + 1, ZMin,
                                                     deltax, nz, ZMaxLayer);
                nzActive = calcnzActive(ZBound_Low, ZBound_High, id, layernumber + 1);
                LocalActiveDomainSize = calcLocalActiveDomainSize(MyXSlices, MyYSlices,
                                                                  nzActive); // Number of active cells on this MPI rank
                // With remelting, also reinitialize temperature views back to zero and resize LayerTimeTempHistory, in
                // preparation for loading the next layer's (layernumber + 1) temperature data from RawData into the
                // temperature views
                if (SimulationType == "S")
                    TempInit_SpotRemelt(layernumber + 1, G, R, SimulationType, id, MyXSlices, MyYSlices, MyXOffset,
                                        MyYOffset, deltax, deltat, ZBound_Low, nz, LocalActiveDomainSize,
                                        LocalDomainSize, CritTimeStep, UndercoolingChange, UndercoolingCurrent,
                                        LayerHeight, FreezingRange, LayerID, NSpotsX, NSpotsY, SpotRadius, SpotOffset,
                                        LayerTimeTempHistory, NumberOfSolidificationEvents, MeltTimeStep,
                                        MaxSolidificationEvents, SolidificationEventCounter);
                else if (SimulationType == "R")
                    TempInit_ReadDataRemelt(layernumber + 1, id, MyXSlices, MyYSlices, nz, LocalActiveDomainSize,
                                            LocalDomainSize, MyXOffset, MyYOffset, deltax, deltat, FreezingRange,
                                            LayerTimeTempHistory, NumberOfSolidificationEvents, MaxSolidificationEvents,
                                            MeltTimeStep, CritTimeStep, UndercoolingChange, UndercoolingCurrent, XMin,
                                            YMin, ZMinLayer, LayerHeight, nzActive, ZBound_Low, FinishTimeStep, LayerID,
                                            FirstValue, LastValue, RawData, SolidificationEventCounter,
                                            TempFilesInSeries);
            }
            else {
                // Determine the bounds of the next layer: Z coordinates span ZBound_Low-ZBound_High, inclusive
                ZBound_Low =
                    calcZBound_Low_NoRemelt(id, MyXSlices, MyYSlices, LocalDomainSize, layernumber + 1, LayerID);
                ZBound_High =
                    calcZBound_High_NoRemelt(id, MyXSlices, MyYSlices, LocalDomainSize, layernumber + 1, LayerID);
                nzActive = calcnzActive(ZBound_Low, ZBound_High, id, layernumber + 1);
                LocalActiveDomainSize = calcLocalActiveDomainSize(MyXSlices, MyYSlices,
                                                                  nzActive); // Number of active cells on this MPI rank
            }

            // Update buffer size
            BufSizeZ = nzActive;

            // Resize and zero all view data relating to the active region from the last layer, in preparation for the
            // next layer
            ZeroResetViews(LocalActiveDomainSize, BufSizeX, BufSizeY, BufSizeZ, DiagonalLength, CritDiagonalLength,
                           DOCenter, DecompositionStrategy, BufferWestSend, BufferEastSend, BufferNorthSend,
                           BufferSouthSend, BufferNorthEastSend, BufferNorthWestSend, BufferSouthEastSend,
                           BufferSouthWestSend, BufferWestRecv, BufferEastRecv, BufferNorthRecv, BufferSouthRecv,
                           BufferNorthEastRecv, BufferNorthWestRecv, BufferSouthEastRecv, BufferSouthWestRecv,
                           SteeringVector);

            MPI_Barrier(MPI_COMM_WORLD);
            if (id == 0)
                std::cout << "Resize executed and new layer set up, GN dimensions are " << BufSizeX << " " << BufSizeY
                          << " " << BufSizeZ << std::endl;

            // If the baseplate was initialized from a substrate grain spacing, initialize powder layer grain structure
            // for the next layer "layernumber + 1" Otherwise, the entire substrate (baseplate + powder) was read from a
            // file, and the powder layers have already been initialized
            if ((!(UseSubstrateFile)) && (!(BaseplateThroughPowder)))
                PowderInit(layernumber + 1, nx, ny, LayerHeight, ZMaxLayer, ZMin, deltax, MyXSlices, MyYSlices,
                           MyXOffset, MyYOffset, id, GrainID, RNGSeed, NextLayer_FirstEpitaxialGrainID, PowderDensity);

            // Initialize active cell data structures and nuclei locations for the next layer "layernumber + 1"
            if (RemeltingYN)
                CellTypeInit_Remelt(MyXSlices, MyYSlices, LocalActiveDomainSize, CellType, CritTimeStep, id,
                                    ZBound_Low);
            else
                CellTypeInit_NoRemelt(
                    layernumber + 1, id, np, DecompositionStrategy, MyXSlices, MyYSlices, MyXOffset, MyYOffset,
                    ZBound_Low, nz, LocalActiveDomainSize, LocalDomainSize, CellType, CritTimeStep, NeighborX,
                    NeighborY, NeighborZ, NGrainOrientations, GrainUnitVector, DiagonalLength, GrainID,
                    CritDiagonalLength, DOCenter, LayerID, BufferWestSend, BufferEastSend, BufferNorthSend,
                    BufferSouthSend, BufferNorthEastSend, BufferNorthWestSend, BufferSouthEastSend, BufferSouthWestSend,
                    BufSizeX, BufSizeY, AtNorthBoundary, AtSouthBoundary, AtEastBoundary, AtWestBoundary);

            // Initialize potential nucleation event data for next layer "layernumber + 1"
            // Views containing nucleation data will be resized to the possible number of nuclei on a given MPI rank for
            // the next layer
            NucleiInit(layernumber + 1, RNGSeed, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nx, ny, nzActive,
                       ZBound_Low, id, NMax, dTN, dTsigma, deltax, NucleiLocation, NucleationTimes_Host, NucleiGrainID,
                       CellType, CritTimeStep, UndercoolingChange, LayerID, PossibleNuclei_ThisRankThisLayer,
                       Nuclei_WholeDomain, AtNorthBoundary, AtSouthBoundary, AtEastBoundary, AtWestBoundary,
                       RemeltingYN, NucleationCounter, MaxSolidificationEvents, NumberOfSolidificationEvents,
                       LayerTimeTempHistory);

            // Update ghost nodes for grain locations and attributes
            MPI_Barrier(MPI_COMM_WORLD);
            if ((np > 1) && (!(RemeltingYN))) {
                if (DecompositionStrategy == 1)
                    GhostNodes1D(-1, id, NeighborRank_North, NeighborRank_South, MyXSlices, MyYSlices, MyXOffset,
                                 MyYOffset, NeighborX, NeighborY, NeighborZ, CellType, DOCenter, GrainID,
                                 GrainUnitVector, DiagonalLength, CritDiagonalLength, NGrainOrientations,
                                 BufferNorthSend, BufferSouthSend, BufferNorthRecv, BufferSouthRecv, BufSizeX, BufSizeY,
                                 BufSizeZ, ZBound_Low);
                else
                    GhostNodes2D(-1, id, NeighborRank_North, NeighborRank_South, NeighborRank_East, NeighborRank_West,
                                 NeighborRank_NorthEast, NeighborRank_NorthWest, NeighborRank_SouthEast,
                                 NeighborRank_SouthWest, MyXSlices, MyYSlices, MyXOffset, MyYOffset, NeighborX,
                                 NeighborY, NeighborZ, CellType, DOCenter, GrainID, GrainUnitVector, DiagonalLength,
                                 CritDiagonalLength, NGrainOrientations, BufferWestSend, BufferEastSend,
                                 BufferNorthSend, BufferSouthSend, BufferNorthEastSend, BufferNorthWestSend,
                                 BufferSouthEastSend, BufferSouthWestSend, BufferWestRecv, BufferEastRecv,
                                 BufferNorthRecv, BufferSouthRecv, BufferNorthEastRecv, BufferNorthWestRecv,
                                 BufferSouthEastRecv, BufferSouthWestRecv, BufSizeX, BufSizeY, BufSizeZ, ZBound_Low);
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

    MPI_Barrier(MPI_COMM_WORLD);
    if (((PrintMisorientation) || (PrintFinalUndercoolingVals) || (PrintFullOutput)) || (PrintDefaultRVE)) {
        if (id == 0)
            std::cout << "Collecting data on rank 0 and printing to files" << std::endl;
        // Host mirrors of CellType and GrainID are not maintained - pass device views and perform copy inside of
        // subroutine
        PrintExaCAData(id, NumberOfLayers - 1, np, nx, ny, nz, MyXSlices, MyYSlices, MyXOffset, MyYOffset,
                       ProcessorsInXDirection, ProcessorsInYDirection, GrainID, CritTimeStep, GrainUnitVector, LayerID,
                       CellType, UndercoolingChange, UndercoolingCurrent, OutputFile, DecompositionStrategy,
                       NGrainOrientations, PathToOutput, 0, PrintMisorientation, PrintFinalUndercoolingVals,
                       PrintFullOutput, false, PrintDefaultRVE, 0, ZBound_Low, nzActive, deltax, XMin, YMin, ZMin,
                       NumberOfLayers, RVESize);
    }
    else {
        if (id == 0)
            std::cout << "No output files to be printed, exiting program" << std::endl;
    }

    double OutTime = MPI_Wtime() - StartOutTime;
    double InitMaxTime, InitMinTime, OutMaxTime, OutMinTime = 0.0;
    double NuclMaxTime, NuclMinTime, CreateSVMinTime, CreateSVMaxTime, CaptureMaxTime, CaptureMinTime, GhostMaxTime,
        GhostMinTime = 0.0;
    MPI_Allreduce(&InitTime, &InitMaxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&InitTime, &InitMinTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&NuclTime, &NuclMaxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&NuclTime, &NuclMinTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&CreateSVTime, &CreateSVMaxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&CreateSVTime, &CreateSVMinTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&CaptureTime, &CaptureMaxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&CaptureTime, &CaptureMinTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&GhostTime, &GhostMaxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&GhostTime, &GhostMinTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&OutTime, &OutMaxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&OutTime, &OutMinTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    PrintExaCALog(id, np, InputFile, SimulationType, DecompositionStrategy, MyXSlices, MyYSlices, MyXOffset, MyYOffset,
                  AConst, BConst, CConst, DConst, FreezingRange, deltax, NMax, dTN, dTsigma, temp_paths,
                  TempFilesInSeries, HT_deltax, RemeltingYN, deltat, NumberOfLayers, LayerHeight, SubstrateFileName,
                  SubstrateGrainSpacing, UseSubstrateFile, G, R, nx, ny, nz, FractSurfaceSitesActive, PathToOutput,
                  NSpotsX, NSpotsY, SpotOffset, SpotRadius, OutputFile, InitTime, RunTime, OutTime, cycle, InitMaxTime,
                  InitMinTime, NuclMaxTime, NuclMinTime, CreateSVMinTime, CreateSVMaxTime, CaptureMaxTime,
                  CaptureMinTime, GhostMaxTime, GhostMinTime, OutMaxTime, OutMinTime);
}
