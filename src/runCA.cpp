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

    int nx, ny, nz, NumberOfLayers, LayerHeight, TempFilesInSeries;
    int NSpotsX, NSpotsY, SpotOffset, SpotRadius, HTtoCAratio, RVESize;
    unsigned int NumberOfTemperatureDataPoints;
    int PrintDebug, TimeSeriesInc;
    bool PrintMisorientation, PrintFinalUndercoolingVals, PrintFullOutput, RemeltingYN, UseSubstrateFile,
        PrintTimeSeries, PrintIdleTimeSeriesFrames, PrintDefaultRVE, BaseplateThroughPowder, LayerwiseTempRead,
        PrintBinary;
    float SubstrateGrainSpacing;
    double HT_deltax, deltax, deltat, FractSurfaceSitesActive, G, R, NMax, dTN, dTsigma, RNGSeed, PowderDensity;
    std::string SubstrateFileName, MaterialFileName, SimulationType, OutputFile, GrainOrientationFile, PathToOutput;
    std::vector<std::string> temp_paths;

    // Read input data
    InputReadFromFile(id, InputFile, SimulationType, deltax, NMax, dTN, dTsigma, OutputFile, GrainOrientationFile,
                      TempFilesInSeries, temp_paths, HT_deltax, RemeltingYN, deltat, NumberOfLayers, LayerHeight,
                      MaterialFileName, SubstrateFileName, SubstrateGrainSpacing, UseSubstrateFile, G, R, nx, ny, nz,
                      FractSurfaceSitesActive, PathToOutput, PrintDebug, PrintMisorientation,
                      PrintFinalUndercoolingVals, PrintFullOutput, NSpotsX, NSpotsY, SpotOffset, SpotRadius,
                      PrintTimeSeries, TimeSeriesInc, PrintIdleTimeSeriesFrames, PrintDefaultRVE, RNGSeed,
                      BaseplateThroughPowder, PowderDensity, RVESize, LayerwiseTempRead, PrintBinary);
    // Read material data.
    InterfacialResponseFunction irf(id, MaterialFileName, deltat, deltax);

    // Variables characterizing local processor grids relative to global domain
    // 1D decomposition in Y: Each MPI rank has a subset consisting of of MyYSlices cells, out of ny cells in Y
    // Each MPI rank's subdomain is offset by MyYOffset cells from the lower bound of the domain in Y
    int MyYSlices, MyYOffset;
    long int LocalDomainSize;
    // Variables characterizing process IDs of neighboring MPI ranks on the grid
    // Positive Y/NegativeY directions are North/South
    int NeighborRank_North, NeighborRank_South;
    // Variables denoting whether or not each MPI rank's grid is at a global domain boundary
    bool AtNorthBoundary, AtSouthBoundary;
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
    // For simulations using input temperature data with remelting: even if only LayerwiseTempRead is true, all files
    // need to be read to determine the domain bounds
    FindXYZBounds(SimulationType, id, deltax, nx, ny, nz, temp_paths, XMin, XMax, YMin, YMax, ZMin, ZMax, LayerHeight,
                  NumberOfLayers, TempFilesInSeries, ZMinLayer, ZMaxLayer, SpotRadius);

    // Ensure that input powder layer init options are compatible with this domain size, if needed for this problem type
    if ((SimulationType == "R") || (SimulationType == "S"))
        checkPowderOverflow(nx, ny, LayerHeight, NumberOfLayers, BaseplateThroughPowder, PowderDensity);

    // Decompose the domain into subdomains on each MPI rank: Calculate MyYSlices and MyYOffset for each rank, where
    // each subdomain contains "MyYSlices" in Y, offset from the full domain origin by "MyYOffset" cells in Y
    DomainDecomposition(id, np, MyYSlices, MyYOffset, NeighborRank_North, NeighborRank_South, nx, ny, nz,
                        LocalDomainSize, AtNorthBoundary, AtSouthBoundary);

    // Read in temperature data from files, stored in "RawData", with the appropriate MPI ranks storing the appropriate
    // data
    if (SimulationType == "R")
        ReadTemperatureData(id, deltax, HT_deltax, HTtoCAratio, MyYSlices, MyYOffset, YMin, temp_paths, NumberOfLayers,
                            TempFilesInSeries, NumberOfTemperatureDataPoints, RawData, FirstValue, LastValue,
                            LayerwiseTempRead, 0);

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
            calcLocalActiveDomainSize(nx, MyYSlices, nzActive); // Number of active cells on this MPI rank
        // Initialize the temperature fields:
        // R: input temperature data from files using reduced/sparse data format (with remelting)
        // S: spot melt array test problem (with remelting)
        if (SimulationType == "R")
            TempInit_ReadDataRemelt(0, id, nx, MyYSlices, nz, LocalActiveDomainSize, LocalDomainSize, MyYOffset, deltax,
                                    deltat, irf.FreezingRange, LayerTimeTempHistory, NumberOfSolidificationEvents,
                                    MaxSolidificationEvents, MeltTimeStep, CritTimeStep, UndercoolingChange,
                                    UndercoolingCurrent, XMin, YMin, ZMinLayer, LayerHeight, nzActive, ZBound_Low,
                                    FinishTimeStep, LayerID, FirstValue, LastValue, RawData, SolidificationEventCounter,
                                    TempFilesInSeries);
        else if (SimulationType == "S")
            TempInit_SpotRemelt(0, G, R, SimulationType, id, nx, MyYSlices, MyYOffset, deltax, deltat, ZBound_Low, nz,
                                LocalActiveDomainSize, LocalDomainSize, CritTimeStep, UndercoolingChange,
                                UndercoolingCurrent, LayerHeight, irf.FreezingRange, LayerID, NSpotsX, NSpotsY,
                                SpotRadius, SpotOffset, LayerTimeTempHistory, NumberOfSolidificationEvents,
                                MeltTimeStep, MaxSolidificationEvents, SolidificationEventCounter);
    }
    else {
        // For simulations without remelting, ZBound_Low/ZBound_High/nzActive/LocalActiveDomainSize are calculated after
        // the multilayer temperature field is initialized, to avoid unneeded inclusion of cells that solidify during
        // the first layer, but remelt during the second layer Initialize the temperature fields R: input temperature
        // data from files using reduced/sparse data format (without remelting) S: spot melt array test problem (without
        // remelting) C: directional/constrained solidification test problem
        if (SimulationType == "R")
            TempInit_ReadDataNoRemelt(id, nx, MyYSlices, MyYOffset, deltax, HTtoCAratio, deltat, nz, LocalDomainSize,
                                      CritTimeStep, UndercoolingChange, XMin, YMin, ZMin, ZMinLayer, ZMaxLayer,
                                      LayerHeight, NumberOfLayers, FinishTimeStep, irf.FreezingRange, LayerID,
                                      FirstValue, LastValue, RawData);
        else if (SimulationType == "S")
            TempInit_SpotNoRemelt(G, R, SimulationType, id, nx, MyYSlices, MyYOffset, deltax, deltat, nz,
                                  LocalDomainSize, CritTimeStep, UndercoolingChange, LayerHeight, NumberOfLayers,
                                  irf.FreezingRange, LayerID, NSpotsX, NSpotsY, SpotRadius, SpotOffset);
        else if (SimulationType == "C")
            TempInit_DirSolidification(G, R, id, nx, MyYSlices, deltax, deltat, nz, LocalDomainSize, CritTimeStep,
                                       UndercoolingChange, LayerID);
        // Assign values to the bounds of the current layer now that the multilayer temperature field is initialized
        ZBound_Low = calcZBound_Low_NoRemelt(id, nx, MyYSlices, LocalDomainSize, 0, LayerID);
        ZBound_High = calcZBound_High_NoRemelt(id, nx, MyYSlices, LocalDomainSize, 0, LayerID);
        nzActive = calcnzActive(ZBound_Low, ZBound_High, id, 0);
        LocalActiveDomainSize =
            calcLocalActiveDomainSize(nx, MyYSlices, nzActive); // Number of active cells on this MPI rank
    }
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
    int BufSizeX = nx;
    int BufSizeZ = nzActive;

    // Send/recv buffers for ghost node data should be initialized with zeros
    Buffer2D BufferSouthSend("BufferSouthSend", BufSizeX * BufSizeZ, 5);
    Buffer2D BufferNorthSend("BufferNorthSend", BufSizeX * BufSizeZ, 5);
    Buffer2D BufferSouthRecv("BufferSouthRecv", BufSizeX * BufSizeZ, 5);
    Buffer2D BufferNorthRecv("BufferNorthRecv", BufSizeX * BufSizeZ, 5);

    // Initialize the grain structure and cell types - for either a constrained solidification problem, using a
    // substrate from a file, or generating a substrate using the existing CA algorithm
    int NextLayer_FirstEpitaxialGrainID;
    if (SimulationType == "C") {
        SubstrateInit_ConstrainedGrowth(id, FractSurfaceSitesActive, MyYSlices, nx, ny, MyYOffset, NeighborX, NeighborY,
                                        NeighborZ, GrainUnitVector, NGrainOrientations, CellType, GrainID,
                                        DiagonalLength, DOCenter, CritDiagonalLength, RNGSeed, np, BufferNorthSend,
                                        BufferSouthSend, BufSizeX, AtNorthBoundary, AtSouthBoundary);
    }
    else {
        if (UseSubstrateFile)
            SubstrateInit_FromFile(SubstrateFileName, nz, nx, MyYSlices, MyYOffset, id, GrainID, nzActive,
                                   BaseplateThroughPowder);
        else
            BaseplateInit_FromGrainSpacing(SubstrateGrainSpacing, nx, ny, ZMinLayer, ZMaxLayer, MyYSlices, MyYOffset,
                                           id, deltax, GrainID, RNGSeed, NextLayer_FirstEpitaxialGrainID, nz,
                                           BaseplateThroughPowder);
        // Separate routine for active cell data structure init for problems other than constrained solidification
        if (RemeltingYN)
            CellTypeInit_Remelt(nx, MyYSlices, LocalActiveDomainSize, CellType, CritTimeStep, id, ZBound_Low);
        else {
            CellTypeInit_NoRemelt(0, id, np, nx, MyYSlices, MyYOffset, ZBound_Low, nz, LocalActiveDomainSize,
                                  LocalDomainSize, CellType, CritTimeStep, NeighborX, NeighborY, NeighborZ,
                                  NGrainOrientations, GrainUnitVector, DiagonalLength, GrainID, CritDiagonalLength,
                                  DOCenter, LayerID, BufferNorthSend, BufferSouthSend, BufSizeX, AtNorthBoundary,
                                  AtSouthBoundary);
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
    NucleiInit(0, RNGSeed, MyYSlices, MyYOffset, nx, ny, nzActive, ZBound_Low, id, NMax, dTN, dTsigma, deltax,
               NucleiLocation, NucleationTimes_Host, NucleiGrainID, CellType, CritTimeStep, UndercoolingChange, LayerID,
               PossibleNuclei_ThisRankThisLayer, Nuclei_WholeDomain, AtNorthBoundary, AtSouthBoundary, RemeltingYN,
               NucleationCounter, MaxSolidificationEvents, NumberOfSolidificationEvents, LayerTimeTempHistory);

    // Steering Vector
    ViewI SteeringVector(Kokkos::ViewAllocateWithoutInitializing("SteeringVector"), LocalActiveDomainSize);
    ViewI_H numSteer_Host(Kokkos::ViewAllocateWithoutInitializing("SteeringVectorSize"), 1);
    numSteer_Host(0) = 0;
    ViewI numSteer = Kokkos::create_mirror_view_and_copy(device_memory_space(), numSteer_Host);

    // Update ghost node data for initial state of simulation - only needed if no remelting, as there are no active
    // cells inititially in simulations that directly model the melting process
    if ((np > 1) && (!(RemeltingYN))) {
        GhostNodes1D(-1, id, NeighborRank_North, NeighborRank_South, nx, MyYSlices, MyYOffset, NeighborX, NeighborY,
                     NeighborZ, CellType, DOCenter, GrainID, GrainUnitVector, DiagonalLength, CritDiagonalLength,
                     NGrainOrientations, BufferNorthSend, BufferSouthSend, BufferNorthRecv, BufferSouthRecv, BufSizeX,
                     BufSizeZ, ZBound_Low);
    }

    // If specified, print initial values in some views for debugging purposes
    double InitTime = MPI_Wtime() - StartInitTime;
    if (id == 0)
        std::cout << "Data initialized: Time spent: " << InitTime << " s" << std::endl;
    if (PrintDebug) {
        // Host mirrors of CellType and GrainID are not maintained - pass device views and perform copy inside of
        // subroutine
        PrintExaCAData(id, -1, np, nx, ny, nz, MyYSlices, MyYOffset, GrainID, CritTimeStep, GrainUnitVector, LayerID,
                       CellType, UndercoolingChange, UndercoolingCurrent, OutputFile, NGrainOrientations, PathToOutput,
                       PrintDebug, false, false, false, false, false, 0, ZBound_Low, nzActive, deltax, XMin, YMin, ZMin,
                       NumberOfLayers, PrintBinary);
        MPI_Barrier(MPI_COMM_WORLD);
        if (id == 0)
            std::cout << "Initialization data file(s) printed" << std::endl;
    }
    int cycle = 0;
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
                PrintExaCAData(id, layernumber, np, nx, ny, nz, MyYSlices, MyYOffset, GrainID, CritTimeStep,
                               GrainUnitVector, LayerID, CellType, UndercoolingChange, UndercoolingCurrent, OutputFile,
                               NGrainOrientations, PathToOutput, 0, false, false, false, true, false,
                               IntermediateFileCounter, ZBound_Low, nzActive, deltax, XMin, YMin, ZMin, NumberOfLayers,
                               PrintBinary);
                IntermediateFileCounter++;
            }
            cycle++;

            // Update cells on GPU - undercooling and diagonal length updates, nucleation
            // Cells with a successful nucleation event are marked and added to a steering vector, later dealt with in
            // CellCapture
            StartNuclTime = MPI_Wtime();
            Nucleation(cycle, SuccessfulNucEvents_ThisRank, NucleationCounter, PossibleNuclei_ThisRankThisLayer,
                       NucleationTimes_Host, NucleiLocation, NucleiGrainID, CellType, GrainID, ZBound_Low, nx,
                       MyYSlices, SteeringVector, numSteer);
            NuclTime += MPI_Wtime() - StartNuclTime;

            // Update cells on GPU - new active cells, solidification of old active cells
            // Cell capture performed in two steps - first, adding cells of interest to a steering vector (different
            // subroutine called with versus without remelting), and second, iterating over the steering vector to
            // perform active cell creation and capture operations
            StartCreateSVTime = MPI_Wtime();
            if (RemeltingYN)
                FillSteeringVector_Remelt(cycle, LocalActiveDomainSize, nx, MyYSlices, NeighborX, NeighborY, NeighborZ,
                                          CritTimeStep, UndercoolingCurrent, UndercoolingChange, CellType, GrainID,
                                          ZBound_Low, nzActive, SteeringVector, numSteer, numSteer_Host, MeltTimeStep,
                                          BufSizeX, AtNorthBoundary, AtSouthBoundary, BufferNorthSend, BufferSouthSend);
            else
                FillSteeringVector_NoRemelt(cycle, LocalActiveDomainSize, nx, MyYSlices, CritTimeStep,
                                            UndercoolingCurrent, UndercoolingChange, CellType, ZBound_Low, layernumber,
                                            LayerID, SteeringVector, numSteer, numSteer_Host);
            CreateSVTime += MPI_Wtime() - StartCreateSVTime;

            StartCaptureTime = MPI_Wtime();
            CellCapture(id, np, cycle, LocalActiveDomainSize, LocalDomainSize, nx, MyYSlices, irf, MyYOffset, NeighborX,
                        NeighborY, NeighborZ, CritTimeStep, UndercoolingCurrent, UndercoolingChange, GrainUnitVector,
                        CritDiagonalLength, DiagonalLength, CellType, DOCenter, GrainID, NGrainOrientations,
                        BufferNorthSend, BufferSouthSend, BufSizeX, ZBound_Low, nzActive, nz, SteeringVector, numSteer,
                        numSteer_Host, AtNorthBoundary, AtSouthBoundary, SolidificationEventCounter, MeltTimeStep,
                        LayerTimeTempHistory, NumberOfSolidificationEvents, RemeltingYN);
            CaptureTime += MPI_Wtime() - StartCaptureTime;

            if (np > 1) {
                // Update ghost nodes
                StartGhostTime = MPI_Wtime();
                GhostNodes1D(cycle, id, NeighborRank_North, NeighborRank_South, nx, MyYSlices, MyYOffset, NeighborX,
                             NeighborY, NeighborZ, CellType, DOCenter, GrainID, GrainUnitVector, DiagonalLength,
                             CritDiagonalLength, NGrainOrientations, BufferNorthSend, BufferSouthSend, BufferNorthRecv,
                             BufferSouthRecv, BufSizeX, BufSizeZ, ZBound_Low);
                GhostTime += MPI_Wtime() - StartGhostTime;
            }

            if (cycle % 1000 == 0) {

                if (RemeltingYN)
                    IntermediateOutputAndCheck_Remelt(
                        id, np, cycle, MyYSlices, MyYOffset, LocalActiveDomainSize, nx, ny, nz, nzActive, deltax, XMin,
                        YMin, ZMin, SuccessfulNucEvents_ThisRank, XSwitch, CellType, CritTimeStep, GrainID,
                        SimulationType, layernumber, NumberOfLayers, ZBound_Low, NGrainOrientations, LayerID,
                        GrainUnitVector, UndercoolingChange, UndercoolingCurrent, PathToOutput, OutputFile,
                        PrintIdleTimeSeriesFrames, TimeSeriesInc, IntermediateFileCounter, NumberOfLayers, MeltTimeStep,
                        PrintBinary);
                else
                    IntermediateOutputAndCheck(id, np, cycle, MyYSlices, MyYOffset, LocalDomainSize,
                                               LocalActiveDomainSize, nx, ny, nz, nzActive, deltax, XMin, YMin, ZMin,
                                               SuccessfulNucEvents_ThisRank, XSwitch, CellType, CritTimeStep, GrainID,
                                               SimulationType, FinishTimeStep, layernumber, NumberOfLayers, ZBound_Low,
                                               NGrainOrientations, LayerID, GrainUnitVector, UndercoolingChange,
                                               UndercoolingCurrent, PathToOutput, OutputFile, PrintIdleTimeSeriesFrames,
                                               TimeSeriesInc, IntermediateFileCounter, NumberOfLayers, PrintBinary);
            }

        } while (XSwitch == 0);
        if (layernumber != NumberOfLayers - 1) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (id == 0)
                std::cout << "Layer " << layernumber << " finished solidification; initializing layer "
                          << layernumber + 1 << std::endl;

            // Reset intermediate file counter to zero if printing video files
            if (PrintTimeSeries)
                IntermediateFileCounter = 0;

            // Determine new active cell domain size and offset from bottom of global domain
            if (RemeltingYN) {
                // Determine the bounds of the next layer: Z coordinates span ZBound_Low-ZBound_High, inclusive
                // If the next layer's temperature data isn't already stored, it should be read
                if ((SimulationType == "R") && (LayerwiseTempRead)) {
                    ReadTemperatureData(id, deltax, HT_deltax, HTtoCAratio, MyYSlices, MyYOffset, YMin, temp_paths,
                                        NumberOfLayers, TempFilesInSeries, NumberOfTemperatureDataPoints, RawData,
                                        FirstValue, LastValue, LayerwiseTempRead, layernumber + 1);
                }
                ZBound_Low =
                    calcZBound_Low_Remelt(SimulationType, LayerHeight, layernumber + 1, ZMinLayer, ZMin, deltax);
                ZBound_High = calcZBound_High_Remelt(SimulationType, SpotRadius, LayerHeight, layernumber + 1, ZMin,
                                                     deltax, nz, ZMaxLayer);
                nzActive = calcnzActive(ZBound_Low, ZBound_High, id, layernumber + 1);
                // Number of active cells on this MPI rank
                LocalActiveDomainSize = calcLocalActiveDomainSize(nx, MyYSlices, nzActive);
                // With remelting, also reinitialize temperature views back to zero and resize LayerTimeTempHistory, in
                // preparation for loading the next layer's (layernumber + 1) temperature data from RawData into the
                // temperature views
                if (SimulationType == "S")
                    TempInit_SpotRemelt(layernumber + 1, G, R, SimulationType, id, nx, MyYSlices, MyYOffset, deltax,
                                        deltat, ZBound_Low, nz, LocalActiveDomainSize, LocalDomainSize, CritTimeStep,
                                        UndercoolingChange, UndercoolingCurrent, LayerHeight, irf.FreezingRange,
                                        LayerID, NSpotsX, NSpotsY, SpotRadius, SpotOffset, LayerTimeTempHistory,
                                        NumberOfSolidificationEvents, MeltTimeStep, MaxSolidificationEvents,
                                        SolidificationEventCounter);
                else if (SimulationType == "R") {
                    TempInit_ReadDataRemelt(
                        layernumber + 1, id, nx, MyYSlices, nz, LocalActiveDomainSize, LocalDomainSize, MyYOffset,
                        deltax, deltat, irf.FreezingRange, LayerTimeTempHistory, NumberOfSolidificationEvents,
                        MaxSolidificationEvents, MeltTimeStep, CritTimeStep, UndercoolingChange, UndercoolingCurrent,
                        XMin, YMin, ZMinLayer, LayerHeight, nzActive, ZBound_Low, FinishTimeStep, LayerID, FirstValue,
                        LastValue, RawData, SolidificationEventCounter, TempFilesInSeries);
                }
            }
            else {
                // Determine the bounds of the next layer: Z coordinates span ZBound_Low-ZBound_High, inclusive
                ZBound_Low = calcZBound_Low_NoRemelt(id, nx, MyYSlices, LocalDomainSize, layernumber + 1, LayerID);
                ZBound_High = calcZBound_High_NoRemelt(id, nx, MyYSlices, LocalDomainSize, layernumber + 1, LayerID);
                nzActive = calcnzActive(ZBound_Low, ZBound_High, id, layernumber + 1);
                // Number of active cells on this MPI rank
                LocalActiveDomainSize = calcLocalActiveDomainSize(nx, MyYSlices, nzActive);
            }

            // Update buffer size
            BufSizeZ = nzActive;

            // Resize and zero all view data relating to the active region from the last layer, in preparation for the
            // next layer
            ZeroResetViews(LocalActiveDomainSize, BufSizeX, BufSizeZ, DiagonalLength, CritDiagonalLength, DOCenter,
                           BufferNorthSend, BufferSouthSend, BufferNorthRecv, BufferSouthRecv, SteeringVector);

            MPI_Barrier(MPI_COMM_WORLD);
            if (id == 0)
                std::cout << "Resize executed and new layer set up, GN dimensions are " << BufSizeX << " " << BufSizeZ
                          << std::endl;

            // If the baseplate was initialized from a substrate grain spacing, initialize powder layer grain structure
            // for the next layer "layernumber + 1" Otherwise, the entire substrate (baseplate + powder) was read from a
            // file, and the powder layers have already been initialized
            if ((!(UseSubstrateFile)) && (!(BaseplateThroughPowder)))
                PowderInit(layernumber + 1, nx, ny, LayerHeight, ZMaxLayer, ZMin, deltax, MyYSlices, MyYOffset, id,
                           GrainID, RNGSeed, NextLayer_FirstEpitaxialGrainID, PowderDensity);

            // Initialize active cell data structures and nuclei locations for the next layer "layernumber + 1"
            if (RemeltingYN)
                CellTypeInit_Remelt(nx, MyYSlices, LocalActiveDomainSize, CellType, CritTimeStep, id, ZBound_Low);
            else
                CellTypeInit_NoRemelt(layernumber + 1, id, np, nx, MyYSlices, MyYOffset, ZBound_Low, nz,
                                      LocalActiveDomainSize, LocalDomainSize, CellType, CritTimeStep, NeighborX,
                                      NeighborY, NeighborZ, NGrainOrientations, GrainUnitVector, DiagonalLength,
                                      GrainID, CritDiagonalLength, DOCenter, LayerID, BufferNorthSend, BufferSouthSend,
                                      BufSizeX, AtNorthBoundary, AtSouthBoundary);

            // Initialize potential nucleation event data for next layer "layernumber + 1"
            // Views containing nucleation data will be resized to the possible number of nuclei on a given MPI rank for
            // the next layer
            NucleiInit(layernumber + 1, RNGSeed, MyYSlices, MyYOffset, nx, ny, nzActive, ZBound_Low, id, NMax, dTN,
                       dTsigma, deltax, NucleiLocation, NucleationTimes_Host, NucleiGrainID, CellType, CritTimeStep,
                       UndercoolingChange, LayerID, PossibleNuclei_ThisRankThisLayer, Nuclei_WholeDomain,
                       AtNorthBoundary, AtSouthBoundary, RemeltingYN, NucleationCounter, MaxSolidificationEvents,
                       NumberOfSolidificationEvents, LayerTimeTempHistory);

            // Update ghost nodes for grain locations and attributes
            MPI_Barrier(MPI_COMM_WORLD);
            if ((np > 1) && (!(RemeltingYN))) {
                GhostNodes1D(-1, id, NeighborRank_North, NeighborRank_South, nx, MyYSlices, MyYOffset, NeighborX,
                             NeighborY, NeighborZ, CellType, DOCenter, GrainID, GrainUnitVector, DiagonalLength,
                             CritDiagonalLength, NGrainOrientations, BufferNorthSend, BufferSouthSend, BufferNorthRecv,
                             BufferSouthRecv, BufSizeX, BufSizeZ, ZBound_Low);
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
        PrintExaCAData(id, NumberOfLayers - 1, np, nx, ny, nz, MyYSlices, MyYOffset, GrainID, CritTimeStep,
                       GrainUnitVector, LayerID, CellType, UndercoolingChange, UndercoolingCurrent, OutputFile,
                       NGrainOrientations, PathToOutput, 0, PrintMisorientation, PrintFinalUndercoolingVals,
                       PrintFullOutput, false, PrintDefaultRVE, 0, ZBound_Low, nzActive, deltax, XMin, YMin, ZMin,
                       NumberOfLayers, PrintBinary, RVESize);
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

    PrintExaCALog(id, np, InputFile, SimulationType, MyYSlices, MyYOffset, irf, deltax, NMax, dTN, dTsigma, temp_paths,
                  TempFilesInSeries, HT_deltax, RemeltingYN, deltat, NumberOfLayers, LayerHeight, SubstrateFileName,
                  SubstrateGrainSpacing, UseSubstrateFile, G, R, nx, ny, nz, FractSurfaceSitesActive, PathToOutput,
                  NSpotsX, NSpotsY, SpotOffset, SpotRadius, OutputFile, InitTime, RunTime, OutTime, cycle, InitMaxTime,
                  InitMinTime, NuclMaxTime, NuclMinTime, CreateSVMinTime, CreateSVMaxTime, CaptureMaxTime,
                  CaptureMinTime, GhostMaxTime, GhostMinTime, OutMaxTime, OutMinTime, XMin, XMax, YMin, YMax, ZMin,
                  ZMax);
}
