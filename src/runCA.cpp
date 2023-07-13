// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "runCA.hpp"

#include "CAcelldata.hpp"
#include "CAghostnodes.hpp"
#include "CAinitialize.hpp"
#include "CAnucleation.hpp"
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
    int NSpotsX, NSpotsY, SpotOffset, SpotRadius, HTtoCAratio;
    bool UseSubstrateFile, BaseplateThroughPowder, LayerwiseTempRead, PowderFirstLayer;
    float SubstrateGrainSpacing;
    double HT_deltax, deltax, deltat, FractSurfaceSitesActive, G, R, NMax, dTN, dTsigma, RNGSeed, PowderActiveFraction;
    std::string SubstrateFileName, MaterialFileName, SimulationType, GrainOrientationFile;
    std::vector<std::string> temp_paths;

    // Data printing structure - contains print options (false by default) and functions
    Print print(np);

    // Read input file - toggle appropriate print options
    InputReadFromFile(id, InputFile, SimulationType, deltax, NMax, dTN, dTsigma, GrainOrientationFile,
                      TempFilesInSeries, temp_paths, HT_deltax, deltat, NumberOfLayers, LayerHeight, MaterialFileName,
                      SubstrateFileName, SubstrateGrainSpacing, UseSubstrateFile, G, R, nx, ny, nz,
                      FractSurfaceSitesActive, NSpotsX, NSpotsY, SpotOffset, SpotRadius, RNGSeed,
                      BaseplateThroughPowder, PowderActiveFraction, LayerwiseTempRead, PowderFirstLayer, print);
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
    double XMin, YMin, ZMin, XMax, YMax, ZMax; // OpenFOAM simulation bounds (if using OpenFOAM data)
    double *ZMinLayer = new double[NumberOfLayers];
    double *ZMaxLayer = new double[NumberOfLayers];
    int *FinishTimeStep = new int[NumberOfLayers];

    // Data structure for storing raw temperature data from file(s)
    // Store data as double - needed for small time steps to resolve local differences in solidification conditions
    // Each data point has 6 values (X, Y, Z coordinates, melting time, liquidus time, and cooling rate)
    // Initial estimate for size
    ViewD_H RawTemperatureData(Kokkos::ViewAllocateWithoutInitializing("RawTemperatureData"), 1000000);

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
        checkPowderOverflow(nx, ny, LayerHeight, NumberOfLayers, BaseplateThroughPowder, PowderActiveFraction);

    // Decompose the domain into subdomains on each MPI rank: Calculate MyYSlices and MyYOffset for each rank, where
    // each subdomain contains "MyYSlices" in Y, offset from the full domain origin by "MyYOffset" cells in Y
    DomainDecomposition(id, np, MyYSlices, MyYOffset, NeighborRank_North, NeighborRank_South, nx, ny, nz,
                        LocalDomainSize, AtNorthBoundary, AtSouthBoundary);

    // Read in temperature data from files, stored in "RawData", with the appropriate MPI ranks storing the appropriate
    // data
    if (SimulationType == "R")
        ReadTemperatureData(id, deltax, HT_deltax, HTtoCAratio, MyYSlices, MyYOffset, YMin, temp_paths, NumberOfLayers,
                            TempFilesInSeries, FirstValue, LastValue, LayerwiseTempRead, 0, RawTemperatureData);

    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Mesh initialized and (if being used), temperature data read" << std::endl;

    // Temperature fields characterized by these variables:
    // Current undercooling is initialized to 0
    ViewI CritTimeStep(Kokkos::ViewAllocateWithoutInitializing("CritTimeStep"), LocalDomainSize);
    ViewF UndercoolingChange(Kokkos::ViewAllocateWithoutInitializing("UndercoolingChange"), LocalDomainSize);
    ViewF UndercoolingCurrent("UndercoolingCurrent", LocalDomainSize);

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
    int ZBound_Low = calcZBound_Low(SimulationType, LayerHeight, 0, ZMinLayer, ZMin, deltax);
    int ZBound_High = calcZBound_High(SimulationType, SpotRadius, LayerHeight, 0, ZMin, deltax, nz, ZMaxLayer);
    int nzActive = calcnzActive(ZBound_Low, ZBound_High, id, 0);
    int LocalActiveDomainSize =
        calcLocalActiveDomainSize(nx, MyYSlices, nzActive); // Number of active cells on this MPI rank
    // Initialize the temperature fields:
    // R: input temperature data from files using reduced/sparse data format
    // S: spot melt array test problem
    // C: directional/constrained solidification test problem
    if (SimulationType == "R")
        TempInit_ReadData(0, id, nx, MyYSlices, nz, LocalActiveDomainSize, LocalDomainSize, MyYOffset, deltax, deltat,
                          irf.FreezingRange, LayerTimeTempHistory, NumberOfSolidificationEvents,
                          MaxSolidificationEvents, MeltTimeStep, CritTimeStep, UndercoolingChange, UndercoolingCurrent,
                          XMin, YMin, ZMinLayer, LayerHeight, nzActive, ZBound_Low, FinishTimeStep, FirstValue,
                          LastValue, RawTemperatureData, SolidificationEventCounter, TempFilesInSeries);
    else if (SimulationType == "S")
        TempInit_Spot(0, G, R, SimulationType, id, nx, MyYSlices, MyYOffset, deltax, deltat, ZBound_Low, nz,
                      LocalActiveDomainSize, LocalDomainSize, CritTimeStep, UndercoolingChange, UndercoolingCurrent,
                      LayerHeight, irf.FreezingRange, NSpotsX, NSpotsY, SpotRadius, SpotOffset, LayerTimeTempHistory,
                      NumberOfSolidificationEvents, MeltTimeStep, MaxSolidificationEvents, SolidificationEventCounter);
    else if (SimulationType == "C")
        TempInit_DirSolidification(G, R, id, nx, MyYSlices, deltax, deltat, nz, LocalDomainSize, CritTimeStep,
                                   UndercoolingChange, NumberOfSolidificationEvents, SolidificationEventCounter,
                                   MeltTimeStep, MaxSolidificationEvents, LayerTimeTempHistory);
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

    // Variables characterizing the active cell region within each rank's grid
    // DiagonalLengths should be initialized to 0, DOCenter/CritDiagonalLength do not need initialization
    ViewF DiagonalLength("DiagonalLength", LocalActiveDomainSize);
    ViewF DOCenter(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * LocalActiveDomainSize);
    ViewF CritDiagonalLength(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength"), 26 * LocalActiveDomainSize);

    // Buffers for ghost node data (fixed size)
    int BufSizeInitialEstimate = 25;
    int BufSize = BufSizeInitialEstimate; // set to initial estimate
    Buffer2D BufferSouthSend(Kokkos::ViewAllocateWithoutInitializing("BufferSouthSend"), BufSize, 8);
    Buffer2D BufferNorthSend(Kokkos::ViewAllocateWithoutInitializing("BufferNorthSend"), BufSize, 8);
    Buffer2D BufferSouthRecv(Kokkos::ViewAllocateWithoutInitializing("BufferSouthRecv"), BufSize, 8);
    Buffer2D BufferNorthRecv(Kokkos::ViewAllocateWithoutInitializing("BufferNorthRecv"), BufSize, 8);
    ViewI SendSizeSouth(Kokkos::ViewAllocateWithoutInitializing("SendSizeSouth"), 1);
    ViewI SendSizeNorth(Kokkos::ViewAllocateWithoutInitializing("SendSizeNorth"), 1);
    ViewI_H SendSizeSouth_Host(Kokkos::ViewAllocateWithoutInitializing("SendSizeSouth_Host"), 1);
    ViewI_H SendSizeNorth_Host(Kokkos::ViewAllocateWithoutInitializing("SendSizeNorth_Host"), 1);
    // Send/recv buffers for ghost node data should be initialized with -1s in the first index as placeholders for empty
    // positions in the buffer, and with send size counts of 0
    ResetSendBuffers(BufSize, BufferNorthSend, BufferSouthSend, SendSizeNorth, SendSizeSouth);

    // Initialize cell types, grain IDs, and layer IDs
    CellData<device_memory_space> cellData(LocalDomainSize, LocalActiveDomainSize, nx, MyYSlices, ZBound_Low);
    if (SimulationType == "C")
        cellData.init_substrate(id, FractSurfaceSitesActive, MyYSlices, nx, ny, MyYOffset, NeighborX, NeighborY,
                                NeighborZ, GrainUnitVector, NGrainOrientations, DiagonalLength, DOCenter,
                                CritDiagonalLength, RNGSeed);
    else
        cellData.init_substrate(SubstrateFileName, UseSubstrateFile, BaseplateThroughPowder, PowderFirstLayer, nx, ny,
                                nz, LayerHeight, LocalActiveDomainSize, ZMaxLayer, ZMin, deltax, MyYSlices, MyYOffset,
                                ZBound_Low, id, RNGSeed, SubstrateGrainSpacing, PowderActiveFraction,
                                NumberOfSolidificationEvents);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Grain struct initialized" << std::endl;

    // Nucleation data structure, containing views of nuclei locations, time steps, and ids, and nucleation event
    // counters - initialized with an estimate on the number of nuclei in the layer Without knowing
    // PossibleNuclei_ThisRankThisLayer yet, initialize nucleation data structures to estimated sizes, resize inside of
    // NucleiInit when the number of nuclei per rank is known
    int EstimatedNuclei_ThisRankThisLayer = NMax * pow(deltax, 3) * LocalActiveDomainSize;
    Nucleation<device_memory_space> nucleation(EstimatedNuclei_ThisRankThisLayer, NMax, deltax);
    // Fill in nucleation data structures, and assign nucleation undercooling values to potential nucleation events
    // Potential nucleation grains are only associated with liquid cells in layer 0 - they will be initialized for each
    // successive layer when layer 0 in complete
    nucleation.placeNuclei(MaxSolidificationEvents, NumberOfSolidificationEvents, LayerTimeTempHistory, RNGSeed, 0, nx,
                           ny, nzActive, dTN, dTsigma, MyYSlices, MyYOffset, ZBound_Low, id, AtNorthBoundary,
                           AtSouthBoundary);

    // Steering Vector
    ViewI SteeringVector(Kokkos::ViewAllocateWithoutInitializing("SteeringVector"), LocalActiveDomainSize);
    ViewI_H numSteer_Host(Kokkos::ViewAllocateWithoutInitializing("SteeringVectorSize"), 1);
    numSteer_Host(0) = 0;
    ViewI numSteer = Kokkos::create_mirror_view_and_copy(device_memory_space(), numSteer_Host);

    // Get the size of buffer data for sending/receiving print information before printing any fields
    print.getSendRecvDataSizes(nx, ny, nz, MyYSlices, MyYOffset, np);
    // If specified, print initial values in some views for debugging purposes
    double InitTime = MPI_Wtime() - StartInitTime;
    if (id == 0)
        std::cout << "Data initialized: Time spent: " << InitTime << " s" << std::endl;
    print.printInitExaCAData(id, np, nx, ny, MyYSlices, nzActive, deltax, XMin, YMin, ZMin, cellData.GrainID_AllLayers,
                             cellData.LayerID_AllLayers, MeltTimeStep, CritTimeStep, UndercoolingChange);
    MPI_Barrier(MPI_COMM_WORLD);
    int cycle = 0;
    double StartRunTime = MPI_Wtime();
    for (int layernumber = 0; layernumber < NumberOfLayers; layernumber++) {

        int XSwitch = 0;
        double LayerTime1 = MPI_Wtime();

        // Loop continues until all liquid cells claimed by solid grains
        do {

            // Start of time step - check and see if intermediate system output is to be printed to files
            if ((print.PrintTimeSeries) && (cycle % print.TimeSeriesInc == 0)) {
                // Print current state of ExaCA simulation (up to and including the current layer's data)
                print.printIntermediateGrainMisorientation(
                    id, np, cycle, nx, ny, nz, MyYSlices, nzActive, deltax, XMin, YMin, ZMin,
                    cellData.GrainID_AllLayers, cellData.LayerID_AllLayers, cellData.CellType_AllLayers,
                    GrainUnitVector, NGrainOrientations, layernumber, ZBound_Low);
            }
            cycle++;

            // Update cells on GPU - undercooling and diagonal length updates, nucleation
            // Cells with a successful nucleation event are marked and added to a steering vector, later dealt with in
            // CellCapture
            StartNuclTime = MPI_Wtime();
            nucleation.nucleate_grain(cycle, cellData, ZBound_Low, nx, MyYSlices, SteeringVector, numSteer);
            NuclTime += MPI_Wtime() - StartNuclTime;

            // Update cells on GPU - new active cells, solidification of old active cells
            // Cell capture performed in two steps - first, adding cells of interest to a steering vector (different
            // subroutine called with versus without remelting), and second, iterating over the steering vector to
            // perform active cell creation and capture operations
            // Constrained/directional solidification problem cells only undergo 1 solidificaton event and are
            // initialized as liquid - the steering vector operation for this problem can be constructed using
            // FillSteeringVector_NoRemelt (a simplified version of FillSteeringVector_Remelt
            StartCreateSVTime = MPI_Wtime();
            if (SimulationType != "C")
                FillSteeringVector_Remelt(cycle, LocalActiveDomainSize, nx, MyYSlices, NeighborX, NeighborY, NeighborZ,
                                          CritTimeStep, UndercoolingCurrent, UndercoolingChange, cellData, ZBound_Low,
                                          nzActive, SteeringVector, numSteer, numSteer_Host, MeltTimeStep,
                                          SolidificationEventCounter, NumberOfSolidificationEvents,
                                          LayerTimeTempHistory);
            else
                FillSteeringVector_NoRemelt(cycle, LocalActiveDomainSize, nx, MyYSlices, CritTimeStep,
                                            UndercoolingCurrent, UndercoolingChange, cellData, ZBound_Low, layernumber,
                                            SteeringVector, numSteer, numSteer_Host);
            CreateSVTime += MPI_Wtime() - StartCreateSVTime;

            StartCaptureTime = MPI_Wtime();
            CellCapture(id, np, cycle, LocalActiveDomainSize, LocalDomainSize, nx, MyYSlices, irf, MyYOffset, NeighborX,
                        NeighborY, NeighborZ, CritTimeStep, UndercoolingCurrent, UndercoolingChange, GrainUnitVector,
                        CritDiagonalLength, DiagonalLength, cellData, DOCenter, NGrainOrientations, BufferNorthSend,
                        BufferSouthSend, SendSizeNorth, SendSizeSouth, ZBound_Low, nzActive, nz, SteeringVector,
                        numSteer, numSteer_Host, AtNorthBoundary, AtSouthBoundary, SolidificationEventCounter,
                        LayerTimeTempHistory, NumberOfSolidificationEvents, BufSize);
            // Count the number of cells' in halo regions where the data did not fit into the send buffers
            // Reduce across all ranks, as the same BufSize should be maintained across all ranks
            // If any rank overflowed its buffer size, resize all buffers to the new size plus 10% padding
            int OldBufSize = BufSize;
            BufSize = ResizeBuffers(BufferNorthSend, BufferSouthSend, BufferNorthRecv, BufferSouthRecv, SendSizeNorth,
                                    SendSizeSouth, SendSizeNorth_Host, SendSizeSouth_Host, OldBufSize);
            if (OldBufSize != BufSize) {
                if (id == 0)
                    std::cout << "Resized number of cells stored in send/recv buffers from " << OldBufSize << " to "
                              << BufSize << std::endl;
                RefillBuffers(nx, nzActive, MyYSlices, ZBound_Low, cellData, BufferNorthSend, BufferSouthSend,
                              SendSizeNorth, SendSizeSouth, AtNorthBoundary, AtSouthBoundary, DOCenter, DiagonalLength,
                              NGrainOrientations, BufSize);
            }
            CaptureTime += MPI_Wtime() - StartCaptureTime;

            if (np > 1) {
                // Update ghost nodes
                StartGhostTime = MPI_Wtime();
                GhostNodes1D(cycle, id, NeighborRank_North, NeighborRank_South, nx, MyYSlices, MyYOffset, NeighborX,
                             NeighborY, NeighborZ, cellData, DOCenter, GrainUnitVector, DiagonalLength,
                             CritDiagonalLength, NGrainOrientations, BufferNorthSend, BufferSouthSend, BufferNorthRecv,
                             BufferSouthRecv, BufSize, ZBound_Low, SendSizeNorth, SendSizeSouth);
                GhostTime += MPI_Wtime() - StartGhostTime;
            }

            if (cycle % 1000 == 0) {
                IntermediateOutputAndCheck(id, np, cycle, MyYSlices, LocalActiveDomainSize, nx, ny, nz, nzActive,
                                           deltax, XMin, YMin, ZMin, nucleation.SuccessfulNucleationCounter, XSwitch,
                                           cellData, CritTimeStep, SimulationType, layernumber, NumberOfLayers,
                                           ZBound_Low, NGrainOrientations, GrainUnitVector, print, MeltTimeStep);
            }

        } while (XSwitch == 0);
        if (layernumber != NumberOfLayers - 1) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (id == 0)
                std::cout << "Layer " << layernumber << " finished solidification; initializing layer "
                          << layernumber + 1 << std::endl;

            // Reset intermediate file counter to zero if printing video files
            if (print.PrintTimeSeries)
                print.resetIntermediateFileCounter();

            // Determine new active cell domain size and offset from bottom of global domain
            ZBound_Low = calcZBound_Low(SimulationType, LayerHeight, layernumber + 1, ZMinLayer, ZMin, deltax);
            ZBound_High =
                calcZBound_High(SimulationType, SpotRadius, LayerHeight, layernumber + 1, ZMin, deltax, nz, ZMaxLayer);
            nzActive = calcnzActive(ZBound_Low, ZBound_High, id, layernumber + 1);
            LocalActiveDomainSize = calcLocalActiveDomainSize(nx, MyYSlices, nzActive);
            // Determine the bounds of the next layer: Z coordinates span ZBound_Low-ZBound_High, inclusive
            // If the next layer's temperature data isn't already stored, it should be read
            if ((SimulationType == "R") && (LayerwiseTempRead)) {
                ReadTemperatureData(id, deltax, HT_deltax, HTtoCAratio, MyYSlices, MyYOffset, YMin, temp_paths,
                                    NumberOfLayers, TempFilesInSeries, FirstValue, LastValue, LayerwiseTempRead,
                                    layernumber + 1, RawTemperatureData);
            }
            // Reinitialize temperature views back to zero and resize LayerTimeTempHistory, in
            // preparation for loading the next layer's (layernumber + 1) temperature data from RawData into the
            // temperature views
            if (SimulationType == "S")
                TempInit_Spot(layernumber + 1, G, R, SimulationType, id, nx, MyYSlices, MyYOffset, deltax, deltat,
                              ZBound_Low, nz, LocalActiveDomainSize, LocalDomainSize, CritTimeStep, UndercoolingChange,
                              UndercoolingCurrent, LayerHeight, irf.FreezingRange, NSpotsX, NSpotsY, SpotRadius,
                              SpotOffset, LayerTimeTempHistory, NumberOfSolidificationEvents, MeltTimeStep,
                              MaxSolidificationEvents, SolidificationEventCounter);
            else if (SimulationType == "R") {
                TempInit_ReadData(layernumber + 1, id, nx, MyYSlices, nz, LocalActiveDomainSize, LocalDomainSize,
                                  MyYOffset, deltax, deltat, irf.FreezingRange, LayerTimeTempHistory,
                                  NumberOfSolidificationEvents, MaxSolidificationEvents, MeltTimeStep, CritTimeStep,
                                  UndercoolingChange, UndercoolingCurrent, XMin, YMin, ZMinLayer, LayerHeight, nzActive,
                                  ZBound_Low, FinishTimeStep, FirstValue, LastValue, RawTemperatureData,
                                  SolidificationEventCounter, TempFilesInSeries);
            }

            // Resize and zero all view data relating to the active region from the last layer, in preparation for the
            // next layer
            ZeroResetViews(LocalActiveDomainSize, DiagonalLength, CritDiagonalLength, DOCenter, SteeringVector);
            MPI_Barrier(MPI_COMM_WORLD);

            // Sets up views, powder layer (if necessary), and cell types for the next layer of a multilayer problem
            cellData.init_next_layer(layernumber + 1, id, nx, ny, MyYSlices, MyYOffset, ZBound_Low, LayerHeight,
                                     LocalActiveDomainSize, UseSubstrateFile, BaseplateThroughPowder, ZMin, ZMaxLayer,
                                     deltax, RNGSeed, PowderActiveFraction, NumberOfSolidificationEvents);

            // Initialize potential nucleation event data for next layer "layernumber + 1"
            // Views containing nucleation data will be resized to the possible number of nuclei on a given MPI rank for
            // the next layer
            nucleation.resetNucleiCounters(); // start counters at 0
            nucleation.placeNuclei(MaxSolidificationEvents, NumberOfSolidificationEvents, LayerTimeTempHistory, RNGSeed,
                                   layernumber + 1, nx, ny, nzActive, dTN, dTsigma, MyYSlices, MyYOffset, ZBound_Low,
                                   id, AtNorthBoundary, AtSouthBoundary);

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
    // Collect and print specified final fields to output files
    print.printFinalExaCAData(id, np, nx, ny, nz, MyYSlices, NumberOfLayers, cellData.LayerID_AllLayers,
                              cellData.CellType_AllLayers, cellData.GrainID_AllLayers, UndercoolingCurrent,
                              UndercoolingChange, MeltTimeStep, CritTimeStep, GrainUnitVector, NGrainOrientations,
                              deltax, XMin, YMin, ZMin);

    // Calculate volume fraction of solidified domain consisting of nucleated grains
    float VolFractionNucleated =
        calcVolFractionNucleated(id, nx, MyYSlices, LocalDomainSize, cellData.LayerID_AllLayers,
                                 cellData.GrainID_AllLayers, AtNorthBoundary, AtSouthBoundary);

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

    // Print the log file with JSON format, timing information to the console
    PrintExaCALog(id, np, InputFile, print.PathToOutput, print.BaseFileName, SimulationType, MyYSlices, MyYOffset, irf,
                  deltax, NMax, dTN, dTsigma, temp_paths, TempFilesInSeries, HT_deltax, deltat, NumberOfLayers,
                  LayerHeight, SubstrateFileName, SubstrateGrainSpacing, UseSubstrateFile, G, R, nx, ny, nz,
                  FractSurfaceSitesActive, NSpotsX, NSpotsY, SpotOffset, SpotRadius, InitTime, RunTime, OutTime, cycle,
                  InitMaxTime, InitMinTime, NuclMaxTime, NuclMinTime, CreateSVMinTime, CreateSVMaxTime, CaptureMaxTime,
                  CaptureMinTime, GhostMaxTime, GhostMinTime, OutMaxTime, OutMinTime, XMin, XMax, YMin, YMax, ZMin,
                  ZMax, GrainOrientationFile, VolFractionNucleated);
    if (id == 0)
        PrintExaCATiming(np, InitTime, RunTime, OutTime, cycle, InitMaxTime, InitMinTime, NuclMaxTime, NuclMinTime,
                         CreateSVMinTime, CreateSVMaxTime, CaptureMaxTime, CaptureMinTime, GhostMaxTime, GhostMinTime,
                         OutMaxTime, OutMinTime);
}
