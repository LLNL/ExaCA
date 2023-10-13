// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "ExaCA.hpp"

#include "mpi.h"

#include <string>
#include <vector>

void RunProgram_Reduced(int id, int np, std::string InputFile) {
    double NuclTime = 0.0, CreateSVTime = 0.0, CaptureTime = 0.0, GhostTime = 0.0;
    double StartNuclTime, StartCreateSVTime, StartCaptureTime, StartGhostTime;
    double StartInitTime = MPI_Wtime();

    // Read input file
    Inputs inputs(id, InputFile);

    // These parameters from the inputs struct will eventually be stored in the grid struct, but for now are temporarily
    // kept in the current scope
    int nx = inputs.domain.nx;
    int ny = inputs.domain.ny;
    int nz = inputs.domain.nz;
    int NumberOfLayers = inputs.domain.NumberOfLayers;
    int LayerHeight = inputs.domain.LayerHeight;
    double deltax = inputs.domain.deltax;
    double HT_deltax = inputs.temperature.HT_deltax;

    // Material response function
    InterfacialResponseFunction irf(id, inputs.MaterialFileName, inputs.domain.deltat, deltax);

    // Variables characterizing local processor grids relative to global domain
    // 1D decomposition in Y: Each MPI rank has a subset consisting of of ny_local cells, out of ny cells in Y
    // Each MPI rank's subdomain is offset by y_offset cells from the lower bound of the domain in Y
    int ny_local, y_offset, DomainSize, DomainSize_AllLayers;
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

    // Intialize neighbor list structures (NeighborX, NeighborY, NeighborZ)
    NeighborListInit(NeighborX, NeighborY, NeighborZ);

    // Obtain the physical XYZ bounds of the domain, using either domain size from the input file, or reading
    // temperature data files and parsing the coordinates
    // For simulations using input temperature data with remelting: even if only LayerwiseTempRead is true, all files
    // need to be read to determine the domain bounds
    FindXYZBounds(id, deltax, nx, ny, nz, XMin, XMax, YMin, YMax, ZMin, ZMax, ZMinLayer, ZMaxLayer, NumberOfLayers,
                  LayerHeight, inputs);

    // Ensure that input powder layer init options are compatible with this domain size, if needed for this problem type
    if ((inputs.SimulationType == "R") || (inputs.SimulationType == "S"))
        checkPowderOverflow(nx, ny, LayerHeight, NumberOfLayers, inputs);

    // Decompose the domain into subdomains on each MPI rank: Calculate ny_local and y_offset for each rank, where
    // each subdomain contains "ny_local" in Y, offset from the full domain origin by "y_offset" cells in Y
    DomainDecomposition(id, np, ny_local, y_offset, NeighborRank_North, NeighborRank_South, nx, ny, nz,
                        DomainSize_AllLayers, AtNorthBoundary, AtSouthBoundary);
    // Bounds of the current layer: Z coordinates span z_layer_bottom-z_layer_top, inclusive
    int z_layer_bottom = calc_z_layer_bottom(inputs.SimulationType, LayerHeight, 0, ZMinLayer, ZMin, deltax);
    int z_layer_top =
        calc_z_layer_top(inputs.SimulationType, inputs.domain.SpotRadius, LayerHeight, 0, ZMin, deltax, nz, ZMaxLayer);
    int nz_layer = calc_nz_layer(z_layer_bottom, z_layer_top, id, 0);
    DomainSize = calcLayerDomainSize(nx, ny_local, nz_layer); // Number of cells in the current layer on this MPI rank
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Mesh initialized: initial domain size is " << nz_layer << " out of " << nz
                  << " total cells in the Z direction" << std::endl;

    // Temperature fields characterized by data in this structure
    Temperature<device_memory_space> temperature(DomainSize, NumberOfLayers, inputs.temperature);
    // Read temperature data if necessary
    if (inputs.SimulationType == "R")
        temperature.readTemperatureData(id, deltax, y_offset, ny_local, YMin, NumberOfLayers, 0);
    // Initialize the temperature fields for the simualtion type of interest
    if ((inputs.SimulationType == "C") || (inputs.SimulationType == "SingleGrain")) {
        if (inputs.temperature.G == 0)
            temperature.initialize(id, DomainSize, inputs.domain.deltat);
        else
            temperature.initialize(id, inputs.SimulationType, nx, ny_local, nz, deltax, DomainSize,
                                   inputs.domain.deltat);
    }
    else if (inputs.SimulationType == "S")
        temperature.initialize(id, nx, ny_local, y_offset, deltax, DomainSize, irf.FreezingRange, inputs,
                               NumberOfLayers);
    else if (inputs.SimulationType == "R")
        temperature.initialize(0, id, nx, ny_local, DomainSize, y_offset, deltax, irf.FreezingRange, XMin, YMin,
                               ZMinLayer, LayerHeight, nz_layer, z_layer_bottom, FinishTimeStep, inputs.domain.deltat);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Done with temperature field initialization" << std::endl;

    int NGrainOrientations = 0; // Number of grain orientations considered in the simulation
    // No initialize size yet, will be resized in OrientationInit
    ViewF GrainUnitVector(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), 0);

    // Initialize grain orientations
    OrientationInit(id, NGrainOrientations, GrainUnitVector, inputs.GrainOrientationFile);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Done with orientation initialization " << std::endl;

    // Variables characterizing the active cell region within each rank's grid
    // DiagonalLengths should be initialized to 0, DOCenter/CritDiagonalLength do not need initialization
    ViewF DiagonalLength("DiagonalLength", DomainSize);
    ViewF DOCenter(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * DomainSize);
    ViewF CritDiagonalLength(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength"), 26 * DomainSize);

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
    CellData<device_memory_space> cellData(DomainSize_AllLayers, DomainSize, nx, ny_local, z_layer_bottom,
                                           inputs.substrate);
    if (inputs.SimulationType == "C")
        cellData.init_substrate(id, ny_local, nx, ny, y_offset, inputs.RNGSeed);
    else if (inputs.SimulationType == "SingleGrain")
        cellData.init_substrate(id, nx, ny, nz, ny_local, y_offset, DomainSize);
    else
        cellData.init_substrate(nx, ny, nz, DomainSize, ZMaxLayer, ZMin, deltax, ny_local, y_offset, z_layer_bottom, id,
                                inputs.RNGSeed, temperature.NumberOfSolidificationEvents);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Grain struct initialized" << std::endl;

    // Nucleation data structure, containing views of nuclei locations, time steps, and ids, and nucleation event
    // counters - initialized with an estimate on the number of nuclei in the layer Without knowing
    // PossibleNuclei_ThisRankThisLayer yet, initialize nucleation data structures to estimated sizes, resize inside of
    // NucleiInit when the number of nuclei per rank is known
    int EstimatedNuclei_ThisRankThisLayer = inputs.nucleation.NMax * pow(deltax, 3) * DomainSize;
    Nucleation<device_memory_space> nucleation(EstimatedNuclei_ThisRankThisLayer, deltax, inputs.nucleation);
    // Fill in nucleation data structures, and assign nucleation undercooling values to potential nucleation events
    // Potential nucleation grains are only associated with liquid cells in layer 0 - they will be initialized for each
    // successive layer when layer 0 in complete
    nucleation.placeNuclei(temperature, inputs.RNGSeed, 0, nx, ny, nz_layer, ny_local, y_offset, z_layer_bottom, id,
                           AtNorthBoundary, AtSouthBoundary);

    // Steering Vector
    ViewI SteeringVector(Kokkos::ViewAllocateWithoutInitializing("SteeringVector"), DomainSize);
    ViewI_H numSteer_Host(Kokkos::ViewAllocateWithoutInitializing("SteeringVectorSize"), 1);
    numSteer_Host(0) = 0;
    ViewI numSteer = Kokkos::create_mirror_view_and_copy(device_memory_space(), numSteer_Host);

    // Initialize printing struct from inputs
    Print print(nx, ny, nz, ny_local, y_offset, np, inputs.print);

    // If specified, print initial values in some views for debugging purposes
    double InitTime = MPI_Wtime() - StartInitTime;
    if (id == 0)
        std::cout << "Data initialized: Time spent: " << InitTime << " s" << std::endl;
    print.printInitExaCAData(id, np, nx, ny, ny_local, nz_layer, DomainSize, deltax, XMin, YMin, ZMin,
                             cellData.GrainID_AllLayers, cellData.LayerID_AllLayers, temperature);
    MPI_Barrier(MPI_COMM_WORLD);
    int cycle = 0;
    double StartRunTime = MPI_Wtime();
    for (int layernumber = 0; layernumber < NumberOfLayers; layernumber++) {

        int XSwitch = 0;
        double LayerTime1 = MPI_Wtime();

        // Loop continues until all liquid cells claimed by solid grains
        do {

            // Start of time step - optional print current state of ExaCA simulation (up to and including the current
            // layer's data)
            print.printIntermediateGrainMisorientation(
                id, np, cycle, nx, ny, nz, ny_local, nz_layer, deltax, XMin, YMin, ZMin, cellData.GrainID_AllLayers,
                cellData.CellType_AllLayers, GrainUnitVector, NGrainOrientations, layernumber, z_layer_bottom);
            cycle++;

            // Update cells on GPU - undercooling and diagonal length updates, nucleation
            // Cells with a successful nucleation event are marked and added to a steering vector, later dealt with in
            // CellCapture
            StartNuclTime = MPI_Wtime();
            nucleation.nucleate_grain(cycle, cellData, z_layer_bottom, nx, ny_local, SteeringVector, numSteer);
            NuclTime += MPI_Wtime() - StartNuclTime;

            // Update cells on GPU - new active cells, solidification of old active cells
            // Cell capture performed in two steps - first, adding cells of interest to a steering vector (different
            // subroutine called with versus without remelting), and second, iterating over the steering vector to
            // perform active cell creation and capture operations
            // Constrained/directional solidification problem cells only undergo 1 solidificaton event and are
            // initialized as liquid - the steering vector operation for this problem can be constructed using
            // FillSteeringVector_NoRemelt (a simplified version of FillSteeringVector_Remelt
            StartCreateSVTime = MPI_Wtime();

            if ((inputs.SimulationType == "C") || (inputs.SimulationType == "SingleGrain"))
                FillSteeringVector_NoRemelt(cycle, DomainSize, temperature, cellData, SteeringVector, numSteer,
                                            numSteer_Host);
            else
                FillSteeringVector_Remelt(cycle, DomainSize, nx, ny_local, NeighborX, NeighborY, NeighborZ, temperature,
                                          cellData, nz_layer, SteeringVector, numSteer, numSteer_Host);
            CreateSVTime += MPI_Wtime() - StartCreateSVTime;

            StartCaptureTime = MPI_Wtime();
            CellCapture(id, np, cycle, nx, ny_local, irf, y_offset, NeighborX, NeighborY, NeighborZ, GrainUnitVector,
                        CritDiagonalLength, DiagonalLength, cellData, temperature, DOCenter, NGrainOrientations,
                        BufferNorthSend, BufferSouthSend, SendSizeNorth, SendSizeSouth, nz_layer, SteeringVector,
                        numSteer, numSteer_Host, AtNorthBoundary, AtSouthBoundary, BufSize);
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
                RefillBuffers(nx, nz_layer, ny_local, cellData, BufferNorthSend, BufferSouthSend, SendSizeNorth,
                              SendSizeSouth, AtNorthBoundary, AtSouthBoundary, DOCenter, DiagonalLength,
                              NGrainOrientations, BufSize);
            }
            CaptureTime += MPI_Wtime() - StartCaptureTime;

            if (np > 1) {
                // Update ghost nodes
                StartGhostTime = MPI_Wtime();
                GhostNodes1D(cycle, id, NeighborRank_North, NeighborRank_South, nx, ny_local, y_offset, NeighborX,
                             NeighborY, NeighborZ, cellData, DOCenter, GrainUnitVector, DiagonalLength,
                             CritDiagonalLength, NGrainOrientations, BufferNorthSend, BufferSouthSend, BufferNorthRecv,
                             BufferSouthRecv, BufSize, SendSizeNorth, SendSizeSouth);
                GhostTime += MPI_Wtime() - StartGhostTime;
            }

            if ((cycle % 1000 == 0) && (inputs.SimulationType != "SingleGrain")) {
                IntermediateOutputAndCheck(id, np, cycle, ny_local, DomainSize, nx, ny, nz, nz_layer, z_layer_bottom,
                                           deltax, XMin, YMin, ZMin, nucleation.SuccessfulNucleationCounter, XSwitch,
                                           cellData, temperature, inputs.SimulationType, layernumber,
                                           NGrainOrientations, GrainUnitVector, print);
            }
            else if (inputs.SimulationType == "SingleGrain") {
                IntermediateOutputAndCheck(id, cycle, ny_local, y_offset, DomainSize, nx, ny, nz, XSwitch,
                                           cellData.CellType_AllLayers);
            }

        } while (XSwitch == 0);
        if (layernumber != NumberOfLayers - 1) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (id == 0)
                std::cout << "Layer " << layernumber << " finished solidification; initializing layer "
                          << layernumber + 1 << std::endl;

            // Reset intermediate file counter to zero if printing video files
            print.resetIntermediateFileCounter();

            // Determine new active cell domain size and offset from bottom of global domain
            z_layer_bottom =
                calc_z_layer_bottom(inputs.SimulationType, LayerHeight, layernumber + 1, ZMinLayer, ZMin, deltax);
            z_layer_top = calc_z_layer_top(inputs.SimulationType, inputs.domain.SpotRadius, LayerHeight,
                                           layernumber + 1, ZMin, deltax, nz, ZMaxLayer);
            nz_layer = calc_nz_layer(z_layer_bottom, z_layer_top, id, layernumber + 1);
            DomainSize = calcLayerDomainSize(nx, ny_local, nz_layer);
            // Determine the bounds of the next layer: Z coordinates span z_layer_bottom-z_layer_top, inclusive
            // For simulation type R, need to initialize new temperature field data for layer "layernumber + 1"
            if (inputs.SimulationType == "R") {
                // If the next layer's temperature data isn't already stored, it should be read
                if (inputs.temperature.LayerwiseTempRead)
                    temperature.readTemperatureData(id, deltax, y_offset, ny_local, YMin, NumberOfLayers,
                                                    layernumber + 1);
                // Initialize next layer's temperature data
                temperature.initialize(layernumber + 1, id, nx, ny_local, DomainSize, y_offset, deltax,
                                       irf.FreezingRange, XMin, YMin, ZMinLayer, LayerHeight, nz_layer, z_layer_bottom,
                                       FinishTimeStep, inputs.domain.deltat);
            }

            // Reset initial undercooling/solidification event counter of all cells to zeros for the next layer,
            // resizing the views to the number of cells associated with the next layer
            temperature.reset_layer_events_undercooling(DomainSize);

            // Resize and zero all view data relating to the active region from the last layer, in preparation for the
            // next layer
            ZeroResetViews(DomainSize, DiagonalLength, CritDiagonalLength, DOCenter, SteeringVector);
            MPI_Barrier(MPI_COMM_WORLD);

            // Sets up views, powder layer (if necessary), and cell types for the next layer of a multilayer problem
            cellData.init_next_layer(layernumber + 1, id, nx, ny, ny_local, y_offset, z_layer_bottom, DomainSize,
                                     inputs.RNGSeed, ZMin, ZMaxLayer, deltax, temperature.NumberOfSolidificationEvents);

            // Initialize potential nucleation event data for next layer "layernumber + 1"
            // Views containing nucleation data will be resized to the possible number of nuclei on a given MPI rank for
            // the next layer
            nucleation.resetNucleiCounters(); // start counters at 0
            nucleation.placeNuclei(temperature, inputs.RNGSeed, layernumber + 1, nx, ny, nz_layer, ny_local, y_offset,
                                   z_layer_bottom, id, AtNorthBoundary, AtSouthBoundary);

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
    print.printFinalExaCAData(id, np, nx, ny, nz, ny_local, NumberOfLayers, DomainSize, cellData.LayerID_AllLayers,
                              cellData.CellType_AllLayers, cellData.GrainID_AllLayers, temperature, GrainUnitVector,
                              NGrainOrientations, deltax, XMin, YMin, ZMin);

    // Calculate volume fraction of solidified domain consisting of nucleated grains
    float VolFractionNucleated = calcVolFractionNucleated(id, nx, ny_local, DomainSize, cellData.LayerID_AllLayers,
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
    inputs.PrintExaCALog(id, np, InputFile, ny_local, y_offset, irf, deltax, HT_deltax, NumberOfLayers, LayerHeight, nx,
                         ny, nz, InitTime, RunTime, OutTime, cycle, InitMaxTime, InitMinTime, NuclMaxTime, NuclMinTime,
                         CreateSVMinTime, CreateSVMaxTime, CaptureMaxTime, CaptureMinTime, GhostMaxTime, GhostMinTime,
                         OutMaxTime, OutMinTime, XMin, XMax, YMin, YMax, ZMin, ZMax, VolFractionNucleated);
    if (id == 0)
        PrintExaCATiming(np, InitTime, RunTime, OutTime, cycle, InitMaxTime, InitMinTime, NuclMaxTime, NuclMinTime,
                         CreateSVMinTime, CreateSVMaxTime, CaptureMaxTime, CaptureMinTime, GhostMaxTime, GhostMinTime,
                         OutMaxTime, OutMinTime);
}
