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
    std::string simulation_type = inputs.SimulationType;

    //  Final time step at which solidification is expected TODO: may not be needed
    int *FinishTimeStep = new int[inputs.domain.NumberOfLayers];

    // Setup local and global grids, decomposing domain
    Grid grid(simulation_type, id, np, inputs.domain.NumberOfLayers, inputs.domain, inputs.temperature);

    // Material response function
    InterfacialResponseFunction irf(id, inputs.MaterialFileName, inputs.domain.deltat, grid.deltax);

    // Ensure that input powder layer init options are compatible with this domain size, if needed for this problem type
    // TODO: Move to input struct and expand to check that inputs are valid for the problem type
    if ((simulation_type == "R") || (simulation_type == "S"))
        checkPowderOverflow(grid.nx, grid.ny, grid.LayerHeight, grid.NumberOfLayers, inputs);

    // Temperature fields characterized by data in this structure
    Temperature<device_memory_space> temperature(grid.DomainSize, grid.NumberOfLayers, inputs.temperature);
    // Read temperature data if necessary
    if (simulation_type == "R")
        temperature.readTemperatureData(id, grid, 0);
    // Initialize the temperature fields for the simualtion type of interest
    if ((simulation_type == "C") || (simulation_type == "SingleGrain")) {
        if (inputs.temperature.G == 0)
            temperature.initialize(id, grid, inputs.domain.deltat);
        else
            temperature.initialize(id, simulation_type, grid, inputs.domain.deltat);
    }
    else if (simulation_type == "S")
        temperature.initialize(id, grid, irf.FreezingRange, inputs);
    else if (simulation_type == "R")
        temperature.initialize(0, id, grid, irf.FreezingRange, FinishTimeStep, inputs.domain.deltat);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Done with temperature field initialization" << std::endl;

    int NGrainOrientations = 0; // Number of grain orientations considered in the simulation
    // No initialize size yet, will be resized in OrientationInit
    ViewF GrainUnitVector(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), 0);

    // Initialize grain orientations
    // TODO: Orientation data to be part of a new struct
    OrientationInit(id, NGrainOrientations, GrainUnitVector, inputs.GrainOrientationFile);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Done with orientation initialization " << std::endl;

    // Variables characterizing the active cell region within each rank's grid
    // DiagonalLengths should be initialized to 0, DOCenter/CritDiagonalLength do not need initialization
    // TODO: Interface parameters/views (active cells and buffers) to be part of a new struct
    ViewF DiagonalLength("DiagonalLength", grid.DomainSize);
    ViewF DOCenter(Kokkos::ViewAllocateWithoutInitializing("DOCenter"), 3 * grid.DomainSize);
    ViewF CritDiagonalLength(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength"), 26 * grid.DomainSize);

    // Buffers for ghost node data (fixed size)
    int BufSizeInitialEstimate = 25;
    int BufSize = BufSizeInitialEstimate; // set to initial estimate
    int BufComponents = 8;
    Buffer2D BufferSouthSend(Kokkos::ViewAllocateWithoutInitializing("BufferSouthSend"), BufSize, BufComponents);
    Buffer2D BufferNorthSend(Kokkos::ViewAllocateWithoutInitializing("BufferNorthSend"), BufSize, BufComponents);
    Buffer2D BufferSouthRecv(Kokkos::ViewAllocateWithoutInitializing("BufferSouthRecv"), BufSize, BufComponents);
    Buffer2D BufferNorthRecv(Kokkos::ViewAllocateWithoutInitializing("BufferNorthRecv"), BufSize, BufComponents);
    ViewI SendSizeSouth(Kokkos::ViewAllocateWithoutInitializing("SendSizeSouth"), 1);
    ViewI SendSizeNorth(Kokkos::ViewAllocateWithoutInitializing("SendSizeNorth"), 1);
    ViewI_H SendSizeSouth_Host(Kokkos::ViewAllocateWithoutInitializing("SendSizeSouth_Host"), 1);
    ViewI_H SendSizeNorth_Host(Kokkos::ViewAllocateWithoutInitializing("SendSizeNorth_Host"), 1);
    // Send/recv buffers for ghost node data should be initialized with -1s in the first index as placeholders for empty
    // positions in the buffer, and with send size counts of 0
    ResetSendBuffers(BufSize, BufferNorthSend, BufferSouthSend, SendSizeNorth, SendSizeSouth);

    // Initialize cell types, grain IDs, and layer IDs
    CellData<device_memory_space> cellData(grid.DomainSize_AllLayers, inputs.substrate);
    if (simulation_type == "C")
        cellData.init_substrate(id, grid, inputs.RNGSeed);
    else if (simulation_type == "SingleGrain")
        cellData.init_substrate(id, grid);
    else
        cellData.init_substrate(id, grid, inputs.RNGSeed, temperature.NumberOfSolidificationEvents);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0)
        std::cout << "Grain struct initialized" << std::endl;

    // Nucleation data structure, containing views of nuclei locations, time steps, and ids, and nucleation event
    // counters - initialized with an estimate on the number of nuclei in the layer Without knowing
    // PossibleNuclei_ThisRankThisLayer yet, initialize nucleation data structures to estimated sizes, resize inside of
    // NucleiInit when the number of nuclei per rank is known
    int EstimatedNuclei_ThisRankThisLayer = inputs.nucleation.NMax * pow(grid.deltax, 3) * grid.DomainSize;
    Nucleation<device_memory_space> nucleation(EstimatedNuclei_ThisRankThisLayer, grid.deltax, inputs.nucleation);
    // Fill in nucleation data structures, and assign nucleation undercooling values to potential nucleation events
    // Potential nucleation grains are only associated with liquid cells in layer 0 - they will be initialized for each
    // successive layer when layer 0 in complete
    nucleation.placeNuclei(temperature, inputs.RNGSeed, 0, grid, id);

    // Steering Vector TODO: part of interface struct in the future
    ViewI SteeringVector(Kokkos::ViewAllocateWithoutInitializing("SteeringVector"), grid.DomainSize);
    ViewI_H numSteer_Host(Kokkos::ViewAllocateWithoutInitializing("SteeringVectorSize"), 1);
    numSteer_Host(0) = 0;
    ViewI numSteer = Kokkos::create_mirror_view_and_copy(device_memory_space(), numSteer_Host);

    // Initialize printing struct from inputs
    Print print(grid, np, inputs.print);

    // If specified, print initial values in some views for debugging purposes
    double InitTime = MPI_Wtime() - StartInitTime;
    if (id == 0)
        std::cout << "Data initialized: Time spent: " << InitTime << " s" << std::endl;
    print.printInitExaCAData(id, np, grid, cellData.GrainID_AllLayers, cellData.LayerID_AllLayers, temperature);
    MPI_Barrier(MPI_COMM_WORLD);
    int cycle = 0;
    double StartRunTime = MPI_Wtime();
    for (int layernumber = 0; layernumber < grid.NumberOfLayers; layernumber++) {

        int XSwitch = 0;
        double LayerTime1 = MPI_Wtime();

        // Loop continues until all liquid cells claimed by solid grains
        do {

            // Start of time step - optional print current state of ExaCA simulation (up to and including the current
            // layer's data)
            print.printIntermediateGrainMisorientation(id, np, cycle, grid, cellData.GrainID_AllLayers,
                                                       cellData.CellType_AllLayers, GrainUnitVector, NGrainOrientations,
                                                       layernumber);
            cycle++;

            // Update cells on GPU - undercooling and diagonal length updates, nucleation
            // Cells with a successful nucleation event are marked and added to a steering vector, later dealt with in
            // CellCapture
            StartNuclTime = MPI_Wtime();
            nucleation.nucleate_grain(cycle, grid, cellData, SteeringVector, numSteer);
            NuclTime += MPI_Wtime() - StartNuclTime;

            // Update cells on GPU - new active cells, solidification of old active cells
            // Cell capture performed in two steps - first, adding cells of interest to a steering vector (different
            // subroutine called with versus without remelting), and second, iterating over the steering vector to
            // perform active cell creation and capture operations
            // Constrained/directional solidification problem cells only undergo 1 solidificaton event and are
            // initialized as liquid - the steering vector operation for this problem can be constructed using
            // FillSteeringVector_NoRemelt (a simplified version of FillSteeringVector_Remelt
            StartCreateSVTime = MPI_Wtime();

            if ((simulation_type == "C") || (simulation_type == "SingleGrain"))
                FillSteeringVector_NoRemelt(cycle, grid, temperature, cellData, SteeringVector, numSteer,
                                            numSteer_Host);
            else
                FillSteeringVector_Remelt(cycle, grid, temperature, cellData, SteeringVector, numSteer, numSteer_Host);
            CreateSVTime += MPI_Wtime() - StartCreateSVTime;

            StartCaptureTime = MPI_Wtime();
            CellCapture(id, np, cycle, grid, irf, GrainUnitVector, CritDiagonalLength, DiagonalLength, cellData,
                        temperature, DOCenter, NGrainOrientations, BufferNorthSend, BufferSouthSend, SendSizeNorth,
                        SendSizeSouth, SteeringVector, numSteer, numSteer_Host, BufSize);
            // Count the number of cells' in halo regions where the data did not fit into the send buffers
            // Reduce across all ranks, as the same BufSize should be maintained across all ranks
            // If any rank overflowed its buffer size, resize all buffers to the new size plus 10% padding
            int OldBufSize = BufSize;
            BufSize = ResizeBuffers(BufferNorthSend, BufferSouthSend, BufferNorthRecv, BufferSouthRecv, SendSizeNorth,
                                    SendSizeSouth, SendSizeNorth_Host, SendSizeSouth_Host, OldBufSize, BufComponents);
            if (OldBufSize != BufSize) {
                if (id == 0)
                    std::cout << "Resized number of cells stored in send/recv buffers from " << OldBufSize << " to "
                              << BufSize << std::endl;
                RefillBuffers(grid, cellData, BufferNorthSend, BufferSouthSend, SendSizeNorth, SendSizeSouth, DOCenter,
                              DiagonalLength, NGrainOrientations, BufSize);
            }
            CaptureTime += MPI_Wtime() - StartCaptureTime;

            if (np > 1) {
                // Update ghost nodes
                StartGhostTime = MPI_Wtime();
                GhostNodes1D(cycle, id, grid, cellData, DOCenter, GrainUnitVector, DiagonalLength, CritDiagonalLength,
                             NGrainOrientations, BufferNorthSend, BufferSouthSend, BufferNorthRecv, BufferSouthRecv,
                             BufSize, SendSizeNorth, SendSizeSouth, BufComponents);
                GhostTime += MPI_Wtime() - StartGhostTime;
            }

            if ((cycle % 1000 == 0) && (simulation_type != "SingleGrain")) {
                IntermediateOutputAndCheck(id, np, cycle, grid, nucleation.SuccessfulNucleationCounter, XSwitch,
                                           cellData, temperature, inputs.SimulationType, layernumber,
                                           NGrainOrientations, GrainUnitVector, print);
            }
            else if (simulation_type == "SingleGrain") {
                IntermediateOutputAndCheck(id, cycle, grid, XSwitch, cellData.CellType_AllLayers);
            }

        } while (XSwitch == 0);
        if (layernumber != grid.NumberOfLayers - 1) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (id == 0)
                std::cout << "Layer " << layernumber << " finished solidification; initializing layer "
                          << layernumber + 1 << std::endl;

            // Reset intermediate file counter to zero if printing video files
            print.resetIntermediateFileCounter();

            // Determine new active cell domain size and offset from bottom of global domain
            grid.init_next_layer(id, simulation_type, layernumber + 1, inputs.domain.SpotRadius);

            // For simulation type R, need to initialize new temperature field data for layer "layernumber + 1"
            // TODO: reorganize these temperature functions calls into a temperature.init_next_layer as done with the
            // substrate
            if (simulation_type == "R") {
                // If the next layer's temperature data isn't already stored, it should be read
                if (inputs.temperature.LayerwiseTempRead)
                    temperature.readTemperatureData(id, grid, layernumber + 1);
                // Initialize next layer's temperature data
                temperature.initialize(layernumber + 1, id, grid, irf.FreezingRange, FinishTimeStep,
                                       inputs.domain.deltat);
            }

            // Reset initial undercooling/solidification event counter of all cells to zeros for the next layer,
            // resizing the views to the number of cells associated with the next layer
            temperature.reset_layer_events_undercooling(grid.DomainSize);

            // Resize and zero all view data relating to the active region from the last layer, in preparation for the
            // next layer
            // TODO: part of the interface struct in the future
            ZeroResetViews(grid.DomainSize, DiagonalLength, CritDiagonalLength, DOCenter, SteeringVector);
            MPI_Barrier(MPI_COMM_WORLD);

            // Sets up views, powder layer (if necessary), and cell types for the next layer of a multilayer problem
            cellData.init_next_layer(layernumber + 1, id, grid, inputs.RNGSeed,
                                     temperature.NumberOfSolidificationEvents);

            // Initialize potential nucleation event data for next layer "layernumber + 1"
            // Views containing nucleation data will be resized to the possible number of nuclei on a given MPI rank for
            // the next layer
            nucleation.resetNucleiCounters(); // start counters at 0
            nucleation.placeNuclei(temperature, inputs.RNGSeed, layernumber + 1, grid, id);

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
    print.printFinalExaCAData(id, np, grid, cellData.LayerID_AllLayers, cellData.CellType_AllLayers,
                              cellData.GrainID_AllLayers, temperature, GrainUnitVector, NGrainOrientations);

    // Calculate volume fraction of solidified domain consisting of nucleated grains
    float VolFractionNucleated = cellData.calcVolFractionNucleated(id, grid);

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
    inputs.PrintExaCALog(id, np, InputFile, grid.ny_local, grid.y_offset, irf, grid.deltax, grid.HT_deltax,
                         grid.NumberOfLayers, grid.LayerHeight, grid.nx, grid.ny, grid.nz, InitTime, RunTime, OutTime,
                         cycle, InitMaxTime, InitMinTime, NuclMaxTime, NuclMinTime, CreateSVMinTime, CreateSVMaxTime,
                         CaptureMaxTime, CaptureMinTime, GhostMaxTime, GhostMinTime, OutMaxTime, OutMinTime, grid.XMin,
                         grid.XMax, grid.YMin, grid.YMax, grid.ZMin, grid.ZMax, VolFractionNucleated);
    if (id == 0)
        PrintExaCATiming(np, InitTime, RunTime, OutTime, cycle, InitMaxTime, InitMinTime, NuclMaxTime, NuclMinTime,
                         CreateSVMinTime, CreateSVMaxTime, CaptureMaxTime, CaptureMinTime, GhostMaxTime, GhostMinTime,
                         OutMaxTime, OutMinTime);
}
