// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_PRINT_HPP
#define EXACA_PRINT_HPP

#include "CAfunctions.hpp"
#include "CAinputs.hpp"
#include "CAinterfacialresponse.hpp"
#include "CAparsefiles.hpp"
#include "CAtemperature.hpp"
#include "CAtypes.hpp"
#include "mpi.h"

#include <Kokkos_Core.hpp>

#include <nlohmann/json.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Write data of type PrintType as ascii or binary, with option to convert between big and small endian binary
template <typename PrintType>
void WriteData(std::ofstream &outstream, PrintType PrintValue, bool PrintBinary, bool SwapEndianYN = false) {
    if (PrintBinary) {
        if (SwapEndianYN)
            SwapEndian(PrintValue);
        int varSize = sizeof(PrintType);
        outstream.write((char *)&PrintValue, varSize);
    }
    else
        outstream << PrintValue << " ";
}

// Struct to hold data printing options and functions
struct Print {

    // If printing the time series, the counter for the number of intermediate files that have been printed for a given
    // layer of a multilayer problem
    int IntermediateFileCounter = 0;

    // Message sizes and data offsets for data send/recieved to/from other ranks- message size different for different
    // ranks
    using view_type_int_host = Kokkos::View<int *, Kokkos::HostSpace>;
    view_type_int_host Recv_y_offset, Recv_ny_local, RBufSize;
    // Y coordinates for a given rank's data being send/loaded into the view of all domain data on rank 0=
    int SendBufStartY, SendBufEndY, SendBufSize;
    // Holds print options from input file
    PrintInputs _inputs;
    // Combined path/file prefix for output files
    std::string PathBaseFileName;

    // Default constructor - options are set in getPrintDataFromFile and copied into this struct
    Print(int nx, int ny, int nz, int ny_local, int y_offset, int np, PrintInputs inputs)
        : Recv_y_offset(view_type_int_host(Kokkos::ViewAllocateWithoutInitializing("Recv_y_offset"), np))
        , Recv_ny_local(view_type_int_host(Kokkos::ViewAllocateWithoutInitializing("Recv_ny_local"), np))
        , RBufSize(view_type_int_host(Kokkos::ViewAllocateWithoutInitializing("RBufSize"), np))
        , _inputs(inputs) {

        // Buffers for sending/receiving data across ranks
        for (int recvrank = 0; recvrank < np; recvrank++) {
            Recv_y_offset(recvrank) = get_yoffset(recvrank, ny, np);
            Recv_ny_local(recvrank) = get_nylocal(recvrank, ny, np);
            RBufSize(recvrank) = nx * Recv_ny_local(recvrank) * nz;
        }

        // Y coordinates for a given rank's data being send/loaded into the view of all domain data on rank 0=
        if (y_offset == 0)
            SendBufStartY = 0;
        else
            SendBufStartY = 1;
        if (ny_local + y_offset == ny)
            SendBufEndY = ny_local;
        else
            SendBufEndY = ny_local - 1;
        SendBufSize = nx * (SendBufEndY - SendBufStartY) * nz;

        PathBaseFileName = _inputs.PathToOutput + _inputs.BaseFileName;
    }

    // Called on rank 0 to collect view data from other ranks, or on other ranks to send data to rank 0
    // MPI datatype corresponding to the view data
    template <typename Collect1DViewTypeDevice>
    auto collectViewData(int id, int np, int nx, int ny, int nz, int ny_local, MPI_Datatype msg_type,
                         Collect1DViewTypeDevice ViewDataThisRank_Device) {
        // View (int or float extracted from 1D View) for 3D data (no initial size, only given size/filled on rank 0)
        using value_type = typename Collect1DViewTypeDevice::value_type;
        Kokkos::View<value_type ***, Kokkos::HostSpace> ViewData_WholeDomain(
            Kokkos::ViewAllocateWithoutInitializing(ViewDataThisRank_Device.label() + "_WholeDomain"), 0, 0, 0);

        // Get the host view type
        using host_view_type = typename Collect1DViewTypeDevice::HostMirror;

        // Copy view data of interest host
        auto ViewDataThisRank = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ViewDataThisRank_Device);
        if (id == 0) {
            // View (short, int, or float) for 3D data given a size
            Kokkos::resize(ViewData_WholeDomain, nz, nx, ny);
            // Place rank 0 data into view for whole domain
            for (int coord_z = 0; coord_z < nz; coord_z++) {
                for (int coord_x = 0; coord_x < nx; coord_x++) {
                    for (int coord_y_local = 0; coord_y_local < ny_local; coord_y_local++) {
                        int index = get1Dindex(coord_x, coord_y_local, coord_z, nx, ny_local);
                        ViewData_WholeDomain(coord_z, coord_x, coord_y_local) = ViewDataThisRank(index);
                    }
                }
            }

            // Recieve values from other ranks - message size different for different ranks
            for (int recvrank = 1; recvrank < np; recvrank++) {
                int RecvBufSize_ThisRank = RBufSize(recvrank);
                host_view_type RecvBufData(Kokkos::ViewAllocateWithoutInitializing("RecvBufData"),
                                           RecvBufSize_ThisRank);
                MPI_Recv(RecvBufData.data(), RecvBufSize_ThisRank, msg_type, recvrank, 0, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
                int DataCounter = 0;
                for (int coord_z = 0; coord_z < nz; coord_z++) {
                    for (int coord_x = 0; coord_x < nx; coord_x++) {
                        for (int coord_y_local = 0; coord_y_local < Recv_ny_local(recvrank); coord_y_local++) {
                            int coord_y_global = coord_y_local + Recv_y_offset(recvrank);
                            ViewData_WholeDomain(coord_z, coord_x, coord_y_global) = RecvBufData(DataCounter);
                            DataCounter++;
                        }
                    }
                }
            }
        }
        else {
            // Send non-ghost node data to rank 0
            int DataCounter = 0;
            host_view_type SendBuf(Kokkos::ViewAllocateWithoutInitializing("SendBuf"), SendBufSize);
            for (int coord_z = 0; coord_z < nz; coord_z++) {
                for (int coord_x = 0; coord_x < nx; coord_x++) {
                    for (int coord_y_local = SendBufStartY; coord_y_local < SendBufEndY; coord_y_local++) {
                        int index = get1Dindex(coord_x, coord_y_local, coord_z, nx, ny_local);
                        SendBuf(DataCounter) = ViewDataThisRank(index);
                        DataCounter++;
                    }
                }
            }
            MPI_Send(SendBuf.data(), SendBufSize, msg_type, 0, 0, MPI_COMM_WORLD);
        }
        return ViewData_WholeDomain;
    }

    // Called on rank 0, prints initial values of selected data structures to Paraview files for cells between the
    // overall simulation bottom and the top of the first layer
    template <typename ViewTypeGrainID, typename ViewTypeLayerID>
    void printInitExaCAData(int id, int np, int nx, int ny, int ny_local, int nz_layer, int DomainSize, double deltax,
                            double XMin, double YMin, double ZMin, ViewTypeGrainID GrainID, ViewTypeLayerID LayerID,
                            Temperature<device_memory_space> &temperature) {

        if ((_inputs.PrintInitGrainID) || (_inputs.PrintInitLayerID) || (_inputs.PrintInitMeltTimeStep) ||
            (_inputs.PrintInitCritTimeStep) || (_inputs.PrintInitUndercoolingChange)) {
            std::string FName = PathBaseFileName + "_init.vtk";
            std::ofstream Grainplot;
            if (id == 0) {
                std::cout << "Printing initial data structures to a vtk file" << std::endl;
                WriteHeader(Grainplot, FName, nx, ny, nz_layer, deltax, XMin, YMin, ZMin);
            }
            if (_inputs.PrintInitGrainID) {
                auto GrainID_WholeDomain = collectViewData(id, np, nx, ny, nz_layer, ny_local, MPI_INT, GrainID);
                printViewData(id, Grainplot, nx, ny, nz_layer, "int", "GrainID", GrainID_WholeDomain);
            }
            if (_inputs.PrintInitLayerID) {
                auto LayerID_WholeDomain = collectViewData(id, np, nx, ny, nz_layer, ny_local, MPI_SHORT, LayerID);
                printViewData(id, Grainplot, nx, ny, nz_layer, "short", "LayerID", LayerID_WholeDomain);
            }
            if (_inputs.PrintInitMeltTimeStep) {
                ViewI MeltTimeStep = temperature.extract_tm_tl_cr_data<ViewI>(0, DomainSize);
                auto MeltTimeStep_WholeDomain =
                    collectViewData(id, np, nx, ny, nz_layer, ny_local, MPI_INT, MeltTimeStep);
                printViewData(id, Grainplot, nx, ny, nz_layer, "int", "MeltTimeStep", MeltTimeStep_WholeDomain);
            }
            if (_inputs.PrintInitCritTimeStep) {
                ViewI CritTimeStep = temperature.extract_tm_tl_cr_data<ViewI>(1, DomainSize);
                auto CritTimeStep_WholeDomain =
                    collectViewData(id, np, nx, ny, nz_layer, ny_local, MPI_INT, CritTimeStep);
                printViewData(id, Grainplot, nx, ny, nz_layer, "int", "CritTimeStep", CritTimeStep_WholeDomain);
            }
            if (_inputs.PrintInitUndercoolingChange) {
                ViewF UndercoolingChange = temperature.extract_tm_tl_cr_data<ViewF>(2, DomainSize, 0);
                auto UndercoolingChange_WholeDomain =
                    collectViewData(id, np, nx, ny, nz_layer, ny_local, MPI_FLOAT, UndercoolingChange);
                printViewData(id, Grainplot, nx, ny, nz_layer, "float", "UndercoolingChange",
                              UndercoolingChange_WholeDomain);
            }
            if (id == 0)
                Grainplot.close();
        }
    }

    // Called on rank 0, prints intermediate values of grain misorientation for all layers up to current layer, also
    // marking which cells are liquid
    template <typename ViewTypeGrainID, typename ViewTypeCellType, typename ViewTypeGrainUnitVector>
    void printIntermediateGrainMisorientation(int id, int np, int cycle, int nx, int ny, int nz, int ny_local,
                                              int nz_layer, double deltax, double XMin, double YMin, double ZMin,
                                              ViewTypeGrainID GrainID, ViewTypeCellType CellType,
                                              ViewTypeGrainUnitVector GrainUnitVector, int NGrainOrientations,
                                              int layernumber, int z_layer_bottom) {

        // Check if option is toggled and whether the output should be printed this time step
        if ((_inputs.PrintTimeSeries) && (cycle % _inputs.TimeSeriesInc == 0)) {
            IntermediateFileCounter++;
            auto GrainID_WholeDomain = collectViewData(id, np, nx, ny, nz, ny_local, MPI_INT, GrainID);
            auto CellType_WholeDomain = collectViewData(id, np, nx, ny, nz, ny_local, MPI_INT, CellType);
            if (id == 0) {
                std::cout << "Intermediate output on time step " << cycle << std::endl;
                auto GrainUnitVector_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainUnitVector);
                printGrainMisorientations(nx, ny, nz, GrainID_WholeDomain, CellType_WholeDomain, GrainUnitVector_Host,
                                          NGrainOrientations, deltax, XMin, YMin, ZMin, true, layernumber,
                                          z_layer_bottom, nz_layer);
            }
        }
    }

    // Check if intermediate output is to be printed during a series of "skipped" time steps (i.e., time steps where no
    // melting or solidification occurs). If it should be printed on the time step, call
    // printIntermediateGrainMisorientation
    template <typename ViewTypeGrainID, typename ViewTypeCellType, typename ViewTypeGrainUnitVector>
    void printIdleIntermediateGrainMisorientation(int id, int np, int cycle, int nx, int ny, int nz, int ny_local,
                                                  int nz_layer, double deltax, double XMin, double YMin, double ZMin,
                                                  ViewTypeGrainID GrainID, ViewTypeCellType CellType,
                                                  ViewTypeGrainUnitVector GrainUnitVector, int NGrainOrientations,
                                                  int layernumber, int z_layer_bottom,
                                                  unsigned long int GlobalNextMeltTimeStep) {

        if (_inputs.PrintIdleTimeSeriesFrames) {
            // Print current state of ExaCA simulation (up to and including the current layer's data) during the skipped
            // time steps, if intermediate output is toggled
            for (unsigned long int cycle_jump = cycle + 1; cycle_jump < GlobalNextMeltTimeStep; cycle_jump++) {
                printIntermediateGrainMisorientation(id, np, cycle_jump, nx, ny, nz, ny_local, nz_layer, deltax, XMin,
                                                     YMin, ZMin, GrainID, CellType, GrainUnitVector, NGrainOrientations,
                                                     layernumber, z_layer_bottom);
            }
        }
    }

    // Prints final values of selected data structures to Paraview files
    template <typename ViewTypeGrainID, typename ViewTypeLayerID, typename ViewTypeCell, typename ViewTypeGrainUnit>
    void printFinalExaCAData(int id, int np, int nx, int ny, int nz, int ny_local, int NumberOfLayers, int DomainSize,
                             ViewTypeLayerID LayerID, ViewTypeCell CellType, ViewTypeGrainID GrainID,
                             Temperature<device_memory_space> &temperature, ViewTypeGrainUnit GrainUnitVector,
                             int NGrainOrientations, double deltax, double XMin, double YMin, double ZMin) {

        if (id == 0)
            std::cout << "Printing final data structures to vtk files" << std::endl;

        if ((_inputs.PrintFinalGrainID) || (_inputs.PrintFinalLayerID) || (_inputs.PrintFinalMisorientation) ||
            (_inputs.PrintDefaultRVE)) {

            // GrainID and LayerID are printed to one vtk file, while Grain misorientation is printed to a separate vtk
            // file GrainID and LayerID are also needed needed if PrintFinalMisorientation = true
            std::string FName = PathBaseFileName + ".vtk";
            std::ofstream Grainplot;
            if (id == 0)
                WriteHeader(Grainplot, FName, nx, ny, nz, deltax, XMin, YMin, ZMin);
            if ((_inputs.PrintFinalGrainID) || (_inputs.PrintFinalLayerID) || (_inputs.PrintFinalMisorientation) ||
                (_inputs.PrintDefaultRVE)) {
                auto GrainID_WholeDomain = collectViewData(id, np, nx, ny, nz, ny_local, MPI_INT, GrainID);
                auto LayerID_WholeDomain = collectViewData(id, np, nx, ny, nz, ny_local, MPI_SHORT, LayerID);
                if (_inputs.PrintFinalGrainID)
                    printViewData(id, Grainplot, nx, ny, nz, "int", "GrainID", GrainID_WholeDomain);
                if (_inputs.PrintFinalLayerID)
                    printViewData(id, Grainplot, nx, ny, nz, "short", "LayerID", LayerID_WholeDomain);
                if ((id == 0) && _inputs.PrintFinalMisorientation) {
                    auto GrainUnitVector_Host =
                        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainUnitVector);
                    // empty view passed as dummy argument (not used in this call to printGrainMisorientation)
                    Kokkos::View<int ***, Kokkos::HostSpace> CellType_WholeDomain(
                        Kokkos::ViewAllocateWithoutInitializing("ViewData_WholeDomain"), 0, 0, 0);
                    printGrainMisorientations(nx, ny, nz, GrainID_WholeDomain, CellType_WholeDomain,
                                              GrainUnitVector_Host, NGrainOrientations, deltax, XMin, YMin, ZMin,
                                              false);
                }
                if ((id == 0) && (_inputs.PrintDefaultRVE))
                    printExaConstitDefaultRVE(nx, ny, nz, LayerID_WholeDomain, GrainID_WholeDomain, deltax,
                                              NumberOfLayers);
            }
            if (id == 0)
                Grainplot.close();
        }

        if ((_inputs.PrintFinalMeltTimeStep) || (_inputs.PrintFinalCritTimeStep) ||
            (_inputs.PrintFinalUndercoolingChange) || (_inputs.PrintFinalCellType) ||
            (_inputs.PrintFinalUndercoolingCurrent)) {
            // Temperature field data is printed to a separate vtk file
            std::string FName = PathBaseFileName + "_final.vtk";
            std::ofstream GrainplotF;
            if (id == 0)
                WriteHeader(GrainplotF, FName, nx, ny, nz, deltax, XMin, YMin, ZMin);
            if (_inputs.PrintFinalMeltTimeStep) {
                ViewI MeltTimeStep = temperature.extract_tm_tl_cr_data<ViewI>(0, DomainSize);
                auto MeltTimeStep_WholeDomain = collectViewData(id, np, nx, ny, nz, ny_local, MPI_INT, MeltTimeStep);
                printViewData(id, GrainplotF, nx, ny, nz, "int", "MeltTimeStep", MeltTimeStep_WholeDomain);
            }
            if (_inputs.PrintFinalCritTimeStep) {
                ViewI CritTimeStep = temperature.extract_tm_tl_cr_data<ViewI>(1, DomainSize);
                auto CritTimeStep_WholeDomain = collectViewData(id, np, nx, ny, nz, ny_local, MPI_INT, CritTimeStep);
                printViewData(id, GrainplotF, nx, ny, nz, "int", "CritTimeStep", CritTimeStep_WholeDomain);
            }
            if (_inputs.PrintFinalUndercoolingChange) {
                ViewF UndercoolingChange = temperature.extract_tm_tl_cr_data<ViewF>(0, DomainSize, 0);
                auto UndercoolingChange_WholeDomain =
                    collectViewData(id, np, nx, ny, nz, ny_local, MPI_FLOAT, UndercoolingChange);
                printViewData(id, GrainplotF, nx, ny, nz, "int", "UndercoolingChange", UndercoolingChange_WholeDomain);
            }
            if (_inputs.PrintFinalUndercoolingCurrent) {
                auto UndercoolingCurrent_WholeDomain =
                    collectViewData(id, np, nx, ny, nz, ny_local, MPI_FLOAT, temperature.UndercoolingCurrent);
                printViewData(id, GrainplotF, nx, ny, nz, "float", "UndercoolingFinal",
                              UndercoolingCurrent_WholeDomain);
            }
            if (_inputs.PrintFinalCellType) {
                auto CellType_WholeDomain = collectViewData(id, np, nx, ny, nz, ny_local, MPI_INT, CellType);
                printViewData(id, GrainplotF, nx, ny, nz, "int", "CellType", CellType_WholeDomain);
            }
            if (id == 0)
                GrainplotF.close();
        }
    }

    // Print Paraview file header data on rank 0, for either binary or ASCII output
    void WriteHeader(std::ofstream &ParaviewOutputStream, std::string FName, int nx, int ny, int nz, double deltax,
                     double XMin, double YMin, double ZMin) {

        if (_inputs.PrintBinary)
            ParaviewOutputStream.open(FName, std::ios::out | std::ios::binary);
        else
            ParaviewOutputStream.open(FName);
        ParaviewOutputStream << "# vtk DataFile Version 3.0" << std::endl;
        ParaviewOutputStream << "vtk output" << std::endl;
        if (_inputs.PrintBinary)
            ParaviewOutputStream << "BINARY" << std::endl;
        else
            ParaviewOutputStream << "ASCII" << std::endl;
        ParaviewOutputStream << "DATASET STRUCTURED_POINTS" << std::endl;
        ParaviewOutputStream << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;
        ParaviewOutputStream << "ORIGIN " << XMin << " " << YMin << " " << ZMin << std::endl;
        ParaviewOutputStream << "SPACING " << deltax << " " << deltax << " " << deltax << std::endl;
        ParaviewOutputStream << std::fixed << "POINT_DATA " << nx * ny * nz << std::endl;
    }

    // Called on rank 0 to write view data to the vtk file
    template <typename Print3DViewType>
    void printViewData(int id, std::ofstream &Grainplot, int nx, int ny, int nz, std::string DataLabel,
                       std::string VarNameLabel, Print3DViewType ViewData_WholeDomain) {
        if (id != 0)
            return;

        // Print data to the vtk file - casting to the appropriate type if necessary
        Grainplot << "SCALARS " << VarNameLabel << " " << DataLabel << " 1" << std::endl;
        Grainplot << "LOOKUP_TABLE default" << std::endl;
        for (int coord_z = 0; coord_z < nz; coord_z++) {
            for (int coord_y_global = 0; coord_y_global < ny; coord_y_global++) {
                for (int coord_x = 0; coord_x < nx; coord_x++) {
                    if (DataLabel == "int") {
                        int writeval = static_cast<int>(ViewData_WholeDomain(coord_z, coord_x, coord_y_global));
                        WriteData(Grainplot, writeval, _inputs.PrintBinary, true);
                    }
                    else if (DataLabel == "short") {
                        short writeval = static_cast<short>(ViewData_WholeDomain(coord_z, coord_x, coord_y_global));
                        WriteData(Grainplot, writeval, _inputs.PrintBinary, true);
                    }
                    else if (DataLabel == "float") {
                        float writeval = static_cast<float>(ViewData_WholeDomain(coord_z, coord_x, coord_y_global));
                        WriteData(Grainplot, writeval, _inputs.PrintBinary, true);
                    }
                }
            }
            // Do not insert newline character if using binary writing, as this will break the binary data read by
            // adding a blank line
            if (!(_inputs.PrintBinary))
                Grainplot << std::endl;
        }
    }

    // Reset the intermediate file counter to 0, performed at the start of a new layer when intermediate output is being
    // printed
    void resetIntermediateFileCounter() {
        if (_inputs.PrintTimeSeries)
            IntermediateFileCounter = 0;
    }

    // On rank 0, print grain misorientation, 0-62 for epitaxial grains and 100-162 for nucleated grains, to a paraview
    // file Optionally add layer label/intermediate frame and print -1 for liquid cells if this is being printed as
    // intermediate output
    template <typename ViewTypeInt3DHost, typename ViewTypeFloat>
    void printGrainMisorientations(int nx, int ny, int nz, ViewTypeInt3DHost GrainID_WholeDomain,
                                   ViewTypeInt3DHost CellType_WholeDomain, ViewTypeFloat GrainUnitVector,
                                   int NGrainOrientations, double deltax, double XMin, double YMin, double ZMin,
                                   bool IntermediatePrint, int layernumber = 0, int ZBound_Low = 0, int nz_layer = 0) {

        // Print grain orientations to file - either all layers, or if in an intermediate state, the layers up to the
        // current one
        std::string FName;
        int ZPrintSize;
        if (IntermediatePrint) {
            FName = PathBaseFileName + "_layer" + std::to_string(layernumber) + "_" +
                    std::to_string(IntermediateFileCounter) + ".vtk";
            ZPrintSize = ZBound_Low + nz_layer;
        }
        else {
            FName = PathBaseFileName + "_Misorientations.vtk";
            ZPrintSize = nz;
        }
        std::cout << "Printing Paraview file of grain misorientations for Z coordinates of 0 through " << ZPrintSize - 1
                  << std::endl;
        std::ofstream GrainplotM;
        WriteHeader(GrainplotM, FName, nx, ny, ZPrintSize, deltax, XMin, YMin, ZMin);
        GrainplotM << "SCALARS Angle_z short 1" << std::endl;
        GrainplotM << "LOOKUP_TABLE default" << std::endl;

        // Get grain <100> misorientation relative to the Z direction for each orientation
        auto GrainMisorientation = MisorientationCalc(NGrainOrientations, GrainUnitVector, 2);
        // For cells that are currently liquid (possible only for intermediate state print, as the final state will only
        // have solid cells), -1 is printed as the misorienatation. Misorientations for grains from the baseplate or
        // powder layer are between 0-62 (degrees, rounded to nearest integer), and cells that are associated with
        // nucleated grains are assigned values between 100-162 to differentiate them. Additionally, 200 is printed as
        // the misorientation for cells in the powder layer that have not been assigned a grain ID.
        for (int k = 0; k < ZPrintSize; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    short IntPrintVal;
                    if (GrainID_WholeDomain(k, i, j) == 0)
                        IntPrintVal = 200;
                    else {
                        int MyOrientation = getGrainOrientation(GrainID_WholeDomain(k, i, j), NGrainOrientations);
                        if (GrainID_WholeDomain(k, i, j) < 0)
                            IntPrintVal = static_cast<short>(std::round(GrainMisorientation(MyOrientation)) + 100);
                        else
                            IntPrintVal = static_cast<short>(std::round(GrainMisorientation(MyOrientation)));
                    }
                    if (IntermediatePrint) {
                        if (CellType_WholeDomain(k, i, j) == Liquid)
                            IntPrintVal = -1;
                    }
                    WriteData(GrainplotM, IntPrintVal, _inputs.PrintBinary, true);
                }
            }
            if (!(_inputs.PrintBinary))
                GrainplotM << std::endl;
        }
        GrainplotM.close();
    }

    // Print a representative volume element (RVE) from this multilayer simulation from the "default" location in the
    // domain. The default location is as close to the center of the domain in X and Y as possible, and as close to the
    // top of the domain while not including the final layer's microstructure. If an RVE size was not specified in the
    // input file, the default size is 0.5 by 0.5 by 0.5 mm
    template <typename ViewTypeShort3DHost, typename ViewTypeInt3DHost>
    void printExaConstitDefaultRVE(int nx, int ny, int nz, ViewTypeShort3DHost LayerID_WholeDomain,
                                   ViewTypeInt3DHost GrainID_WholeDomain, double deltax, int NumberOfLayers) {

        // Determine the lower and upper Y bounds of the RVE
        int RVE_XLow = std::floor(nx / 2) - std::floor(_inputs.RVESize / 2);
        int RVE_XHigh = RVE_XLow + _inputs.RVESize - 1;
        int RVE_YLow = std::floor(ny / 2) - std::floor(_inputs.RVESize / 2);
        int RVE_YHigh = RVE_YLow + _inputs.RVESize - 1;

        // Make sure the RVE fits in the simulation domain in X and Y
        if ((RVE_XLow < 0) || (RVE_XHigh > nx - 1) || (RVE_YLow < 0) || (RVE_YHigh > ny - 1)) {
            std::cout << "WARNING: Simulation domain is too small to obtain default RVE data (should be at least "
                      << _inputs.RVESize << " cells in the X and Y directions" << std::endl;
            if (RVE_XLow < 0)
                RVE_XLow = 0;
            if (RVE_XHigh > nx - 1)
                RVE_XHigh = nx - 1;
            if (RVE_YLow < 0)
                RVE_YLow = 0;
            if (RVE_YHigh > ny - 1)
                RVE_YHigh = ny - 1;
        }

        // Determine the upper Z bound of the RVE - largest Z for which all cells in the RVE do not contain LayerID
        // values of the last layer
        int RVE_ZHigh = nz - 1;
        for (int k = nz - 1; k >= 0; k--) {
            [&] {
                for (int j = RVE_YLow; j <= RVE_YHigh; j++) {
                    for (int i = RVE_XLow; i <= RVE_XHigh; i++) {
                        if (LayerID_WholeDomain(k, i, j) == (NumberOfLayers - 1))
                            return; // check next lowest value for k
                    }
                }
                RVE_ZHigh = k;
                k = 0; // leave loop
            }();
        }

        // Determine the lower Z bound of the RVE, and make sure the RVE fits in the simulation domain in X and Y
        int RVE_ZLow = RVE_ZHigh - _inputs.RVESize + 1;
        if (RVE_ZLow < 0) {
            std::cout << "WARNING: Simulation domain is too small to obtain default RVE data (should be at least "
                      << _inputs.RVESize << " cells in the Z direction, more layers are required" << std::endl;
            RVE_ZLow = 0;
            RVE_ZHigh = nz - 1;
        }

        // Print RVE data to file
        std::string FName = PathBaseFileName + "_ExaConstit.csv";
        std::cout << "Default size RVE with X coordinates " << RVE_XLow << "," << RVE_XHigh << "; Y coordinates "
                  << RVE_YLow << "," << RVE_YHigh << "; Z coordinates " << RVE_ZLow << "," << RVE_ZHigh
                  << " being printed to file " << FName << " for ExaConstit" << std::endl;
        std::ofstream GrainplotE;
        GrainplotE.open(FName);
        GrainplotE << "Coordinates are in CA units, 1 cell = " << deltax << " m. Data is cell-centered. Origin at "
                   << RVE_XLow << "," << RVE_YLow << "," << RVE_ZLow << " , domain size is " << RVE_XHigh - RVE_XLow + 1
                   << " by " << RVE_YHigh - RVE_YLow + 1 << " by " << RVE_ZHigh - RVE_ZLow + 1 << " cells" << std::endl;
        GrainplotE << "X coord, Y coord, Z coord, Grain ID" << std::endl;
        for (int k = RVE_ZLow; k <= RVE_ZHigh; k++) {
            for (int j = RVE_YLow; j <= RVE_YHigh; j++) {
                for (int i = RVE_XLow; i <= RVE_XHigh; i++) {
                    GrainplotE << i << "," << j << "," << k << "," << GrainID_WholeDomain(k, i, j) << std::endl;
                }
            }
        }
        GrainplotE.close();
    }
};

void PrintExaCATiming(int np, double InitTime, double RunTime, double OutTime, int cycle, double InitMaxTime,
                      double InitMinTime, double NuclMaxTime, double NuclMinTime, double CreateSVMinTime,
                      double CreateSVMaxTime, double CaptureMaxTime, double CaptureMinTime, double GhostMaxTime,
                      double GhostMinTime, double OutMaxTime, double OutMinTime);
#endif
