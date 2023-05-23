// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAprint.hpp"
#include "CAconfig.hpp"
#include "CAfunctions.hpp"
#include "CAparsefiles.hpp"

#include "mpi.h"

// Print Paraview file header data, for either binary or ASCII output
void WriteHeader(std::ofstream &ParaviewOutputStream, std::string FName, bool PrintBinary, int nx, int ny, int nz,
                 double deltax, double XMin, double YMin, double ZMin) {

    if (PrintBinary)
        ParaviewOutputStream.open(FName, std::ios::out | std::ios::binary);
    else
        ParaviewOutputStream.open(FName);
    ParaviewOutputStream << "# vtk DataFile Version 3.0" << std::endl;
    ParaviewOutputStream << "vtk output" << std::endl;
    if (PrintBinary)
        ParaviewOutputStream << "BINARY" << std::endl;
    else
        ParaviewOutputStream << "ASCII" << std::endl;
    ParaviewOutputStream << "DATASET STRUCTURED_POINTS" << std::endl;
    ParaviewOutputStream << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;
    ParaviewOutputStream << "ORIGIN " << XMin << " " << YMin << " " << ZMin << std::endl;
    ParaviewOutputStream << "SPACING " << deltax << " " << deltax << " " << deltax << std::endl;
    ParaviewOutputStream << std::fixed << "POINT_DATA " << nx * ny * nz << std::endl;
}
//*****************************************************************************/
// On rank 0, collect data for one single int view
void CollectIntField(ViewI3D_H IntVar_WholeDomain, ViewI_H IntVar, int nz, int nx, int MyYSlices, int np,
                     ViewI_H RecvYOffset, ViewI_H RecvYSlices, ViewI_H RBufSize) {

    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                IntVar_WholeDomain(k, i, j) = IntVar(k * nx * MyYSlices + i * MyYSlices + j);
            }
        }
    }

    // Recieve values from other ranks - message size different for different ranks
    for (int p = 1; p < np; p++) {
        int RecvBufSize_ThisRank = RBufSize(p);
        ViewI_H RecvBufIntVar(Kokkos::ViewAllocateWithoutInitializing("RecvBufIntVar"), RecvBufSize_ThisRank);
        MPI_Recv(RecvBufIntVar.data(), RecvBufSize_ThisRank, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int DataCounter = 0;
        for (int k = 0; k < nz; k++) {
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < RecvYSlices(p); j++) {
                    IntVar_WholeDomain(k, i, j + RecvYOffset(p)) = RecvBufIntVar(DataCounter);
                    DataCounter++;
                }
            }
        }
    }
}

// On rank 0, collect data for one single float view
void CollectFloatField(ViewF3D_H FloatVar_WholeDomain, ViewF_H FloatVar, int nz, int nx, int MyYSlices, int np,
                       ViewI_H RecvYOffset, ViewI_H RecvYSlices, ViewI_H RBufSize) {

    // Set float variable to 0 for whole domain and place values for rank 0
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                FloatVar_WholeDomain(k, i, j) = FloatVar(k * nx * MyYSlices + i * MyYSlices + j);
            }
        }
    }

    // Recieve values from other ranks - message size different for different ranks
    for (int p = 1; p < np; p++) {
        int RecvBufSize_ThisRank = RBufSize(p);
        ViewF_H RecvBufFloatVar(Kokkos::ViewAllocateWithoutInitializing("RecvBufFloatVar"), RecvBufSize_ThisRank);
        MPI_Recv(RecvBufFloatVar.data(), RecvBufSize_ThisRank, MPI_FLOAT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int DataCounter = 0;
        for (int k = 0; k < nz; k++) {
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < RecvYSlices(p); j++) {
                    FloatVar_WholeDomain(k, i, j + RecvYOffset(p)) = RecvBufFloatVar(DataCounter);
                    DataCounter++;
                }
            }
        }
    }
}

// On rank 0, collect data for one single bool array (and convert it to integer 0s and 1s for MPI)
void CollectBoolField(ViewI3D_H IntVar_WholeDomain, bool *BoolVar, int nz, int nx, int MyYSlices, int np,
                      ViewI_H RecvYOffset, ViewI_H RecvYSlices, ViewI_H RBufSize) {

    // Resize bool variable for whole domain and place values for rank 0
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                if (BoolVar[k * nx * MyYSlices + i * MyYSlices + j])
                    IntVar_WholeDomain(k, i, j) = 1;
            }
        }
    }

    // Recieve values from other ranks - message size different for different ranks
    for (int p = 1; p < np; p++) {
        int RecvBufSize_ThisRank = RBufSize(p);
        ViewI_H RecvBufIntVar(Kokkos::ViewAllocateWithoutInitializing("RecvBufIntVar"), RecvBufSize_ThisRank);
        MPI_Recv(RecvBufIntVar.data(), RecvBufSize_ThisRank, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int DataCounter = 0;
        for (int k = 0; k < nz; k++) {
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < RecvYSlices[p]; j++) {
                    IntVar_WholeDomain(k, i, j + RecvYOffset(p)) = RecvBufIntVar(DataCounter);
                    DataCounter++;
                }
            }
        }
    }
}

//*****************************************************************************/
// On rank > 0, send data for an integer view to rank 0
void SendIntField(ViewI_H VarToSend, int nz, int nx, int MyYSlices, int SendBufSize, int SendBufStartY,
                  int SendBufEndY) {

    // Send non-ghost node data to rank 0
    int DataCounter = 0;
    ViewI_H SendBuf(Kokkos::ViewAllocateWithoutInitializing("SendBuf"), SendBufSize);
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = SendBufStartY; j < SendBufEndY; j++) {
                SendBuf(DataCounter) = VarToSend(k * nx * MyYSlices + i * MyYSlices + j);
                DataCounter++;
            }
        }
    }
    MPI_Send(SendBuf.data(), SendBufSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
}

// On rank > 0, send data for an float view to rank 0
void SendFloatField(ViewF_H VarToSend, int nz, int nx, int MyYSlices, int SendBufSize, int SendBufStartY,
                    int SendBufEndY) {

    // Send non-ghost node data to rank 0
    int DataCounter = 0;
    ViewF_H SendBuf(Kokkos::ViewAllocateWithoutInitializing("SendBuf"), SendBufSize);
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = SendBufStartY; j < SendBufEndY; j++) {
                SendBuf(DataCounter) = VarToSend(k * nx * MyYSlices + i * MyYSlices + j);
                DataCounter++;
            }
        }
    }
    MPI_Send(SendBuf.data(), SendBufSize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
}

// On rank > 0, send data for a bool array (converted into integers for MPI) to rank 0
void SendBoolField(bool *VarToSend, int nz, int nx, int MyYSlices, int SendBufSize, int SendBufStartY,
                   int SendBufEndY) {

    // Send non-ghost node data to rank 0
    int DataCounter = 0;
    ViewI_H SendBuf(Kokkos::ViewAllocateWithoutInitializing("SendBuf"), SendBufSize);
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            for (int j = SendBufStartY; j < SendBufEndY; j++) {
                if (VarToSend[k * nx * MyYSlices + i * MyYSlices + j])
                    SendBuf(DataCounter) = 1;
                else
                    SendBuf(DataCounter) = 0;
                DataCounter++;
            }
        }
    }
    MPI_Send(SendBuf.data(), SendBufSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
}

//*****************************************************************************/
// Prints values of selected data structures to Paraview files
void PrintExaCAData(int id, int layernumber, int np, int nx, int ny, int nz, int MyYSlices, int MyYOffset,
                    ViewI GrainID, ViewI CritTimeStep, ViewF GrainUnitVector, ViewI LayerID, ViewI CellType,
                    ViewF UndercoolingChange, ViewF UndercoolingCurrent, std::string BaseFileName,
                    int NGrainOrientations, std::string PathToOutput, int PrintDebug, bool PrintMisorientation,
                    bool PrintFinalUndercoolingVals, bool PrintFullOutput, bool PrintTimeSeries, bool PrintDefaultRVE,
                    int IntermediateFileCounter, int ZBound_Low, int nzActive, double deltax, double XMin, double YMin,
                    double ZMin, int NumberOfLayers, bool PrintBinary, int RVESize) {

    if (id == 0) {
        // Message sizes and data offsets for data recieved from other ranks- message size different for different ranks
        ViewI_H RecvYOffset(Kokkos::ViewAllocateWithoutInitializing("RecvYOffset"), np);
        ViewI_H RecvYSlices(Kokkos::ViewAllocateWithoutInitializing("RecvYSlices"), np);
        ViewI_H RBufSize(Kokkos::ViewAllocateWithoutInitializing("RBufSize"), np);

        for (int p = 1; p < np; p++) {
            RecvYOffset(p) = YOffsetCalc(p, ny, np);
            RecvYSlices(p) = YMPSlicesCalc(p, ny, np);
            RBufSize(p) = nx * RecvYSlices(p) * nz;
        }
        // Create variables for each possible data structure being collected on rank 0
        // If PrintDebug = 0, we aren't printing any debug files after initialization, but we are printing either the
        // misorientations (requiring GrainID) or we are printing the full end-of-run dataset (GrainID, LayerID)
        // If PrintDebug = 1, we are plotting CellType, LayerID, and CritTimeStep If PrintDebug = 2,
        // we are plotting all possible data structures
        // These views are allocated as size 0 x 0 x 0, and resized to nz by nx by ny, or size 0 by 0 by 0 if being
        // collected/printed on rank 0 (to avoid allocating memory unncessarily)
        // Utilize the fact that kokkos automatically initializes all values in these views to 0
        ViewI3D_H GrainID_WholeDomain(Kokkos::ViewAllocateWithoutInitializing("GrainID_WholeDomain"), 0, 0, 0);
        ViewI3D_H LayerID_WholeDomain(Kokkos::ViewAllocateWithoutInitializing("LayerID_WholeDomain"), 0, 0, 0);
        ViewI3D_H CritTimeStep_WholeDomain(
            Kokkos::ViewAllocateWithoutInitializing("CritTimeStep_WholeDomain_WholeDomain"), 0, 0, 0);
        ViewI3D_H CellType_WholeDomain(Kokkos::ViewAllocateWithoutInitializing("CellType_WholeDomain"), 0, 0, 0);
        ViewF3D_H UndercoolingChange_WholeDomain(
            Kokkos::ViewAllocateWithoutInitializing("UndercoolingChange_WholeDomain"), 0, 0, 0);
        ViewF3D_H UndercoolingCurrent_WholeDomain(
            Kokkos::ViewAllocateWithoutInitializing("UndercoolingCurrent_WholeDomain"), 0, 0, 0);

        // Collect all data on rank 0, for all data structures of interest
        ViewI_H GrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);
        Kokkos::resize(GrainID_WholeDomain, nz, nx, ny);
        CollectIntField(GrainID_WholeDomain, GrainID_Host, nz, nx, MyYSlices, np, RecvYOffset, RecvYSlices, RBufSize);
        ViewI_H LayerID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), LayerID);
        Kokkos::resize(LayerID_WholeDomain, nz, nx, ny);
        CollectIntField(LayerID_WholeDomain, LayerID_Host, nz, nx, MyYSlices, np, RecvYOffset, RecvYSlices, RBufSize);
        if ((PrintDebug > 0) || (PrintTimeSeries)) {
            ViewI_H CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
            Kokkos::resize(CellType_WholeDomain, nz, nx, ny);
            CollectIntField(CellType_WholeDomain, CellType_Host, nz, nx, MyYSlices, np, RecvYOffset, RecvYSlices,
                            RBufSize);
        }
        if (PrintDebug > 0) {
            ViewI_H CritTimeStep_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CritTimeStep);
            Kokkos::resize(CritTimeStep_WholeDomain, nz, nx, ny);
            CollectIntField(CritTimeStep_WholeDomain, CritTimeStep_Host, nz, nx, MyYSlices, np, RecvYOffset,
                            RecvYSlices, RBufSize);
        }
        if (PrintDebug == 2) {
            ViewF_H UndercoolingChange_Host =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), UndercoolingChange);
            Kokkos::resize(UndercoolingChange_WholeDomain, nz, nx, ny);
            CollectFloatField(UndercoolingChange_WholeDomain, UndercoolingChange_Host, nz, nx, MyYSlices, np,
                              RecvYOffset, RecvYSlices, RBufSize);
        }
        if ((PrintDebug == 2) || (PrintFinalUndercoolingVals)) {
            ViewF_H UndercoolingCurrent_Host =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), UndercoolingCurrent);
            Kokkos::resize(UndercoolingCurrent_WholeDomain, nz, nx, ny);
            CollectFloatField(UndercoolingCurrent_WholeDomain, UndercoolingCurrent_Host, nz, nx, MyYSlices, np,
                              RecvYOffset, RecvYSlices, RBufSize);
        }
        if ((PrintMisorientation) || (PrintTimeSeries)) {
            ViewF_H GrainUnitVector_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainUnitVector);
            if (PrintMisorientation)
                PrintGrainMisorientations(BaseFileName, PathToOutput, nx, ny, nz, LayerID_WholeDomain,
                                          GrainID_WholeDomain, GrainUnitVector_Host, NGrainOrientations, deltax, XMin,
                                          YMin, ZMin, PrintBinary);
            if (PrintTimeSeries)
                PrintIntermediateExaCAState(IntermediateFileCounter, layernumber, BaseFileName, PathToOutput,
                                            ZBound_Low, nzActive, nx, ny, GrainID_WholeDomain, CellType_WholeDomain,
                                            GrainUnitVector_Host, NGrainOrientations, deltax, XMin, YMin, ZMin,
                                            PrintBinary);
        }
        if (PrintFinalUndercoolingVals)
            PrintFinalUndercooling(BaseFileName, PathToOutput, nx, ny, nz, UndercoolingCurrent_WholeDomain,
                                   LayerID_WholeDomain, deltax, XMin, YMin, ZMin, PrintBinary);
        if (PrintDefaultRVE)
            PrintExaConstitDefaultRVE(BaseFileName, PathToOutput, nx, ny, nz, LayerID_WholeDomain, GrainID_WholeDomain,
                                      deltax, NumberOfLayers, RVESize);
        if ((PrintFullOutput) || (PrintDebug > 0))
            PrintCAFields(nx, ny, nz, GrainID_WholeDomain, LayerID_WholeDomain, CritTimeStep_WholeDomain,
                          CellType_WholeDomain, UndercoolingChange_WholeDomain, UndercoolingCurrent_WholeDomain,
                          PathToOutput, BaseFileName, PrintDebug, PrintFullOutput, deltax, XMin, YMin, ZMin,
                          PrintBinary);
    }
    else {

        // No ghost nodes in sent data:
        int SendBufStartY, SendBufEndY;
        if (MyYOffset == 0)
            SendBufStartY = 0;
        else
            SendBufStartY = 1;
        if (MyYSlices + MyYOffset == ny)
            SendBufEndY = MyYSlices;
        else
            SendBufEndY = MyYSlices - 1;

        int SendBufSize = nx * (SendBufEndY - SendBufStartY) * nz;

        // Send Grain ID and Layer ID data to rank 0 for all print options
        ViewI_H GrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);
        SendIntField(GrainID_Host, nz, nx, MyYSlices, SendBufSize, SendBufStartY, SendBufEndY);
        ViewI_H LayerID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), LayerID);
        SendIntField(LayerID_Host, nz, nx, MyYSlices, SendBufSize, SendBufStartY, SendBufEndY);
        if ((PrintDebug > 0) || (PrintTimeSeries)) {
            ViewI_H CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
            SendIntField(CellType_Host, nz, nx, MyYSlices, SendBufSize, SendBufStartY, SendBufEndY);
        }
        if (PrintDebug > 0) {
            ViewI_H CritTimeStep_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CritTimeStep);
            SendIntField(CritTimeStep_Host, nz, nx, MyYSlices, SendBufSize, SendBufStartY, SendBufEndY);
        }
        if (PrintDebug == 2) {
            ViewF_H UndercoolingChange_Host =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), UndercoolingChange);
            SendFloatField(UndercoolingChange_Host, nz, nx, MyYSlices, SendBufSize, SendBufStartY, SendBufEndY);
        }
        if ((PrintDebug == 2) || (PrintFinalUndercoolingVals)) {
            ViewF_H UndercoolingCurrent_Host =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), UndercoolingCurrent);
            SendFloatField(UndercoolingCurrent_Host, nz, nx, MyYSlices, SendBufSize, SendBufStartY, SendBufEndY);
        }
    }
}

//*****************************************************************************/
// Print specified fields to a paraview file
void PrintCAFields(int nx, int ny, int nz, ViewI3D_H GrainID_WholeDomain, ViewI3D_H LayerID_WholeDomain,
                   ViewI3D_H CritTimeStep_WholeDomain, ViewI3D_H CellType_WholeDomain,
                   ViewF3D_H UndercoolingChange_WholeDomain, ViewF3D_H UndercoolingCurrent_WholeDomain,
                   std::string PathToOutput, std::string BaseFileName, int PrintDebug, bool PrintFullOutput,
                   double deltax, double XMin, double YMin, double ZMin, bool PrintBinary) {

    std::string FName;
    if (PrintFullOutput) {
        FName = PathToOutput + BaseFileName + ".vtk";
        std::cout << "Printing layer ID, grain ID, and melted 0/1 data to a vtk file " << FName
                  << " for post-processing" << std::endl;
    }
    else {
        FName = PathToOutput + BaseFileName + "_debug.vtk";
        std::cout << "Printing initial data structures to a vtk file " << FName << " for debugging" << std::endl;
    }
    std::ofstream Grainplot;
    WriteHeader(Grainplot, FName, PrintBinary, nx, ny, nz, deltax, XMin, YMin, ZMin);
    // Print Layer ID data
    Grainplot << "SCALARS LayerID short 1" << std::endl;
    Grainplot << "LOOKUP_TABLE default" << std::endl;
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                short ShortLayerID = static_cast<short>(LayerID_WholeDomain(k, i, j));
                WriteData(Grainplot, ShortLayerID, PrintBinary, true);
            }
        }
        // Do not insert newline character if using binary writing, as this will break the binary data read by adding a
        // blank line
        if (!(PrintBinary))
            Grainplot << std::endl;
    }
    if ((PrintFullOutput) || (PrintDebug == 2)) {
        // Print Grain ID data
        Grainplot << "SCALARS GrainID int 1" << std::endl;
        Grainplot << "LOOKUP_TABLE default" << std::endl;
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    WriteData(Grainplot, GrainID_WholeDomain(k, i, j), PrintBinary, true);
                }
            }
            if (!(PrintBinary))
                Grainplot << std::endl;
        }
    }
    if (PrintDebug > 0) {
        // Print CritTimeStep data
        Grainplot << "SCALARS CritTimeStep int 1" << std::endl;
        Grainplot << "LOOKUP_TABLE default" << std::endl;
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    WriteData(Grainplot, CritTimeStep_WholeDomain(k, i, j), PrintBinary, true);
                }
            }
            if (!(PrintBinary))
                Grainplot << std::endl;
        }
        // Print CellType data
        Grainplot << "SCALARS CellType int 1" << std::endl;
        Grainplot << "LOOKUP_TABLE default" << std::endl;
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    WriteData(Grainplot, CellType_WholeDomain(k, i, j), PrintBinary, true);
                }
            }
            if (!(PrintBinary))
                Grainplot << std::endl;
        }
    }
    if (PrintDebug == 2) {
        // Print Undercooling change data
        Grainplot << "SCALARS UndercoolingChange float 1" << std::endl;
        Grainplot << "LOOKUP_TABLE default" << std::endl;
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    WriteData(Grainplot, UndercoolingChange_WholeDomain(k, i, j), PrintBinary, true);
                }
            }
            if (!(PrintBinary))
                Grainplot << std::endl;
        }
        // Print Undercooling current data
        Grainplot << "SCALARS UndercoolingCurrent float 1" << std::endl;
        Grainplot << "LOOKUP_TABLE default" << std::endl;
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    WriteData(Grainplot, UndercoolingCurrent_WholeDomain(k, i, j), PrintBinary, true);
                }
            }
            if (!(PrintBinary))
                Grainplot << std::endl;
        }
    }
    Grainplot.close();
}
//*****************************************************************************/
// Print grain misorientation, 0-62 for epitaxial grains and 100-162 for nucleated grains, to a paraview file
void PrintGrainMisorientations(std::string BaseFileName, std::string PathToOutput, int nx, int ny, int nz,
                               ViewI3D_H LayerID_WholeDomain, ViewI3D_H GrainID_WholeDomain, ViewF_H GrainUnitVector,
                               int NGrainOrientations, double deltax, double XMin, double YMin, double ZMin,
                               bool PrintBinary) {

    std::string FName = PathToOutput + BaseFileName + "_Misorientations.vtk";
    std::cout << "Printing Paraview file of grain misorientations" << std::endl;
    // Print grain orientations to file
    std::ofstream GrainplotM;
    WriteHeader(GrainplotM, FName, PrintBinary, nx, ny, nz, deltax, XMin, YMin, ZMin);
    GrainplotM << "SCALARS Angle_z unsigned_short 1" << std::endl;
    GrainplotM << "LOOKUP_TABLE default" << std::endl;

    int NucleatedGrainCells = 0;
    int MeltedCells = 0;
    // Get grain misorientations relative to the Z direction for each orientation
    ViewF_H GrainMisorientation = MisorientationCalc(NGrainOrientations, GrainUnitVector, 2);
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                unsigned short IntPrintVal;
                if (LayerID_WholeDomain(k, i, j) == -1)
                    IntPrintVal = 200;
                else {
                    MeltedCells++;
                    int MyOrientation = getGrainOrientation(GrainID_WholeDomain(k, i, j), NGrainOrientations);
                    if (GrainID_WholeDomain(k, i, j) < 0) {
                        IntPrintVal = static_cast<unsigned short>(std::round(GrainMisorientation(MyOrientation)) + 100);
                        NucleatedGrainCells++;
                    }
                    else
                        IntPrintVal = static_cast<unsigned short>(std::round(GrainMisorientation(MyOrientation)));
                }
                WriteData(GrainplotM, IntPrintVal, PrintBinary, true);
            }
        }
        if (!(PrintBinary))
            GrainplotM << std::endl;
    }
    GrainplotM.close();
    std::cout << "Volume fraction of solidified portion of domain claimed by nucleated grains: "
              << (float)(NucleatedGrainCells) / (float)(MeltedCells) << std::endl;
}

//*****************************************************************************/
// Print the undercooling at which each cell solidified (went from active to solid type), for all cells that underwent
// solidification, to a paraview file
void PrintFinalUndercooling(std::string BaseFileName, std::string PathToOutput, int nx, int ny, int nz,
                            ViewF3D_H UndercoolingCurrent_WholeDomain, ViewI3D_H LayerID_WholeDomain, double deltax,
                            double XMin, double YMin, double ZMin, bool PrintBinary) {

    std::string FName = PathToOutput + BaseFileName + "_FinalUndercooling.vtk";
    std::cout << "Printing Paraview file of final undercooling values" << std::endl;
    // Print undercooling to file
    std::ofstream UndercoolingPlot;
    WriteHeader(UndercoolingPlot, FName, PrintBinary, nx, ny, nz, deltax, XMin, YMin, ZMin);
    UndercoolingPlot << "SCALARS UndercoolingFinal float 1" << std::endl;
    UndercoolingPlot << "LOOKUP_TABLE default" << std::endl;
    // Print 0s for the undercooling for cells that did not undergo solidification
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                float FloatPrintVal;
                if (LayerID_WholeDomain(k, i, j) == -1)
                    FloatPrintVal = 0.0;
                else
                    FloatPrintVal = UndercoolingCurrent_WholeDomain(k, i, j);
                WriteData(UndercoolingPlot, FloatPrintVal, PrintBinary, true);
            }
        }
        if (!(PrintBinary))
            UndercoolingPlot << std::endl;
    }
    UndercoolingPlot.close();
}

//*****************************************************************************/
// Print a representative volume element (RVE) from this multilayer simulation from the "default" location in the
// domain. The default location is as close to the center of the domain in X and Y as possible, and as close to the top
// of the domain while not including the final layer's microstructure. If an RVE size was not specified in the input
// file, the default size is 0.5 by 0.5 by 0.5 mm
void PrintExaConstitDefaultRVE(std::string BaseFileName, std::string PathToOutput, int nx, int ny, int nz,
                               ViewI3D_H LayerID_WholeDomain, ViewI3D_H GrainID_WholeDomain, double deltax,
                               int NumberOfLayers, int RVESize) {

    // Determine the lower and upper Y bounds of the RVE
    int RVE_XLow = std::floor(nx / 2) - std::floor(RVESize / 2);
    int RVE_XHigh = RVE_XLow + RVESize - 1;
    int RVE_YLow = std::floor(ny / 2) - std::floor(RVESize / 2);
    int RVE_YHigh = RVE_YLow + RVESize - 1;

    // Make sure the RVE fits in the simulation domain in X and Y
    if ((RVE_XLow < 0) || (RVE_XHigh > nx - 1) || (RVE_YLow < 0) || (RVE_YHigh > ny - 1)) {
        std::cout << "WARNING: Simulation domain is too small to obtain default RVE data (should be at least "
                  << RVESize << " cells in the X and Y directions" << std::endl;
        if (RVE_XLow < 0)
            RVE_XLow = 0;
        if (RVE_XHigh > nx - 1)
            RVE_XHigh = nx - 1;
        if (RVE_YLow < 0)
            RVE_YLow = 0;
        if (RVE_YHigh > ny - 1)
            RVE_YHigh = ny - 1;
    }

    // Determine the upper Z bound of the RVE - largest Z for which all cells in the RVE do not contain LayerID values
    // of the last layer
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
    int RVE_ZLow = RVE_ZHigh - RVESize + 1;
    if (RVE_ZLow < 0) {
        std::cout << "WARNING: Simulation domain is too small to obtain default RVE data (should be at least "
                  << RVESize << " cells in the Z direction, more layers are required" << std::endl;
        RVE_ZLow = 0;
        RVE_ZHigh = nz - 1;
    }

    // Print RVE data to file
    std::string FName = PathToOutput + BaseFileName + "_ExaConstit.csv";
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

//*****************************************************************************/
// Print a log file for this ExaCA run in json file format, containing information about the run parameters used
// from the input file as well as the decomposition scheme
void PrintExaCALog(int id, int np, std::string InputFile, std::string SimulationType, int MyYSlices, int MyYOffset,
                   InterfacialResponseFunction irf, double deltax, double NMax, double dTN, double dTsigma,
                   std::vector<std::string> temp_paths, int TempFilesInSeries, double HT_deltax, bool RemeltingYN,
                   double deltat, int NumberOfLayers, int LayerHeight, std::string SubstrateFileName,
                   double SubstrateGrainSpacing, bool SubstrateFile, double G, double R, int nx, int ny, int nz,
                   double FractSurfaceSitesActive, std::string PathToOutput, int NSpotsX, int NSpotsY, int SpotOffset,
                   int SpotRadius, std::string BaseFileName, double InitTime, double RunTime, double OutTime, int cycle,
                   double InitMaxTime, double InitMinTime, double NuclMaxTime, double NuclMinTime,
                   double CreateSVMinTime, double CreateSVMaxTime, double CaptureMaxTime, double CaptureMinTime,
                   double GhostMaxTime, double GhostMinTime, double OutMaxTime, double OutMinTime, double XMin,
                   double XMax, double YMin, double YMax, double ZMin, double ZMax, std::string GrainOrientationFile) {

    int *YSlices = new int[np];
    int *YOffset = new int[np];
    MPI_Gather(&MyYSlices, 1, MPI_INT, YSlices, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&MyYOffset, 1, MPI_INT, YOffset, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (id == 0) {
        std::string FName = PathToOutput + BaseFileName + ".json";
        std::cout << "Printing ExaCA log file" << std::endl;
        std::ofstream ExaCALog;
        ExaCALog.open(FName);
        ExaCALog << "{" << std::endl;
        ExaCALog << "   \"ExaCAVersion\": \"" << version() << "\", " << std::endl;
        ExaCALog << "   \"ExaCACommitHash\": \"" << gitCommitHash() << "\", " << std::endl;
        ExaCALog << "   \"KokkosVersion\": \"" << kokkosVersion() << "\", " << std::endl;
        ExaCALog << "   \"InputFile\": \"" << InputFile << "\", " << std::endl;
        ExaCALog << "   \"TimeStepOfOutput\": " << cycle << "," << std::endl;
        ExaCALog << "   \"SimulationType\": \"" << SimulationType << "\"," << std::endl;
        ExaCALog << "   \"GrainOrientationFile\": \"" << GrainOrientationFile << "\"," << std::endl;
        ExaCALog << "   \"Domain\": {" << std::endl;
        ExaCALog << "      \"Nx\": " << nx << "," << std::endl;
        ExaCALog << "      \"Ny\": " << ny << "," << std::endl;
        ExaCALog << "      \"Nz\": " << nz << "," << std::endl;
        ExaCALog << "      \"CellSize\": " << deltax << "," << std::endl;
        ExaCALog << "      \"TimeStep\": " << deltat << "," << std::endl;
        ExaCALog << "      \"XBounds\": [" << XMin << "," << XMax << "]," << std::endl;
        ExaCALog << "      \"YBounds\": [" << YMin << "," << YMax << "]," << std::endl;
        ExaCALog << "      \"ZBounds\": [" << ZMin << "," << ZMax << "]";
        if (SimulationType != "C") {
            ExaCALog << "," << std::endl;
            ExaCALog << "      \"NumberOfLayers\": " << NumberOfLayers << "," << std::endl;
            ExaCALog << "      \"LayerOffset\": " << LayerHeight << "," << std::endl;
            ExaCALog << "      \"Remelting\": " << RemeltingYN;
            if (SimulationType == "S") {
                ExaCALog << "," << std::endl;
                ExaCALog << "      \"NSpotsX\": " << NSpotsX << "," << std::endl;
                ExaCALog << "      \"NSpotsY\": " << NSpotsY << "," << std::endl;
                ExaCALog << "      \"RSpots\": " << SpotRadius << "," << std::endl;
                ExaCALog << "      \"SpotOffset\": " << SpotOffset << std::endl;
            }
            else
                ExaCALog << std::endl;
        }
        else
            ExaCALog << std::endl;
        ExaCALog << "   }," << std::endl;
        ExaCALog << "   \"Nucleation\": {" << std::endl;
        ExaCALog << "      \"Density\": " << NMax << "," << std::endl;
        ExaCALog << "      \"MeanUndercooling\": " << dTN << "," << std::endl;
        ExaCALog << "      \"StDevUndercooling\": " << dTsigma << std::endl;
        ExaCALog << "   }," << std::endl;
        ExaCALog << "   \"TemperatureData\": {" << std::endl;
        if (SimulationType == "R") {
            ExaCALog << "       \"TemperatureFiles\": [";
            for (int i = 0; i < TempFilesInSeries - 1; i++) {
                ExaCALog << "\"" << temp_paths[i] << "\", ";
            }
            ExaCALog << "\"" << temp_paths[TempFilesInSeries - 1] << "\"]," << std::endl;
            ExaCALog << "       \"HeatTransferCellSize\": " << HT_deltax << std::endl;
        }
        else {
            ExaCALog << "      \"G\": " << G << "," << std::endl;
            ExaCALog << "      \"R\": " << R << std::endl;
        }
        ExaCALog << "   }," << std::endl;
        ExaCALog << "   \"Substrate\": {" << std::endl;
        if (SimulationType == "C")
            ExaCALog << "       \"FractionSurfaceSitesActive\": " << FractSurfaceSitesActive << std::endl;
        else {
            if (SubstrateFile)
                ExaCALog << "       \"SubstrateFilename\": " << SubstrateFileName << std::endl;
            else
                ExaCALog << "       \"MeanSize\": " << SubstrateGrainSpacing << std::endl;
        }
        ExaCALog << "   }," << std::endl;
        ExaCALog << irf.print() << std::endl;
        ExaCALog << "   \"NumberMPIRanks\": " << np << "," << std::endl;
        ExaCALog << "   \"Decomposition\": {" << std::endl;
        ExaCALog << "       \"SubdomainYSize\": [";
        for (int i = 0; i < np - 1; i++)
            ExaCALog << YSlices[i] << ",";
        ExaCALog << YSlices[np - 1] << "]," << std::endl;
        ExaCALog << "       \"SubdomainYOffset\": [";
        for (int i = 0; i < np - 1; i++)
            ExaCALog << YOffset[i] << ",";
        ExaCALog << YOffset[np - 1] << "]" << std::endl;
        ExaCALog << "   }," << std::endl;
        ExaCALog << "   \"Timing\": {" << std::endl;
        ExaCALog << "       \"Runtime\": " << InitTime + RunTime + OutTime << "," << std::endl;
        ExaCALog << "       \"InitRunOutputBreakdown\": [" << InitTime << "," << RunTime << "," << OutTime << "],"
                 << std::endl;
        ExaCALog << "       \"MaxMinInitTime\": [" << InitMaxTime << "," << InitMinTime << "]," << std::endl;
        ExaCALog << "       \"MaxMinNucleationTime\": [" << NuclMaxTime << "," << NuclMinTime << "]," << std::endl;
        ExaCALog << "       \"MaxMinSteeringVectorCreationTime\": [" << CreateSVMaxTime << "," << CreateSVMinTime
                 << "]," << std::endl;
        ExaCALog << "       \"MaxMinCellCaptureTime\": [" << CaptureMaxTime << "," << CaptureMinTime << "],"
                 << std::endl;
        ExaCALog << "       \"MaxMinGhostExchangeTime\": [" << GhostMaxTime << "," << GhostMinTime << "]," << std::endl;
        ExaCALog << "       \"MaxMinOutputTime\": [" << OutMaxTime << "," << OutMinTime << "]" << std::endl;
        ExaCALog << "   }" << std::endl;
        ExaCALog << "}" << std::endl;
        ExaCALog.close();
    }
}

// Print timing info to console
void PrintExaCATiming(int np, double InitTime, double RunTime, double OutTime, int cycle, double InitMaxTime,
                      double InitMinTime, double NuclMaxTime, double NuclMinTime, double CreateSVMinTime,
                      double CreateSVMaxTime, double CaptureMaxTime, double CaptureMinTime, double GhostMaxTime,
                      double GhostMinTime, double OutMaxTime, double OutMinTime) {

    std::cout << "===================================================================================" << std::endl;
    std::cout << "Having run with = " << np << " processors" << std::endl;
    std::cout << "Output written at cycle = " << cycle << std::endl;
    std::cout << "Total time = " << InitTime + RunTime + OutTime << std::endl;
    std::cout << "Time spent initializing data = " << InitTime << " s" << std::endl;
    std::cout << "Time spent performing CA calculations = " << RunTime << " s" << std::endl;
    std::cout << "Time spent collecting and printing output data = " << OutTime << " s\n" << std::endl;

    std::cout << "Max/min rank time initializing data  = " << InitMaxTime << " / " << InitMinTime << " s" << std::endl;
    std::cout << "Max/min rank time in CA nucleation   = " << NuclMaxTime << " / " << NuclMinTime << " s" << std::endl;
    std::cout << "Max/min rank time in CA steering vector creation = " << CreateSVMaxTime << " / " << CreateSVMinTime
              << " s" << std::endl;
    std::cout << "Max/min rank time in CA cell capture = " << CaptureMaxTime << " / " << CaptureMinTime << " s"
              << std::endl;
    std::cout << "Max/min rank time in CA ghosting     = " << GhostMaxTime << " / " << GhostMinTime << " s"
              << std::endl;
    std::cout << "Max/min rank time exporting data     = " << OutMaxTime << " / " << OutMinTime << " s\n" << std::endl;

    std::cout << "===================================================================================" << std::endl;
}

// Prints grain misorientation and cell type at some intermediate time during the simulation
// Files are named with the ending "_[FrameNumber].vtk", allowing Paraview to read these in as a series to make videos
// Cells that are liquid are given the value -1
// Epitaxial grains are colored 0-62 (by misorientation with +Z direction)
// Nucleated grains are colored 100-162 (by misorientation with +Z direction plus 100)
void PrintIntermediateExaCAState(int IntermediateFileCounter, int layernumber, std::string BaseFileName,
                                 std::string PathToOutput, int ZBound_Low, int nzActive, int nx, int ny,
                                 ViewI3D_H GrainID_WholeDomain, ViewI3D_H CellType_WholeDomain, ViewF_H GrainUnitVector,
                                 int NGrainOrientations, double deltax, double XMin, double YMin, double ZMin,
                                 bool PrintBinary) {

    std::string FName = PathToOutput + BaseFileName + "_layer" + std::to_string(layernumber) + "_" +
                        std::to_string(IntermediateFileCounter) + ".vtk";

    // Size of output in Z will depend on the current layer bounds - always start from Z = 0, but only collect data up
    // to Z = ZBound_Low + nzActive
    int ZPrintSize = ZBound_Low + nzActive;
    std::cout << "Printing file " << FName << " : Z coordinates of 0 through " << ZPrintSize - 1 << std::endl;
    ViewF_H GrainMisorientation(Kokkos::ViewAllocateWithoutInitializing("GrainMisorientation"), NGrainOrientations);
    for (int n = 0; n < NGrainOrientations; n++) {
        // Find the smallest possible misorientation between the domain +Z direction, and this grain orientations' 6
        // possible 001 directions (where 62.7 degrees is the largest possible misorientation between two 001 directions
        // for a cubic crystal system)
        float AngleZmin = 62.7;
        for (int ll = 0; ll < 3; ll++) {
            float AngleZ = std::abs((180 / M_PI) * std::acos(GrainUnitVector(9 * n + 3 * ll + 2)));
            if (AngleZ < AngleZmin) {
                AngleZmin = AngleZ;
            }
        }
        GrainMisorientation(n) = AngleZmin;
    }
    std::ofstream GrainplotM;
    WriteHeader(GrainplotM, FName, PrintBinary, nx, ny, ZPrintSize, deltax, XMin, YMin, ZMin);
    GrainplotM << "SCALARS Angle_z float 1" << std::endl;
    GrainplotM << "LOOKUP_TABLE default" << std::endl;
    for (int k = 0; k < ZPrintSize; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                float FloatPrintVal;
                if (CellType_WholeDomain(k, i, j) != Liquid) {
                    int MyOrientation = getGrainOrientation(GrainID_WholeDomain(k, i, j), NGrainOrientations);
                    if (GrainID_WholeDomain(k, i, j) < 0)
                        FloatPrintVal = GrainMisorientation(MyOrientation) + 100.0;
                    else
                        FloatPrintVal = GrainMisorientation(MyOrientation);
                }
                else
                    FloatPrintVal = -1.0;
                WriteData(GrainplotM, FloatPrintVal, PrintBinary, true);
            }
        }
        if (!(PrintBinary))
            GrainplotM << std::endl;
    }
    GrainplotM.close();
}

std::string version() { return ExaCA_VERSION; }

std::string gitCommitHash() { return ExaCA_GIT_COMMIT_HASH; }

std::string kokkosVersion() { return ExaCA_Kokkos_VERSION_STRING; }
