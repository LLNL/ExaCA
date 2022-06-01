// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAprint.hpp"

#include "CAfunctions.hpp"

#include "mpi.h"

#include <cmath>
#include <fstream>

//*****************************************************************************/
//*****************************************************************************/
// On rank 0, collect data for one single int view
void CollectIntField(ViewI3D_H IntVar_WholeDomain, ViewI_H IntVar, int nz, int MyXSlices, int MyYSlices, int np,
                     ViewI_H RecvXOffset, ViewI_H RecvYOffset, ViewI_H RecvXSlices, ViewI_H RecvYSlices,
                     ViewI_H RBufSize) {

    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < MyXSlices; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                IntVar_WholeDomain(k, i, j) = IntVar(k * MyXSlices * MyYSlices + i * MyYSlices + j);
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
            for (int i = 0; i < RecvXSlices(p); i++) {
                for (int j = 0; j < RecvYSlices(p); j++) {
                    IntVar_WholeDomain(k, i + RecvXOffset(p), j + RecvYOffset(p)) = RecvBufIntVar(DataCounter);
                    DataCounter++;
                }
            }
        }
    }
}

// On rank 0, collect data for one single float view
void CollectFloatField(ViewF3D_H FloatVar_WholeDomain, ViewF_H FloatVar, int nz, int MyXSlices, int MyYSlices, int np,
                       ViewI_H RecvXOffset, ViewI_H RecvYOffset, ViewI_H RecvXSlices, ViewI_H RecvYSlices,
                       ViewI_H RBufSize) {

    // Set float variable to 0 for whole domain and place values for rank 0
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < MyXSlices; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                FloatVar_WholeDomain(k, i, j) = FloatVar(k * MyXSlices * MyYSlices + i * MyYSlices + j);
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
            for (int i = 0; i < RecvXSlices(p); i++) {
                for (int j = 0; j < RecvYSlices(p); j++) {
                    FloatVar_WholeDomain(k, i + RecvXOffset(p), j + RecvYOffset(p)) = RecvBufFloatVar(DataCounter);
                    DataCounter++;
                }
            }
        }
    }
}

// On rank 0, collect data for one single bool array (and convert it to integer 0s and 1s for MPI)
void CollectBoolField(ViewI3D_H IntVar_WholeDomain, bool *BoolVar, int nz, int MyXSlices, int MyYSlices, int np,
                      ViewI_H RecvXOffset, ViewI_H RecvYOffset, ViewI_H RecvXSlices, ViewI_H RecvYSlices,
                      ViewI_H RBufSize) {

    // Resize bool variable for whole domain and place values for rank 0
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < MyXSlices; i++) {
            for (int j = 0; j < MyYSlices; j++) {
                if (BoolVar[k * MyXSlices * MyYSlices + i * MyYSlices + j])
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
            for (int i = 0; i < RecvXSlices[p]; i++) {
                for (int j = 0; j < RecvYSlices[p]; j++) {
                    IntVar_WholeDomain(k, i + RecvXOffset(p), j + RecvYOffset(p)) = RecvBufIntVar(DataCounter);
                    DataCounter++;
                }
            }
        }
    }
}

//*****************************************************************************/
// On rank > 0, send data for an integer view to rank 0
void SendIntField(ViewI_H VarToSend, int nz, int MyXSlices, int MyYSlices, int SendBufSize, int SendBufStartX,
                  int SendBufEndX, int SendBufStartY, int SendBufEndY) {

    // Send non-ghost node data to rank 0
    int DataCounter = 0;
    ViewI_H SendBuf(Kokkos::ViewAllocateWithoutInitializing("SendBuf"), SendBufSize);
    for (int k = 0; k < nz; k++) {
        for (int i = SendBufStartX; i < SendBufEndX; i++) {
            for (int j = SendBufStartY; j < SendBufEndY; j++) {
                SendBuf(DataCounter) = VarToSend(k * MyXSlices * MyYSlices + i * MyYSlices + j);
                DataCounter++;
            }
        }
    }
    MPI_Send(SendBuf.data(), SendBufSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
}

// On rank > 0, send data for an float view to rank 0
void SendFloatField(ViewF_H VarToSend, int nz, int MyXSlices, int MyYSlices, int SendBufSize, int SendBufStartX,
                    int SendBufEndX, int SendBufStartY, int SendBufEndY) {

    // Send non-ghost node data to rank 0
    int DataCounter = 0;
    ViewF_H SendBuf(Kokkos::ViewAllocateWithoutInitializing("SendBuf"), SendBufSize);
    for (int k = 0; k < nz; k++) {
        for (int i = SendBufStartX; i < SendBufEndX; i++) {
            for (int j = SendBufStartY; j < SendBufEndY; j++) {
                SendBuf(DataCounter) = VarToSend(k * MyXSlices * MyYSlices + i * MyYSlices + j);
                DataCounter++;
            }
        }
    }
    MPI_Send(SendBuf.data(), SendBufSize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
}

// On rank > 0, send data for a bool array (converted into integers for MPI) to rank 0
void SendBoolField(bool *VarToSend, int nz, int MyXSlices, int MyYSlices, int SendBufSize, int SendBufStartX,
                   int SendBufEndX, int SendBufStartY, int SendBufEndY) {

    // Send non-ghost node data to rank 0
    int DataCounter = 0;
    ViewI_H SendBuf(Kokkos::ViewAllocateWithoutInitializing("SendBuf"), SendBufSize);
    for (int k = 0; k < nz; k++) {
        for (int i = SendBufStartX; i < SendBufEndX; i++) {
            for (int j = SendBufStartY; j < SendBufEndY; j++) {
                if (VarToSend[k * MyXSlices * MyYSlices + i * MyYSlices + j])
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
void PrintExaCAData(int id, int layernumber, int np, int nx, int ny, int nz, int MyXSlices, int MyYSlices,
                    int MyXOffset, int MyYOffset, int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI GrainID,
                    ViewI CritTimeStep, ViewF GrainUnitVector, ViewI LayerID, ViewI CellType, ViewF UndercoolingChange,
                    ViewF UndercoolingCurrent, std::string BaseFileName, int DecompositionStrategy,
                    int NGrainOrientations, bool *Melted, std::string PathToOutput, int PrintDebug,
                    bool PrintMisorientation, bool PrintFinalUndercoolingVals, bool PrintFullOutput,
                    bool PrintTimeSeries, bool PrintDefaultRVE, int IntermediateFileCounter, int ZBound_Low,
                    int nzActive, double deltax, float XMin, float YMin, float ZMin, int NumberOfLayers) {

    if (id == 0) {
        // Message sizes and data offsets for data recieved from other ranks- message size different for different ranks
        ViewI_H RecvXOffset(Kokkos::ViewAllocateWithoutInitializing("RecvXOffset"), np);
        ViewI_H RecvYOffset(Kokkos::ViewAllocateWithoutInitializing("RecvYOffset"), np);
        ViewI_H RecvXSlices(Kokkos::ViewAllocateWithoutInitializing("RecvXSlices"), np);
        ViewI_H RecvYSlices(Kokkos::ViewAllocateWithoutInitializing("RecvYSlices"), np);
        ViewI_H RBufSize(Kokkos::ViewAllocateWithoutInitializing("RBufSize"), np);

        for (int p = 1; p < np; p++) {
            RecvXOffset(p) = XOffsetCalc(p, nx, ProcessorsInXDirection, ProcessorsInYDirection, DecompositionStrategy);
            RecvXSlices(p) =
                XMPSlicesCalc(p, nx, ProcessorsInXDirection, ProcessorsInYDirection, DecompositionStrategy);
            RecvYOffset(p) = YOffsetCalc(p, ny, ProcessorsInYDirection, np, DecompositionStrategy);
            RecvYSlices(p) = YMPSlicesCalc(p, ny, ProcessorsInYDirection, np, DecompositionStrategy);

            RBufSize(p) = RecvXSlices(p) * RecvYSlices(p) * nz;
        }
        // Create variables for each possible data structure being collected on rank 0
        // If PrintDebug = 0, we aren't printing any debug files after initialization, but we are printing either the
        // misorientations (requiring Melted and GrainID) or we are printing the full end-of-run dataset (Melted,
        // GrainID, LayerID) If PrintDebug = 1, we are plotting CellType, LayerID, and CritTimeStep If PrintDebug = 2,
        // we are plotting all possible data structures
        // These views are allocated as size 0 x 0 x 0, and resized to nz by nx by ny, or size 0 by 0 by 0 if being
        // collected/printed on rank 0 (to avoid allocating memory unncessarily)
        // Utilize the fact that kokkos automatically initializes all values in these views to 0
        ViewI3D_H GrainID_WholeDomain(Kokkos::ViewAllocateWithoutInitializing("GrainID_WholeDomain"), 0, 0, 0);
        ViewI3D_H LayerID_WholeDomain(Kokkos::ViewAllocateWithoutInitializing("LayerID_WholeDomain"), 0, 0, 0);
        ViewI3D_H Melted_WholeDomain(Kokkos::ViewAllocateWithoutInitializing("Melted_WholeDomain"), 0, 0, 0);
        ViewI3D_H CritTimeStep_WholeDomain(
            Kokkos::ViewAllocateWithoutInitializing("CritTimeStep_WholeDomain_WholeDomain"), 0, 0, 0);
        ViewI3D_H CellType_WholeDomain(Kokkos::ViewAllocateWithoutInitializing("CellType_WholeDomain"), 0, 0, 0);
        ViewF3D_H UndercoolingChange_WholeDomain(
            Kokkos::ViewAllocateWithoutInitializing("UndercoolingChange_WholeDomain"), 0, 0, 0);
        ViewF3D_H UndercoolingCurrent_WholeDomain(
            Kokkos::ViewAllocateWithoutInitializing("UndercoolingCurrent_WholeDomain"), 0, 0, 0);

        if ((PrintTimeSeries) || (PrintMisorientation) || (PrintDebug == 2) || (PrintFullOutput) || (PrintDefaultRVE)) {
            // Collect all data on rank 0, for all data structures of interest
            ViewI_H GrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);
            Kokkos::resize(GrainID_WholeDomain, nz, nx, ny);
            CollectIntField(GrainID_WholeDomain, GrainID_Host, nz, MyXSlices, MyYSlices, np, RecvXOffset, RecvYOffset,
                            RecvXSlices, RecvYSlices, RBufSize);
        }
        if ((PrintMisorientation) || (PrintFullOutput) || (PrintDebug == 2) || (PrintFinalUndercoolingVals)) {
            Kokkos::resize(Melted_WholeDomain, nz, nx, ny);
            CollectBoolField(Melted_WholeDomain, Melted, nz, MyXSlices, MyYSlices, np, RecvXOffset, RecvYOffset,
                             RecvXSlices, RecvYSlices, RBufSize);
        }
        if ((PrintFullOutput) || (PrintDebug > 0) || (PrintDefaultRVE)) {
            ViewI_H LayerID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), LayerID);
            Kokkos::resize(LayerID_WholeDomain, nz, nx, ny);
            CollectIntField(LayerID_WholeDomain, LayerID_Host, nz, MyXSlices, MyYSlices, np, RecvXOffset, RecvYOffset,
                            RecvXSlices, RecvYSlices, RBufSize);
        }
        if ((PrintDebug > 0) || (PrintTimeSeries)) {
            ViewI_H CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
            Kokkos::resize(CellType_WholeDomain, nz, nx, ny);
            CollectIntField(CellType_WholeDomain, CellType_Host, nz, MyXSlices, MyYSlices, np, RecvXOffset, RecvYOffset,
                            RecvXSlices, RecvYSlices, RBufSize);
        }
        if (PrintDebug > 0) {
            ViewI_H CritTimeStep_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CritTimeStep);
            Kokkos::resize(CritTimeStep_WholeDomain, nz, nx, ny);
            CollectIntField(CritTimeStep_WholeDomain, CritTimeStep_Host, nz, MyXSlices, MyYSlices, np, RecvXOffset,
                            RecvYOffset, RecvXSlices, RecvYSlices, RBufSize);
        }
        if (PrintDebug == 2) {
            ViewF_H UndercoolingChange_Host =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), UndercoolingChange);
            Kokkos::resize(UndercoolingChange_WholeDomain, nz, nx, ny);
            CollectFloatField(UndercoolingChange_WholeDomain, UndercoolingChange_Host, nz, MyXSlices, MyYSlices, np,
                              RecvXOffset, RecvYOffset, RecvXSlices, RecvYSlices, RBufSize);
        }
        if ((PrintDebug == 2) || (PrintFinalUndercoolingVals)) {
            ViewF_H UndercoolingCurrent_Host =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), UndercoolingCurrent);
            Kokkos::resize(UndercoolingCurrent_WholeDomain, nz, nx, ny);
            CollectFloatField(UndercoolingCurrent_WholeDomain, UndercoolingCurrent_Host, nz, MyXSlices, MyYSlices, np,
                              RecvXOffset, RecvYOffset, RecvXSlices, RecvYSlices, RBufSize);
        }
        if ((PrintMisorientation) || (PrintTimeSeries)) {
            ViewF_H GrainUnitVector_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainUnitVector);
            if (PrintMisorientation)
                PrintGrainMisorientations(BaseFileName, PathToOutput, nx, ny, nz, Melted_WholeDomain,
                                          GrainID_WholeDomain, GrainUnitVector_Host, NGrainOrientations, deltax, XMin,
                                          YMin, ZMin);
            if (PrintTimeSeries)
                PrintIntermediateExaCAState(IntermediateFileCounter, layernumber, BaseFileName, PathToOutput,
                                            ZBound_Low, nzActive, nx, ny, GrainID_WholeDomain, CellType_WholeDomain,
                                            GrainUnitVector_Host, NGrainOrientations, deltax, XMin, YMin, ZMin);
        }
        if (PrintFinalUndercoolingVals)
            PrintFinalUndercooling(BaseFileName, PathToOutput, nx, ny, nz, Melted_WholeDomain,
                                   UndercoolingCurrent_WholeDomain, deltax, XMin, YMin, ZMin);
        if (PrintDefaultRVE)
            PrintExaConstitDefaultRVE(BaseFileName, PathToOutput, nx, ny, nz, LayerID_WholeDomain, GrainID_WholeDomain,
                                      deltax, NumberOfLayers);
        if ((PrintFullOutput) || (PrintDebug > 0))
            PrintCAFields(nx, ny, nz, GrainID_WholeDomain, LayerID_WholeDomain, CritTimeStep_WholeDomain,
                          CellType_WholeDomain, UndercoolingChange_WholeDomain, UndercoolingCurrent_WholeDomain,
                          Melted_WholeDomain, PathToOutput, BaseFileName, PrintDebug, PrintFullOutput, deltax, XMin,
                          YMin, ZMin);
    }
    else {

        // No ghost nodes in sent data:
        int SendBufStartX, SendBufEndX, SendBufStartY, SendBufEndY;
        if (MyYOffset == 0)
            SendBufStartY = 0;
        else
            SendBufStartY = 1;
        if (MyYSlices + MyYOffset == ny)
            SendBufEndY = MyYSlices;
        else
            SendBufEndY = MyYSlices - 1;
        if (MyXOffset == 0)
            SendBufStartX = 0;
        else
            SendBufStartX = 1;
        if (MyXSlices + MyXOffset == nx)
            SendBufEndX = MyXSlices;
        else
            SendBufEndX = MyXSlices - 1;

        int SendBufSize = (SendBufEndX - SendBufStartX) * (SendBufEndY - SendBufStartY) * nz;

        // Collect Melted/Grain ID data on rank 0
        if ((PrintTimeSeries) || (PrintMisorientation) || (PrintDebug == 2) || (PrintFullOutput) || (PrintDefaultRVE)) {
            ViewI_H GrainID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), GrainID);
            SendIntField(GrainID_Host, nz, MyXSlices, MyYSlices, SendBufSize, SendBufStartX, SendBufEndX, SendBufStartY,
                         SendBufEndY);
        }
        if ((PrintMisorientation) || (PrintFullOutput) || (PrintDebug == 2) || (PrintFinalUndercoolingVals))
            SendBoolField(Melted, nz, MyXSlices, MyYSlices, SendBufSize, SendBufStartX, SendBufEndX, SendBufStartY,
                          SendBufEndY);
        if ((PrintFullOutput) || (PrintDebug > 0) || (PrintDefaultRVE)) {
            ViewI_H LayerID_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), LayerID);
            SendIntField(LayerID_Host, nz, MyXSlices, MyYSlices, SendBufSize, SendBufStartX, SendBufEndX, SendBufStartY,
                         SendBufEndY);
        }
        if ((PrintDebug > 0) || (PrintTimeSeries)) {
            ViewI_H CellType_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CellType);
            SendIntField(CellType_Host, nz, MyXSlices, MyYSlices, SendBufSize, SendBufStartX, SendBufEndX,
                         SendBufStartY, SendBufEndY);
        }
        if (PrintDebug > 0) {
            ViewI_H CritTimeStep_Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), CritTimeStep);
            SendIntField(CritTimeStep_Host, nz, MyXSlices, MyYSlices, SendBufSize, SendBufStartX, SendBufEndX,
                         SendBufStartY, SendBufEndY);
        }
        if (PrintDebug == 2) {
            ViewF_H UndercoolingChange_Host =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), UndercoolingChange);
            SendFloatField(UndercoolingChange_Host, nz, MyXSlices, MyYSlices, SendBufSize, SendBufStartX, SendBufEndX,
                           SendBufStartY, SendBufEndY);
        }
        if ((PrintDebug == 2) || (PrintFinalUndercoolingVals)) {
            ViewF_H UndercoolingCurrent_Host =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), UndercoolingCurrent);
            SendFloatField(UndercoolingCurrent_Host, nz, MyXSlices, MyYSlices, SendBufSize, SendBufStartX, SendBufEndX,
                           SendBufStartY, SendBufEndY);
        }
    }
}

//*****************************************************************************/
// Print specified fields to a paraview file
void PrintCAFields(int nx, int ny, int nz, ViewI3D_H GrainID_WholeDomain, ViewI3D_H LayerID_WholeDomain,
                   ViewI3D_H CritTimeStep_WholeDomain, ViewI3D_H CellType_WholeDomain,
                   ViewF3D_H UndercoolingChange_WholeDomain, ViewF3D_H UndercoolingCurrent_WholeDomain,
                   ViewI3D_H Melted_WholeDomain, std::string PathToOutput, std::string BaseFileName, int PrintDebug,
                   bool PrintFullOutput, double deltax, float XMin, float YMin, float ZMin) {

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
    Grainplot.open(FName);
    Grainplot << "# vtk DataFile Version 3.0" << std::endl;
    Grainplot << "vtk output" << std::endl;
    Grainplot << "ASCII" << std::endl;
    Grainplot << "DATASET STRUCTURED_POINTS" << std::endl;
    Grainplot << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;
    Grainplot << "ORIGIN " << XMin << " " << YMin << " " << ZMin << std::endl;
    Grainplot << "SPACING " << deltax << " " << deltax << " " << deltax << std::endl;
    Grainplot << std::fixed << "POINT_DATA " << nx * ny * nz << std::endl;
    // Print Layer ID data
    Grainplot << "SCALARS LayerID int 1" << std::endl;
    Grainplot << "LOOKUP_TABLE default" << std::endl;
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                Grainplot << LayerID_WholeDomain(k, i, j) << " ";
            }
        }
        Grainplot << std::endl;
    }
    if ((PrintFullOutput) || (PrintDebug == 2)) {
        // Print Grain ID data
        Grainplot << "SCALARS GrainID int 1" << std::endl;
        Grainplot << "LOOKUP_TABLE default" << std::endl;
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    Grainplot << GrainID_WholeDomain(k, i, j) << " ";
                }
            }
            Grainplot << std::endl;
        }
        // Print Melted data
        Grainplot << "SCALARS Melted int 1" << std::endl;
        Grainplot << "LOOKUP_TABLE default" << std::endl;
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    Grainplot << Melted_WholeDomain(k, i, j) << " ";
                }
            }
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
                    Grainplot << CritTimeStep_WholeDomain(k, i, j) << " ";
                }
            }
            Grainplot << std::endl;
        }
        // Print CellType data
        Grainplot << "SCALARS CellType int 1" << std::endl;
        Grainplot << "LOOKUP_TABLE default" << std::endl;
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    Grainplot << CellType_WholeDomain(k, i, j) << " ";
                }
            }
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
                    Grainplot << UndercoolingChange_WholeDomain(k, i, j) << " ";
                }
            }
            Grainplot << std::endl;
        }
        // Print Undercooling current data
        Grainplot << "SCALARS Melted float 1" << std::endl;
        Grainplot << "LOOKUP_TABLE default" << std::endl;
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    Grainplot << UndercoolingCurrent_WholeDomain(k, i, j) << " ";
                }
            }
            Grainplot << std::endl;
        }
    }
    Grainplot.close();
}
//*****************************************************************************/
// Print grain misorientation, 0-62 for epitaxial grains and 100-162 for nucleated grains, to a paraview file
void PrintGrainMisorientations(std::string BaseFileName, std::string PathToOutput, int nx, int ny, int nz,
                               ViewI3D_H Melted_WholeDomain, ViewI3D_H GrainID_WholeDomain, ViewF_H GrainUnitVector,
                               int NGrainOrientations, double deltax, float XMin, float YMin, float ZMin) {

    std::string FName = PathToOutput + BaseFileName + "_Misorientations.vtk";
    std::cout << "Printing Paraview file of grain misorientations" << std::endl;
    // Print grain orientations to file
    std::ofstream GrainplotM;
    // New vtk file - don't print wall cells at +/- X/Y boundaries, nor at -Z boundary
    GrainplotM.open(FName);
    GrainplotM << "# vtk DataFile Version 3.0" << std::endl;
    GrainplotM << "vtk output" << std::endl;
    GrainplotM << "ASCII" << std::endl;
    GrainplotM << "DATASET STRUCTURED_POINTS" << std::endl;
    GrainplotM << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;
    GrainplotM << "ORIGIN " << XMin << " " << YMin << " " << ZMin << std::endl;
    GrainplotM << "SPACING " << deltax << " " << deltax << " " << deltax << std::endl;
    GrainplotM << std::fixed << "POINT_DATA " << nx * ny * nz << std::endl;
    GrainplotM << "SCALARS Angle_z int 1" << std::endl;
    GrainplotM << "LOOKUP_TABLE default" << std::endl;

    int NucleatedGrainCells = 0;
    int MeltedCells = 0;
    ViewI_H GrainMisorientation_Round(Kokkos::ViewAllocateWithoutInitializing("GrainMisorientation"),
                                      NGrainOrientations);
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
        GrainMisorientation_Round(n) = round(AngleZmin);
    }
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                if (Melted_WholeDomain(k, i, j) == 0)
                    GrainplotM << 200 << " ";
                else {
                    MeltedCells++;
                    int RoundedAngle;
                    int MyOrientation = getGrainOrientation(GrainID_WholeDomain(k, i, j), NGrainOrientations);
                    if (GrainID_WholeDomain(k, i, j) < 0) {
                        RoundedAngle = GrainMisorientation_Round(MyOrientation) + 100;
                        NucleatedGrainCells++;
                    }
                    else
                        RoundedAngle = GrainMisorientation_Round(MyOrientation);
                    GrainplotM << RoundedAngle << " ";
                }
            }
        }
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
                            ViewI3D_H Melted_WholeDomain, ViewF3D_H UndercoolingCurrent_WholeDomain, double deltax,
                            float XMin, float YMin, float ZMin) {

    std::string FName = PathToOutput + BaseFileName + "_FinalUndercooling.vtk";
    std::cout << "Printing Paraview file of final undercooling values" << std::endl;
    // Print undercooling to file
    std::ofstream UndercoolingPlot;
    UndercoolingPlot.open(FName);
    UndercoolingPlot << "# vtk DataFile Version 3.0" << std::endl;
    UndercoolingPlot << "vtk output" << std::endl;
    UndercoolingPlot << "ASCII" << std::endl;
    UndercoolingPlot << "DATASET STRUCTURED_POINTS" << std::endl;
    UndercoolingPlot << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;
    UndercoolingPlot << "ORIGIN " << XMin << " " << YMin << " " << ZMin << std::endl;
    UndercoolingPlot << "SPACING " << deltax << " " << deltax << " " << deltax << std::endl;
    UndercoolingPlot << std::fixed << "POINT_DATA " << nx * ny * nz << std::endl;
    UndercoolingPlot << "SCALARS UndercoolingFinal float 1" << std::endl;
    UndercoolingPlot << "LOOKUP_TABLE default" << std::endl;
    // Print 0s for the undercooling for cells that did not undergo solidification
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                if (Melted_WholeDomain(k, i, j) == 0)
                    UndercoolingPlot << 0.0 << " ";
                else
                    UndercoolingPlot << UndercoolingCurrent_WholeDomain(k, i, j) << " ";
            }
        }
        UndercoolingPlot << std::endl;
    }
    UndercoolingPlot.close();
}

//*****************************************************************************/
// Print the "default" representative volume element from this multilayer simulation, not including the final layer's
// microstructure, and centered in the simulation domain in X and Y Default RVE size is 0.5 by 0.5 by 0.5 mm
void PrintExaConstitDefaultRVE(std::string BaseFileName, std::string PathToOutput, int nx, int ny, int nz,
                               ViewI3D_H LayerID_WholeDomain, ViewI3D_H GrainID_WholeDomain, double deltax,
                               int NumberOfLayers) {

    // Determine the size, in CA cells of the RVE
    long int RVESize = std::lrint(0.0005 / deltax);

    // Determine the lower and upper Y bounds of the RVE
    long int RVE_XLow = std::floor(nx / 2) - std::floor(RVESize / 2);
    long int RVE_XHigh = RVE_XLow + RVESize - 1;
    long int RVE_YLow = std::floor(ny / 2) - std::floor(RVESize / 2);
    long int RVE_YHigh = RVE_YLow + RVESize - 1;

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
    long int RVE_ZHigh = nz - 1;
    for (int k = nz - 1; k >= 0; k--) {
        [&] {
            for (int i = RVE_XLow; i <= RVE_XHigh; i++) {
                for (int j = RVE_YLow; j <= RVE_YHigh; j++) {
                    if (LayerID_WholeDomain(k, i, j) == (NumberOfLayers - 1))
                        return; // check next lowest value for k
                }
            }
            RVE_ZHigh = k;
            k = 0; // leave loop
        }();
    }

    // Determine the lower Z bound of the RVE, and make sure the RVE fits in the simulation domain in X and Y
    long int RVE_ZLow = RVE_ZHigh - RVESize + 1;
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
        for (int i = RVE_XLow; i <= RVE_XHigh; i++) {
            for (int j = RVE_YLow; j <= RVE_YHigh; j++) {
                GrainplotE << i << "," << j << "," << k << "," << GrainID_WholeDomain(k, i, j) << std::endl;
            }
        }
    }
    GrainplotE.close();
}

//*****************************************************************************/
// Print a log file for this ExaCA run, containing information about the run parameters used
// from the input file as well as the decomposition scheme
void PrintExaCALog(int id, int np, std::string InputFile, std::string SimulationType, int DecompositionStrategy,
                   int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, double AConst, double BConst,
                   double CConst, double DConst, double FreezingRange, double deltax, double NMax, double dTN,
                   double dTsigma, std::string tempfile, int TempFilesInSeries, double HT_deltax, bool RemeltingYN,
                   double deltat, int NumberOfLayers, int LayerHeight, std::string SubstrateFileName,
                   double SubstrateGrainSpacing, bool SubstrateFile, double G, double R, int nx, int ny, int nz,
                   double FractSurfaceSitesActive, std::string PathToOutput, int NSpotsX, int NSpotsY, int SpotOffset,
                   int SpotRadius, double FractPowderSitesActive, std::string BaseFileName, double InitTime,
                   double RunTime, double OutTime, int cycle, double InitMaxTime, double InitMinTime,
                   double NuclMaxTime, double NuclMinTime, double CaptureMaxTime, double CaptureMinTime,
                   double GhostMaxTime, double GhostMinTime, double OutMaxTime, double OutMinTime) {

    int *XSlices = new int[np];
    int *YSlices = new int[np];
    int *XOffset = new int[np];
    int *YOffset = new int[np];
    MPI_Gather(&MyXSlices, 1, MPI_INT, XSlices, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&MyYSlices, 1, MPI_INT, YSlices, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&MyXOffset, 1, MPI_INT, XOffset, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&MyYOffset, 1, MPI_INT, YOffset, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (id == 0) {

        std::string FName = PathToOutput + BaseFileName + ".log";
        std::cout << "Printing ExaCA log file" << std::endl;
        std::ofstream ExaCALog;
        ExaCALog.open(FName);
        ExaCALog << "log file for a simulation run with input file " << InputFile << "  run on " << np
                 << " MPI ranks, output written at cycle " << cycle << std::endl;
        ExaCALog << "This simulation took " << InitTime + RunTime + OutTime
                 << " seconds to run, with the init/run/output breakdown as " << InitTime << "/" << RunTime << "/"
                 << OutTime << std::endl;
        ExaCALog << "This simulation was type: " << SimulationType << std::endl;
        ExaCALog << "Domain size in x: " << nx << std::endl;
        ExaCALog << "Domain size in y: " << ny << std::endl;
        ExaCALog << "Domain size in z: " << nz << std::endl;
        ExaCALog << "Cell size: " << deltax << " microns" << std::endl;
        ExaCALog << "Time step: " << deltat << " microseconds" << std::endl;
        ExaCALog << "Nucleation density was " << NMax << " m^-3 , mean nucleation undercooling was " << dTN
                 << " K, and standard deviation of nucleation undercooling was " << dTsigma << " K" << std::endl;
        ExaCALog << "Interfacial response function parameters used were " << AConst << " , " << BConst << " , "
                 << CConst << " , " << DConst << " , and the alloy freezing range was " << FreezingRange << std::endl;
        if (SimulationType == "C") {
            ExaCALog << "The thermal gradient was " << G << " K/m, and the cooling rate " << R << " K/s" << std::endl;
            ExaCALog << "The fraction of surface sites active was " << FractSurfaceSitesActive << std::endl;
        }
        else if ((SimulationType == "S") || (SimulationType == "SM")) {
            if (RemeltingYN)
                ExaCALog << "Remelting was included" << std::endl;
            else
                ExaCALog << "Remelting was not included" << std::endl;
            ExaCALog << "A total of " << NSpotsX << " in X and " << NSpotsY << " in Y were considered" << std::endl;
            ExaCALog << "The spots were offset by " << SpotOffset << " microns, and had radii of " << SpotRadius
                     << " microns" << std::endl;
            ExaCALog << "This pattern was repeated for " << NumberOfLayers << " layers" << std::endl;
            ExaCALog << "The fraction of active cells at the powder interface was " << FractPowderSitesActive
                     << std::endl;
        }
        else {
            ExaCALog << NumberOfLayers << " layers were simulated, with an offset of " << LayerHeight << " cells"
                     << std::endl;
            if (RemeltingYN)
                ExaCALog << "Remelting was included" << std::endl;
            else
                ExaCALog << "Remelting was not included" << std::endl;
            if (SubstrateFile)
                ExaCALog << "The substrate file was " << SubstrateFileName << std::endl;
            else
                ExaCALog << "The mean substrate grain size was " << SubstrateGrainSpacing << " microns" << std::endl;
            ExaCALog << "The " << TempFilesInSeries << " temperature file(s) used had the name " << tempfile
                     << std::endl;
            ExaCALog << "The temperature data resolution was " << HT_deltax << " microns" << std::endl;
            ExaCALog << "The fraction of active cells at the powder interface was " << FractPowderSitesActive
                     << std::endl;
        }
        ExaCALog << "The decomposition scheme used was: " << DecompositionStrategy << std::endl;
        for (int i = 0; i < np; i++) {
            ExaCALog << "Rank " << i << " contained " << XSlices[i] << " cells in x , " << YSlices[i]
                     << " cells in y; subdomain was offset by " << XOffset[i] << " in x , " << YOffset[i] << " in y"
                     << std::endl;
        }
        ExaCALog << "Max/min rank time initializing data  = " << InitMaxTime << " / " << InitMinTime << " s"
                 << std::endl;
        ExaCALog << "Max/min rank time in CA nucleation   = " << NuclMaxTime << " / " << NuclMinTime << " s"
                 << std::endl;
        ExaCALog << "Max/min rank time in CA cell capture = " << CaptureMaxTime << " / " << CaptureMinTime << " s"
                 << std::endl;
        ExaCALog << "Max/min rank time in CA ghosting     = " << GhostMaxTime << " / " << GhostMinTime << " s"
                 << std::endl;
        ExaCALog << "Max/min rank time exporting data     = " << OutMaxTime << " / " << OutMinTime << " s\n"
                 << std::endl;
        ExaCALog.close();
        // Also print this log information to the console
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

// Prints grain misorientation and cell type at some intermediate time during the simulation
// Files are named with the ending "_[FrameNumber].vtk", allowing Paraview to read these in as a series to make videos
// Cells that are liquid are given the value -1
// Epitaxial grains are colored 0-62 (by misorientation with +Z direction)
// Nucleated grains are colored 100-162 (by misorientation with +Z direction plus 100)
void PrintIntermediateExaCAState(int IntermediateFileCounter, int layernumber, std::string BaseFileName,
                                 std::string PathToOutput, int ZBound_Low, int nzActive, int nx, int ny,
                                 ViewI3D_H GrainID_WholeDomain, ViewI3D_H CellType_WholeDomain, ViewF_H GrainUnitVector,
                                 int NGrainOrientations, double deltax, float XMin, float YMin, float ZMin) {

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

    GrainplotM.open(FName);
    GrainplotM << "# vtk DataFile Version 3.0" << std::endl;
    GrainplotM << "vtk output" << std::endl;
    GrainplotM << "ASCII" << std::endl;
    GrainplotM << "DATASET STRUCTURED_POINTS" << std::endl;
    GrainplotM << "DIMENSIONS " << nx << " " << ny << " " << ZPrintSize << std::endl;
    GrainplotM << "ORIGIN " << XMin << " " << YMin << " " << ZMin << std::endl;
    GrainplotM << "SPACING " << deltax << " " << deltax << " " << deltax << std::endl;
    GrainplotM << std::fixed << "POINT_DATA " << nx * ny * ZPrintSize << std::endl;
    GrainplotM << "SCALARS Angle_z float 1" << std::endl;
    GrainplotM << "LOOKUP_TABLE default" << std::endl;
    for (int k = 0; k < ZPrintSize; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                // Liquid cells are assigned values of -1, epitaxial grains values between 0-63 (depending on
                // misorientation with Z direction), nucleated grains values between 100-163 (depending on
                // misorientation with Z direction, plus an offset of 100 to differentiate them from the epitaxial
                // grains)
                if (CellType_WholeDomain(k, i, j) != Liquid) {
                    // GrainID = 0 sites in the solid correspond to "inactive" parts of a powder layer
                    // i.e., they don't have an associated orientation. Assign these a unique value of 200
                    if (GrainID_WholeDomain(k, i, j) == 0)
                        GrainplotM << 200.0 << " ";
                    else {
                        int MyOrientation = getGrainOrientation(GrainID_WholeDomain(k, i, j), NGrainOrientations);
                        if (GrainID_WholeDomain(k, i, j) < 0)
                            GrainplotM << GrainMisorientation(MyOrientation) + 100.0 << " ";
                        else
                            GrainplotM << GrainMisorientation(MyOrientation) << " ";
                    }
                }
                else
                    GrainplotM << -1.0 << " ";
            }
        }
        GrainplotM << std::endl;
    }
    GrainplotM.close();
}
