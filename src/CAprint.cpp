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
        for (int i = 1; i < MyXSlices - 1; i++) {
            for (int j = 1; j < MyYSlices - 1; j++) {
                IntVar_WholeDomain(k, i - 1, j - 1) = IntVar(k * MyXSlices * MyYSlices + i * MyYSlices + j);
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
        for (int i = 1; i < MyXSlices - 1; i++) {
            for (int j = 1; j < MyYSlices - 1; j++) {
                FloatVar_WholeDomain(k, i - 1, j - 1) = FloatVar(k * MyXSlices * MyYSlices + i * MyYSlices + j);
            }
        }
    }

    // Recieve values from other ranks - message size different for different ranks
    for (int p = 1; p < np; p++) {
        int RecvBufSize_ThisRank = RBufSize(p);
        ViewI_H RecvBufFloatVar(Kokkos::ViewAllocateWithoutInitializing("RecvBufFloatVar"), RecvBufSize_ThisRank);
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
        for (int i = 1; i < MyXSlices - 1; i++) {
            for (int j = 1; j < MyYSlices - 1; j++) {
                if (BoolVar[k * MyXSlices * MyYSlices + i * MyYSlices + j])
                    IntVar_WholeDomain(k, i - 1, j - 1) = 1;
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
void SendIntField(ViewI_H VarToSend, int nz, int MyXSlices, int MyYSlices, int SendBufSize) {

    // Send non-ghost node data to rank 0
    int DataCounter = 0;
    ViewI_H SendBuf(Kokkos::ViewAllocateWithoutInitializing("SendBuf"), SendBufSize);
    for (int k = 0; k < nz; k++) {
        for (int i = 1; i < MyXSlices - 1; i++) {
            for (int j = 1; j < MyYSlices - 1; j++) {
                SendBuf(DataCounter) = VarToSend(k * MyXSlices * MyYSlices + i * MyYSlices + j);
                DataCounter++;
            }
        }
    }
    MPI_Send(SendBuf.data(), SendBufSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
}

// On rank > 0, send data for an float view to rank 0
void SendFloatField(ViewF_H VarToSend, int nz, int MyXSlices, int MyYSlices, int SendBufSize) {

    // Send non-ghost node data to rank 0
    int DataCounter = 0;
    ViewF_H SendBuf(Kokkos::ViewAllocateWithoutInitializing("SendBuf"), SendBufSize);
    for (int k = 0; k < nz; k++) {
        for (int i = 1; i < MyXSlices - 1; i++) {
            for (int j = 1; j < MyYSlices - 1; j++) {
                SendBuf(DataCounter) = VarToSend(k * MyXSlices * MyYSlices + i * MyYSlices + j);
                DataCounter++;
            }
        }
    }
    MPI_Send(SendBuf.data(), SendBufSize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
}

// On rank > 0, send data for a bool array (converted into integers for MPI) to rank 0
void SendBoolField(bool *VarToSend, int nz, int MyXSlices, int MyYSlices, int SendBufSize) {

    // Send non-ghost node data to rank 0
    int DataCounter = 0;
    ViewI_H SendBuf(Kokkos::ViewAllocateWithoutInitializing("SendBuf"), SendBufSize);
    for (int k = 0; k < nz; k++) {
        for (int i = 1; i < MyXSlices - 1; i++) {
            for (int j = 1; j < MyYSlices - 1; j++) {
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
                    int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI_H GrainID, ViewI_H GrainOrientation,
                    ViewI_H CritTimeStep, ViewF_H GrainUnitVector, ViewI_H LayerID, ViewI_H CellType,
                    ViewF_H UndercoolingChange, ViewF_H UndercoolingCurrent, std::string BaseFileName,
                    int DecompositionStrategy, int NGrainOrientations, bool *Melted, std::string PathToOutput,
                    int PrintDebug, bool PrintMisorientation, bool PrintFullOutput, bool PrintTimeSeries,
                    int IntermediateFileCounter, int ZBound_Low, int nzActive, double deltax, float XMin, float YMin,
                    float ZMin) {

    // Collect all data on rank 0, for all data structures of interest
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

            // No ghost nodes in send and recieved data
            RecvYSlices(p) = RecvYSlices(p) - 2;
            RecvXSlices(p) = RecvXSlices(p) - 2;

            // Readjust X and Y offsets of incoming data based on the lack of ghost nodes
            RecvYOffset(p)++;
            RecvXOffset(p)++;

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

        if ((PrintTimeSeries) || (PrintMisorientation) || (PrintDebug == 2) || (PrintFullOutput)) {
            Kokkos::resize(GrainID_WholeDomain, nz, nx, ny);
            CollectIntField(GrainID_WholeDomain, GrainID, nz, MyXSlices, MyYSlices, np, RecvXOffset, RecvYOffset,
                            RecvXSlices, RecvYSlices, RBufSize);
        }
        if ((PrintMisorientation) || (PrintFullOutput) || (PrintDebug == 2)) {
            Kokkos::resize(Melted_WholeDomain, nz, nx, ny);
            CollectBoolField(Melted_WholeDomain, Melted, nz, MyXSlices, MyYSlices, np, RecvXOffset, RecvYOffset,
                             RecvXSlices, RecvYSlices, RBufSize);
        }
        if ((PrintFullOutput) || (PrintDebug > 0)) {
            Kokkos::resize(LayerID_WholeDomain, nz, nx, ny);
            CollectIntField(LayerID_WholeDomain, LayerID, nz, MyXSlices, MyYSlices, np, RecvXOffset, RecvYOffset,
                            RecvXSlices, RecvYSlices, RBufSize);
        }
        if ((PrintDebug > 0) || (PrintTimeSeries)) {
            Kokkos::resize(CellType_WholeDomain, nz, nx, ny);
            CollectIntField(CellType_WholeDomain, CellType, nz, MyXSlices, MyYSlices, np, RecvXOffset, RecvYOffset,
                            RecvXSlices, RecvYSlices, RBufSize);
        }
        if (PrintDebug > 0) {
            Kokkos::resize(CritTimeStep_WholeDomain, nz, nx, ny);
            CollectIntField(CritTimeStep_WholeDomain, CritTimeStep, nz, MyXSlices, MyYSlices, np, RecvXOffset,
                            RecvYOffset, RecvXSlices, RecvYSlices, RBufSize);
        }
        if (PrintDebug == 2) {
            Kokkos::resize(UndercoolingChange_WholeDomain, nz, nx, ny);
            Kokkos::resize(UndercoolingCurrent_WholeDomain, nz, nx, ny);
            CollectFloatField(UndercoolingChange_WholeDomain, UndercoolingChange, nz, MyXSlices, MyYSlices, np,
                              RecvXOffset, RecvYOffset, RecvXSlices, RecvYSlices, RBufSize);
            CollectFloatField(UndercoolingCurrent_WholeDomain, UndercoolingCurrent, nz, MyXSlices, MyYSlices, np,
                              RecvXOffset, RecvYOffset, RecvXSlices, RecvYSlices, RBufSize);
        }

        if (PrintMisorientation)
            PrintGrainMisorientations(BaseFileName, PathToOutput, nx, ny, nz, Melted_WholeDomain, GrainID_WholeDomain,
                                      GrainOrientation, GrainUnitVector, NGrainOrientations, deltax, XMin, YMin, ZMin);
        if ((PrintFullOutput) || (PrintDebug > 0))
            PrintCAFields(nx, ny, nz, GrainID_WholeDomain, LayerID_WholeDomain, CritTimeStep_WholeDomain,
                          CellType_WholeDomain, UndercoolingChange_WholeDomain, UndercoolingCurrent_WholeDomain,
                          Melted_WholeDomain, PathToOutput, BaseFileName, PrintDebug, PrintFullOutput, deltax, XMin,
                          YMin, ZMin);
        if (PrintTimeSeries)
            PrintIntermediateExaCAState(IntermediateFileCounter, layernumber, BaseFileName, PathToOutput, ZBound_Low,
                                        nzActive, nx, ny, GrainID_WholeDomain, CellType_WholeDomain, GrainOrientation,
                                        GrainUnitVector, NGrainOrientations, deltax, XMin, YMin, ZMin);
    }
    else {
        int SendBufSize = (MyXSlices - 2) * (MyYSlices - 2) * nz;

        // Collect Melted/Grain ID data on rank 0
        if ((PrintTimeSeries) || (PrintMisorientation) || (PrintDebug == 2) || (PrintFullOutput))
            SendIntField(GrainID, nz, MyXSlices, MyYSlices, SendBufSize);
        if ((PrintMisorientation) || (PrintFullOutput) || (PrintDebug == 2))
            SendBoolField(Melted, nz, MyXSlices, MyYSlices, SendBufSize);
        if ((PrintFullOutput) || (PrintDebug > 0))
            SendIntField(LayerID, nz, MyXSlices, MyYSlices, SendBufSize);
        if ((PrintDebug > 0) || (PrintTimeSeries)) {
            SendIntField(CellType, nz, MyXSlices, MyYSlices, SendBufSize);
        }
        if (PrintDebug > 0) {
            SendIntField(CritTimeStep, nz, MyXSlices, MyYSlices, SendBufSize);
        }
        if (PrintDebug == 2) {
            SendFloatField(UndercoolingChange, nz, MyXSlices, MyYSlices, SendBufSize);
            SendFloatField(UndercoolingCurrent, nz, MyXSlices, MyYSlices, SendBufSize);
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
    // The CA domain is technically 2 cells bigger than the temperature field in all dimensions (except the top surface)
    // The start of the data in X, Y, Z is actually the XMin, YMin, ZMin from the temperature file minus 1 cell size
    Grainplot << "ORIGIN " << XMin - deltax << " " << YMin - deltax << " " << ZMin - deltax << std::endl;
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
                               ViewI3D_H Melted_WholeDomain, ViewI3D_H GrainID_WholeDomain, ViewI_H GrainOrientation,
                               ViewF_H GrainUnitVector, int NGrainOrientations, double deltax, float XMin, float YMin,
                               float ZMin) {

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
    GrainplotM << "DIMENSIONS " << nx - 2 << " " << ny - 2 << " " << nz - 1 << std::endl;
    // The CA domain is technically 2 cells bigger than the temperature field in all dimensions (except the top surface)
    // The start of the data in X, Y, Z is actually the XMin, YMin, ZMin from the temperature file minus 1 cell size
    GrainplotM << "ORIGIN " << XMin - deltax << " " << YMin - deltax << " " << ZMin - deltax << std::endl;
    GrainplotM << "SPACING " << deltax << " " << deltax << " " << deltax << std::endl;
    GrainplotM << std::fixed << "POINT_DATA " << (nx - 2) * (ny - 2) * (nz - 1) << std::endl;
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
    for (int k = 1; k < nz; k++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                if (Melted_WholeDomain(k, i, j) == 0)
                    GrainplotM << 200 << " ";
                else {
                    MeltedCells++;
                    int RoundedAngle;
                    int MyOrientation =
                        GrainOrientation(((abs(GrainID_WholeDomain(k, i, j)) - 1) % NGrainOrientations));
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
// Print a log file for this ExaCA run, containing information about the run parameters used
// from the input file as well as the decomposition scheme
void PrintExaCALog(int id, int np, std::string InputFile, std::string SimulationType, int DecompositionStrategy,
                   int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, double AConst, double BConst,
                   double CConst, double DConst, double FreezingRange, double deltax, double NMax, double dTN,
                   double dTsigma, std::string tempfile, int TempFilesInSeries, double HT_deltax, bool RemeltingYN,
                   double deltat, int NumberOfLayers, int LayerHeight, std::string SubstrateFileName,
                   double SubstrateGrainSpacing, bool SubstrateFile, double G, double R, int nx, int ny, int nz,
                   double FractSurfaceSitesActive, std::string PathToOutput, int NSpotsX, int NSpotsY, int SpotOffset,
                   int SpotRadius, std::string BaseFileName, double InitTime, double RunTime, double OutTime, int cycle,
                   double InitMaxTime, double InitMinTime, double NuclMaxTime, double NuclMinTime,
                   double CaptureMaxTime, double CaptureMinTime, double GhostMaxTime, double GhostMinTime,
                   double OutMaxTime, double OutMinTime) {

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
        else if (SimulationType == "S") {
            ExaCALog << "A total of " << NSpotsX << " in X and " << NSpotsY << " in Y were considered" << std::endl;
            ExaCALog << "The spots were offset by " << SpotOffset << " microns, and had radii of " << SpotRadius
                     << " microns" << std::endl;
            ExaCALog << "This pattern was repeated for " << NumberOfLayers << " layers" << std::endl;
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
            if (SimulationType == "R") {
                ExaCALog << "The " << TempFilesInSeries << " temperature file(s) used had the name " << tempfile
                         << std::endl;
                ExaCALog << "The temperature data resolution was " << HT_deltax << " microns" << std::endl;
            }
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
                                 ViewI3D_H GrainID_WholeDomain, ViewI3D_H CellType_WholeDomain,
                                 ViewI_H GrainOrientation, ViewF_H GrainUnitVector, int NGrainOrientations,
                                 double deltax, float XMin, float YMin, float ZMin) {

    std::string FName = PathToOutput + BaseFileName + "_layer" + std::to_string(layernumber) + "_" +
                        std::to_string(IntermediateFileCounter) + ".vtk";

    // Size of output in Z will depend on the current layer bounds - always start from Z = 1, but only collect data up
    // to Z = ZBound_Low + nzActive - 1
    int ZPrintSize = ZBound_Low + nzActive - 1;
    std::cout << "Printing file " << FName << " : Z coordinates of 1 through " << ZPrintSize << std::endl;
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
    GrainplotM << "DIMENSIONS " << nx - 2 << " " << ny - 2 << " " << ZPrintSize << std::endl;
    // The CA domain is technically 2 cells bigger than the temperature field in all dimensions (except the top surface)
    // The start of the data in X, Y, Z is actually the XMin, YMin, ZMin from the temperature file minus 1 cell size
    GrainplotM << "ORIGIN " << XMin - deltax << " " << YMin - deltax << " " << ZMin - deltax << std::endl;
    GrainplotM << "SPACING " << deltax << " " << deltax << " " << deltax << std::endl;
    GrainplotM << std::fixed << "POINT_DATA " << (nx - 2) * (ny - 2) * ZPrintSize << std::endl;
    GrainplotM << "SCALARS Angle_z float 1" << std::endl;
    GrainplotM << "LOOKUP_TABLE default" << std::endl;
    for (int k = 1; k <= ZPrintSize; k++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                if (CellType_WholeDomain(k, i, j) != Liquid) {
                    if (GrainID_WholeDomain(k, i, j) == 0)
                        std::cout << i << " " << j << " " << k << " " << CellType_WholeDomain(k, i, j) << std::endl;
                    int MyOrientation =
                        GrainOrientation(((abs(GrainID_WholeDomain(k, i, j)) - 1) % NGrainOrientations));
                    if (GrainID_WholeDomain(k, i, j) < 0)
                        GrainplotM << GrainMisorientation(MyOrientation) + 100.0 << " ";
                    else
                        GrainplotM << GrainMisorientation(MyOrientation) << " ";
                }
                else
                    GrainplotM << -1.0 << " ";
            }
        }
        GrainplotM << std::endl;
    }
    GrainplotM.close();
}
