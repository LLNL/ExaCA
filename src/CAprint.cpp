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
void CollectIntField(std::vector<std::vector<std::vector<int>>> &IntVar_WholeDomain, ViewI_H IntVar, int nx, int ny,
                     int nz, int MyXSlices, int MyYSlices, int np, int *RecvXOffset, int *RecvYOffset, int *RecvXSlices,
                     int *RecvYSlices, int *RBufSize) {

    // Resize int variable for whole domain and place values for rank 0
    for (int k = 0; k < nz; k++) {
        std::vector<std::vector<int>> IntVar_JK;
        for (int i = 0; i < nx; i++) {
            std::vector<int> IntVar_K;
            for (int j = 0; j < ny; j++) {
                IntVar_K.push_back(0);
            }
            IntVar_JK.push_back(IntVar_K);
        }
        IntVar_WholeDomain.push_back(IntVar_JK);
    }
    for (int k = 1; k < nz - 1; k++) {
        for (int i = 1; i < MyXSlices - 1; i++) {
            for (int j = 1; j < MyYSlices - 1; j++) {
                IntVar_WholeDomain[k][i - 1][j - 1] = IntVar(k * MyXSlices * MyYSlices + i * MyYSlices + j);
            }
        }
    }

    // Recieve values from other ranks - message size different for different ranks
    for (int p = 1; p < np; p++) {
        int RecvBufSize_ThisRank = RBufSize[p];
        int *RecvBufIntVar = new int[RecvBufSize_ThisRank];
        MPI_Recv(RecvBufIntVar, RecvBufSize_ThisRank, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int DataCounter = 0;
        for (int k = 1; k < nz - 1; k++) {
            for (int i = 0; i < RecvXSlices[p]; i++) {
                for (int j = 0; j < RecvYSlices[p]; j++) {
                    IntVar_WholeDomain[k][i + RecvXOffset[p]][j + RecvYOffset[p]] = RecvBufIntVar[DataCounter];
                    DataCounter++;
                }
            }
        }
    }
}

// On rank 0, collect data for one single float view
void CollectFloatField(std::vector<std::vector<std::vector<float>>> &FloatVar_WholeDomain, ViewF_H FloatVar, int nx,
                       int ny, int nz, int MyXSlices, int MyYSlices, int np, int *RecvXOffset, int *RecvYOffset,
                       int *RecvXSlices, int *RecvYSlices, int *RBufSize) {

    // Resize float variable for whole domain and place values for rank 0
    for (int k = 0; k < nz; k++) {
        std::vector<std::vector<float>> FloatVar_JK;
        for (int i = 0; i < nx; i++) {
            std::vector<float> FloatVar_K;
            for (int j = 0; j < ny; j++) {
                FloatVar_K.push_back(0);
            }
            FloatVar_JK.push_back(FloatVar_K);
        }
        FloatVar_WholeDomain.push_back(FloatVar_JK);
    }
    for (int k = 1; k < nz - 1; k++) {
        for (int i = 1; i < MyXSlices - 1; i++) {
            for (int j = 1; j < MyYSlices - 1; j++) {
                FloatVar_WholeDomain[k][i - 1][j - 1] = FloatVar(k * MyXSlices * MyYSlices + i * MyYSlices + j);
            }
        }
    }

    // Recieve values from other ranks - message size different for different ranks
    for (int p = 1; p < np; p++) {
        int RecvBufSize_ThisRank = RBufSize[p];
        float *RecvBufFloatVar = new float[RecvBufSize_ThisRank];
        MPI_Recv(RecvBufFloatVar, RecvBufSize_ThisRank, MPI_FLOAT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int DataCounter = 0;
        for (int k = 1; k < nz - 1; k++) {
            for (int i = 0; i < RecvXSlices[p]; i++) {
                for (int j = 0; j < RecvYSlices[p]; j++) {
                    FloatVar_WholeDomain[k][i + RecvXOffset[p]][j + RecvYOffset[p]] = RecvBufFloatVar[DataCounter];
                    DataCounter++;
                }
            }
        }
    }
}

// On rank 0, collect data for one single bool array (and convert it to integer 0s and 1s for MPI)
void CollectBoolField(std::vector<std::vector<std::vector<int>>> &IntVar_WholeDomain, bool *BoolVar, int nx, int ny,
                      int nz, int MyXSlices, int MyYSlices, int np, int *RecvXOffset, int *RecvYOffset,
                      int *RecvXSlices, int *RecvYSlices, int *RBufSize) {

    // Resize bool variable for whole domain and place values for rank 0
    for (int k = 0; k < nz; k++) {
        std::vector<std::vector<int>> IntVar_JK;
        for (int i = 0; i < nx; i++) {
            std::vector<int> IntVar_K;
            for (int j = 0; j < ny; j++) {
                IntVar_K.push_back(0);
            }
            IntVar_JK.push_back(IntVar_K);
        }
        IntVar_WholeDomain.push_back(IntVar_JK);
    }
    for (int k = 1; k < nz - 1; k++) {
        for (int i = 1; i < MyXSlices - 1; i++) {
            for (int j = 1; j < MyYSlices - 1; j++) {
                if (BoolVar[k * MyXSlices * MyYSlices + i * MyYSlices + j])
                    IntVar_WholeDomain[k][i - 1][j - 1] = 1;
                else
                    IntVar_WholeDomain[k][i - 1][j - 1] = 0;
            }
        }
    }

    // Recieve values from other ranks - message size different for different ranks
    for (int p = 1; p < np; p++) {
        int RecvBufSize_ThisRank = RBufSize[p];
        int *RecvBufIntVar = new int[RecvBufSize_ThisRank];
        MPI_Recv(RecvBufIntVar, RecvBufSize_ThisRank, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int DataCounter = 0;
        for (int k = 1; k < nz - 1; k++) {
            for (int i = 0; i < RecvXSlices[p]; i++) {
                for (int j = 0; j < RecvYSlices[p]; j++) {
                    IntVar_WholeDomain[k][i + RecvXOffset[p]][j + RecvYOffset[p]] = RecvBufIntVar[DataCounter];
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
    int *SendBuf = new int[SendBufSize];
    for (int k = 1; k < nz - 1; k++) {
        for (int i = 1; i < MyXSlices - 1; i++) {
            for (int j = 1; j < MyYSlices - 1; j++) {
                SendBuf[DataCounter] = VarToSend(k * MyXSlices * MyYSlices + i * MyYSlices + j);
                DataCounter++;
            }
        }
    }
    MPI_Send(SendBuf, SendBufSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
}

// On rank > 0, send data for an float view to rank 0
void SendFloatField(ViewF_H VarToSend, int nz, int MyXSlices, int MyYSlices, int SendBufSize) {

    // Send non-ghost node data to rank 0
    int DataCounter = 0;
    float *SendBuf = new float[SendBufSize];
    for (int k = 1; k < nz - 1; k++) {
        for (int i = 1; i < MyXSlices - 1; i++) {
            for (int j = 1; j < MyYSlices - 1; j++) {
                SendBuf[DataCounter] = VarToSend(k * MyXSlices * MyYSlices + i * MyYSlices + j);
                DataCounter++;
            }
        }
    }
    MPI_Send(SendBuf, SendBufSize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
}

// On rank > 0, send data for a bool array (converted into integers for MPI) to rank 0
void SendBoolField(bool *VarToSend, int nz, int MyXSlices, int MyYSlices, int SendBufSize) {

    // Send non-ghost node data to rank 0
    int DataCounter = 0;
    int *SendBuf = new int[SendBufSize];
    for (int k = 1; k < nz - 1; k++) {
        for (int i = 1; i < MyXSlices - 1; i++) {
            for (int j = 1; j < MyYSlices - 1; j++) {
                if (VarToSend[k * MyXSlices * MyYSlices + i * MyYSlices + j])
                    SendBuf[DataCounter] = 1;
                else
                    SendBuf[DataCounter] = 0;
                DataCounter++;
            }
        }
    }
    MPI_Send(SendBuf, SendBufSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
}

//*****************************************************************************/
// Prints values of selected data structures to Paraview files
void PrintExaCAData(int id, int np, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection,
                    int ProcessorsInYDirection, ViewI_H GrainID, ViewI_H GrainOrientation, ViewI_H CritTimeStep,
                    ViewF_H GrainUnitVector, ViewI_H LayerID, ViewI_H CellType, ViewF_H UndercoolingChange,
                    ViewF_H UndercoolingCurrent, std::string BaseFileName, int DecompositionStrategy,
                    int NGrainOrientations, bool *Melted, std::string PathToOutput, int PrintDebug,
                    bool PrintMisorientation, bool PrintFullOutput) {

    // Collect all data on rank 0, for all data structures of interest
    if (id == 0) {
        // Message sizes and data offsets for data recieved from other ranks- message size different for different ranks
        int *RecvXOffset = new int[np];
        int *RecvYOffset = new int[np];
        int *RecvXSlices = new int[np];
        int *RecvYSlices = new int[np];
        int *RBufSize = new int[np];

        for (int p = 1; p < np; p++) {
            RecvXOffset[p] = XOffsetCalc(p, nx, ProcessorsInXDirection, ProcessorsInYDirection, DecompositionStrategy);
            RecvXSlices[p] =
                XMPSlicesCalc(p, nx, ProcessorsInXDirection, ProcessorsInYDirection, DecompositionStrategy);

            RecvYOffset[p] = YOffsetCalc(p, ny, ProcessorsInYDirection, np, DecompositionStrategy);
            RecvYSlices[p] = YMPSlicesCalc(p, ny, ProcessorsInYDirection, np, DecompositionStrategy);

            // No ghost nodes in send and recieved data
            RecvYSlices[p] = RecvYSlices[p] - 2;
            RecvXSlices[p] = RecvXSlices[p] - 2;

            // Readjust X and Y offsets of incoming data based on the lack of ghost nodes
            RecvYOffset[p]++;
            RecvXOffset[p]++;

            RBufSize[p] = RecvXSlices[p] * RecvYSlices[p] * (nz - 2);
        }
        // Create variables for each possible data structure being collected on rank 0
        std::vector<std::vector<std::vector<int>>> GrainID_WholeDomain, LayerID_WholeDomain, CellType_WholeDomain,
            CritTimeStep_WholeDomain;
        std::vector<std::vector<std::vector<float>>> UndercoolingChange_WholeDomain, UndercoolingCurrent_WholeDomain;
        std::vector<std::vector<std::vector<int>>> Melted_WholeDomain;

        // If PrintDebug = 0, we aren't printing any debug files after initialization, but we are printing either the
        // misorientations (requiring Melted and GrainID) or we are printing the full end-of-run dataset (Melted,
        // GrainID, LayerID) If PrintDebug = 1, we are plotting CellType, LayerID, and CritTimeStep If PrintDebug = 2,
        // we are plotting all possible data structures

        if ((PrintMisorientation) || (PrintDebug == 2) || (PrintFullOutput))
            CollectIntField(GrainID_WholeDomain, GrainID, nx, ny, nz, MyXSlices, MyYSlices, np, RecvXOffset,
                            RecvYOffset, RecvXSlices, RecvYSlices, RBufSize);
        if ((PrintMisorientation) || (PrintFullOutput) || (PrintDebug == 2))
            CollectBoolField(Melted_WholeDomain, Melted, nx, ny, nz, MyXSlices, MyYSlices, np, RecvXOffset, RecvYOffset,
                             RecvXSlices, RecvYSlices, RBufSize);
        if ((PrintFullOutput) || (PrintDebug > 0))
            CollectIntField(LayerID_WholeDomain, LayerID, nx, ny, nz, MyXSlices, MyYSlices, np, RecvXOffset,
                            RecvYOffset, RecvXSlices, RecvYSlices, RBufSize);
        if (PrintDebug > 0) {
            CollectIntField(CellType_WholeDomain, CellType, nx, ny, nz, MyXSlices, MyYSlices, np, RecvXOffset,
                            RecvYOffset, RecvXSlices, RecvYSlices, RBufSize);
            CollectIntField(CritTimeStep_WholeDomain, CritTimeStep, nx, ny, nz, MyXSlices, MyYSlices, np, RecvXOffset,
                            RecvYOffset, RecvXSlices, RecvYSlices, RBufSize);
        }
        if (PrintDebug == 2) {
            CollectFloatField(UndercoolingChange_WholeDomain, UndercoolingChange, nx, ny, nz, MyXSlices, MyYSlices, np,
                              RecvXOffset, RecvYOffset, RecvXSlices, RecvYSlices, RBufSize);
            CollectFloatField(UndercoolingCurrent_WholeDomain, UndercoolingCurrent, nx, ny, nz, MyXSlices, MyYSlices,
                              np, RecvXOffset, RecvYOffset, RecvXSlices, RecvYSlices, RBufSize);
        }

        if (PrintMisorientation)
            PrintParaview(BaseFileName, PathToOutput, nx, ny, nz, Melted_WholeDomain, GrainID_WholeDomain,
                          GrainOrientation, GrainUnitVector, NGrainOrientations);
        PrintParaviewGeneric(nx, ny, nz, GrainID_WholeDomain, LayerID_WholeDomain, CritTimeStep_WholeDomain,
                             CellType_WholeDomain, UndercoolingChange_WholeDomain, UndercoolingCurrent_WholeDomain,
                             Melted_WholeDomain, PathToOutput, BaseFileName, PrintDebug, PrintFullOutput);
    }
    else {
        int SendBufSize = (MyXSlices - 2) * (MyYSlices - 2) * (nz - 2);

        // Collect Melted/Grain ID data on rank 0
        if ((PrintMisorientation) || (PrintDebug == 2) || (PrintFullOutput))
            SendIntField(GrainID, nz, MyXSlices, MyYSlices, SendBufSize);
        if ((PrintMisorientation) || (PrintFullOutput) || (PrintDebug == 2))
            SendBoolField(Melted, nz, MyXSlices, MyYSlices, SendBufSize);
        if ((PrintFullOutput) || (PrintDebug > 0))
            SendIntField(LayerID, nz, MyXSlices, MyYSlices, SendBufSize);
        if (PrintDebug > 0) {
            SendIntField(CellType, nz, MyXSlices, MyYSlices, SendBufSize);
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
void PrintParaviewGeneric(int nx, int ny, int nz, std::vector<std::vector<std::vector<int>>> GrainID_WholeDomain,
                          std::vector<std::vector<std::vector<int>>> LayerID_WholeDomain,
                          std::vector<std::vector<std::vector<int>>> CritTimeStep_WholeDomain,
                          std::vector<std::vector<std::vector<int>>> CellType_WholeDomain,
                          std::vector<std::vector<std::vector<float>>> UndercoolingChange_WholeDomain,
                          std::vector<std::vector<std::vector<float>>> UndercoolingCurrent_WholeDomain,
                          std::vector<std::vector<std::vector<int>>> Melted_WholeDomain, std::string PathToOutput,
                          std::string BaseFileName, int PrintDebug, bool PrintFullOutput) {

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
    Grainplot << "ORIGIN 0 0 0" << std::endl;
    Grainplot << "SPACING 1 1 1" << std::endl;
    Grainplot << std::fixed << "POINT_DATA " << nx * ny * nz << std::endl;
    // Print Layer ID data
    Grainplot << "SCALARS LayerID int 1" << std::endl;
    Grainplot << "LOOKUP_TABLE default" << std::endl;
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                Grainplot << LayerID_WholeDomain[k][i][j] << " ";
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
                    Grainplot << GrainID_WholeDomain[k][i][j] << " ";
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
                    Grainplot << Melted_WholeDomain[k][i][j] << " ";
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
                    Grainplot << CritTimeStep_WholeDomain[k][i][j] << " ";
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
                    Grainplot << CellType_WholeDomain[k][i][j] << " ";
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
                    Grainplot << UndercoolingChange_WholeDomain[k][i][j] << " ";
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
                    Grainplot << UndercoolingCurrent_WholeDomain[k][i][j] << " ";
                }
            }
            Grainplot << std::endl;
        }
    }
    Grainplot.close();
}
//*****************************************************************************/
// Print grain misorientation, 0-62 for epitaxial grains and 100-162 for nucleated grains, to a paraview file
void PrintParaview(std::string BaseFileName, std::string PathToOutput, int nx, int ny, int nz,
                   std::vector<std::vector<std::vector<int>>> Melted_WholeDomain,
                   std::vector<std::vector<std::vector<int>>> GrainID_WholeDomain, ViewI_H GrainOrientation,
                   ViewF_H GrainUnitVector, int NGrainOrientations) {

    std::string FName = PathToOutput + BaseFileName + "_Misorientations.vtk";
    std::cout << "Printing Paraview file of grain misorientations" << std::endl;
    // Print grain orientations to file
    std::ofstream GrainplotM;
    // New vtk file
    GrainplotM.open(FName);
    GrainplotM << "# vtk DataFile Version 3.0" << std::endl;
    GrainplotM << "vtk output" << std::endl;
    GrainplotM << "ASCII" << std::endl;
    GrainplotM << "DATASET STRUCTURED_POINTS" << std::endl;
    GrainplotM << "DIMENSIONS " << nx - 2 << " " << ny - 2 << " " << nz - 2 << std::endl;
    GrainplotM << "ORIGIN 0 0 0" << std::endl;
    GrainplotM << "SPACING 1 1 1" << std::endl;
    GrainplotM << std::fixed << "POINT_DATA " << (nx - 2) * (ny - 2) * (nz - 2) << std::endl;
    GrainplotM << "SCALARS Angle_z int 1" << std::endl;
    GrainplotM << "LOOKUP_TABLE default" << std::endl;

    int NucleatedGrainCells = 0;
    int MeltedCells = 0;
    for (int k = 1; k < nz - 1; k++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                if (Melted_WholeDomain[k][i][j] == 0)
                    GrainplotM << 200 << " ";
                else {
                    MeltedCells++;
                    int RoundedAngle;
                    int MyOrientation =
                        GrainOrientation(((abs(GrainID_WholeDomain[k][i][j]) - 1) % NGrainOrientations));
                    double AngleZmin = 62.7;
                    for (int ll = 0; ll < 3; ll++) {
                        double AngleZ = std::abs((180 / M_PI) * acos(GrainUnitVector(9 * MyOrientation + 3 * ll + 2)));
                        if (AngleZ < AngleZmin) {
                            AngleZmin = AngleZ;
                        }
                    }
                    if (GrainID_WholeDomain[k][i][j] < 0) {
                        RoundedAngle = round(AngleZmin) + 100;
                        NucleatedGrainCells++;
                    }
                    else
                        RoundedAngle = round(AngleZmin);
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
                   int SpotRadius, std::string BaseFileName, double InitTime, double RunTime, double OutTime) {

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
        ExaCALog << "log file for a simulation run with input file " << InputFile << "  run on " << np << " MPI ranks"
                 << std::endl;
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
        ExaCALog.close();
    }
}
