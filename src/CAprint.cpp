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
// Prints values of grain orientation for all cells to files
void CollectGrainData(int id, int np, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection,
                      int ProcessorsInYDirection, ViewI_H GrainID, ViewI_H GrainOrientation, ViewF_H GrainUnitVector,
                      std::string BaseFileName, int DecompositionStrategy, int NGrainOrientations, bool *Melted,
                      std::string PathToOutput, bool FilesToPrint[4], double deltax) {

    if (id == 0) {
        // Create GrainID variable for entire domain, place GrainIDs for rank 0
        std::vector<std::vector<std::vector<int>>> GrainID_WholeDomain;
        std::vector<std::vector<std::vector<bool>>> Melted_WholeDomain;
        for (int k = 0; k < nz; k++) {
            std::vector<std::vector<int>> GrainID_JK;
            std::vector<std::vector<bool>> Melted_JK;
            for (int i = 0; i < nx; i++) {
                std::vector<int> GrainID_K;
                std::vector<bool> Melted_K;
                for (int j = 0; j < ny; j++) {
                    GrainID_K.push_back(0);
                    Melted_K.push_back(false);
                }
                GrainID_JK.push_back(GrainID_K);
                Melted_JK.push_back(Melted_K);
            }
            GrainID_WholeDomain.push_back(GrainID_JK);
            Melted_WholeDomain.push_back(Melted_JK);
        }
        for (int k = 1; k < nz - 1; k++) {
            for (int i = 1; i < MyXSlices - 1; i++) {
                for (int j = 1; j < MyYSlices - 1; j++) {
                    GrainID_WholeDomain[k][i - 1][j - 1] = GrainID(k * MyXSlices * MyYSlices + i * MyYSlices + j);
                    Melted_WholeDomain[k][i - 1][j - 1] = Melted[k * MyXSlices * MyYSlices + i * MyYSlices + j];
                }
            }
        }

        // Recieve Grain ID value from other ranks
        // Message size different for different ranks
        for (int p = 1; p < np; p++) {
            int RecvXOffset = XOffsetCalc(p, nx, ProcessorsInXDirection, ProcessorsInYDirection, DecompositionStrategy);
            int RecvXSlices =
                XMPSlicesCalc(p, nx, ProcessorsInXDirection, ProcessorsInYDirection, DecompositionStrategy);

            int RecvYOffset = YOffsetCalc(p, ny, ProcessorsInYDirection, np, DecompositionStrategy);
            int RecvYSlices = YMPSlicesCalc(p, ny, ProcessorsInYDirection, np, DecompositionStrategy);

            // No ghost nodes in send and recieved data
            RecvYSlices = RecvYSlices - 2;
            RecvXSlices = RecvXSlices - 2;

            // Readjust X and Y offsets of incoming data based on the lack of ghost nodes
            RecvYOffset++;
            RecvXOffset++;

            int RBufSize = (RecvXSlices) * (RecvYSlices) * (nz - 2);
            int *RecvBufGID = new int[RBufSize];
            int *RecvBufM = new int[RBufSize];
            int DataCounter = 0;
            MPI_Recv(RecvBufGID, RBufSize, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(RecvBufM, RBufSize, MPI_INT, p, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (int k = 1; k < nz - 1; k++) {
                for (int i = 0; i < RecvXSlices; i++) {
                    for (int j = 0; j < RecvYSlices; j++) {
                        GrainID_WholeDomain[k][i + RecvXOffset][j + RecvYOffset] = RecvBufGID[DataCounter];
                        if (RecvBufM[DataCounter] == 1)
                            Melted_WholeDomain[k][i + RecvXOffset][j + RecvYOffset] = true;
                        else
                            Melted_WholeDomain[k][i + RecvXOffset][j + RecvYOffset] = false;
                        DataCounter++;
                    }
                }
            }
        }

        std::string FName = PathToOutput + BaseFileName;

        if (FilesToPrint[0])
            PrintOrientations(FName, nx, ny, nz, Melted_WholeDomain, GrainID_WholeDomain, NGrainOrientations);
        if (FilesToPrint[1])
            PrintGrainIDs(FName, nx, ny, nz, Melted_WholeDomain, GrainID_WholeDomain);
        if (FilesToPrint[2])
            PrintGrainIDsForExaConstit(FName, nx, ny, nz, GrainID_WholeDomain, deltax);
        if (FilesToPrint[3])
            PrintParaview(FName, nx, ny, nz, Melted_WholeDomain, GrainID_WholeDomain, GrainOrientation, GrainUnitVector,
                          NGrainOrientations);
        if (FilesToPrint[4])
            PrintGrainAreas(FName, deltax, nx, ny, nz, GrainID_WholeDomain);
        if (FilesToPrint[5])
            PrintWeightedGrainAreas(FName, deltax, nx, ny, nz, GrainID_WholeDomain);
    }
    else {

        // Send non-ghost node data to rank 0
        int SBufSize = (MyXSlices - 2) * (MyYSlices - 2) * (nz - 2);
        int *SendBufGID = new int[SBufSize];
        int *SendBufM = new int[SBufSize];
        int DataCounter = 0;

        for (int k = 1; k < nz - 1; k++) {
            for (int i = 1; i < MyXSlices - 1; i++) {
                for (int j = 1; j < MyYSlices - 1; j++) {
                    SendBufGID[DataCounter] = GrainID(k * MyXSlices * MyYSlices + i * MyYSlices + j);
                    if (Melted[k * MyXSlices * MyYSlices + i * MyYSlices + j])
                        SendBufM[DataCounter] = 1;
                    else
                        SendBufM[DataCounter] = 0;
                    DataCounter++;
                }
            }
        }

        MPI_Send(SendBufGID, SBufSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(SendBufM, SBufSize, MPI_INT, 0, 1, MPI_COMM_WORLD);
    }
}

//*****************************************************************************/
// Prints values of critical undercooling for all cells to files
void PrintTempValues(int id, int np, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection,
                     int ProcessorsInYDirection, ViewI_H CritTimeStep, int DecompositionStrategy,
                     std::string PathToOutput) {

    // Critical time step for solidification start printed to file first
    if (id == 0) {
        // Create GrainID variable for entire domain, place GrainIDs for rank 0
        std::vector<std::vector<std::vector<int>>> CritTimeStep_WholeDomain;
        for (int k = 0; k < nz; k++) {
            std::vector<std::vector<int>> CTS_JK;
            for (int i = 0; i < nx; i++) {
                std::vector<int> CTS_K;
                for (int j = 0; j < ny; j++) {
                    CTS_K.push_back(0);
                }
                CTS_JK.push_back(CTS_K);
            }
            CritTimeStep_WholeDomain.push_back(CTS_JK);
        }
        for (int k = 1; k < nz - 1; k++) {
            for (int i = 1; i < MyXSlices - 1; i++) {
                for (int j = 1; j < MyYSlices - 1; j++) {
                    CritTimeStep_WholeDomain[k][i - 1][j - 1] =
                        CritTimeStep(k * MyXSlices * MyYSlices + i * MyYSlices + j);
                }
            }
        }

        // Recieve values from other ranks
        // Message size different for different ranks
        for (int p = 1; p < np; p++) {
            int RecvXOffset = XOffsetCalc(p, nx, ProcessorsInXDirection, ProcessorsInYDirection, DecompositionStrategy);
            int RecvXSlices =
                XMPSlicesCalc(p, nx, ProcessorsInXDirection, ProcessorsInYDirection, DecompositionStrategy);

            int RecvYOffset = YOffsetCalc(p, ny, ProcessorsInYDirection, np, DecompositionStrategy);
            int RecvYSlices = YMPSlicesCalc(p, ny, ProcessorsInYDirection, np, DecompositionStrategy);

            // No ghost nodes in send and recieved data
            RecvYSlices = RecvYSlices - 2;
            RecvXSlices = RecvXSlices - 2;

            // Readjust X and Y offsets of incoming data based on the lack of ghost nodes
            RecvYOffset++;
            RecvXOffset++;

            int RBufSize = RecvXSlices * RecvYSlices * (nz - 2);
            int *RecvBufGID = new int[RBufSize];
            MPI_Recv(RecvBufGID, RBufSize, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            int DataCounter = 0;
            for (int k = 1; k < nz - 1; k++) {
                for (int i = 0; i < RecvXSlices; i++) {
                    for (int j = 0; j < RecvYSlices; j++) {
                        CritTimeStep_WholeDomain[k][i + RecvXOffset][j + RecvYOffset] = RecvBufGID[DataCounter];
                        DataCounter++;
                    }
                }
            }
        }

        // Print grain orientations to file
        std::string FName = PathToOutput + "CriticalTimeStep.vtk";
        std::ofstream CTSplot;
        // New vtk file
        CTSplot.open(FName);
        CTSplot << "# vtk DataFile Version 3.0" << std::endl;
        CTSplot << "vtk output" << std::endl;
        CTSplot << "ASCII" << std::endl;
        CTSplot << "DATASET STRUCTURED_POINTS" << std::endl;
        CTSplot << "DIMENSIONS " << nx - 2 << " " << ny - 2 << " " << nz - 2 << std::endl;
        CTSplot << "ORIGIN 0 0 0" << std::endl;
        CTSplot << "SPACING 1 1 1" << std::endl;
        CTSplot << std::fixed << "POINT_DATA " << (nx - 2) * (ny - 2) * (nz - 2) << std::endl;
        CTSplot << "SCALARS Angle_z float 1" << std::endl;
        CTSplot << "LOOKUP_TABLE default" << std::endl;

        for (int k = 1; k < nz - 1; k++) {
            for (int j = 1; j < ny - 1; j++) {
                for (int i = 1; i < nx - 1; i++) {
                    CTSplot << 0.75 * pow(10, -6) * (float)(CritTimeStep_WholeDomain[k][i][j]) << " ";
                }
            }
            CTSplot << std::endl;
        }
        CTSplot.close();
    }
    else {
        // Send non-ghost node data to rank 0
        int SBufSize = (MyXSlices - 2) * (MyYSlices - 2) * (nz - 2);
        int *SendBufGID = new int[SBufSize];
        int DataCounter = 0;
        for (int k = 1; k < nz - 1; k++) {
            for (int i = 1; i < MyXSlices - 1; i++) {
                for (int j = 1; j < MyYSlices - 1; j++) {
                    SendBufGID[DataCounter] = CritTimeStep(k * MyXSlices * MyYSlices + i * MyYSlices + j);
                    DataCounter++;
                }
            }
        }
        MPI_Send(SendBufGID, SBufSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}

//*****************************************************************************/
void PrintCT(int id, int np, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection,
             int ProcessorsInYDirection, ViewI_H CellType, std::string BaseFileName, int DecompositionStrategy) {

    std::string FName = BaseFileName + "_CT.vtk";
    std::ofstream Grainplot;
    if (id == 0) {
        // Create GrainID variable for entire domain, place GrainIDs for rank 0
        std::vector<std::vector<std::vector<int>>> CellType_WholeDomain;
        for (int k = 0; k < nz; k++) {
            std::vector<std::vector<int>> CellType_JK;
            for (int i = 0; i < nx; i++) {
                std::vector<int> CellType_K;
                for (int j = 0; j < ny; j++) {
                    CellType_K.push_back(Wall);
                }
                CellType_JK.push_back(CellType_K);
            }
            CellType_WholeDomain.push_back(CellType_JK);
        }
        for (int k = 1; k < nz - 1; k++) {
            for (int i = 1; i < MyXSlices - 1; i++) {
                for (int j = 1; j < MyYSlices - 1; j++) {
                    CellType_WholeDomain[k][i - 1][j - 1] = CellType(k * MyXSlices * MyYSlices + i * MyYSlices + j);
                }
            }
        }

        std::cout << "Domain created: Rank 0 placed slices 0 through" << MyYSlices - 2 << std::endl;
        // Recieve Grain ID value from other ranks
        // Message size different for different ranks
        for (int p = 1; p < np; p++) {
            int RecvXOffset = XOffsetCalc(p, nx, ProcessorsInXDirection, ProcessorsInYDirection, DecompositionStrategy);
            int RecvXSlices =
                XMPSlicesCalc(p, nx, ProcessorsInXDirection, ProcessorsInYDirection, DecompositionStrategy);

            int RecvYOffset = YOffsetCalc(p, ny, ProcessorsInYDirection, np, DecompositionStrategy);
            int RecvYSlices = YMPSlicesCalc(p, ny, ProcessorsInYDirection, np, DecompositionStrategy);

            // No ghost nodes in send and recieved data
            RecvYSlices = RecvYSlices - 2;
            RecvXSlices = RecvXSlices - 2;

            // Readjust X and Y offsets of incoming data based on the lack of ghost nodes
            RecvYOffset++;
            RecvXOffset++;

            int RBufSize = RecvXSlices * RecvYSlices * (nz - 2);
            int *RecvBufGID = new int[RBufSize];
            MPI_Recv(RecvBufGID, RBufSize, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            int DataCounter = 0;

            for (int k = 1; k < nz - 1; k++) {
                for (int i = 0; i < RecvXSlices; i++) {
                    for (int j = 0; j < RecvYSlices; j++) {
                        CellType_WholeDomain[k][i + RecvXOffset][j + RecvYOffset] = RecvBufGID[DataCounter];
                        DataCounter++;
                    }
                }
            }
        }
        std::cout << "Opening file" << std::endl;
        // Print grain orientations to file
        std::ofstream Grainplot;
        // New vtk file
        Grainplot.open(FName);
        Grainplot << "# vtk DataFile Version 3.0" << std::endl;
        Grainplot << "vtk output" << std::endl;
        Grainplot << "ASCII" << std::endl;
        Grainplot << "DATASET STRUCTURED_POINTS" << std::endl;
        Grainplot << "DIMENSIONS " << nx - 2 << " " << ny - 2 << " " << nz - 2 << std::endl;
        Grainplot << "ORIGIN 0 0 0" << std::endl;
        Grainplot << "SPACING 1 1 1" << std::endl;
        Grainplot << std::fixed << "POINT_DATA " << (nx - 2) * (ny - 2) * (nz - 2) << std::endl;
        Grainplot << "SCALARS Angle_z int 1" << std::endl;
        Grainplot << "LOOKUP_TABLE default" << std::endl;

        for (int k = 1; k < nz - 1; k++) {
            for (int j = 1; j < ny - 1; j++) {
                for (int i = 1; i < nx - 1; i++) {
                    if (CellType_WholeDomain[k][i][j] == Solid)
                        Grainplot << 0 << " ";
                    else if (CellType_WholeDomain[k][i][j] == Liquid)
                        Grainplot << 2 << " ";
                    else if (CellType_WholeDomain[k][i][j] == Active)
                        Grainplot << 1 << " ";
                    else
                        Grainplot << 3 << " ";
                }
            }
            Grainplot << std::endl;
        }
        Grainplot.close();
    }
    else {

        // Send non-ghost node data to rank 0
        int SBufSize = (MyXSlices - 2) * (MyYSlices - 2) * (nz - 2);
        int *SendBufGID = new int[SBufSize];
        int DataCounter = 0;
        for (int k = 1; k < nz - 1; k++) {
            for (int i = 1; i < MyXSlices - 1; i++) {
                for (int j = 1; j < MyYSlices - 1; j++) {
                    SendBufGID[DataCounter] = CellType(k * MyXSlices * MyYSlices + i * MyYSlices + j);
                    DataCounter++;
                }
            }
        }
        MPI_Send(SendBufGID, SBufSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}

//*****************************************************************************/
void PrintOrientations(std::string FName, int nx, int ny, int nz,
                       std::vector<std::vector<std::vector<bool>>> Melted_WholeDomain,
                       std::vector<std::vector<std::vector<int>>> GrainID_WholeDomain, int NGrainOrientations) {

    // Histogram of orientations for texture determination
    std::cout << "Printing file of grain orientation distribution, out of " << NGrainOrientations
              << " possible orientations for each grain" << std::endl;
    std::ofstream Grainplot0;
    std::string FName0 = FName + "_Orientations.csv";
    Grainplot0.open(FName0);
    ViewI_H GOHistogram("grain_orientations", NGrainOrientations);
    std::cout << "Histogram made" << std::endl;
    // frequency data on grain ids
    for (int k = 1; k < nz - 1; k++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                if (Melted_WholeDomain[k][i][j]) {
                    int GOVal = (abs(GrainID_WholeDomain[k][i][j]) - 1) % NGrainOrientations;
                    GOHistogram(GOVal)++;
                }
            }
        }
    }
    std::cout << "Histogram filled" << std::endl;
    for (int i = 0; i < NGrainOrientations; i++) {
        Grainplot0 << GOHistogram[i] << std::endl;
    }
    Grainplot0.close();
}

//*****************************************************************************/
void PrintGrainIDs(std::string FName, int nx, int ny, int nz,
                   std::vector<std::vector<std::vector<bool>>> Melted_WholeDomain,
                   std::vector<std::vector<std::vector<int>>> GrainID_WholeDomain) {

    std::string FName1 = FName + ".csv";
    std::cout << "Printing file of Grain ID values" << std::endl;
    std::ofstream Grainplot1;
    Grainplot1.open(FName1);
    Grainplot1 << "Outermost dimension is nz = " << nz - 2 << std::endl;
    Grainplot1 << "ny = " << ny - 2 << std::endl;
    Grainplot1 << "Innermost dimension is nx = " << nx - 2 << std::endl;
    for (int k = 1; k < nz - 1; k++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                if (!Melted_WholeDomain[k][i][j])
                    Grainplot1 << 0 << std::endl;
                else
                    Grainplot1 << GrainID_WholeDomain[k][i][j] << std::endl;
            }
        }
    }
    Grainplot1.close();
}

//*****************************************************************************/
void PrintParaview(std::string FName, int nx, int ny, int nz,
                   std::vector<std::vector<std::vector<bool>>> Melted_WholeDomain,
                   std::vector<std::vector<std::vector<int>>> GrainID_WholeDomain, ViewI_H GrainOrientation,
                   ViewF_H GrainUnitVector, int NGrainOrientations) {

    std::string FName2 = FName + ".vtk";
    std::cout << "Printing Paraview file of grain misorientations" << std::endl;
    // Print grain orientations to file
    std::ofstream Grainplot2;
    // New vtk file
    Grainplot2.open(FName2);
    Grainplot2 << "# vtk DataFile Version 3.0" << std::endl;
    Grainplot2 << "vtk output" << std::endl;
    Grainplot2 << "ASCII" << std::endl;
    Grainplot2 << "DATASET STRUCTURED_POINTS" << std::endl;
    Grainplot2 << "DIMENSIONS " << nx - 2 << " " << ny - 2 << " " << nz - 2 << std::endl;
    Grainplot2 << "ORIGIN 0 0 0" << std::endl;
    Grainplot2 << "SPACING 1 1 1" << std::endl;
    Grainplot2 << std::fixed << "POINT_DATA " << (nx - 2) * (ny - 2) * (nz - 2) << std::endl;
    Grainplot2 << "SCALARS Angle_z int 1" << std::endl;
    Grainplot2 << "LOOKUP_TABLE default" << std::endl;

    int NVC = 0;
    for (int k = 1; k < nz - 1; k++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                if (!Melted_WholeDomain[k][i][j])
                    Grainplot2 << 200 << " ";
                else {
                    int RoundedAngle;
                    int MyOrientation =
                        GrainOrientation(((abs(GrainID_WholeDomain[k][i][j]) - 1) % NGrainOrientations));
                    double AngleZmin = 62.7;
                    for (int ll = 0; ll < 3; ll++) {
                        double AngleZ = abs((180 / M_PI) * acos(GrainUnitVector(9 * MyOrientation + 3 * ll + 2)));
                        if (AngleZ < AngleZmin) {
                            AngleZmin = AngleZ;
                        }
                    }
                    if (GrainID_WholeDomain[k][i][j] < 0) {
                        RoundedAngle = round(AngleZmin) + 100;
                        NVC++;
                    }
                    else
                        RoundedAngle = round(AngleZmin);
                    Grainplot2 << RoundedAngle << " ";
                }
            }
        }
        Grainplot2 << std::endl;
    }
    Grainplot2.close();
    std::cout << "Volume fraction of domain claimed by nucleated grains: "
              << (float)(NVC) / (float)((nx - 2) * (ny - 2) * (nz - 2)) << std::endl;
}

//*****************************************************************************/
void PrintGrainIDsForExaConstit(std::string FName, int nx, int ny, int nz,
                                std::vector<std::vector<std::vector<int>>> GrainID_WholeDomain, double deltax) {

    std::string FName1 = FName + "_ExaConstit.csv";
    std::cout << "Printing file of X, Y, Z, GrainID values for ExaConstit" << std::endl;

    std::ofstream Grainplot1;
    Grainplot1.open(FName1);
    Grainplot1 << "Coordinates are in CA units, 1 cell = " << deltax
               << " microns. Data is cell-centered. Origin at 0,0,0, domain size is " << nx - 2 << " by " << ny - 2
               << " by " << nz - 2 << " cells" << std::endl;
    Grainplot1 << "X coord, Y coord, Z coord, Grain ID" << std::endl;
    for (int k = 1; k < nz - 1; k++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                Grainplot1 << i - 1 << "," << j - 1 << "," << k - 1 << "," << GrainID_WholeDomain[k][i][j] << std::endl;
            }
        }
    }
    Grainplot1.close();
}

//*****************************************************************************/
void PrintGrainAreas(std::string FName, double deltax, int nx, int ny, int nz,
                     std::vector<std::vector<std::vector<int>>> GrainID_WholeDomain) {

    std::string FName1 = FName + "_GrainAreas.csv";
    std::cout << "Printing file of grain area values (in square microns) for all Z coordinates" << std::endl;
    std::ofstream Grainplot1;
    Grainplot1.open(FName1);
    for (int k = 1; k <= nz - 2; k++) {
        std::vector<int> GIDVals_ThisLayer;
        for (int j = 1; j <= ny - 2; j++) {
            for (int i = 1; i <= nx - 2; i++) {
                GIDVals_ThisLayer.push_back(GrainID_WholeDomain[k][i][j]);
            }
        }
        std::vector<int>::iterator it;
        sort(GIDVals_ThisLayer.begin(), GIDVals_ThisLayer.end());
        it = std::unique(GIDVals_ThisLayer.begin(), GIDVals_ThisLayer.end());
        int CellsThisLayer = GIDVals_ThisLayer.size();
        GIDVals_ThisLayer.resize(std::distance(GIDVals_ThisLayer.begin(), it));
        int GrainsThisLayer = GIDVals_ThisLayer.size();
        double MeanGrainAreaThisLayer = (double)(CellsThisLayer) / (double)(GrainsThisLayer);
        Grainplot1 << MeanGrainAreaThisLayer * deltax * deltax / pow(10, -12) << std::endl;
    }
    Grainplot1.close();
}

//*****************************************************************************/
void PrintWeightedGrainAreas(std::string FName, double deltax, int nx, int ny, int nz,
                             std::vector<std::vector<std::vector<int>>> GrainID_WholeDomain) {

    std::string FName1 = FName + "_WeightedGrainAreas.csv";
    std::cout << "Printing file of Grain ID values (in square microns) for every 10th Z coordinate" << std::endl;
    std::ofstream Grainplot1;
    Grainplot1.open(FName1);
    for (int k = 1; k <= nz - 2; k++) {
        if ((k + 2) % 12 == 0) {
            std::vector<int> GIDVals_ThisLayer;
            for (int j = 1; j <= ny - 2; j++) {
                for (int i = 1; i <= nx - 2; i++) {
                    GIDVals_ThisLayer.push_back(GrainID_WholeDomain[k][i][j]);
                }
            }
            std::vector<int> GIDAllVals_ThisLayer;
            GIDAllVals_ThisLayer = GIDVals_ThisLayer;
            std::vector<int>::iterator it;
            sort(GIDVals_ThisLayer.begin(), GIDVals_ThisLayer.end());
            it = std::unique(GIDVals_ThisLayer.begin(), GIDVals_ThisLayer.end());
            int CellsThisLayer = GIDVals_ThisLayer.size();
            GIDVals_ThisLayer.resize(std::distance(GIDVals_ThisLayer.begin(), it));
            int GrainsThisLayer = GIDVals_ThisLayer.size();
            long int AreaXArea = 0;
            for (int l = 0; l < GrainsThisLayer; l++) {
                long int MyGrainArea = 0;
                for (int ll = 0; ll < CellsThisLayer; ll++) {
                    if (GIDVals_ThisLayer[l] == GIDAllVals_ThisLayer[ll])
                        MyGrainArea++;
                }
                AreaXArea += MyGrainArea * MyGrainArea;
            }
            double WeightedArea = ((double)(AreaXArea) / (double)((nx - 2) * (ny - 2)));
            Grainplot1 << WeightedArea * deltax * deltax / pow(10, -12) << std::endl;
        }
    }
    Grainplot1.close();
}
