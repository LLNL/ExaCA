#include "header.h"
using namespace std;

//*****************************************************************************/
// Initial placement of data in ghost nodes
void GhostNodesInit_GPU(int, int, int DecompositionStrategy, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int ZBound_Low, int nzActive, int LocalActiveDomainSize, int NGrainOrientations, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewF GrainUnitVector, ViewI GrainOrientation, ViewI GrainID, ViewI CellType, ViewF DOCenter, ViewF DiagonalLength, ViewF CritDiagonalLength, ViewI Locks) {

    // Fill buffers with ghost node data following initialization of data on GPUs
    // Similar to the calls to GhostNodes1D/GhostNodes2D, but the information sent/received in the halo regions is different
    // Need to send and receive cell type data, grain id data from other ranks
    Buffer3D BufferA("BufferA",MyXSlices,nzActive,2);
    Buffer3D BufferB("BufferB",MyXSlices,nzActive,2);
    Buffer3D BufferC("BufferC",MyYSlices,nzActive,5);
    Buffer3D BufferD("BufferD",MyYSlices,nzActive,2);
    Buffer2D BufferE("BufferE",nzActive,2);
    Buffer2D BufferF("BufferF",nzActive,2);
    Buffer2D BufferG("BufferG",nzActive,2);
    Buffer2D BufferH("BufferH",nzActive,2);
    Buffer3D BufferAR("BufferAR",MyXSlices,nzActive,2);
    Buffer3D BufferBR("BufferBR",MyXSlices,nzActive,2);
    Buffer3D BufferCR("BufferCR",MyYSlices,nzActive,2);
    Buffer3D BufferDR("BufferDR",MyYSlices,nzActive,2);
    Buffer2D BufferER("BufferER",nzActive,2);
    Buffer2D BufferFR("BufferFR",nzActive,2);
    Buffer2D BufferGR("BufferGR",nzActive,2);
    Buffer2D BufferHR("BufferHR",nzActive,2);

    // Load send buffers
    Kokkos::parallel_for ("GNInit",LocalActiveDomainSize, KOKKOS_LAMBDA (const long int& D3D1ConvPosition) {
        int RankZ = D3D1ConvPosition/(MyXSlices*MyYSlices);
        int Rem = D3D1ConvPosition % (MyXSlices*MyYSlices);
        int RankX = Rem/MyYSlices;
        int RankY = Rem % MyYSlices;
        int GlobalZ = RankZ + ZBound_Low;
        int D3D1ConvPositionGlobal = GlobalZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
        int GhostGID = GrainID(D3D1ConvPositionGlobal);
        int GhostCT = CellType(D3D1ConvPositionGlobal);
        // Collect data for the ghost nodes:
        if (RankY == 1) {
            BufferA(RankX,RankZ,0) = GhostGID;
            BufferA(RankX,RankZ,1) = GhostCT;
        }
        else if (RankY == MyYSlices-2) {
            BufferB(RankX,RankZ,0) = GhostGID;
            BufferB(RankX,RankZ,1) = GhostCT;
        }
        if (DecompositionStrategy > 1) {
            if (RankX == MyXSlices-2) {
                BufferC(RankY,RankZ,0) = GhostGID;
                BufferC(RankY,RankZ,1) = GhostCT;
            }
            else if (RankX == 1) {
                BufferD(RankY,RankZ,0) = GhostGID;
                BufferD(RankY,RankZ,1) = GhostCT;
            }
            if ((RankY == 1)&&(RankX == MyXSlices-2)) {
                BufferE(RankZ,0) = GhostGID;
                BufferE(RankZ,1) = GhostCT;
            }
            else if ((RankY == 1)&&(RankX == MyXSlices-2)) {
                BufferG(RankZ,0) = GhostGID;
                BufferG(RankZ,1) = GhostCT;
            }
            else if ((RankY == MyYSlices-2)&&(RankX == MyXSlices-2)) {
                BufferF(RankZ,0) = GhostGID;
                BufferF(RankZ,1) = GhostCT;
            }
            else if ((RankY == MyYSlices-2)&&(RankX == 1)) {
                BufferH(RankZ,0) = GhostGID;
                BufferH(RankZ,1) = GhostCT;
            }
        }
    });

    int NumberOfSends;
    if (DecompositionStrategy > 1) NumberOfSends = 8;
    else NumberOfSends = 2;
    std::vector<MPI_Request> SendRequests(NumberOfSends, MPI_REQUEST_NULL);
    std::vector<MPI_Request> RecvRequests(NumberOfSends, MPI_REQUEST_NULL);
    
    // Send data to each other rank (MPI_Isend)
    MPI_Isend(BufferA.data(),2*MyXSlices*nzActive,MPI_FLOAT,MyLeft,0,MPI_COMM_WORLD,&SendRequests[0]);
    MPI_Isend(BufferB.data(),2*MyXSlices*nzActive,MPI_FLOAT,MyRight,0,MPI_COMM_WORLD,&SendRequests[1]);
    if (DecompositionStrategy > 1) {
        MPI_Isend(BufferC.data(),2*MyYSlices*nzActive,MPI_FLOAT,MyIn,0,MPI_COMM_WORLD,&SendRequests[2]);
        MPI_Isend(BufferD.data(),2*MyYSlices*nzActive,MPI_FLOAT,MyOut,0,MPI_COMM_WORLD,&SendRequests[3]);
        MPI_Isend(BufferE.data(),2*nzActive,MPI_FLOAT,MyLeftIn,0,MPI_COMM_WORLD,&SendRequests[4]);
        MPI_Isend(BufferF.data(),2*nzActive,MPI_FLOAT,MyRightOut,0,MPI_COMM_WORLD,&SendRequests[5]);
        MPI_Isend(BufferG.data(),2*nzActive,MPI_FLOAT,MyRightIn,0,MPI_COMM_WORLD,&SendRequests[6]);
        MPI_Isend(BufferH.data(),2*nzActive,MPI_FLOAT,MyLeftOut,0,MPI_COMM_WORLD,&SendRequests[7]);
    }
    
    // Receive buffers for all neighbors (MPI_Irecv)
    MPI_Irecv(BufferAR.data(),2*MyXSlices*nzActive,MPI_FLOAT,MyLeft,0,MPI_COMM_WORLD,&RecvRequests[0]);
    MPI_Irecv(BufferBR.data(),2*MyXSlices*nzActive,MPI_FLOAT,MyRight,0,MPI_COMM_WORLD,&RecvRequests[1]);
    if (DecompositionStrategy > 1) {
        MPI_Irecv(BufferCR.data(),2*MyYSlices*nzActive,MPI_FLOAT,MyIn,0,MPI_COMM_WORLD,&RecvRequests[2]);
        MPI_Irecv(BufferDR.data(),2*MyYSlices*nzActive,MPI_FLOAT,MyOut,0,MPI_COMM_WORLD,&RecvRequests[3]);
        MPI_Irecv(BufferER.data(),2*nzActive,MPI_FLOAT,MyLeftIn,0,MPI_COMM_WORLD,&RecvRequests[4]);
        MPI_Irecv(BufferFR.data(),2*nzActive,MPI_FLOAT,MyRightOut,0,MPI_COMM_WORLD,&RecvRequests[5]);
        MPI_Irecv(BufferGR.data(),2*nzActive,MPI_FLOAT,MyRightIn,0,MPI_COMM_WORLD,&RecvRequests[6]);
        MPI_Irecv(BufferHR.data(),2*nzActive,MPI_FLOAT,MyLeftOut,0,MPI_COMM_WORLD,&RecvRequests[7]);
    }

    // unpack in any order
    bool unpack_complete = false;
    while (!unpack_complete) {
        // Get the next buffer to unpack from rank "unpack_index"
        int unpack_index = MPI_UNDEFINED;
        MPI_Waitany(NumberOfSends, RecvRequests.data(), &unpack_index, MPI_STATUS_IGNORE);
        
        // If there are no more buffers to unpack, leave the while loop
        if (MPI_UNDEFINED == unpack_index) {
            unpack_complete = true;
        }
        // Otherwise unpack the next buffer.
        else {
            Kokkos::parallel_for("BufferUnpack",nzActive, KOKKOS_LAMBDA (const int& RankZ) {
                
                int RankX, RankY;
                int GlobalZ = RankZ + ZBound_Low;

                // Which rank was the data received from?
                if ((unpack_index == 0)||(unpack_index == 1)) {
                    // Data receieved from MyLeft or MyRight
                    if (unpack_index == 0) RankY = 0;
                    else RankY = MyYSlices-1;
                    // If 1D decomposition, unpack all X coordinates from buffer
                    int XIterationRangeStart = 0;
                    int XIterationRangeEnd = MyXSlices;
                    // If 2D decomposition, unpack all X coordinates except the two on the ends from the buffer (those are filled using corner ghost nodes)
                    if (DecompositionStrategy > 1) {
                        XIterationRangeStart = 1;
                        XIterationRangeEnd = MyXSlices-1;
                    }
                    for (int RankX=XIterationRangeStart; RankX<XIterationRangeEnd; RankX++) {
                        int GlobalCellLocation = GlobalZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                        if (unpack_index == 0) {
                            GrainID(GlobalCellLocation) = BufferAR(RankX,RankZ,0);
                            CellType(GlobalCellLocation) = BufferAR(RankX,RankZ,1);
                        }
                        else {
                            GrainID(GlobalCellLocation) = BufferBR(RankX,RankZ,0);
                            CellType(GlobalCellLocation) = BufferBR(RankX,RankZ,1);
                        }
                        if (CellType(GlobalCellLocation) == Active) {
                            // Mark cell for later placement of octahedra, critical diagonal length calculations
                            CellType(GlobalCellLocation) = Ghost1;
                        }
                    }
                }

                else if ((unpack_index == 2)||(unpack_index == 3)) {
                    // Data received from MyIn or MyOut
                    if (unpack_index == 2) RankX = MyXSlices-1;
                    else RankX = 0;
                    // This is a 2D decomposition, only unpack Y coordinates other than the two on the ends from the buffer (those are filled using corner ghost nodes)
                    for (int RankY=1; RankY<MyYSlices-1; RankY++) {
                        int GlobalCellLocation = GlobalZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                        if (unpack_index == 2) {
                            GrainID(GlobalCellLocation) = BufferCR(RankY,RankZ,0);
                            CellType(GlobalCellLocation) = BufferCR(RankY,RankZ,1);
                        }
                        else {
                            GrainID(GlobalCellLocation) = BufferDR(RankY,RankZ,0);
                            CellType(GlobalCellLocation) = BufferDR(RankY,RankZ,1);
                        }
                        if (CellType(GlobalCellLocation) == Active) {
                            // Mark cell for later placement of octahedra, critical diagonal length calculations
                            CellType(GlobalCellLocation) = Ghost1;
                        }
                    }
                }
                else {
                    if (unpack_index == 4) {
                        // Data received from MyLeftIn
                        RankX = MyXSlices-1;
                        RankY = 0;
                    }
                    else if (unpack_index == 5) {
                        // Data received from MyRightOut
                        RankX = MyXSlices-1;
                        RankY = MyYSlices-1;
                    }
                    else if (unpack_index == 6) {
                        // Data received from MyRightIn
                        RankX = 0;
                        RankY = MyYSlices-1;
                    }
                    else {
                        // Data received from MyLeftOut
                        RankX = 0;
                        RankY = 0;
                    }
                    int GlobalCellLocation = GlobalZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                    if (unpack_index == 4) {
                        GrainID(GlobalCellLocation) = BufferER(RankZ,0);
                        CellType(GlobalCellLocation) = BufferER(RankZ,1);
                    }
                    else if (unpack_index == 5) {
                        GrainID(GlobalCellLocation) = BufferFR(RankZ,0);
                        CellType(GlobalCellLocation) = BufferFR(RankZ,1);
                    }
                    else if (unpack_index == 6) {
                        GrainID(GlobalCellLocation) = BufferGR(RankZ,0);
                        CellType(GlobalCellLocation) = BufferGR(RankZ,1);
                    }
                    else {
                        GrainID(GlobalCellLocation) = BufferHR(RankZ,0);
                        CellType(GlobalCellLocation) = BufferHR(RankZ,1);
                    }
                    if (CellType(GlobalCellLocation) == Active) {
                        // Mark cell for later placement of octahedra, critical diagonal length calculations
                        CellType(GlobalCellLocation) = Ghost1;
                    }
                }
            });
        }
    }
    // Wait on send requests
    MPI_Waitall(NumberOfSends, SendRequests.data(), MPI_STATUSES_IGNORE );

    // Octahedra for active cells that were marked in the previous loop, initially have centers aligned with cell centers, Diagonal lengths of 0.01
    // Also ensure locks values match up with cell types that may have been changed through ghost node update (Wall/Solid/Active Locks = 0, Liquid/LiqSol Locks = 1)
    Kokkos::parallel_for ("GNInit",LocalActiveDomainSize, KOKKOS_LAMBDA (const int& CellLocation) {
        int RankZ = CellLocation/(MyXSlices*MyYSlices);
        int Rem = CellLocation % (MyXSlices*MyYSlices);
        int RankX = Rem/MyYSlices;
        int RankY = Rem % MyYSlices;
        int GlobalZ = RankZ + ZBound_Low;
        int GlobalCellLocation = GlobalZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
        if ((CellType(GlobalCellLocation) == Wall)||(CellType(GlobalCellLocation) == Solid)||(CellType(GlobalCellLocation) == Active)) {
            Locks(CellLocation) = 0;
        }
        else if ((CellType(GlobalCellLocation) == Liquid)||(CellType(GlobalCellLocation) == LiqSol)) {
            Locks(CellLocation) = 1;
        }
        if (CellType(GlobalCellLocation) == Ghost1) {
            CellType(GlobalCellLocation) = Active;
            Locks(CellLocation) = 0;
            DOCenter((long int)(3)*CellLocation) = RankX + MyXOffset + 0.5;
            DOCenter((long int)(3)*CellLocation+(long int)(1)) = RankY + MyYOffset + 0.5;
            DOCenter((long int)(3)*CellLocation+(long int)(2)) = GlobalZ + 0.5;
            int MyOrientation = GrainOrientation(((abs(GrainID(GlobalCellLocation)) - 1) % NGrainOrientations));
            DiagonalLength(CellLocation) = 0.01;
            // Global coordinates of cell center
            double xp = RankX + MyXOffset + 0.5;
            double yp = RankY + MyYOffset + 0.5;
            double zp = GlobalZ + 0.5;
            // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
            for (int n=0; n<26; n++)  {
                
                int MyNeighborX = RankX + NeighborX(n);
                int MyNeighborY = RankY + NeighborY(n);
                int MyNeighborZ = RankZ + NeighborZ(n);
                long int NeighborPosition = MyNeighborZ*MyXSlices*MyYSlices + MyNeighborX*MyYSlices + MyNeighborY;

                if (NeighborPosition == CellLocation) {
                    // Do not calculate critical diagonal length req'd for the newly captured cell to capture the original
                    CritDiagonalLength((long int)(26)*NeighborPosition+(long int)(n)) = 10000000.0;
                }
                
                // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the global coordinates of the center of a neighbor cell
                double x0 = xp + NeighborX(n) - (RankX + MyXOffset + 0.5);
                double y0 = yp + NeighborY(n) - (RankY + MyYOffset + 0.5);
                double z0 = zp + NeighborZ(n) - (GlobalZ + 0.5);
                // mag0 is the magnitude of (x0,y0,z0)
                double mag0 = pow(pow(x0,2.0) + pow(y0,2.0) + pow(z0,2.0),0.5);
                
                // Calculate unit vectors for the octahedron that intersect the new cell center
                double Diag1X, Diag1Y, Diag1Z, Diag2X, Diag2Y, Diag2Z, Diag3X, Diag3Y, Diag3Z;
                double Angle1 = (GrainUnitVector(9*MyOrientation)*x0 + GrainUnitVector(9*MyOrientation + 1)*y0 + GrainUnitVector(9*MyOrientation + 2)*z0)/mag0;
                if (Angle1 < 0) {
                    Diag1X = GrainUnitVector(9*MyOrientation);
                    Diag1Y = GrainUnitVector(9*MyOrientation + 1);
                    Diag1Z = GrainUnitVector(9*MyOrientation + 2);
                }
                else {
                    Diag1X = -GrainUnitVector(9*MyOrientation);
                    Diag1Y = -GrainUnitVector(9*MyOrientation + 1);
                    Diag1Z = -GrainUnitVector(9*MyOrientation + 2);
                }

                double Angle2 = (GrainUnitVector(9*MyOrientation + 3)*x0 + GrainUnitVector(9*MyOrientation + 4)*y0 + GrainUnitVector(9*MyOrientation + 5)*z0)/mag0;
                if (Angle2 < 0) {
                    Diag2X = GrainUnitVector(9*MyOrientation + 3);
                    Diag2Y = GrainUnitVector(9*MyOrientation + 4);
                    Diag2Z = GrainUnitVector(9*MyOrientation + 5);
                }
                else {
                    Diag2X = -GrainUnitVector(9*MyOrientation + 3);
                    Diag2Y = -GrainUnitVector(9*MyOrientation + 4);
                    Diag2Z = -GrainUnitVector(9*MyOrientation + 5);
                }

                double Angle3 = (GrainUnitVector(9*MyOrientation + 6)*x0 + GrainUnitVector(9*MyOrientation + 7)*y0 + GrainUnitVector(9*MyOrientation + 8)*z0)/mag0;
                if (Angle3 < 0) {
                    Diag3X = GrainUnitVector(9*MyOrientation + 6);
                    Diag3Y = GrainUnitVector(9*MyOrientation + 7);
                    Diag3Z = GrainUnitVector(9*MyOrientation + 8);
                }
                else {
                    Diag3X = -GrainUnitVector(9*MyOrientation + 6);
                    Diag3Y = -GrainUnitVector(9*MyOrientation + 7);
                    Diag3Z = -GrainUnitVector(9*MyOrientation + 8);
                }
                
                double U1[3], U2[3], UU[3], Norm[3];
                U1[0] = Diag2X - Diag1X;
                U1[1] = Diag2Y - Diag1Y;
                U1[2] = Diag2Z - Diag1Z;
                U2[0] = Diag3X - Diag1X;
                U2[1] = Diag3Y - Diag1Y;
                U2[2] = Diag3Z - Diag1Z;
                UU[0] = U1[1]*U2[2] - U1[2]*U2[1];
                UU[1] = U1[2]*U2[0] - U1[0]*U2[2];
                UU[2] = U1[0]*U2[1] - U1[1]*U2[0];
                double NDem = sqrt(UU[0]*UU[0] + UU[1]*UU[1] + UU[2]*UU[2]);
                Norm[0] = UU[0]/NDem;
                Norm[1] = UU[1]/NDem;
                Norm[2] = UU[2]/NDem;
                // normal to capturing plane
                double normx = Norm[0];
                double normy = Norm[1];
                double normz = Norm[2];
                double ParaT = (normx*x0+normy*y0+normz*z0)/(normx*Diag1X+normy*Diag1Y+normz*Diag1Z);
                float CDLVal = pow(pow(ParaT*Diag1X,2.0) + pow(ParaT*Diag1Y,2.0) + pow(ParaT*Diag1Z,2.0),0.5);
                CritDiagonalLength((long int)(26)*CellLocation+(long int)(n)) = CDLVal;
            }
        }
    });

}

//*****************************************************************************/
// 2D domain decomposition: update ghost nodes with new cell data from Nucleation and CellCapture routines
void GhostNodes2D_GPU(int, int, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewI CellType, ViewF DOCenter, ViewI GrainID, ViewF GrainUnitVector, ViewI GrainOrientation, ViewF DiagonalLength, ViewF CritDiagonalLength, int NGrainOrientations, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferC, Buffer2D BufferD, Buffer2D BufferE, Buffer2D BufferF, Buffer2D BufferG, Buffer2D BufferH, Buffer2D BufferAR, Buffer2D BufferBR, Buffer2D BufferCR, Buffer2D BufferDR, Buffer2D BufferER, Buffer2D BufferFR, Buffer2D BufferGR, Buffer2D BufferHR, int BufSizeX, int BufSizeY, int BufSizeZ, ViewI Locks, int ZBound_Low) {

    Kokkos::fence();
    MPI_Barrier(MPI_COMM_WORLD);

    // Allocate requests.
    std::vector<MPI_Request> SendRequests(8, MPI_REQUEST_NULL);
    std::vector<MPI_Request> RecvRequests(8, MPI_REQUEST_NULL);
    
    // Send data to each other rank (MPI_Isend)
    MPI_Isend(BufferA.data(),5*BufSizeX*BufSizeZ,MPI_FLOAT,MyLeft,0,MPI_COMM_WORLD,&SendRequests[0]);
    MPI_Isend(BufferB.data(),5*BufSizeX*BufSizeZ,MPI_FLOAT,MyRight,0,MPI_COMM_WORLD,&SendRequests[1]);
    MPI_Isend(BufferC.data(),5*BufSizeY*BufSizeZ,MPI_FLOAT,MyIn,0,MPI_COMM_WORLD,&SendRequests[2]);
    MPI_Isend(BufferD.data(),5*BufSizeY*BufSizeZ,MPI_FLOAT,MyOut,0,MPI_COMM_WORLD,&SendRequests[3]);
    MPI_Isend(BufferE.data(),5*BufSizeZ,MPI_FLOAT,MyLeftIn,0,MPI_COMM_WORLD,&SendRequests[4]);
    MPI_Isend(BufferF.data(),5*BufSizeZ,MPI_FLOAT,MyRightOut,0,MPI_COMM_WORLD,&SendRequests[5]);
    MPI_Isend(BufferG.data(),5*BufSizeZ,MPI_FLOAT,MyRightIn,0,MPI_COMM_WORLD,&SendRequests[6]);
    MPI_Isend(BufferH.data(),5*BufSizeZ,MPI_FLOAT,MyLeftOut,0,MPI_COMM_WORLD,&SendRequests[7]);
    
    // Receive buffers for all neighbors (MPI_Irecv)
    MPI_Irecv(BufferAR.data(),5*BufSizeX*BufSizeZ,MPI_FLOAT,MyLeft,0,MPI_COMM_WORLD,&RecvRequests[0]);
    MPI_Irecv(BufferBR.data(),5*BufSizeX*BufSizeZ,MPI_FLOAT,MyRight,0,MPI_COMM_WORLD,&RecvRequests[1]);
    MPI_Irecv(BufferCR.data(),5*BufSizeY*BufSizeZ,MPI_FLOAT,MyIn,0,MPI_COMM_WORLD,&RecvRequests[2]);
    MPI_Irecv(BufferDR.data(),5*BufSizeY*BufSizeZ,MPI_FLOAT,MyOut,0,MPI_COMM_WORLD,&RecvRequests[3]);
    MPI_Irecv(BufferER.data(),5*BufSizeZ,MPI_FLOAT,MyLeftIn,0,MPI_COMM_WORLD,&RecvRequests[4]);
    MPI_Irecv(BufferFR.data(),5*BufSizeZ,MPI_FLOAT,MyRightOut,0,MPI_COMM_WORLD,&RecvRequests[5]);
    MPI_Irecv(BufferGR.data(),5*BufSizeZ,MPI_FLOAT,MyRightIn,0,MPI_COMM_WORLD,&RecvRequests[6]);
    MPI_Irecv(BufferHR.data(),5*BufSizeZ,MPI_FLOAT,MyLeftOut,0,MPI_COMM_WORLD,&RecvRequests[7]);
    
    // unpack in any order
    bool unpack_complete = false;
    while (!unpack_complete) {
        // Get the next buffer to unpack from rank "unpack_index"
        int unpack_index = MPI_UNDEFINED;
        MPI_Waitany(8, RecvRequests.data(), &unpack_index, MPI_STATUS_IGNORE);
        
        // If there are no more buffers to unpack, leave the while loop
        if (MPI_UNDEFINED == unpack_index) {
            unpack_complete = true;
        }
        // Otherwise unpack the next buffer.
        else {
            int RecvBufSize;
            if (unpack_index <= 1) {
                RecvBufSize = BufSizeX*BufSizeZ;
            }
            else if (unpack_index >= 4) {
                RecvBufSize = BufSizeZ;
            }
            else {
                RecvBufSize = BufSizeY*BufSizeZ;
            }
            Kokkos::parallel_for("BufferUnpack",RecvBufSize, KOKKOS_LAMBDA (const int& BufPosition) {
                int RankX, RankY, RankZ, NewGrainID;
                long int CellLocation;
                bool Place = false;
                float DOCenterX, DOCenterY, DOCenterZ, NewDiagonalLength;
                // Which rank was the data received from?
                if (unpack_index == 0) {
                    // Data receieved from MyLeft
                    // Adjust X Position by +1 since X = 0 is not included in this buffer
                    RankX = BufPosition % BufSizeX + 1;
                    RankY = 0;
                    RankZ = BufPosition/BufSizeX;
                    CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
                    if ((BufferAR(BufPosition,4) > 0)&&(DiagonalLength(CellLocation) == 0)) {
                        Place = true;
                        NewGrainID = (int)(BufferAR(BufPosition,0));
                        DOCenterX = BufferAR(BufPosition,1);
                        DOCenterY = BufferAR(BufPosition,2);
                        DOCenterZ = BufferAR(BufPosition,3);
                        NewDiagonalLength = BufferAR(BufPosition,4);
                    }
                }
                else if (unpack_index == 1) {
                    // Data received from MyRight
                    // Adjust X Position by +1 since X = 0 is not included in this buffer
                    RankX = BufPosition % BufSizeX + 1;
                    RankY = MyYSlices-1;
                    RankZ = BufPosition/BufSizeX;
                    CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
                    if ((BufferBR(BufPosition,4) > 0)&&(DiagonalLength(CellLocation) == 0)) {
                        Place = true;
                        NewGrainID = (int)(BufferBR(BufPosition,0));
                        DOCenterX = BufferBR(BufPosition,1);
                        DOCenterY = BufferBR(BufPosition,2);
                        DOCenterZ = BufferBR(BufPosition,3);
                        NewDiagonalLength = BufferBR(BufPosition,4);
                    }
                }
                else if (unpack_index == 2) {
                    // Data received from MyIn
                    // Adjust Y Position by +1 since Y = 0 is not included in this buffer
                    RankX = MyXSlices-1;
                    RankY = BufPosition % BufSizeY + 1;
                    RankZ = BufPosition/BufSizeY;
                    CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
                    if ((BufferCR(BufPosition,4) > 0)&&(DiagonalLength(CellLocation) == 0)) {
                        Place = true;
                        NewGrainID = (int)(BufferCR(BufPosition,0));
                        DOCenterX = BufferCR(BufPosition,1);
                        DOCenterY = BufferCR(BufPosition,2);
                        DOCenterZ = BufferCR(BufPosition,3);
                        NewDiagonalLength = BufferCR(BufPosition,4);
                    }
                }
                else if (unpack_index == 3) {
                    // Data received from MyOut
                    // Adjust Y Position by +1 since Y = 0 is not included in this buffer
                    RankX = 0;
                    RankY = BufPosition % BufSizeY + 1;
                    RankZ = BufPosition/BufSizeY;
                    CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
                    if ((BufferDR(BufPosition,4) > 0)&&(DiagonalLength(CellLocation) == 0)) {
                        Place = true;
                        NewGrainID = (int)(BufferDR(BufPosition,0));
                        DOCenterX = BufferDR(BufPosition,1);
                        DOCenterY = BufferDR(BufPosition,2);
                        DOCenterZ = BufferDR(BufPosition,3);
                        NewDiagonalLength = BufferDR(BufPosition,4);
                    }
                }
                else if (unpack_index == 4) {
                    // Data received from MyLeftIn
                    RankX = MyXSlices-1;
                    RankY = 0;
                    RankZ = BufPosition;
                    CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
                    if ((BufferER(BufPosition,4) > 0)&&(DiagonalLength(CellLocation) == 0)) {
                        Place = true;
                        NewGrainID = (int)(BufferER(BufPosition,0));
                        DOCenterX = BufferER(BufPosition,1);
                        DOCenterY = BufferER(BufPosition,2);
                        DOCenterZ = BufferER(BufPosition,3);
                        NewDiagonalLength = BufferER(BufPosition,4);
                    }
                }
                else if (unpack_index == 5) {
                    // Data received from MyRightOut
                    RankX = MyXSlices-1;
                    RankY = MyYSlices-1;
                    RankZ = BufPosition;
                    CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
                    if ((BufferFR(BufPosition,4) > 0)&&(DiagonalLength(CellLocation) == 0)) {
                        Place = true;
                        NewGrainID = (int)(BufferFR(BufPosition,0));
                        DOCenterX = BufferFR(BufPosition,1);
                        DOCenterY = BufferFR(BufPosition,2);
                        DOCenterZ = BufferFR(BufPosition,3);
                        NewDiagonalLength = BufferFR(BufPosition,4);
                    }
                }
                else if (unpack_index == 6) {
                    // Data received from MyRightIn
                    RankX = 0;
                    RankY = MyYSlices-1;
                    RankZ = BufPosition;
                    CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
                    if ((BufferGR(BufPosition,4) > 0)&&(DiagonalLength(CellLocation) == 0)) {
                        Place = true;
                        NewGrainID = (int)(BufferGR(BufPosition,0));
                        DOCenterX = BufferGR(BufPosition,1);
                        DOCenterY = BufferGR(BufPosition,2);
                        DOCenterZ = BufferGR(BufPosition,3);
                        NewDiagonalLength = BufferGR(BufPosition,4);
                    }
                }
                else if (unpack_index == 7) {
                    // Data received from MyLeftOut
                    RankX = 0;
                    RankY = 0;
                    RankZ = BufPosition;
                    CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
                    if ((BufferHR(BufPosition,4) > 0)&&(DiagonalLength(CellLocation) == 0)) {
                        Place = true;
                        NewGrainID = (int)(BufferHR(BufPosition,0));
                        DOCenterX = BufferHR(BufPosition,1);
                        DOCenterY = BufferHR(BufPosition,2);
                        DOCenterZ = BufferHR(BufPosition,3);
                        NewDiagonalLength = BufferHR(BufPosition,4);
                    }
                }

                if (Place) {
                    // Update this ghost node cell's information with data from other rank
                    int GlobalZ = RankZ + ZBound_Low;
                    int GlobalCellLocation = GlobalZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                    CellType(GlobalCellLocation) = Active;
                    Locks(CellLocation) = 0;
                    GrainID(GlobalCellLocation) = NewGrainID;
                    DOCenter((long int)(3)*CellLocation) = DOCenterX;
                    DOCenter((long int)(3)*CellLocation+(long int)(1)) = DOCenterY;
                    DOCenter((long int)(3)*CellLocation+(long int)(2)) = DOCenterZ;
                    int MyOrientation = GrainOrientation(((abs(GrainID(GlobalCellLocation)) - 1) % NGrainOrientations));
                    DiagonalLength(CellLocation) = NewDiagonalLength;
                    // Global coordinates of cell center
                    double xp = RankX + MyXOffset + 0.5;
                    double yp = RankY + MyYOffset + 0.5;
                    double zp = GlobalZ + 0.5;
                    // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
                    for (int n=0; n<26; n++)  {
                        
                        int MyNeighborX = RankX + NeighborX(n);
                        int MyNeighborY = RankY + NeighborY(n);
                        int MyNeighborZ = RankZ + NeighborZ(n);
                        long int NeighborPosition = MyNeighborZ*MyXSlices*MyYSlices + MyNeighborX*MyYSlices + MyNeighborY;

                        if (NeighborPosition == CellLocation) {
                            // Do not calculate critical diagonal length req'd for the newly captured cell to capture the original
                            CritDiagonalLength((long int)(26)*NeighborPosition+(long int)(n)) = 10000000.0;
                        }
                        
                        // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the global coordinates of the center of a neighbor cell
                        double x0 = xp + NeighborX(n) - DOCenterX;
                        double y0 = yp + NeighborY(n) - DOCenterY;
                        double z0 = zp + NeighborZ(n) - DOCenterZ;
                        // mag0 is the magnitude of (x0,y0,z0)
                        double mag0 = pow(pow(x0,2.0) + pow(y0,2.0) + pow(z0,2.0),0.5);
                        
                        // Calculate unit vectors for the octahedron that intersect the new cell center
                        double Diag1X, Diag1Y, Diag1Z, Diag2X, Diag2Y, Diag2Z, Diag3X, Diag3Y, Diag3Z;
                        double Angle1 = (GrainUnitVector(9*MyOrientation)*x0 + GrainUnitVector(9*MyOrientation + 1)*y0 + GrainUnitVector(9*MyOrientation + 2)*z0)/mag0;
                        if (Angle1 < 0) {
                            Diag1X = GrainUnitVector(9*MyOrientation);
                            Diag1Y = GrainUnitVector(9*MyOrientation + 1);
                            Diag1Z = GrainUnitVector(9*MyOrientation + 2);
                        }
                        else {
                            Diag1X = -GrainUnitVector(9*MyOrientation);
                            Diag1Y = -GrainUnitVector(9*MyOrientation + 1);
                            Diag1Z = -GrainUnitVector(9*MyOrientation + 2);
                        }

                        double Angle2 = (GrainUnitVector(9*MyOrientation + 3)*x0 + GrainUnitVector(9*MyOrientation + 4)*y0 + GrainUnitVector(9*MyOrientation + 5)*z0)/mag0;
                        if (Angle2 < 0) {
                            Diag2X = GrainUnitVector(9*MyOrientation + 3);
                            Diag2Y = GrainUnitVector(9*MyOrientation + 4);
                            Diag2Z = GrainUnitVector(9*MyOrientation + 5);
                        }
                        else {
                            Diag2X = -GrainUnitVector(9*MyOrientation + 3);
                            Diag2Y = -GrainUnitVector(9*MyOrientation + 4);
                            Diag2Z = -GrainUnitVector(9*MyOrientation + 5);
                        }

                        double Angle3 = (GrainUnitVector(9*MyOrientation + 6)*x0 + GrainUnitVector(9*MyOrientation + 7)*y0 + GrainUnitVector(9*MyOrientation + 8)*z0)/mag0;
                        if (Angle3 < 0) {
                            Diag3X = GrainUnitVector(9*MyOrientation + 6);
                            Diag3Y = GrainUnitVector(9*MyOrientation + 7);
                            Diag3Z = GrainUnitVector(9*MyOrientation + 8);
                        }
                        else {
                            Diag3X = -GrainUnitVector(9*MyOrientation + 6);
                            Diag3Y = -GrainUnitVector(9*MyOrientation + 7);
                            Diag3Z = -GrainUnitVector(9*MyOrientation + 8);
                        }
                        
                        double U1[3], U2[3], UU[3], Norm[3];
                        U1[0] = Diag2X - Diag1X;
                        U1[1] = Diag2Y - Diag1Y;
                        U1[2] = Diag2Z - Diag1Z;
                        U2[0] = Diag3X - Diag1X;
                        U2[1] = Diag3Y - Diag1Y;
                        U2[2] = Diag3Z - Diag1Z;
                        UU[0] = U1[1]*U2[2] - U1[2]*U2[1];
                        UU[1] = U1[2]*U2[0] - U1[0]*U2[2];
                        UU[2] = U1[0]*U2[1] - U1[1]*U2[0];
                        double NDem = sqrt(UU[0]*UU[0] + UU[1]*UU[1] + UU[2]*UU[2]);
                        Norm[0] = UU[0]/NDem;
                        Norm[1] = UU[1]/NDem;
                        Norm[2] = UU[2]/NDem;
                        // normal to capturing plane
                        double normx = Norm[0];
                        double normy = Norm[1];
                        double normz = Norm[2];
                        double ParaT = (normx*x0+normy*y0+normz*z0)/(normx*Diag1X+normy*Diag1Y+normz*Diag1Z);
                        float CDLVal = pow(pow(ParaT*Diag1X,2.0) + pow(ParaT*Diag1Y,2.0) + pow(ParaT*Diag1Z,2.0),0.5);
                        CritDiagonalLength((long int)(26)*CellLocation+(long int)(n)) = CDLVal;
                    }
                }
            });
        }
    }
    // Wait on send requests
    MPI_Waitall(8, SendRequests.data(), MPI_STATUSES_IGNORE );
}

//*****************************************************************************/
// 1D domain decomposition: update ghost nodes with new cell data from Nucleation and CellCapture routines
void GhostNodes1D_GPU(int, int, int MyLeft, int MyRight, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewI CellType, ViewF DOCenter, ViewI GrainID, ViewF GrainUnitVector, ViewI GrainOrientation, ViewF DiagonalLength, ViewF CritDiagonalLength, int NGrainOrientations, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferAR, Buffer2D BufferBR, int BufSizeX, int, int BufSizeZ, ViewI Locks, int ZBound_Low) {
    
    std::vector<MPI_Request> SendRequests(2, MPI_REQUEST_NULL);
    std::vector<MPI_Request> RecvRequests(2, MPI_REQUEST_NULL);
    
    // Send data to each other rank (MPI_Isend)
    MPI_Isend(BufferA.data(),5*BufSizeX*BufSizeZ,MPI_FLOAT,MyLeft,0,MPI_COMM_WORLD,&SendRequests[0]);
    MPI_Isend(BufferB.data(),5*BufSizeX*BufSizeZ,MPI_FLOAT,MyRight,0,MPI_COMM_WORLD,&SendRequests[1]);
    
    // Receive buffers for all neighbors (MPI_Irecv)
    MPI_Irecv(BufferAR.data(),5*BufSizeX*BufSizeZ,MPI_FLOAT,MyLeft,0,MPI_COMM_WORLD,&RecvRequests[0]);
    MPI_Irecv(BufferBR.data(),5*BufSizeX*BufSizeZ,MPI_FLOAT,MyRight,0,MPI_COMM_WORLD,&RecvRequests[1]);
    
    // unpack in any order
    bool unpack_complete = false;
    while (!unpack_complete) {
        // Get the next buffer to unpack from rank "unpack_index"
        int unpack_index = MPI_UNDEFINED;
        MPI_Waitany(2, RecvRequests.data(), &unpack_index, MPI_STATUS_IGNORE);
        // If there are no more buffers to unpack, leave the while loop
        if (MPI_UNDEFINED == unpack_index) {
            unpack_complete = true;
        }
        // Otherwise unpack the next buffer.
        else {
            int RecvBufSize = BufSizeX*BufSizeZ;
            Kokkos::parallel_for("BufferUnpack",RecvBufSize, KOKKOS_LAMBDA (const int& BufPosition) {
                int RankX, RankY, RankZ, NewGrainID;
                long int CellLocation;
                float DOCenterX, DOCenterY, DOCenterZ, NewDiagonalLength;
                bool Place = false;
                RankZ = BufPosition/BufSizeX;
                RankX = BufPosition % BufSizeX;
                // Which rank was the data received from?
                if (unpack_index == 0) {
                    // Data receieved from MyLeft
                    RankY = 0;
                    CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
                    if ((BufferAR(BufPosition,4) > 0)&&(DiagonalLength(CellLocation) == 0)) {
                        Place = true;
                        NewGrainID = (int)(BufferAR(BufPosition,0));
                        DOCenterX = BufferAR(BufPosition,1);
                        DOCenterY = BufferAR(BufPosition,2);
                        DOCenterZ = BufferAR(BufPosition,3);
                        NewDiagonalLength = BufferAR(BufPosition,4);
                    }
                }
                else if (unpack_index == 1) {
                    // Data received from MyRight
                    RankY = MyYSlices-1;
                    CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
                    if ((BufferBR(BufPosition,4) > 0)&&(DiagonalLength(CellLocation) == 0)) {
                        Place = true;
                        NewGrainID = (int)(BufferBR(BufPosition,0));
                        DOCenterX = BufferBR(BufPosition,1);
                        DOCenterY = BufferBR(BufPosition,2);
                        DOCenterZ = BufferBR(BufPosition,3);
                        NewDiagonalLength = BufferBR(BufPosition,4);
                    }
                }
                if (Place) {
                    int GlobalZ = RankZ + ZBound_Low;
                    int GlobalCellLocation = GlobalZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                    
                    // Update this ghost node cell's information with data from other rank
                    GrainID(GlobalCellLocation) = NewGrainID;
                    Locks(CellLocation) = 0;
                    DOCenter((long int)(3)*CellLocation) = DOCenterX;
                    DOCenter((long int)(3)*CellLocation+(long int)(1)) = DOCenterY;
                    DOCenter((long int)(3)*CellLocation+(long int)(2)) = DOCenterZ;
                    int MyOrientation = GrainOrientation(((abs(GrainID(GlobalCellLocation)) - 1) % NGrainOrientations));
                    DiagonalLength(CellLocation) = NewDiagonalLength;
                    // Global coordinates of cell center
                    double xp = RankX + MyXOffset + 0.5;
                    double yp = RankY + MyYOffset + 0.5;
                    double zp = GlobalZ + 0.5;
                    // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
                    for (int n=0; n<26; n++)  {
                        
                        int MyNeighborX = RankX + NeighborX(n);
                        int MyNeighborY = RankY + NeighborY(n);
                        int MyNeighborZ = RankZ + NeighborZ(n);
                        long int NeighborPosition = MyNeighborZ*MyXSlices*MyYSlices + MyNeighborX*MyYSlices + MyNeighborY;
                        if (NeighborPosition == CellLocation) {
                            // Do not calculate critical diagonal length req'd for the newly captured cell to capture the original
                            CritDiagonalLength((long int)(26)*NeighborPosition+(long int)(n)) = 10000000.0;
                        }
                        
                        // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the global coordinates of the center of a neighbor cell
                        double x0 = xp + NeighborX(n) - DOCenterX;
                        double y0 = yp + NeighborY(n) - DOCenterY;
                        double z0 = zp + NeighborZ(n) - DOCenterZ;
                        // mag0 is the magnitude of (x0,y0,z0)
                        double mag0 = pow(pow(x0,2.0) + pow(y0,2.0) + pow(z0,2.0),0.5);
                        
                        // Calculate unit vectors for the octahedron that intersect the new cell center
                        double Diag1X, Diag1Y, Diag1Z, Diag2X, Diag2Y, Diag2Z, Diag3X, Diag3Y, Diag3Z;
                        double Angle1 = (GrainUnitVector(9*MyOrientation)*x0 + GrainUnitVector(9*MyOrientation + 1)*y0 + GrainUnitVector(9*MyOrientation + 2)*z0)/mag0;
                        if (Angle1 < 0) {
                            Diag1X = GrainUnitVector(9*MyOrientation);
                            Diag1Y = GrainUnitVector(9*MyOrientation + 1);
                            Diag1Z = GrainUnitVector(9*MyOrientation + 2);
                        }
                        else {
                            Diag1X = -GrainUnitVector(9*MyOrientation);
                            Diag1Y = -GrainUnitVector(9*MyOrientation + 1);
                            Diag1Z = -GrainUnitVector(9*MyOrientation + 2);
                        }

                        double Angle2 = (GrainUnitVector(9*MyOrientation + 3)*x0 + GrainUnitVector(9*MyOrientation + 4)*y0 + GrainUnitVector(9*MyOrientation + 5)*z0)/mag0;
                        if (Angle2 < 0) {
                            Diag2X = GrainUnitVector(9*MyOrientation + 3);
                            Diag2Y = GrainUnitVector(9*MyOrientation + 4);
                            Diag2Z = GrainUnitVector(9*MyOrientation + 5);
                        }
                        else {
                            Diag2X = -GrainUnitVector(9*MyOrientation + 3);
                            Diag2Y = -GrainUnitVector(9*MyOrientation + 4);
                            Diag2Z = -GrainUnitVector(9*MyOrientation + 5);
                        }

                        double Angle3 = (GrainUnitVector(9*MyOrientation + 6)*x0 + GrainUnitVector(9*MyOrientation + 7)*y0 + GrainUnitVector(9*MyOrientation + 8)*z0)/mag0;
                        if (Angle3 < 0) {
                            Diag3X = GrainUnitVector(9*MyOrientation + 6);
                            Diag3Y = GrainUnitVector(9*MyOrientation + 7);
                            Diag3Z = GrainUnitVector(9*MyOrientation + 8);
                        }
                        else {
                            Diag3X = -GrainUnitVector(9*MyOrientation + 6);
                            Diag3Y = -GrainUnitVector(9*MyOrientation + 7);
                            Diag3Z = -GrainUnitVector(9*MyOrientation + 8);
                        }
                        double U1[3], U2[3], UU[3], Norm[3];
                        U1[0] = Diag2X - Diag1X;
                        U1[1] = Diag2Y - Diag1Y;
                        U1[2] = Diag2Z - Diag1Z;
                        U2[0] = Diag3X - Diag1X;
                        U2[1] = Diag3Y - Diag1Y;
                        U2[2] = Diag3Z - Diag1Z;
                        UU[0] = U1[1]*U2[2] - U1[2]*U2[1];
                        UU[1] = U1[2]*U2[0] - U1[0]*U2[2];
                        UU[2] = U1[0]*U2[1] - U1[1]*U2[0];
                        double NDem = sqrt(UU[0]*UU[0] + UU[1]*UU[1] + UU[2]*UU[2]);
                        Norm[0] = UU[0]/NDem;
                        Norm[1] = UU[1]/NDem;
                        Norm[2] = UU[2]/NDem;
                        // normal to capturing plane
                        double normx = Norm[0];
                        double normy = Norm[1];
                        double normz = Norm[2];
                        double ParaT = (normx*x0+normy*y0+normz*z0)/(normx*Diag1X+normy*Diag1Y+normz*Diag1Z);
                        float CDLVal = pow(pow(ParaT*Diag1X,2.0) + pow(ParaT*Diag1Y,2.0) + pow(ParaT*Diag1Z,2.0),0.5);
                        CritDiagonalLength((long int)(26)*CellLocation+(long int)(n)) = CDLVal;
                    }
                    CellType(GlobalCellLocation) = Active;
                }
            });
        }
    }
    
    // Wait on send requests
    MPI_Waitall(2, SendRequests.data(), MPI_STATUSES_IGNORE );
}
