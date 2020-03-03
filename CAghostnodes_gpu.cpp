#include "header.h"
using namespace std;

// 2D domain decomposition: update ghost nodes with new cell data from Nucleation and CellCapture routines
void GhostNodes2D_GPU(int cycle, int id, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nz, int NeighborX[26], int NeighborY[26], int NeighborZ[26], ViewI CellType, ViewF DOCenter, ViewI GrainID, float* GrainUnitVector, ViewI TriangleIndex, int* GrainOrientation, ViewF DiagonalLength, ViewF CritDiagonalLength, int NGrainOrientations, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferC, Buffer2D BufferD, Buffer2D BufferE, Buffer2D BufferF, Buffer2D BufferG, Buffer2D BufferH, Buffer2D BufferAR, Buffer2D BufferBR, Buffer2D BufferCR, Buffer2D BufferDR, Buffer2D BufferER, Buffer2D BufferFR, Buffer2D BufferGR, Buffer2D BufferHR, int BufSizeX, int BufSizeY, int BufSizeZ, ViewI Locks) {

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
                int RankX, RankY, RankZ, NewGrainID, CellLocation;
                bool Place = false;
                float DOCenterX, DOCenterY, DOCenterZ, NewDiagonalLength;
                // Which rank was the data received from?
                if (unpack_index == 0) {
                    // Data receieved from MyLeft
                    // Adjust X Position by +1 since X = 0 is not included in this buffer
                    RankX = BufPosition % BufSizeX + 1;
                    RankY = 0;
                    RankZ = floor(BufPosition/BufSizeX);
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
                    RankZ = floor(BufPosition/BufSizeX);
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
                    RankZ = floor(BufPosition/BufSizeY);
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
                    RankZ = floor(BufPosition/BufSizeY);
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
                    //printf("unpack, X/Y/Z, GID, Center, DL %d %d %d %d %d %f %f %f %f \n",unpack_index,CellLocation,RankX,RankY,RankZ,DOCenterX,DOCenterY,DOCenterZ,NewDiagonalLength);
                    CellType(CellLocation) = Active;
                    Locks(CellLocation) = 0;
                    GrainID(CellLocation) = NewGrainID;
                    DOCenter(3*CellLocation) = DOCenterX;
                    DOCenter(3*CellLocation+1) = DOCenterY;
                    DOCenter(3*CellLocation+2) = DOCenterZ;
                    int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
                    DiagonalLength(CellLocation) = NewDiagonalLength;
                    // Global coordinates of cell center
                    double xp = RankX + MyXOffset + 0.5;
                    double yp = RankY + MyYOffset + 0.5;
                    double zp = RankZ + 0.5;
                    // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
                    for (int n=0; n<26; n++)  {
                        
                        int MyNeighborX = RankX + NeighborX[n];
                        int MyNeighborY = RankY + NeighborY[n];
                        int MyNeighborZ = RankZ + NeighborZ[n];
                        int NeighborPosition = MyNeighborZ*MyXSlices*MyYSlices + MyNeighborX*MyYSlices + MyNeighborY;
                        if (NeighborPosition == CellLocation) {
                            // Do not calculate critical diagonal length req'd for the newly captured cell to capture the original
                            TriangleIndex(78*NeighborPosition + 3*n) = 6;
                            TriangleIndex(78*NeighborPosition + 3*n + 1) = 6;
                            TriangleIndex(78*NeighborPosition + 3*n + 2) = 6;
                            CritDiagonalLength(26*NeighborPosition+n) = 10000000.0;
                        }
                        
                        // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the global coordinates of the center of a neighbor cell
                        double x0 = xp + NeighborX[n] - DOCenterX;
                        double y0 = yp + NeighborY[n] - DOCenterY;
                        double z0 = zp + NeighborZ[n] - DOCenterZ;
                        // mag0 is the magnitude of (x0,y0,z0)
                        double mag0 = pow(pow(x0,2) + pow(y0,2) + pow(z0,2),0.5);
                        
                        // Calculate angles between the octahedron diagonal directions and the vector x0,y0,z0
                        double AnglesA[6];
                        for (int aa=0; aa<6; aa++) {
                            double xd = GrainUnitVector[18*MyOrientation + 3*aa];
                            double yd = GrainUnitVector[18*MyOrientation + 3*aa + 1];
                            double zd = GrainUnitVector[18*MyOrientation + 3*aa + 2];
                            AnglesA[aa] = (xd*x0 + yd*y0 + zd*z0)/mag0;
                        }
                        
                        int index1, index2, index3;
                        index1 = 0;
                        for (int ii=1; ii<6; ii++) {
                            if (AnglesA[index1] < AnglesA[ii]) {
                                index1 = ii;
                            }
                        }
                        
                        TriangleIndex(78*CellLocation + 3*n) = index1;
                        // First diagonal of the capturing face is that which makes the smallest (?) angle with x0,y0,z0
                        double Diag1X = GrainUnitVector[18*MyOrientation + 3*index1 + 0];
                        double Diag1Y = GrainUnitVector[18*MyOrientation + 3*index1 + 1];
                        double Diag1Z = GrainUnitVector[18*MyOrientation + 3*index1 + 2];
                        AnglesA[index1] = -1;
                        if (index1 % 2 == 0) AnglesA[index1+1] = -1;
                        if (index1 % 2 == 1) AnglesA[index1-1] = -1;
                        
                        double MaxA = AnglesA[0];
                        for (int ii=1; ii<6; ii++) {
                            if (MaxA < AnglesA[ii]) {
                                MaxA = AnglesA[ii];
                            }
                        }
                        
                        
                        double Diag2X, Diag2Y, Diag2Z, Diag3X, Diag3Y, Diag3Z;
                        if (MaxA == 0) {
                            // Special case- other diagonals are all perpendicular to the first one (e.g. the octahedron corner captures the new cell center)
                            // manually assign other diagonals (one of the 4 possible "capturing" faces)
                            if ((index1 == 0)||(index1 == 1)) {
                                index2 = 2;
                                index3 = 4;
                            }
                            else if ((index1 == 2)||(index1 == 3)) {
                                index2 = 0;
                                index3 = 4;
                            }
                            else if ((index1 == 4)||(index1 == 5)) {
                                index2 = 0;
                                index3 = 2;
                            }
                        }
                        else {
                            // 2nd and 3rd closest diagonals to x0,y0,z0
                            // note that if these are the same length, this means that the octahedron edge captures the new cell center
                            // in this case, either of the 2 possible "capturing" faces will work
                            index2 = 0;
                            for (int ii=1; ii<6; ii++) {
                                if (AnglesA[index2] < AnglesA[ii]) {
                                    index2 = ii;
                                }
                            }
                            AnglesA[index2] = -1;
                            if (index2 % 2 == 0) AnglesA[index2+1] = -1;
                            if (index2 % 2 == 1) AnglesA[index2-1] = -1;
                            // index3 = MaxIndex(AnglesA);
                            index3 = 0;
                            for (int ii=1; ii<6; ii++) {
                                if (AnglesA[index3] < AnglesA[ii]) {
                                    index3 = ii;
                                }
                            }
                            
                        }
                        TriangleIndex(78*CellLocation + 3*n + 1) = index2;
                        Diag2X = GrainUnitVector[18*MyOrientation + 3*index2 + 0];
                        Diag2Y = GrainUnitVector[18*MyOrientation + 3*index2 + 1];
                        Diag2Z = GrainUnitVector[18*MyOrientation + 3*index2 + 2];
                        TriangleIndex(78*CellLocation + 3*n + 2) = index3;
                        Diag3X = GrainUnitVector[18*MyOrientation + 3*index3 + 0];
                        Diag3Y = GrainUnitVector[18*MyOrientation + 3*index3 + 1];
                        Diag3Z = GrainUnitVector[18*MyOrientation + 3*index3 + 2];
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
                        float CDLVal = pow(pow(ParaT*Diag1X,2) + pow(ParaT*Diag1Y,2) + pow(ParaT*Diag1Z,2),0.5);
                        //                                if ((normx*Diag1X+normy*Diag1Y+normz*Diag1Z) == 0.0) {
                        //                                    printf("Captured cell : %d %d %d %f %d %d %d %f %f %f",MyNeighborX,MyNeighborY,MyNeighborZ,mag0,index1,index2,index3,normx,normy,normz);
                        //                                }
                        CritDiagonalLength(26*CellLocation+n) = CDLVal;
                        //printf("Unpack %d CLOC %f %f %f OLOC %f %f %f TI = %d %d %d DiagonalN %d CDL Val %f \n",unpack_index,RankX+MyXOffset+0.5,RankY+MyYOffset+0.5,RankZ+0.5,DOCenterX,DOCenterY,DOCenterZ,TriangleIndex(78*CellLocation + 3*n),TriangleIndex(78*CellLocation + 3*n + 1),TriangleIndex(78*CellLocation + 3*n+2),n,CDLVal);
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
void GhostNodes1D_GPU(int cycle, int id, int MyLeft, int MyRight, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nz, int NeighborX[26], int NeighborY[26], int NeighborZ[26], ViewI CellType, ViewF DOCenter, ViewI GrainID, float* GrainUnitVector, ViewI TriangleIndex, int* GrainOrientation, ViewF DiagonalLength, ViewF CritDiagonalLength, int NGrainOrientations, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferAR, Buffer2D BufferBR, int BufSizeX, int BufSizeY, int BufSizeZ, ViewI Locks) {
    
    Kokkos::fence();
    MPI_Barrier(MPI_COMM_WORLD);
    
    std::vector<MPI_Request> SendRequests(2, MPI_REQUEST_NULL);
    std::vector<MPI_Request> RecvRequests(2, MPI_REQUEST_NULL);
    
    // Send data to each other rank (MPI_Isend)
    MPI_Isend(BufferA.data(),5*BufSizeX*BufSizeZ,MPI_FLOAT,MyLeft,0,MPI_COMM_WORLD,&SendRequests[0]);
    MPI_Isend(BufferB.data(),5*BufSizeX*BufSizeZ,MPI_FLOAT,MyRight,0,MPI_COMM_WORLD,&SendRequests[1]);
    
    // Receive buffers for all neighbors (MPI_Irecv)
    MPI_Irecv(BufferAR.data(),5*BufSizeX*BufSizeZ,MPI_FLOAT,MyLeft,0,MPI_COMM_WORLD,&RecvRequests[0]);
    MPI_Irecv(BufferBR.data(),5*BufSizeX*BufSizeZ,MPI_FLOAT,MyRight,0,MPI_COMM_WORLD,&RecvRequests[1]);
    
    //MPI_Barrier(MPI_COMM_WORLD);
    //if (id == 0) cout << "Data transfer complete" << endl;
    
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
                int RankX, RankY, RankZ, NewGrainID, CellLocation;
                float DOCenterX, DOCenterY, DOCenterZ, NewDiagonalLength;
                bool Place = false;
                RankZ = floor(BufPosition/BufSizeX);
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
                    //if (BufferAR(BufPosition,4) > 0) printf("NewDL for A is %f \n",BufferAR(BufPosition,4));
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
                    //if (BufferBR(BufPosition,4) > 0) printf("NewDL for B is %f \n",BufferBR(BufPosition,4));
                }
                if (Place) {
                    //printf("ID, RankX, RankZ, SizeX, SizeZ %d %d %d %d %d \n",id,RankX,RankZ,BufSizeX,BufSizeZ);
                    //printf("ID, X/Y/Z, GID, Center, DL %d %d %d %d %d %f %f %f %f \n",id,CellLocation,RankX,RankY,RankZ,DOCenterX,DOCenterY,DOCenterZ,NewDiagonalLength);
                    // Update this ghost node cell's information with data from other rank
                    GrainID(CellLocation) = NewGrainID;
                    Locks(CellLocation) = 0;
                    DOCenter(3*CellLocation) = DOCenterX;
                    DOCenter(3*CellLocation+1) = DOCenterY;
                    DOCenter(3*CellLocation+2) = DOCenterZ;
                    int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
                    DiagonalLength(CellLocation) = NewDiagonalLength;
                    // Global coordinates of cell center
                    double xp = RankX + MyXOffset + 0.5;
                    double yp = RankY + MyYOffset + 0.5;
                    double zp = RankZ + 0.5;
                    // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
                    for (int n=0; n<26; n++)  {
                        
                        int MyNeighborX = RankX + NeighborX[n];
                        int MyNeighborY = RankY + NeighborY[n];
                        int MyNeighborZ = RankZ + NeighborZ[n];
                        int NeighborPosition = MyNeighborZ*MyXSlices*MyYSlices + MyNeighborX*MyYSlices + MyNeighborY;
                        if (NeighborPosition == CellLocation) {
                            // Do not calculate critical diagonal length req'd for the newly captured cell to capture the original
                            TriangleIndex(78*NeighborPosition + 3*n) = 6;
                            TriangleIndex(78*NeighborPosition + 3*n + 1) = 6;
                            TriangleIndex(78*NeighborPosition + 3*n + 2) = 6;
                            CritDiagonalLength(26*NeighborPosition+n) = 10000000.0;
                        }
                        
                        // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the global coordinates of the center of a neighbor cell
                        double x0 = xp + NeighborX[n] - DOCenterX;
                        double y0 = yp + NeighborY[n] - DOCenterY;
                        double z0 = zp + NeighborZ[n] - DOCenterZ;
                        // mag0 is the magnitude of (x0,y0,z0)
                        double mag0 = pow(pow(x0,2) + pow(y0,2) + pow(z0,2),0.5);
                        
                        // Calculate angles between the octahedron diagonal directions and the vector x0,y0,z0
                        double AnglesA[6];
                        for (int aa=0; aa<6; aa++) {
                            double xd = GrainUnitVector[18*MyOrientation + 3*aa];
                            double yd = GrainUnitVector[18*MyOrientation + 3*aa + 1];
                            double zd = GrainUnitVector[18*MyOrientation + 3*aa + 2];
                            AnglesA[aa] = (xd*x0 + yd*y0 + zd*z0)/mag0;
                        }
                        
                        int index1, index2, index3;
                        index1 = 0;
                        for (int ii=1; ii<6; ii++) {
                            if (AnglesA[index1] < AnglesA[ii]) {
                                index1 = ii;
                            }
                        }
                        
                        TriangleIndex(78*CellLocation + 3*n) = index1;
                        // First diagonal of the capturing face is that which makes the smallest (?) angle with x0,y0,z0
                        double Diag1X = GrainUnitVector[18*MyOrientation + 3*index1 + 0];
                        double Diag1Y = GrainUnitVector[18*MyOrientation + 3*index1 + 1];
                        double Diag1Z = GrainUnitVector[18*MyOrientation + 3*index1 + 2];
                        AnglesA[index1] = -1;
                        if (index1 % 2 == 0) AnglesA[index1+1] = -1;
                        if (index1 % 2 == 1) AnglesA[index1-1] = -1;
                        
                        double MaxA = AnglesA[0];
                        for (int ii=1; ii<6; ii++) {
                            if (MaxA < AnglesA[ii]) {
                                MaxA = AnglesA[ii];
                            }
                        }
                        
                        
                        double Diag2X, Diag2Y, Diag2Z, Diag3X, Diag3Y, Diag3Z;
                        if (MaxA == 0) {
                            // Special case- other diagonals are all perpendicular to the first one (e.g. the octahedron corner captures the new cell center)
                            // manually assign other diagonals (one of the 4 possible "capturing" faces)
                            if ((index1 == 0)||(index1 == 1)) {
                                index2 = 2;
                                index3 = 4;
                            }
                            else if ((index1 == 2)||(index1 == 3)) {
                                index2 = 0;
                                index3 = 4;
                            }
                            else if ((index1 == 4)||(index1 == 5)) {
                                index2 = 0;
                                index3 = 2;
                            }
                        }
                        else {
                            // 2nd and 3rd closest diagonals to x0,y0,z0
                            // note that if these are the same length, this means that the octahedron edge captures the new cell center
                            // in this case, either of the 2 possible "capturing" faces will work
                            index2 = 0;
                            for (int ii=1; ii<6; ii++) {
                                if (AnglesA[index2] < AnglesA[ii]) {
                                    index2 = ii;
                                }
                            }
                            AnglesA[index2] = -1;
                            if (index2 % 2 == 0) AnglesA[index2+1] = -1;
                            if (index2 % 2 == 1) AnglesA[index2-1] = -1;
                            // index3 = MaxIndex(AnglesA);
                            index3 = 0;
                            for (int ii=1; ii<6; ii++) {
                                if (AnglesA[index3] < AnglesA[ii]) {
                                    index3 = ii;
                                }
                            }
                            
                        }
                        TriangleIndex(78*CellLocation + 3*n + 1) = index2;
                        Diag2X = GrainUnitVector[18*MyOrientation + 3*index2 + 0];
                        Diag2Y = GrainUnitVector[18*MyOrientation + 3*index2 + 1];
                        Diag2Z = GrainUnitVector[18*MyOrientation + 3*index2 + 2];
                        TriangleIndex(78*CellLocation + 3*n + 2) = index3;
                        Diag3X = GrainUnitVector[18*MyOrientation + 3*index3 + 0];
                        Diag3Y = GrainUnitVector[18*MyOrientation + 3*index3 + 1];
                        Diag3Z = GrainUnitVector[18*MyOrientation + 3*index3 + 2];
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
                        float CDLVal = pow(pow(ParaT*Diag1X,2) + pow(ParaT*Diag1Y,2) + pow(ParaT*Diag1Z,2),0.5);
                        //                                if ((normx*Diag1X+normy*Diag1Y+normz*Diag1Z) == 0.0) {
                        //                                    printf("Captured cell : %d %d %d %f %d %d %d %f %f %f",MyNeighborX,MyNeighborY,MyNeighborZ,mag0,index1,index2,index3,normx,normy,normz);
                        //                                }
                        CritDiagonalLength(26*CellLocation+n) = CDLVal;
                        //if (id == 1) printf("Cell at CYCLE %d LOC %d %d %d ORIENTATION %d DL = %f TI = %d %d %d DiagonalN %d CDL Val %f \n",cycle,RankX,RankY,RankZ,GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)],DiagonalLength(CellLocation),TriangleIndex(78*CellLocation + 3*n),TriangleIndex(78*CellLocation + 3*n + 1),TriangleIndex(78*CellLocation + 3*n+2),n,CDLVal);
                    }
                    CellType(CellLocation) = Active;
                }
            });
        }
    }
    // Wait on send requests
    MPI_Waitall(2, SendRequests.data(), MPI_STATUSES_IGNORE );
    
}
