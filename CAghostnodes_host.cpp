#include "header.h"
using namespace std;

// 2D domain decomposition: update ghost nodes with new cell data from Nucleation and CellCapture routines
void GhostNodes2D(int cycle, int id, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nz, int NeighborX[26], int NeighborY[26], int NeighborZ[26], ViewI::HostMirror CellType, ViewF::HostMirror DOCenter, ViewI::HostMirror GrainID, float* GrainUnitVector, ViewI::HostMirror TriangleIndex, int* GrainOrientation, ViewF::HostMirror DiagonalLength, ViewF::HostMirror CritDiagonalLength, int NGrainOrientations) {

    
    int ACount = 0;
    int BCount = 0;
    int CCount = 0;
    int DCount = 0;
    int ECount = 0;
    int FCount = 0;
    int GCount = 0;
    int HCount = 0;
    
    int LocalDomainSize = nz*MyXSlices*MyYSlices;
    //Kokkos::parallel_for("GhostNodesCount",LocalDomainSize, KOKKOS_LAMBDA (const int& D3D1ConvPosition) {
    for (int D3D1ConvPosition=0; D3D1ConvPosition<LocalDomainSize; D3D1ConvPosition++) {
        if ((CellType(D3D1ConvPosition) == Ghost1)||(CellType(D3D1ConvPosition) == Ghost2)||(CellType(D3D1ConvPosition) == Ghost3)) {
            // This cell is at (RankZ,RankX,RankY)
            int RankZ = floor(D3D1ConvPosition/(MyXSlices*MyYSlices));
            int Rem = D3D1ConvPosition % (MyXSlices*MyYSlices);
            int RankX = floor(Rem/MyYSlices);
            int RankY = Rem % MyYSlices;
            if (RankY == 1) {
                // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                if (RankX == MyXSlices-2) {
                    ECount++;
                    CCount++;
                    ACount++;
                }
                else if (RankX == 1) {
                    GCount++;
                    DCount++;
                    ACount++;
                }
                else if ((RankX > 1)&&(RankX < MyXSlices-2)) {
                    ACount++;
                }
            }
            else if (RankY == MyYSlices-2) {
                // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                if (RankX == MyXSlices-2) {
                    FCount++;
                    CCount++;
                    BCount++;
                }
                else if (RankX == 1) {
                    HCount++;
                    DCount++;
                    BCount++;
                }
                else if ((RankX > 1)&&(RankX < MyXSlices-2)) {
                    BCount++;
                }
            }
            else if ((RankX == 1)&&(RankY > 1)&&(RankY < MyYSlices-2)) {
                DCount++;
            }
            else if ((RankX == MyXSlices-2)&&(RankY > 1)&&(RankY < MyYSlices-2)) {
                CCount++;
            }
        }
    }

    // Determine whether or not ghost node information transfer needs to take place
    int ARCount, BRCount, CRCount, DRCount, ERCount, FRCount, GRCount, HRCount;

    // Send BCount, Recieve ARCount (send to the right, recieve on the left)
    MPI_Sendrecv(&BCount,1,MPI_INT,MyRight,0,&ARCount,1,MPI_INT,MyLeft,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Send ACount, Recieve BRCount (send to the left, recieve on the right)
    MPI_Sendrecv(&ACount,1,MPI_INT,MyLeft,0,&BRCount,1,MPI_INT,MyRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Send CCount, Recieve DRCount (send into the plane, recieve out of the plane)
    MPI_Sendrecv(&CCount,1,MPI_INT,MyIn,0,&DRCount,1,MPI_INT,MyOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Send DCount, Recieve CRCount (send out of the plane, recieve into the plane)
    MPI_Sendrecv(&DCount,1,MPI_INT,MyOut,0,&CRCount,1,MPI_INT,MyIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);


    // Send HCount, Recieve ERCount
    MPI_Sendrecv(&HCount,1,MPI_INT,MyRightOut,0,&ERCount,1,MPI_INT,MyLeftIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Send ECount, Recieve HRCount
    MPI_Sendrecv(&ECount,1,MPI_INT,MyLeftIn,0,&HRCount,1,MPI_INT,MyRightOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Send GCount, Recieve FRCount
    MPI_Sendrecv(&GCount,1,MPI_INT,MyLeftOut,0,&FRCount,1,MPI_INT,MyRightIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Send FCount, Recieve GRCount
    MPI_Sendrecv(&FCount,1,MPI_INT,MyRightIn,0,&GRCount,1,MPI_INT,MyLeftOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    if (MyLeft == MPI_PROC_NULL) ARCount = 0;
    if (MyRight == MPI_PROC_NULL) BRCount = 0;
    if (MyIn == MPI_PROC_NULL) CRCount = 0;
    if (MyOut == MPI_PROC_NULL) DRCount = 0;
    if (MyLeftIn == MPI_PROC_NULL) ERCount = 0;
    if (MyRightIn == MPI_PROC_NULL) FRCount = 0;
    if (MyLeftOut == MPI_PROC_NULL) GRCount = 0;
    if (MyRightOut == MPI_PROC_NULL) HRCount = 0;


    // Buffers for recieving ghost node data
    float GhostNodesAR[7*ARCount];
    float GhostNodesBR[7*BRCount];
    float GhostNodesCR[7*CRCount];
    float GhostNodesDR[7*DRCount];
    float GhostNodesER[6*ERCount];
    float GhostNodesFR[6*FRCount];
    float GhostNodesGR[6*GRCount];
    float GhostNodesHR[6*HRCount];


    // Collect ghost node data and send to other ranks- left and right
    if (ACount > 0) {
        float GhostNodesA[7*ACount];
        int BufACount = 0;
        for (int RankZ=1; RankZ<nz-1; RankZ++) {
            for (int RankX=1; RankX<MyXSlices-1; RankX++) {
                int CellPosition = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + 1;
                if ((CellType(CellPosition) == Ghost1)||(CellType(CellPosition) == Ghost2)||(CellType(CellPosition) == Ghost3)) {
                    // This cell just became active, store its information to send to id = id - 1
                    GhostNodesA[7*BufACount] = (float)(RankX);
                    GhostNodesA[7*BufACount+1] = (float)(RankZ);
                    GhostNodesA[7*BufACount+2] = (float)(GrainID(CellPosition));
                    GhostNodesA[7*BufACount+3] = DOCenter(3*CellPosition);
                    GhostNodesA[7*BufACount+4] = DOCenter(3*CellPosition+1);
                    GhostNodesA[7*BufACount+5] = DOCenter(3*CellPosition+2);
                    GhostNodesA[7*BufACount+6] = DiagonalLength(CellPosition);
                    if (CellType(CellPosition) == Ghost3) CellType(CellPosition) = Ghost2;
                    else if (CellType(CellPosition) == Ghost2) CellType(CellPosition) = Ghost1;
                    else if (CellType(CellPosition) == Ghost1) CellType(CellPosition) = Active;
                    BufACount++;
                }
            }
        }
        if (BRCount == 0) {
            // Sending data to id = id - 1 only
            MPI_Send(&GhostNodesA,ACount*7,MPI_FLOAT,MyLeft,0,MPI_COMM_WORLD);
        }
        else {
            // Sending data to id = id - 1 and recieving data from id = id + 1
            MPI_Sendrecv(&GhostNodesA,ACount*7,MPI_FLOAT,MyLeft,0,&GhostNodesBR,BRCount*7,MPI_FLOAT,MyRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
       // cout << "ID = " << id << " ACount = " << ACount << " ABuf = " << BufACount << endl;
    }
    else if (BRCount > 0) {
        // Recieving data from id = id + 1 only
        MPI_Recv(&GhostNodesBR,BRCount*7,MPI_FLOAT,MyRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }

    if (BCount > 0) {
        float GhostNodesB[7*BCount];
        int BufBCount = 0;
        for (int RankZ=1; RankZ<nz-1; RankZ++) {
            for (int RankX=1; RankX<MyXSlices-1; RankX++) {
                int CellPosition = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + MyYSlices-2;
                if ((CellType(CellPosition) == Ghost1)||(CellType(CellPosition) == Ghost2)||(CellType(CellPosition) == Ghost3)) {
                    // This cell just became active, store its information to send to id = id - 1
                    GhostNodesB[7*BufBCount] = (float)(RankX);
                    GhostNodesB[7*BufBCount+1] = (float)(RankZ);
                    GhostNodesB[7*BufBCount+2] = (float)(GrainID(CellPosition));
                    GhostNodesB[7*BufBCount+3] = DOCenter(3*CellPosition);
                    GhostNodesB[7*BufBCount+4] = DOCenter(3*CellPosition+1);
                    GhostNodesB[7*BufBCount+5] = DOCenter(3*CellPosition+2);
                    GhostNodesB[7*BufBCount+6] = DiagonalLength(CellPosition);
                    if (CellType(CellPosition) == Ghost3) CellType(CellPosition) = Ghost2;
                    else if (CellType(CellPosition) == Ghost2) CellType(CellPosition) = Ghost1;
                    else if (CellType(CellPosition) == Ghost1) CellType(CellPosition) = Active;
                    BufBCount++;
                    //cout << "CT after B: " << CellType[RankZ][RankX][MyYSlices-2] << endl;
                }
            }
        }
      //  cout << "ID = " << id << " BCount = " << BCount << " BBuf = " << BufBCount << endl;
        if (ARCount == 0) {
            // Sending data to id = id + 1 only
            MPI_Send(&GhostNodesB,BCount*7,MPI_FLOAT,MyRight,1,MPI_COMM_WORLD);
        }
        else {
            // Sending data to id = id + 1 and recieving data from id = id - 1
            MPI_Sendrecv(&GhostNodesB,BCount*7,MPI_FLOAT,MyRight,1,&GhostNodesAR,ARCount*7,MPI_FLOAT,MyLeft,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    }
    else if (ARCount > 0) {
        // Recieving data from id = id - 1 only
        MPI_Recv(&GhostNodesAR,ARCount*7,MPI_FLOAT,MyLeft,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }

    // Collect ghost node data and send to other ranks- in and out
    if (CCount > 0) {
        float GhostNodesC[7*CCount];
        int BufCCount = 0;
        for (int RankZ=1; RankZ<nz-1; RankZ++) {
            for (int RankY=1; RankY<MyYSlices-1; RankY++) {
                int CellLocation = RankZ*MyXSlices*MyYSlices +(MyXSlices-2)*MyYSlices + RankY;
                if ((CellType(CellLocation) == Ghost3)||(CellType(CellLocation) == Ghost2) ||(CellType(CellLocation) == Ghost1)) {
                    // This cell just became active, store its information
                    GhostNodesC[7*BufCCount] = (float)(RankY);
                    GhostNodesC[7*BufCCount+1] = (float)(RankZ);
                    GhostNodesC[7*BufCCount+2] = (float)(GrainID(CellLocation));
                    GhostNodesC[7*BufCCount+3] = DOCenter(3*CellLocation);
                    GhostNodesC[7*BufCCount+4] = DOCenter(3*CellLocation+1);
                    GhostNodesC[7*BufCCount+5] = DOCenter(3*CellLocation+2);
                    GhostNodesC[7*BufCCount+6] = DiagonalLength(CellLocation);
                    if (CellType(CellLocation) == Ghost3) CellType(CellLocation) = Ghost2;
                    else if (CellType(CellLocation) == Ghost2) CellType(CellLocation) = Ghost1;
                    else if (CellType(CellLocation) == Ghost1) CellType(CellLocation) = Active;
                    BufCCount++;
                }
            }
        }

        if (DRCount == 0) {
            // Sending data only
            MPI_Send(&GhostNodesC,CCount*7,MPI_FLOAT,MyIn,0,MPI_COMM_WORLD);
        }
        else {
            // Sending data and recieving data
            MPI_Sendrecv(&GhostNodesC,CCount*7,MPI_FLOAT,MyIn,0,&GhostNodesDR,DRCount*7,MPI_FLOAT,MyOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
       // cout << "ID = " << id << " CCount = " << CCount << " CBuf = " << BufCCount << endl;
    }
    else if (DRCount > 0) {
        // Recieving data only
        MPI_Recv(&GhostNodesDR,DRCount*7,MPI_FLOAT,MyOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }

    if (DCount > 0) {
        float GhostNodesD[7*DCount];
        int BufDCount = 0;
        for (int RankZ=1; RankZ<nz-1; RankZ++) {
            for (int RankY=1; RankY<MyYSlices-1; RankY++) {
                int CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices + RankY;
                if ((CellType(CellLocation) == Ghost3)||(CellType(CellLocation) == Ghost2)||(CellType(CellLocation) == Ghost1)) {
                    // This cell just became active, store its information
                    GhostNodesD[7*BufDCount] = (float)(RankY);
                    GhostNodesD[7*BufDCount+1] = (float)(RankZ);
                    GhostNodesD[7*BufDCount+2] = (float)(GrainID(CellLocation));
                    GhostNodesD[7*BufDCount+3] = DOCenter(3*CellLocation);
                    GhostNodesD[7*BufDCount+4] = DOCenter(3*CellLocation+1);
                    GhostNodesD[7*BufDCount+5] = DOCenter(3*CellLocation+2);
                    GhostNodesD[7*BufDCount+6] = DiagonalLength(CellLocation);
                    if (CellType(CellLocation) == Ghost3) CellType(CellLocation) = Ghost2;
                    else if (CellType(CellLocation) == Ghost2) CellType(CellLocation) = Ghost1;
                    else if (CellType(CellLocation) == Ghost1) CellType(CellLocation) = Active;
                    BufDCount++;
                }
            }
        }
        if (CRCount == 0) {
            // Sending data only
            MPI_Send(&GhostNodesD,DCount*7,MPI_FLOAT,MyOut,1,MPI_COMM_WORLD);
        }
        else {
            // Sending data and recieving data
            MPI_Sendrecv(&GhostNodesD,DCount*7,MPI_FLOAT,MyOut,1,&GhostNodesCR,CRCount*7,MPI_FLOAT,MyIn,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    //    cout << "ID = " << id << " DCount = " << DCount << " DBuf = " << BufDCount << endl;
    }
    else if (CRCount > 0) {
        // Recieving data only
        MPI_Recv(&GhostNodesCR,CRCount*7,MPI_FLOAT,MyIn,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }

    // Collect ghost node data and send to other ranks- MyLeftIn and MyRightOut
    if (ECount > 0) {
        float GhostNodesE[6*ECount];
        int BufECount = 0;
        for (int RankZ=1; RankZ<nz-1; RankZ++) {
            int CellLocation = RankZ*MyXSlices*MyYSlices + (MyXSlices-2)*(MyYSlices) + 1;
            if ((CellType(CellLocation) == Ghost3)||(CellType(CellLocation) == Ghost2)||(CellType(CellLocation) == Ghost1)) {
                // This cell just became active, store its information
                GhostNodesE[6*BufECount] = (float)(RankZ);
                GhostNodesE[6*BufECount+1] = (float)(GrainID(CellLocation));
                GhostNodesE[6*BufECount+2] = DOCenter(3*CellLocation);
                GhostNodesE[6*BufECount+3] = DOCenter(3*CellLocation+1);
                GhostNodesE[6*BufECount+4] = DOCenter(3*CellLocation+2);
                GhostNodesE[6*BufECount+5] = DiagonalLength(CellLocation);
                if (CellType(CellLocation) == Ghost3) CellType(CellLocation) = Ghost2;
                else if (CellType(CellLocation) == Ghost2) CellType(CellLocation) = Ghost1;
                else if (CellType(CellLocation) == Ghost1) CellType(CellLocation) = Active;
                BufECount++;
            }
        }
        if (HRCount == 0) {
            // Sending data only
            MPI_Send(&GhostNodesE,ECount*6,MPI_FLOAT,MyLeftIn,0,MPI_COMM_WORLD);
        }
        else {
            // Sending data and recieving data
            MPI_Sendrecv(&GhostNodesE,ECount*6,MPI_FLOAT,MyLeftIn,0,&GhostNodesHR,HRCount*6,MPI_FLOAT,MyRightOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
      //  cout << "ID = " << id << " ECount = " << ECount << " EBuf = " << BufECount << endl;
    }
    else if (HRCount > 0) {
        // Recieving data only
        MPI_Recv(&GhostNodesHR,HRCount*6,MPI_FLOAT,MyRightOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }

    if (HCount > 0) {
        float GhostNodesH[6*HCount];
        int BufHCount = 0;
        for (int RankZ=1; RankZ<nz-1; RankZ++) {
            int CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices + MyYSlices-2;
            if ((CellType(CellLocation) == Ghost3)||(CellType(CellLocation) == Ghost2)||(CellType(CellLocation) == Ghost1)) {
                // This cell just became active, store its information
                GhostNodesH[6*BufHCount] = (float)(RankZ);
                 GhostNodesH[6*BufHCount+1] = (float)(GrainID(CellLocation));
                GhostNodesH[6*BufHCount+2] = DOCenter(3*CellLocation);
                GhostNodesH[6*BufHCount+3] = DOCenter(3*CellLocation+1);
                GhostNodesH[6*BufHCount+4] = DOCenter(3*CellLocation+2);
                GhostNodesH[6*BufHCount+5] = DiagonalLength(CellLocation);
                if (CellType[CellLocation] == Ghost3) CellType(CellLocation) = Ghost2;
                else if (CellType[CellLocation] == Ghost2) CellType(CellLocation) = Ghost1;
                else if (CellType[CellLocation] == Ghost1) CellType(CellLocation) = Active;
                BufHCount++;
            }
        }
        if (ERCount == 0) {
            // Sending data only
            MPI_Send(&GhostNodesH,HCount*6,MPI_FLOAT,MyRightOut,0,MPI_COMM_WORLD);
        }
        else {
            // Sending data and recieving data
            MPI_Sendrecv(&GhostNodesH,HCount*6,MPI_FLOAT,MyRightOut,0,&GhostNodesER,ERCount*6,MPI_FLOAT,MyLeftIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    }
    else if (ERCount > 0) {
        // Recieving data only
        MPI_Recv(&GhostNodesER,ERCount*6,MPI_FLOAT,MyLeftIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }

    // Collect ghost node data and send to other ranks- MyRightIn and MyLeftOut
    if (FCount > 0) {
        float GhostNodesF[6*FCount];
        int BufFCount = 0;
        for (int RankZ=1; RankZ<nz-1; RankZ++) {
            int CellLocation = RankZ*MyXSlices*MyYSlices + (MyXSlices-2)*MyYSlices + MyYSlices-2;
            if ((CellType(CellLocation) == Ghost3)||(CellType(CellLocation) == Ghost2)||(CellType(CellLocation) == Ghost1)) {
                // This cell just became active, store its information
                GhostNodesF[6*BufFCount] = (float)(RankZ);
                GhostNodesF[6*BufFCount+1] = (float)(GrainID(CellLocation));
                GhostNodesF[6*BufFCount+2] = DOCenter(3*CellLocation);
                GhostNodesF[6*BufFCount+3] = DOCenter(3*CellLocation+1);
                GhostNodesF[6*BufFCount+4] = DOCenter(3*CellLocation+2);
                GhostNodesF[6*BufFCount+5] = DiagonalLength(CellLocation);
                if (CellType(CellLocation) == Ghost3) CellType(CellLocation) = Ghost2;
                else if (CellType(CellLocation) == Ghost2) CellType(CellLocation) = Ghost1;
                else if (CellType(CellLocation) == Ghost1) CellType(CellLocation) = Active;
                BufFCount++;
            }
        }
        if (GRCount == 0) {
            // Sending data only
            MPI_Send(&GhostNodesF,FCount*6,MPI_FLOAT,MyRightIn,1,MPI_COMM_WORLD);
        }
        else {
            // Sending data and recieving data
            MPI_Sendrecv(&GhostNodesF,FCount*6,MPI_FLOAT,MyRightIn,1,&GhostNodesGR,GRCount*6,MPI_FLOAT,MyLeftOut,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    }
    else if (GRCount > 0) {
        // Recieving data only
        MPI_Recv(&GhostNodesGR,GRCount*6,MPI_FLOAT,MyLeftOut,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }

    if (GCount > 0) {
        float GhostNodesG[6*GCount];
        int BufGCount = 0;
        for (int RankZ=1; RankZ<nz-1; RankZ++) {
            int CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices + 1;
            if ((CellType(CellLocation) == Ghost3)||(CellType(CellLocation) == Ghost2)||(CellType(CellLocation) == Ghost1)) {
                // This cell just became active, store its information
                GhostNodesG[6*BufGCount] = (float)(RankZ);
                GhostNodesG[6*BufGCount+1] = (float)(GrainID(CellLocation));
                GhostNodesG[6*BufGCount+2] = DOCenter(3*CellLocation);
                GhostNodesG[6*BufGCount+3] = DOCenter(3*CellLocation+1);
                GhostNodesG[6*BufGCount+4] = DOCenter(3*CellLocation+2);
                GhostNodesG[6*BufGCount+5] = DiagonalLength(CellLocation);
                if (CellType(CellLocation) == Ghost3) CellType(CellLocation) = Ghost2;
                if (CellType(CellLocation) == Ghost2) CellType(CellLocation) = Ghost1;
                if (CellType(CellLocation) == Ghost1) CellType(CellLocation) = Active;
                BufGCount++;
            }
        }
        if (FRCount == 0) {
            // Sending data only
            MPI_Send(&GhostNodesG,GCount*6,MPI_FLOAT,MyLeftOut,1,MPI_COMM_WORLD);
        }
        else {
            // Sending data and recieving data
            MPI_Sendrecv(&GhostNodesG,GCount*6,MPI_FLOAT,MyLeftOut,1,&GhostNodesFR,FRCount*6,MPI_FLOAT,MyRightIn,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
       // cout << "ID = " << id << " GCount = " << GCount << " GBuf = " << BufGCount << endl;
    }
    else if (FRCount > 0) {
        // Recieving data only
        MPI_Recv(&GhostNodesFR,FRCount*6,MPI_FLOAT,MyRightIn,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }

     //   if ((id == 169)&&(cycle == 11998)) cout << " data collection end" << endl;
    // Place ghost node data recieved from the left (if needed)
    if (ARCount > 0) {
        for (int i=0; i<ARCount; i++) {
//                cout << "ID = " << id << " " << i << " of " << ARCount-1 << endl;
            // GhostNodesAR: X coord, Z coord, GrainID, DOCenterX/Y/Z, DiagLength;
            int RankX = (int)(GhostNodesAR[7*i]);
            int RankY = 0;
            int RankZ = (int)(GhostNodesAR[7*i+1]);
            int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
            if ((cycle == 0)||((CellType(CellLocation) == Liquid)||(CellType(CellLocation) == Delayed)||(CellType(CellLocation) == LiqSol))) {
                //cout << "ID = " << id << "A" << endl;
                int GlobalX = RankX + MyXOffset;
                int GlobalY = RankY + MyYOffset;
                // Update this ghost node cell's information with data from other rank
                CellType(CellLocation) = Active;
                GrainID(CellLocation) = (int)(GhostNodesAR[7*i+2]);
                double cx = GhostNodesAR[7*i+3];
                double cy = GhostNodesAR[7*i+4];
                double cz = GhostNodesAR[7*i+5];
                DOCenter(3*CellLocation) = cx;
                DOCenter(3*CellLocation+1) = cy;
                DOCenter(3*CellLocation+2) = cz;
                int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
                DiagonalLength(CellLocation) = GhostNodesAR[7*i+6];
                double xp = GlobalX + 0.5;
                double yp = GlobalY + 0.5;
                double zp = RankZ + 0.5;
                // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
                CritDiagLengthCalc(xp, yp, zp, MyOrientation,RankX,RankY,RankZ,CellLocation,cx,cy,cz,NeighborX,NeighborY,NeighborZ, GrainUnitVector,TriangleIndex,CritDiagonalLength);

            }
        }
    }
//        MPI_Barrier(MPI_COMM_WORLD);
//        if (id == 0) cout << " A " << endl;
    //   if ((id == 169)&&(cycle == 11998)) cout << " data collection end A" << endl;
    // Place ghost node data recieved from the right (if needed)
    if (BRCount > 0) {
        for (int i=0; i<BRCount; i++) {
            // GhostNodesBR: X coord, Z coord, GrainID, DOCenterX/Y/Z, DiagLength;
            int RankX = (int)(GhostNodesBR[7*i]);
            int RankY = MyYSlices-1;
            int RankZ = (int)(GhostNodesBR[7*i+1]);
            int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
            if ((cycle == 0)||((CellType(CellLocation) == Liquid)||(CellType(CellLocation) == Delayed)||(CellType(CellLocation) == LiqSol))) {
                int GlobalX = MyXOffset + RankX;
                int GlobalY = MyYOffset + RankY;
                // Update this ghost node cell's information with data from other rank
                CellType(CellLocation) = Active;
                GrainID(CellLocation) = (int)(GhostNodesBR[7*i+2]);
                double cx = GhostNodesBR[7*i+3];
                double cy = GhostNodesBR[7*i+4];
                double cz = GhostNodesBR[7*i+5];
                //if (id == 2) cout << "RECIEVING CYCLE " << cycle << " Center " << cx << " " << cy << " " << cz << endl;
                DOCenter(3*CellLocation) = cx;
                DOCenter(3*CellLocation+1) = cy;
                DOCenter(3*CellLocation+2) = cz;
                int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
                DiagonalLength(CellLocation) = GhostNodesBR[7*i+6];
                double xp = GlobalX + 0.5;
                double yp = GlobalY + 0.5;
                double zp = RankZ + 0.5;
                // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
                CritDiagLengthCalc(xp, yp, zp, MyOrientation,RankX,RankY,RankZ,CellLocation,cx,cy,cz,NeighborX,NeighborY,NeighborZ, GrainUnitVector,TriangleIndex,CritDiagonalLength);
            }
        }
    }
//        MPI_Barrier(MPI_COMM_WORLD);
//           if (id == 0) cout << " data collection end B" << endl;
    // Place ghost node data recieved from in plane (if needed)
    if (CRCount > 0) {
        //cout << "CYCLE = " << cycle << " ID = " << id << " CRCount = " << CRCount << endl;
        for (int i=0; i<CRCount; i++) {
            //cout << "i = " << i << " GNCR = " << GhostNodesCR[0] << endl;
            // GhostNodesCR: Y coord, Z coord, GrainID, DOCenterX/Y/Z, DiagLength;
            int RankX = MyXSlices-1;
            int RankY = (int)(GhostNodesCR[7*i]);
            int RankZ = (int)(GhostNodesCR[7*i+1]);
            int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
            //cout << "ID = " << id << " C checking " << RankX << " " << RankY << " " << RankZ << endl;
            if ((cycle == 0)||((CellType(CellLocation) == Liquid)||(CellType(CellLocation) == Delayed)||(CellType(CellLocation) == LiqSol))) {
                //cout << "ID = " << id << "C" << endl;
                int GlobalX = MyXOffset + RankX;
                int GlobalY = MyYOffset + RankY;
                // Update this ghost node cell's information with data from other rank
                CellType(CellLocation) = Active;
                GrainID(CellLocation) = (int)(GhostNodesCR[7*i+2]);
                double cx = GhostNodesCR[7*i+3];
                double cy = GhostNodesCR[7*i+4];
                double cz = GhostNodesCR[7*i+5];
                DOCenter(3*CellLocation) = cx;
                DOCenter(3*CellLocation+1) = cy;
                DOCenter(3*CellLocation+2) = cz;
                int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
                DiagonalLength(CellLocation) = GhostNodesCR[7*i+6];
                double xp = GlobalX + 0.5;
                double yp = GlobalY + 0.5;
                double zp = RankZ + 0.5;
                // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
                CritDiagLengthCalc(xp, yp, zp, MyOrientation,RankX,RankY,RankZ,CellLocation,cx,cy,cz,NeighborX,NeighborY,NeighborZ, GrainUnitVector,TriangleIndex,CritDiagonalLength);
            }
        }
    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (id == 0) cout << "End data collection C" << endl;
//    }
    // Place ghost node data recieved from out of plane (if needed)
    if (DRCount > 0) {
        for (int i=0; i<DRCount; i++) {
            // GhostNodesDR: Y coord, Z coord, GrainID, DOCenterX/Y/Z, DiagLength;
            int RankX = 0;
            int RankY = (int)(GhostNodesDR[7*i]);
            int RankZ = (int)(GhostNodesDR[7*i+1]);
            int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
            //cout << "ID = " << id << " D checking " << RankX << " " << RankY << " " << RankZ << endl;
            if ((cycle == 0)||((CellType(CellLocation) == Liquid)||(CellType(CellLocation) == Delayed)||(CellType(CellLocation) == LiqSol))) {
                int GlobalX = MyXOffset + RankX;
                int GlobalY = MyYOffset + RankY;
                // Update this ghost node cell's information with data from other rank
                CellType(CellLocation) = Active;
                GrainID(CellLocation) = (int)(GhostNodesDR[7*i+2]);
                double cx = GhostNodesDR[7*i+3];
                double cy = GhostNodesDR[7*i+4];
                double cz = GhostNodesDR[7*i+5];
                DOCenter(3*CellLocation) = cx;
                DOCenter(3*CellLocation+1) = cy;
                DOCenter(3*CellLocation+2) = cz;
                int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
                DiagonalLength(CellLocation) = GhostNodesDR[7*i+6];
                double xp = GlobalX + 0.5;
                double yp = GlobalY + 0.5;
                double zp = RankZ + 0.5;
                // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
                CritDiagLengthCalc(xp, yp, zp, MyOrientation,RankX,RankY,RankZ,CellLocation,cx,cy,cz,NeighborX,NeighborY,NeighborZ, GrainUnitVector,TriangleIndex,CritDiagonalLength);
            }
        }
    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (id == 0) cout << "End data collection D" << endl;
    // Place ghost node data recieved from left and out of plane (if needed)
    if (ERCount > 0) {
        for (int i=0; i<ERCount; i++) {
            // GhostNodesER: Z coord, GrainID, DOCenterX/Y/Z, DiagLength;
            int RankX = MyXSlices-1;
            int RankY = 0;
            int RankZ = (int)(GhostNodesER[6*i]);
            int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
            if ((cycle == 0)||((CellType(CellLocation) == Liquid)||(CellType(CellLocation) == Delayed)||(CellType(CellLocation) == LiqSol))) {
                int GlobalX = MyXOffset + RankX;
                int GlobalY = MyYOffset + RankY;
                // Update this ghost node cell's information with data from other rank
                CellType(CellLocation) = Active;
                GrainID(CellLocation) = (int)(GhostNodesER[6*i+1]);
                double cx = GhostNodesER[6*i+2];
                double cy = GhostNodesER[6*i+3];
                double cz = GhostNodesER[6*i+4];
                DOCenter(3*CellLocation) = cx;
                DOCenter(3*CellLocation+1) = cy;
                DOCenter(3*CellLocation+2) = cz;
                int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
                DiagonalLength(CellLocation) = GhostNodesER[6*i+5];
                double xp = GlobalX + 0.5;
                double yp = GlobalY + 0.5;
                double zp = RankZ + 0.5;
                // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
                CritDiagLengthCalc(xp, yp, zp, MyOrientation,RankX,RankY,RankZ,CellLocation,cx,cy,cz,NeighborX,NeighborY,NeighborZ, GrainUnitVector,TriangleIndex,CritDiagonalLength);
            }
        }
    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (id == 0) cout << "End data collection E" << endl;
    // Place ghost node data recieved from right and out of plane (if needed)
    if (FRCount > 0) {
        for (int i=0; i<FRCount; i++) {
            // GhostNodesFR: Z coord, GrainID, DOCenterX/Y/Z, DiagLength;
            int RankX = MyXSlices-1;
            int RankY = MyYSlices-1;
            int RankZ = (int)(GhostNodesFR[6*i]);
            int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
            if ((cycle == 0)||((CellType(CellLocation) == Liquid)||(CellType(CellLocation) == Delayed)||(CellType(CellLocation) == LiqSol))) {
                int GlobalX = MyXOffset + RankX;
                int GlobalY = MyYOffset + RankY;
                // Update this ghost node cell's information with data from other rank
                CellType(CellLocation) = Active;
                GrainID(CellLocation) = (int)(GhostNodesFR[6*i+1]);
                double cx = GhostNodesFR[6*i+2];
                double cy = GhostNodesFR[6*i+3];
                double cz = GhostNodesFR[6*i+4];
                DOCenter(3*CellLocation) = cx;
                DOCenter(3*CellLocation+1) = cy;
                DOCenter(3*CellLocation+2) = cz;
                int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
                DiagonalLength(CellLocation) = GhostNodesFR[6*i+5];
                double xp = GlobalX + 0.5;
                double yp = GlobalY + 0.5;
                double zp = RankZ + 0.5;
                // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
                CritDiagLengthCalc(xp, yp, zp, MyOrientation,RankX,RankY,RankZ,CellLocation,cx,cy,cz,NeighborX,NeighborY,NeighborZ, GrainUnitVector,TriangleIndex,CritDiagonalLength);
            }
        }
    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (id == 0) cout << "End data collection F" << endl;
    // Place ghost node data recieved from left and into plane (if needed)
    if (GRCount > 0) {
        for (int i=0; i<GRCount; i++) {
            // GhostNodesFR: Z coord, GrainID, DOCenterX/Y/Z, DiagLength;
            int RankX = 0;
            int RankY = 0;
            int RankZ = (int)(GhostNodesGR[6*i]);
            int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
            if ((cycle == 0)||((CellType(CellLocation) == Liquid)||(CellType(CellLocation) == Delayed)||(CellType(CellLocation) == LiqSol))) {
                int GlobalX = MyXOffset + RankX;
                int GlobalY = MyYOffset + RankY;
                // Update this ghost node cell's information with data from other rank
                CellType(CellLocation) = Active;
                GrainID(CellLocation) = (int)(GhostNodesGR[6*i+1]);
                double cx = GhostNodesGR[6*i+2];
                double cy = GhostNodesGR[6*i+3];
                double cz = GhostNodesGR[6*i+4];
                DOCenter(3*CellLocation) = cx;
                DOCenter(3*CellLocation+1) = cy;
                DOCenter(3*CellLocation+2) = cz;
                int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
                DiagonalLength(CellLocation) = GhostNodesGR[6*i+5];
                double xp = GlobalX + 0.5;
                double yp = GlobalY + 0.5;
                double zp = RankZ + 0.5;
                // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
                CritDiagLengthCalc(xp, yp, zp, MyOrientation,RankX,RankY,RankZ,CellLocation,cx,cy,cz,NeighborX,NeighborY,NeighborZ, GrainUnitVector,TriangleIndex,CritDiagonalLength);
            }
        }
    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (id == 0) cout << "End data collection G" << endl;
     //  if ((id == 169)&&(cycle == 11998)) cout << " data collection end G" << endl;
    // Place ghost node data recieved from right and into plane (if needed)
    if (HRCount > 0) {
        for (int i=0; i<HRCount; i++) {
            // GhostNodesFR: Z coord, GrainID, DOCenterX/Y/Z, DiagLength;
            int RankX = 0;
            int RankY = MyYSlices-1;
            int RankZ = (int)(GhostNodesHR[6*i]);
            int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
            if ((cycle == 0)||((CellType(CellLocation) == Liquid)||(CellType(CellLocation) == Delayed)||(CellType(CellLocation) == LiqSol))) {
                int GlobalX = MyXOffset + RankX;
                int GlobalY = MyYOffset + RankY;
                // Update this ghost node cell's information with data from other rank
                CellType(CellLocation) = Active;
                GrainID(CellLocation) = (int)(GhostNodesHR[6*i+1]);
                double cx = GhostNodesHR[6*i+2];
                double cy = GhostNodesHR[6*i+3];
                double cz = GhostNodesHR[6*i+4];
                DOCenter(3*CellLocation) = cx;
                DOCenter(3*CellLocation+1) = cy;
                DOCenter(3*CellLocation+2) = cz;
                int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
                DiagonalLength(CellLocation) = GhostNodesHR[6*i+5];
                double xp = GlobalX + 0.5;
                double yp = GlobalY + 0.5;
                double zp = RankZ + 0.5;
                // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
                CritDiagLengthCalc(xp, yp, zp, MyOrientation,RankX,RankY,RankZ,CellLocation,cx,cy,cz,NeighborX,NeighborY,NeighborZ, GrainUnitVector,TriangleIndex,CritDiagonalLength);
            }
        }
    }
//        MPI_Barrier(MPI_COMM_WORLD);
//        cout << " Done w ghost nodes " << endl;
//        for (int k=0; k<nz; k++)   {
//            for(int i=0; i<MyXSlices; i++)  {
//                for(int j=0; j<MyYSlices; j++)  {
//                    if ((CellType[k][i][j] == Ghost1)||(CellType[k][i][j] == Ghost2)||(CellType[k][i][j] == Ghost3))
//                        cout << "CellType = " << CellType[k][i][j] << endl;
//                }
//            }
//        }
    
}
//*****************************************************************************/

// 1D domain decomposition: update ghost nodes with new cell data from Nucleation and CellCapture routines
void GhostNodes1D(int cycle, int id, int MyLeft, int MyRight, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nz, int NeighborX[26], int NeighborY[26], int NeighborZ[26], ViewI::HostMirror CellType, ViewF::HostMirror DOCenter, ViewI::HostMirror GrainID, float* GrainUnitVector, ViewI::HostMirror TriangleIndex, int* GrainOrientation, ViewF::HostMirror DiagonalLength, ViewF::HostMirror CritDiagonalLength, int NGrainOrientations) {
    
        // Determine whether or not ghost node information transfer needs to take place
    int ACount = 0;
    int BCount = 0;
    int ARCount, BRCount;
    
    int LocalDomainSize = nz*MyXSlices*MyYSlices;
    for (int D3D1ConvPosition=0; D3D1ConvPosition<LocalDomainSize; D3D1ConvPosition++) {
        if (CellType(D3D1ConvPosition) == Ghost1) {
            // Find Y coordinate of this cell and add to appropriate counter
            int RankZ = floor(D3D1ConvPosition/(MyXSlices*MyYSlices));
            int Rem = D3D1ConvPosition % (MyXSlices*MyYSlices);
            int RankY = Rem % MyYSlices;
            if (RankY == 1) {
                ACount++;
            }
            else if (RankY == MyYSlices-2) {
                BCount++;
            }
        }
    }
    
        // Send BCount, Recieve ARCount (send to the right, recieve on the left)
        MPI_Sendrecv(&BCount,1,MPI_INT,MyRight,0,&ARCount,1,MPI_INT,MyLeft,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        // Send ACount, Recieve BRCount (send to the left, recieve on the right)
        MPI_Sendrecv(&ACount,1,MPI_INT,MyLeft,0,&BRCount,1,MPI_INT,MyRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    
        if (MyLeft == MPI_PROC_NULL) ARCount = 0;
        if (MyRight == MPI_PROC_NULL) BRCount = 0;
    
        // Buffers for recieving ghost node data
        float GhostNodesAR[7*ARCount];
        float GhostNodesBR[7*BRCount];
    
    
        // Collect ghost node data and send to other ranks
        if (ACount > 0) {
            float GhostNodesA[7*ACount];
            int BufACount = 0;
            for (int RankZ=1; RankZ<nz-1; RankZ++) {
                for (int RankX=1; RankX<MyXSlices-1; RankX++) {
                    int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + 1;
                    if (CellType(CellLocation) == Ghost1) {
                        // This cell just became active, store its information to send to id = id - 1
                        GhostNodesA[7*BufACount] = (float)(RankX);
                        // cout << "A " << BufACount << " " << GhostNodesA[7*BufACount] << endl;
                        GhostNodesA[7*BufACount+1] = (float)(RankZ);
                        GhostNodesA[7*BufACount+2] = (float)(GrainID(CellLocation));
                        GhostNodesA[7*BufACount+3] = DOCenter(3*CellLocation);
                        GhostNodesA[7*BufACount+4] = DOCenter(3*CellLocation+1);
                        GhostNodesA[7*BufACount+5] = DOCenter(3*CellLocation+2);
                        GhostNodesA[7*BufACount+6] = DiagonalLength(CellLocation);
                        CellType(CellLocation) = Active;
                        BufACount++;
                        //cout << "ID = " << id << " Q A at " << RankX << " " << RankZ << endl;
                    }
                }
            }
            //cout << "ID " << id << " has ACount = " << ACount << " and is sending " << BufACount << endl;
            if (BRCount == 0) {
                // Sending data to id = id - 1 only
                MPI_Send(&GhostNodesA,ACount*7,MPI_FLOAT,MyLeft,0,MPI_COMM_WORLD);
            }
            else {
                // Sending data to id = id - 1 and recieving data from id = id + 1
                MPI_Sendrecv(&GhostNodesA,ACount*7,MPI_FLOAT,MyLeft,0,&GhostNodesBR,BRCount*7,MPI_FLOAT,MyRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
    
        }
        else if (BRCount > 0) {
            // Recieving data from id = id + 1 only
            MPI_Recv(&GhostNodesBR,BRCount*7,MPI_FLOAT,MyRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    
        if (BCount > 0) {
            float GhostNodesB[7*BCount];
            int BufBCount = 0;
            for (int RankZ=1; RankZ<nz-1; RankZ++) {
                for (int RankX=1; RankX<MyXSlices-1; RankX++) {
                    int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + MyYSlices-2;
                    if (CellType(CellLocation) == Ghost1) {
    
                        // This cell just became active, store its information to send to id = id - 1
                        GhostNodesB[7*BufBCount] = (float)(RankX);
                        GhostNodesB[7*BufBCount+1] = (float)(RankZ);
                        GhostNodesB[7*BufBCount+2] = (float)(GrainID(CellLocation));
                        GhostNodesB[7*BufBCount+3] = DOCenter(3*CellLocation);
                        GhostNodesB[7*BufBCount+4] = DOCenter(3*CellLocation+1);
                        GhostNodesB[7*BufBCount+5] = DOCenter(3*CellLocation+2);
                        GhostNodesB[7*BufBCount+6] = DiagonalLength(CellLocation);
                        CellType(CellLocation) = Active;
                        BufBCount++;
                        //cout << "ID = " << id << " Q B at " << RankX << " " << RankZ << endl;
                    }
                }
            }
            //cout << "ID " << id << " has BCount = " << BCount << " and is sending " << BufBCount << endl;
            if (ARCount == 0) {
                // Sending data to id = id + 1 only
                MPI_Send(&GhostNodesB,BCount*7,MPI_FLOAT,MyRight,1,MPI_COMM_WORLD);
            }
            else {
                // Sending data to id = id + 1 and recieving data from id = id - 1
                MPI_Sendrecv(&GhostNodesB,BCount*7,MPI_FLOAT,MyRight,1,&GhostNodesAR,ARCount*7,MPI_FLOAT,MyLeft,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        else if (ARCount > 0) {
            // Recieving data from id = id - 1 only
            MPI_Recv(&GhostNodesAR,ARCount*7,MPI_FLOAT,MyLeft,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    
        //   MPI_Barrier(MPI_COMM_WORLD);
        //   if (id == 0) cout << "Cycle = " << cycle << endl;
        // if (id == 2) cout << "Rank 2 recieved " << ARCount << " elements from rank 1, starting with " << GhostNodesAR[0] << endl;
        // if (id == 2) cout << "Rank 2 recieved " << BRCount << " elements from rank 3, starting with " << GhostNodesBR[0] << endl;
        //if (id == 0) cout << "sending and recieving done " << cycle << endl;
        // Place ghost node data recieved from the left (if needed)
        if (ARCount > 0) {
            //cout << "id = " << id << " ARC = " << ARCount << endl;
            for (int i=0; i<ARCount; i++) {
                // GhostNodesAR: X coord, Z coord, GrainID, DOCenterX/Y/Z, DiagLength;
                int RankX = (int)(GhostNodesAR[7*i]);
                int RankY = 0;
                int RankZ = (int)(GhostNodesAR[7*i+1]);
                int CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
                if ((cycle == 0)||((CellType(CellLocation) == Liquid)||(CellType(CellLocation) == Delayed)||(CellType(CellLocation) == LiqSol))) {
                    int GlobalX = RankX + MyXOffset;
                    int GlobalY = MyYOffset;
                    // Update this ghost node cell's information with data from other rank
                    CellType(CellLocation) = Active;
                    GrainID(CellLocation) = (int)(GhostNodesAR[7*i+2]);
                    double cx = GhostNodesAR[7*i+3];
                    double cy = GhostNodesAR[7*i+4];
                    double cz = GhostNodesAR[7*i+5];
                    DOCenter(3*CellLocation) = cx;
                    DOCenter(3*CellLocation+1) = cy;
                    DOCenter(3*CellLocation+2) = cz;
                    int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
                    DiagonalLength(CellLocation) = GhostNodesAR[7*i+6];
                    double xp = GlobalX + 0.5;
                    double yp = GlobalY + 0.5;
                    double zp = RankZ + 0.5;
                    // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
                    CritDiagLengthCalc(xp, yp, zp, MyOrientation,RankX,RankY,RankZ,CellLocation,cx,cy,cz,NeighborX,NeighborY,NeighborZ, GrainUnitVector,TriangleIndex,CritDiagonalLength);
                }
            }
        }
    
        // Place ghost node data recieved from the right (if needed)
        if (BRCount > 0) {
            //cout << "id = " << id << " BRC = " << BRCount << endl;
            for (int i=0; i<BRCount; i++) {
                // GhostNodesAR: X coord, Z coord, GrainID, DOCenterX/Y/Z, DiagLength;
                int RankX = (int)(GhostNodesBR[7*i]);
                int RankY = MyYSlices-1;
                int RankZ = (int)(GhostNodesBR[7*i+1]);
                int CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
                if ((cycle == 0)||((CellType(CellLocation) == Liquid)||(CellType(CellLocation) == Delayed)||(CellType(CellLocation) == LiqSol))) {
                    //if (cycle == 1003) cout << "Cell at " << RankX << " " << RankY << " " <<RankZ << " now active" << endl;
                    int GlobalX = MyXOffset + RankX;
                    int GlobalY = MyYOffset + RankY;
                    // Update this ghost node cell's information with data from other rank
                    CellType(CellLocation) = Active;
                    GrainID(CellLocation) = (int)(GhostNodesBR[7*i+2]);
                    double cx = GhostNodesBR[7*i+3];
                    double cy = GhostNodesBR[7*i+4];
                    double cz = GhostNodesBR[7*i+5];
                    //if (id == 2) cout << "RECIEVING CYCLE " << cycle << " Center " << cx << " " << cy << " " << cz << endl;
                    DOCenter(3*CellLocation) = cx;
                    DOCenter(3*CellLocation+1) = cy;
                    DOCenter(3*CellLocation+2) = cz;
                    int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
                    DiagonalLength(CellLocation) = GhostNodesBR[7*i+6];
                    double xp = GlobalX + 0.5;
                    double yp = GlobalY + 0.5;
                    double zp = RankZ + 0.5;
                    // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
                    CritDiagLengthCalc(xp, yp, zp, MyOrientation,RankX,RankY,RankZ,CellLocation,cx,cy,cz,NeighborX,NeighborY,NeighborZ, GrainUnitVector,TriangleIndex,CritDiagonalLength);
                    //                for (int l=0; l<26; l++) {
                    //                    cout << "l = " << l << " CDL = " << CritDiagonalLength[BoxZ][RankX][RankY][l] << endl;
                    //                }
                }
                //            else if (CellType[RankZ][RankX][RankY] == 'Q') {
                //                cout << "ID " << id << " competition for cell " << RankX << " " << RankY << " " << RankZ << endl;
                //                //CellType[RankZ][RankX][RankY] = Active;
                //            }
            }
        }
        //  MPI_Barrier(MPI_COMM_WORLD);
        //  if (id == 0) cout << "Cycle = " << cycle << endl;
    
}
