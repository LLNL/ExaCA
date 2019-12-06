#include "header.h"
using namespace std;

// 2D domain decomposition: update ghost nodes with new cell data from Nucleation and CellCapture routines
void GhostNodes2D_GPU(int cycle, int id, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nz, int NeighborX[26], int NeighborY[26], int NeighborZ[26], ViewI CellType, ViewF DOCenter, ViewI GrainID, float* GrainUnitVector, ViewI TriangleIndex, int* GrainOrientation, ViewF DiagonalLength, ViewF CritDiagonalLength, int NGrainOrientations, ViewBufCounts ACount, ViewBufCounts BCount, ViewBufCounts CCount, ViewBufCounts DCount, ViewBufCounts ECount, ViewBufCounts FCount, ViewBufCounts GCount, ViewBufCounts HCount) {

    Kokkos::fence();
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Determine whether or not ghost node information transfer needs to take place
    ViewBufCounts ARCount("ARCount",1);
    ViewBufCounts BRCount("BRCount",1);
    ViewBufCounts CRCount("CRCount",1);
    ViewBufCounts DRCount("DRCount",1);
    ViewBufCounts ERCount("ERCount",1);
    ViewBufCounts FRCount("FRCount",1);
    ViewBufCounts GRCount("GRCount",1);
    ViewBufCounts HRCount("HRCount",1);
    
    // Send BCount, Recieve ARCount (send to the right, recieve on the left)
    MPI_Sendrecv(BCount.data(),1,MPI_INT,MyRight,0,ARCount.data(),1,MPI_INT,MyLeft,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Send ACount, Recieve BRCount (send to the left, recieve on the right)
    MPI_Sendrecv(ACount.data(),1,MPI_INT,MyLeft,0,BRCount.data(),1,MPI_INT,MyRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Send CCount, Recieve DRCount (send into the plane, recieve out of the plane)
    MPI_Sendrecv(CCount.data(),1,MPI_INT,MyIn,0,DRCount.data(),1,MPI_INT,MyOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Send DCount, Recieve CRCount (send out of the plane, recieve into the plane)
    MPI_Sendrecv(DCount.data(),1,MPI_INT,MyOut,0,CRCount.data(),1,MPI_INT,MyIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Send HCount, Recieve ERCount
    MPI_Sendrecv(HCount.data(),1,MPI_INT,MyRightOut,0,ERCount.data(),1,MPI_INT,MyLeftIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Send ECount, Recieve HRCount
    MPI_Sendrecv(ECount.data(),1,MPI_INT,MyLeftIn,0,HRCount.data(),1,MPI_INT,MyRightOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Send GCount, Recieve FRCount
    MPI_Sendrecv(GCount.data(),1,MPI_INT,MyLeftOut,0,FRCount.data(),1,MPI_INT,MyRightIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Send FCount, Recieve GRCount
    MPI_Sendrecv(FCount.data(),1,MPI_INT,MyRightIn,0,GRCount.data(),1,MPI_INT,MyLeftOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Buffers for recieving ghost node data - as views
    ViewBuf GhostNodesAR("GhostNodesAR",7*ARCount(0));
    ViewBuf GhostNodesBR("GhostNodesBR",7*BRCount(0));
    ViewBuf GhostNodesCR("GhostNodesCR",7*CRCount(0));
    ViewBuf GhostNodesDR("GhostNodesDR",7*DRCount(0));
    ViewBuf GhostNodesER("GhostNodesER",7*ERCount(0));
    ViewBuf GhostNodesFR("GhostNodesFR",7*FRCount(0));
    ViewBuf GhostNodesGR("GhostNodesGR",7*GRCount(0));
    ViewBuf GhostNodesHR("GhostNodesHR",7*HRCount(0));

    // Collect ghost node data and send to other ranks- left and right
    ViewBuf GhostNodesA("GhostNodesA",7*ACount(0));
    Kokkos::parallel_for("BufALoop",ACount(0), KOKKOS_LAMBDA (const int& LocalCount) {
        int Flag = 0;
        for (int RankZ=1; RankZ<nz-1; RankZ++) {
            for (int RankX=1; RankX<MyXSlices-1; RankX++) {
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + 1;
                int GhostCheck3 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost3,Ghost2);
                int GhostCheck2 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost2,Ghost1);
                int GhostCheck1 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost1,Active);
                if ((GhostCheck3 == Ghost3)||(GhostCheck2 == Ghost2)||(GhostCheck1 == Ghost1)) {
                    GhostNodesA(7*LocalCount) = (float)(RankX);
                    GhostNodesA(7*LocalCount+1) = (float)(RankZ);
                    GhostNodesA(7*LocalCount+2) = (float)(GrainID(CellLocation));
                    GhostNodesA(7*LocalCount+3) = DOCenter(3*CellLocation);
                    GhostNodesA(7*LocalCount+4) = DOCenter(3*CellLocation+1);
                    GhostNodesA(7*LocalCount+5) = DOCenter(3*CellLocation+2);
                    GhostNodesA(7*LocalCount+6) = DiagonalLength(CellLocation);
                    Flag = 1;
                    break;
                }
            }
            if (Flag == 1) break;
        }
    });
    ViewBuf GhostNodesB("GhostNodesB",7*BCount(0));
    Kokkos::parallel_for("BufBLoop",BCount(0), KOKKOS_LAMBDA (const int& LocalCount) {
        int Flag = 0;
        for (int RankZ=1; RankZ<nz-1; RankZ++) {
            for (int RankX=1; RankX<MyXSlices-1; RankX++) {
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + MyYSlices-2;
                int GhostCheck3 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost3,Ghost2);
                int GhostCheck2 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost2,Ghost1);
                int GhostCheck1 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost1,Active);
                if ((GhostCheck3 == Ghost3)||(GhostCheck2 == Ghost2)||(GhostCheck1 == Ghost1)) {
                    GhostNodesB(7*LocalCount) = (float)(RankX);
                    GhostNodesB(7*LocalCount+1) = (float)(RankZ);
                    GhostNodesB(7*LocalCount+2) = (float)(GrainID(CellLocation));
                    GhostNodesB(7*LocalCount+3) = DOCenter(3*CellLocation);
                    GhostNodesB(7*LocalCount+4) = DOCenter(3*CellLocation+1);
                    GhostNodesB(7*LocalCount+5) = DOCenter(3*CellLocation+2);
                    GhostNodesB(7*LocalCount+6) = DiagonalLength(CellLocation);
                    Flag = 1;
                    break;
                }
            }
            if (Flag == 1) break;
        }
    });
    // Collect ghost node data to send to other ranks - in and out of plane
    ViewBuf GhostNodesC("GhostNodesC",7*CCount(0));
    Kokkos::parallel_for("BufCLoop",CCount(0), KOKKOS_LAMBDA (const int& LocalCount) {
        int Flag = 0;
        for (int RankZ=1; RankZ<nz-1; RankZ++) {
            for (int RankY=1; RankY<MyYSlices-1; RankY++) {
                int CellLocation = RankZ*MyXSlices*MyYSlices +(MyXSlices-2)*MyYSlices + RankY;
                int GhostCheck3 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost3,Ghost2);
                int GhostCheck2 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost2,Ghost1);
                int GhostCheck1 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost1,Active);
                if ((GhostCheck3 == Ghost3)||(GhostCheck2 == Ghost2)||(GhostCheck1 == Ghost1)) {
                    GhostNodesC(7*LocalCount) = (float)(RankY);
                    GhostNodesC(7*LocalCount+1) = (float)(RankZ);
                    GhostNodesC(7*LocalCount+2) = (float)(GrainID(CellLocation));
                    GhostNodesC(7*LocalCount+3) = DOCenter(3*CellLocation);
                    GhostNodesC(7*LocalCount+4) = DOCenter(3*CellLocation+1);
                    GhostNodesC(7*LocalCount+5) = DOCenter(3*CellLocation+2);
                    GhostNodesC(7*LocalCount+6) = DiagonalLength(CellLocation);
                    Flag = 1;
                    break;
                }
            }
            if (Flag == 1) break;
        }
    });
    ViewBuf GhostNodesD("GhostNodesD",7*DCount(0));
    Kokkos::parallel_for("BufDLoop",DCount(0), KOKKOS_LAMBDA (const int& LocalCount) {
        int Flag = 0;
        for (int RankZ=1; RankZ<nz-1; RankZ++) {
            for (int RankY=1; RankY<MyYSlices-1; RankY++) {
                int CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices + RankY;
                int GhostCheck3 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost3,Ghost2);
                int GhostCheck2 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost2,Ghost1);
                int GhostCheck1 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost1,Active);
                if ((GhostCheck3 == Ghost3)||(GhostCheck2 == Ghost2)||(GhostCheck1 == Ghost1)) {
                    GhostNodesD(7*LocalCount) = (float)(RankY);
                    GhostNodesD(7*LocalCount+1) = (float)(RankZ);
                    GhostNodesD(7*LocalCount+2) = (float)(GrainID(CellLocation));
                    GhostNodesD(7*LocalCount+3) = DOCenter(3*CellLocation);
                    GhostNodesD(7*LocalCount+4) = DOCenter(3*CellLocation+1);
                    GhostNodesD(7*LocalCount+5) = DOCenter(3*CellLocation+2);
                    GhostNodesD(7*LocalCount+6) = DiagonalLength(CellLocation);
                    Flag = 1;
                    break;
                }
            }
            if (Flag == 1) break;
        }
    });
    // Collect ghost node data to send to other ranks - corners
    ViewBuf GhostNodesE("GhostNodesE",6*ECount(0));
    Kokkos::parallel_for("BufELoop",ECount(0), KOKKOS_LAMBDA (const int& LocalCount) {
        for (int RankZ=1; RankZ<nz-1; RankZ++) {
            int CellLocation = RankZ*MyXSlices*MyYSlices + (MyXSlices-2)*(MyYSlices) + 1;
            int GhostCheck3 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost3,Ghost2);
            int GhostCheck2 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost2,Ghost1);
            int GhostCheck1 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost1,Active);
            if ((GhostCheck3 == Ghost3)||(GhostCheck2 == Ghost2)||(GhostCheck1 == Ghost1)) {
                GhostNodesE(6*LocalCount) = (float)(RankZ);
                GhostNodesE(6*LocalCount+1) = (float)(GrainID(CellLocation));
                GhostNodesE(6*LocalCount+2) = DOCenter(3*CellLocation);
                GhostNodesE(6*LocalCount+3) = DOCenter(3*CellLocation+1);
                GhostNodesE(6*LocalCount+4) = DOCenter(3*CellLocation+2);
                GhostNodesE(6*LocalCount+5) = DiagonalLength(CellLocation);
                break;
            }
        }
    });
    ViewBuf GhostNodesF("GhostNodesF",6*FCount(0));
    Kokkos::parallel_for("BufFLoop",FCount(0), KOKKOS_LAMBDA (const int& LocalCount) {
        for (int RankZ=1; RankZ<nz-1; RankZ++) {
            int CellLocation = RankZ*MyXSlices*MyYSlices + (MyXSlices-2)*MyYSlices + MyYSlices-2;
            int GhostCheck3 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost3,Ghost2);
            int GhostCheck2 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost2,Ghost1);
            int GhostCheck1 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost1,Active);
            if ((GhostCheck3 == Ghost3)||(GhostCheck2 == Ghost2)||(GhostCheck1 == Ghost1)) {
                GhostNodesF(6*LocalCount) = (float)(RankZ);
                GhostNodesF(6*LocalCount+1) = (float)(GrainID(CellLocation));
                GhostNodesF(6*LocalCount+2) = DOCenter(3*CellLocation);
                GhostNodesF(6*LocalCount+3) = DOCenter(3*CellLocation+1);
                GhostNodesF(6*LocalCount+4) = DOCenter(3*CellLocation+2);
                GhostNodesF(6*LocalCount+5) = DiagonalLength(CellLocation);
                break;
            }
        }
    });
    ViewBuf GhostNodesG("GhostNodesG",6*GCount(0));
    Kokkos::parallel_for("BufGLoop",GCount(0), KOKKOS_LAMBDA (const int& LocalCount) {
        for (int RankZ=1; RankZ<nz-1; RankZ++) {
            int CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices + 1;
            int GhostCheck3 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost3,Ghost2);
            int GhostCheck2 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost2,Ghost1);
            int GhostCheck1 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost1,Active);
            if ((GhostCheck3 == Ghost3)||(GhostCheck2 == Ghost2)||(GhostCheck1 == Ghost1)) {
                GhostNodesG(6*LocalCount) = (float)(RankZ);
                GhostNodesG(6*LocalCount+1) = (float)(GrainID(CellLocation));
                GhostNodesG(6*LocalCount+2) = DOCenter(3*CellLocation);
                GhostNodesG(6*LocalCount+3) = DOCenter(3*CellLocation+1);
                GhostNodesG(6*LocalCount+4) = DOCenter(3*CellLocation+2);
                GhostNodesG(6*LocalCount+5) = DiagonalLength(CellLocation);
                break;
            }
        }
    });
    ViewBuf GhostNodesH("GhostNodesH",6*HCount(0));
    Kokkos::parallel_for("BufHLoop",HCount(0), KOKKOS_LAMBDA (const int& LocalCount) {
        for (int RankZ=1; RankZ<nz-1; RankZ++) {
            int CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices + MyYSlices-2;
            int GhostCheck3 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost3,Ghost2);
            int GhostCheck2 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost2,Ghost1);
            int GhostCheck1 = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost1,Active);
            if ((GhostCheck3 == Ghost3)||(GhostCheck2 == Ghost2)||(GhostCheck1 == Ghost1)) {
                GhostNodesH(6*LocalCount) = (float)(RankZ);
                GhostNodesH(6*LocalCount+1) = (float)(GrainID(CellLocation));
                GhostNodesH(6*LocalCount+2) = DOCenter(3*CellLocation);
                GhostNodesH(6*LocalCount+3) = DOCenter(3*CellLocation+1);
                GhostNodesH(6*LocalCount+4) = DOCenter(3*CellLocation+2);
                GhostNodesH(6*LocalCount+5) = DiagonalLength(CellLocation);
                break;
            }
        }
    });

   // Sending/receiving data to left/right
   MPI_Sendrecv(GhostNodesA.data(),7*ACount(0),MPI_FLOAT,MyLeft,0,GhostNodesBR.data(),7*BRCount(0),MPI_FLOAT,MyRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
   MPI_Sendrecv(GhostNodesB.data(),7*BCount(0),MPI_FLOAT,MyRight,1,GhostNodesAR.data(),7*ARCount(0),MPI_FLOAT,MyLeft,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

   // Sending data and recieving data in/out of plane
   MPI_Sendrecv(GhostNodesC.data(),7*CCount(0),MPI_FLOAT,MyIn,0,GhostNodesDR.data(),7*DRCount(0),MPI_FLOAT,MyOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
   MPI_Sendrecv(GhostNodesD.data(),7*DCount(0),MPI_FLOAT,MyOut,1,GhostNodesCR.data(),7*CRCount(0),MPI_FLOAT,MyIn,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

   // Sending/receiving first pair of corners
   MPI_Sendrecv(GhostNodesE.data(),6*ECount(0),MPI_FLOAT,MyLeftIn,0,GhostNodesHR.data(),6*HRCount(0),MPI_FLOAT,MyRightOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
   MPI_Sendrecv(GhostNodesH.data(),6*HCount(0),MPI_FLOAT,MyRightOut,0,GhostNodesER.data(),6*ERCount(0),MPI_FLOAT,MyLeftIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

   // Second pair of corners
   MPI_Sendrecv(GhostNodesF.data(),6*FCount(0),MPI_FLOAT,MyRightIn,1,GhostNodesGR.data(),6*GRCount(0),MPI_FLOAT,MyLeftOut,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
   MPI_Sendrecv(GhostNodesG.data(),6*GCount(0),MPI_FLOAT,MyLeftOut,1,GhostNodesFR.data(),6*FRCount(0),MPI_FLOAT,MyRightIn,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    
    // Unpack buffers
    Kokkos::parallel_for("UpdateALoop",ARCount(0), KOKKOS_LAMBDA (const int& i) {
        int RankX = (int)(GhostNodesAR(7*i));
        int RankY = 0;
        int RankZ = (int)(GhostNodesAR(7*i+1));
        int CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
        //printf("CellLocation %d CellType %d \n",CellLocation,CellType(CellLocation));
        if ((cycle == 0)||((CellType(CellLocation) == Liquid)||(CellType(CellLocation) == Delayed)||(CellType(CellLocation) == LiqSol))) {
            //int GlobalX = RankX + MyXOffset;
            //int GlobalY = MyYOffset;
            // Update this ghost node cell's information with data from other rank
            CellType(CellLocation) = Active;
            GrainID(CellLocation) = (int)(GhostNodesAR(7*i+2));
            double cx = GhostNodesAR(7*i+3);
            double cy = GhostNodesAR(7*i+4);
            double cz = GhostNodesAR(7*i+5);
            DOCenter(3*CellLocation) = cx;
            DOCenter(3*CellLocation+1) = cy;
            DOCenter(3*CellLocation+2) = cz;
            int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
            DiagonalLength(CellLocation) = GhostNodesAR(7*i+6);
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
                double x0 = xp + NeighborX[n] - cx;
                double y0 = yp + NeighborY[n] - cy;
                double z0 = zp + NeighborZ[n] - cz;
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
        }
    });
    
    Kokkos::parallel_for("UpdateBLoop",BRCount(0), KOKKOS_LAMBDA (const int& i) {
        // GhostNodesBR: X coord, Z coord, GrainID, DOCenterX/Y/Z, DiagLength;
        int RankX = (int)(GhostNodesBR(7*i));
        int RankY = MyYSlices-1;
        int RankZ = (int)(GhostNodesBR(7*i+1));
        int CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
        if ((cycle == 0)||((CellType(CellLocation) == Liquid)||(CellType(CellLocation) == Delayed)||(CellType(CellLocation) == LiqSol))) {
            //if (cycle == 1003) cout << "Cell at " << RankX << " " << RankY << " " <<RankZ << " now active" << endl;
            //int GlobalX = MyXOffset + RankX;
            //int GlobalY = MyYOffset + RankY;
            // Update this ghost node cell's information with data from other rank
            CellType(CellLocation) = Active;
            GrainID(CellLocation) = (int)(GhostNodesBR(7*i+2));
            double cx = GhostNodesBR(7*i+3);
            double cy = GhostNodesBR(7*i+4);
            double cz = GhostNodesBR(7*i+5);
            //if (id == 2) cout << "RECIEVING CYCLE " << cycle << " Center " << cx << " " << cy << " " << cz << endl;
            DOCenter(3*CellLocation) = cx;
            DOCenter(3*CellLocation+1) = cy;
            DOCenter(3*CellLocation+2) = cz;
            int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
            DiagonalLength(CellLocation) = GhostNodesBR(7*i+6);
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
                double x0 = xp + NeighborX[n] - cx;
                double y0 = yp + NeighborY[n] - cy;
                double z0 = zp + NeighborZ[n] - cz;
                // mag0 is the magnitude of (x0,y0,z0)
                
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
                //                                if (CDLVal == 0.0) printf("Zero CDLVal : %d %d %d %d %f %d %d %d %f %f %f",MyNeighborX,MyNeighborY,MyNeighborZ,n,mag0,index1,index2,index3,normx,normy,normz);
            }
            //                for (int l=0; l<26; l++) {
            //                    cout << "l = " << l << " CDL = " << CritDiagonalLength[BoxZ][RankX][RankY][l] << endl;
            //                }
        }
        //            else if (CellType[RankZ][RankX][RankY] == 'Q') {
        //                cout << "ID " << id << " competition for cell " << RankX << " " << RankY << " " << RankZ << endl;
        //                //CellType[RankZ][RankX][RankY] = Active;
        //            }
    });
    
    Kokkos::parallel_for("UpdateCLoop",CRCount(0), KOKKOS_LAMBDA (const int& i) {
        // GhostNodesCR: Y coord, Z coord, GrainID, DOCenterX/Y/Z, DiagLength;
        int RankX = MyXSlices-2;
        int RankY = (int)(GhostNodesCR(7*i));
        int RankZ = (int)(GhostNodesCR(7*i+1));
        int CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
        if ((cycle == 0)||((CellType(CellLocation) == Liquid)||(CellType(CellLocation) == Delayed)||(CellType(CellLocation) == LiqSol))) {
            //if (cycle == 1003) cout << "Cell at " << RankX << " " << RankY << " " <<RankZ << " now active" << endl;
            //int GlobalX = MyXOffset + RankX;
            //int GlobalY = MyYOffset + RankY;
            // Update this ghost node cell's information with data from other rank
            CellType(CellLocation) = Active;
            GrainID(CellLocation) = (int)(GhostNodesCR(7*i+2));
            double cx = GhostNodesCR(7*i+3);
            double cy = GhostNodesCR(7*i+4);
            double cz = GhostNodesCR(7*i+5);
            //if (id == 2) cout << "RECIEVING CYCLE " << cycle << " Center " << cx << " " << cy << " " << cz << endl;
            DOCenter(3*CellLocation) = cx;
            DOCenter(3*CellLocation+1) = cy;
            DOCenter(3*CellLocation+2) = cz;
            int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
            DiagonalLength(CellLocation) = GhostNodesCR(7*i+6);
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
                double x0 = xp + NeighborX[n] - cx;
                double y0 = yp + NeighborY[n] - cy;
                double z0 = zp + NeighborZ[n] - cz;
                // mag0 is the magnitude of (x0,y0,z0)
                
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
                //                                if (CDLVal == 0.0) printf("Zero CDLVal : %d %d %d %d %f %d %d %d %f %f %f",MyNeighborX,MyNeighborY,MyNeighborZ,n,mag0,index1,index2,index3,normx,normy,normz);
            }
            //                for (int l=0; l<26; l++) {
            //                    cout << "l = " << l << " CDL = " << CritDiagonalLength[BoxZ][RankX][RankY][l] << endl;
            //                }
        }
        //            else if (CellType[RankZ][RankX][RankY] == 'Q') {
        //                cout << "ID " << id << " competition for cell " << RankX << " " << RankY << " " << RankZ << endl;
        //                //CellType[RankZ][RankX][RankY] = Active;
        //            }
    });
    
    Kokkos::parallel_for("UpdateDLoop",DRCount(0), KOKKOS_LAMBDA (const int& i) {
        // GhostNodesDR: Y coord, Z coord, GrainID, DOCenterX/Y/Z, DiagLength;
        int RankX = 1;
        int RankY = (int)(GhostNodesDR(7*i));
        int RankZ = (int)(GhostNodesDR(7*i+1));
        int CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
        if ((cycle == 0)||((CellType(CellLocation) == Liquid)||(CellType(CellLocation) == Delayed)||(CellType(CellLocation) == LiqSol))) {
            //if (cycle == 1003) cout << "Cell at " << RankX << " " << RankY << " " <<RankZ << " now active" << endl;
            //int GlobalX = MyXOffset + RankX;
            //int GlobalY = MyYOffset + RankY;
            // Update this ghost node cell's information with data from other rank
            CellType(CellLocation) = Active;
            GrainID(CellLocation) = (int)(GhostNodesDR(7*i+2));
            double cx = GhostNodesDR(7*i+3);
            double cy = GhostNodesDR(7*i+4);
            double cz = GhostNodesDR(7*i+5);
            //if (id == 2) cout << "RECIEVING CYCLE " << cycle << " Center " << cx << " " << cy << " " << cz << endl;
            DOCenter(3*CellLocation) = cx;
            DOCenter(3*CellLocation+1) = cy;
            DOCenter(3*CellLocation+2) = cz;
            int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
            DiagonalLength(CellLocation) = GhostNodesDR(7*i+6);
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
                double x0 = xp + NeighborX[n] - cx;
                double y0 = yp + NeighborY[n] - cy;
                double z0 = zp + NeighborZ[n] - cz;
                // mag0 is the magnitude of (x0,y0,z0)
                
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
                //                                if (CDLVal == 0.0) printf("Zero CDLVal : %d %d %d %d %f %d %d %d %f %f %f",MyNeighborX,MyNeighborY,MyNeighborZ,n,mag0,index1,index2,index3,normx,normy,normz);
            }
            //                for (int l=0; l<26; l++) {
            //                    cout << "l = " << l << " CDL = " << CritDiagonalLength[BoxZ][RankX][RankY][l] << endl;
            //                }
        }
        //            else if (CellType[RankZ][RankX][RankY] == 'Q') {
        //                cout << "ID " << id << " competition for cell " << RankX << " " << RankY << " " << RankZ << endl;
        //                //CellType[RankZ][RankX][RankY] = Active;
        //            }
    });
    
    Kokkos::parallel_for("UpdateELoop",ERCount(0), KOKKOS_LAMBDA (const int& i) {
        // GhostNodesER: Z coord, GrainID, DOCenterX/Y/Z, DiagLength;
        int RankX = MyXSlices-2;
        int RankY = 1;
        int RankZ = (int)(GhostNodesER(6*i));
        int CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
        if ((cycle == 0)||((CellType(CellLocation) == Liquid)||(CellType(CellLocation) == Delayed)||(CellType(CellLocation) == LiqSol))) {
            //if (cycle == 1003) cout << "Cell at " << RankX << " " << RankY << " " <<RankZ << " now active" << endl;
            //int GlobalX = MyXOffset + RankX;
            //int GlobalY = MyYOffset + RankY;
            // Update this ghost node cell's information with data from other rank
            CellType(CellLocation) = Active;
            GrainID(CellLocation) = (int)(GhostNodesER(6*i+1));
            double cx = GhostNodesER(6*i+2);
            double cy = GhostNodesER(6*i+3);
            double cz = GhostNodesER(6*i+4);
            //if (id == 2) cout << "RECIEVING CYCLE " << cycle << " Center " << cx << " " << cy << " " << cz << endl;
            DOCenter(3*CellLocation) = cx;
            DOCenter(3*CellLocation+1) = cy;
            DOCenter(3*CellLocation+2) = cz;
            int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
            DiagonalLength(CellLocation) = GhostNodesER(6*i+5);
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
                double x0 = xp + NeighborX[n] - cx;
                double y0 = yp + NeighborY[n] - cy;
                double z0 = zp + NeighborZ[n] - cz;
                // mag0 is the magnitude of (x0,y0,z0)
                
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
                //                                if (CDLVal == 0.0) printf("Zero CDLVal : %d %d %d %d %f %d %d %d %f %f %f",MyNeighborX,MyNeighborY,MyNeighborZ,n,mag0,index1,index2,index3,normx,normy,normz);
            }
            //                for (int l=0; l<26; l++) {
            //                    cout << "l = " << l << " CDL = " << CritDiagonalLength[BoxZ][RankX][RankY][l] << endl;
            //                }
        }
        //            else if (CellType[RankZ][RankX][RankY] == 'Q') {
        //                cout << "ID " << id << " competition for cell " << RankX << " " << RankY << " " << RankZ << endl;
        //                //CellType[RankZ][RankX][RankY] = Active;
        //            }
    });
    Kokkos::parallel_for("UpdateFLoop",FRCount(0), KOKKOS_LAMBDA (const int& i) {
        // GhostNodesFR: Z coord, GrainID, DOCenterX/Y/Z, DiagLength;
        int RankX = MyXSlices-2;
        int RankY = MyYSlices-2;
        int RankZ = (int)(GhostNodesFR(6*i));
        int CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
        if ((cycle == 0)||((CellType(CellLocation) == Liquid)||(CellType(CellLocation) == Delayed)||(CellType(CellLocation) == LiqSol))) {
            //if (cycle == 1003) cout << "Cell at " << RankX << " " << RankY << " " <<RankZ << " now active" << endl;
            //int GlobalX = MyXOffset + RankX;
            //int GlobalY = MyYOffset + RankY;
            // Update this ghost node cell's information with data from other rank
            CellType(CellLocation) = Active;
            GrainID(CellLocation) = (int)(GhostNodesFR(6*i+1));
            double cx = GhostNodesFR(6*i+2);
            double cy = GhostNodesFR(6*i+3);
            double cz = GhostNodesFR(6*i+4);
            //if (id == 2) cout << "RECIEVING CYCLE " << cycle << " Center " << cx << " " << cy << " " << cz << endl;
            DOCenter(3*CellLocation) = cx;
            DOCenter(3*CellLocation+1) = cy;
            DOCenter(3*CellLocation+2) = cz;
            int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
            DiagonalLength(CellLocation) = GhostNodesFR(6*i+5);
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
                double x0 = xp + NeighborX[n] - cx;
                double y0 = yp + NeighborY[n] - cy;
                double z0 = zp + NeighborZ[n] - cz;
                // mag0 is the magnitude of (x0,y0,z0)
                
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
                //                                if (CDLVal == 0.0) printf("Zero CDLVal : %d %d %d %d %f %d %d %d %f %f %f",MyNeighborX,MyNeighborY,MyNeighborZ,n,mag0,index1,index2,index3,normx,normy,normz);
            }
            //                for (int l=0; l<26; l++) {
            //                    cout << "l = " << l << " CDL = " << CritDiagonalLength[BoxZ][RankX][RankY][l] << endl;
            //                }
        }
        //            else if (CellType[RankZ][RankX][RankY] == 'Q') {
        //                cout << "ID " << id << " competition for cell " << RankX << " " << RankY << " " << RankZ << endl;
        //                //CellType[RankZ][RankX][RankY] = Active;
        //            }
    });
    Kokkos::parallel_for("UpdateGLoop",GRCount(0), KOKKOS_LAMBDA (const int& i) {
        // GhostNodesGR: Z coord, GrainID, DOCenterX/Y/Z, DiagLength;
        int RankX = 1;
        int RankY = 1;
        int RankZ = (int)(GhostNodesGR(6*i));
        int CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
        if ((cycle == 0)||((CellType(CellLocation) == Liquid)||(CellType(CellLocation) == Delayed)||(CellType(CellLocation) == LiqSol))) {
            //if (cycle == 1003) cout << "Cell at " << RankX << " " << RankY << " " <<RankZ << " now active" << endl;
            //int GlobalX = MyXOffset + RankX;
            //int GlobalY = MyYOffset + RankY;
            // Update this ghost node cell's information with data from other rank
            CellType(CellLocation) = Active;
            GrainID(CellLocation) = (int)(GhostNodesGR(6*i+1));
            double cx = GhostNodesGR(6*i+2);
            double cy = GhostNodesGR(6*i+3);
            double cz = GhostNodesGR(6*i+4);
            //if (id == 2) cout << "RECIEVING CYCLE " << cycle << " Center " << cx << " " << cy << " " << cz << endl;
            DOCenter(3*CellLocation) = cx;
            DOCenter(3*CellLocation+1) = cy;
            DOCenter(3*CellLocation+2) = cz;
            int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
            DiagonalLength(CellLocation) = GhostNodesGR(6*i+5);
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
                double x0 = xp + NeighborX[n] - cx;
                double y0 = yp + NeighborY[n] - cy;
                double z0 = zp + NeighborZ[n] - cz;
                // mag0 is the magnitude of (x0,y0,z0)
                
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
                //                                if (CDLVal == 0.0) printf("Zero CDLVal : %d %d %d %d %f %d %d %d %f %f %f",MyNeighborX,MyNeighborY,MyNeighborZ,n,mag0,index1,index2,index3,normx,normy,normz);
            }
        }
    });
    Kokkos::parallel_for("UpdateHLoop",HRCount(0), KOKKOS_LAMBDA (const int& i) {
        // GhostNodesHR: Z coord, GrainID, DOCenterX/Y/Z, DiagLength;
        int RankX = 1;
        int RankY = MyYSlices-2;
        int RankZ = (int)(GhostNodesHR(6*i));
        int CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
        if ((cycle == 0)||((CellType(CellLocation) == Liquid)||(CellType(CellLocation) == Delayed)||(CellType(CellLocation) == LiqSol))) {
            //if (cycle == 1003) cout << "Cell at " << RankX << " " << RankY << " " <<RankZ << " now active" << endl;
            //int GlobalX = MyXOffset + RankX;
            //int GlobalY = MyYOffset + RankY;
            // Update this ghost node cell's information with data from other rank
            CellType(CellLocation) = Active;
            GrainID(CellLocation) = (int)(GhostNodesHR(6*i+1));
            double cx = GhostNodesHR(6*i+2);
            double cy = GhostNodesHR(6*i+3);
            double cz = GhostNodesHR(6*i+4);
            //if (id == 2) cout << "RECIEVING CYCLE " << cycle << " Center " << cx << " " << cy << " " << cz << endl;
            DOCenter(3*CellLocation) = cx;
            DOCenter(3*CellLocation+1) = cy;
            DOCenter(3*CellLocation+2) = cz;
            int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
            DiagonalLength(CellLocation) = GhostNodesHR(6*i+5);
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
                double x0 = xp + NeighborX[n] - cx;
                double y0 = yp + NeighborY[n] - cy;
                double z0 = zp + NeighborZ[n] - cz;
                // mag0 is the magnitude of (x0,y0,z0)
                
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
                //                                if (CDLVal == 0.0) printf("Zero CDLVal : %d %d %d %d %f %d %d %d %f %f %f",MyNeighborX,MyNeighborY,MyNeighborZ,n,mag0,index1,index2,index3,normx,normy,normz);
            }
            //                for (int l=0; l<26; l++) {
            //                    cout << "l = " << l << " CDL = " << CritDiagonalLength[BoxZ][RankX][RankY][l] << endl;
            //                }
        }
        //            else if (CellType[RankZ][RankX][RankY] == 'Q') {
        //                cout << "ID " << id << " competition for cell " << RankX << " " << RankY << " " << RankZ << endl;
        //                //CellType[RankZ][RankX][RankY] = Active;
        //            }
    });
}
//*****************************************************************************/

// 1D domain decomposition: update ghost nodes with new cell data from Nucleation and CellCapture routines
void GhostNodes1D_GPU(int cycle, int id, int MyLeft, int MyRight, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nz, int NeighborX[26], int NeighborY[26], int NeighborZ[26], ViewI CellType, ViewF DOCenter, ViewI GrainID, float* GrainUnitVector, ViewI TriangleIndex, int* GrainOrientation, ViewF DiagonalLength, ViewF CritDiagonalLength, int NGrainOrientations, ViewBufCounts ACount, ViewBufCounts BCount) {
    
    Kokkos::fence();
    MPI_Barrier(MPI_COMM_WORLD);
    //cout << "Start" << id << endl;
    
    
    // Determine whether or not ghost node information transfer needs to take place
    ViewBufCounts ARCount("ARCount",1);
    ViewBufCounts BRCount("BRCount",1);

    // Send BCount, Recieve ARCount (send to the right, recieve on the left)
    MPI_Sendrecv(BCount.data(),1,MPI_INT,MyRight,0,ARCount.data(),1,MPI_INT,MyLeft,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    // Send ACount, Recieve BRCount (send to the left, recieve on the right)
    MPI_Sendrecv(ACount.data(),1,MPI_INT,MyLeft,0,BRCount.data(),1,MPI_INT,MyRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Buffers for recieving ghost node data
    ViewBuf GhostNodesAR("GhostNodesAR",7*ARCount(0));
    ViewBuf GhostNodesBR("GhostNodesBR",7*BRCount(0));

    // Collect ghost node data and send to other ranks
    ViewBuf GhostNodesA("GhostNodesA",7*ACount(0));

    //MPI_Barrier(MPI_COMM_WORLD);
    //cout << "A" << id << endl;


    // Load buffer
    Kokkos::parallel_for("BufALoop",ACount(0), KOKKOS_LAMBDA (const int& LocalCount) {
        int Flag = 0;
        for (int RankZ=1; RankZ<nz-1; RankZ++) {
            for (int RankX=1; RankX<MyXSlices-1; RankX++) {
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + 1;
                int GhostCheck = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost1,Active);
                if (GhostCheck == Ghost1) {
                    // This cell just became active, store its information to send to id = id + 1
                    //if (id == 0) printf("Ghost cell at %d ; Count is %d out of %d total \n",CellLocation,LocalCount,ACount(0));
                    GhostNodesA(7*LocalCount) = (float)(RankX);
                    GhostNodesA(7*LocalCount+1) = (float)(RankZ);
                    GhostNodesA(7*LocalCount+2) = (float)(GrainID(CellLocation));
                    GhostNodesA(7*LocalCount+3) = DOCenter(3*CellLocation);
                    GhostNodesA(7*LocalCount+4) = DOCenter(3*CellLocation+1);
                    GhostNodesA(7*LocalCount+5) = DOCenter(3*CellLocation+2);
                    GhostNodesA(7*LocalCount+6) = DiagonalLength(CellLocation);
                    Flag = 1;
                    break;
                }
            }
            if (Flag == 1) break;
        }
    });
    
    //MPI_Barrier(MPI_COMM_WORLD);
    //cout << "B" << id << endl;
    
    // Collect ghost node data and send to other ranks
    ViewBuf GhostNodesB("GhostNodesB",7*BCount(0));
    // Load buffer
    Kokkos::parallel_for("BufBLoop",BCount(0), KOKKOS_LAMBDA (const int& LocalCount) {
        int Flag = 0;
        for (int RankZ=1; RankZ<nz-1; RankZ++) {
            for (int RankX=1; RankX<MyXSlices-1; RankX++) {
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + MyYSlices-2;
                int GhostCheck = Kokkos::atomic_compare_exchange(&CellType(CellLocation),Ghost1,Active);
                if (GhostCheck == Ghost1) {
                    // This cell just became active, store its information to send to id = id - 1
                    //if (id == 0) printf("Ghost cell at %d ; Count is %d out of %d total \n",CellLocation,LocalCount,BCount(0));
                    //int RankY = MyYSlices-2;
                    GhostNodesB(7*LocalCount) = (float)(RankX);
                    GhostNodesB(7*LocalCount+1) = (float)(RankZ);
                    GhostNodesB(7*LocalCount+2) = (float)(GrainID(CellLocation));
                    GhostNodesB(7*LocalCount+3) = DOCenter(3*CellLocation);
                    GhostNodesB(7*LocalCount+4) = DOCenter(3*CellLocation+1);
                    GhostNodesB(7*LocalCount+5) = DOCenter(3*CellLocation+2);
                    GhostNodesB(7*LocalCount+6) = DiagonalLength(CellLocation);
                    for (int n=0; n<26; n++) {
                        //if (id == 0) printf("SEND CELL at CYCLE %d LOC %d %d %d ORIENTATION %d DL = %f TI = %d %d %d DiagonalN %d CDL Val %f \n",cycle,RankX,RankY,RankZ,GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)],DiagonalLength(CellLocation),TriangleIndex(78*CellLocation + 3*n),TriangleIndex(78*CellLocation + 3*n + 1),TriangleIndex(78*CellLocation + 3*n + 2),n,CritDiagonalLength(26*CellLocation+n));
                    }
                    Flag = 1;
                    break;
                }
            }
            if (Flag == 1) break;
        }
    });
    
    
    Kokkos::fence();
    MPI_Barrier(MPI_COMM_WORLD);
    //cout << "C" << id << endl;
//    Kokkos::parallel_for("BufBLoop",BCount(0), KOKKOS_LAMBDA (const int& i) {
//        if (id == 0) printf("Rank 0 sending to Rank 1 a total of %d data points : points %d through %d are %f %f %f %f %f %f %f \n",7*BCount(0),7*i,7*i+6,GhostNodesB(7*i),GhostNodesB(7*i+1),GhostNodesB(7*i+2),GhostNodesB(7*i+3),GhostNodesB(7*i+4),GhostNodesB(7*i+5),GhostNodesB(7*i+6));
//    });
    
    // Send and recieve data
    MPI_Sendrecv(GhostNodesA.data(),7*ACount(0),MPI_FLOAT,MyLeft,0,GhostNodesBR.data(),7*BRCount(0),MPI_FLOAT,MyRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Sendrecv(GhostNodesB.data(),7*BCount(0),MPI_FLOAT,MyRight,1,GhostNodesAR.data(),7*ARCount(0),MPI_FLOAT,MyLeft,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //MPI_Barrier(MPI_COMM_WORLD);
    //cout << "D" << id << endl;
//    Kokkos::parallel_for("UpdateALoop",ARCount(0), KOKKOS_LAMBDA (const int& i) {
//        if (id == 1) printf("Rank 1 received from Rank 0 a total of %d data points : points %d through %d are %f %f %f %f %f %f %f \n",7*ARCount(0),7*i,7*i+6,GhostNodesAR(7*i),GhostNodesAR(7*i+1),GhostNodesAR(7*i+2),GhostNodesAR(7*i+3),GhostNodesAR(7*i+4),GhostNodesAR(7*i+5),GhostNodesAR(7*i+6));
//    });

    Kokkos::parallel_for("UpdateALoop",ARCount(0), KOKKOS_LAMBDA (const int& i) {
        int RankX = (int)(GhostNodesAR(7*i));
        int RankY = 0;
        int RankZ = (int)(GhostNodesAR(7*i+1));
        int CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
        //printf("CellLocation %d CellType %d \n",CellLocation,CellType(CellLocation));
        if ((cycle == 0)||((CellType(CellLocation) == Liquid)||(CellType(CellLocation) == Delayed)||(CellType(CellLocation) == LiqSol))) {
            //int GlobalX = RankX + MyXOffset;
            //int GlobalY = MyYOffset;
            // Update this ghost node cell's information with data from other rank
            CellType(CellLocation) = Active;
            GrainID(CellLocation) = (int)(GhostNodesAR(7*i+2));
            double cx = GhostNodesAR(7*i+3);
            double cy = GhostNodesAR(7*i+4);
            double cz = GhostNodesAR(7*i+5);
            DOCenter(3*CellLocation) = cx;
            DOCenter(3*CellLocation+1) = cy;
            DOCenter(3*CellLocation+2) = cz;
            int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
            DiagonalLength(CellLocation) = GhostNodesAR(7*i+6);
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
                double x0 = xp + NeighborX[n] - cx;
                double y0 = yp + NeighborY[n] - cy;
                double z0 = zp + NeighborZ[n] - cz;
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
        }
    });
    //MPI_Barrier(MPI_COMM_WORLD);
    //cout << "E" << id << endl;
    Kokkos::parallel_for("UpdateBLoop",BRCount(0), KOKKOS_LAMBDA (const int& i) {
        // GhostNodesBR: X coord, Z coord, GrainID, DOCenterX/Y/Z, DiagLength;
        int RankX = (int)(GhostNodesBR(7*i));
        int RankY = MyYSlices-1;
        int RankZ = (int)(GhostNodesBR(7*i+1));
        int CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
        if ((cycle == 0)||((CellType(CellLocation) == Liquid)||(CellType(CellLocation) == Delayed)||(CellType(CellLocation) == LiqSol))) {
            //if (cycle == 1003) cout << "Cell at " << RankX << " " << RankY << " " <<RankZ << " now active" << endl;
            //int GlobalX = MyXOffset + RankX;
            //int GlobalY = MyYOffset + RankY;
            // Update this ghost node cell's information with data from other rank
            CellType(CellLocation) = Active;
            GrainID(CellLocation) = (int)(GhostNodesBR(7*i+2));
            double cx = GhostNodesBR(7*i+3);
            double cy = GhostNodesBR(7*i+4);
            double cz = GhostNodesBR(7*i+5);
            //if (id == 2) cout << "RECIEVING CYCLE " << cycle << " Center " << cx << " " << cy << " " << cz << endl;
            DOCenter(3*CellLocation) = cx;
            DOCenter(3*CellLocation+1) = cy;
            DOCenter(3*CellLocation+2) = cz;
            int MyOrientation = GrainOrientation[((abs(GrainID(CellLocation)) - 1) % NGrainOrientations)];
            DiagonalLength(CellLocation) = GhostNodesBR(7*i+6);
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
                double x0 = xp + NeighborX[n] - cx;
                double y0 = yp + NeighborY[n] - cy;
                double z0 = zp + NeighborZ[n] - cz;
                // mag0 is the magnitude of (x0,y0,z0)

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
                //                                if (CDLVal == 0.0) printf("Zero CDLVal : %d %d %d %d %f %d %d %d %f %f %f",MyNeighborX,MyNeighborY,MyNeighborZ,n,mag0,index1,index2,index3,normx,normy,normz);
            }
            //                for (int l=0; l<26; l++) {
            //                    cout << "l = " << l << " CDL = " << CritDiagonalLength[BoxZ][RankX][RankY][l] << endl;
            //                }
        }
        //            else if (CellType[RankZ][RankX][RankY] == 'Q') {
        //                cout << "ID " << id << " competition for cell " << RankX << " " << RankY << " " << RankZ << endl;
        //                //CellType[RankZ][RankX][RankY] = Active;
        //            }
    });
    //MPI_Barrier(MPI_COMM_WORLD);
   // cout << "End" << id << endl;
}
