#include "header.h"
using namespace std;


void TemperatureUpdate(int id, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, const int nz, int cycle, int &nn, double AConst, double BConst, double CConst, ViewI CritTimeStep, ViewC CellType, ViewF UndercoolingCurrent, ViewF UndercoolingChange, vector <int> NucLocI, vector <int> NucLocJ, vector <int> NucLocK, vector <int> NucleationTimes, vector <float> NucleationUndercooling, ViewI GrainID, int* GrainOrientation, ViewF DOCenter, int NeighborX[26], int NeighborY[26], int NeighborZ[26], float* GrainUnitVector, ViewI TriangleIndex, ViewF CritDiagonalLength, ViewF DiagonalLength, int NGrainOrientations) {

    int LocalDomainSize = nz*MyXSlices*MyYSlices;
    // Loop through CellType - any delay cells that are now liquid (below the liquidus)?
    // For the liquid and interface cells - update their local undercooling
    int NucleationThisDT = 0;
    Kokkos::parallel_reduce("TempUpdateLoop",LocalDomainSize, KOKKOS_LAMBDA (const int& D3D1ConvPosition, int &update) {
    //for (int D3D1ConvPosition=0; D3D1ConvPosition<LocalDomainSize; D3D1ConvPosition++) {
        int RankZ = floor(D3D1ConvPosition/(MyXSlices*MyYSlices));
        int Rem = D3D1ConvPosition % (MyXSlices*MyYSlices);
        int RankX = floor(Rem/MyYSlices);
        int RankY = Rem % MyYSlices;
        if (cycle > CritTimeStep(D3D1ConvPosition)) {

            if (CellType(D3D1ConvPosition) == 'D') {
                // This delayed liquid cell is at the liquidus and ready for potential solidification
                CellType(D3D1ConvPosition) = 'L';
                //if (id == 14) cout << RankZ << " liquid at cycle " << cycle << endl;
                UndercoolingCurrent(D3D1ConvPosition) = 0.0;
            }
            else if (CellType(D3D1ConvPosition) == 'N') {
                // Update local undercooling - linear cooling between solidus and liquidus
                UndercoolingCurrent(D3D1ConvPosition) += UndercoolingChange(D3D1ConvPosition);
                // What is the critical time step for this site to nucleate?
//                for (int nnumber=0; nnumber<NucLocI.size(); nnumber++) {
//                    if ((RankX == NucLocI[nnumber])&&(RankY == NucLocJ[nnumber])&&(RankZ == NucLocK[nnumber])) {
//                        if (cycle >= NucleationTimes[nnumber]) {
//
//                            // This undercooled liquid cell is now a nuclei
//                            update++;
//                            //cout << "New grain cycle " << cycle << " Rank " << id << " ID = " << GrainID[RankZ][RankX][RankY] << " Orientation " << GrainOrientation[GrainID[RankZ][RankX][RankY]] << endl;
//                            UndercoolingCurrent[D3D1ConvPosition] = NucleationUndercooling[nnumber];
//                            int GlobalX = RankX + MyXOffset;
//                            int GlobalY = RankY + MyYOffset;
//                            int MyGrainID = GrainID[D3D1ConvPosition];
//                            //cout << "N " << MyGrainID << endl;
//                            NewGrain(MyXSlices,MyYSlices,nz,RankX,RankY,RankZ,MyGrainID,GlobalX,GlobalY,CellType,GrainID,DiagonalLength,DOCenter);
//                            // The orientation for the new grain will depend on its Grain ID
//                            int MyOrientation = GrainOrientation[((abs(MyGrainID) - 1) % NGrainOrientations)];
//                            // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
//                            // (xp,yp,zp) is the new cell's center on the global grid
//                            double xp = GlobalX + 0.5;
//                            double yp = GlobalY + 0.5;
//                            double zp = RankZ + 0.5;
//                        CritDiagLengthCalc(xp,yp,zp,MyOrientation,RankX,RankY,RankZ,D3D1ConvPosition,DOCenter[3*D3D1ConvPosition],DOCenter[3*D3D1ConvPosition+1],DOCenter[3*D3D1ConvPosition+2],NeighborX,NeighborY,NeighborZ, GrainUnitVector,TriangleIndex,CritDiagonalLength);
//                        }
//                    }
//                }
            }
            else if (CellType(D3D1ConvPosition) == 'L') {
                // Update local undercooling - linear cooling between solidus and liquidus
                UndercoolingCurrent(D3D1ConvPosition) += UndercoolingChange(D3D1ConvPosition);
            }
            else if (CellType(D3D1ConvPosition) == 'A') {
                // Update local undercooling - linear cooling between solidus and liquidus
                UndercoolingCurrent(D3D1ConvPosition) += UndercoolingChange(D3D1ConvPosition);
                // Update local diagonal length of active cell
                double LocU = UndercoolingCurrent(D3D1ConvPosition);
                LocU = min(210.0,LocU);
                double V = AConst*pow(LocU,3) + BConst*pow(LocU,2) + CConst*LocU;
                V = max(0.0,V);
                DiagonalLength(D3D1ConvPosition) += min(0.045,V); // Max amount the diagonal can grow per time step
            }
        }
    },NucleationThisDT);
    nn += NucleationThisDT;

}


// Decentered octahedron algorithm for the capture of new interface cells by grains
//void CellCapture(int id, int cycle, int DecompositionStrategy, int &ACount, int &BCount, int &CCount, int &DCount, int &ECount, int &FCount, int &GCount, int &HCount, int MyXSlices, int MyYSlices, const int nz, int MyXOffset, int MyYOffset, int ItList[9][26], int NeighborX[26], int NeighborY[26], int NeighborZ[26], float* GrainUnitVector, ViewI::HostMirror  TriangleIndex, ViewF::HostMirror  CritDiagonalLength, ViewF::HostMirror  DiagonalLength, int* GrainOrientation, ViewC::HostMirror  CellType, ViewF::HostMirror  DOCenter, ViewI::HostMirror  GrainID, int NGrainOrientations, ViewF::HostMirror  UndercoolingCurrent) {
void CellCapture(int id, int cycle, int DecompositionStrategy, int &ACount, int &BCount, int &CCount, int &DCount, int &ECount, int &FCount, int &GCount, int &HCount, int MyXSlices, int MyYSlices, const int nz, int MyXOffset, int MyYOffset, int ItList[9][26], int NeighborX[26], int NeighborY[26], int NeighborZ[26], float* GrainUnitVector, ViewI TriangleIndex, ViewF CritDiagonalLength, ViewF DiagonalLength, int* GrainOrientation, ViewC CellType, ViewF DOCenter, ViewI GrainID, int NGrainOrientations, ViewF UndercoolingCurrent) {


    int LocalDomainSize = nz*MyXSlices*MyYSlices;
    // Cell capture - parallel loop over all type 'A' cells
    Kokkos::parallel_for("CellupdateLoop",LocalDomainSize, KOKKOS_LAMBDA (const int& D3D1ConvPosition) {
    //for (int D3D1ConvPosition=0; D3D1ConvPosition<LocalDomainSize; D3D1ConvPosition++) {
        int RankZ = floor(D3D1ConvPosition/(MyXSlices*MyYSlices));
        int Rem = D3D1ConvPosition % (MyXSlices*MyYSlices);
        int RankX = floor(Rem/MyYSlices);
        int RankY = Rem % MyYSlices;
        if (CellType(D3D1ConvPosition) == 'A') {

            // Cycle through all neigboring cells on this processor to see if they have been captured
            // Cells in ghost nodes cannot capture cells on other processors
            int LCount = 0;
            // Which neighbors should be iterated over?
            int ItBounds;
            // If X and Y coordinates are not on edges, Case 0: iteratation over neighbors 0-25 possible
            // If Y coordinate is on lower edge, Case 1: iteration over only neighbors 9-25 possible
            // If Y coordinate is on upper edge, Case 2: iteration over only neighbors 0-16 possible
            // If X coordinate is on lower edge, Case 3: iteration over only neighbors 0,1,3,4,6,8,9,10,11,13,15,17,18,20,21,22,24
            // If X coordinate is on upper edge, Case 4: iteration over only neighbors 0,2,3,4,5,7,9,10,12,14,16,17,19,20,21,23,25
            // If X/Y coordinates are on lower edge, Case 5: iteration over only neighbors 9,10,11,13,15,17,18,20,21,22,24
            // If X coordinate is on upper edge/Y on lower edge, Case 6:
            // If X coordinate is on lower edge/Y on upper edge, Case 7:
            // If X/Y coordinates are on upper edge, Case 8:
            if (RankY == 0) {
                if (RankX == 0) {
                    ItBounds = 5;
                }
                else if (RankX == MyXSlices-1) {
                    ItBounds = 6;
                }
                else {
                    ItBounds = 1;
                }
            }
            else if (RankY == MyYSlices-1) {
                if (RankX == 0) {
                    ItBounds = 7;
                }
                else if (RankX == MyXSlices-1) {
                    ItBounds = 8;
                }
                else {
                    ItBounds = 2;
                }
            }
            else {
                if (RankX == 0) {
                    ItBounds = 3;
                }
                else if (RankX == MyXSlices-1) {
                    ItBounds = 4;
                }
                else {
                    ItBounds = 0;
                }
            }
            int NListLength;
            if (ItBounds == 0) {
                NListLength = 26;
            }
            else if (ItBounds > 4) {
                NListLength = 11;
            }
            else {
                NListLength = 17;
            }
            // "ll" corresponds to the specific position on the list of neighboring cells
            for (int ll=0; ll<NListLength; ll++) {
                // "l" correpsponds to the specific neighboring cell
                int l = ItList[ItBounds][ll];
                // Local coordinates of adjacent cell center
                int MyNeighborX = RankX + NeighborX[l];
                int MyNeighborY = RankY + NeighborY[l];
                int MyNeighborZ = RankZ + NeighborZ[l];
                int NeighborD3D1ConvPosition = MyNeighborZ*MyXSlices*MyYSlices + MyNeighborX*MyYSlices + MyNeighborY;
                if ((CellType(NeighborD3D1ConvPosition) == 'L')||(CellType(NeighborD3D1ConvPosition) == 'N')||(CellType(NeighborD3D1ConvPosition) == 'D')) {
                    LCount = 1;
                    if (DiagonalLength(D3D1ConvPosition) >= CritDiagonalLength(26*D3D1ConvPosition+l)) {
                        
                        // atomic operation - new cell's state known across GPU
                        Kokkos::atomic_store(CellType(NeighborD3D1ConvPosition),'A');
                        //CellType(NeighborD3D1ConvPosition) = 'A';
                        
                        int GlobalX = RankX + MyXOffset;
                        int GlobalY = RankY + MyYOffset;
                        int h = GrainID(D3D1ConvPosition);
                        int MyOrientation = GrainOrientation[((abs(h) - 1) % NGrainOrientations)];

                        // The new cell is captured by this cell's growing octahedron (Grain "h")
                        GrainID(NeighborD3D1ConvPosition) = h;
                        // (cxold, cyold, czold) are the coordiantes of this decentered octahedron
                        double cxold = DOCenter(3*D3D1ConvPosition);
                        double cyold = DOCenter(3*D3D1ConvPosition+1);
                        double czold = DOCenter(3*D3D1ConvPosition+2);

                        // (xp,yp,zp) are the global coordinates of the new cell's center
                        double xp = GlobalX + NeighborX[l] + 0.5;
                        double yp = GlobalY + NeighborY[l] + 0.5;
                        double zp = RankZ + NeighborZ[l] + 0.5;
                        // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the image of the center of the new  cell
                        double x0 = xp - cxold;
                        double y0 = yp - cyold;
                        double z0 = zp - czold;

                        int index1 = TriangleIndex(78*D3D1ConvPosition + 3*l);
                        int index2 = TriangleIndex(78*D3D1ConvPosition + 3*l + 1);
                        int index3 = TriangleIndex(78*D3D1ConvPosition + 3*l + 2);

                        double Diag1X = GrainUnitVector[18*MyOrientation + 3*index1];
                        double Diag1Y = GrainUnitVector[18*MyOrientation + 3*index1 + 1];
                        double Diag1Z = GrainUnitVector[18*MyOrientation + 3*index1 + 2];

                        double Diag2X = GrainUnitVector[18*MyOrientation + 3*index2];
                        double Diag2Y = GrainUnitVector[18*MyOrientation + 3*index2 + 1];
                        double Diag2Z = GrainUnitVector[18*MyOrientation + 3*index2 + 2];
                        double Diag3X = GrainUnitVector[18*MyOrientation + 3*index3];
                        double Diag3Y = GrainUnitVector[18*MyOrientation + 3*index3 + 1];
                        double Diag3Z = GrainUnitVector[18*MyOrientation + 3*index3 + 2];

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
                        double TriangleX1 = cxold+ParaT*Diag1X;
                        double TriangleX2 = cxold+ParaT*Diag2X;
                        double TriangleX3 = cxold+ParaT*Diag3X;
                        double TriangleY1 = cyold+ParaT*Diag1Y;
                        double TriangleY2 = cyold+ParaT*Diag2Y;
                        double TriangleY3 = cyold+ParaT*Diag3Y;
                        double TriangleZ1 = czold+ParaT*Diag1Z;
                        double TriangleZ2 = czold+ParaT*Diag2Z;
                        double TriangleZ3 = czold+ParaT*Diag3Z;

                        // Determine which of the 3 corners of the capturing face is closest to the captured cell center
                        double Disttocorner0 = pow(pow(TriangleX1-xp,2) + pow(TriangleY1-yp,2) + pow(TriangleZ1-zp,2),0.5);
                        double Disttocorner1 = pow(pow(TriangleX2-xp,2) + pow(TriangleY2-yp,2) + pow(TriangleZ2-zp,2),0.5);
                        double Disttocorner2 = pow(pow(TriangleX3-xp,2) + pow(TriangleY3-yp,2) + pow(TriangleZ3-zp,2),0.5);

                        int mindisttocornerindex;
                        double mindisttocorner, xc, yc, zc;

                        if (Disttocorner0 < Disttocorner1) {
                            if (Disttocorner2 < Disttocorner0) {
                                mindisttocornerindex = 2;
                                mindisttocorner = Disttocorner2;
                                xc = TriangleX3;
                                yc = TriangleY3;
                                zc = TriangleZ3;
                            }
                            else {
                                mindisttocornerindex = 0;
                                mindisttocorner = Disttocorner0;
                                xc = TriangleX1;
                                yc = TriangleY1;
                                zc = TriangleZ1;
                            }
                        }
                        else {
                            if (Disttocorner2 < Disttocorner1) {
                                mindisttocornerindex = 2;
                                mindisttocorner = Disttocorner2;
                                xc = TriangleX3;
                                yc = TriangleY3;
                                zc = TriangleZ3;
                            }
                            else   {
                                mindisttocornerindex = 1;
                                mindisttocorner = Disttocorner1;
                                xc = TriangleX2;
                                yc = TriangleY2;
                                zc = TriangleZ2;
                            }
                        }

                        double x1, y1, z1, x2, y2, z2;
                        if (mindisttocornerindex == 0) {
                            x1 = TriangleX2;
                            y1 = TriangleY2;
                            z1 = TriangleZ2;
                            x2 = TriangleX3;
                            y2 = TriangleY3;
                            z2 = TriangleZ3;
                        }

                        if (mindisttocornerindex == 1) {
                            x1 = TriangleX1;
                            y1 = TriangleY1;
                            z1 = TriangleZ1;
                            x2 = TriangleX3;
                            y2 = TriangleY3;
                            z2 = TriangleZ3;
                        }

                        if (mindisttocornerindex == 2) {
                            x1 = TriangleX1;
                            y1 = TriangleY1;
                            z1 = TriangleZ1;
                            x2 = TriangleX2;
                            y2 = TriangleY2;
                            z2 = TriangleZ2;
                        }

                        double D1 = pow(pow(xp-x2,2) + pow(yp-y2,2) + pow(zp-z2,2),0.5);
                        double D2 = pow(pow(xc-x2,2) + pow(yc-y2,2) + pow(zc-z2,2),0.5);
                        double D3 = pow(pow(xp-x1,2) + pow(yp-y1,2) + pow(zp-z1,2),0.5);
                        double D4 = pow(pow(xc-x1,2) + pow(yc-y1,2) + pow(zc-z1,2),0.5);

                        double I1, I2, J1, J2;
                        // If minimum distance to corner = 0, the octahedron corner captured the new cell center
                        if (mindisttocorner == 0) {
                            I1 = 0;
                            I2 = D2;
                            J1 = 0;
                            J2 = D4;
                        }
                        else {
                            I1 = D1*((xp-x2)*(xc-x2) + (yp-y2)*(yc-y2) + (zp-z2)*(zc-z2))/(D1*D2);
                            I2 = D2 - I1;
                            J1 = D3*((xp-x1)*(xc-x1) + (yp-y1)*(yc-y1) + (zp-z1)*(zc-z1))/(D3*D4);
                            J2 = D4 - J1;
                        }
                        double L12 = 0.5*(min(I1,sqrt(3)) + min(I2,sqrt(3)));
                        double L13 = 0.5*(min(J1,sqrt(3)) + min(J2,sqrt(3)));
                        double NewODiagL = sqrt(2)*max(L12,L13); // half diagonal length of new octahedron

                        DiagonalLength(NeighborD3D1ConvPosition) = NewODiagL;
                        //if ((id == 1)&&(cycle == 6403)) cout << "Cell calculated New diag" << endl;
                        // Calculate coordinates of new decentered octahedron center
                        double CaptDiag[3], CaptDiagUV[3];
                        CaptDiag[0] = xc - cxold;
                        CaptDiag[1] = yc - cyold;
                        CaptDiag[2] = zc - czold;

                        NDem = sqrt(CaptDiag[0]*CaptDiag[0] + CaptDiag[1]*CaptDiag[1] + CaptDiag[2]*CaptDiag[2]);
                        CaptDiagUV[0] = CaptDiag[0]/NDem;
                        CaptDiagUV[1] = CaptDiag[1]/NDem;
                        CaptDiagUV[2] = CaptDiag[2]/NDem;
                        // (cx, cy, cz) are the coordiantes of the new active cell's decentered octahedron
                        double cx = xc - NewODiagL*CaptDiagUV[0];
                        double cy = yc - NewODiagL*CaptDiagUV[1];
                        double cz = zc - NewODiagL*CaptDiagUV[2];

                        DOCenter(3*NeighborD3D1ConvPosition) = cx;
                        DOCenter(3*NeighborD3D1ConvPosition+1) = cy;
                        DOCenter(3*NeighborD3D1ConvPosition+2) = cz;

                        // Calculate critical diagonal lengths for the new active cell located at (xp,yp,zp) on the local grid
                        // For each neighbor (l=0 to 25), calculate which octahedron face leads to cell capture
                        // Calculate critical octahedron diagonal length to activate each nearest neighbor, as well as the coordinates of the triangle vertices on the capturing face
                        for (int n=0; n<26; n++)  {
                            // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the image of the center of a neighbor cell
                            double x0 = xp + NeighborX[n] - cx;
                            double y0 = yp + NeighborY[n] - cy;
                            double z0 = zp + NeighborZ[n] - cz;
                            // mag0 is the magnitude of (x0,y0,z0)
                            double mag0 = pow(pow(x0,2) + pow(y0,2) + pow(z0,2),0.5);
                            if (mag0 == 0) {
                                mag0 = max(0.0,0.00001);
                                x0 = 0.00001;
                            }
                            
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
                            

                            TriangleIndex(78*NeighborD3D1ConvPosition + 3*n) = index1;
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
                            
                            TriangleIndex(78*NeighborD3D1ConvPosition + 3*n + 1) = index2;
                            Diag2X = GrainUnitVector[18*MyOrientation + 3*index2 + 0];
                            Diag2Y = GrainUnitVector[18*MyOrientation + 3*index2 + 1];
                            Diag2Z = GrainUnitVector[18*MyOrientation + 3*index2 + 2];
                            TriangleIndex(78*NeighborD3D1ConvPosition + 3*n + 2) = index3;
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
                            if ((normx*Diag1X+normy*Diag1Y+normz*Diag1Z) == 0) CDLVal = 0.0;
                            CritDiagonalLength(26*NeighborD3D1ConvPosition+n)  = CDLVal;
                        }

                        // Collect data for the ghost nodes:
                        if (DecompositionStrategy == 1) {
                            if (MyNeighborY == 1) {
                                CellType(NeighborD3D1ConvPosition) = '1';
                            }
                            else if (MyNeighborY == MyYSlices-2) {
                                CellType(NeighborD3D1ConvPosition) = '1';
                            }
                        }
                        else {
                            if (MyNeighborY == 1) {
                                // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                                if (MyNeighborX == MyXSlices-2) {
                                    CellType(NeighborD3D1ConvPosition) = '3';
                                }
                                else if (MyNeighborX == 1) {
                                    CellType(NeighborD3D1ConvPosition) = '3';
                                }
                                else if ((MyNeighborX > 1)&&(MyNeighborX < MyXSlices-2)) {
                                    // This is being sent to MyLeft
                                    CellType(NeighborD3D1ConvPosition) = '1';
                                }
                            }
                            else if (MyNeighborY == MyYSlices-2) {
                                // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                                if (MyNeighborX == MyXSlices-2) {
                                    CellType(NeighborD3D1ConvPosition) = '3';
                                }
                                else if (MyNeighborX == 1) {
                                    CellType(NeighborD3D1ConvPosition) = '3';
                                }
                                else if ((MyNeighborX > 1)&&(MyNeighborX < MyXSlices-2)) {
                                    CellType(NeighborD3D1ConvPosition) = '1';
                                }
                            }
                            else if ((MyNeighborX == 1)&&(MyNeighborY > 1)&&(MyNeighborY < MyYSlices-2)) {
                                CellType(NeighborD3D1ConvPosition) = '1';
                                //if (id == 0) cout << "RANK 0 LISTED " << MyNeighborX << " " << MyNeighborY << " " << MyNeighborZ << endl;
                            }
                            else if ((MyNeighborX == MyXSlices-2)&&(MyNeighborY > 1)&&(MyNeighborY < MyYSlices-2)) {
                                CellType(NeighborD3D1ConvPosition) = '1';
                            }
                        } // End if statement for ghost node marking
                    } // End if statement for capture event
                } // End if statement for neighbors of type L and N
            } // End loop over all neighbors of this active cell
            if (LCount == 0) {
                // This active cell has no more neighboring cells to be captured, becomes solid
                CellType(D3D1ConvPosition) = 'S';
            }
        } // end "if" loop for Active type cells
    }); // end Kokkos parallel for loop over all cells on this MPI rank

}


// Prints intermediate code output to stdout, checks to see if solidification is complete
 void IntermediateOutputAndCheck(int id, int cycle, int MyXSlices, int MyYSlices, int nz, int nn, int &XSwitch, ViewC::HostMirror CellType) {

    signed long int Type_sumL = 0;
    signed long int Type_sumD = 0;
    signed long int Type_sumA = 0;
    signed long int Global_sumL, Global_sumD, Global_sumA, Global_nn;
    Global_nn = 0;
    for (int k=1; k<nz-1; k++)   {
        for(int i=1; i<MyXSlices-1; i++)  {
            for(int j=1; j<MyYSlices-1; j++)  {
                int D3D1ConvPosition = k*MyXSlices*MyYSlices + i*MyYSlices + j;
                if (CellType(D3D1ConvPosition) == 'D') {
                    Type_sumD += 1;
                }
                else if (CellType(D3D1ConvPosition) == 'L') {
                    Type_sumL += 1;
                }
                else if (CellType(D3D1ConvPosition) == 'A') {
                    Type_sumA += 1;
                }
            }
        }
    }
//     if (id == 57) {
//         cout << "======================================================" << endl;
//         cout << "On rank 57" << endl;
//         cout << "Superheated liquid cells remaining = " << Type_sumD << endl;
//         cout << "Inactive interface cells remaining = " << Type_sumI << endl;
//         cout << "Undercooled liquid cells remaining = " << Type_sumL << endl;
//         cout << "Active interface cells remaining = " << Type_sumA << endl;
//         //        cout << "Liquid cells remaining = " << Global_sumL + Global_sumD << endl;
//         //        cout << "Number of nucleation events = " << Global_nn << endl;
//         cout << "======================================================" << endl;
//     }
    MPI_Reduce(&Type_sumL,&Global_sumL,1,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&Type_sumD,&Global_sumD,1,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    //MPI_Reduce(&Type_sumI,&Global_sumI,1,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&Type_sumA,&Global_sumA,1,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&nn,&Global_nn,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    if (id == 0) {
        cout << "======================================================" << endl;
        cout << "cycle = " << cycle << endl;
        cout << "Superheated liquid cells remaining = " << Global_sumD << endl;
        //cout << "Inactive interface cells remaining = " << Global_sumI << endl;
        cout << "Undercooled liquid cells remaining = " << Global_sumL << endl;
        cout << "Active interface cells remaining = " << Global_sumA << endl;
//        cout << "Liquid cells remaining = " << Global_sumL + Global_sumD << endl;
        cout << "Number of nucleation events = " << Global_nn << endl;
        cout << "======================================================" << endl;
        if (Global_sumL + Global_sumD == 0) {
            XSwitch = 1;
        }
    }

    MPI_Bcast(&XSwitch,1,MPI_INT,0,MPI_COMM_WORLD);

}
