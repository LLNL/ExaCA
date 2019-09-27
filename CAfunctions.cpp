#include "header.h"
using namespace std;

double CrossP1(double TestVec1[3], double TestVec2[3]);
double CrossP2(double TestVec1[3], double TestVec2[3]);
double CrossP3(double TestVec1[3], double TestVec2[3]);
int FindItBounds(int RankX, int RankY, int MyXSlices, int MyYSlices);
int MaxIndex(double TestVec3[6]);
int XMPSlicesCalc(int p, int nx, int ProcessorsInXDirection, int ProcessorsInYDirection, int np, int DecompositionStrategy);
int XOffsetCalc(int p, int nx, int ProcessorsInXDirection, int ProcessorsInYDirection, int np, int DecompositionStrategy);
int YMPSlicesCalc(int p, int ny, int ProcessorsInYDirection, int np, int DecompositionStrategy);
int YOffsetCalc(int p, int ny, int ProcessorsInYDirection, int np, int DecompositionStrategy);
double MaxVal(double TestVec3[6],int NVals);
void NewGrain(int nx, int ny, int nz, int RankX, int RankY, int RankZ, int MyGrainID, int GlobalX, int GlobalY, ViewC::HostMirror CellType, ViewI::HostMirror GrainID, ViewF::HostMirror DiagonalLength, ViewF::HostMirror DOCenter);
void CritDiagLengthCalc(double xp, double yp, double zp, int MyOrientation, int NewCellX, int NewCellY, int NewCellZ, int CellLocation, double cx, double cy, double cz, int NeighborX[26], int NeighborY[26], int NeighborZ[26], float* GrainUnitVector, ViewI::HostMirror TriangleIndex, ViewF::HostMirror CritDiagonalLength);
void InitialDecomposition(int DecompositionStrategy, int nx, int ny, int &ProcessorsInXDirection, int &ProcessorsInYDirection, int id, int np, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset,int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn, int &MyLeftOut, int &MyRightIn, int &MyRightOut);
void XYLimitCalc(int &LLX, int &LLY, int &ULX, int &ULY, int MyXSlices, int MyYSlices, int MyLeft, int MyRight, int MyIn, int MyOut);

/*************************** FUNCTIONS CALLED THROUGH MAIN SUBROUTINES ***************************/

// Return the three terms of the cross product between two column vectors of length three
double CrossP1(double TestVec1[3], double TestVec2[3]) {
    double CrossPOut;
    CrossPOut = TestVec1[1]*TestVec2[2] - TestVec1[2]*TestVec2[1];
    return CrossPOut;
}
double CrossP2(double TestVec1[3], double TestVec2[3]) {
    double CrossPOut;
    CrossPOut = TestVec1[2]*TestVec2[0]-TestVec1[0]*TestVec2[2];
    return CrossPOut;
}
double CrossP3(double TestVec1[3], double TestVec2[3]) {
    double CrossPOut;
    CrossPOut = TestVec1[0]*TestVec2[1] - TestVec1[1]*TestVec2[0];
    return CrossPOut;
}

int XMPSlicesCalc(int p, int nx, int ProcessorsInXDirection, int ProcessorsInYDirection, int np, int DecompositionStrategy) {
    int XRemoteMPSlices;
    if (DecompositionStrategy == 1) {
        XRemoteMPSlices = nx;
    }
    else if (DecompositionStrategy >= 2) {
        int XSlicesPerP = nx/ProcessorsInXDirection;
        int XRemainder = nx % ProcessorsInXDirection;
        if (XRemainder == 0) {
            XRemoteMPSlices = XSlicesPerP;
        }
        else {
            int XPosition = p/ProcessorsInYDirection;
            if (XPosition < XRemainder) {
                XRemoteMPSlices = XSlicesPerP + 1;
            }
            else {
                XRemoteMPSlices = XSlicesPerP;
            }
        }
    }
    // Add "ghost nodes" for other processors
    XRemoteMPSlices = XRemoteMPSlices+2;
    return XRemoteMPSlices;
}

int XOffsetCalc(int p, int nx, int ProcessorsInXDirection, int ProcessorsInYDirection, int np, int DecompositionStrategy) {
    int RemoteXOffset;
    if (DecompositionStrategy == 1) {
        RemoteXOffset = 0;
    }
    else if (DecompositionStrategy >= 2) {
        int XSlicesPerP = nx/ProcessorsInXDirection;
        int XRemainder = nx % ProcessorsInXDirection;
        int XPosition = p/ProcessorsInYDirection;
        if (XRemainder == 0) {
            RemoteXOffset = XPosition*XSlicesPerP;
        }
        else {
            if (XRemainder > XPosition) {
                RemoteXOffset = XPosition*(XSlicesPerP + 1);
            }
            else {
                RemoteXOffset = (XSlicesPerP + 1)*(XRemainder-1) + XSlicesPerP + 1 + (XPosition-XRemainder)*XSlicesPerP;
            }
        }
    }
    // Account for "ghost nodes" for other processors
    RemoteXOffset--;
    // cout << " My ID is " << p << " and my X Offset is " << RemoteXOffset << endl;
    return RemoteXOffset;
}

int YMPSlicesCalc(int p, int ny, int ProcessorsInYDirection, int np, int DecompositionStrategy) {
    int YRemoteMPSlices;
    if (DecompositionStrategy == 1) {
        int YSlicesPerP = ny/np;
        int YRemainder = ny % np;
        if (YRemainder == 0) {
            YRemoteMPSlices = YSlicesPerP;
        }
        else {
            if (YRemainder > p) {
                YRemoteMPSlices = YSlicesPerP + 1;
            }
            else {
                YRemoteMPSlices = YSlicesPerP;
            }
        }
    }
    else if (DecompositionStrategy >= 2) {
        int YSlicesPerP = ny/ProcessorsInYDirection;
        int YRemainder = ny % ProcessorsInYDirection;
        if (YRemainder == 0) {
            YRemoteMPSlices = YSlicesPerP;
        }
        else {
            int YPosition = p % ProcessorsInYDirection;
            if (YPosition < YRemainder) {
                YRemoteMPSlices = YSlicesPerP + 1;
            }
            else {
                YRemoteMPSlices = YSlicesPerP;
            }
        }
    }
    // Add "ghost nodes" for other processors
    YRemoteMPSlices = YRemoteMPSlices+2;

    return YRemoteMPSlices;
}

int YOffsetCalc(int p, int ny, int ProcessorsInYDirection, int np, int DecompositionStrategy) {
    int RemoteYOffset;
    if (DecompositionStrategy == 1) {
        int YSlicesPerP = ny/np;
        int YRemainder = ny % np;
        if (YRemainder == 0) {
            RemoteYOffset = p*YSlicesPerP;
        }
        else {
            if (YRemainder > p) {
                RemoteYOffset = p*(YSlicesPerP + 1);
            }
            else {
                RemoteYOffset = (YSlicesPerP + 1)*(YRemainder-1) + YSlicesPerP + 1 + (p-YRemainder)*YSlicesPerP;
            }
        }
    }
    else if (DecompositionStrategy >= 2) {
        int YSlicesPerP = ny/ProcessorsInYDirection;
        int YRemainder = ny % ProcessorsInYDirection;
        int YPosition = p % ProcessorsInYDirection;
        if (YRemainder == 0) {
            RemoteYOffset = YPosition*YSlicesPerP;
        }
        else {
            if (YRemainder > YPosition) {
                RemoteYOffset = YPosition*(YSlicesPerP + 1);
            }
            else {
                RemoteYOffset = (YSlicesPerP + 1)*(YRemainder-1) + YSlicesPerP + 1 + (YPosition-YRemainder)*YSlicesPerP;
            }
        }
    }
    // Account for "ghost nodes" for other processors
    RemoteYOffset--;
    // cout << " My ID is " << p << " and my Y Offset is " << RemoteYOffset << endl;
    return RemoteYOffset;
}

int FindItBounds(int RankX, int RankY, int MyXSlices, int MyYSlices) {
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
    return ItBounds;
}

// Given a list of 6 values, return the index of the largest one
int MaxIndex(double TestVec3[6]) {
    int MaxIT = 0;
    for (int i=1; i<6; i++) {
        if (TestVec3[MaxIT] < TestVec3[i]) {
            MaxIT = i;
        }
    }
    return MaxIT;
}

// Given a list of 6 values, return the largest one
double MaxVal(double TestVec3[6], int NVals) {
    double MaxIT = TestVec3[0];
    for (int i=1; i<NVals; i++) {
        if (MaxIT < TestVec3[i]) {
            MaxIT = TestVec3[i];
        }
    }
    return MaxIT;
}

// Given a CA cell location and grain ID, create a new active cell and octahedron from the liquid
void NewGrain(int MyXSlices, int MyYSlices, int nz, int RankX, int RankY, int RankZ, int MyGrainID, int GlobalX, int GlobalY, ViewC::HostMirror CellType, ViewI::HostMirror GrainID, ViewF::HostMirror DiagonalLength, ViewF::HostMirror DOCenter) {
    int CellLocation = RankZ*MyXSlices*MyYSlices + MyYSlices*RankX + RankY;
    CellType(CellLocation) = 'A';
    GrainID(CellLocation) = MyGrainID;
    DiagonalLength(CellLocation) = 0.01;
    DOCenter(3*CellLocation) = GlobalX + 0.5;
    DOCenter(3*CellLocation+1) = GlobalY + 0.5;
    DOCenter(3*CellLocation+2) = RankZ + 0.5;
}

void CritDiagLengthCalc(double xp, double yp, double zp, int MyOrientation, int NewCellX, int NewCellY, int NewCellZ, int CellLocation, double cx, double cy, double cz, int NeighborX[26], int NeighborY[26], int NeighborZ[26], float* GrainUnitVector, ViewI::HostMirror TriangleIndex, ViewF::HostMirror CritDiagonalLength) {
    
    // Calculate critical diagonal lengths for the new active cell located at (xp,yp,zp) on the local grid
    // For each neighbor (l=0 to 25), calculate which octahedron face leads to cell capture
    // Calculate critical octahedron diagonal length to activate each nearest neighbor, as well as the coordinates of the triangle vertices on the capturing face
    
    // if ((id == 0)&&(cycle == 588)) cout << " Cell at " << RankX << " " << RankY << " " << RankZ << " has been captured by a grain with orientation " << MyOrientation << endl;
    //    if ((id == 0)&&(cycle > 0)) cout << " cycle " << cycle << " Orientation captured " << MyOrientation << endl;
    for (int l=0; l<26; l++)  {
        // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the image of the center of a neighbor cell
        double x0 = xp + NeighborX[l] - cx;
        double y0 = yp + NeighborY[l] - cy;
        double z0 = zp + NeighborZ[l] - cz;
        // mag0 is the magnitude of (x0,y0,z0)
        double mag0 = pow(pow(x0,2) + pow(y0,2) + pow(z0,2),0.5);
        if (mag0 == 0) {
            mag0 = max(0.0,0.00001);
            x0 = 0.00001;
        }
        
        // Calculate angles between the octahedron diagonal directions and the vector x0,y0,z0
        double AnglesA[6];
        for (int ll=0; ll<6; ll++) {
            double xd = GrainUnitVector[18*MyOrientation + 3*ll];
            double yd = GrainUnitVector[18*MyOrientation + 3*ll + 1];
            double zd = GrainUnitVector[18*MyOrientation + 3*ll + 2];
            AnglesA[ll] = (xd*x0 + yd*y0 + zd*z0)/mag0;
        }
        int index1, index2, index3;
        index1 = MaxIndex(AnglesA);
        TriangleIndex(78*CellLocation + 3*l) = index1;
        // First diagonal of the capturing face is that which makes the smallest (?) angle with x0,y0,z0
        double Diag1X = GrainUnitVector[18*MyOrientation + 3*index1 + 0];
        double Diag1Y = GrainUnitVector[18*MyOrientation + 3*index1 + 1];
        double Diag1Z = GrainUnitVector[18*MyOrientation + 3*index1 + 2];
        AnglesA[index1] = -1;
        if (index1 % 2 == 0) AnglesA[index1+1] = -1;
        if (index1 % 2 == 1) AnglesA[index1-1] = -1;
        double MaxA = MaxVal(AnglesA,6);
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
            index2 = MaxIndex(AnglesA);
            AnglesA[index2] = -1;
            if (index2 % 2 == 0) AnglesA[index2+1] = -1;
            if (index2 % 2 == 1) AnglesA[index2-1] = -1;
            index3 = MaxIndex(AnglesA);
        }

        TriangleIndex(78*CellLocation + 3*l + 1) = index2;
        Diag2X = GrainUnitVector[18*MyOrientation + 3*index2 + 0];
        Diag2Y = GrainUnitVector[18*MyOrientation + 3*index2 + 1];
        Diag2Z = GrainUnitVector[18*MyOrientation + 3*index2 + 2];
        TriangleIndex(78*CellLocation + 3*l + 2) = index3;
        // if (id == 1) cout << "RankX/RankY/BoxZ are " << RankX << " " << RankY << " " << BoxZ << " indexes are " << index1 << " " << index2 << " " << index3 << endl;
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
        UU[0] = CrossP1(U1,U2);
        UU[1] = CrossP2(U1,U2);
        UU[2] = CrossP3(U1,U2);
        double NDem = sqrt(UU[0]*UU[0] + UU[1]*UU[1] + UU[2]*UU[2]);
        Norm[0] = UU[0]/NDem;
        Norm[1] = UU[1]/NDem;
        Norm[2] = UU[2]/NDem;
        // normal to capturing plane
        double normx = Norm[0];
        double normy = Norm[1];
        double normz = Norm[2];
        double ParaT = (normx*x0+normy*y0+normz*z0)/(normx*Diag1X+normy*Diag1Y+normz*Diag1Z);
        //        if (isnan(pow(pow(ParaT*Diag1X,2) + pow(ParaT*Diag1Y,2) + pow(ParaT*Diag1Z,2),0.5))) {
        //            cout << "Cycle " << cycle << " id = " << id << " NaN CDL for cell " << NewCellX << " " << NewCellY << " " << BoxZ << endl;
        //        }
        float CDLVal = pow(pow(ParaT*Diag1X,2) + pow(ParaT*Diag1Y,2) + pow(ParaT*Diag1Z,2),0.5);
        if (std::isnan(CDLVal)) CDLVal = 0.0;
        CritDiagonalLength(26*CellLocation+l)  = CDLVal;

    }
}

/*******************************************************************************************************************/

// Determine the mapping of processors to grid data
void InitialDecomposition(int DecompositionStrategy, int nx, int ny, int &ProcessorsInXDirection, int &ProcessorsInYDirection, int id, int np, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset,int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn, int &MyLeftOut, int &MyRightIn, int &MyRightOut) {
    
    if (DecompositionStrategy == 1) {
    OneDim:
        ProcessorsInXDirection = 1;
        ProcessorsInYDirection = np;
        MyLeft = id-1;
        MyRight = id+1;
        if (id == 0) MyLeft = MPI_PROC_NULL;
        if (id == np-1) MyRight = MPI_PROC_NULL;
    }
    else if (DecompositionStrategy == 2) {
        if (np % 2 != 0) {
            if (id == 0) cout << "Warning: Number of processors not divisible by two- defaulting to 1D decomposition" << endl;
            goto OneDim;
        }
        ProcessorsInXDirection = 2;
        ProcessorsInYDirection = np/2;
        MyLeft = id-1;
        MyRight = id+1;
        MyIn = id+ProcessorsInYDirection;
        MyOut = id-ProcessorsInYDirection;
        MyLeftIn = MyIn-1;
        MyRightIn = MyIn+1;;
        MyLeftOut = MyOut-1;
        MyRightOut = MyOut+1;
        if (id % ProcessorsInYDirection == 0)     MyLeft = MPI_PROC_NULL;
        if (id % ProcessorsInYDirection == ProcessorsInYDirection-1)   MyRight = MPI_PROC_NULL;
        if (id < ProcessorsInYDirection)        MyOut = MPI_PROC_NULL;
        if (id >= (np-ProcessorsInYDirection))      MyIn = MPI_PROC_NULL;
        if ((MyLeft == MPI_PROC_NULL)||(MyIn == MPI_PROC_NULL)) MyLeftIn = MPI_PROC_NULL;
        if ((MyRight == MPI_PROC_NULL)||(MyIn == MPI_PROC_NULL)) MyRightIn = MPI_PROC_NULL;
        if ((MyLeft == MPI_PROC_NULL)||(MyOut == MPI_PROC_NULL)) MyLeftOut = MPI_PROC_NULL;
        if ((MyRight == MPI_PROC_NULL)||(MyOut == MPI_PROC_NULL)) MyRightOut = MPI_PROC_NULL;
    }
    else if (DecompositionStrategy == 3) {
        if (np % 2 != 0) {
            if (id == 0) cout << "Warning: Number of processors not divisible by two- defaulting to 1D decomposition" << endl;
            goto OneDim;
        }
        int PX, PY;
        double Square = sqrt(np);
        if (Square == round(Square)) {
            PX = Square;
            PY = Square;
        }
        else {
            int RoundSquare = floor(Square);
            int NotSquare = 1;
            while (NotSquare == 1) {
                if (np % RoundSquare == 0) {
                    PY = np/RoundSquare;
                    PX = RoundSquare;
                    //cout << PX << " " << PY << endl;
                    NotSquare = 0;
                }
                else {
                    RoundSquare--;
                }
            }
        }
        if (((nx > ny)&&(PY > PX))||((ny > nx)&&(PX > PY))) {
            ProcessorsInXDirection = PY;
            ProcessorsInYDirection = PX;
        }
        else {
            ProcessorsInXDirection = PX;
            ProcessorsInYDirection = PY;
        }
        if ((ProcessorsInXDirection == 2)&&(id == 0)) cout << "Warning: 2D decomposition pattern defaulting to 2DA" << endl;
        if (id == 0) cout << " Processors in X: " << ProcessorsInXDirection << " Processors in Y: " << ProcessorsInYDirection << endl;
        MyLeft = id-1;
        MyRight = id+1;
        MyIn = id+ProcessorsInYDirection;
        MyOut = id-ProcessorsInYDirection;
        MyLeftIn = MyIn-1;
        MyRightIn = MyIn+1;;
        MyLeftOut = MyOut-1;
        MyRightOut = MyOut+1;
        if (id % ProcessorsInYDirection == 0)     MyLeft = MPI_PROC_NULL;
        if (id % ProcessorsInYDirection == ProcessorsInYDirection-1)   MyRight = MPI_PROC_NULL;
        if (id < ProcessorsInYDirection)        MyOut = MPI_PROC_NULL;
        if (id >= (np-ProcessorsInYDirection))      MyIn = MPI_PROC_NULL;
        if ((MyLeft == MPI_PROC_NULL)||(MyIn == MPI_PROC_NULL)) MyLeftIn = MPI_PROC_NULL;
        if ((MyRight == MPI_PROC_NULL)||(MyIn == MPI_PROC_NULL)) MyRightIn = MPI_PROC_NULL;
        if ((MyLeft == MPI_PROC_NULL)||(MyOut == MPI_PROC_NULL)) MyLeftOut = MPI_PROC_NULL;
        if ((MyRight == MPI_PROC_NULL)||(MyOut == MPI_PROC_NULL)) MyRightOut = MPI_PROC_NULL;
    }
    // if (id == 0) cout << "PX = " << ProcessorsInXDirection << " PY = " << ProcessorsInYDirection << endl;
    //cout << "ID = " << id << " MyLeft = " << MyLeft << " MyRight = "  << MyRight << endl;
}


/*******************************************************************************************************************/

// Determine the scan limits for each rank in X and Y directions

void XYLimitCalc(int &LLX, int &LLY, int &ULX, int &ULY, int MyXSlices, int MyYSlices, int MyLeft, int MyRight, int MyIn, int MyOut) {
    LLX = 0;
    LLY = 0;
    ULX = MyXSlices-1;
    ULY = MyYSlices-1;
    if (MyLeft == MPI_PROC_NULL)  LLY = 1;
    if (MyRight == MPI_PROC_NULL) ULY = MyYSlices-2;
    if (MyIn == MPI_PROC_NULL) ULX = MyXSlices-2;
    if (MyOut == MPI_PROC_NULL) LLX = 1;
}
