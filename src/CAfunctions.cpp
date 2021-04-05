#include "header.h"
using namespace std;

/*************************** FUNCTIONS CALLED THROUGH MAIN SUBROUTINES ***************************/

//*****************************************************************************/
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

//*****************************************************************************/
int XMPSlicesCalc(int p, int nx, int ProcessorsInXDirection, int ProcessorsInYDirection, int DecompositionStrategy) {
    int XRemoteMPSlices = 0;
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


//*****************************************************************************/
int XOffsetCalc(int p, int nx, int ProcessorsInXDirection, int ProcessorsInYDirection, int DecompositionStrategy) {
    int RemoteXOffset = 0;
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
    return RemoteXOffset;
}

//*****************************************************************************/
int YMPSlicesCalc(int p, int ny, int ProcessorsInYDirection, int np, int DecompositionStrategy) {
    int YRemoteMPSlices = 0;
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

//*****************************************************************************/
int YOffsetCalc(int p, int ny, int ProcessorsInYDirection, int np, int DecompositionStrategy) {
    int RemoteYOffset = 0;
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
    return RemoteYOffset;
}

//*****************************************************************************/
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

//*****************************************************************************/
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

//*****************************************************************************/
// Determine the mapping of processors to grid data
void InitialDecomposition(int &DecompositionStrategy, int nx, int ny, int &ProcessorsInXDirection, int &ProcessorsInYDirection, int id, int np, int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn, int &MyLeftOut, int &MyRightIn, int &MyRightOut) {
    
    if (DecompositionStrategy == 1) {
    OneDim:
        ProcessorsInXDirection = 1;
        ProcessorsInYDirection = np;
        if (np > 1) {
            MyLeft = id-1;
            MyRight = id+1;
            if (id == 0) MyLeft = MPI_PROC_NULL;
            if (id == np-1) MyRight = MPI_PROC_NULL;
            MyIn = MPI_PROC_NULL;
            MyRightIn = MPI_PROC_NULL;
            MyOut = MPI_PROC_NULL;
            MyRightOut = MPI_PROC_NULL;
            MyLeftIn = MPI_PROC_NULL;
            MyLeftOut = MPI_PROC_NULL;
        }
        else {
            // No MPI communication
            MyLeft = MPI_PROC_NULL;
            MyRight = MPI_PROC_NULL;
            MyIn = MPI_PROC_NULL;
            MyRightIn = MPI_PROC_NULL;
            MyOut = MPI_PROC_NULL;
            MyRightOut = MPI_PROC_NULL;
            MyLeftIn = MPI_PROC_NULL;
            MyLeftOut = MPI_PROC_NULL;
        }
    }
    else if (DecompositionStrategy == 2) {
        if (np % 2 != 0) {
            if (id == 0) cout << "Warning: Number of processors not divisible by two- defaulting to 1D decomposition" << endl;
            DecompositionStrategy = 1;
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
        
        // Is the number of ranks a perfect square or divisible by 2?
        bool DivisibleBy2, PerfectSquare;
        if (np % 2 != 0) {
            DivisibleBy2 = false;
        }
        else {
            DivisibleBy2 = true;
        }
        double Square = sqrt(np);
        if (Square == round(Square)) {
            PerfectSquare = true;
        }
        else {
            PerfectSquare = false;
        }
        
        if (((!PerfectSquare)&&(!DivisibleBy2))||(np == 2)) {
            if (id == 0) cout << "Note: Number of processors is either 2, not divisible by two, or not divisible by itself- defaulting to 1D decomposition" << endl;
            DecompositionStrategy = 1;
            goto OneDim;
        }
        int PX, PY;

        if (PerfectSquare) {
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
        if ((ProcessorsInXDirection == 2)&&(id == 0)) cout << "Note: Decomposition pattern 3 is equivalent to 2" << endl;
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
}


//*****************************************************************************/
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