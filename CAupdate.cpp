#include "header.h"
using namespace std;

namespace sample {  // namespace helps with name resolution in reduction identity
    template< class ScalarType, int N >
    struct array_type {
        ScalarType the_array[N];
        
        KOKKOS_INLINE_FUNCTION   // Default constructor - Initialize to 0's
        array_type() {
            for (int i = 0; i < N; i++ ) { the_array[i] = 0; }
        }
        KOKKOS_INLINE_FUNCTION   // Copy Constructor
        array_type(const array_type & rhs) {
            for (int i = 0; i < N; i++ ){
                the_array[i] = rhs.the_array[i];
            }
        }
        KOKKOS_INLINE_FUNCTION   // add operator
        array_type& operator += (const array_type& src) {
            for ( int i = 0; i < N; i++ ) {
                the_array[i]+=src.the_array[i];
            }
            return *this;
        }
        KOKKOS_INLINE_FUNCTION   // volatile add operator
        void operator += (const volatile array_type& src) volatile {
            for ( int i = 0; i < N; i++ ) {
                the_array[i]+=src.the_array[i];
            }
        }
    };
    typedef array_type<int,8> ValueType;  // used to simplify code below
}
namespace Kokkos { //reduction identity must be defined in Kokkos namespace
    template<>
    struct reduction_identity< sample::ValueType > {
        KOKKOS_FORCEINLINE_FUNCTION static sample::ValueType sum() {
            return sample::ValueType();
        }
    };
}
    
void Nucleation(int id, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, const int nz, int cycle, int &nn, ViewI CritTimeStep, ViewI CellType, ViewF UndercoolingCurrent, ViewF UndercoolingChange, ViewI NucleiLocations, ViewI NucleationTimes, ViewI GrainID, int* GrainOrientation, ViewF DOCenter, int NeighborX[26], int NeighborY[26], int NeighborZ[26], float* GrainUnitVector, ViewF CritDiagonalLength, ViewF DiagonalLength, int NGrainOrientations, int PossibleNuclei_ThisRank, ViewI Locks) {
    
    // Loop through local list of nucleation events - has the time step exceeded the time step for nucleation at the sites?
    int NucleationThisDT = 0;
    Kokkos::parallel_reduce("NucleiUpdateLoop",PossibleNuclei_ThisRank, KOKKOS_LAMBDA (const int& NucCounter, int &update) {
        //printf("Nucleation Check: rank, cycle, i %d %d %d \n",id,cycle,NucCounter);
        if ((cycle >= NucleationTimes(NucCounter))&&(CellType(NucleiLocations(NucCounter)) == LiqSol)) {
            Locks(NucleiLocations(NucCounter)) = 1;
            //printf("Nucleation: rank, cycle, NucTime %d %d %d",id,cycle,NucleationTimes(NucCounter));
            long int D3D1ConvPosition = NucleiLocations(NucCounter);
            // This undercooled liquid cell is now a nuclei (add to count if it isn't in the ghost nodes, to avoid double counting)
            int RankZ = floor(D3D1ConvPosition/(MyXSlices*MyYSlices));
            int Rem = D3D1ConvPosition % (MyXSlices*MyYSlices);
            int RankX = floor(Rem/MyYSlices);
            int RankY = Rem % MyYSlices;
            if ((RankX > 0)&&(RankX < MyXSlices-1)&&(RankY > 0)&&(RankY < MyYSlices-1)) update++;
            int GlobalX = RankX + MyXOffset;
            int GlobalY = RankY + MyYOffset;
            int MyGrainID = GrainID[D3D1ConvPosition];
            CellType(D3D1ConvPosition) = Active;
            GrainID(D3D1ConvPosition) = MyGrainID;
            DiagonalLength(D3D1ConvPosition) = 0.01;
            long int DOX = (long int)(3)*D3D1ConvPosition;
            long int DOY = (long int)(3)*D3D1ConvPosition+(long int)(1);
            long int DOZ = (long int)(3)*D3D1ConvPosition+(long int)(2);
            DOCenter(DOX) = GlobalX + 0.5;
            DOCenter(DOY) = GlobalY + 0.5;
            DOCenter(DOZ) = RankZ + 0.5;
            // The orientation for the new grain will depend on its Grain ID (nucleated grains have negative GrainID values)
            int MyOrientation = GrainOrientation[((abs(MyGrainID) - 1) % NGrainOrientations)];
            // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
            for (int n=0; n<26; n++)  {

//                int MyNeighborX = RankX + NeighborX[n];
//                int MyNeighborY = RankY + NeighborY[n];
//                int MyNeighborZ = RankZ + NeighborZ[n];

                // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the image of the center of a neighbor cell
                double x0 = NeighborX[n];
                double y0 = NeighborY[n];
                double z0 = NeighborZ[n];

                // mag0 is the magnitude of (x0,y0,z0)
                double mag0 = pow(pow(x0,2) + pow(y0,2) + pow(z0,2),0.5);

                // Calculate unit vectors for the octahedron that intersect the new cell center
                double Diag1X, Diag1Y, Diag1Z, Diag2X, Diag2Y, Diag2Z, Diag3X, Diag3Y, Diag3Z;
                double Angle1 = (GrainUnitVector[18*MyOrientation]*x0 + GrainUnitVector[18*MyOrientation + 1]*y0 + GrainUnitVector[18*MyOrientation + 2]*z0)/mag0;
                if (Angle1 < 0) {
                    Diag1X = GrainUnitVector[18*MyOrientation];
                    Diag1Y = GrainUnitVector[18*MyOrientation + 1];
                    Diag1Z = GrainUnitVector[18*MyOrientation + 2];
                }
                else {
                    Diag1X = GrainUnitVector[18*MyOrientation + 3];
                    Diag1Y = GrainUnitVector[18*MyOrientation + 4];
                    Diag1Z = GrainUnitVector[18*MyOrientation + 5];
                }
                
                double Angle2 = (GrainUnitVector[18*MyOrientation + 6]*x0 + GrainUnitVector[18*MyOrientation + 7]*y0 + GrainUnitVector[18*MyOrientation + 8]*z0)/mag0;
                if (Angle2 < 0) {
                    Diag1X = GrainUnitVector[18*MyOrientation + 6];
                    Diag1Y = GrainUnitVector[18*MyOrientation + 7];
                    Diag1Z = GrainUnitVector[18*MyOrientation + 8];
                }
                else {
                    Diag1X = GrainUnitVector[18*MyOrientation + 9];
                    Diag1Y = GrainUnitVector[18*MyOrientation + 10];
                    Diag1Z = GrainUnitVector[18*MyOrientation + 11];
                }
                
                double Angle3 = (GrainUnitVector[18*MyOrientation + 12]*x0 + GrainUnitVector[18*MyOrientation + 13]*y0 + GrainUnitVector[18*MyOrientation + 14]*z0)/mag0;
                if (Angle3 < 0) {
                    Diag2X = GrainUnitVector[18*MyOrientation + 12];
                    Diag2Y = GrainUnitVector[18*MyOrientation + 13];
                    Diag2Z = GrainUnitVector[18*MyOrientation + 14];
                }
                else {
                    Diag3X = GrainUnitVector[18*MyOrientation + 15];
                    Diag3Y = GrainUnitVector[18*MyOrientation + 16];
                    Diag3Z = GrainUnitVector[18*MyOrientation + 17];
                }
                
                // Calculate angles between the octahedron diagonal directions and the vector x0,y0,z0
//                double AnglesA[6];
//                for (int aa=0; aa<6; aa++) {
//                    double xd = GrainUnitVector[18*MyOrientation + 3*aa];
//                    double yd = GrainUnitVector[18*MyOrientation + 3*aa + 1];
//                    double zd = GrainUnitVector[18*MyOrientation + 3*aa + 2];
//                    AnglesA[aa] = (xd*x0 + yd*y0 + zd*z0)/mag0;
//                }
//                int index1, index2, index3;
//                index1 = 0;
//                for (int ii=1; ii<6; ii++) {
//                    if (AnglesA[index1] < AnglesA[ii]) {
//                        index1 = ii;
//                    }
//                }
//
//                //TriangleIndex(78*D3D1ConvPosition + 3*n) = index1;
//                // First diagonal of the capturing face is that which makes the smallest (?) angle with x0,y0,z0
//                double Diag1X = GrainUnitVector[18*MyOrientation + 3*index1 + 0];
//                double Diag1Y = GrainUnitVector[18*MyOrientation + 3*index1 + 1];
//                double Diag1Z = GrainUnitVector[18*MyOrientation + 3*index1 + 2];
//                AnglesA[index1] = -1;
//                if (index1 % 2 == 0) AnglesA[index1+1] = -1;
//                if (index1 % 2 == 1) AnglesA[index1-1] = -1;
//
//                double MaxA = AnglesA[0];
//                for (int ii=1; ii<6; ii++) {
//                    if (MaxA < AnglesA[ii]) {
//                        MaxA = AnglesA[ii];
//                    }
//                }
//
//
//                double Diag2X, Diag2Y, Diag2Z, Diag3X, Diag3Y, Diag3Z;
//                if (MaxA == 0) {
//                    // Special case- other diagonals are all perpendicular to the first one (e.g. the octahedron corner captures the new cell center)
//                    // manually assign other diagonals (one of the 4 possible "capturing" faces)
//                    if ((index1 == 0)||(index1 == 1)) {
//                        index2 = 2;
//                        index3 = 4;
//                    }
//                    else if ((index1 == 2)||(index1 == 3)) {
//                        index2 = 0;
//                        index3 = 4;
//                    }
//                    else if ((index1 == 4)||(index1 == 5)) {
//                        index2 = 0;
//                        index3 = 2;
//                    }
//                }
//                else {
//                    // 2nd and 3rd closest diagonals to x0,y0,z0
//                    // note that if these are the same length, this means that the octahedron edge captures the new cell center
//                    // in this case, either of the 2 possible "capturing" faces will work
//                    index2 = 0;
//                    for (int ii=1; ii<6; ii++) {
//                        if (AnglesA[index2] < AnglesA[ii]) {
//                            index2 = ii;
//                        }
//                    }
//                    AnglesA[index2] = -1;
//                    if (index2 % 2 == 0) AnglesA[index2+1] = -1;
//                    if (index2 % 2 == 1) AnglesA[index2-1] = -1;
//                    // index3 = MaxIndex(AnglesA);
//                    index3 = 0;
//                    for (int ii=1; ii<6; ii++) {
//                        if (AnglesA[index3] < AnglesA[ii]) {
//                            index3 = ii;
//                        }
//                    }
//
//                }
//                //TriangleIndex(78*D3D1ConvPosition + 3*n + 1) = index2;
//                Diag2X = GrainUnitVector[18*MyOrientation + 3*index2 + 0];
//                Diag2Y = GrainUnitVector[18*MyOrientation + 3*index2 + 1];
//                Diag2Z = GrainUnitVector[18*MyOrientation + 3*index2 + 2];
//                //TriangleIndex(78*D3D1ConvPosition + 3*n + 2) = index3;
//                Diag3X = GrainUnitVector[18*MyOrientation + 3*index3 + 0];
//                Diag3Y = GrainUnitVector[18*MyOrientation + 3*index3 + 1];
//                Diag3Z = GrainUnitVector[18*MyOrientation + 3*index3 + 2];
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
                long int CDLIndex = (long int)(26)*D3D1ConvPosition + (long int)(n);
                CritDiagonalLength(CDLIndex) = CDLVal;
                //printf("ID = %d Orient = %d GID = %d Diag %d CDL %f \n",id,MyOrientation,MyGrainID,n,CDLVal);
                //printf("CDLVal : %d %d %d %d %f %d %d %d %f %f %f",MyNeighborX,MyNeighborY,MyNeighborZ,n,mag0,index1,index2,index3,normx,normy,normz);
            }
        }
    },NucleationThisDT);
    nn += NucleationThisDT;

//    MPI_Barrier(MPI_COMM_WORLD);
//    if (id == 0) cout << "End nucleation " << cycle << endl;
    
}
    
// Decentered octahedron algorithm for the capture of new interface cells by grains
void CellCapture(int id, int np, int cycle, int DecompositionStrategy, int LocalDomainSize, int MyXSlices, int MyYSlices, const int nz, double AConst, double BConst, double CConst, double DConst, int MyXOffset, int MyYOffset, int ItList[9][26], int NeighborX[26], int NeighborY[26], int NeighborZ[26], ViewI CritTimeStep, ViewF UndercoolingCurrent, ViewF UndercoolingChange, float* GrainUnitVector, ViewF CritDiagonalLength, ViewF DiagonalLength, int* GrainOrientation, ViewI CellType, ViewF DOCenter, ViewI GrainID, int NGrainOrientations, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferC, Buffer2D BufferD, Buffer2D BufferE, Buffer2D BufferF, Buffer2D BufferG, Buffer2D BufferH, int BufSizeX, int BufSizeY, ViewI Locks) {
    
    // Cell capture - parallel reduce loop over all type Active cells, counting number of ghost node cells that need to be accounted for
    Kokkos::parallel_for ("CellCapture",LocalDomainSize, KOKKOS_LAMBDA (const long int& D3D1ConvPosition) {
        
        // Test
//        if ((UndercoolingCurrent(D3D1ConvPosition) < 0)&&(cycle > CritTimeStep(D3D1ConvPosition)+1)&&(CellType(D3D1ConvPosition) != Wall)) {
//            printf("Superheated cell %d %d %d %f \n",id,D3D1ConvPosition,CellType(D3D1ConvPosition),UndercoolingCurrent(D3D1ConvPosition));
//        }
        // Cells of interest for the CA
        if ((CellType(D3D1ConvPosition) != Solid)&&(cycle > CritTimeStep(D3D1ConvPosition))) {

            if (CellType(D3D1ConvPosition) == Delayed) {
                // This delayed liquid cell is at the liquidus and ready for potential solidification
                //if (cycle == 2) printf("Liquid cell %d %d \n",id,D3D1ConvPosition);
                CellType(D3D1ConvPosition) = Liquid;
                UndercoolingCurrent(D3D1ConvPosition) = 0.0;
            }
            else if ((CellType(D3D1ConvPosition) == LiqSol)||(CellType(D3D1ConvPosition) == Liquid)) {
                // Update local undercooling - linear cooling between solidus and liquidus
                UndercoolingCurrent(D3D1ConvPosition) += UndercoolingChange(D3D1ConvPosition);
            }
            else if (CellType(D3D1ConvPosition) == Active) {

                int RankZ = floor(D3D1ConvPosition/(MyXSlices*MyYSlices));
                int Rem = D3D1ConvPosition % (MyXSlices*MyYSlices);
                int RankX = floor(Rem/MyYSlices);
                int RankY = Rem % MyYSlices;
                
                // Update local undercooling - linear cooling between solidus and liquidus
                UndercoolingCurrent(D3D1ConvPosition) += UndercoolingChange(D3D1ConvPosition);
                // Update local diagonal length of active cell
                double LocU = UndercoolingCurrent(D3D1ConvPosition);
                LocU = min(210.0,LocU);
                double V = AConst*pow(LocU,3) + BConst*pow(LocU,2) + CConst*LocU + DConst;
                V = max(0.0,V);
                DiagonalLength(D3D1ConvPosition) += min(0.045,V); // Max amount the diagonal can grow per time step
                //if (cycle >= 20000) printf("Active cell rank %d with Undercooling %f and UC %f and DL %f \n",id,LocU,UndercoolingChange(D3D1ConvPosition),V);
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
                
                
                //printf("Rank %d cell %d %d %d \n",id,RankX,RankY,RankZ);
                // "ll" corresponds to the specific position on the list of neighboring cells
                for (int ll=0; ll<NListLength; ll++) {
                    // "l" correpsponds to the specific neighboring cell
                    int l = ItList[ItBounds][ll];
                    // Local coordinates of adjacent cell center
                    int MyNeighborX = RankX + NeighborX[l];
                    int MyNeighborY = RankY + NeighborY[l];
                    int MyNeighborZ = RankZ + NeighborZ[l];
                    long int NeighborD3D1ConvPosition = MyNeighborZ*MyXSlices*MyYSlices + MyNeighborX*MyYSlices + MyNeighborY;
                    //                printf("Old cell at %d %d %d , CDL for direction %d is %f : \n", RankX,RankY,RankZ,ll,CritDiagonalLength(26*D3D1ConvPosition+l));
                    if ((CellType(NeighborD3D1ConvPosition) == Liquid)||(CellType(NeighborD3D1ConvPosition) == LiqSol)||(CellType(NeighborD3D1ConvPosition) == Delayed)) LCount = 1;
                    // Capture of cell located at "NeighborD3D1ConvPosition" if this condition is satisfied
                    if ((DiagonalLength(D3D1ConvPosition) >= CritDiagonalLength(26*D3D1ConvPosition+l))&&(Locks(NeighborD3D1ConvPosition) == 1)) {
                        // Use of atomic_compare_exchange (https://github.com/kokkos/kokkos/wiki/Kokkos%3A%3Aatomic_compare_exchange)
                        // old_val = atomic_compare_exchange(ptr_to_value,comparison_value, new_value);
                        // Atomicly sets the value at the address given by ptr_to_value to new_value if the current value at ptr_to_value is equal to comparison_value
                        // Returns the previously stored value at the address independent on whether the exchange has happened.
                        
                        // If this cell's value for "Locks" is still 1, replace it with 0 and return a value of 1
                        // If this cell's value for "Locks" has been changed to 0, return a value of 0
                        int OldLocksValue = Kokkos::atomic_compare_exchange(&Locks(NeighborD3D1ConvPosition),1,0);
                        
                        // If OldLocksValue is 0, this capture event already happened
                        // Only proceed if OldLocksValue is 1
                        if (OldLocksValue == 1) {
                            //printf("Old cell was at %d %d %d , CDL was %f and DL was %f :New cell at %d %d %d captured \n",RankX,RankY,RankZ,CritDiagonalLength(26*D3D1ConvPosition+l),DiagonalLength(D3D1ConvPosition),MyNeighborX,MyNeighborY,MyNeighborZ);
                            int GlobalX = RankX + MyXOffset;
                            int GlobalY = RankY + MyYOffset;
                            int h = GrainID(D3D1ConvPosition);
                            int MyOrientation = GrainOrientation[((abs(h) - 1) % NGrainOrientations)];
                            
                            // The new cell is captured by this cell's growing octahedron (Grain "h")
                            GrainID(NeighborD3D1ConvPosition) = h;
                            // (cxold, cyold, czold) are the coordiantes of this decentered octahedron
                            double cxold = DOCenter((long int)(3)*D3D1ConvPosition);
                            double cyold = DOCenter((long int)(3)*D3D1ConvPosition+(long int)(1));
                            double czold = DOCenter((long int)(3)*D3D1ConvPosition+(long int)(2));
                            
                            // (xp,yp,zp) are the global coordinates of the new cell's center
                            double xp = GlobalX + NeighborX[l] + 0.5;
                            double yp = GlobalY + NeighborY[l] + 0.5;
                            double zp = RankZ + NeighborZ[l] + 0.5;
                            
                            // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the image of the center of the new cell
                            double x0 = xp - cxold;
                            double y0 = yp - cyold;
                            double z0 = zp - czold;
                            
                            // mag0 is the magnitude of (x0,y0,z0)
                            double mag0 = pow(pow(x0,2) + pow(y0,2) + pow(z0,2),0.5);
                            
                            // Calculate unit vectors for the octahedron that intersect the new cell center
                            double Diag1X, Diag1Y, Diag1Z, Diag2X, Diag2Y, Diag2Z, Diag3X, Diag3Y, Diag3Z;
                            double Angle1 = (GrainUnitVector[9*MyOrientation]*x0 + GrainUnitVector[9*MyOrientation + 1]*y0 + GrainUnitVector[9*MyOrientation + 2]*z0)/mag0;
                            if (Angle1 < 0) {
                                Diag1X = GrainUnitVector[9*MyOrientation];
                                Diag1Y = GrainUnitVector[9*MyOrientation + 1];
                                Diag1Z = GrainUnitVector[9*MyOrientation + 2];
                            }
                            else {
                                Diag1X = -GrainUnitVector[9*MyOrientation];
                                Diag1Y = -GrainUnitVector[9*MyOrientation + 1];
                                Diag1Z = -GrainUnitVector[9*MyOrientation + 2];
                            }
                            
                            double Angle2 = (GrainUnitVector[9*MyOrientation + 3]*x0 + GrainUnitVector[9*MyOrientation + 4]*y0 + GrainUnitVector[9*MyOrientation + 5]*z0)/mag0;
                            if (Angle2 < 0) {
                                Diag2X = GrainUnitVector[9*MyOrientation + 3];
                                Diag2Y = GrainUnitVector[9*MyOrientation + 4];
                                Diag2Z = GrainUnitVector[9*MyOrientation + 5];
                            }
                            else {
                                Diag2X = -GrainUnitVector[9*MyOrientation + 3];
                                Diag2Y = -GrainUnitVector[9*MyOrientation + 4];
                                Diag2Z = -GrainUnitVector[9*MyOrientation + 5];
                            }
                            
                            double Angle3 = (GrainUnitVector[9*MyOrientation + 6]*x0 + GrainUnitVector[9*MyOrientation + 7]*y0 + GrainUnitVector[9*MyOrientation + 8]*z0)/mag0;
                            if (Angle3 < 0) {
                                Diag3X = GrainUnitVector[9*MyOrientation + 6];
                                Diag3Y = GrainUnitVector[9*MyOrientation + 7];
                                Diag3Z = GrainUnitVector[9*MyOrientation + 8];
                            }
                            else {
                                Diag3X = -GrainUnitVector[9*MyOrientation + 6];
                                Diag3Y = -GrainUnitVector[9*MyOrientation + 7];
                                Diag3Z = -GrainUnitVector[9*MyOrientation + 8];
                            }
                            
                            // Calculate angles between the octahedron diagonal directions and the vector x0,y0,z0
//                            double AnglesA[6];
//                            int index1, index2, index3;
//                            for (int aa=0; aa<6; aa++) {
//                                double xd = GrainUnitVector[18*MyOrientation + 3*aa];
//                                double yd = GrainUnitVector[18*MyOrientation + 3*aa + 1];
//                                double zd = GrainUnitVector[18*MyOrientation + 3*aa + 2];
//                                AnglesA[aa] = (xd*x0 + yd*y0 + zd*z0)/mag0;
//                            }
//                            index1 = 0;
//                            for (int ii=1; ii<6; ii++) {
//                                if (AnglesA[index1] < AnglesA[ii]) {
//                                    index1 = ii;
//                                }
//                            }
//                            AnglesA[index1] = -1;
//                            if (index1 % 2 == 0) AnglesA[index1+1] = -1;
//                            if (index1 % 2 == 1) AnglesA[index1-1] = -1;
//
//                            double MaxA = AnglesA[0];
//                            for (int ii=1; ii<6; ii++) {
//                                if (MaxA < AnglesA[ii]) {
//                                    MaxA = AnglesA[ii];
//                                }
//                            }
//                            if (MaxA == 0) {
//                                // Special case- other diagonals are all perpendicular to the first one (e.g. the octahedron corner captures the new cell center)
//                                // manually assign other diagonals (one of the 4 possible "capturing" faces)
//                                if ((index1 == 0)||(index1 == 1)) {
//                                    index2 = 2;
//                                    index3 = 4;
//                                }
//                                else if ((index1 == 2)||(index1 == 3)) {
//                                    index2 = 0;
//                                    index3 = 4;
//                                }
//                                else if ((index1 == 4)||(index1 == 5)) {
//                                    index2 = 0;
//                                    index3 = 2;
//                                }
//                            }
//                            else {
//                                // 2nd and 3rd closest diagonals to x0,y0,z0
//                                // note that if these are the same length, this means that the octahedron edge captures the new cell center
//                                // in this case, either of the 2 possible "capturing" faces will work
//                                index2 = 0;
//                                for (int ii=1; ii<6; ii++) {
//                                    if (AnglesA[index2] < AnglesA[ii]) {
//                                        index2 = ii;
//                                    }
//                                }
//                                AnglesA[index2] = -1;
//                                if (index2 % 2 == 0) AnglesA[index2+1] = -1;
//                                if (index2 % 2 == 1) AnglesA[index2-1] = -1;
//                                // index3 = MaxIndex(AnglesA);
//                                index3 = 0;
//                                for (int ii=1; ii<6; ii++) {
//                                    if (AnglesA[index3] < AnglesA[ii]) {
//                                        index3 = ii;
//                                    }
//                                }
//
//                            }
//                            double Diag1X = GrainUnitVector[18*MyOrientation + 3*index1];
//                            double Diag1Y = GrainUnitVector[18*MyOrientation + 3*index1 + 1];
//                            double Diag1Z = GrainUnitVector[18*MyOrientation + 3*index1 + 2];
//
//                            double Diag2X = GrainUnitVector[18*MyOrientation + 3*index2];
//                            double Diag2Y = GrainUnitVector[18*MyOrientation + 3*index2 + 1];
//                            double Diag2Z = GrainUnitVector[18*MyOrientation + 3*index2 + 2];
//                            double Diag3X = GrainUnitVector[18*MyOrientation + 3*index3];
//                            double Diag3Y = GrainUnitVector[18*MyOrientation + 3*index3 + 1];
//                            double Diag3Z = GrainUnitVector[18*MyOrientation + 3*index3 + 2];
                            
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
                            
                            DOCenter((long int)(3)*NeighborD3D1ConvPosition) = cx;
                            DOCenter((long int)(3)*NeighborD3D1ConvPosition+(long int)(1)) = cy;
                            DOCenter((long int)(3)*NeighborD3D1ConvPosition+(long int)(2)) = cz;
                            // Calculate critical diagonal lengths for the new active cell located at (xp,yp,zp) on the local grid
                            // For each neighbor (l=0 to 25), calculate which octahedron face leads to cell capture
                            // Calculate critical octahedron diagonal length to activate each nearest neighbor, as well as the coordinates of the triangle vertices on the capturing face
                            for (int n=0; n<26; n++)  {
                                int NeighborOfNeighborX = MyNeighborX + NeighborX[n];
                                int NeighborOfNeighborY = MyNeighborY + NeighborY[n];
                                int NeighborOfNeighborZ = MyNeighborZ + NeighborZ[n];
                                long int NeighborOfNeighborPosition = NeighborOfNeighborZ*MyXSlices*MyYSlices + NeighborOfNeighborX*MyYSlices + NeighborOfNeighborY;
                                if (NeighborOfNeighborPosition == D3D1ConvPosition) {
                                    // Do not calculate critical diagonal length req'd for the newly captured cell to capture the original
                                    //TriangleIndex(78*NeighborD3D1ConvPosition + 3*n) = 6;
                                    //TriangleIndex(78*NeighborD3D1ConvPosition + 3*n + 1) = 6;
                                    //TriangleIndex(78*NeighborD3D1ConvPosition + 3*n + 2) = 6;
                                    CritDiagonalLength((long int)(26)*NeighborD3D1ConvPosition+(long int)(n)) = 10000000.0;
                                }
                                else {
                                    // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the image of the center of a neighbor cell
                                    double x0 = xp + NeighborX[n] - cx;
                                    double y0 = yp + NeighborY[n] - cy;
                                    double z0 = zp + NeighborZ[n] - cz;
                                    // mag0 is the magnitude of (x0,y0,z0)
                                    double mag0 = pow(pow(x0,2) + pow(y0,2) + pow(z0,2),0.5);
                                    
                                    // Calculate unit vectors for the octahedron that intersect the new cell center
                                    double Diag1X, Diag1Y, Diag1Z, Diag2X, Diag2Y, Diag2Z, Diag3X, Diag3Y, Diag3Z;
                                    double Angle1 = (GrainUnitVector[9*MyOrientation]*x0 + GrainUnitVector[9*MyOrientation + 1]*y0 + GrainUnitVector[9*MyOrientation + 2]*z0)/mag0;
                                    if (Angle1 < 0) {
                                        Diag1X = GrainUnitVector[9*MyOrientation];
                                        Diag1Y = GrainUnitVector[9*MyOrientation + 1];
                                        Diag1Z = GrainUnitVector[9*MyOrientation + 2];
                                    }
                                    else {
                                        Diag1X = -GrainUnitVector[9*MyOrientation];
                                        Diag1Y = -GrainUnitVector[9*MyOrientation + 1];
                                        Diag1Z = -GrainUnitVector[9*MyOrientation + 2];
                                    }

                                    double Angle2 = (GrainUnitVector[9*MyOrientation + 3]*x0 + GrainUnitVector[9*MyOrientation + 4]*y0 + GrainUnitVector[9*MyOrientation + 5]*z0)/mag0;
                                    if (Angle2 < 0) {
                                        Diag2X = GrainUnitVector[9*MyOrientation + 3];
                                        Diag2Y = GrainUnitVector[9*MyOrientation + 4];
                                        Diag2Z = GrainUnitVector[9*MyOrientation + 5];
                                    }
                                    else {
                                        Diag2X = -GrainUnitVector[9*MyOrientation + 3];
                                        Diag2Y = -GrainUnitVector[9*MyOrientation + 4];
                                        Diag2Z = -GrainUnitVector[9*MyOrientation + 5];
                                    }

                                    double Angle3 = (GrainUnitVector[9*MyOrientation + 6]*x0 + GrainUnitVector[9*MyOrientation + 7]*y0 + GrainUnitVector[9*MyOrientation + 8]*z0)/mag0;
                                    if (Angle3 < 0) {
                                        Diag3X = GrainUnitVector[9*MyOrientation + 6];
                                        Diag3Y = GrainUnitVector[9*MyOrientation + 7];
                                        Diag3Z = GrainUnitVector[9*MyOrientation + 8];
                                    }
                                    else {
                                        Diag3X = -GrainUnitVector[9*MyOrientation + 6];
                                        Diag3Y = -GrainUnitVector[9*MyOrientation + 7];
                                        Diag3Z = -GrainUnitVector[9*MyOrientation + 8];
                                    }
                                    // Calculate angles between the octahedron diagonal directions and the vector x0,y0,z0

//                                    double AnglesA[6];
//                                    for (int aa=0; aa<6; aa++) {
//                                        double xd = GrainUnitVector[18*MyOrientation + 3*aa];
//                                        double yd = GrainUnitVector[18*MyOrientation + 3*aa + 1];
//                                        double zd = GrainUnitVector[18*MyOrientation + 3*aa + 2];
//                                        AnglesA[aa] = (xd*x0 + yd*y0 + zd*z0)/mag0;
//                                    }
//
//                                    int index1N, index2N, index3N;
//                                    index1N = 0;
//                                    for (int ii=1; ii<6; ii++) {
//                                        if (AnglesA[index1N] < AnglesA[ii]) {
//                                            index1N = ii;
//                                        }
//                                    }
//
//                                    //TriangleIndex(78*NeighborD3D1ConvPosition + 3*n) = index1N;
//
//                                    // First diagonal of the capturing face is that which makes the smallest (?) angle with x0,y0,z0
//                                    double Diag1X = GrainUnitVector[18*MyOrientation + 3*index1N + 0];
//                                    double Diag1Y = GrainUnitVector[18*MyOrientation + 3*index1N + 1];
//                                    double Diag1Z = GrainUnitVector[18*MyOrientation + 3*index1N + 2];
//                                    AnglesA[index1N] = -1;
//                                    if (index1N % 2 == 0) AnglesA[index1N+1] = -1;
//                                    if (index1N % 2 == 1) AnglesA[index1N-1] = -1;
//
//                                    double MaxA = AnglesA[0];
//                                    for (int ii=1; ii<6; ii++) {
//                                        if (MaxA < AnglesA[ii]) {
//                                            MaxA = AnglesA[ii];
//                                        }
//                                    }
//
//
//                                    double Diag2X, Diag2Y, Diag2Z, Diag3X, Diag3Y, Diag3Z;
//                                    if (MaxA == 0) {
//                                        // Special case- other diagonals are all perpendicular to the first one (e.g. the octahedron corner captures the new cell center)
//                                        // manually assign other diagonals (one of the 4 possible "capturing" faces)
//                                        if ((index1N == 0)||(index1N == 1)) {
//                                            index2 = 2;
//                                            index3 = 4;
//                                        }
//                                        else if ((index1N == 2)||(index1N == 3)) {
//                                            index2 = 0;
//                                            index3 = 4;
//                                        }
//                                        else if ((index1N == 4)||(index1N == 5)) {
//                                            index2 = 0;
//                                            index3 = 2;
//                                        }
//                                    }
//                                    else {
//                                        // 2nd and 3rd closest diagonals to x0,y0,z0
//                                        // note that if these are the same length, this means that the octahedron edge captures the new cell center
//                                        // in this case, either of the 2 possible "capturing" faces will work
//                                        index2N = 0;
//                                        for (int ii=1; ii<6; ii++) {
//                                            if (AnglesA[index2N] < AnglesA[ii]) {
//                                                index2N = ii;
//                                            }
//                                        }
//                                        AnglesA[index2N] = -1;
//                                        if (index2N % 2 == 0) AnglesA[index2N+1] = -1;
//                                        if (index2N % 2 == 1) AnglesA[index2N-1] = -1;
//                                        // index3 = MaxIndex(AnglesA);
//                                        index3N = 0;
//                                        for (int ii=1; ii<6; ii++) {
//                                            if (AnglesA[index3N] < AnglesA[ii]) {
//                                                index3N = ii;
//                                            }
//                                        }
//
//                                    }
//                                    //TriangleIndex(78*NeighborD3D1ConvPosition + 3*n + 1) = index2N;
//                                    Diag2X = GrainUnitVector[18*MyOrientation + 3*index2N + 0];
//                                    Diag2Y = GrainUnitVector[18*MyOrientation + 3*index2N + 1];
//                                    Diag2Z = GrainUnitVector[18*MyOrientation + 3*index2N + 2];
//                                    //TriangleIndex(78*NeighborD3D1ConvPosition + 3*n + 2) = index3N;
//                                    Diag3X = GrainUnitVector[18*MyOrientation + 3*index3N + 0];
//                                    Diag3Y = GrainUnitVector[18*MyOrientation + 3*index3N + 1];
//                                    Diag3Z = GrainUnitVector[18*MyOrientation + 3*index3N + 2];
//                                    
                                    
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
                                    CritDiagonalLength((long int)(26)*NeighborD3D1ConvPosition+(long int)(n)) = CDLVal;
                                    //                                if (CDLVal == 0.0) printf("Zero CDLVal : %d %d %d %d %f %d %d %d %f %f %f",MyNeighborX,MyNeighborY,MyNeighborZ,n,mag0,index1,index2,index3,normx,normy,normz);
                                }
                            }
                            
                            if (np > 0) {
                                
                                float GhostGID = (float)(GrainID(NeighborD3D1ConvPosition));
                                float GhostDOCX = DOCenter((long int)(3)*NeighborD3D1ConvPosition);
                                float GhostDOCY = DOCenter((long int)(3)*NeighborD3D1ConvPosition+(long int)(1));
                                float GhostDOCZ = DOCenter((long int)(3)*NeighborD3D1ConvPosition+(long int)(2));
                                float GhostDL = DiagonalLength(NeighborD3D1ConvPosition);
                                // Collect data for the ghost nodes:
                                if (DecompositionStrategy == 1) {
                                    int GNPosition = MyNeighborZ*BufSizeX + MyNeighborX;
                                    if (MyNeighborY == 1) {
                                        BufferA(GNPosition,0) = GhostGID;
                                        BufferA(GNPosition,1) = GhostDOCX;
                                        BufferA(GNPosition,2) = GhostDOCY;
                                        BufferA(GNPosition,3) = GhostDOCZ;
                                        BufferA(GNPosition,4) = GhostDL;
                                        //printf("NewDL for A is %f \n",BufferA(GNPosition,4));
                                    }
                                    else if (MyNeighborY == MyYSlices-2) {
                                        int GNPosition = MyNeighborZ*BufSizeX + MyNeighborX;
                                        BufferB(GNPosition,0) = GhostGID;
                                        BufferB(GNPosition,1) = GhostDOCX;
                                        BufferB(GNPosition,2) = GhostDOCY;
                                        BufferB(GNPosition,3) = GhostDOCZ;
                                        BufferB(GNPosition,4) = GhostDL;
                                        //printf("NewDL for B is %f \n",BufferB(GNPosition,4));
                                    }
                                }
                                else {
                                    if (MyNeighborY == 1) {
                                        // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                                        if (MyNeighborX == MyXSlices-2) {
                                            // To MyLeft (BufferA)
                                            int GNPosition = MyNeighborZ*BufSizeX + MyNeighborX-1;
                                            BufferA(GNPosition,0) = GhostGID;
                                            BufferA(GNPosition,1) = GhostDOCX;
                                            BufferA(GNPosition,2) = GhostDOCY;
                                            BufferA(GNPosition,3) = GhostDOCZ;
                                            BufferA(GNPosition,4) = GhostDL;
                                            // To MyOut (BufferC)
                                            GNPosition = MyNeighborZ*BufSizeY + MyNeighborY-1;
                                            BufferC(GNPosition,0) = GhostGID;
                                            BufferC(GNPosition,1) = GhostDOCX;
                                            BufferC(GNPosition,2) = GhostDOCY;
                                            BufferC(GNPosition,3) = GhostDOCZ;
                                            BufferC(GNPosition,4) = GhostDL;
                                            // To MyLeftOut (BufferE)
                                            GNPosition = MyNeighborZ;
                                            BufferE(GNPosition,0) = GhostGID;
                                            BufferE(GNPosition,1) = GhostDOCX;
                                            BufferE(GNPosition,2) = GhostDOCY;
                                            BufferE(GNPosition,3) = GhostDOCZ;
                                            BufferE(GNPosition,4) = GhostDL;
                                        }
                                        else if (MyNeighborX == 1) {
                                            int GNPosition = MyNeighborZ*BufSizeX + MyNeighborX-1;
                                            BufferA(GNPosition,0) = GhostGID;
                                            BufferA(GNPosition,1) = GhostDOCX;
                                            BufferA(GNPosition,2) = GhostDOCY;
                                            BufferA(GNPosition,3) = GhostDOCZ;
                                            BufferA(GNPosition,4) = GhostDL;
                                            GNPosition = MyNeighborZ*BufSizeY + MyNeighborY-1;
                                            BufferD(GNPosition,0) = GhostGID;
                                            BufferD(GNPosition,1) = GhostDOCX;
                                            BufferD(GNPosition,2) = GhostDOCY;
                                            BufferD(GNPosition,3) = GhostDOCZ;
                                            BufferD(GNPosition,4) = GhostDL;
                                            GNPosition = MyNeighborZ;
                                            BufferG(GNPosition,0) = GhostGID;
                                            BufferG(GNPosition,1) = GhostDOCX;
                                            BufferG(GNPosition,2) = GhostDOCY;
                                            BufferG(GNPosition,3) = GhostDOCZ;
                                            BufferG(GNPosition,4) = GhostDL;
                                        }
                                        else if ((MyNeighborX > 1)&&(MyNeighborX < MyXSlices-2)) {
                                            // This is being sent to MyLeft
                                            int GNPosition = MyNeighborZ*BufSizeX + MyNeighborX-1;
                                            BufferA(GNPosition,0) = GhostGID;
                                            BufferA(GNPosition,1) = GhostDOCX;
                                            BufferA(GNPosition,2) = GhostDOCY;
                                            BufferA(GNPosition,3) = GhostDOCZ;
                                            BufferA(GNPosition,4) = GhostDL;
                                        }
                                    }
                                    else if (MyNeighborY == MyYSlices-2) {
                                        // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                                        if (MyNeighborX == MyXSlices-2) {
                                            int GNPosition = MyNeighborZ*BufSizeX + MyNeighborX-1;
                                            BufferB(GNPosition,0) = GhostGID;
                                            BufferB(GNPosition,1) = GhostDOCX;
                                            BufferB(GNPosition,2) = GhostDOCY;
                                            BufferB(GNPosition,3) = GhostDOCZ;
                                            BufferB(GNPosition,4) = GhostDL;
                                            GNPosition = MyNeighborZ*BufSizeY + MyNeighborY-1;
                                            BufferC(GNPosition,0) = GhostGID;
                                            BufferC(GNPosition,1) = GhostDOCX;
                                            BufferC(GNPosition,2) = GhostDOCY;
                                            BufferC(GNPosition,3) = GhostDOCZ;
                                            BufferC(GNPosition,4) = GhostDL;
                                            GNPosition = MyNeighborZ;
                                            BufferF(GNPosition,0) = GhostGID;
                                            BufferF(GNPosition,1) = GhostDOCX;
                                            BufferF(GNPosition,2) = GhostDOCY;
                                            BufferF(GNPosition,3) = GhostDOCZ;
                                            BufferF(GNPosition,4) = GhostDL;
                                        }
                                        else if (MyNeighborX == 1) {
                                            int GNPosition = MyNeighborZ*BufSizeX + MyNeighborX-1;
                                            BufferB(GNPosition,0) = GhostGID;
                                            BufferB(GNPosition,1) = GhostDOCX;
                                            BufferB(GNPosition,2) = GhostDOCY;
                                            BufferB(GNPosition,3) = GhostDOCZ;
                                            BufferB(GNPosition,4) = GhostDL;
                                            GNPosition = MyNeighborZ*BufSizeY + MyNeighborY-1;
                                            BufferD(GNPosition,0) = GhostGID;
                                            BufferD(GNPosition,1) = GhostDOCX;
                                            BufferD(GNPosition,2) = GhostDOCY;
                                            BufferD(GNPosition,3) = GhostDOCZ;
                                            BufferD(GNPosition,4) = GhostDL;
                                            GNPosition = MyNeighborZ;
                                            BufferH(GNPosition,0) = GhostGID;
                                            BufferH(GNPosition,1) = GhostDOCX;
                                            BufferH(GNPosition,2) = GhostDOCY;
                                            BufferH(GNPosition,3) = GhostDOCZ;
                                            BufferH(GNPosition,4) = GhostDL;
                                        }
                                        else if ((MyNeighborX > 1)&&(MyNeighborX < MyXSlices-2)) {
                                            int GNPosition = MyNeighborZ*BufSizeX + MyNeighborX-1;
                                            BufferB(GNPosition,0) = GhostGID;
                                            BufferB(GNPosition,1) = GhostDOCX;
                                            BufferB(GNPosition,2) = GhostDOCY;
                                            BufferB(GNPosition,3) = GhostDOCZ;
                                            BufferB(GNPosition,4) = GhostDL;
                                        }
                                    }
                                    else if ((MyNeighborX == 1)&&(MyNeighborY > 1)&&(MyNeighborY < MyYSlices-2)) {
                                        int GNPosition = MyNeighborZ*BufSizeY + MyNeighborY-1;
                                        BufferD(GNPosition,0) = GhostGID;
                                        BufferD(GNPosition,1) = GhostDOCX;
                                        BufferD(GNPosition,2) = GhostDOCY;
                                        BufferD(GNPosition,3) = GhostDOCZ;
                                        BufferD(GNPosition,4) = GhostDL;
                                        //if (id == 0) cout << "RANK 0 LISTED " << MyNeighborX << " " << MyNeighborY << " " << MyNeighborZ << endl;
                                    }
                                    else if ((MyNeighborX == MyXSlices-2)&&(MyNeighborY > 1)&&(MyNeighborY < MyYSlices-2)) {
                                        int GNPosition = MyNeighborZ*BufSizeY + MyNeighborY-1;
                                        BufferC(GNPosition,0) = GhostGID;
                                        BufferC(GNPosition,1) = GhostDOCX;
                                        BufferC(GNPosition,2) = GhostDOCY;
                                        BufferC(GNPosition,3) = GhostDOCZ;
                                        BufferC(GNPosition,4) = GhostDL;
                                    }
                                } // End if statement for ghost node marking
                            } // End if statement for serial/parallel code
                            // Only update the new cell's type once Critical Diagonal Length, Triangle Index, and Diagonal Length values have been assigned to it
                            // Avoids the race condition in which the new cell is activated, and another thread acts on the new active cell before
                            // the cell's new critical diagonal length/triangle index/diagonal length values are assigned
                            CellType(NeighborD3D1ConvPosition) = Active;
                        } // End if statement within locked capture loop
                    } // End if statement for outer capture loop
                } // End loop over all neighbors of this active cell
                if (LCount == 0) {
                    // This active cell has no more neighboring cells to be captured, becomes solid
                    CellType(D3D1ConvPosition) = Solid;
                }
            } // end "if" loop for Active type cells
        } // end "if" loop over cells of interest
    });

    Kokkos::fence();
    // Fix corrupted lock values
    Kokkos::parallel_for ("CellCapture",LocalDomainSize, KOKKOS_LAMBDA (const long int& D3D1ConvPosition) {
        if ((CellType(D3D1ConvPosition) == Liquid)||(CellType(D3D1ConvPosition) == LiqSol)||(CellType(D3D1ConvPosition) == Delayed)) {
            Kokkos::atomic_compare_exchange(&Locks(D3D1ConvPosition),0,1);
        }
    });
    
}

// Prints intermediate code output to stdout, checks to see if solidification is complete
void IntermediateOutputAndCheck(int id, int &cycle, int LocalDomainSize, int nn, int &XSwitch, ViewI CellType, ViewI CritTimeStep, string TemperatureDataType) {
    
    sample::ValueType CellTypeStorage;
    Kokkos::parallel_reduce (  LocalDomainSize, KOKKOS_LAMBDA (const int& D3D1ConvPosition, sample::ValueType & upd) {
            if (CellType(D3D1ConvPosition) == Delayed) upd.the_array[0] += 1;
            else if (CellType(D3D1ConvPosition) == Liquid)  upd.the_array[1] += 1;
            else if (CellType(D3D1ConvPosition) == Active)     upd.the_array[2] += 1;
            else if (CellType(D3D1ConvPosition) == Solid)      upd.the_array[3] += 1;
    }, Kokkos::Sum<sample::ValueType>(CellTypeStorage) );
    signed long int Global_nn = 0;
    signed long int Type_sumD = CellTypeStorage.the_array[0];
    signed long int Type_sumL = CellTypeStorage.the_array[1];
    signed long int Type_sumA = CellTypeStorage.the_array[2];
    signed long int Type_sumS = CellTypeStorage.the_array[3];

    signed long int Global_sumL, Global_sumD, Global_sumA, Global_sumS;
    MPI_Reduce(&Type_sumD,&Global_sumD,1,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&Type_sumL,&Global_sumL,1,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&Type_sumA,&Global_sumA,1,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&Type_sumS,&Global_sumS,1,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&nn,&Global_nn,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    
    if (id == 0) {
        cout << "======================================================" << endl;
        cout << "cycle = " << cycle << endl;
        cout << "Superheated liquid cells remaining = " << Global_sumD << endl;
        cout << "Undercooled liquid cells remaining = " << Global_sumL << endl;
        cout << "Active interface cells remaining = " << Global_sumA << endl;
        cout << "Solid cells remaining = " << Global_sumS << endl;
        cout << "Number of nucleation events = " << Global_nn << endl;
        cout << "======================================================" << endl;
        if (Global_sumL + Global_sumD == 0) {
            XSwitch = 1;
        }
    }
    MPI_Bcast(&XSwitch,1,MPI_INT,0,MPI_COMM_WORLD);
    if (XSwitch == 0) {
        if (TemperatureDataType == "R") {
            MPI_Bcast(&Global_sumL,1,MPI_LONG_LONG,0,MPI_COMM_WORLD);
            if (Global_sumL == 0) {
                // Check when the next superheated cells go below the liquidus
                int NextCTS;
                Kokkos::parallel_reduce("CellCapture",LocalDomainSize, KOKKOS_LAMBDA (const int& D3D1ConvPosition, int &tempv) {
                    if (CellType(D3D1ConvPosition) == Delayed) {
                        //if (CritTimeStep(D3D1ConvPosition) > 20000000) printf("CTS on rank %d is %d \n",id,CritTimeStep(D3D1ConvPosition));
                        if (CritTimeStep(D3D1ConvPosition) < tempv) tempv = CritTimeStep(D3D1ConvPosition);
                    }
                },Kokkos::Min<int>(NextCTS));
                //cout << "ID " << id << " NextCTS = " << NextCTS << endl;
                signed long int GlobalNextCTS;
                MPI_Allreduce(&NextCTS,&GlobalNextCTS,1,MPI_LONG_LONG,MPI_MIN,MPI_COMM_WORLD);
                if ((GlobalNextCTS - cycle) > 5000) {
                    // Jump to next time step when solidification starts again
                    cycle = GlobalNextCTS-1;
                    if (id == 0) cout << "Jumping to cycle " << cycle+1 << endl;
                }
            }
        }
    }
}
