#include "header.h"
using namespace std;

///*****************************************************************************/
/**
  Prints values of grain orientation for all cells to files
*/
void CollectGrainData(int id, int np, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI::HostMirror GrainID, ViewI::HostMirror GrainOrientation, ViewF::HostMirror GrainUnitVector, string BaseFileName, int DecompositionStrategy, int NGrainOrientations, bool* Melted, string PathToOutput, bool FilesToPrint[4], double deltax) {

    if (id == 0) {
        // Create GrainID variable for entire domain, place GrainIDs for rank 0
        vector <vector <vector <int> > > GrainID_WholeDomain;
        vector <vector <vector <bool> > > Melted_WholeDomain;
        for (int k=0; k<nz; k++) {
            vector <vector <int> > GrainID_JK;
            vector <vector <bool> > Melted_JK;
            for (int i=0; i<nx; i++) {
                vector <int> GrainID_K;
                vector <bool> Melted_K;
                for (int j=0; j<ny; j++) {
                    GrainID_K.push_back(0);
                    Melted_K.push_back(false);
                }
                GrainID_JK.push_back(GrainID_K);
                Melted_JK.push_back(Melted_K);
            }
            GrainID_WholeDomain.push_back(GrainID_JK);
            Melted_WholeDomain.push_back(Melted_JK);
        }
        for (int k=1; k<nz-1; k++) {
            for (int i=1; i<MyXSlices-1; i++) {
                for (int j=1; j<MyYSlices-1; j++) {
                    GrainID_WholeDomain[k][i-1][j-1] = GrainID(k*MyXSlices*MyYSlices + i*MyYSlices + j);
                    Melted_WholeDomain[k][i-1][j-1] = Melted[k*MyXSlices*MyYSlices+i*MyYSlices+j];
                }
            }
        }
        
//      cout << "RANK 0 PLACED AT X = 1 through " << MyXSlices-2 << " , Y = 1 through " << MyYSlices-2 << endl;
        // Recieve Grain ID value from other ranks
        // Message size different for different ranks
        for (int p=1; p<np; p++) {
           // cout << "FEEDING " << p << " " << nx << " " << ProcessorsInXDirection << " " << np << " " << DecompositionStrategy << endl;
            int RecvXOffset = XOffsetCalc(p,nx,ProcessorsInXDirection,ProcessorsInYDirection,np,DecompositionStrategy);
            int RecvXSlices = XMPSlicesCalc(p,nx,ProcessorsInXDirection,ProcessorsInYDirection,np,DecompositionStrategy);

            int RecvYOffset = YOffsetCalc(p,ny,ProcessorsInYDirection,np,DecompositionStrategy);
            int RecvYSlices = YMPSlicesCalc(p,ny,ProcessorsInYDirection,np,DecompositionStrategy);

            // No ghost nodes in send and recieved data
            //int XPosition = p/ProcessorsInYDirection;
            //int YPosition = p % ProcessorsInYDirection;
            RecvYSlices = RecvYSlices-2;
            RecvXSlices = RecvXSlices-2;
            
            // Readjust X and Y offsets of incoming data based on the lack of ghost nodes
            RecvYOffset++;
            RecvXOffset++;

            int RBufSize = (RecvXSlices)*(RecvYSlices)*(nz-2);
            int *RecvBufGID = new int[RBufSize];
            int *RecvBufM = new int[RBufSize];
            int DataCounter = 0;
            MPI_Recv(RecvBufGID,RBufSize,MPI_INT,p,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(RecvBufM,RBufSize,MPI_INT,p,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

            for (int k=1; k<nz-1; k++) {
                for (int i=0; i<RecvXSlices; i++) {
                    for (int j=0; j<RecvYSlices; j++) {
                        GrainID_WholeDomain[k][i+RecvXOffset][j+RecvYOffset] = RecvBufGID[DataCounter];
                        if (RecvBufM[DataCounter] == 1) Melted_WholeDomain[k][i+RecvXOffset][j+RecvYOffset] = true;
                        else Melted_WholeDomain[k][i+RecvXOffset][j+RecvYOffset] = false;
                        DataCounter++;
                    }
                }
            }
            //cout << "DataCounter = " << DataCounter << " of " << RBufSize << endl;
        }
        
        string FName = PathToOutput + BaseFileName;

        if (FilesToPrint[0]) PrintOrientations(FName, nx, ny, nz, Melted_WholeDomain, GrainID_WholeDomain, NGrainOrientations);
        if (FilesToPrint[1]) PrintGrainIDs(FName, nx, ny, nz, Melted_WholeDomain, GrainID_WholeDomain);
        if (FilesToPrint[2]) PrintGrainIDsForExaConstit(FName, nx, ny, nz, GrainID_WholeDomain,deltax);
        if (FilesToPrint[3]) PrintParaview(FName, nx, ny, nz, Melted_WholeDomain, GrainID_WholeDomain, GrainOrientation, GrainUnitVector, NGrainOrientations);
        if (FilesToPrint[4]) PrintGrainAreas(FName, deltax, nx, ny, nz, GrainID_WholeDomain);
        if (FilesToPrint[5]) PrintWeightedGrainAreas(FName, deltax, nx, ny, nz, GrainID_WholeDomain);
        
    }
    else {

        // Send non-ghost node data to rank 0
        //cout << "Rank " << id << endl;
        int SBufSize = (MyXSlices-2)*(MyYSlices-2)*(nz-2);
        int *SendBufGID = new int[SBufSize];
        int *SendBufM = new int [SBufSize];
        //if (id == 4) cout << "Cell type is " << GrainID((nz-2)*MyXSlices*MyYSlices + 233*MyYSlices + 1) << endl;
        int DataCounter = 0;

        for (int k=1; k<nz-1; k++) {
            for (int i=1; i<MyXSlices-1; i++) {
                for (int j=1; j<MyYSlices-1; j++) {
                    SendBufGID[DataCounter] = GrainID(k*MyXSlices*MyYSlices+i*MyYSlices+j);
                    if (Melted[k*MyXSlices*MyYSlices+i*MyYSlices+j]) SendBufM[DataCounter] = 1;
                    else SendBufM[DataCounter] = 0;
                    DataCounter++;
                }
            }
        }
        
        //cout << "Rank " << id << " sending " <<  SBufSize << " data points" << endl;
        MPI_Send(SendBufGID,SBufSize,MPI_INT,0,0,MPI_COMM_WORLD);
        MPI_Send(SendBufM,SBufSize,MPI_INT,0,1,MPI_COMM_WORLD);
    }

}

///*****************************************************************************/
/**
 Prints values of critical undercooling for all cells to files
 */
void PrintTempValues(int id, int np, int nx, int ny, int nz, int MyXSlices,int MyYSlices,int ProcessorsInXDirection,int ProcessorsInYDirection, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange, int DecompositionStrategy, string PathToOutput) {
    
    // Critical time step for solidification start printed to file first
    if (id == 0) {
        // Create GrainID variable for entire domain, place GrainIDs for rank 0
        vector <vector <vector <int> > > CritTimeStep_WholeDomain;
        for (int k=0; k<nz; k++) {
            vector <vector <int> > CTS_JK;
            for (int i=0; i<nx; i++) {
                vector <int> CTS_K;
                for (int j=0; j<ny; j++) {
                    CTS_K.push_back(0);
                }
                CTS_JK.push_back(CTS_K);
            }
            CritTimeStep_WholeDomain.push_back(CTS_JK);
        }
        //cout << "X" << endl;
        for (int k=1; k<nz-1; k++) {
            for (int i=1; i<MyXSlices-1; i++) {
                for (int j=1; j<MyYSlices-1; j++) {
                    CritTimeStep_WholeDomain[k][i-1][j-1] = CritTimeStep(k*MyXSlices*MyYSlices + i*MyYSlices + j);
                }
            }
        }
        //cout << "Y" << endl;
        // Recieve values from other ranks
        // Message size different for different ranks
        for (int p=1; p<np; p++) {
            //cout << "FEEDING " << p << " " << nx << " " << ProcessorsInXDirection << " " << np << " " << DecompositionStrategy << endl;
            int RecvXOffset = XOffsetCalc(p,nx,ProcessorsInXDirection,ProcessorsInYDirection,np,DecompositionStrategy);
            int RecvXSlices = XMPSlicesCalc(p,nx,ProcessorsInXDirection,ProcessorsInYDirection,np,DecompositionStrategy);
            
            int RecvYOffset = YOffsetCalc(p,ny,ProcessorsInYDirection,np,DecompositionStrategy);
            int RecvYSlices = YMPSlicesCalc(p,ny,ProcessorsInYDirection,np,DecompositionStrategy);
            
            // No ghost nodes in send and recieved data
            RecvYSlices = RecvYSlices-2;
            RecvXSlices = RecvXSlices-2;
            
            // Readjust X and Y offsets of incoming data based on the lack of ghost nodes
            RecvYOffset++;
            RecvXOffset++;

            int RBufSize = RecvXSlices*RecvYSlices*(nz-2);
            //cout << "RBufSize = " << RBufSize << endl;
            int *RecvBufGID = new int[RBufSize];
            //cout << "Rank 0 ready to recieve " << RBufSize << " int data from rank " << p << endl;
            MPI_Recv(RecvBufGID,RBufSize,MPI_INT,p,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            //cout << "Recieved data from rank " << p << endl;
            int DataCounter = 0;
            for (int k=1; k<nz-1; k++) {
                for (int i=0; i<RecvXSlices; i++) {
                    for (int j=0; j<RecvYSlices; j++) {
                        CritTimeStep_WholeDomain[k][i+RecvXOffset][j+RecvYOffset] = RecvBufGID[DataCounter];
                        DataCounter++;
                    }
                }
            }
            //cout << "DataCounter = " << DataCounter << " of " << RBufSize << endl;
            //cout << "Placed data in Rows " << RecvXOffset << " through " << RecvXOffset+RecvXSlices-1 << " and in cols " << RecvYOffset << " through " << RecvYOffset+RecvYSlices-1 << endl;
        }
        
        // Print grain orientations to file
        string FName = PathToOutput + "CriticalTimeStep.vtk";
        std::ofstream CTSplot;
//        // New vtk file
        CTSplot.open(FName);
        CTSplot << "# vtk DataFile Version 3.0" << endl;
        CTSplot << "vtk output" << endl;
        CTSplot << "ASCII" << endl;
        CTSplot << "DATASET STRUCTURED_POINTS" << endl;
        CTSplot << "DIMENSIONS " << nx-2 << " " << ny-2 << " " << nz-2 << endl;
        CTSplot << "ORIGIN 0 0 0" << endl;
        CTSplot << "SPACING 1 1 1" << endl;
        CTSplot << fixed << "POINT_DATA " << (nx-2)*(ny-2)*(nz-2) << endl;
        CTSplot << "SCALARS Angle_z float 1" << endl;
        CTSplot << "LOOKUP_TABLE default" << endl;
        
        for (int k=1; k<nz-1; k++) {
            for (int j=1; j<ny-1; j++) {
                for (int i=1; i<nx-1; i++) {
                    CTSplot << 0.75*pow(10,-6)*(float)(CritTimeStep_WholeDomain[k][i][j]) << " ";
                }
            }
            CTSplot << endl;
        }
        CTSplot.close();

        
    }
    else {
        // Send non-ghost node data to rank 0
        //cout << "Rank " << id << endl;
        int SBufSize = (MyXSlices-2)*(MyYSlices-2)*(nz-2);
        int *SendBufGID = new int[SBufSize];
        //cout << "Rank " << id << endl;
        int DataCounter = 0;
        // cout << "Rank " << id << " sending " << nx*nz*(MyYSlices-2) << " data points" << endl;
        for (int k=1; k<nz-1; k++) {
            for (int i=1; i<MyXSlices-1; i++) {
                for (int j=1; j<MyYSlices-1; j++) {
                    SendBufGID[DataCounter] = CritTimeStep(k*MyXSlices*MyYSlices + i*MyYSlices + j);
                    DataCounter++;
                }
            }
        }
        //cout << "Rank " << id << " sending " << SBufSize << " int data points" << endl;
        MPI_Send(SendBufGID,SBufSize,MPI_INT,0,0,MPI_COMM_WORLD);
        //cout << "RANK " << id << " SENT DATA" << endl;
    }

}

///*****************************************************************************/

void PrintCT(int id, int np, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI::HostMirror CellType, string BaseFileName, int DecompositionStrategy) {

    string FName = BaseFileName + "_CT.vtk";
    std::ofstream Grainplot;
    if (id == 0) {
        // Create GrainID variable for entire domain, place GrainIDs for rank 0
        vector <vector <vector <int> > > CellType_WholeDomain;
        for (int k=0; k<nz; k++) {
            vector <vector <int> > CellType_JK;
            for (int i=0; i<nx; i++) {
                vector <int> CellType_K;
                for (int j=0; j<ny; j++) {
                    CellType_K.push_back(Wall);
                }
                CellType_JK.push_back(CellType_K);
            }
            CellType_WholeDomain.push_back(CellType_JK);
        }
        for (int k=1; k<nz-1; k++) {
            for (int i=1; i<MyXSlices-1; i++) {
                for (int j=1; j<MyYSlices-1; j++) {
                    CellType_WholeDomain[k][i-1][j-1] = CellType(k*MyXSlices*MyYSlices + i*MyYSlices + j);
                }
            }
        }
        
        cout << "Domain created: Rank 0 placed slices 0 through" << MyYSlices-2 << endl;
        // Recieve Grain ID value from other ranks
        // Message size different for different ranks
        for (int p=1; p<np; p++) {
            // cout << "FEEDING " << p << " " << nx << " " << ProcessorsInXDirection << " " << np << " " << DecompositionStrategy << endl;
            int RecvXOffset = XOffsetCalc(p,nx,ProcessorsInXDirection,ProcessorsInYDirection,np,DecompositionStrategy);
            int RecvXSlices = XMPSlicesCalc(p,nx,ProcessorsInXDirection,ProcessorsInYDirection,np,DecompositionStrategy);
            
            int RecvYOffset = YOffsetCalc(p,ny,ProcessorsInYDirection,np,DecompositionStrategy);
            int RecvYSlices = YMPSlicesCalc(p,ny,ProcessorsInYDirection,np,DecompositionStrategy);
            
            // No ghost nodes in send and recieved data
//            int XPosition = p/ProcessorsInYDirection;
//            int YPosition = p % ProcessorsInYDirection;
            RecvYSlices = RecvYSlices-2;
            RecvXSlices = RecvXSlices-2;
            
            // Readjust X and Y offsets of incoming data based on the lack of ghost nodes
            RecvYOffset++;
            RecvXOffset++;
                //}
                        //cout << "Allocating memory" << endl;
            int RBufSize = RecvXSlices*RecvYSlices*(nz-2);
            int *RecvBufGID = new int[RBufSize];
            //cout << "Rank 0 recieving from " << p << " a total of " << RecvYSlices << " slices; Y from " << RecvYOffset << " to " << RecvYOffset+RecvYSlices-1 << " and " << RBufSize << " data points total" << endl;
            MPI_Recv(RecvBufGID,RBufSize,MPI_INT,p,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            //cout << "Placing data in Rows " << RecvXOffset << " through " << RecvXOffset+RecvXSlices-1 << " and in cols " << RecvYOffset << " through " << RecvYOffset+RecvYSlices-1 << endl;
            int DataCounter = 0;
            for (int k=1; k<nz-1; k++) {
                for (int i=0; i<RecvXSlices; i++) {
                    for (int j=0; j<RecvYSlices; j++) {
                        CellType_WholeDomain[k][i+RecvXOffset][j+RecvYOffset] = RecvBufGID[DataCounter];
                        DataCounter++;
                    }
                }
            }
            //cout << "DataCounter = " << DataCounter << " of " << RBufSize << endl;
        }
        cout << "Opening file" << endl;
        // Print grain orientations to file
        std::ofstream Grainplot;
        // New vtk file
        Grainplot.open(FName);
        Grainplot << "# vtk DataFile Version 3.0" << endl;
        Grainplot << "vtk output" << endl;
        Grainplot << "ASCII" << endl;
        Grainplot << "DATASET STRUCTURED_POINTS" << endl;
        Grainplot << "DIMENSIONS " << nx-2 << " " << ny-2 << " " << nz-2 << endl;
        Grainplot << "ORIGIN 0 0 0" << endl;
        Grainplot << "SPACING 1 1 1" << endl;
        Grainplot << fixed << "POINT_DATA " << (nx-2)*(ny-2)*(nz-2) << endl;
        Grainplot << "SCALARS Angle_z int 1" << endl;
        Grainplot << "LOOKUP_TABLE default" << endl;
        
        for (int k=1; k<nz-1; k++) {
            for (int j=1; j<ny-1; j++) {
                for (int i=1; i<nx-1; i++) {    
                    if (CellType_WholeDomain[k][i][j] == Solid) Grainplot << 0 << " ";
                    else if (CellType_WholeDomain[k][i][j] == Liquid) Grainplot << 2 << " ";
                    else if (CellType_WholeDomain[k][i][j] == Active) Grainplot << 1 << " ";
                    else  Grainplot << 3 << " ";
                }
            }
            Grainplot << endl;
        }
        Grainplot.close();
    }
    else {

        // Send non-ghost node data to rank 0
        //cout << "Rank " << id << endl;
        int SBufSize = (MyXSlices-2)*(MyYSlices-2)*(nz-2);
        int *SendBufGID = new int[SBufSize];
        //cout << "Rank " << id << " Sending " << MyYSlices-2 << " slices " << endl; //; Y from " << MyYOffset+1 << " to " << MyYOffset+MyYSlices-2 << endl;
        int DataCounter = 0;
        for (int k=1; k<nz-1; k++) {
            for (int i=1; i<MyXSlices-1; i++) {
                for (int j=1; j<MyYSlices-1; j++) {
                    SendBufGID[DataCounter] = CellType(k*MyXSlices*MyYSlices + i*MyYSlices + j);
                    DataCounter++;
                }
            }
        }
       // cout << "Rank " << id << " sending " << SBufSize << " data points" << endl;
        MPI_Send(SendBufGID,SBufSize,MPI_INT,0,0,MPI_COMM_WORLD);
       // cout << "RANK " << id << " SENT DATA" << endl;
    }
    
}

void PrintOrientations(string FName, int nx, int ny, int nz, vector <vector <vector <bool> > > Melted_WholeDomain, vector <vector <vector <int> > > GrainID_WholeDomain, int NGrainOrientations) {

    // Histogram of orientations for texture determination
    cout << "Printing file of grain orientation distribution, out of " << NGrainOrientations << " possible orientations for each grain" << endl;
    ofstream Grainplot0;
    string FName0 = FName + "_Orientations.csv";
    Grainplot0.open(FName0);
    int GOHistogram[NGrainOrientations];
    for (int i=0; i<NGrainOrientations; i++) {
        GOHistogram[i] = 0;
    }
    cout << "Histogram made" << endl;
    // frequency data on grain ids
    for (int k=1; k<nz-1; k++) {
        for (int j=1; j<ny-1; j++) {
            for (int i=1; i<nx-1; i++) {
                if (Melted_WholeDomain[k][i][j]) {
                    //cout << GrainID_WholeDomain[k][i][j] << endl;
                    int GOVal = (abs(GrainID_WholeDomain[k][i][j]) - 1) % NGrainOrientations;
                    GOHistogram[GOVal]++;
                }
            }
        }
    }
    cout << "Histogram filled" << endl;
    for (int i=0; i<NGrainOrientations; i++) {
        Grainplot0 << GOHistogram[i] << endl;
    }
    Grainplot0.close();
}

void PrintGrainIDs(string FName, int nx, int ny, int nz, vector <vector <vector <bool> > > Melted_WholeDomain, vector <vector <vector <int> > > GrainID_WholeDomain) {

    string FName1 = FName + ".csv";
    cout << "Printing file of Grain ID values" << endl;
    std::ofstream Grainplot1;
    Grainplot1.open(FName1);
    Grainplot1 << "Outermost dimension is nz = " << nz-2 << endl;
    Grainplot1 << "ny = " << ny-2 << endl;
    Grainplot1 << "Innermost dimension is nx = " << nx-2 << endl;
    for (int k=1; k<nz-1; k++) {
        for (int j=1; j<ny-1; j++) {
            for (int i=1; i<nx-1; i++) {
                if (!Melted_WholeDomain[k][i][j]) Grainplot1 << 0 << endl;
                else Grainplot1 << GrainID_WholeDomain[k][i][j] << endl;
            }
        }
    }
    Grainplot1.close();
}

void PrintParaview(string FName, int nx, int ny, int nz, vector <vector <vector <bool> > > Melted_WholeDomain, vector <vector <vector <int> > > GrainID_WholeDomain, ViewI::HostMirror GrainOrientation, ViewF::HostMirror GrainUnitVector, int NGrainOrientations) {

    string FName2 = FName + ".vtk";
    cout << "Printing Paraview file of grain misorientations" << endl;
    // Print grain orientations to file
    std::ofstream Grainplot2;
    // New vtk file
    Grainplot2.open(FName2);
    Grainplot2 << "# vtk DataFile Version 3.0" << endl;
    Grainplot2 << "vtk output" << endl;
    Grainplot2 << "ASCII" << endl;
    Grainplot2 << "DATASET STRUCTURED_POINTS" << endl;
    Grainplot2 << "DIMENSIONS " << nx-2 << " " << ny-2 << " " << nz-2 << endl;
    Grainplot2 << "ORIGIN 0 0 0" << endl;
    Grainplot2 << "SPACING 1 1 1" << endl;
    Grainplot2 << fixed << "POINT_DATA " << (nx-2)*(ny-2)*(nz-2) << endl;
    Grainplot2 << "SCALARS Angle_z int 1" << endl;
    Grainplot2 << "LOOKUP_TABLE default" << endl;

    int NVC = 0;
    for (int k=1; k<nz-1; k++) {
        for (int j=1; j<ny-1; j++) {
            for (int i=1; i<nx-1; i++) {
                if (!Melted_WholeDomain[k][i][j]) Grainplot2 << 200 << " ";
                else {
                    int RoundedAngle;
                    int MyOrientation = GrainOrientation(((abs(GrainID_WholeDomain[k][i][j]) - 1) % NGrainOrientations));
                    double AngleZmin = 62.7;
                    for(int ll=0; ll<3; ll++) {
                        double AngleZ = abs((180/M_PI)*acos(GrainUnitVector(9*MyOrientation + 3*ll + 2)));
                        if (AngleZ < AngleZmin) {
                            AngleZmin = AngleZ;
                        }
                    }
                    if (GrainID_WholeDomain[k][i][j] < 0) {
                       RoundedAngle = round(AngleZmin) + 100;
                       NVC++;
                    }
                    else RoundedAngle = round(AngleZmin);
                    Grainplot2 << RoundedAngle << " ";
                }
            }
        }
        Grainplot2 << endl;
    }
    Grainplot2.close();
    cout << "Volume fraction of domain claimed by nucleated grains: " << (float)(NVC)/(float)((nx-2)*(ny-2)*(nz-2)) << endl;
}

void PrintGrainIDsForExaConstit(string FName, int nx, int ny, int nz, vector <vector <vector <int> > > GrainID_WholeDomain, double deltax) {
    
    string FName1 = FName + "_ExaConstit.csv";
    cout << "Printing file of X, Y, Z, GrainID values for ExaConstit" << endl;

    std::ofstream Grainplot1;
    Grainplot1.open(FName1);
    Grainplot1 << "Coordinates are in CA units, 1 cell = " << deltax << " microns. Data is cell-centered. Origin at 0,0,0, domain size is " << nx-2 << " by " << ny-2 << " by " << nz-2 << " cells" << endl;
    Grainplot1 << "X coord, Y coord, Z coord, Grain ID" << endl;
    for (int k=1; k<nz-1; k++) {
        for (int j=1; j<ny-1; j++) {
            for (int i=1; i<nx-1; i++) {
                Grainplot1 << i-1 << "," << j-1 << "," << k-1 << "," << GrainID_WholeDomain[k][i][j] << endl;
            }
        }
    }
    Grainplot1.close();
}

void PrintGrainAreas(string FName, double deltax, int nx, int ny, int nz, vector <vector <vector <int> > > GrainID_WholeDomain) {

    string FName1 = FName + "_GrainAreas.csv";
    cout << "Printing file of grain area values (in square microns) for all Z coordinates" << endl;
    std::ofstream Grainplot1;
    Grainplot1.open(FName1);
    for (int k=1; k<=nz-2; k++) {
        vector <int> GIDVals_ThisLayer;
        for (int j=1; j<=ny-2; j++) {
           for (int i=1; i<=nx-2; i++) {
                GIDVals_ThisLayer.push_back(GrainID_WholeDomain[k][i][j]);
            }
        }
        std::vector<int>::iterator it;
        sort(GIDVals_ThisLayer.begin(),GIDVals_ThisLayer.end());
        it = std::unique (GIDVals_ThisLayer.begin(), GIDVals_ThisLayer.end());
        int CellsThisLayer = GIDVals_ThisLayer.size();
        GIDVals_ThisLayer.resize( std::distance(GIDVals_ThisLayer.begin(),it) );
        int GrainsThisLayer = GIDVals_ThisLayer.size();
        double MeanGrainAreaThisLayer = (double)(CellsThisLayer) / (double)(GrainsThisLayer);
        Grainplot1 << MeanGrainAreaThisLayer*deltax*deltax/pow(10,-12) << endl;
    }
    Grainplot1.close();
}

void PrintWeightedGrainAreas(string FName, double deltax, int nx, int ny, int nz, vector <vector <vector <int> > > GrainID_WholeDomain) {

    string FName1 = FName + "_WeightedGrainAreas.csv";
    cout << "Printing file of Grain ID values (in square microns) for every 10th Z coordinate" << endl;
    std::ofstream Grainplot1;
    Grainplot1.open(FName1);
     for (int k=1; k<=nz-2; k++) {
       if ((k+2) % 12 == 0) {
         vector <int> GIDVals_ThisLayer;
         for (int j=1; j<=ny-2; j++) {
            for (int i=1; i<=nx-2; i++) {
                GIDVals_ThisLayer.push_back(GrainID_WholeDomain[k][i][j]);
             }
         }
         vector<int> GIDAllVals_ThisLayer;
         GIDAllVals_ThisLayer = GIDVals_ThisLayer;
         std::vector<int>::iterator it;
         sort(GIDVals_ThisLayer.begin(),GIDVals_ThisLayer.end());
         it = std::unique (GIDVals_ThisLayer.begin(), GIDVals_ThisLayer.end());
         int CellsThisLayer = GIDVals_ThisLayer.size();
         GIDVals_ThisLayer.resize( std::distance(GIDVals_ThisLayer.begin(),it) );
         int GrainsThisLayer = GIDVals_ThisLayer.size();
         long int AreaXArea = 0;
         for (int l=0; l<GrainsThisLayer; l++) {
             long int MyGrainArea = 0;
             for (int ll=0; ll<CellsThisLayer; ll++) {
                 if (GIDVals_ThisLayer[l] == GIDAllVals_ThisLayer[ll]) MyGrainArea++;
             }
             AreaXArea += MyGrainArea*MyGrainArea;
         }
         double WeightedArea = ((double)(AreaXArea) / (double)((nx-2)*(ny-2)));
         Grainplot1 << WeightedArea*deltax*deltax/pow(10,-12) << endl;
       }
    }
    Grainplot1.close();
    
}
