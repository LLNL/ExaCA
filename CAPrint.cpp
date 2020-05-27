#include "header.h"
using namespace std;

///*****************************************************************************/
/**
  Prints values of grain orientation for all cells to files
*/
void PrintValues(int id, int np, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI::HostMirror GrainID, int* GrainOrientation, float* GrainUnitVector, string BaseFileName, int DecompositionStrategy, int NGrainOrientations, bool* Melted) {

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
        
        // File 1: Histogram of orientations for texture determination
        ofstream Grainplot0;
        string FName0 = BaseFileName + "_Orientations.csv";
        Grainplot0.open(FName0);
        int GOHistogram[NGrainOrientations];
        for (int i=0; i<NGrainOrientations; i++) {
            GOHistogram[i] = 0;
        }
        // frequency data on grain ids
        for (int k=2; k<nz-1; k++) {
            for (int j=1; j<ny-1; j++) {
                for (int i=1; i<nx-1; i++) {
                    if (Melted_WholeDomain[k][i][j]) {
                        int GOVal = (abs(GrainID_WholeDomain[k][i][j]) - 1) % NGrainOrientations;
                        GOHistogram[GOVal]++;
                    }
                }
            }
        }
        for (int i=0; i<NGrainOrientations; i++) {
            Grainplot0 << GOHistogram[i] << endl;
        }
        Grainplot0.close();
        
//        // "/gpfs/alpine/world-shared/mat190/rolchigo/" +
        string FName1 = BaseFileName + ".csv";
        cout << "Opening file " << FName1 << endl;
        std::ofstream Grainplot1;
        Grainplot1.open(FName1);
        Grainplot1 << "Outermost dimension is nz = " << nz-2 << endl;
        Grainplot1 << "ny = " << ny-2 << endl;
        Grainplot1 << "Innermost dimension is nx = " << nx-2 << endl;
        for (int k=1; k<nz-2; k++) {
            for (int j=1; j<ny-1; j++) {
                for (int i=1; i<nx-1; i++) {
                    Grainplot1 << GrainID_WholeDomain[k][i][j] << endl;
                }
            }
        }
        Grainplot1.close();
        

        string FName2 = BaseFileName + ".vtk";
        cout << "Opening file " << FName2 << endl;
        // Print grain orientations to file
        std::ofstream Grainplot2;
        // New vtk file
        Grainplot2.open(FName2);
        Grainplot2 << "# vtk DataFile Version 3.0" << endl;
        Grainplot2 << "vtk output" << endl;
        Grainplot2 << "ASCII" << endl;
        Grainplot2 << "DATASET STRUCTURED_POINTS" << endl;
        Grainplot2 << "DIMENSIONS " << nx << " " << ny << " " << nz-2 << endl;
        Grainplot2 << "ORIGIN 0 0 0" << endl;
        Grainplot2 << "SPACING 1 1 1" << endl;
        Grainplot2 << fixed << "POINT_DATA " << (nx)*(ny)*(nz-2) << endl;
        Grainplot2 << "SCALARS Angle_z int 1" << endl;
        Grainplot2 << "LOOKUP_TABLE default" << endl;
        
        int NVC = 0;
        for (int k=1; k<nz-1; k++) {
            for (int j=0; j<ny; j++) {
                for (int i=0; i<nx; i++) {
                    if (!Melted_WholeDomain[k][i][j]) Grainplot2 << 200 << " ";
                    else {
                        //if (GrainID_WholeDomain[k][i][j] == -1) {
                        //    Grainplot2 << 300 << " ";
                        //}
                        //else {
                            int RoundedAngle;
                            int MyOrientation = GrainOrientation[((abs(GrainID_WholeDomain[k][i][j]) - 1) % NGrainOrientations)];
                            double AngleZmin = 54.7;
                            for(int ll=0; ll<3; ll++) {
                                double AngleZ = abs((180/M_PI)*acos(GrainUnitVector[9*MyOrientation + 3*ll + 2]));
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
                        //}
                    }
                }
            }
            Grainplot2 << endl;
        }
        Grainplot2.close();
        cout << "Volume fraction of domain claimed by nucleated grains: " << (float)(NVC)/(float)((nx*ny*(nz-2))) << endl;
        
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
 Prints values of critical undercooling and critical solidification time step for all cells to files
 */
void PrintTempValues(int id, int np, int nx, int ny, int nz, int MyXSlices,int MyYSlices,int ProcessorsInXDirection,int ProcessorsInYDirection, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange, int DecompositionStrategy) {
    
    
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
        string FName = "CriticalTimeStep.vtk";
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
                    CTSplot << CritTimeStep_WholeDomain[k][i][j] << " ";
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
/**
 Prints values of grain orientation for all cells to files - multiple layers of results
 */
void PrintValuesMultilayer(int ZSizeStore, int id, int np, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI::HostMirror GrainID, int* GrainID_Stored, int* GrainOrientation, float* GrainUnitVector, string BaseFileName, int DecompositionStrategy, int NGrainOrientations, bool* Melted, bool* MeltedStored) {
    // vector <vector <vector <int> > > &GrainID

    //    MPI_Barrier(MPI_COMM_WORLD);
    //    cout << "RANK " << id << endl;
    if (id == 0) {

        // Create GrainID variable for entire domain, place GrainIDs for rank 0
        vector <vector <vector <int> > > GrainID_WholeDomain;
        vector <vector <vector <bool> > > Melted_WholeDomain;
        for (int k=0; k<ZSizeStore+nz-2; k++) {
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
        cout << "Stored size: " << ZSizeStore << endl;
        for (int k=0; k<ZSizeStore; k++)  {
            for(int i=1; i<MyXSlices-1; i++) {
                for(int j=1; j<MyYSlices-1; j++) {
                    GrainID_WholeDomain[k][i-1][j-1] = GrainID_Stored[k*MyXSlices*MyYSlices+i*MyYSlices+j];
                    Melted_WholeDomain[k][i-1][j-1] = MeltedStored[k*MyXSlices*MyYSlices+i*MyYSlices+j];
                }
            }
        }
        cout << "Active size: " << nz-2 << endl;
        for (int k=1; k<nz-1; k++) {
            for (int i=1; i<MyXSlices-1; i++) {
                for (int j=1; j<MyYSlices-1; j++) {
                    GrainID_WholeDomain[k-1+ZSizeStore][i-1][j-1] = GrainID(k*MyXSlices*MyYSlices + i*MyYSlices + j);
                    Melted_WholeDomain[k-1+ZSizeStore][i-1][j-1] = Melted[k*MyXSlices*MyYSlices + i*MyYSlices + j];
                }
            }
        }
        //cout << "Z" << endl;
        // Recieve Grain ID value from other ranks
        // Message size different for different ranks
        for (int p=1; p<np; p++) {
            // cout << "FEEDING " << p << " " << nx << " " << ProcessorsInXDirection << " " << np << " " << DecompositionStrategy << endl;
            int RecvXOffset = XOffsetCalc(p,nx,ProcessorsInXDirection,ProcessorsInYDirection,np,DecompositionStrategy);
            int RecvXSlices = XMPSlicesCalc(p,nx,ProcessorsInXDirection,ProcessorsInYDirection,np,DecompositionStrategy);
            
            int RecvYOffset = YOffsetCalc(p,ny,ProcessorsInYDirection,np,DecompositionStrategy);
            int RecvYSlices = YMPSlicesCalc(p,ny,ProcessorsInYDirection,np,DecompositionStrategy);
            
            // No ghost nodes in send and recieved data
            RecvYSlices = RecvYSlices-2;
            RecvXSlices = RecvXSlices-2;
            RecvXOffset++;
            RecvYOffset++;

            int RBufSize = (RecvXSlices)*(RecvYSlices)*(ZSizeStore+nz-2);
            int *RecvBufGID = new int[RBufSize];
            int *RecvBufM = new int[RBufSize];
            //cout << "RBuf size for rank " << p << " : " << RBufSize << endl;
            int DataCounter = 0;
            //cout << "Rank 0 recieved " << RBufSize << " data points" << endl;
            MPI_Recv(RecvBufGID,RBufSize,MPI_INT,p,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(RecvBufM,RBufSize,MPI_INT,p,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            for (int k=0; k<ZSizeStore; k++)  {
                for(int i=0; i<RecvXSlices; i++) {
                    for(int j=0; j<RecvYSlices; j++) {
                        GrainID_WholeDomain[k][i+RecvXOffset][j+RecvYOffset] = RecvBufGID[DataCounter];
                        if (RecvBufM[DataCounter] == 1) Melted_WholeDomain[k][i+RecvXOffset][j+RecvYOffset] = true;
                        else Melted_WholeDomain[k][i+RecvXOffset][j+RecvYOffset] = false;
                        DataCounter++;
                    }
                }
            }

            for (int k=1; k<nz-1; k++) {
                for (int i=0; i<RecvXSlices; i++) {
                    for (int j=0; j<RecvYSlices; j++) {
                        GrainID_WholeDomain[k-1+ZSizeStore][i+RecvXOffset][j+RecvYOffset] = RecvBufGID[DataCounter];
                        if (RecvBufM[DataCounter] == 1) Melted_WholeDomain[k-1+ZSizeStore][i+RecvXOffset][j+RecvYOffset] = true;
                        else Melted_WholeDomain[k-1+ZSizeStore][i+RecvXOffset][j+RecvYOffset] = false;
                        DataCounter++;
                    }
                }
            }


        }

        //string SummitExt = "/gpfs/alpine/mat190/world-shared/rolchigo";
        // File 1: Histogram of orientations for texture determination
        ofstream Grainplot1;
        string FName1 = BaseFileName  + "_MultilayerOrientations.csv";
        Grainplot1.open(FName1);
        int GOHistogram[NGrainOrientations];
        for (int i=0; i<NGrainOrientations; i++) {
            GOHistogram[i] = 0;
        }
       //  frequency data on grain ids
        for (int k=0; k<ZSizeStore+nz-2; k++) {
            for (int j=1; j<ny-1; j++) {
                for (int i=1; i<nx-1; i++) {
                    //int D3D1ConvPosition = k*MyXSlices*MyYSlices + i*MyYSlices + j;
                    if (Melted_WholeDomain[k][i][j]) {
                        int GOVal = (abs(GrainID_WholeDomain[k][i][j]) - 1) % NGrainOrientations;
                        GOHistogram[GOVal]++;
                    }
                }
            }
        }
        for (int i=0; i<NGrainOrientations; i++) {
            Grainplot1 << GOHistogram[i] << endl;
        }
        Grainplot1.close();
        
        // File 2: GrainID values for grain size/aspect ratio stats
        string FName2 = BaseFileName  + "_MultilayerGrainIDs.csv";
        ofstream Grainplot2;
        Grainplot2.open(FName2);
        Grainplot2 << "Outermost dimension is nz = " << ZSizeStore+nz-2 << endl;
        Grainplot2 << "ny = " << ny-2 << endl;
        Grainplot2 << "Innermost dimension is nx = " << nx-2 << endl;
        for (int k=0; k<ZSizeStore+nz-2; k++) {
            for (int j=1; j<ny-1; j++) {
                for (int i=1; i<nx-1; i++) {
                    int D3D1ConvPosition = k*MyXSlices*MyYSlices + i*MyYSlices + j;
                    if (Melted_WholeDomain[k][i][j]) {
                      Grainplot2 << GrainID_WholeDomain[k][i][j] << endl;
                    }
                    else {
                        Grainplot2 << "0" << endl;
                    }
                }
            }
        }
        Grainplot2.close();
        
        // File 3: Grain orientations relative to vertical direction for paraview
        string FName3 = BaseFileName  + "_Multilayer.vtk";
        std::ofstream Grainplot3;
        // New vtk file
        Grainplot3.open(FName3);
        Grainplot3 << "# vtk DataFile Version 3.0" << endl;
        Grainplot3 << "vtk output" << endl;
        Grainplot3 << "ASCII" << endl;
        Grainplot3 << "DATASET STRUCTURED_POINTS" << endl;
        Grainplot3 << "DIMENSIONS " << nx << " " << ny << " " << ZSizeStore+nz-2 << endl;
        Grainplot3 << "ORIGIN 0 0 0" << endl;
        Grainplot3 << "SPACING 1 1 1" << endl;
        Grainplot3 << fixed << "POINT_DATA " << (nx)*(ny)*(ZSizeStore+nz-2) << endl;
        Grainplot3 << "SCALARS Angle_z int 1" << endl;
        Grainplot3 << "LOOKUP_TABLE default" << endl;
        
        int NVC = 0;
        for (int k=0; k<ZSizeStore+nz-2; k++) {
            for (int j=0; j<ny; j++) {
                for (int i=0; i<nx; i++) {
                    if (!Melted_WholeDomain[k][i][j]) Grainplot3 << 200 << " ";
                    else {
                        //if (GrainID_WholeDomain[k][i][j] == -1) {
                        //    Grainplot3 << 300 << " ";
                        //}
                        //else {
                            int RoundedAngle;
                            //cout << GrainID_WholeDomain[k][i][j] << endl;
                            int MyOrientation = GrainOrientation[((abs(GrainID_WholeDomain[k][i][j]) - 1) % NGrainOrientations)];
                            //cout << GrainID_WholeDomain[k][i][j] << endl;
                            double AngleZmin = 54.7;
                            for(int ll=0; ll<3; ll++) {
                                double AngleZ = abs((180/M_PI)*acos(GrainUnitVector[9*MyOrientation + 3*ll + 2]));
                                if (AngleZ < AngleZmin) {
                                    AngleZmin = AngleZ;
                                }
                            }
                            if (GrainID_WholeDomain[k][i][j] < 0) {
                                RoundedAngle = round(AngleZmin) + 100;
                                NVC++;
                            }
                            else {
                                RoundedAngle = round(AngleZmin);
                            }
                            Grainplot3 << RoundedAngle << " ";
                        //}
                    }
                }
            }
            Grainplot3 << endl;
        }
        Grainplot3.close();
        //cout << "Fraction of cells claimed by nucleated grains: " << (double)(NVC)/(double)((ZSizeStore+nz-2)*(nx-2)*(ny-2)) << endl;
    }
    else {
        // Send non-ghost node data to rank 0
        int SBufSize = (MyXSlices-2)*(MyYSlices-2)*(ZSizeStore+nz-2);
        int *SendBufGID = new int[SBufSize];
        int *SendBufM = new int[SBufSize];
        //cout << "Rank " << id << " SBuf size : " << SBufSize << endl;
        int DataCounter = 0;
        //cout << "Rank " << id << " sending " << SBufSize << " data points" << endl;
        
        

        for (int k=0; k<ZSizeStore; k++)  {
            for(int i=1; i<MyXSlices-1; i++) {
                for(int j=1; j<MyYSlices-1; j++) {
                    SendBufGID[DataCounter] = GrainID_Stored[k*MyXSlices*MyYSlices+i*MyYSlices+j];
                    if (MeltedStored[k*MyXSlices*MyYSlices+i*MyYSlices+j]) SendBufM[DataCounter] = 1;
                    else SendBufM[DataCounter] = 0;
                    DataCounter++;
                }
            }
        }

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
        //cout << "RANK " << id << " SENT DATA" << endl;
    }
    
}

///*****************************************************************************/

void PrintCT(int id, int np, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI::HostMirror CellType, string BaseFileName, int DecompositionStrategy) {
    // vector <vector <vector <int> > > &GrainID
    string FName = BaseFileName + "_CT.vtk";
    std::ofstream Grainplot;
    //    MPI_Barrier(MPI_COMM_WORLD);
    //    cout << "RANK " << id << endl;
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
                    CellType_WholeDomain[k][i][j] = CellType(k*MyXSlices*MyYSlices + i*MyYSlices + j);
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
                    else if (CellType_WholeDomain[k][i][j] == Wall) Grainplot << 10 << " ";
                    else if (CellType_WholeDomain[k][i][j] == Liquid) Grainplot << 2 << " ";
                    else if (CellType_WholeDomain[k][i][j] == Active) Grainplot << 1 << " ";
                    else  Grainplot << 3 << " ";
                }
            }
            Grainplot << endl;
        }
        Grainplot.close();
        
        //        std::ofstream Grainplot2;
        //        // New csv file of centerline
        //        Grainplot2.open("Centerline.csv");
        //        int k = nz-3;
        //        for (int j=1; j<ny-1; j++) {
        //            for (int i=1; i<nx-1; i++) {
        //                if (GrainID_WholeDomain[k][i][j] == 0) {
        //                    Grainplot2 << 60 << " ";
        //                }
        //                else {
        //                    if (GrainID_WholeDomain[k][i][j] != 1) {
        //                        cout << i << " " << j << " " << GrainID_WholeDomain[k][i][j] << endl;
        //                    }
        //                    int RoundedAngle;
        //                    //cout << GrainID_WholeDomain[k][i][j] << endl;
        //                    int MyOrientation = GrainOrientation[GrainID_WholeDomain[k][i][j]];
        //                    //cout << GrainID_WholeDomain[k][i][j] << endl;
        //                    double AngleZmin = 54.7;
        //                    for(int ll=0; ll<6; ll++) {
        //                        double AngleZ = (180/M_PI)*acos(GrainUnitVector[18*MyOrientation + 3*ll + 2]);
        //                        if (AngleZ < AngleZmin) {
        //                            AngleZmin = AngleZ;
        //                        }
        //                    }
        //                    RoundedAngle = round(AngleZmin);
        //                    Grainplot2 << RoundedAngle;
        //
        //                }
        //                if (i < nx-2) Grainplot2 << ",";
        //                else Grainplot2 << endl;
        //            }
        //        }
        //        Grainplot2.close();
        
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

