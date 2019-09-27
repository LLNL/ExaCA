#include "header.h"
using namespace std;

///*****************************************************************************/
/**
  Prints values of grain orientation for all cells to files
*/
void PrintValues(int id, int np,int nx, int ny,int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI::HostMirror GrainID, int* GrainOrientation,float* GrainUnitVector,string BaseFileName, int DecompositionStrategy, int NGrainOrientations) {

    if (id == 0) {
        // Create GrainID variable for entire domain, place GrainIDs for rank 0
        vector <vector <vector <int> > > GrainID_WholeDomain;
        for (int k=0; k<nz; k++) {
            vector <vector <int> > GrainID_JK;
            for (int i=0; i<nx; i++) {
                vector <int> GrainID_K;
                for (int j=0; j<ny; j++) {
                    GrainID_K.push_back(0);
                }
                GrainID_JK.push_back(GrainID_K);
            }
            GrainID_WholeDomain.push_back(GrainID_JK);
        }
        for (int k=1; k<nz-1; k++) {
            for (int i=1; i<MyXSlices-1; i++) {
                for (int j=1; j<MyYSlices-1; j++) {
                    GrainID_WholeDomain[k][i-1][j-1] = GrainID[k*MyXSlices*MyYSlices + i*MyYSlices + j];
                    //if (GrainID_WholeDomain[k][i-1][j-1] < 0) cout << "N 0" << endl;
                }
            }
        }

        //cout << "Domain created" << endl;
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
            //cout << "Allocating memory" << endl;
            int RBufSize = RecvXSlices*RecvYSlices*(nz-2);
            int *RecvBufGID = new int[RBufSize];
            //cout << "Rank 0 ready to recieve " << RBufSize << " data from rank " << p << endl;
            MPI_Recv(RecvBufGID,RBufSize,MPI_INT,p,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            //cout << "Placing data in Rows " << RecvXOffset << " through " << RecvXOffset+RecvXSlices-1 << " and in cols " << RecvYOffset << " through " << RecvYOffset+RecvYSlices-1 << endl;
            int DataCounter = 0;
            for (int k=1; k<nz-1; k++) {
                for (int i=0; i<RecvXSlices; i++) {
                    for (int j=0; j<RecvYSlices; j++) {
                        GrainID_WholeDomain[k][i+RecvXOffset][j+RecvYOffset] = RecvBufGID[DataCounter];
                        //if (GrainID_WholeDomain[k][i+RecvXOffset][j+RecvYOffset] < 0) cout << "NR" << endl;
                        DataCounter++;
                    }
                }
            }
            //cout << "DataCounter = " << DataCounter << " of " << RBufSize << endl;
        }
//        cout << "Opening file" << endl;
//        string FName1 = BaseFileName + ".csv";
//        std::ofstream Grainplot1;
//        Grainplot1.open(FName1);
//        for (int k=1; k<nz-1; k++) {
//            for (int j=1; j<ny-1; j++) {
//                for (int i=1; i<nx-1; i++) {
//                    if (GrainID_WholeDomain[k][i][j] != 0) {
//                        int MyOrientation = GrainOrientation[((abs(GrainID_WholeDomain[k][i][j]) - 1) % NGrainOrientations)];
//                        double AngleZmin = 54.7;
//                        for(int ll=0; ll<6; ll++) {
//                            double AngleZ = (180/M_PI)*acos(GrainUnitVector[18*MyOrientation + 3*ll + 2]);
//                            if (AngleZ < AngleZmin) {
//                                AngleZmin = AngleZ;
//                            }
//                        }
//                        Grainplot1 << i << "," << j << "," << k << "," << GrainID_WholeDomain[k][i][j] << "," << MyOrientation << endl;
//                    }
//                }
//            }
//            Grainplot1 << endl;
//        }
//        Grainplot1.close();
        
        string FName2 = BaseFileName + ".vtk";
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
                    if (GrainID_WholeDomain[k][i][j] == 0) {
                        Grainplot2 << 200 << " ";
                    }
                    else {
                        int RoundedAngle;
                        //cout << GrainID_WholeDomain[k][i][j] << endl;
                        int MyOrientation = GrainOrientation[((abs(GrainID_WholeDomain[k][i][j]) - 1) % NGrainOrientations)];
                        //cout << GrainID_WholeDomain[k][i][j] << endl;
                        double AngleZmin = 54.7;
                        for(int ll=0; ll<6; ll++) {
                            double AngleZ = (180/M_PI)*acos(GrainUnitVector[18*MyOrientation + 3*ll + 2]);
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
        cout << "Number of cells claimed by nucleated grains: " << NVC << endl;
        
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
        //cout << "Rank " << id << endl;
        int DataCounter = 0;
        for (int k=1; k<nz-1; k++) {
            for (int i=1; i<MyXSlices-1; i++) {
                for (int j=1; j<MyYSlices-1; j++) {
                    SendBufGID[DataCounter] = GrainID[k*MyXSlices*MyYSlices + i*MyYSlices + j];
//                    if (GrainID[k*MyXSlices*MyYSlices + i*MyYSlices + j] < 0) cout << "N " << np << endl;
//                    if (SendBufGID[DataCounter] < 0) cout << "S " << np << endl;
                    DataCounter++;
                }
            }
        }
        //cout << "Rank " << id << " sending " << SBufSize << " data points" << endl;
        MPI_Send(SendBufGID,SBufSize,MPI_INT,0,0,MPI_COMM_WORLD);
//        //cout << "RANK " << id << " SENT DATA" << endl;
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
        for (int k=1; k<nz-1; k++) {
            for (int i=1; i<MyXSlices-1; i++) {
                for (int j=1; j<MyYSlices-1; j++) {
                    CritTimeStep_WholeDomain[k][i-1][j-1] = CritTimeStep[k*MyXSlices*MyYSlices + i*MyYSlices + j];
                }
            }
        }
        // Recieve values from other ranks
        // Message size different for different ranks
        for (int p=1; p<np; p++) {
            //cout << "FEEDING " << p << " " << nx << " " << ProcessorsInXDirection << " " << np << " " << DecompositionStrategy << endl;
            int RecvXOffset = XOffsetCalc(p,nx,ProcessorsInXDirection,ProcessorsInYDirection,np,DecompositionStrategy);
            int RecvXSlices = XMPSlicesCalc(p,nx,ProcessorsInXDirection,ProcessorsInYDirection,np,DecompositionStrategy);

            int RecvYOffset = YOffsetCalc(p,ny,ProcessorsInYDirection,np,DecompositionStrategy);
            int RecvYSlices = YMPSlicesCalc(p,ny,ProcessorsInYDirection,np,DecompositionStrategy);

            // No ghost nodes in sent and recieved data
            RecvYSlices = RecvYSlices-2;
            RecvXSlices = RecvXSlices-2;
            RecvXOffset++;
            RecvYOffset++;
            int RBufSize = RecvXSlices*RecvYSlices*(nz-2);
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
                    SendBufGID[DataCounter] = CritTimeStep[k*MyXSlices*MyYSlices + i*MyYSlices + j];
                    DataCounter++;
                }
            }
        }
        //cout << "Rank " << id << " sending " << SBufSize << " int data points" << endl;
        MPI_Send(SendBufGID,SBufSize,MPI_INT,0,0,MPI_COMM_WORLD);
        //cout << "RANK " << id << " SENT DATA" << endl;
    }
    
//    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
//    // Cooling rate printed to file second
//    if (id == 0) {
//        // Create GrainID variable for entire domain, place GrainIDs for rank 0
//        vector <vector <vector <float> > > UndercoolingChange_WholeDomain;
//        for (int k=0; k<nz; k++) {
//            vector <vector <float> > UC_JK;
//            for (int i=0; i<nx; i++) {
//                vector <float> UC_K;
//                for (int j=0; j<ny; j++) {
//                    UC_K.push_back(0);
//                }
//                UC_JK.push_back(UC_K);
//            }
//            UndercoolingChange_WholeDomain.push_back(UC_JK);
//        }
//        for (int k=1; k<nz-1; k++) {
//            for (int i=1; i<MyXSlices-1; i++) {
//                for (int j=1; j<MyYSlices-1; j++) {
//                    UndercoolingChange_WholeDomain[k][i-1][j-1] = UndercoolingChange[k*MyXSlices*MyYSlices + i*MyYSlices + j];
//                }
//            }
//        }
//        // Recieve Grain ID value from other ranks
//        // Message size different for different ranks
//        for (int p=1; p<np; p++) {
//            // cout << "FEEDING " << p << " " << nx << " " << ProcessorsInXDirection << " " << np << " " << DecompositionStrategy << endl;
//            int RecvXOffset = XOffsetCalc(p,nx,ProcessorsInXDirection,ProcessorsInYDirection,np,DecompositionStrategy);
//            int RecvXSlices = XMPSlicesCalc(p,nx,ProcessorsInXDirection,ProcessorsInYDirection,np,DecompositionStrategy);
//            
//            int RecvYOffset = YOffsetCalc(p,ny,ProcessorsInYDirection,np,DecompositionStrategy);
//            int RecvYSlices = YMPSlicesCalc(p,ny,ProcessorsInYDirection,np,DecompositionStrategy);
//            
//            // No ghost nodes in send and recieved data
//            RecvYSlices = RecvYSlices-2;
//            RecvXSlices = RecvXSlices-2;
//            RecvXOffset++;
//            RecvYOffset++;
//            //cout << "Allocating memory" << endl;
//            int RBufSize = RecvXSlices*RecvYSlices*(nz-2);
//            int *RecvBufGID = new int[RBufSize];
//            //cout << "Rank 0 ready to recieve " << RecvXSlices << " " << RecvYSlices << " data from rank " << p << endl;
//            MPI_Recv(RecvBufGID,RBufSize,MPI_FLOAT,p,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//            //            //cout << "Recieved data from rank " << p << endl;
//            int DataCounter = 0;
//            for (int k=1; k<nz-1; k++) {
//                for (int i=0; i<RecvXSlices; i++) {
//                    for (int j=0; j<RecvYSlices; j++) {
//                        UndercoolingChange_WholeDomain[k][i+RecvXOffset][j+RecvYOffset] = RecvBufGID[DataCounter];
//                        DataCounter++;
//                    }
//                }
//            }
//            //cout << "DataCounter = " << DataCounter << " of " << RBufSize << endl;
//            //cout << "Placed data in Rows " << RecvXOffset << " through " << RecvXOffset+RecvXSlices-1 << " and in cols " << RecvYOffset << " through " << RecvYOffset+RecvYSlices-1 << endl;
//        }
//        
//        // Print grain orientations to file
//        string FName = "UndercoolingChange.vtk";
//        std::ofstream UCplot;
//        // New vtk file
//        UCplot.open(FName);
//        UCplot << "# vtk DataFile Version 3.0" << endl;
//        UCplot << "vtk output" << endl;
//        UCplot << "ASCII" << endl;
//        UCplot << "DATASET STRUCTURED_POINTS" << endl;
//        UCplot << "DIMENSIONS " << nx-2 << " " << ny-2 << " " << nz-2 << endl;
//        UCplot << "ORIGIN 0 0 0" << endl;
//        UCplot << "SPACING 1 1 1" << endl;
//        UCplot << fixed << "POINT_DATA " << (nx-2)*(ny-2)*(nz-2) << endl;
//        UCplot << "SCALARS Angle_z float 1" << endl;
//        UCplot << "LOOKUP_TABLE default" << endl;
//        
//        for (int k=1; k<nz-1; k++) {
//            for (int j=1; j<ny-1; j++) {
//                for (int i=1; i<nx-1; i++) {
//                    UCplot <<UndercoolingChange_WholeDomain[k][i][j] << " ";
//                }
//            }
//            UCplot << endl;
//        }
//        UCplot.close();
//    }
//    else {
//        // Send non-ghost node data to rank 0
//        //cout << "Rank " << id << endl;
//        int SBufSize = (MyXSlices-2)*(MyYSlices-2)*(nz-2);
//        int *SendBufGID = new int[SBufSize];
//        //cout << "Rank " << id << endl;
//        int DataCounter = 0;
//        // cout << "Rank " << id << " sending " << nx*nz*(MyYSlices-2) << " data points" << endl;
//        for (int k=1; k<nz-1; k++) {
//            for (int i=1; i<MyXSlices-1; i++) {
//                for (int j=1; j<MyYSlices-1; j++) {
//                    SendBufGID[DataCounter] = UndercoolingChange[k*MyXSlices*MyYSlices + i*MyYSlices + j];
//                    DataCounter++;
//                }
//            }
//        }
//        //cout << "Rank " << id << " sending " << MyXSlices-2 << " " << MyYSlices-2 << " data points" << endl;
//        MPI_Send(SendBufGID,SBufSize,MPI_FLOAT,0,0,MPI_COMM_WORLD);
//        //cout << "RANK " << id << " SENT DATA" << endl;
//    }
//
}
///*****************************************************************************/
/**
 Prints values of grain orientation for all cells to files - multiple layers of results
 */
void PrintValuesMultilayer(int NumberOfLayers, int LayerHeight, int id, int np, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI::HostMirror GrainID, int* GrainID_Stored, int* GrainOrientation, float* GrainUnitVector, string BaseFileName, int DecompositionStrategy, int NGrainOrientations) {
    // vector <vector <vector <int> > > &GrainID

    //    MPI_Barrier(MPI_COMM_WORLD);
    //    cout << "RANK " << id << endl;
    if (id == 0) {

        // Create GrainID variable for entire domain, place GrainIDs for rank 0
        vector <vector <vector <int> > > GrainID_WholeDomain;
        for (int k=0; k<(NumberOfLayers-1)*LayerHeight+nz-2; k++) {
            vector <vector <int> > GrainID_JK;
            for (int i=0; i<nx; i++) {
                vector <int> GrainID_K;
                for (int j=0; j<ny; j++) {
                    GrainID_K.push_back(0);
                }
                GrainID_JK.push_back(GrainID_K);
            }
            GrainID_WholeDomain.push_back(GrainID_JK);
        }

        int LastLayerOffset = (NumberOfLayers-1)*LayerHeight;
        for (int k=0; k<LastLayerOffset; k++)  {
            for(int i=1; i<MyXSlices-1; i++) {
                for(int j=1; j<MyYSlices-1; j++) {
                    GrainID_WholeDomain[k][i-1][j-1] = GrainID_Stored[k*MyXSlices*MyYSlices+i*MyYSlices+j];
                }
            }
        }

        for (int k=1; k<nz-1; k++) {
            for (int i=1; i<MyXSlices-1; i++) {
                for (int j=1; j<MyYSlices-1; j++) {
                    GrainID_WholeDomain[k-1+LastLayerOffset][i-1][j-1] = GrainID[k*MyXSlices*MyYSlices + i*MyYSlices + j];
                }
            }
        }

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

            int RBufSize = (RecvXSlices)*(RecvYSlices)*(LastLayerOffset+nz-2);
            int *RecvBufGID = new int[RBufSize];
            //cout << "RBuf size for rank " << p << " : " << RBufSize << endl;
            int DataCounter = 0;
            //cout << "Rank 0 recieved " << RBufSize << " data points" << endl;
            MPI_Recv(RecvBufGID,RBufSize,MPI_INT,p,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

            for (int k=0; k<LastLayerOffset; k++)  {
                for(int i=0; i<RecvXSlices; i++) {
                    for(int j=0; j<RecvYSlices; j++) {
                        GrainID_WholeDomain[k][i+RecvXOffset][j+RecvYOffset] = RecvBufGID[DataCounter];
                        DataCounter++;
                    }
                }
            }

            for (int k=1; k<nz-1; k++) {
                for (int i=0; i<RecvXSlices; i++) {
                    for (int j=0; j<RecvYSlices; j++) {
                        GrainID_WholeDomain[k-1+LastLayerOffset][i+RecvXOffset][j+RecvYOffset] = RecvBufGID[DataCounter];
                        DataCounter++;
                    }
                }
            }


        }

//        // Print grain orientations to file
//        string FName = BaseFileName + "_MultilayerG.vtk";
//        std::ofstream Grainplot;
//        // New vtk file
//        Grainplot.open(FName);
//        Grainplot << "# vtk DataFile Version 3.0" << endl;
//        Grainplot << "vtk output" << endl;
//        Grainplot << "ASCII" << endl;
//        Grainplot << "DATASET STRUCTURED_POINTS" << endl;
//        Grainplot << "DIMENSIONS " << nx-2 << " " << ny-2 << " " << LastLayerOffset+nz-2 << endl;
//        Grainplot << "ORIGIN 0 0 0" << endl;
//        Grainplot << "SPACING 1 1 1" << endl;
//        Grainplot << fixed << "POINT_DATA " << (nx-2)*(ny-2)*(LastLayerOffset+nz-2) << endl;
//        Grainplot << "SCALARS Angle_z int 1" << endl;
//        Grainplot << "LOOKUP_TABLE default" << endl;
//        for (int k=0; k<LastLayerOffset+nz-2; k++) {
//            for (int j=1; j<ny-1; j++) {
//                for (int i=1; i<nx-1; i++) {
//
//                    Grainplot << GrainID_WholeDomain[k][i][j] << " ";
//                }
//            }
//            Grainplot << endl;
//        }
//        Grainplot.close();
        
        cout << "Opening file 1" << endl;
        string FName1 = BaseFileName + "_MultilayerO.csv";
        std::ofstream Grainplot1;
        Grainplot1.open(FName1);
        for (int k=0; k<LastLayerOffset+nz-2; k++) {
            for (int j=1; j<ny-1; j++) {
                for (int i=1; i<nx-1; i++) {
                    if (GrainID_WholeDomain[k][i][j] != 0) {
                        int MyOrientation = GrainOrientation[((abs(GrainID_WholeDomain[k][i][j]) - 1) % NGrainOrientations)];
                        double AngleZmin = 54.7;
                        for(int ll=0; ll<6; ll++) {
                            double AngleZ = (180/M_PI)*acos(GrainUnitVector[18*MyOrientation + 3*ll + 2]);
                            if (AngleZ < AngleZmin) {
                                AngleZmin = AngleZ;
                            }
                        }
                        Grainplot1 << i << "," << j << "," << k << "," << GrainID_WholeDomain[k][i][j] << "," << MyOrientation << endl;
                    }
                }
            }
            Grainplot1 << endl;
        }
        Grainplot1.close();
        
        // Print grain orientations to file
        string FName = BaseFileName + "_MultilayerO.vtk";
        std::ofstream Grainplot2;
        // New vtk file
        Grainplot2.open(FName);
        Grainplot2 << "# vtk DataFile Version 3.0" << endl;
        Grainplot2 << "vtk output" << endl;
        Grainplot2 << "ASCII" << endl;
        Grainplot2 << "DATASET STRUCTURED_POINTS" << endl;
        Grainplot2 << "DIMENSIONS " << nx-2 << " " << ny-2 << " " << LastLayerOffset+nz-2 << endl;
        Grainplot2 << "ORIGIN 0 0 0" << endl;
        Grainplot2 << "SPACING 1 1 1" << endl;
        Grainplot2 << fixed << "POINT_DATA " << (nx-2)*(ny-2)*(LastLayerOffset+nz-2) << endl;
        Grainplot2 << "SCALARS Angle_z int 1" << endl;
        Grainplot2 << "LOOKUP_TABLE default" << endl;

        int NVC = 0;
        for (int k=0; k<LastLayerOffset+nz-2; k++) {
            for (int j=1; j<ny-1; j++) {
                for (int i=1; i<nx-1; i++) {
                    if (GrainID_WholeDomain[k][i][j] == 0) {
                        Grainplot2 << 200 << " ";
                    }
                    else {
                        int RoundedAngle;
                        //cout << GrainID_WholeDomain[k][i][j] << endl;
                        int MyOrientation = GrainOrientation[((abs(GrainID_WholeDomain[k][i][j]) - 1) % NGrainOrientations)];
                        //cout << GrainID_WholeDomain[k][i][j] << endl;
                        double AngleZmin = 54.7;
                        for(int ll=0; ll<6; ll++) {
                            double AngleZ = (180/M_PI)*acos(GrainUnitVector[18*MyOrientation + 3*ll + 2]);
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
        cout << "Fraction of cells claimed by nucleated grains: " << (double)(NVC)/(double)((LastLayerOffset+nz-2)*(nx-2)*(ny-2)) << endl;
//        if ((NumberOfLayers-1)*LayerHeight+nz-2-300 > 0) {
//
//
//            string FNameS = BaseFileName + "_Selection.vtk";
//            std::ofstream GrainplotS;
//            // New vtk file
//            GrainplotS.open(FNameS);
//            GrainplotS << "# vtk DataFile Version 3.0" << endl;
//            GrainplotS << "vtk output" << endl;
//            GrainplotS << "ASCII" << endl;
//            GrainplotS << "DATASET STRUCTURED_POINTS" << endl;
//            GrainplotS << "DIMENSIONS " << 100 << " " << 100 << " " << 300 << endl;
//            GrainplotS << "ORIGIN 0 0 0" << endl;
//            GrainplotS << "SPACING 1 1 1" << endl;
//            GrainplotS << fixed << "POINT_DATA " << 100*100*300 << endl;
//            GrainplotS << "SCALARS Angle_z int 1" << endl;
//            GrainplotS << "LOOKUP_TABLE default" << endl;
//
//
//            string FName2 = BaseFileName + "_Grains.csv";
//            std::ofstream GPlot;
//            // New csv file
//            GPlot.open(FName2);
//
//            vector <int> GDR;
//            for (int k=(NumberOfLayers-1)*LayerHeight+nz-2-300; k<(NumberOfLayers-1)*LayerHeight+nz-2; k++) {
//                for (int j=150; j<250; j++) {
//                    for (int i=600; i<700; i++) {
//                        GPlot << i << "," << j << "," << k << "," << GrainID_WholeDomain[k][i][j] << endl;
//                        GDR.push_back(GrainID_WholeDomain[k][i][j]);
//                        if (GrainID_WholeDomain[k][i][j] == 0) {
//                            GrainplotS << 60 << " ";
//                        }
//                        else {
//                            int RoundedAngle;
//                            int MyOrientation = GrainOrientation[GrainID_WholeDomain[k][i][j]];
//                            double AngleZmin = 54.7;
//                            for(int ll=0; ll<6; ll++) {
//                                double AngleZ = (180/M_PI)*acos(GrainUnitVector[18*MyOrientation + 3*ll + 2]);
//                                if (AngleZ < AngleZmin) {
//                                    AngleZmin = AngleZ;
//                                }
//                            }
//                            RoundedAngle = round(AngleZmin);
//                            GrainplotS << RoundedAngle << " ";
//                        }
//                    }
//                }
//                GrainplotS << endl;
//            }
//            GrainplotS.close();
//            GPlot.close();
//
//            vector <int> GDR_Unique;
//            GDR_Unique.push_back(GDR[0]);
//            for (int i=1; i<GDR.size(); i++) {
//                int U = 0;
//                // Check to see if element "i" of the list matches and of the unique GrainIDs in this region...
//                for (int j=0; j<GDR_Unique.size(); j++) {
//                    if (GDR[i] == GDR_Unique[j]) {
//                        U = 1;
//                        goto NotUnique;
//                    }
//                }
//                // ...and if it doesn't, add it to the list of unique GrainIDs
//                NotUnique:
//                if (U == 0) GDR_Unique.push_back(GDR[i]);
//            }
//            std::sort (GDR_Unique.begin(), GDR_Unique.end());
//            string FNameO = BaseFileName + "_Orientations.csv";
//            std::ofstream GrainplotO;
//            // New file for orientations of grains within the selection
//            GrainplotO.open(FNameO);
//            for (int j=0; j<GDR_Unique.size(); j++) {
//                GrainplotO << GDR_Unique[j] << "," << GrainOrientation[GDR_Unique[j]] << endl;
//            }
//            GrainplotO.close();
//        }
//
    }
    else {
        // Send non-ghost node data to rank 0
        int LastLayerOffset = (NumberOfLayers-1)*LayerHeight;
        int SBufSize = (MyXSlices-2)*(MyYSlices-2)*(LastLayerOffset+nz-2);
        int *SendBufGID = new int[SBufSize];
        //cout << "Rank " << id << " SBuf size : " << SBufSize << endl;
        int DataCounter = 0;
        //cout << "Rank " << id << " sending " << SBufSize << " data points" << endl;
        
        

        for (int k=0; k<LastLayerOffset; k++)  {
            for(int i=1; i<MyXSlices-1; i++) {
                for(int j=1; j<MyYSlices-1; j++) {
                    SendBufGID[DataCounter] = GrainID_Stored[k*MyXSlices*MyYSlices+i*MyYSlices+j];
                    DataCounter++;
                }
            }
        }

        for (int k=1; k<nz-1; k++) {
            for (int i=1; i<MyXSlices-1; i++) {
                for (int j=1; j<MyYSlices-1; j++) {
                    SendBufGID[DataCounter] = GrainID[k*MyXSlices*MyYSlices+i*MyYSlices+j];
                    DataCounter++;
                }
            }
        }
        
        //cout << "Rank " << id << " sending " <<  SBufSize << " data points" << endl;
        MPI_Send(SendBufGID,SBufSize,MPI_INT,0,0,MPI_COMM_WORLD);
        //cout << "RANK " << id << " SENT DATA" << endl;
    }
    
}

///*****************************************************************************/

void PrintCT(int id, int np,int nx, int ny,int nz, int MyXSlices, int MyYSlices, int ProcessorsInXDirection, int ProcessorsInYDirection, ViewI::HostMirror CellType, string BaseFileName, int DecompositionStrategy) {
    // vector <vector <vector <int> > > &GrainID
    string FName = BaseFileName + "_CT.vtk";
    std::ofstream Grainplot;
    //    MPI_Barrier(MPI_COMM_WORLD);
    //    cout << "RANK " << id << endl;
    if (id == 0) {
        // Create GrainID variable for entire domain, place GrainIDs for rank 0
        vector <vector <vector <char> > > GrainID_WholeDomain;
        for (int k=0; k<nz; k++) {
            vector <vector <char> > GrainID_JK;
            for (int i=0; i<nx; i++) {
                vector <char> GrainID_K;
                for (int j=0; j<ny; j++) {
                    GrainID_K.push_back(0);
                }
                GrainID_JK.push_back(GrainID_K);
            }
            GrainID_WholeDomain.push_back(GrainID_JK);
        }
        for (int k=1; k<nz-1; k++) {
            for (int i=1; i<MyXSlices-1; i++) {
                for (int j=1; j<MyYSlices-1; j++) {
                    GrainID_WholeDomain[k][i-1][j-1] = CellType[k*MyXSlices*MyYSlices + i*MyYSlices + j];
                }
            }
        }
        
        //cout << "Domain created" << endl;
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
            //cout << "Allocating memory" << endl;
            int RBufSize = RecvXSlices*RecvYSlices*(nz-2);
            char *RecvBufGID = new char[RBufSize];
            //cout << "Rank 0 ready to recieve " << RBufSize << " data from rank " << p << endl;
            MPI_Recv(RecvBufGID,RBufSize,MPI_CHAR,p,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            //cout << "Placing data in Rows " << RecvXOffset << " through " << RecvXOffset+RecvXSlices-1 << " and in cols " << RecvYOffset << " through " << RecvYOffset+RecvYSlices-1 << endl;
            int DataCounter = 0;
            for (int k=1; k<nz-1; k++) {
                for (int i=0; i<RecvXSlices; i++) {
                    for (int j=0; j<RecvYSlices; j++) {
                        GrainID_WholeDomain[k][i+RecvXOffset][j+RecvYOffset] = RecvBufGID[DataCounter];
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
                    if ((GrainID_WholeDomain[k][i][j] == 'S')||(GrainID_WholeDomain[k][i][j] == 'W')) Grainplot << 0 << " ";
                    else if (GrainID_WholeDomain[k][i][j] == 'L') Grainplot << 2 << " ";
                    else if (GrainID_WholeDomain[k][i][j] == 'A') Grainplot << 1 << " ";
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
        char *SendBufGID = new char[SBufSize];
        //cout << "Rank " << id << endl;
        int DataCounter = 0;
        for (int k=1; k<nz-1; k++) {
            for (int i=1; i<MyXSlices-1; i++) {
                for (int j=1; j<MyYSlices-1; j++) {
                    SendBufGID[DataCounter] = CellType[k*MyXSlices*MyYSlices + i*MyYSlices + j];
                    DataCounter++;
                }
            }
        }
        //cout << "Rank " << id << " sending " << SBufSize << " data points" << endl;
        MPI_Send(SendBufGID,SBufSize,MPI_CHAR,0,0,MPI_COMM_WORLD);
        //cout << "RANK " << id << " SENT DATA" << endl;
    }
    
}

