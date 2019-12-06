#include "header.h"
using namespace std;

void RunProgram_Bulky(int id, int np, int ierr, int DecompositionStrategy, double deltax, double AConst, double BConst, double CConst, double NMax, double dTN, double dTsigma, string BaseFileName, string GrainOrientationFile) {
    
//    double InitTime = MPI_Wtime();
//    bool UseSubstrateFile;
//    int NRatio, nx, ny, nz, NumberOfLayers, LayerHeight;
//    double HT_deltax, deltat, DomUndercooling, InitialGrainWidth, G, R;
//    string BaseFileName, GrainOrientationFile, TemperatureDataType, AFile1, AFile2, AFile3, tempfile;
//
//    if (TemperatureDataType == "A") {
//        AInputRead(AFile1, AFile2, AFile3, NRatio, NumberOfLayers, LayerHeight, UseSubstrateFile);
//    }
//    else if (TemperatureDataType == "C") {
//        CInputRead(G, R, nx, ny, nz, NRatio, NumberOfLayers);
//        UseSubstrateFile = false;
//    }
//    else if (TemperatureDataType == "U") {
//        UInputRead(G, R, DomUndercooling, nx, ny, nz, NumberOfLayers);
//        UseSubstrateFile = false;
//    }
//    else {
//        RInputRead(tempfile, HT_deltax, deltat, NumberOfLayers, LayerHeight, BaseFileName, UseSubstrateFile);
//    }
//    
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (id == 0) cout << "Input file read" << endl;
//    
//    // Grid decomposition
//    int ProcessorsInXDirection, ProcessorsInYDirection;
//    // Variables characterizing local processor grids relative to global domain
//    int MyXSlices, MyXOffset, MyYSlices, MyYOffset, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyLeftOut, MyRightIn, MyRightOut;
//    // Neighbor lists for cells
//    int NeighborX[26], NeighborY[26], NeighborZ[26], ItList[9][26];
//    float XMin, YMin, ZMin, XMax, YMax, ZMax; // OpenFOAM simulation bounds (if using OpenFOAM data)
//
//    // Initialization of the grid and decomposition, along with deltax and deltat
//    ParallelMeshInit(G, R, DecompositionStrategy, NeighborX, NeighborY, NeighborZ, ItList, TemperatureDataType, ierr, id, np, MyXSlices, MyYSlices, MyXOffset, MyYOffset, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyLeftOut, MyRightIn, MyRightOut, deltax, HT_deltax, deltat, nx, ny, nz, ProcessorsInXDirection, ProcessorsInYDirection, AFile3, tempfile, XMin, XMax, YMin, YMax, ZMin, ZMax, DomUndercooling, NRatio);
//    int LocalDomainSize = MyXSlices*MyYSlices*nz; // Number of cells on this MPI rank
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (id == 0) cout << "Mesh initialized" << endl;
//    
//    // Temperature fields characterized by these variables:
////    int* CritTimeStep = new int[LocalDomainSize];
////    float* UndercoolingChange = new float[LocalDomainSize];
////    float* UndercoolingCurrent = new float[LocalDomainSize];
//     ViewI CritTimeStep_G("CritTimeStep_G", LocalDomainSize);
//     ViewF UndercoolingChange_G("UndercoolingChange_G",LocalDomainSize);
//     ViewF UndercoolingCurrent_G("UndercoolingCurrent_G",LocalDomainSize);
//     ViewI::HostMirror CritTimeStep_H = Kokkos::create_mirror_view( CritTimeStep_G );
//     ViewF::HostMirror UndercoolingChange_H = Kokkos::create_mirror_view( UndercoolingChange_G );
//     ViewF::HostMirror UndercoolingCurrent_H = Kokkos::create_mirror_view( UndercoolingCurrent_G );
//    
//
//    // CA cell variables
////    int* GrainID = new int[LocalDomainSize];
////    char* CellType = new char[LocalDomainSize];
//    ViewI GrainID_G("GrainID_G",LocalDomainSize);
//    ViewI CellType_G("CellType_G",LocalDomainSize);
//    ViewI::HostMirror GrainID_H = Kokkos::create_mirror_view( GrainID_G );
//    ViewI::HostMirror CellType_H = Kokkos::create_mirror_view( CellType_G );
//    bool* Melted = new bool[LocalDomainSize];
//    
//    // Initialize the temperature fields
//    TempInit(G, R, DecompositionStrategy,NeighborX, NeighborY, NeighborZ, ItList, TemperatureDataType, ierr, id, np, MyXSlices, MyYSlices, MyXOffset, MyYOffset, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyLeftOut, MyRightIn, MyRightOut, deltax, HT_deltax, deltat, nx, ny, nz, ProcessorsInXDirection, ProcessorsInYDirection,  CritTimeStep_H, UndercoolingChange_H, UndercoolingCurrent_H, DomUndercooling, tempfile, AFile1, AFile2, XMin, XMax, YMin, YMax, ZMin, ZMax, Melted);
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (id == 0) cout << "Done with temperature field initialization " << endl;
//    
//    int NGrainOrientations = 10000; // Number of grain orientations considered in the simulation
//    float* GrainUnitVector = new float[18*NGrainOrientations];
//    int* GrainOrientation = new int[NGrainOrientations];
//    // Initialize grain orientations
//    OrientationInit(id, NGrainOrientations, GrainOrientation, GrainUnitVector, GrainOrientationFile);
//    
//    // Variables characterizing the active cell region within each rank's grid
////    float* DiagonalLength = new float[LocalDomainSize];
////    float* CritDiagonalLength = new float[26*LocalDomainSize];
////    float* DOCenter = new float[3*LocalDomainSize];
////    int* TriangleIndex = new int[26*3*LocalDomainSize];
//    ViewF DiagonalLength_G("DiagonalLength_G",LocalDomainSize);
//    ViewF CritDiagonalLength_G("CritDiagonalLength_G",26*LocalDomainSize);
//    ViewF DOCenter_G("DOCenter_G",3*LocalDomainSize);
//    ViewI TriangleIndex_G("TriangleIndex_G",26*3*LocalDomainSize);
//    ViewF::HostMirror DiagonalLength_H = Kokkos::create_mirror_view( DiagonalLength_G );
//    ViewF::HostMirror CritDiagonalLength_H = Kokkos::create_mirror_view( CritDiagonalLength_G );
//    ViewF::HostMirror DOCenter_H = Kokkos::create_mirror_view( DOCenter_G );
//    ViewI::HostMirror TriangleIndex_H = Kokkos::create_mirror_view( TriangleIndex_G );
//    
//    // Initialize the grain structure
//    int PossibleNuclei_ThisRank, NextLayer_FirstNucleatedGrainID;
//    GrainInit(UseSubstrateFile, InitialGrainWidth, NGrainOrientations, DecompositionStrategy, ProcessorsInXDirection, ProcessorsInYDirection, nx, ny, nz, MyXSlices, MyYSlices, MyXOffset, MyYOffset, id, np, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, ItList, NeighborX, NeighborY, NeighborZ, GrainOrientation, GrainUnitVector, DiagonalLength_H, CellType_H, TriangleIndex_H, GrainID_H, CritDiagonalLength_H, DOCenter_H, CritTimeStep_H, UndercoolingChange_H, Melted, deltax, NMax, NextLayer_FirstNucleatedGrainID, PossibleNuclei_ThisRank);
//    
//    // Update ghost nodes for grain locations and attributes
//    if (DecompositionStrategy == 1) GhostNodes1D(0, id, MyLeft, MyRight, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, NeighborX, NeighborY, NeighborZ, CellType_H, DOCenter_H, GrainID_H, GrainUnitVector,TriangleIndex_H, GrainOrientation, DiagonalLength_H, CritDiagonalLength_H, NGrainOrientations);
//    else GhostNodes2D(0, id, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, NeighborX, NeighborY, NeighborZ, CellType_H, DOCenter_H, GrainID_H, GrainUnitVector, TriangleIndex_H, GrainOrientation, DiagonalLength_H, CritDiagonalLength_H, NGrainOrientations);
//    
//    ViewI NucleationTimes_G("NucleationTimes_G",PossibleNuclei_ThisRank);
//    ViewI NucleiLocation_G("NucleiLocation_G",PossibleNuclei_ThisRank);
//    ViewI::HostMirror NucleationTimes_H = Kokkos::create_mirror_view( NucleationTimes_G );
//    ViewI::HostMirror NucleiLocation_H = Kokkos::create_mirror_view( NucleiLocation_G );
//    
//    // Update nuclei on ghost nodes, fill in nucleation data structures, and assign nucleation undercooling values to potential nucleation events
//    NucleiInit(DecompositionStrategy, MyXSlices, MyYSlices, id, dTN, dTsigma, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, PossibleNuclei_ThisRank, NucleiLocation_H, NucleationTimes_H, GrainOrientation, CellType_H, GrainID_H, CritTimeStep_H, UndercoolingChange_H);
//    
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (id == 0) cout << "Done with grain initialization " << endl;
//    
//    // Normalize solidification parameters
//    AConst = AConst*deltat/deltax;
//    BConst = BConst*deltat/deltax;
//    CConst = CConst*deltat/deltax;
//
//    double InitTime2 = MPI_Wtime();
//    
//    int* GrainID_Stored = new int[MyXSlices*MyYSlices*((NumberOfLayers-1)*LayerHeight)];
//    bool* MeltedStored = new bool[MyXSlices*MyYSlices*((NumberOfLayers-1)*LayerHeight)];
//    int cycle;
//    
//    // Copy view data to GPU
//    Kokkos::deep_copy( GrainID_G, GrainID_H );
//    Kokkos::deep_copy( CellType_G, CellType_H );
//    Kokkos::deep_copy( DiagonalLength_G, DiagonalLength_H );
//    Kokkos::deep_copy( CritDiagonalLength_G, CritDiagonalLength_H );
//    Kokkos::deep_copy( DOCenter_G, DOCenter_H );
//    Kokkos::deep_copy( TriangleIndex_G, TriangleIndex_H );
//    Kokkos::deep_copy( CritTimeStep_G, CritTimeStep_H );
//    Kokkos::deep_copy( UndercoolingChange_G, UndercoolingChange_H );
//    Kokkos::deep_copy( UndercoolingCurrent_G, UndercoolingCurrent_H );
//    Kokkos::deep_copy( NucleiLocation_G, NucleiLocation_H );
//    Kokkos::deep_copy( NucleationTimes_G, NucleationTimes_H );
//    
//    if (id == 0) cout << "Time spent initializing data = " << InitTime2-InitTime << " s" << endl;
//
//    
//    double TimeA = 0;
//    double TimeB = 0;
//    double TimeC = 0;
//    
//    for (int layernumber=0; layernumber<NumberOfLayers; layernumber++) {
//        cycle = 0;
//        int nn = 0; // Counter for the number of nucleation events
//        int XSwitch = 0;
//
//        double Time1, Time2, Time3, Time4;
//        // Loop continues until all liquid cells claimed by solid grains
//        do {
//            cycle++;
//
//            MPI_Barrier(MPI_COMM_WORLD);
//            if (id == 0) {
//                Time1 = MPI_Wtime();
//            }
//            
//            // Update cells on GPU - undercooling and diagonal length updates, nucleation
//            TemperatureUpdate(id, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, cycle, nn, AConst, BConst, CConst, CritTimeStep_G, CellType_G, UndercoolingCurrent_G, UndercoolingChange_G, NucleiLocation_G, NucleationTimes_G, GrainID_G, GrainOrientation, DOCenter_G, NeighborX,  NeighborY, NeighborZ, GrainUnitVector, TriangleIndex_G, CritDiagonalLength_G, DiagonalLength_G, NGrainOrientations, PossibleNuclei_ThisRank);
//            
//            MPI_Barrier(MPI_COMM_WORLD);
//            if (id == 0) {
//                Time2 = MPI_Wtime();
//                TimeA += (Time2-Time1);
//            }
//
//            // Update cells on GPU - new active cells, solidification of old active cells
//            CellCapture(id, cycle, DecompositionStrategy, MyXSlices, MyYSlices, nz, MyXOffset, MyYOffset, ItList, NeighborX, NeighborY, NeighborZ, GrainUnitVector, TriangleIndex_G, CritDiagonalLength_G, DiagonalLength_G, GrainOrientation, CellType_G, DOCenter_G, GrainID_G, NGrainOrientations);
//            
//            MPI_Barrier(MPI_COMM_WORLD);
//            if (id == 0) {
//                Time3 = MPI_Wtime();
//                TimeB += (Time3-Time2);
//            }
//            
//            // Use GPU vals on CPU
//            Kokkos::deep_copy( GrainID_H, GrainID_G );
//            Kokkos::deep_copy( CellType_H, CellType_G );
//            Kokkos::deep_copy( DiagonalLength_H, DiagonalLength_G );
//            Kokkos::deep_copy( CritDiagonalLength_H, CritDiagonalLength_G );
//            Kokkos::deep_copy( DOCenter_H, DOCenter_G );
//            Kokkos::deep_copy( TriangleIndex_H, TriangleIndex_G );
//            
//            // Update ghost nodes on host
//            if (DecompositionStrategy == 1) GhostNodes1D(cycle, id, MyLeft, MyRight, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, NeighborX, NeighborY, NeighborZ, CellType_H, DOCenter_H,GrainID_H, GrainUnitVector,TriangleIndex_H, GrainOrientation, DiagonalLength_H, CritDiagonalLength_H, NGrainOrientations);
//            else GhostNodes2D(cycle, id, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, NeighborX, NeighborY, NeighborZ, CellType_H, DOCenter_H, GrainID_H, GrainUnitVector, TriangleIndex_H, GrainOrientation, DiagonalLength_H, CritDiagonalLength_H, NGrainOrientations);
//
//            // Copy cell state and octahedron attribute data from CPU to GPU
//            Kokkos::deep_copy( GrainID_G, GrainID_H );
//            Kokkos::deep_copy( CellType_G, CellType_H );
//            Kokkos::deep_copy( DiagonalLength_G, DiagonalLength_H );
//            Kokkos::deep_copy( CritDiagonalLength_G, CritDiagonalLength_H );
//            Kokkos::deep_copy( DOCenter_G, DOCenter_H );
//            Kokkos::deep_copy( TriangleIndex_G, TriangleIndex_H );
//            
//            MPI_Barrier(MPI_COMM_WORLD);
//            if (id == 0) {
//                Time4 = MPI_Wtime();
//                TimeC += (Time4-Time3);
//            }
//
//            if (cycle % 500 == 0) {
//                IntermediateOutputAndCheck(id, cycle, MyXSlices, MyYSlices, nz, nn, XSwitch, CellType_H);
//            }
//
//            
//        } while(XSwitch == 0);
//
//
//        if (layernumber != NumberOfLayers-1) {
//            if (id == 0) cout << "Done with layer " << layernumber+1 << " of " << NumberOfLayers << endl;
//            LayerSetup(UseSubstrateFile, layernumber, LayerHeight, InitialGrainWidth, NGrainOrientations, DecompositionStrategy, ProcessorsInXDirection, ProcessorsInYDirection, nx, ny, nz, TemperatureDataType, MyXSlices, MyYSlices, MyXOffset, MyYOffset, id, np, ierr, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyLeftOut, MyRightIn, MyRightOut, ItList, NeighborX, NeighborY, NeighborZ, GrainOrientation, GrainUnitVector, DiagonalLength_H, CellType_H, TriangleIndex_H, GrainID_H, GrainID_Stored, CritDiagonalLength_H, DOCenter_H, CritTimeStep_H, UndercoolingChange_H, UndercoolingCurrent_H, tempfile, deltax, deltat, NMax, dTN, dTsigma, AFile1, AFile2, XMin, XMax, YMin, YMax, ZMin, ZMax, NextLayer_FirstNucleatedGrainID, Melted, MeltedStored, PossibleNuclei_ThisRank);
//            if (id == 0) cout << " Initialized ghost nodes for layer " << layernumber+2 << endl;
//            if (DecompositionStrategy == 1) GhostNodes1D(0, id,MyLeft, MyRight, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, NeighborX, NeighborY, NeighborZ, CellType_H, DOCenter_H ,GrainID_H , GrainUnitVector,TriangleIndex_H, GrainOrientation, DiagonalLength_H, CritDiagonalLength_H, NGrainOrientations);
//            else GhostNodes2D(0, id, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, NeighborX, NeighborY, NeighborZ, CellType_H, DOCenter_H, GrainID_H, GrainUnitVector, TriangleIndex_H, GrainOrientation, DiagonalLength_H, CritDiagonalLength_H, NGrainOrientations);
//             if (id == 0) cout << " Initialized nuclei for layer " << layernumber+2 << endl;
//            
//            //Kokkos::resize(a, 200,50);
//            
//            // Change nuclei back from active type, now that ghost nodes have been updated
//            NucleiInit(DecompositionStrategy, MyXSlices, MyYSlices, id, dTN, dTsigma, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, PossibleNuclei_ThisRank, NucleiLocation_H, NucleationTimes_H, GrainOrientation, CellType_H, GrainID_H, CritTimeStep_H, UndercoolingChange_H);
//            MPI_Barrier(MPI_COMM_WORLD);
//            if (id == 0) cout << "Starting layer " << layernumber+2 << endl;
//        }
//
//    }
//
//    double InitTime3 = MPI_Wtime();
//   // PrintTempValues(id,np,nx,ny,nz,MyXSlices,MyYSlices,ProcessorsInXDirection,ProcessorsInYDirection,CritTimeStep,UndercoolingChange,DecompositionStrategy);
//    
//    
//    if (NumberOfLayers == 1) {
//    PrintValues(id,np,nx,ny,nz,MyXSlices,MyYSlices,ProcessorsInXDirection,ProcessorsInYDirection,GrainID_H,GrainOrientation,GrainUnitVector,BaseFileName,DecompositionStrategy,NGrainOrientations,Melted);
//    }
//    else {
//        PrintValuesMultilayer(NumberOfLayers, LayerHeight, id,np,nx,ny,nz,MyXSlices,MyYSlices,ProcessorsInXDirection,ProcessorsInYDirection,GrainID_H,GrainID_Stored,GrainOrientation,GrainUnitVector,BaseFileName,DecompositionStrategy,NGrainOrientations,Melted,MeltedStored);
//    }
//    double InitTime4 = MPI_Wtime();
//    if (id == 0) {
//        cout << "===================================================================================" << endl;
//        cout << "Having run with = " << np << " processors" << endl;
//        cout << "Output written at cycle = " << cycle << endl;
//        cout << "Time to solution = " << InitTime3-InitTime << " s" << endl;
//        cout << "Time spent in TemperatureUpdate = " << TimeA << " s" << endl;
//        cout << "Time spent in CellCapture = " << TimeB << " s" << endl;
//        cout << "Time spent in GhostNodes = " << TimeC << " s" << endl;
//        cout << "Time spent collecting and printing output data = " << InitTime4-InitTime3 << " s" << endl;
//        cout << "===================================================================================" << endl;
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    cout << "ID = " << id << " ready to exit" << endl;
//    double GlobalT;
//    MPI_Reduce(&CoolTime,&GlobalT,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//    if (id == 0) cout << "Time spent in cooling = " << GlobalT << " s" << endl;
//    MPI_Reduce(&NucTime,&GlobalT,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//    if (id == 0) cout << "Time spent in nucleation = " << GlobalT << " s" << endl;
//    MPI_Reduce(&GrowTime,&GlobalT,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//    if (id == 0) cout << "Time spent in growth = " << GlobalT << " s" << endl;
//    MPI_Reduce(&CaptTime,&GlobalT,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//    if (id == 0) cout << "Time spent in capture = " << GlobalT << " s" << endl;
//    MPI_Reduce(&DeactTime,&GlobalT,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//    if (id == 0) cout << "Time spent in deactivation = " << GlobalT << " s" << endl;
    
}
