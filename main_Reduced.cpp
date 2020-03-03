#include "header.h"
using namespace std;

void RunProgram_Reduced(int id, int np, int ierr, int DecompositionStrategy, double deltax, double AConst, double BConst, double CConst, double DConst, double NMax, double dTN, double dTsigma, string BaseFileName, string GrainOrientationFile, string TemperatureDataType) {
    
    
    double InitTime = MPI_Wtime();
    int nx, ny, nz, NumberOfLayers, LayerHeight, TempFilesInSeries;
    double HT_deltax, deltat, FractSurfaceSitesActive, G, R;
    string SubstrateFileName, tempfile, TemperatureDataSource;
    bool LayerwiseTemperature;
    
    int DecompositionStrategy_A = DecompositionStrategy;
    
    if (TemperatureDataType == "C") {
        // Constrained alloy solidification test problem: fixed thermal gradient/solidification velocity
        CInputRead(G, R, nx, ny, nz, deltax, deltat, FractSurfaceSitesActive);
        NumberOfLayers = 1;
        LayerHeight = nz;
    }
    else {
        RInputRead(tempfile, HT_deltax, deltat, NumberOfLayers, LayerHeight, BaseFileName, TemperatureDataSource,SubstrateFileName,LayerwiseTemperature,TempFilesInSeries);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0) cout << "Input files read" << endl;
    
    // Grid decomposition
    int ProcessorsInXDirection, ProcessorsInYDirection;
    // Variables characterizing local processor grids relative to global domain
    int MyXSlices, MyXOffset, MyYSlices, MyYOffset, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyLeftOut, MyRightIn, MyRightOut;
    // Neighbor lists for cells
    int NeighborX[26], NeighborY[26], NeighborZ[26], ItList[9][26];
    float XMin, YMin, ZMin, XMax, YMax, ZMax; // OpenFOAM simulation bounds (if using OpenFOAM data)
    
    // Initialization of the grid and decomposition, along with deltax and deltat
    ParallelMeshInit(G, R, DecompositionStrategy_A, NeighborX, NeighborY, NeighborZ, ItList, TemperatureDataType, ierr, id, np, MyXSlices, MyYSlices, MyXOffset, MyYOffset, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyLeftOut, MyRightIn, MyRightOut, deltax, HT_deltax, deltat, nx, ny, nz, ProcessorsInXDirection, ProcessorsInYDirection, tempfile, XMin, XMax, YMin, YMax, ZMin, ZMax, TemperatureDataSource,LayerwiseTemperature);
    long int LocalDomainSize = MyXSlices*MyYSlices*nz; // Number of cells on this MPI rank
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0) cout << "Mesh initialized" << endl;
    
    // Temperature fields characterized by these variables:
     ViewI CritTimeStep_G("CritTimeStep_G", LocalDomainSize);
     ViewF UndercoolingChange_G("UndercoolingChange_G",LocalDomainSize);
     ViewF UndercoolingCurrent_G("UndercoolingCurrent_G",LocalDomainSize);
     ViewI::HostMirror CritTimeStep_H = Kokkos::create_mirror_view( CritTimeStep_G );
     ViewF::HostMirror UndercoolingChange_H = Kokkos::create_mirror_view( UndercoolingChange_G );
     ViewF::HostMirror UndercoolingCurrent_H = Kokkos::create_mirror_view( UndercoolingCurrent_G );
     bool* Melted = new bool[LocalDomainSize];
    
    // Initialize the temperature fields
    TempInit(-1, TempFilesInSeries, G, R, DecompositionStrategy_A,NeighborX, NeighborY, NeighborZ, ItList, TemperatureDataType, ierr, id, np, MyXSlices, MyYSlices, MyXOffset, MyYOffset, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyLeftOut, MyRightIn, MyRightOut, deltax, HT_deltax, deltat, nx, ny, nz, ProcessorsInXDirection, ProcessorsInYDirection,  CritTimeStep_H, UndercoolingChange_H, UndercoolingCurrent_H, tempfile, XMin, XMax, YMin, YMax, ZMin, ZMax, Melted, TemperatureDataSource, LayerwiseTemperature);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0) cout << "Done with temperature field initialization " << endl;
    
    int NGrainOrientations = 10000; // Number of grain orientations considered in the simulation
    float* GrainUnitVector = new float[18*NGrainOrientations];
    int* GrainOrientation = new int[NGrainOrientations];
    
    // Initialize grain orientations
    OrientationInit(id, NGrainOrientations, GrainOrientation, GrainUnitVector, GrainOrientationFile);
    
    // CA cell variables
    ViewI GrainID_G("GrainID_G",LocalDomainSize);
    ViewI CellType_G("CellType_G",LocalDomainSize);
    ViewI::HostMirror GrainID_H = Kokkos::create_mirror_view( GrainID_G );
    ViewI::HostMirror CellType_H = Kokkos::create_mirror_view( CellType_G );
    
    // Variables characterizing the active cell region within each rank's grid
    ViewF DiagonalLength_G("DiagonalLength_G",LocalDomainSize);
    ViewF CritDiagonalLength_G("CritDiagonalLength_G",26*LocalDomainSize);
    ViewF DOCenter_G("DOCenter_G",3*LocalDomainSize);
    ViewI TriangleIndex_G("TriangleIndex_G",26*3*LocalDomainSize);
    ViewF::HostMirror DiagonalLength_H = Kokkos::create_mirror_view( DiagonalLength_G );
    ViewF::HostMirror CritDiagonalLength_H = Kokkos::create_mirror_view( CritDiagonalLength_G );
    ViewF::HostMirror DOCenter_H = Kokkos::create_mirror_view( DOCenter_G );
    ViewI::HostMirror TriangleIndex_H = Kokkos::create_mirror_view( TriangleIndex_G );

    //if (id == 0) cout << "125*DMS = " << 125*LocalDomainSize << " " << " Local Domain Size = " << LocalDomainSize << " = " << (125*LocalDomainSize)/125 << endl;
    // Initialize the grain structure
    int PossibleNuclei_ThisRank, NextLayer_FirstNucleatedGrainID;
    GrainInit(-1, LayerHeight, TemperatureDataType, SubstrateFileName, FractSurfaceSitesActive, NGrainOrientations, DecompositionStrategy_A, ProcessorsInXDirection, ProcessorsInYDirection, nx, ny, nz, MyXSlices, MyYSlices, MyXOffset, MyYOffset, id, np, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, ItList, NeighborX, NeighborY, NeighborZ, GrainOrientation, GrainUnitVector, DiagonalLength_H, CellType_H, TriangleIndex_H, GrainID_H, CritDiagonalLength_H, DOCenter_H, CritTimeStep_H, UndercoolingChange_H, Melted, deltax, NMax, NextLayer_FirstNucleatedGrainID, PossibleNuclei_ThisRank);
   // MPI_Barrier(MPI_COMM_WORLD);
   // if (id == 0) cout << "GInit done" << endl;
    // Update ghost nodes for grain locations and attributes
    if (np > 1) {
        if (DecompositionStrategy_A == 1) GhostNodes1D(0, id, MyLeft, MyRight, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, NeighborX, NeighborY, NeighborZ, CellType_H, DOCenter_H, GrainID_H, GrainUnitVector,TriangleIndex_H, GrainOrientation, DiagonalLength_H, CritDiagonalLength_H, NGrainOrientations);
        else GhostNodes2D(0, id, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, NeighborX, NeighborY, NeighborZ, CellType_H, DOCenter_H, GrainID_H, GrainUnitVector, TriangleIndex_H, GrainOrientation, DiagonalLength_H, CritDiagonalLength_H, NGrainOrientations);
    }
  //  cout << "GN done, nuclei on rank " << id << " = " << PossibleNuclei_ThisRank << endl;
    ViewI NucleationTimes_G("NucleationTimes_G",PossibleNuclei_ThisRank);
    ViewI NucleiLocation_G("NucleiLocation_G",PossibleNuclei_ThisRank);
    ViewI::HostMirror NucleationTimes_H = Kokkos::create_mirror_view( NucleationTimes_G );
    ViewI::HostMirror NucleiLocation_H = Kokkos::create_mirror_view( NucleiLocation_G );
    
    // Update nuclei on ghost nodes, fill in nucleation data structures, and assign nucleation undercooling values to potential nucleation events
    NucleiInit(DecompositionStrategy_A, MyXSlices, MyYSlices, nz, id, dTN, dTsigma, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, PossibleNuclei_ThisRank, NucleiLocation_H, NucleationTimes_H, GrainOrientation, CellType_H, GrainID_H, CritTimeStep_H, UndercoolingChange_H);
    
   // MPI_Barrier(MPI_COMM_WORLD);
  //  if (id == 0) cout << "Done with grain initialization " << endl;
    
    // Normalize solidification parameters
    AConst = AConst*deltat/deltax;
    BConst = BConst*deltat/deltax;
    CConst = CConst*deltat/deltax;

    double InitTime2 = MPI_Wtime();
    
    int* GrainID_Stored = new int[MyXSlices*MyYSlices*((NumberOfLayers-1)*LayerHeight)];
    bool* MeltedStored = new bool[MyXSlices*MyYSlices*((NumberOfLayers-1)*LayerHeight)];
    int cycle;
    
    // Buffers for ghost node data (fixed size)
    int BufSizeX, BufSizeY, BufSizeZ;
    if (DecompositionStrategy_A == 1) {
        BufSizeX = MyXSlices;
        BufSizeY = 0;
        BufSizeZ = nz;
    }
    else {
        BufSizeX = MyXSlices-2;
        BufSizeY = MyYSlices-2;
        BufSizeZ = nz;
    }
    Buffer2D BufferA("BufferA",BufSizeX*BufSizeZ,5);
    Buffer2D BufferB("BufferB",BufSizeX*BufSizeZ,5);
    Buffer2D BufferC("BufferC",BufSizeY*BufSizeZ,5);
    Buffer2D BufferD("BufferD",BufSizeY*BufSizeZ,5);
    Buffer2D BufferE("BufferE",BufSizeZ,5);
    Buffer2D BufferF("BufferF",BufSizeZ,5);
    Buffer2D BufferG("BufferG",BufSizeZ,5);
    Buffer2D BufferH("BufferH",BufSizeZ,5);
    Buffer2D BufferAR("BufferAR",BufSizeX*BufSizeZ,5);
    Buffer2D BufferBR("BufferBR",BufSizeX*BufSizeZ,5);
    Buffer2D BufferCR("BufferCR",BufSizeY*BufSizeZ,5);
    Buffer2D BufferDR("BufferDR",BufSizeY*BufSizeZ,5);
    Buffer2D BufferER("BufferER",BufSizeZ,5);
    Buffer2D BufferFR("BufferFR",BufSizeZ,5);
    Buffer2D BufferGR("BufferGR",BufSizeZ,5);
    Buffer2D BufferHR("BufferHR",BufSizeZ,5);
    // Copy view data to GPU
    Kokkos::deep_copy( GrainID_G, GrainID_H );
    Kokkos::deep_copy( CellType_G, CellType_H );
    Kokkos::deep_copy( DiagonalLength_G, DiagonalLength_H );
    Kokkos::deep_copy( CritDiagonalLength_G, CritDiagonalLength_H );
    Kokkos::deep_copy( DOCenter_G, DOCenter_H );
    Kokkos::deep_copy( TriangleIndex_G, TriangleIndex_H );
    Kokkos::deep_copy( CritTimeStep_G, CritTimeStep_H );
    Kokkos::deep_copy( UndercoolingChange_G, UndercoolingChange_H );
    Kokkos::deep_copy( UndercoolingCurrent_G, UndercoolingCurrent_H );
    Kokkos::deep_copy( NucleiLocation_G, NucleiLocation_H );
    Kokkos::deep_copy( NucleationTimes_G, NucleationTimes_H );
    
    // Locks for cell capture
    // 0 = cannot be captured, 1 = can be capured
    ViewI Locks("Locks",LocalDomainSize);
    Kokkos::parallel_for("LockInit",LocalDomainSize, KOKKOS_LAMBDA (const int& i) {
        if ((CellType_G(i) == Delayed)||(CellType_G(i) == LiqSol)||(CellType_G(i) == Liquid)) Locks(i) = 1;
    });
    
    if (id == 0) cout << "Time spent initializing data = " << InitTime2-InitTime << " s" << endl;

    
    double TimeA = 0;
    double TimeB = 0;
    double TimeC = 0;
    
    for (int layernumber=0; layernumber<NumberOfLayers; layernumber++) {

        cycle = 0;
        int nn = 0; // Counter for the number of nucleation events
        int XSwitch = 0;

        double Time1, Time2, Time3, Time4;
        // Loop continues until all liquid cells claimed by solid grains
        do {
            cycle++;

//            MPI_Barrier(MPI_COMM_WORLD);
//            if (id == 0) {
//                Time1 = MPI_Wtime();
//            }
//            MPI_Barrier(MPI_COMM_WORLD);
            
            //if (id == 0) cout << " CYCLE " << cycle << endl;
            // Update cells on GPU - undercooling and diagonal length updates, nucleation
            Nucleation(id, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, cycle, nn, CritTimeStep_G, CellType_G, UndercoolingCurrent_G, UndercoolingChange_G, NucleiLocation_G, NucleationTimes_G, GrainID_G, GrainOrientation, DOCenter_G, NeighborX,  NeighborY, NeighborZ, GrainUnitVector, TriangleIndex_G, CritDiagonalLength_G, DiagonalLength_G, NGrainOrientations, PossibleNuclei_ThisRank);

//            if (id == 0) cout << "Cycle = " << cycle << endl;
//            MPI_Barrier(MPI_COMM_WORLD);
//            if (id == 0) {
//                Time2 = MPI_Wtime();
//                TimeA += (Time2-Time1);
//            }
//            MPI_Barrier(MPI_COMM_WORLD);
            
            //if (id == 0) cout << " CYCLE " << cycle << endl;
            // Update cells on GPU - new active cells, solidification of old active cells
            CellCapture(id, np, cycle, DecompositionStrategy_A, LocalDomainSize, MyXSlices, MyYSlices, nz, AConst, BConst, CConst, DConst, MyXOffset, MyYOffset, ItList, NeighborX, NeighborY, NeighborZ, CritTimeStep_G, UndercoolingCurrent_G, UndercoolingChange_G,  GrainUnitVector, TriangleIndex_G, CritDiagonalLength_G, DiagonalLength_G, GrainOrientation, CellType_G, DOCenter_G, GrainID_G, NGrainOrientations, BufferA, BufferB, BufferC, BufferD, BufferE, BufferF, BufferG, BufferH, BufSizeX, BufSizeY, Locks);
            
            //if (id == 0) cout << " CYCLE " << cycle << " starting ghost nodes" << endl;
//            MPI_Barrier(MPI_COMM_WORLD);
//            if (id == 0) {
//                Time3 = MPI_Wtime();
//                TimeB += (Time3-Time2);
//            }
//            MPI_Barrier(MPI_COMM_WORLD);

            if (np > 1) {
            // Update ghost nodes
                if (DecompositionStrategy_A == 1) GhostNodes1D_GPU(cycle, id, MyLeft, MyRight, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, NeighborX, NeighborY, NeighborZ, CellType_G, DOCenter_G,GrainID_G, GrainUnitVector,TriangleIndex_G, GrainOrientation, DiagonalLength_G, CritDiagonalLength_G, NGrainOrientations, BufferA, BufferB, BufferAR, BufferBR, BufSizeX,  BufSizeY, BufSizeZ, Locks);
                else GhostNodes2D_GPU(cycle, id, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, NeighborX, NeighborY, NeighborZ, CellType_G, DOCenter_G, GrainID_G, GrainUnitVector, TriangleIndex_G, GrainOrientation, DiagonalLength_G, CritDiagonalLength_G, NGrainOrientations, BufferA, BufferB, BufferC, BufferD, BufferE, BufferF, BufferG, BufferH, BufferAR, BufferBR, BufferCR, BufferDR, BufferER, BufferFR, BufferGR, BufferHR, BufSizeX, BufSizeY, BufSizeZ, Locks);
            }

//            MPI_Barrier(MPI_COMM_WORLD);
//            if (id == 0) {
//                Time4 = MPI_Wtime();
//                TimeC += (Time4-Time3);
//            }
//            MPI_Barrier(MPI_COMM_WORLD);
            
            if (cycle % 1000 == 0) {
                IntermediateOutputAndCheck(id, cycle, LocalDomainSize, nn, XSwitch, CellType_G, CritTimeStep_G, TemperatureDataType);
            }
            //if (cycle == 21000) XSwitch = 1;

        } while(XSwitch == 0);


        if (layernumber != NumberOfLayers-1) {

            // Use GPU vals on CPU for GrainID
            Kokkos::deep_copy( GrainID_H, GrainID_G );

            // Store CPU GrainID values from last layer, reset other variables on CPU for re-initialization
            LayerSetup(SubstrateFileName, layernumber, LayerHeight, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, LocalDomainSize, id, np, DiagonalLength_H, CellType_H, TriangleIndex_H, GrainID_H, GrainID_Stored, CritDiagonalLength_H, DOCenter_H, CritTimeStep_H, UndercoolingChange_H, UndercoolingCurrent_H, Melted, MeltedStored, LayerwiseTemperature, BufferA, BufferB, BufferC, BufferD, BufferE, BufferF, BufferG, BufferH, BufferAR, BufferBR, BufferCR, BufferDR, BufferER, BufferFR, BufferGR, BufferHR, BufSizeX, BufSizeY, BufSizeZ);
            
            // Initialize the temperature fields again, if temperature field is different each layer
            if (LayerwiseTemperature) TempInit(layernumber, TempFilesInSeries, G, R, DecompositionStrategy_A,NeighborX, NeighborY, NeighborZ, ItList, TemperatureDataType, ierr, id, np, MyXSlices, MyYSlices, MyXOffset, MyYOffset, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyLeftOut, MyRightIn, MyRightOut, deltax, HT_deltax, deltat, nx, ny, nz, ProcessorsInXDirection, ProcessorsInYDirection,  CritTimeStep_H, UndercoolingChange_H, UndercoolingCurrent_H, tempfile, XMin, XMax, YMin, YMax, ZMin, ZMax, Melted, TemperatureDataSource, LayerwiseTemperature);
            GrainInit(layernumber, LayerHeight, TemperatureDataType, SubstrateFileName, FractSurfaceSitesActive, NGrainOrientations, DecompositionStrategy_A, ProcessorsInXDirection, ProcessorsInYDirection, nx, ny, nz, MyXSlices, MyYSlices, MyXOffset, MyYOffset, id, np, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, ItList, NeighborX, NeighborY, NeighborZ, GrainOrientation, GrainUnitVector, DiagonalLength_H, CellType_H, TriangleIndex_H, GrainID_H, CritDiagonalLength_H, DOCenter_H, CritTimeStep_H, UndercoolingChange_H, Melted, deltax, NMax, NextLayer_FirstNucleatedGrainID, PossibleNuclei_ThisRank);

            // Resize data structures for nucleated grains (on CPU) - this new layer may have a different number of nucleation events
            Kokkos::resize(NucleationTimes_H,PossibleNuclei_ThisRank);
            Kokkos::resize(NucleiLocation_H,PossibleNuclei_ThisRank);
            // Resize data structures for nucleated grains (on GPU) - this new layer may have a different number of nucleation events
            Kokkos::resize(NucleationTimes_G,PossibleNuclei_ThisRank);
            Kokkos::resize(NucleiLocation_G,PossibleNuclei_ThisRank);
            
            if (id == 0) cout << " Initialized ghost nodes for layer " << layernumber+2 << endl;
            // Update ghost nodes for grain locations and attributes
            if (np > 1) {
                if (DecompositionStrategy_A == 1) GhostNodes1D(0, id, MyLeft, MyRight, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, NeighborX, NeighborY, NeighborZ, CellType_H, DOCenter_H, GrainID_H, GrainUnitVector,TriangleIndex_H, GrainOrientation, DiagonalLength_H, CritDiagonalLength_H, NGrainOrientations);
                else GhostNodes2D(0, id, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, NeighborX, NeighborY, NeighborZ, CellType_H, DOCenter_H, GrainID_H, GrainUnitVector, TriangleIndex_H, GrainOrientation, DiagonalLength_H, CritDiagonalLength_H, NGrainOrientations);
            }
             if (id == 0) cout << " Initialized nuclei for layer " << layernumber+2 << endl;

            // Change nuclei back from active type, now that ghost nodes have been updated
            NucleiInit(DecompositionStrategy_A, MyXSlices, MyYSlices, nz, id, dTN, dTsigma, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, PossibleNuclei_ThisRank, NucleiLocation_H, NucleationTimes_H, GrainOrientation, CellType_H, GrainID_H, CritTimeStep_H, UndercoolingChange_H);

            // Copy cell state and octahedron attribute data from CPU to GPU
            Kokkos::deep_copy( GrainID_G, GrainID_H );
            Kokkos::deep_copy( CellType_G, CellType_H );
            Kokkos::deep_copy( DiagonalLength_G, DiagonalLength_H );
            Kokkos::deep_copy( CritDiagonalLength_G, CritDiagonalLength_H );
            Kokkos::deep_copy( DOCenter_G, DOCenter_H );
            Kokkos::deep_copy( TriangleIndex_G, TriangleIndex_H );
            Kokkos::deep_copy( CritTimeStep_G, CritTimeStep_H );
            Kokkos::deep_copy( UndercoolingChange_G, UndercoolingChange_H );
            Kokkos::deep_copy( UndercoolingCurrent_G, UndercoolingCurrent_H );
            // Deep copy new layer nucleation data to GPU
            Kokkos::deep_copy( NucleationTimes_G, NucleationTimes_H );
            Kokkos::deep_copy( NucleiLocation_G, NucleiLocation_H );
            
            // Reinitialize locks
            Kokkos::parallel_for("LockInit",LocalDomainSize, KOKKOS_LAMBDA (const int& i) {
                if ((CellType_G(i) == Delayed)||(CellType_G(i) == LiqSol)||(CellType_G(i) == Liquid)) Locks(i) = 1;
            });
            
            MPI_Barrier(MPI_COMM_WORLD);
            if (id == 0) cout << "Starting layer " << layernumber+2 << endl;
        }

    }

    double InitTime3 = MPI_Wtime();
    
    // Copy GPU results for GrainID back to CPU for printing to file(s)
    Kokkos::deep_copy( GrainID_H, GrainID_G );
    Kokkos::deep_copy( CellType_H, CellType_G );
    //PrintCT( id,  np,  nx,  ny,  nz,  MyXSlices,  MyYSlices,  ProcessorsInXDirection,  ProcessorsInYDirection,  CellType_H,  BaseFileName,  DecompositionStrategy_A);
    //PrintTempValues(id,np,nx,ny,nz, MyXSlices, MyYSlices, ProcessorsInXDirection, ProcessorsInYDirection, CritTimeStep, UndercoolingChange, DecompositionStrategy_A);
    
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (id == 0) cout << "Printing to files" << endl;
    if (NumberOfLayers == 1) {
    PrintValues(id,np,nx,ny,nz,MyXSlices,MyYSlices, MyXOffset, MyYOffset, ProcessorsInXDirection, ProcessorsInYDirection, GrainID_H,GrainOrientation,GrainUnitVector,BaseFileName,DecompositionStrategy_A,NGrainOrientations,Melted);
    }
    else {
        PrintValuesMultilayer(NumberOfLayers, LayerHeight, id,np,nx,ny,nz,MyXSlices,MyYSlices,ProcessorsInXDirection,ProcessorsInYDirection,GrainID_H,GrainID_Stored,GrainOrientation,GrainUnitVector,BaseFileName,DecompositionStrategy_A,NGrainOrientations,Melted,MeltedStored);
    }
    double InitTime4 = MPI_Wtime();
    if (id == 0) {
        cout << "===================================================================================" << endl;
        cout << "Having run with = " << np << " processors" << endl;
        cout << "Output written at cycle = " << cycle << endl;
        cout << "Total time = " << InitTime4 - InitTime << endl;
//        cout << "Time spent initializing data = " << InitTime2-InitTime << " s" << endl;
//        cout << "Time spent in Nucleation = " << TimeA << " s" << endl;
//        cout << "Time spent in CellCapture = " << TimeB << " s" << endl;
//        cout << "Time spent in GhostNodes = " << TimeC << " s" << endl;
//        cout << "Time spent collecting and printing output data = " << InitTime4-InitTime3 << " s" << endl;
        cout << "===================================================================================" << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
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
