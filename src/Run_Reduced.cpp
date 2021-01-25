#include "header.h"
using namespace std;

void RunProgram_Reduced(int id, int np, int ierr, string InputFile) {
    
    double NuclTime = 0.0, CaptureTime = 0.0, GhostTime = 0.0;
    double StartNuclTime, StartCaptureTime, StartGhostTime;
    double StartTime = MPI_Wtime();
    
    int nx, ny, nz, DecompositionStrategy, NumberOfLayers, LayerHeight, TempFilesInSeries, NumberOfTruchasRanks;
    bool TruchasMultilayer = false; // If reading from Truchas multilayer data, this is changed to true in the input file (previously, this variable was called "BurstBuffer")
    bool PrintFilesYN;
    bool FilesToPrint[6] = {0}; // Which specific files to print are specified in the input file
    double HT_deltax, deltax, deltat, FractSurfaceSitesActive, G, R, AConst, BConst, CConst, DConst, FreezingRange, NMax, dTN, dTsigma;
    string SubstrateFileName, tempfile, SimulationType, OutputFile, GrainOrientationFile, TemperatureDataSource, ExtraWalls, PathToOutput;
    
    // Read input data
    InputReadFromFile(id, InputFile, SimulationType, DecompositionStrategy, AConst, BConst, CConst, DConst, FreezingRange, deltax, NMax, dTN, dTsigma, OutputFile, GrainOrientationFile, tempfile, TempFilesInSeries, TruchasMultilayer, ExtraWalls, HT_deltax, TemperatureDataSource, deltat, NumberOfLayers, LayerHeight, SubstrateFileName, G, R, nx, ny, nz, FractSurfaceSitesActive,PathToOutput,NumberOfTruchasRanks,FilesToPrint,PrintFilesYN);
    
    // Grid decomposition
    int ProcessorsInXDirection, ProcessorsInYDirection;
    // Variables characterizing local processor grids relative to global domain
    int MyXSlices, MyXOffset, MyYSlices, MyYOffset, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyLeftOut, MyRightIn, MyRightOut;
    // Neighbor lists for cells
    ViewI_H NeighborX_H(Kokkos::ViewAllocateWithoutInitializing("NeighborX"),26);
    ViewI_H NeighborY_H(Kokkos::ViewAllocateWithoutInitializing("NeighborY"),26);
    ViewI_H NeighborZ_H(Kokkos::ViewAllocateWithoutInitializing("NeighborZ"),26);
    ViewI2D_H ItList_H(Kokkos::ViewAllocateWithoutInitializing("ItList"),9,26);
    float XMin, YMin, ZMin, XMax, YMax, ZMax; // OpenFOAM simulation bounds (if using OpenFOAM data)
    float* ZMinLayer = new float[NumberOfLayers];
    float* ZMaxLayer = new float[NumberOfLayers];
    int* FinishTimeStep = new int[NumberOfLayers];
    
    // Temporary data structure for storing temperature data from file(s)
    // Initial estimate for size
    vector <float> RawData(1000000);
    // Contains "NumberOfLayers" values corresponding to the location within "RawData" of the first data element in each layer
    int* FirstValue = new int[NumberOfLayers];

    // Initialization of the grid and decomposition, along with deltax and deltat
    // Read in temperature data
    ParallelMeshInit(DecompositionStrategy, NeighborX_H, NeighborY_H, NeighborZ_H, ItList_H, SimulationType, id, np, MyXSlices, MyYSlices, MyXOffset, MyYOffset, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyLeftOut, MyRightIn, MyRightOut, deltax, nx, ny, nz, ProcessorsInXDirection, ProcessorsInYDirection, tempfile, XMin, XMax, YMin, YMax, ZMin, ZMax, TemperatureDataSource, LayerHeight, NumberOfLayers, TempFilesInSeries, ZMinLayer, ZMaxLayer, FirstValue, RawData, TruchasMultilayer,NumberOfTruchasRanks);
    
    long int LocalDomainSize = MyXSlices*MyYSlices*nz; // Number of cells on this MPI rank

    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0) cout << "Mesh initialized" << endl;
    
    // Temperature fields characterized by these variables:
     ViewI_H CritTimeStep_H(Kokkos::ViewAllocateWithoutInitializing("CritTimeStep"), LocalDomainSize);
     ViewI_H LayerID_H(Kokkos::ViewAllocateWithoutInitializing("LayerID"), LocalDomainSize);
     ViewF_H UndercoolingChange_H(Kokkos::ViewAllocateWithoutInitializing("UndercoolingChange"), LocalDomainSize);
     ViewF_H UndercoolingCurrent_H(Kokkos::ViewAllocateWithoutInitializing("UndercoolingCurrent"), LocalDomainSize);
     bool* Melted = new bool[LocalDomainSize];

    // By default, the active domain bounds are the same as the global domain bounds
    // For multilayer problems, this is not the case and ZBound_High and nzActive will be adjusted in TempInit to account for only the first layer of solidification
    int ZBound_Low;
    int ZBound_High;
    int nzActive;
    
    // Initialize the temperature fields
    TempInit(-1, TempFilesInSeries, G, R, SimulationType, ierr, id, MyXSlices, MyYSlices, MyXOffset, MyYOffset, deltax, HT_deltax, deltat, nx, ny, nz,  CritTimeStep_H, UndercoolingChange_H, UndercoolingCurrent_H, XMin, YMin, ZMin, Melted, ZMinLayer, ZMaxLayer, LayerHeight, NumberOfLayers, nzActive, ZBound_Low, ZBound_High, FinishTimeStep, FreezingRange, LayerID_H, FirstValue, RawData, TruchasMultilayer);
    // Delete temporary data structure for temperature data read
    RawData.clear();
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0) cout << "Done with temperature field initialization, active domain size is " << nzActive << " out of " << nz << " cells in the Z direction" << endl;

    int LocalActiveDomainSize = MyXSlices*MyYSlices*nzActive; // Number of active cells on this MPI rank
    
    // PrintTempValues(id,np,nx,ny,nz, MyXSlices, MyYSlices, ProcessorsInXDirection, ProcessorsInYDirection, CritTimeStep_H, UndercoolingChange_H, DecompositionStrategy,PathToOutput);
    
    int NGrainOrientations = 10000; // Number of grain orientations considered in the simulation
    ViewF_H GrainUnitVector_H(Kokkos::ViewAllocateWithoutInitializing("GrainUnitVector"), 9*NGrainOrientations);
    ViewI_H GrainOrientation_H(Kokkos::ViewAllocateWithoutInitializing("GrainOrientation"), NGrainOrientations);
    
    // Initialize grain orientations
    OrientationInit(id, NGrainOrientations, GrainOrientation_H, GrainUnitVector_H, GrainOrientationFile);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0) cout << "Done with orientation initialization " << endl;
    
    // CA cell variables
    ViewI_H GrainID_H(Kokkos::ViewAllocateWithoutInitializing("GrainID"),LocalDomainSize);
    ViewI_H CellType_H(Kokkos::ViewAllocateWithoutInitializing("CellType"),LocalDomainSize);

    // Variables characterizing the active cell region within each rank's grid
    ViewF_H DiagonalLength_H(Kokkos::ViewAllocateWithoutInitializing("DiagonalLength"),LocalActiveDomainSize);
    ViewF_H CritDiagonalLength_H(Kokkos::ViewAllocateWithoutInitializing("CritDiagonalLength"),26*LocalActiveDomainSize);
    ViewF_H DOCenter_H(Kokkos::ViewAllocateWithoutInitializing("DOCenter"),3*LocalActiveDomainSize);

    // Initialize the grain structure
    int PossibleNuclei_ThisRank, NextLayer_FirstNucleatedGrainID;
    
    GrainInit(-1, SimulationType, SubstrateFileName, FractSurfaceSitesActive, NGrainOrientations, DecompositionStrategy, nx, ny, nz, MyXSlices, MyYSlices, MyXOffset, MyYOffset, id, np, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, ItList_H, NeighborX_H, NeighborY_H, NeighborZ_H, GrainOrientation_H, GrainUnitVector_H, DiagonalLength_H, CellType_H, GrainID_H, CritDiagonalLength_H, DOCenter_H, CritTimeStep_H, UndercoolingChange_H, Melted, deltax, NMax, NextLayer_FirstNucleatedGrainID, PossibleNuclei_ThisRank, ZBound_High, ZBound_Low, ExtraWalls);
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0) cout << "Grain struct initialized" << endl;

    ViewI_H NucleationTimes_H(Kokkos::ViewAllocateWithoutInitializing("NucleationTimes"),PossibleNuclei_ThisRank);
    ViewI_H NucleiLocation_H(Kokkos::ViewAllocateWithoutInitializing("NucleiLocation"),PossibleNuclei_ThisRank);
    
    // Update nuclei on ghost nodes, fill in nucleation data structures, and assign nucleation undercooling values to potential nucleation events
    if (id == 0) cout << " Possible nucleation events (rank: # events): " << endl;
    NucleiInit(DecompositionStrategy, MyXSlices, MyYSlices, nz, id, dTN, dTsigma, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, NucleiLocation_H, NucleationTimes_H, CellType_H, GrainID_H, CritTimeStep_H, UndercoolingChange_H);

    MPI_Barrier(MPI_COMM_WORLD);
    
    // Normalize solidification parameters
    AConst = AConst*deltat/deltax;
    BConst = BConst*deltat/deltax;
    CConst = CConst*deltat/deltax;
    int cycle;
    
    // Buffers for ghost node data (fixed size)
    int BufSizeX, BufSizeY, BufSizeZ;
    if (DecompositionStrategy == 1) {
        BufSizeX = MyXSlices;
        BufSizeY = 0;
        BufSizeZ = nzActive;
    }
    else {
        BufSizeX = MyXSlices-2;
        BufSizeY = MyYSlices-2;
        BufSizeZ = nzActive;
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
    using memory_space = Kokkos::DefaultExecutionSpace::memory_space;
    ViewI GrainID_G = Kokkos::create_mirror_view_and_copy( memory_space(), GrainID_H );
    ViewI CellType_G = Kokkos::create_mirror_view_and_copy( memory_space(), CellType_H );
    ViewF DiagonalLength_G = Kokkos::create_mirror_view_and_copy( memory_space(), DiagonalLength_H );
    ViewF CritDiagonalLength_G = Kokkos::create_mirror_view_and_copy( memory_space(), CritDiagonalLength_H );
    ViewF DOCenter_G = Kokkos::create_mirror_view_and_copy( memory_space(), DOCenter_H );
    ViewI CritTimeStep_G = Kokkos::create_mirror_view_and_copy( memory_space(), CritTimeStep_H );
    ViewI LayerID_G = Kokkos::create_mirror_view_and_copy( memory_space(), LayerID_H );
    ViewF UndercoolingChange_G = Kokkos::create_mirror_view_and_copy( memory_space(), UndercoolingChange_H );
    ViewF UndercoolingCurrent_G = Kokkos::create_mirror_view_and_copy( memory_space(), UndercoolingCurrent_H );
    ViewI NucleiLocation_G = Kokkos::create_mirror_view_and_copy( memory_space(), NucleiLocation_H );
    ViewI NucleationTimes_G = Kokkos::create_mirror_view_and_copy( memory_space(), NucleationTimes_H );
    ViewI NeighborX_G = Kokkos::create_mirror_view_and_copy( memory_space(), NeighborX_H );
    ViewI NeighborY_G = Kokkos::create_mirror_view_and_copy( memory_space(), NeighborY_H );
    ViewI NeighborZ_G = Kokkos::create_mirror_view_and_copy( memory_space(), NeighborZ_H );
    ViewI2D ItList_G = Kokkos::create_mirror_view_and_copy( memory_space(), ItList_H );
    ViewI GrainOrientation_G = Kokkos::create_mirror_view_and_copy( memory_space(), GrainOrientation_H );
    ViewF GrainUnitVector_G = Kokkos::create_mirror_view_and_copy( memory_space(), GrainUnitVector_H );
    
    // Locks for cell capture
    // 0 = cannot be captured, 1 = can be capured
    ViewI Locks(Kokkos::ViewAllocateWithoutInitializing("Locks"),LocalActiveDomainSize);
    Kokkos::parallel_for("LockInit",LocalActiveDomainSize, KOKKOS_LAMBDA (const int& D3D1ConvPosition) {
        int RankZ = D3D1ConvPosition/(MyXSlices*MyYSlices);
        int Rem = D3D1ConvPosition % (MyXSlices*MyYSlices);
        int RankX = Rem/MyYSlices;
        int RankY = Rem % MyYSlices;
        int GlobalZ = ZBound_Low + RankZ;
        int GlobalD3D1ConvPosition = GlobalZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
        if ((CellType_G(GlobalD3D1ConvPosition) == Delayed)||(CellType_G(GlobalD3D1ConvPosition) == LiqSol)||(CellType_G(GlobalD3D1ConvPosition) == Liquid)) Locks(D3D1ConvPosition) = 1;
        else Locks(D3D1ConvPosition) = 0;
    });
    
    if (np > 1) {
        // Ghost nodes for initial microstructure state
        GhostNodesInit_GPU(DecompositionStrategy, MyXSlices, MyYSlices, GrainID_G, DOCenter_G, DiagonalLength_G, BufferA, BufferB, BufferC, BufferD, BufferE, BufferF, BufferG, BufferH, BufSizeX, BufSizeY, LocalActiveDomainSize, ZBound_Low);
        if (DecompositionStrategy == 1) GhostNodes1D_GPU(0, id, MyLeft, MyRight, MyXSlices, MyYSlices, MyXOffset, MyYOffset, NeighborX_G, NeighborY_G, NeighborZ_G, CellType_G, DOCenter_G,GrainID_G, GrainUnitVector_G, GrainOrientation_G, DiagonalLength_G, CritDiagonalLength_G, NGrainOrientations, BufferA, BufferB, BufferAR, BufferBR, BufSizeX,  BufSizeY, BufSizeZ, Locks, ZBound_Low);
        else GhostNodes2D_GPU(0, id, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, MyXSlices, MyYSlices, MyXOffset, MyYOffset, NeighborX_G, NeighborY_G, NeighborZ_G, CellType_G, DOCenter_G, GrainID_G, GrainUnitVector_G, GrainOrientation_G, DiagonalLength_G, CritDiagonalLength_G, NGrainOrientations, BufferA, BufferB, BufferC, BufferD, BufferE, BufferF, BufferG, BufferH, BufferAR, BufferBR, BufferCR, BufferDR, BufferER, BufferFR, BufferGR, BufferHR, BufSizeX, BufSizeY, BufSizeZ, Locks, ZBound_Low);
    }

    double InitTime = MPI_Wtime() - StartTime;
    if (id == 0) cout << "\nData initialized: Time spent: " << InitTime << " s" << endl;

    cycle = 0;
    
    for (int layernumber=0; layernumber<NumberOfLayers; layernumber++) {
        
        int nn = 0; // Counter for the number of nucleation events
        int XSwitch = 0;
        double LayerTime1 = MPI_Wtime();
        
        // Loop continues until all liquid cells claimed by solid grains
        do {
            cycle++;
            // Update cells on GPU - undercooling and diagonal length updates, nucleation
            StartNuclTime = MPI_Wtime();
            Nucleation(id, MyXSlices, MyYSlices, MyXOffset, MyYOffset, cycle, nn, CellType_G, NucleiLocation_G, NucleationTimes_G, GrainID_G, GrainOrientation_G, DOCenter_G, NeighborX_G,  NeighborY_G, NeighborZ_G, GrainUnitVector_G, CritDiagonalLength_G, DiagonalLength_G, NGrainOrientations, PossibleNuclei_ThisRank, Locks, ZBound_Low, layernumber, LayerID_G);
            NuclTime += MPI_Wtime() - StartNuclTime;

            // Update cells on GPU - new active cells, solidification of old active cells
            StartCaptureTime = MPI_Wtime();
            CellCapture(id, np, cycle, DecompositionStrategy, LocalActiveDomainSize, MyXSlices, MyYSlices, AConst, BConst, CConst, DConst, MyXOffset, MyYOffset, ItList_G, NeighborX_G, NeighborY_G, NeighborZ_G, CritTimeStep_G, UndercoolingCurrent_G, UndercoolingChange_G,  GrainUnitVector_G, CritDiagonalLength_G, DiagonalLength_G, GrainOrientation_G, CellType_G, DOCenter_G, GrainID_G, NGrainOrientations, BufferA, BufferB, BufferC, BufferD, BufferE, BufferF, BufferG, BufferH, BufSizeX, BufSizeY, Locks, ZBound_Low, nzActive, layernumber, LayerID_G);
            CaptureTime += MPI_Wtime() - StartCaptureTime;

            if (np > 1) {
                // Update ghost nodes
                StartGhostTime = MPI_Wtime();
                if (DecompositionStrategy == 1) GhostNodes1D_GPU(cycle, id, MyLeft, MyRight, MyXSlices, MyYSlices, MyXOffset, MyYOffset, NeighborX_G, NeighborY_G, NeighborZ_G, CellType_G, DOCenter_G,GrainID_G, GrainUnitVector_G, GrainOrientation_G, DiagonalLength_G, CritDiagonalLength_G, NGrainOrientations, BufferA, BufferB, BufferAR, BufferBR, BufSizeX,  BufSizeY, BufSizeZ, Locks, ZBound_Low);
                else GhostNodes2D_GPU(cycle, id, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, MyXSlices, MyYSlices, MyXOffset, MyYOffset, NeighborX_G, NeighborY_G, NeighborZ_G, CellType_G, DOCenter_G, GrainID_G, GrainUnitVector_G, GrainOrientation_G, DiagonalLength_G, CritDiagonalLength_G, NGrainOrientations, BufferA, BufferB, BufferC, BufferD, BufferE, BufferF, BufferG, BufferH, BufferAR, BufferBR, BufferCR, BufferDR, BufferER, BufferFR, BufferGR, BufferHR, BufSizeX, BufSizeY, BufSizeZ, Locks, ZBound_Low);
                GhostTime += MPI_Wtime() - StartGhostTime;
            }

            if (cycle % 1000 == 0) {
                IntermediateOutputAndCheck(id, cycle, MyXSlices, MyYSlices, LocalDomainSize, LocalActiveDomainSize, nn, XSwitch, CellType_G, CritTimeStep_G, SimulationType, FinishTimeStep, layernumber, NumberOfLayers, ZBound_Low, LayerID_G );
            }
            
        } while(XSwitch == 0);

        if (layernumber != NumberOfLayers-1) {
             // Determine new active cell domain size and offset from bottom of global domain
             int ZShift;
             DomainShiftAndResize(id, MyXSlices, MyYSlices, ZShift, ZBound_Low, ZBound_High, nzActive, LocalDomainSize, LocalActiveDomainSize, BufSizeZ, LayerHeight, CellType_G, layernumber, LayerID_G);

             // Resize active cell data structures
             Kokkos::resize(DiagonalLength_G,LocalActiveDomainSize);
             Kokkos::resize(DOCenter_G,3*LocalActiveDomainSize);
             Kokkos::resize(CritDiagonalLength_G,26*LocalActiveDomainSize);
             Kokkos::resize(Locks,LocalActiveDomainSize);
            
             Kokkos::resize(BufferA,BufSizeX*BufSizeZ,5);
             Kokkos::resize(BufferB,BufSizeX*BufSizeZ,5);
             Kokkos::resize(BufferC,BufSizeY*BufSizeZ,5);
             Kokkos::resize(BufferD,BufSizeY*BufSizeZ,5);
             Kokkos::resize(BufferE,BufSizeZ,5);
             Kokkos::resize(BufferF,BufSizeZ,5);
             Kokkos::resize(BufferG,BufSizeZ,5);
             Kokkos::resize(BufferH,BufSizeZ,5);
            
             Kokkos::resize(BufferAR,BufSizeX*BufSizeZ,5);
             Kokkos::resize(BufferBR,BufSizeX*BufSizeZ,5);
             Kokkos::resize(BufferCR,BufSizeY*BufSizeZ,5);
             Kokkos::resize(BufferDR,BufSizeY*BufSizeZ,5);
             Kokkos::resize(BufferER,BufSizeZ,5);
             Kokkos::resize(BufferFR,BufSizeZ,5);
             Kokkos::resize(BufferGR,BufSizeZ,5);
             Kokkos::resize(BufferHR,BufSizeZ,5);
            
             MPI_Barrier(MPI_COMM_WORLD);
             if (id == 0) cout << "Resize executed" << endl;
            
             // Update active cell data structures for simulation of next layer
             LayerSetup(MyXSlices, MyYSlices, MyXOffset, MyYOffset, LocalActiveDomainSize, GrainOrientation_G, NGrainOrientations, GrainUnitVector_G, NeighborX_G, NeighborY_G, NeighborZ_G, DiagonalLength_G, CellType_G, GrainID_G, CritDiagonalLength_G, DOCenter_G, BufferA, BufferB, BufferC, BufferD, BufferE, BufferF, BufferG, BufferH, BufferAR, BufferBR, BufferCR, BufferDR, BufferER, BufferFR, BufferGR, BufferHR, BufSizeX, BufSizeY, BufSizeZ, ZBound_Low, Locks);

            if (id == 0) cout << "New layer setup, GN dimensions are " << BufSizeX << " " << BufSizeY << " " << BufSizeZ << endl;
            // Update ghost nodes for grain locations and attributes
            GhostNodesInit_GPU(DecompositionStrategy, MyXSlices, MyYSlices, GrainID_G, DOCenter_G, DiagonalLength_G, BufferA, BufferB, BufferC, BufferD, BufferE, BufferF, BufferG, BufferH, BufSizeX, BufSizeY, LocalActiveDomainSize, ZBound_Low);
            MPI_Barrier(MPI_COMM_WORLD);
            if (id == 0) cout << "New layer ghost nodes initialized" << endl;
             if (np > 1) {
             // Update ghost nodes
                 if (DecompositionStrategy == 1) GhostNodes1D_GPU(cycle, id, MyLeft, MyRight, MyXSlices, MyYSlices, MyXOffset, MyYOffset, NeighborX_G, NeighborY_G, NeighborZ_G, CellType_G, DOCenter_G,GrainID_G, GrainUnitVector_G, GrainOrientation_G, DiagonalLength_G, CritDiagonalLength_G, NGrainOrientations, BufferA, BufferB, BufferAR, BufferBR, BufSizeX,  BufSizeY, BufSizeZ, Locks, ZBound_Low);
                 else GhostNodes2D_GPU(cycle, id, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, MyXSlices, MyYSlices, MyXOffset, MyYOffset, NeighborX_G, NeighborY_G, NeighborZ_G, CellType_G, DOCenter_G, GrainID_G, GrainUnitVector_G, GrainOrientation_G, DiagonalLength_G, CritDiagonalLength_G, NGrainOrientations, BufferA, BufferB, BufferC, BufferD, BufferE, BufferF, BufferG, BufferH, BufferAR, BufferBR, BufferCR, BufferDR, BufferER, BufferFR, BufferGR, BufferHR, BufSizeX, BufSizeY, BufSizeZ, Locks, ZBound_Low);
             }
             XSwitch = 0;
             MPI_Barrier(MPI_COMM_WORLD);
             double LayerTime2 = MPI_Wtime();
             cycle = 0;
             if (id == 0) cout << "Time for layer number " << layernumber << " was " << LayerTime2-LayerTime1 << " s, starting layer " << layernumber+1 << endl;
         }
         else {
            MPI_Barrier(MPI_COMM_WORLD);
            double LayerTime2 = MPI_Wtime();
            if (id == 0) cout << "Time for final layer was " << LayerTime2-LayerTime1 << " s" << endl;
         }
    }

    double RunTime = MPI_Wtime() - InitTime;

    // Copy GPU results for GrainID back to CPU for printing to file(s)
    Kokkos::deep_copy( GrainID_H, GrainID_G );
    Kokkos::deep_copy( CellType_H, CellType_G );

    MPI_Barrier(MPI_COMM_WORLD);
    if (PrintFilesYN) {
        if (id == 0) cout << "Collecting data on rank 0 and printing to files" << endl;
        CollectGrainData(id,np,nx,ny,nz,MyXSlices,MyYSlices, ProcessorsInXDirection, ProcessorsInYDirection, GrainID_H,GrainOrientation_H,GrainUnitVector_H,OutputFile,DecompositionStrategy,NGrainOrientations,Melted,PathToOutput,FilesToPrint,deltax);
    }
    else {
        if (id == 0) cout << "No output files to be printed, exiting program" << endl;
    }
    
    double OutTime = MPI_Wtime() - RunTime - InitTime;
    double InitMaxTime, InitMinTime, OutMaxTime, OutMinTime = 0.0;
    double NuclMaxTime, NuclMinTime, CaptureMaxTime, CaptureMinTime, GhostMaxTime, GhostMinTime = 0.0;
    MPI_Allreduce( &InitTime, &InitMaxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
    MPI_Allreduce( &InitTime, &InitMinTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
    MPI_Allreduce( &NuclTime, &NuclMaxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
    MPI_Allreduce( &NuclTime, &NuclMinTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
    MPI_Allreduce( &CaptureTime, &CaptureMaxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
    MPI_Allreduce( &CaptureTime, &CaptureMinTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
    MPI_Allreduce( &GhostTime, &GhostMaxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
    MPI_Allreduce( &GhostTime, &GhostMinTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
    MPI_Allreduce( &OutTime, &OutMaxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
    MPI_Allreduce( &OutTime, &OutMinTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );

    if (id == 0) {
        cout << "===================================================================================" << endl;
        cout << "Having run with = " << np << " processors" << endl;
        cout << "Output written at cycle = " << cycle << endl;
        cout << "Total time = " << InitTime + RunTime + OutTime << endl;
        cout << "Time spent initializing data = " << InitTime << " s" << endl;
        cout << "Time spent performing CA calculations = " << RunTime << " s" << endl;
        cout << "Time spent collecting and printing output data = " << OutTime << " s\n" << endl;

        cout << "Max/min rank time initializing data  = " << InitMaxTime << " / " << InitMinTime <<" s" << endl;
        cout << "Max/min rank time in CA nucleation   = " << NuclMaxTime << " / " << NuclMinTime <<" s" << endl;
        cout << "Max/min rank time in CA cell capture = " << CaptureMaxTime << " / " << CaptureMinTime << " s" << endl;
        cout << "Max/min rank time in CA ghosting     = " << GhostMaxTime << " / " << GhostMinTime << " s" << endl;
        cout << "Max/min rank time exporting data     = " << OutMaxTime << " / " << OutMinTime << " s" << endl << endl;

        cout << "===================================================================================" << endl;
    }
}
