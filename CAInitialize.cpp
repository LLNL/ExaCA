#include "header.h"
using namespace std;
// Initializes input parameters, mesh, temperature field, and grain structures for CA simulations

void MasterInputRead(int &DecompositionStrategy, double &deltax, double &AConst, double &BConst, double &CConst, double &DConst, double &NMax, double &dTN, double &dTsigma, string &BaseFileName, string &GrainOrientationFile, string &TemperatureDataType) {
    
    ifstream InputData;
    string Colon = ":";
    InputData.open("MasterInputs.txt");
    bool SkippingLines = true;
    while(SkippingLines) {
        string dummyline;
        getline(InputData,dummyline);
        if (dummyline == "*****")
            SkippingLines = false;
    }
    string ValueRead;
    // Simulation Type
    getline(InputData,ValueRead);
    std::size_t found = ValueRead.find(Colon);
    std::string str2 = ValueRead.substr(found+1,string::npos);
    TemperatureDataType = str2;
    // Decomposition strategy
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    DecompositionStrategy = stoi(str2,nullptr,10);
    // Interfacial response function A
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    AConst = atof(str2.c_str());
    // Interfacial response function B
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    BConst = atof(str2.c_str());
    // Interfacial response function C
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    CConst = atof(str2.c_str());
    // Interfacial response function D
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    DConst = atof(str2.c_str());
    // CA cell size
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    deltax = atof(str2.c_str())*pow(10,-6);
    // Nucleation density
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    double NRead = atof(str2.c_str());
    NMax = NRead*pow(10,12);
    // Mean nucleation undercooling
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    dTN = atof(str2.c_str());
    // Standard deviation of nucleation undercooling
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    dTsigma = atof(str2.c_str());
    string Quote = "'";
    // Output base file name
    getline(InputData,ValueRead);
    std::size_t found1 = ValueRead.find_first_of(Quote);
    std::size_t found2 = ValueRead.find_last_of(Quote);
    BaseFileName = ValueRead.substr(found1+1,found2-found1-1);
    // File of grain orientations
    getline(InputData,ValueRead);
    found1 = ValueRead.find_first_of(Quote);
    found2 = ValueRead.find_last_of(Quote);
    GrainOrientationFile = ValueRead.substr(found1+1,found2-found1-1);
    InputData.close();

}


void CInputRead(double &G, double &R, int &nx, int &ny, int &nz, double deltax, double &deltat, double &FractSurfaceSitesActive) {
    
    ifstream InputData;
    string Colon = ":";
    InputData.open("Inputs_ConstrainedSolidification.txt");
    bool SkippingLines = true;
    while(SkippingLines) {
        string dummyline;
        getline(InputData,dummyline);
        if (dummyline == "*****")
            SkippingLines = false;
    }
    string ValueRead;
    // Thermal gradient
    getline(InputData,ValueRead);
    std::size_t found = ValueRead.find(Colon);
    std::string str2 = ValueRead.substr(found+1,string::npos);
    G = atof(str2.c_str());
    // Cooling rate
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    R = atof(str2.c_str());
    // Domain size in x
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    nx = stoi(str2,nullptr,10);
    nx = nx+2;
    // Domain size in y
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    ny = stoi(str2,nullptr,10);
    ny = ny+2;
    // Domain size in z
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    nz = stoi(str2,nullptr,10);
    nz = nz+2;
    // delta t (using ratio between cooling rate R, thermal gradient G, and cell size delta x
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    int NRatio = stoi(str2,nullptr,10);
    deltat = deltax/(NRatio*(R/G));
    // Fraction of bottom surface sites active
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    FractSurfaceSitesActive = atof(str2.c_str());
    InputData.close();
}


void RInputRead(string &tempfile, double &HT_deltax, double &deltat, int &NumberOfLayers, int &LayerOffset, string &BaseFileName, string &TemperatureDataSource, string &SubstrateFileName, bool &LayerwiseTemeperature, int &TempFilesInSeries) {
    
    ifstream InputData;
    string Colon = ":";
    string Quote = "'";
    InputData.open("Inputs_ReducedTemperature.txt");
    bool SkippingLines = true;
    while(SkippingLines) {
        string dummyline;
        getline(InputData,dummyline);
        if (dummyline == "*****")
            SkippingLines = false;
    }
    string ValueRead;
    // File containing temperature data
    getline(InputData,ValueRead);
    std::size_t found1 = ValueRead.find_first_of(Quote);
    std::size_t found2 = ValueRead.find_last_of(Quote);
    tempfile = ValueRead.substr(found1+1,found2-found1-1);
    // Heat transport mesh size
    getline(InputData,ValueRead);
    std::size_t found = ValueRead.find(Colon);
    string str2 = ValueRead.substr(found+1,string::npos);
    HT_deltax = atof(str2.c_str())*pow(10,-6);
    // time step (s)
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    deltat = atof(str2.c_str())*pow(10,-6);
    // Number of layers
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    NumberOfLayers = stoi(str2,nullptr,10);
    // Layer height
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    LayerOffset = stoi(str2,nullptr,10);
    // Usage of burst buffer for parameters
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    int BurstBuffer = stoi(str2,nullptr,10);
    if (BurstBuffer != 0) {
        // tempfile, BaseFileName are altered
        char* InPath;
        char* InFileName;
        char* OutPath;
        InPath = getenv ("PATH_TO_INPUT");
        InFileName = getenv ("TRUCHAS_INPUT_FILENAME");
        OutPath = getenv ("PATH_TO_OUTPUT");
        string InPathS(InPath);
        string InFileNameS(InFileName);
        string OutPathS(OutPath);
        OutPathS = OutPathS + "/";
        tempfile.clear();
        tempfile = InPathS + "/";
        tempfile += InFileNameS;
        BaseFileName = OutPathS + BaseFileName;
    }
    // OpenFOAM or Truchas as temperature data source?
    getline(InputData,ValueRead);
    found = ValueRead.find(Quote);
    TemperatureDataSource = ValueRead.substr(found+1,1);
    //cout << "T Data Source = " << TemperatureDataSource << endl;
    // Name of substrate file
    getline(InputData,ValueRead);
    found1 = ValueRead.find_first_of(Quote);
    found2 = ValueRead.find_last_of(Quote);
    SubstrateFileName = ValueRead.substr(found1+1,found2-found1-1);
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    if (str2 == "Y") LayerwiseTemeperature = true;
    else LayerwiseTemeperature = false;
    // Temperature files in series (used if LayerwiseTemperature = true)
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    TempFilesInSeries = stoi(str2,nullptr,10);
    InputData.close();
}

void ParallelMeshInit(double &G, double &R, int &DecompositionStrategy, int (&NeighborX)[26], int (&NeighborY)[26], int (&NeighborZ)[26], int (&ItList)[9][26], string TemperatureDataType, int ierr, int id, int np, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset,int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn, int &MyLeftOut, int &MyRightIn, int &MyRightOut, double &deltax, double HT_deltax, double &deltat, int &nx, int &ny, int &nz, int &ProcessorsInXDirection, int &ProcessorsInYDirection, string tempfile, float &XMin, float &XMax, float &YMin, float &YMax, float &ZMin, float &ZMax, string TemperatureDataSource, bool LayerwiseTemeperature) {
        
    // Assignment of neighbors around a cell "X" is as follows (in order of closest to furthest from cell "X")
    // Neighbors 0 through 8 are in the -Y direction
    // Neighbors 9 through 16 are in the XY plane with cell X
    // Neighbors 17 through 25 are in the +Y direction

    NeighborX[0] = 0; NeighborY[0] = -1; NeighborZ[0] = 0;
    NeighborX[1] = 1; NeighborY[1] = -1; NeighborZ[1] = 0;
    NeighborX[2] = -1; NeighborY[2] = -1; NeighborZ[2] = 0;
    NeighborX[3] = 0; NeighborY[3] = -1; NeighborZ[3] = 1;
    NeighborX[4] = 0; NeighborY[4] = -1; NeighborZ[4] = -1;
    NeighborX[5] = -1; NeighborY[5] = -1; NeighborZ[5] = 1;
    NeighborX[6] = 1; NeighborY[6] = -1; NeighborZ[6] = 1;
    NeighborX[7] = -1; NeighborY[7] = -1; NeighborZ[7] = -1;
    NeighborX[8] = 1; NeighborY[8] = -1; NeighborZ[8] = -1;
    
    NeighborX[9] = 0; NeighborY[9] = 0; NeighborZ[9] = 1;
    NeighborX[10] = 0; NeighborY[10] = 0; NeighborZ[10] = -1;
    NeighborX[11] = 1; NeighborY[11] = 0; NeighborZ[11] = 1;
    NeighborX[12] = -1; NeighborY[12] = 0; NeighborZ[12] = 1;
    NeighborX[13] = 1; NeighborY[13] = 0; NeighborZ[13] = -1;
    NeighborX[14] = -1; NeighborY[14] = 0; NeighborZ[14] = -1;
    NeighborX[15] = 1; NeighborY[15] = 0; NeighborZ[15] = 0;
    NeighborX[16] = -1; NeighborY[16] = 0; NeighborZ[16] = 0;
    
    NeighborX[17] = 0; NeighborY[17] = 1; NeighborZ[17] = 0;
    NeighborX[18] = 1; NeighborY[18] = 1; NeighborZ[18] = 0;
    NeighborX[19] = -1; NeighborY[19] = 1; NeighborZ[19] = 0;
    NeighborX[20] = 0; NeighborY[20] = 1; NeighborZ[20] = 1;
    NeighborX[21] = 0; NeighborY[21] = 1; NeighborZ[21] = -1;
    NeighborX[22] = 1; NeighborY[22] = 1; NeighborZ[22] = 1;
    NeighborX[23] = -1; NeighborY[23] = 1; NeighborZ[23] = 1;
    NeighborX[24] = 1; NeighborY[24] = 1; NeighborZ[24] = -1;
    NeighborX[25] = -1; NeighborY[25] = 1; NeighborZ[25] = -1;
    
//    // Opposite neighbors
//    Opp_Neighbor[0] = 17;
//    Opp_Neighbor[1] = 19;
//    Opp_Neighbor[2] = 18;
//    Opp_Neighbor[3] = 21;
//    Opp_Neighbor[4] = 20;
//    Opp_Neighbor[5] = 24;
//    Opp_Neighbor[6] = 25;
//    Opp_Neighbor[7] = 22;
//    Opp_Neighbor[8] = 23;
//
//    Opp_Neighbor[9] = 10;
//    Opp_Neighbor[10] = 9;
//    Opp_Neighbor[11] = 14;
//    Opp_Neighbor[12] = 13;
//    Opp_Neighbor[13] = 12;
//    Opp_Neighbor[14] = 11;
//    Opp_Neighbor[15] = 16;
//    Opp_Neighbor[16] = 15;
//
//    Opp_Neighbor[17] = 0;
//    Opp_Neighbor[18] = 2;
//    Opp_Neighbor[19] = 1;
//    Opp_Neighbor[20] = 4;
//    Opp_Neighbor[21] = 3;
//    Opp_Neighbor[22] = 7;
//    Opp_Neighbor[23] = 8;
//    Opp_Neighbor[24] = 5;
//    Opp_Neighbor[25] = 6;

    
    // If X and Y coordinates are not on edges, Case 0: iteratation over neighbors 0-25 possible
    for (int i=0; i<=25; i++) {
        ItList[0][i] = i;
    }
    // If Y coordinate is on lower edge, Case 1: iteration over only neighbors 9-25 possible
    int Pos = 0;
    for (int i=0; i<=25; i++) {
        if (NeighborY[i] != -1) {
            ItList[1][Pos] = i;
            Pos++;
        }
    }
    // If Y coordinate is on upper edge, Case 2: iteration over only neighbors 0-16 possible
    Pos = 0;
    for (int i=0; i<=25; i++) {
        if (NeighborY[i] != 1) {
            ItList[2][Pos] = i;
            Pos++;
        }
    }
    // If X coordinate is on lower edge, Case 3: iteration over only neighbors 0,1,3,4,6,8,9,10,11,13,15,17,18,20,21,22,24
    Pos = 0;
    for (int i=0; i<=25; i++) {
        if (NeighborX[i] != -1) {
            ItList[3][Pos] = i;
            Pos++;
        }
    }
    // If X coordinate is on upper edge, Case 4: iteration over only neighbors 0,2,3,4,5,7,9,10,12,14,16,17,19,20,21,23,25
    Pos = 0;
    for (int i=0; i<=25; i++) {
        if (NeighborX[i] != 1) {
            ItList[4][Pos] = i;
            Pos++;
        }
    }
    // If X/Y coordinates are on lower edge, Case 5: iteration over only neighbors 9,10,11,13,15,17,18,20,21,22,24
    Pos = 0;
    for (int i=0; i<=25; i++) {
        if ((NeighborX[i] != -1)&&(NeighborY[i] != -1)) {
            ItList[5][Pos] = i;
            Pos++;
        }
    }
    // If X coordinate is on upper edge/Y on lower edge, Case 6:
    Pos = 0;
    for (int i=0; i<=25; i++) {
        if ((NeighborX[i] != 1)&&(NeighborY[i] != -1)) {
            ItList[6][Pos] = i;
            Pos++;
        }
    }
    // If X coordinate is on lower edge/Y on upper edge, Case 7:
    Pos = 0;
    for (int i=0; i<=25; i++) {
        if ((NeighborX[i] != -1)&&(NeighborY[i] != 1)) {
            ItList[7][Pos] = i;
            Pos++;
        }
    }
    // If X/Y coordinates are on upper edge, Case 8:
    Pos = 0;
    for (int i=0; i<=25; i++) {
        if ((NeighborX[i] != 1)&&(NeighborY[i] != 1)) {
            ItList[8][Pos] = i;
            Pos++;
        }
    }

    // Mesh initialization for each solidification problem
    
    if (TemperatureDataType == "R") {
        // Determine mesh size needed based on OpenFOAM/Truchas/Analytical model data
        // Read geometry data from OpenFOAM/Truchas
        int nx_HT, ny_HT, nz_HT; // OpenFOAM/Truchas mesh limits

        XMin = 1000000.0;
        YMin = 1000000.0;
        ZMin = 1000000.0;
        XMax = -1000000.0;
        YMax = -1000000.0;
        ZMax = -1000000.0;
        ifstream Geom;
        // Find bounds of temperature data - if LayerwiseTemperature is true, all files must be read
        bool FindingBounds = true;
        int LayerReadCount = 1;
        string tempfile_thislayer;
        if (LayerwiseTemeperature) tempfile_thislayer = "1" + tempfile;
        else tempfile_thislayer = tempfile;
        Geom.open(tempfile_thislayer);
        if (id == 0) cout << "Layer file 1 is " << tempfile_thislayer << endl;
        while (FindingBounds) {
            
            while (!Geom.eof()) {
                string s;
                getline(Geom,s);
                bool ReadingLine = true;
                int i = 0;
                int j = 0;
                int FirstChar = 0;
                int LastChar;
                
                while (ReadingLine) {
                    
                    if (s.empty()) break;
                    char C = s.at(i);
                    // If this character is not a space, convert from string
                    if (isblank(C)) {
                        
                        LastChar = i;
                        string NewDataS = s.substr(FirstChar,LastChar-FirstChar);
                        float MeshData = atof(NewDataS.c_str());
                        bool Searching = true;
                        while (Searching) {
                            i++;
                            C = s.at(i);
                            if (!(isblank(C))) {
                                FirstChar = i;
                                Searching = false;
                            }
                        }
                        
                        if (j == 0) {
                            if (XMin > MeshData) XMin = MeshData;
                            if (XMax < MeshData) XMax = MeshData;
                            //cout << XMax << endl;
                        }
                        else if (j == 1) {
                            //if (id == 0) cout << "Y " << MeshData << endl;
                            if (YMin > MeshData) YMin = MeshData;
                            if (YMax < MeshData) YMax = MeshData;
                        }
                        else if (j == 2) {
                            //if (id == 0) cout << "Z " << MeshData << endl;
                            if (ZMin > MeshData) ZMin = MeshData;
                            if (ZMax < MeshData) ZMax = MeshData;
                            ReadingLine = false;
                        }
                        j++;
                    }
                    // advance to the next character
                    i++;
                }
            }
            Geom.close();
            if (!LayerwiseTemeperature) FindingBounds = false;
            else {
                LayerReadCount++;
                string NextLayerFileS = to_string(LayerReadCount);
                tempfile_thislayer = NextLayerFileS + tempfile;
                Geom.open(tempfile_thislayer);
                if (!Geom.good()) {
                    FindingBounds = false;
                    Geom.close();
                }
            }
            if (id == 0) cout << "Layer file " << LayerReadCount << " is " << tempfile_thislayer << endl;
        }

        // Is the input in m (OpenFOAM) or mm (Truchas)?
        double UnitConversion;
        if (TemperatureDataSource == "O") UnitConversion = 1;
        else UnitConversion = 1000;
        XMin = XMin/UnitConversion;
        XMax = XMax/UnitConversion;
        YMin = YMin/UnitConversion;
        YMax = YMax/UnitConversion;
        ZMin = ZMin/UnitConversion;
        ZMax = ZMax/UnitConversion;

        // OpenFOAM/Truchas nodes containing temperature data
        nx_HT = round((XMax-XMin)/HT_deltax) + 1;
        ny_HT = round((YMax-YMin)/HT_deltax) + 1;
        nz_HT = round((ZMax-ZMin)/HT_deltax) + 1;
        
        // Ratio of temperature-CA grid sizes
        int HTratio = round(HT_deltax/deltax);
        
        // CA nodes in each direction (+2 for wall cells at the boundaries) (+2 for solid cells at X/Y boundaries, +1 for solid cells at lower Z boundary)
        nx = (nx_HT - 1)*HTratio + 1 + 2 + 2;
        ny = (ny_HT - 1)*HTratio + 1 + 2 + 2;
        nz = (nz_HT - 1)*HTratio + 1 + 2 + 1;
    
        if (id == 0) {
            if (TemperatureDataSource == "O") cout << "CA model of OpenFOAM heat transport problem solidification" << endl;
            else cout << "CA model of Truchas heat transport problem solidification" << endl;
            cout << "Cell size (m): " << deltax << endl;
            cout << "Time step (s): " << deltat << endl;
            cout << "Domain size: " << nx << " by " << ny << " by " << nz << endl;
            cout << "X Limits of domain: " << XMin << " and " << XMax << endl;
            cout << "Y Limits of domain: " << YMin << " and " << YMax << endl;
            cout << "Z Limits of domain: " << ZMin << " and " << ZMax << endl;
            cout << "================================================================" << endl;
        }
    }
    else {
        if (id == 0) {
            cout << "================================================================" << endl;
            cout << "CA model of fixed temperature field solidification" << endl;
            cout << "Cell size (m): " << deltax << endl;
            cout << "Time step (s): " << deltat << endl;
            cout << "Domain size: " << nx << " by " << ny << " by " << nz << endl;
            cout << "Thermal gradient: " << G << " K/m, Cooling rate: " << R << " K/s" << endl;
            cout << "================================================================" << endl;
        }
        
        
    }

    InitialDecomposition(DecompositionStrategy, nx, ny, ProcessorsInXDirection, ProcessorsInYDirection, id, np, MyXSlices, MyYSlices, MyXOffset, MyYOffset, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyLeftOut, MyRightIn, MyRightOut);

    MyXOffset = XOffsetCalc(id,nx,ProcessorsInXDirection,ProcessorsInYDirection,np,DecompositionStrategy);
    MyXSlices = XMPSlicesCalc(id,nx,ProcessorsInXDirection,ProcessorsInYDirection,np,DecompositionStrategy);

    MyYOffset = YOffsetCalc(id,ny,ProcessorsInYDirection,np,DecompositionStrategy);
    MyYSlices = YMPSlicesCalc(id,ny,ProcessorsInYDirection,np,DecompositionStrategy);
   // cout << "ID = " << id << " X RANGE = " << MyXOffset << " TO = " << MyXOffset+MyXSlices-1 << " ; YRANGE = " << MyYOffset << " TO = " << MyYOffset+MyYSlices-1 << endl;
}

void TempInit(int layernumber, int TempFilesInSeries, double G, double R, int DecompositionStrategy, int (&NeighborX)[26], int (&NeighborY)[26], int (&NeighborZ)[26], int (&ItList)[9][26], string TemperatureDataType, int ierr, int id, int np, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset,int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn, int &MyLeftOut, int &MyRightIn, int &MyRightOut, double deltax, double HT_deltax, double deltat, int &nx, int &ny, int &nz, int &ProcessorsInXDirection, int &ProcessorsInYDirection, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange, ViewF::HostMirror UndercoolingCurrent, string tempfile, float XMin, float XMax, float YMin, float YMax, float ZMin, float ZMax, bool* Melted, string TemperatureDataSource, bool LayerwiseTemeperature) {

    if (TemperatureDataType == "C") {
        
        // Contrained solidification test problem
        // Initialize temperature field in Z direction with thermal gradient G set in input file
        for (int k=0; k<nz; k++) {
            for (int i=0; i<MyXSlices; i++) {
                for (int j=0; j<MyYSlices; j++) {
                    UndercoolingCurrent(k*MyXSlices*MyYSlices + i*MyYSlices + j) = 0;
                    UndercoolingChange(k*MyXSlices*MyYSlices + i*MyYSlices + j) = R*deltat;
                    CritTimeStep(k*MyXSlices*MyYSlices + i*MyYSlices + j) = (int)(((k-1)*G*deltax)/(R*deltat));
                    int GlobalX = i + MyXOffset;
                    int GlobalY = j + MyYOffset;
                    if ((GlobalX > -1)&&(GlobalX < nx)&&(GlobalY > -1)&&(GlobalY < ny)&&(k > 0)&&(k < nz-1)) {
                        Melted[k*MyXSlices*MyYSlices+i*MyYSlices+j] = true;
                    }
                    else {
                        Melted[k*MyXSlices*MyYSlices+i*MyYSlices+j] = false;
                    }
                }
            }
        }
   }
   else {
            
        // Temperature data read
        int HTtoCAratio = HT_deltax/deltax; // OpenFOAM/CA cell size ratio

       // With row/col 0 being wall cells and row/col 1 being solid cells outside of the melted area, the domain starts at row/col 2
       // As the wall cells are not part of the physical domain (solid cells at row/col 1 are defined as X = Y = 0,
       // the melted region domain data starts at X = Y = deltax, with data points at X or Y = deltax + N*HT_deltax through X or Y = nx-3 or ny-3
       
       // The X and Y bounds are the region (for this MPI rank) of the physical domain that needs to be read
       // Extends past the actual spatial extent of the local domain for purposes of interpolating from HT_deltax to deltax
       int LowerXBound, LowerYBound, UpperXBound, UpperYBound;
       if (MyXOffset <= 2) LowerXBound = 2;
       else LowerXBound = MyXOffset - ((MyXOffset-2) % HTtoCAratio);
       if (MyYOffset <= 2) LowerYBound = 2;
       else LowerYBound = MyYOffset - ((MyYOffset-2) % HTtoCAratio);
       
       if (MyXOffset+MyXSlices-1 >= nx-3) UpperXBound = nx-3;
       else UpperXBound = MyXOffset + MyXSlices - 1 + HTtoCAratio - ((MyXOffset + (MyXSlices-1) - 2) % HTtoCAratio);
       if (MyYOffset+MyYSlices-1 >= ny-3) UpperYBound = ny-3;
       else UpperYBound = MyYOffset + MyYSlices - 1 + HTtoCAratio - ((MyYOffset + (MyYSlices-1) - 2) % HTtoCAratio);

       vector <vector <vector <double> > > CritTS, CritTL;
       for (int k=0; k<nz; k++) {
           vector <vector <double> > TemperatureXX;
           for (int i=LowerXBound; i<=UpperXBound; i++) {
               vector <double> TemperatureX;
               for (int j=LowerYBound; j<=UpperYBound; j++) {
                   TemperatureX.push_back(-1.0);
               }
               TemperatureXX.push_back(TemperatureX);
           }
           CritTS.push_back(TemperatureXX);
           CritTL.push_back(TemperatureXX);
       }

       double UnitConversion;
       if (TemperatureDataSource == "O") UnitConversion = 1;
       else UnitConversion = 1000;
       
       
       ifstream Geom;
       string thislayertempfile;
       
       if (LayerwiseTemeperature) {
           if (layernumber == -1) {
               // No layers printed yet
               thislayertempfile = "1" + tempfile;
               // No sites have melted yet
               for (int i=0; i<MyXSlices*MyYSlices*nz; i++) {
                   Melted[i] = false;
               }
           }
           else {
               // "layernumber" is the layer that was just printed
               int NextLayerFile = layernumber + 2;
               string NextLayerFileS;
               NextLayerFile = NextLayerFile % TempFilesInSeries;
               if (NextLayerFile == 0) NextLayerFile = TempFilesInSeries;
               NextLayerFileS = to_string(NextLayerFile);
               thislayertempfile = NextLayerFileS + tempfile;
               ifstream TemperatureFileRead;
               TemperatureFileRead.open(thislayertempfile);
               if (!TemperatureFileRead.good()) {
                   if (id == 0) cout << "Temperature file " << thislayertempfile << " not found" << endl;
                   TemperatureFileRead.close();
               }
           }
       }
       else thislayertempfile = tempfile;
       
       if (id == 0) cout << "Layer " << layernumber+1 << " temperature file is: " << thislayertempfile << endl;
       //cout << "ID = " << id << " X Bounds are " << LowerXBound << "/" << UpperXBound << " , Y Bounds are " << LowerYBound << "/" << UpperYBound << endl;
       Geom.open(thislayertempfile);
       //if (id == 0) {
       while (!Geom.eof()) {
           string s;
           getline(Geom,s);
           bool ReadingLine = true;
           int i = 0;
           int j = 0;
           int XInt, YInt, ZInt;
           double MyX, MyY, MyZ, FirstTime, LastTime;
           
           int FirstChar = 0;
           int LastChar;
           while (ReadingLine) {
               if (s.empty()) break;
               char C = s.at(i);
               // If this character is not a space, convert from string
               if ((isblank(C))||(i == s.length()-1)) {
                   string NewDataS;
                   if (j < 4) {
                       LastChar = i;
                       NewDataS = s.substr(FirstChar,LastChar-FirstChar);
                       if (j < 3) {
                           bool Searching = true;
                           while (Searching) {
                               //cout << "N searching: Char " << i << " is " << s.at(i) << " of line " << s << endl;
                               i++;
                               C = s.at(i);
                               if (!(isblank(C))) {
                                   FirstChar = i;
                                   Searching = false;
                               }
                           }
                       }
                   }
                   else {
                       i = s.length()-1;
                       bool Searching = true;
                       while (Searching) {
                           //cout << "Space searching: Char " << i << " is " << s.at(i) << " of line " << s << endl;
                           C = s.at(i);
                           if (isblank(C)) {
                               FirstChar = i;
                               Searching = false;
                           }
                           i--;
                       }
                       //cout << "Chars " << FirstChar << " + " << s.length()-FirstChar << endl;
                       NewDataS = s.substr(FirstChar,s.length()-FirstChar);
                       //cout << NewDataS << endl;
                   }
                   if (j == 0) {
                       MyX = atof(NewDataS.c_str())/UnitConversion; // Only divide by 1000 if from Truchas
                       XInt = round((MyX-XMin)/deltax) + 2; // + HTtoCAratio;
                       //if (id == 0) cout << XInt << endl;
                   }
                   else if (j == 1) {
                       //if (id == 0) cout << "Y " << MeshData << endl;
                       MyY = atof(NewDataS.c_str())/UnitConversion; // Only divide by 1000 if from Truchas
                       YInt = round((MyY-YMin)/deltax) + 2;// + HTtoCAratio;
                       //if (id == 0) cout << YInt << endl;
                   }
                   else if (j == 2) {
                       //if (id == 0) cout << "Z " << MeshData << endl;
                       MyZ = atof(NewDataS.c_str())/UnitConversion; // Only divide by 1000 if from Truchas
                       ZInt = round((MyZ-ZMin)/deltax) + 2;
                       //if (id == 0) cout << ZInt << " " << nz << endl;
                   }
                   else if (j == 3) {
                       // Determine if this X and Y fall in this rank's range of interest
                       FirstTime = stod(NewDataS.c_str());
                       //cout << XInt-LowerXBound << " " << YInt-LowerYBound << " " << ZInt << " " << FirstTime/UnitConversion << " " << LastTime/UnitConversion << " of " << UpperXBound-LowerXBound+1 << " " << UpperYBound-LowerYBound+1 << " " << nz-1 << endl;
                       if ((XInt >= LowerXBound)&&(XInt <= UpperXBound)&&(YInt >= LowerYBound)&&(YInt <= UpperYBound)) {
                           CritTL[ZInt][XInt-LowerXBound][YInt-LowerYBound] = FirstTime/UnitConversion; // Truchas in ms
                       }

                   }
                   else if (j == 4) {
                       LastTime = stod(NewDataS.c_str());
                       if ((XInt >= LowerXBound)&&(XInt <= UpperXBound)&&(YInt >= LowerYBound)&&(YInt <= UpperYBound)) {
                           CritTS[ZInt][XInt-LowerXBound][YInt-LowerYBound] = LastTime/UnitConversion; // Truchas in ms
                       }
                       ReadingLine = false;
                       //if (FirstTime > LastTime) cout << "FT = " << FirstTime << " LT = " << LastTime << endl;
                   }
                   j++;
               }
               // advance to the next character
               i++;
           }
       }
      // }
       Geom.close();
       
       MPI_Barrier(MPI_COMM_WORLD);
       if (id == 0) cout << "Data read " << endl;
       // Data interpolation between heat transport and CA grids, if necessary
       if (HTtoCAratio != 1) {
           for (int k=2; k<=nz-2; k++) {
               
               int LowZ = k - ((k-2) % HTtoCAratio);
               int HighZ = LowZ + HTtoCAratio;
               double FHighZ = (double)(k - LowZ)/(double)(HTtoCAratio);
               double FLowZ = 1.0 - FHighZ;
               if (HighZ >= nz-2) HighZ = LowZ;
               
               for (int i=0; i<=UpperXBound-LowerXBound; i++) {
                   int LowX =  i - (i % HTtoCAratio);
                   int HighX = LowX + HTtoCAratio;
                   double FHighX = (double)(i - LowX)/(double)(HTtoCAratio);
                   double FLowX = 1.0 - FHighX;
                   if (HighX >= UpperXBound-LowerXBound) HighX = UpperXBound-LowerXBound;
                   
                   for (int j=0; j<=UpperYBound-LowerYBound; j++) {
                       int LowY = j - (j % HTtoCAratio);
                       int HighY = LowY + HTtoCAratio;
                       double FHighY = (float)(j - LowY)/(float)(HTtoCAratio);
                       double FLowY = 1.0 - FHighY;
                       if (HighY >= UpperYBound-LowerYBound) HighY = UpperYBound-LowerYBound;
                       //if (id == 0) cout << "i = " << i << " j = " << j << " k = " << k << " Bounds: " << LowX << " " << HighX << " " << LowY << " " << HighY << " " << LowZ << " " << HighZ << endl;
                       double Pt1 = CritTL[LowZ][LowX][LowY];
                       double Pt2 = CritTL[LowZ][HighX][LowY];
                       double Pt12 = FLowX*Pt1 + FHighX*Pt2;
                       double Pt3 = CritTL[LowZ][LowX][HighY];
                       double Pt4 = CritTL[LowZ][HighX][HighY];
                       double Pt34 = FLowX*Pt3 + FHighX*Pt4;
                       double Pt1234 = Pt12*FLowY + Pt34*FHighY;
                       double Pt5 = CritTL[HighZ][LowX][LowY];
                       double Pt6 = CritTL[HighZ][HighX][LowY];
                       double Pt56 = FLowX*Pt5 + FHighX*Pt6;
                       double Pt7 = CritTL[HighZ][LowX][HighY];
                       double Pt8 = CritTL[HighZ][HighX][HighY];
                       double Pt78 = FLowX*Pt7 + FHighX*Pt8;
                       double Pt5678 = Pt56*FLowY + Pt78*FHighY;
                       if ((Pt1 > 0)&&(Pt2 > 0)&&(Pt3 > 0)&&(Pt4 > 0)&&(Pt5 > 0)&&(Pt6 > 0)&&(Pt7 > 0)&&(Pt8 > 0)) {
                           CritTL[k][i][j] = Pt1234*FLowZ + Pt5678*FHighZ;
                           // cout << CritTL[k][i][j] << endl;
                       }
                       Pt1 = CritTS[LowZ][LowX][LowY];
                       Pt2 = CritTS[LowZ][HighX][LowY];
                       Pt12 = FLowX*Pt1 + FHighX*Pt2;
                       Pt3 = CritTS[LowZ][LowX][HighY];
                       Pt4 = CritTS[LowZ][HighX][HighY];
                       Pt34 = FLowX*Pt3 + FHighX*Pt4;
                       Pt1234 = Pt12*FLowY + Pt34*FHighY;
                       Pt5 = CritTS[HighZ][LowX][LowY];
                       Pt6 = CritTS[HighZ][HighX][LowY];
                       Pt56 = FLowX*Pt5 + FHighX*Pt6;
                       Pt7 = CritTS[HighZ][LowX][HighY];
                       Pt8 = CritTS[HighZ][HighX][HighY];
                       Pt78 = FLowX*Pt7 + FHighX*Pt8;
                       Pt5678 = Pt56*FLowY + Pt78*FHighY;
                       if ((Pt1 > 0)&&(Pt2 > 0)&&(Pt3 > 0)&&(Pt4 > 0)&&(Pt5 > 0)&&(Pt6 > 0)&&(Pt7 > 0)&&(Pt8 > 0)) {
                           CritTS[k][i][j] = Pt1234*FLowZ + Pt5678*FHighZ;
                       }
                       //if (k == 56) cout << LowZ << " " << HighZ << " " << LowX << " " << HighX << " " << LowY << " " << HighY << endl;
                   }
               }
           }
       }
       // }
       MPI_Barrier(MPI_COMM_WORLD);
       if (id == 0) cout << "Interpolation done" << endl;
       
       //  cout << "ID = " << id << " XStart = " << MyXOffset << " XEnd = " << MyXOffset+MyXSlices-1 << " YStart = " << MyYOffset << " YEnd = " << MyYOffset+MyYSlices-1 << endl;
       
       // Convert CritTL, CritTS matrices into CritTimeStep and UndercoolingChange (change in undercooling with time step)
       int MCTS = 0;
       int MinCTS = 10000000;
       for (int k=0; k<nz; k++) {
           for (int ii=LowerXBound; ii<=UpperXBound; ii++) {
               for (int jj=LowerYBound; jj<=UpperYBound; jj++) {
                   if ((ii >= MyXOffset)&&(ii < MyXOffset+MyXSlices)&&(jj >= MyYOffset)&&(jj < MyYOffset+MyYSlices)) {
                       int Adj_i = ii - MyXOffset;
                       int Adj_j = jj - MyYOffset;
                       double CTLiq = CritTL[k][ii-LowerXBound][jj-LowerYBound];
                       int Coord3D1D = k*MyXSlices*MyYSlices + Adj_i*MyYSlices + Adj_j;
                       if (CTLiq > 0)  {
                           Melted[Coord3D1D] = true;
                           //if (id == 0) cout << "True" << endl;
                           CritTimeStep(Coord3D1D) = round(CTLiq/deltat);
                           //cout << 200*k << endl;
                           UndercoolingChange(Coord3D1D) = (1610.0-1420.0)*deltat/(CritTS[k][ii-LowerXBound][jj-LowerYBound] - CTLiq);
                           // cout << CritTimeStep[Coord3D1D] << " " << UndercoolingChange[Coord3D1D] << endl;
                           if (CritTimeStep(Coord3D1D) > MCTS) MCTS = CritTimeStep(Coord3D1D);
                           if (CritTimeStep(Coord3D1D) < MinCTS) MinCTS = CritTimeStep(Coord3D1D);
                           //if (id == 1) cout << CritTimeStep[Coord3D1D] << " " << UndercoolingChange[Coord3D1D] << endl;
                           if (UndercoolingChange(Coord3D1D) > 0.25) {
                               UndercoolingChange(Coord3D1D) = 0.25;
                           }
                           UndercoolingCurrent(Coord3D1D) = 0;
                       }
                       else {
                           CritTimeStep(Coord3D1D) = 0;
                           UndercoolingChange(Coord3D1D) = 0.0;
                           UndercoolingCurrent(Coord3D1D) = 0.0;
                       }
                   }
               }
           }
       }
       
        int GMax, GMin;
       //cout << "My min ts = " << MinCTS << endl;
        MPI_Reduce(&MCTS, &GMax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&MinCTS, &GMin, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Bcast(&GMin,1,MPI_INT,0,MPI_COMM_WORLD);
        //cout << "Global min ts = " << GMin << endl;
        for (int k=0; k<nz; k++) {
           for (int i=0; i<MyXSlices; i++) {
              for (int j=0; j<MyYSlices; j++) {
                  int Coord3D1D = k*MyXSlices*MyYSlices + i*MyYSlices + j;
                  if (CritTimeStep(Coord3D1D) != 0) {
                      CritTimeStep(Coord3D1D) = CritTimeStep(Coord3D1D) - GMin;
                      //if (CritTimeStep(Coord3D1D) < 0) cout << "id = " << id << " CTS = " << CritTimeStep(Coord3D1D) << endl;
                  }
              }
           }
        }
        if (id == 0) cout << "Min Crit time step is = " << GMin << endl;
        if (id == 0) cout << "Max Crit time step is = " << GMax << endl;
    }


}
////*****************************************************************************/

// Initialize grain orientations and unit vectors

void OrientationInit(int id, int NGrainOrientations, int* GrainOrientation, float* GrainUnitVector, string GrainOrientationFile) {
    
    // Read file of grain orientations
    ifstream O;
    O.open(GrainOrientationFile);
    
    // Line 1 is the number of orientation values to read (if not specified already)
    string ValueRead;
    getline(O,ValueRead);
    
    // Populate data structure for grain unit vectors
    for(int i=0; i<NGrainOrientations; i++) {
        string s;
        if (!getline(O,s)) break;
        istringstream ss(s);
        int Comp = 0;
        int UVNumber = 0;
        while (ss) { // This is the 3 grain orientation angles
            string s;
            if (!getline(ss,s,',')) break;
            float ReadGO = atof(s.c_str());
            // X,Y,Z of a single unit vector
            GrainUnitVector[18*i + 3*UVNumber + Comp] = ReadGO;
            GrainUnitVector[18*i + 3*(UVNumber+1) + Comp] = -ReadGO;
            Comp++;
            if (Comp > 2) {
                Comp = 0;
                UVNumber = UVNumber + 2;
            }
        }
    }
    O.close();
    
    // The grain orientations that correspond to each Grain ID must be the same across all ranks
    // Shuffle list of "NGrainOrientation" orientations
    int* GrainOrientation_master = new int[NGrainOrientations];
    if (id == 0) {
        for (int h=0; h<NGrainOrientations; h++) {
            //double R = ((double) rand() / (RAND_MAX));
            GrainOrientation_master[h] = h; //round((NGrainOrientations-1)*R);
        }
    }
    MPI_Bcast(&GrainOrientation_master[0], NGrainOrientations, MPI_INT, 0, MPI_COMM_WORLD);
    for (int h=0; h<NGrainOrientations; h++) {
        GrainOrientation[h] = GrainOrientation_master[h];
    }
    
}


////*****************************************************************************/
////*****************************************************************************/
//
///*
// Initializes cell types where the substrate comes from a file
//*/

void GrainInit(int layernumber, int LayerHeight, string TemperatureDataType, string SubstrateFileName, double FractSurfaceSitesActive, int NGrainOrientations, int DecompositionStrategy, int ProcessorsInXDirection, int ProcessorsInYDirection, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int id, int np, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, int ItList[9][26], int NeighborX[26], int NeighborY[26], int NeighborZ[26], int* GrainOrientation, float* GrainUnitVector, ViewF::HostMirror DiagonalLength, ViewI::HostMirror CellType, ViewI::HostMirror TriangleIndex, ViewI::HostMirror GrainID, ViewF::HostMirror CritDiagonalLength, ViewF::HostMirror DOCenter, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange, bool* Melted, double deltax, double NMax, int &NextLayer_FirstNucleatedGrainID, int &PossibleNuclei_ThisRank) {
    
    mt19937_64 gen(id*(layernumber+1));
    uniform_real_distribution<double> dis(0.0, 1.0);
    
    // Convert initial grain spacing to a grain density
    double BulkProb = NMax*deltax*deltax*deltax;
    if (id == 0) cout << "Fraction of melt pool sites to potentially be activated: " << BulkProb << endl;
    
    // Counter for the number of active cells
    int SubstrateActCells_ThisRank = 0;
    PossibleNuclei_ThisRank = 0;
    
    // Wall cells at global domain boundaries
    // Other cells are either regions that will melt, or part of the substrate
    for (int k=0; k<nz; k++)  {
        for(int i=0; i<MyXSlices; i++) {
            for(int j=0; j<MyYSlices; j++) {
                int GlobalX = i + MyXOffset;
                int GlobalY = j + MyYOffset;
                if ((GlobalX == -1)||(GlobalX == nx)||(GlobalY == -1)||(GlobalY == ny)||(k == 0)||(k == nz-1)) {
                    int CAGridLocation = k*MyXSlices*MyYSlices + i*MyYSlices + j;
                    CellType(CAGridLocation) = Wall;
                    GrainID(CAGridLocation) = 0;
                }
            }
        }
    }
    
    if (TemperatureDataType == "C") {
        
        // Constrained solidification test problem - side surfaces are walls, liquid domain
        for (int k=0; k<nz; k++)  {
            for(int i=0; i<MyXSlices; i++) {
                for(int j=0; j<MyYSlices; j++) {
                    int GlobalX = i + MyXOffset;
                    int GlobalY = j + MyYOffset;
                    int CAGridLocation = k*MyXSlices*MyYSlices + i*MyYSlices + j;
                    GrainID(CAGridLocation) = 0;
                    if ((GlobalX != -1)&&(GlobalX != nx)&&(GlobalY != -1)&&(GlobalY != ny)&&(k != 0)&&(k != nz-1)) {
                        Melted[CAGridLocation] = true;
                        CellType(CAGridLocation) = Delayed;
                    }
                }
            }
        }

        // Other cells may be substrate or heterogenous solid sites
        for (int k=1; k<nz-1; k++)  {
            for(int i=1; i<MyXSlices-1; i++) {
                for(int j=1; j<MyYSlices-1; j++) {
                    int CAGridLocation = k*MyXSlices*MyYSlices + i*MyYSlices + j;
                    double R = dis(gen);
                    if (k == 1) {
                        // Randomly locate substrate grain seeds
                        if (R < FractSurfaceSitesActive) {
                            SubstrateActCells_ThisRank++;
                            CellType(CAGridLocation) = Active;
                        }
                    }
                    else {
                        // Randomly locate bulk site seeds
                        if (R < BulkProb) {
                            PossibleNuclei_ThisRank++;
                            CellType(k*MyXSlices*MyYSlices+i*MyYSlices+j) = LiqSol;
                            // GrainID for these are assigned later
                        }
                    }
                }
            }
        }
        
        // Assign grain IDs to bottom surface grains
        int FirstEpitaxialGrainID = 1;
        if (np > 1) {
            // Grains for epitaxial growth - determine GrainIDs on each MPI rank
            if (id == 0) {
                //cout << " ID 0 sending" << endl;
                int SBuf = FirstEpitaxialGrainID+SubstrateActCells_ThisRank;
                MPI_Send(&SBuf,1,MPI_INT,1,0,MPI_COMM_WORLD);
            }
            else if (id == np-1) {
                //cout << " ID " << np-1 << " recieving" << endl;
                int RBuf;
                MPI_Recv(&RBuf,1,MPI_INT,np-2,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                FirstEpitaxialGrainID = RBuf;
            }
            else {
                int RBuf;
                //cout << " ID " << id << " sending/recieving" << endl;
                MPI_Recv(&RBuf,1,MPI_INT,id-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                FirstEpitaxialGrainID = RBuf;
                int SBuf = RBuf + SubstrateActCells_ThisRank;
                MPI_Send(&SBuf,1,MPI_INT,id+1,0,MPI_COMM_WORLD);
            }
        }
        for (int i=0; i<MyXSlices*MyYSlices*nz; i++)  {
            if (CellType(i) == Active)  {
                GrainID(i) = FirstEpitaxialGrainID;
                FirstEpitaxialGrainID++;
            }
        }
        
    }
    else {
        // Assign GrainID values to cells that are part of the substrate
        // Cells that border the melted region are type active, others are type solid
        ifstream Substrate;
        Substrate.open(SubstrateFileName);
        
        int Substrate_LowX = MyXOffset;
        int Substrate_HighX = MyXOffset+MyXSlices;
        int Substrate_LowY = MyYOffset;
        int Substrate_HighY = MyYOffset+MyYSlices;
        int nxS, nyS, nzS;
        string s;
        getline(Substrate,s);
        std::size_t found = s.find("=");
        string str = s.substr(found+1,s.length()-1);
        nzS = stoi(str,nullptr,10);
        getline(Substrate,s);
        found = s.find("=");
        str = s.substr(found+1,s.length()-1);
        nyS = stoi(str,nullptr,10);
        getline(Substrate,s);
        found = s.find("=");
        str = s.substr(found+1,s.length()-1);
        nxS = stoi(str,nullptr,10);
        
        int FileLBound, FileHBound;
        if (layernumber == -1) {
            FileLBound = 0;
            FileHBound = nz-3;
        }
        else {
            FileLBound = nz - 2 + (layernumber)*LayerHeight;
            FileHBound = nz - 2 + (layernumber)*LayerHeight + (LayerHeight-1);
        }

        if (id == 0) {
            if (layernumber == -1) cout << "Initial substrate: reading layers 0 through " << nz-3 << " from file " << SubstrateFileName << endl;
            else cout << "Reading substrate layers " << FileLBound << " through "  << FileHBound << " from file" << endl;
        }
        
        // Assign GrainID values to cells that are part of the substrate
        // Cells that border the melted region are type active, others are type solid
        for (int k=0; k<nzS; k++) {
            if (k < FileLBound) {
                for (int j=0; j<nyS; j++) {
                    for (int i=0; i<nxS; i++) {
                        string GIDVal;
                        getline(Substrate,GIDVal);
                    }
                }
            }
            else if ((k >= FileLBound)&&(k <= FileHBound)) {
                for (int j=0; j<nyS; j++) {
                    for (int i=0; i<nxS; i++) {
                        string GIDVal;
                        getline(Substrate,GIDVal);
                        if ((i >= Substrate_LowX)&&(i < Substrate_HighX)&&(j >= Substrate_LowY)&&(j < Substrate_HighY)) {
                            int CAGridLocation;
                            if (layernumber == -1) CAGridLocation = (k+1)*MyXSlices*MyYSlices + (i-MyXOffset)*MyYSlices + (j-MyYOffset);
                            else CAGridLocation = (nz-2-LayerHeight+1 + (k-FileLBound))*MyXSlices*MyYSlices + (i-MyXOffset)*MyYSlices + (j-MyYOffset);
                            if (CritTimeStep(CAGridLocation) == 0) {
                                GrainID(CAGridLocation) = stoi(GIDVal,nullptr,10);
                                //cout << stoi(GIDVal,nullptr,10) << endl;
                            }
                            else {
                                GrainID(CAGridLocation) = 0;
                            }
                        }
                    }
                }
            }
            else {
                break;
            }
        }
        Substrate.close();
        
//        // Extra set of wall cells around edges for spot melt problem
//        for (int k=0; k<nz; k++)  {
//            for(int i=0; i<MyXSlices; i++) {
//                for(int j=0; j<MyYSlices; j++) {
//                    int GlobalX = i + MyXOffset;
//                    int GlobalY = j + MyYOffset;
//                    if ((GlobalX == 0)||(GlobalX == nx-1)||(GlobalY == 0)||(GlobalY == ny-1)||(GlobalX == 1)||(GlobalX == nx-2)||(GlobalY == 1)||(GlobalY == ny-2)) {
//                        int CAGridLocation = k*MyXSlices*MyYSlices + i*MyYSlices + j;
//                        CellType(CAGridLocation) = Wall;
//                        GrainID(CAGridLocation) = 0;
//                    }
//                }
//            }
//        }
    
        // Count number of active cells are at the solid-liquid boundary, as well as the number of nucleation events that may potentially occur
        for (int k=1; k<nz-1; k++) {
            for (int j=0; j<MyYSlices; j++) {
                for (int i=0; i<MyXSlices; i++) {
                    int CAGridLocation = k*MyXSlices*MyYSlices + i*MyYSlices + j;
                    if (CellType(CAGridLocation) != Wall) {
                        if (CritTimeStep(CAGridLocation) != 0) {
                           // This is a liquid ("Delayed") cell or a nuclei, if not in a ghost node and the RNG places one at this site
                           double R = dis(gen);
                           if (R < BulkProb) {
                               if ((i != 0)&&(i != MyXSlices-1)&&(j != 0)&&(j != MyYSlices-1)) {
                                   PossibleNuclei_ThisRank++;
                                   CellType(CAGridLocation) = LiqSol;
                               }
                           }
                           else {
                               CellType(CAGridLocation) = Delayed;
                               GrainID(CAGridLocation) = layernumber+1;
                           }
                        }
                        else {
                            // This is a solid or active cell, depending on whether it is located at the interface of the liquid
                            // Check to see if this site is actually at the solid-liquid interface
                            int LCount = 0;
                            // Which neighbors should be iterated over?
                            int ItBounds = FindItBounds(i,j,MyXSlices,MyYSlices);
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
                                int MyNeighborX = i + NeighborX[l];
                                int MyNeighborY = j + NeighborY[l];
                                int MyNeighborZ = k + NeighborZ[l];
                                int NeighborD3D1ConvPosition = MyNeighborZ*MyXSlices*MyYSlices + MyNeighborX*MyYSlices + MyNeighborY;
                                if (CritTimeStep(NeighborD3D1ConvPosition) > 0) {
                                    LCount++;
                                }
                            }
                            if (LCount == 0) {
                                // Not at the interface
                                CellType(CAGridLocation) = Solid;
                            }
                            else {
                                // At the interface
                                CellType(CAGridLocation) = Active;
                                UndercoolingChange(CAGridLocation) = 0.1;
                                SubstrateActCells_ThisRank++;
                            }
                        }
                    }
                }
            }
        }
    }
    int TotalSubstrateActCells, TotalNucleatedGrains;
    MPI_Reduce(&SubstrateActCells_ThisRank,&TotalSubstrateActCells,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&PossibleNuclei_ThisRank,&TotalNucleatedGrains,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    if (id == 0) cout << "Number of substrate active cells: " << TotalSubstrateActCells << endl;
    if (id == 0) cout << "Number of potential nucleated grains: " << TotalNucleatedGrains << endl;
    int FirstNucleatedGID_Rank0;
    if (layernumber == -1) {
        FirstNucleatedGID_Rank0 = -1;
    }
    else {
        FirstNucleatedGID_Rank0 = NextLayer_FirstNucleatedGrainID;
    }
    int MyFirstNGrainID; // First GrainID for nuclei on this rank, for this layer
    if (np > 1) {
        // Assign GrainIDs for nucleated grains (negative values)
        // Grains for nucleated growth
        if (id == 0) {
            int SBuf = FirstNucleatedGID_Rank0-PossibleNuclei_ThisRank;
            MPI_Send(&SBuf,1,MPI_INT,1,0,MPI_COMM_WORLD);
            MyFirstNGrainID = FirstNucleatedGID_Rank0;
        }
        else if (id == np-1) {
            int RBuf;
            MPI_Recv(&RBuf,1,MPI_INT,np-2,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            MyFirstNGrainID = RBuf;
            NextLayer_FirstNucleatedGrainID = MyFirstNGrainID - PossibleNuclei_ThisRank;
        }
        else {
            int RBuf;
            MPI_Recv(&RBuf,1,MPI_INT,id-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            MyFirstNGrainID = RBuf;
            int SBuf = MyFirstNGrainID - PossibleNuclei_ThisRank;
            MPI_Send(&SBuf,1,MPI_INT,id+1,0,MPI_COMM_WORLD);
        }
        MPI_Bcast(&NextLayer_FirstNucleatedGrainID, 1, MPI_INT, np-1, MPI_COMM_WORLD);
    }
    else {
        // No communication among ranks
        NextLayer_FirstNucleatedGrainID = FirstNucleatedGID_Rank0-PossibleNuclei_ThisRank;
        MyFirstNGrainID = FirstNucleatedGID_Rank0;
    }

    // Assign Grain IDs to nucleated grains
    // Set up active cell data structures for appropriate cells
    int ANCount = 0;
    int BNCount = 0;
    int CNCount = 0;
    int DNCount = 0;
    int ENCount = 0;
    int FNCount = 0;
    int GNCount = 0;
    int HNCount = 0;
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (id == 0) cout << "A" << endl;
    // Set up active cell octahedra for growth, mark cell data to be communicated across ranks in ghost nodes
    // Assign GrainIDs to nuclei sites
    int NCounter = MyFirstNGrainID;
    
    //if (id == 0) {
    for (int RankZ=1; RankZ<nz-1; RankZ++) {
        for (int RankX=0; RankX<MyXSlices; RankX++) {
            for (int RankY=0; RankY<MyYSlices; RankY++) {
                int D3D1ConvPosition = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                if (CellType(D3D1ConvPosition) == Active) {
                    //cout << "Coords" << RankX << " " << RankY << " " << RankZ << endl;
                    int GlobalX = RankX + MyXOffset;
                    int GlobalY = RankY + MyYOffset;
                    int MyGrainID = GrainID(D3D1ConvPosition);
                    NewGrain(MyXSlices,MyYSlices,nz,RankX,RankY,RankZ,MyGrainID,GlobalX,GlobalY,CellType,GrainID,DiagonalLength,DOCenter);
                    //cout << "Grain of ID " << MyGrainID << " made" << endl;
                    // The orientation for the new grain will depend on its Grain ID
                    int MyOrientation = GrainOrientation[((abs(MyGrainID) - 1) % NGrainOrientations)];
                    // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
                    // (xp,yp,zp) is the new cell's center on the global grid
                    double xp = GlobalX + 0.5;
                    double yp = GlobalY + 0.5;
                    double zp = RankZ + 0.5;
                CritDiagLengthCalc(xp,yp,zp,MyOrientation,RankX,RankY,RankZ,D3D1ConvPosition,DOCenter(3*D3D1ConvPosition),DOCenter(3*D3D1ConvPosition+1),DOCenter(3*D3D1ConvPosition+2),NeighborX,NeighborY,NeighborZ, GrainUnitVector,TriangleIndex,CritDiagonalLength);
                    //cout << "CDL calculated" << endl;
                    if (np > 1) {
                        if (DecompositionStrategy == 1) {
                            if (RankY == 1) {
                                CellType(D3D1ConvPosition) = Ghost1;
                            }
                            else if (RankY == MyYSlices-2) {
                                CellType(D3D1ConvPosition) = Ghost1;
                            }
                        }
                        else {
                            if (RankY == 1) {
                                // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                                if (RankX == MyXSlices-2) {
                                    CellType(D3D1ConvPosition) = Ghost3;
                                }
                                else if (RankX == 1) {
                                    CellType(D3D1ConvPosition) = Ghost3;
                                }
                                else if ((RankX > 1)&&(RankX < MyXSlices-2)) {
                                    // This is being sent to MyLeft
                                    CellType(D3D1ConvPosition) = Ghost1;
                                }
                            }
                            else if (RankY == MyYSlices-2) {
                                // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                                if (RankX == MyXSlices-2) {
                                    CellType(D3D1ConvPosition) = Ghost3;
                                }
                                else if (RankX == 1) {
                                    CellType(D3D1ConvPosition) = Ghost3;
                                }
                                else if ((RankX > 1)&&(RankX < MyXSlices-2)) {
                                    CellType(D3D1ConvPosition) = Ghost1;
                                }
                            }
                            else if ((RankX == 1)&&(RankY > 1)&&(RankY < MyYSlices-2)) {
                                CellType(D3D1ConvPosition) = Ghost1;
                                //if (id == 0) cout << "RANK 0 LISTED " << MyNeighborX << " " << MyNeighborY << " " << MyNeighborZ << endl;
                            }
                            else if ((RankX == MyXSlices-2)&&(RankY > 1)&&(RankY < MyYSlices-2)) {
                                CellType(D3D1ConvPosition) = Ghost1;
                            }
                        }
                    }
                }
                else if (CellType(D3D1ConvPosition) == LiqSol) {
                    // Mark and count the number of nucleation events to be sent to other ranks
                    GrainID(RankZ*MyXSlices*MyYSlices+RankX*MyYSlices+RankY) = NCounter;
                    NCounter--;
                    if (np > 1) {
                        if (DecompositionStrategy == 1) {
                            if (RankY == 1) {
                                ANCount++;
                            }
                            else if (RankY == MyYSlices-2) {
                                BNCount++;
                            }
                        }
                        else {
                            if (RankY == 1) {
                                // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                                if (RankX == MyXSlices-2) {
                                    ENCount++;
                                    CNCount++;
                                    ANCount++;
                                }
                                else if (RankX == 1) {
                                    GNCount++;
                                    DNCount++;
                                    ANCount++;
                                }
                                else if ((RankX > 1)&&(RankX < MyXSlices-2)) {
                                    // This is being sent to MyLeft
                                    ANCount++;
                                }
                            }
                            else if (RankY == MyYSlices-2) {
                                // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                                if (RankX == MyXSlices-2) {
                                    FNCount++;
                                    CNCount++;
                                    BNCount++;
                                }
                                else if (RankX == 1) {
                                    HNCount++;
                                    DNCount++;
                                    BNCount++;
                                }
                                else if ((RankX > 1)&&(RankX < MyXSlices-2)) {
                                    BNCount++;
                                }
                            }
                            else if ((RankX == 1)&&(RankY > 1)&&(RankY < MyYSlices-2)) {
                                DNCount++;
                            }
                            else if ((RankX == MyXSlices-2)&&(RankY > 1)&&(RankY < MyYSlices-2)) {
                                CNCount++;
                            }
                        }
                    }
                }
            }
        }
    }
//    //}
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (id == 0) cout << "B" << endl;
    // Send/recieve number of nuclei in the ghost nodes (ARNCount-HRNCount) so that each rank knows the total number of nuclei in it's domain
    int ARNCount = 0;
    int BRNCount = 0;
    int CRNCount = 0;
    int DRNCount = 0;
    int ERNCount = 0;
    int FRNCount = 0;
    int GRNCount = 0;
    int HRNCount = 0;
    
    // Send BNCount, Recieve ARNCount (send to the right, recieve on the left)
    MPI_Sendrecv(&BNCount,1,MPI_INT,MyRight,0,&ARNCount,1,MPI_INT,MyLeft,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    
    // Send ANCount, Recieve BRNCount (send to the left, recieve on the right)
    MPI_Sendrecv(&ANCount,1,MPI_INT,MyLeft,0,&BRNCount,1,MPI_INT,MyRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    if (DecompositionStrategy != 1) {
        // Send CNCount, Recieve DRNCount (send into the plane, recieve out of the plane)
        MPI_Sendrecv(&CNCount,1,MPI_INT,MyIn,0,&DRNCount,1,MPI_INT,MyOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        
        // Send DNCount, Recieve CRNCount (send out of the plane, recieve into the plane)
        MPI_Sendrecv(&DNCount,1,MPI_INT,MyOut,0,&CRNCount,1,MPI_INT,MyIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        
        
        // Send HNCount, Recieve ERNCount
        MPI_Sendrecv(&HNCount,1,MPI_INT,MyRightOut,0,&ERNCount,1,MPI_INT,MyLeftIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        
        // Send ENCount, Recieve HRNCount
        MPI_Sendrecv(&ENCount,1,MPI_INT,MyLeftIn,0,&HRNCount,1,MPI_INT,MyRightOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        
        // Send GNCount, Recieve FRNCount
        MPI_Sendrecv(&GNCount,1,MPI_INT,MyLeftOut,0,&FRNCount,1,MPI_INT,MyRightIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        
        // Send FNCount, Recieve GRNCount
        MPI_Sendrecv(&FNCount,1,MPI_INT,MyRightIn,0,&GRNCount,1,MPI_INT,MyLeftOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    
    if (MyLeft == MPI_PROC_NULL) ARNCount = 0;
    if (MyRight == MPI_PROC_NULL) BRNCount = 0;
    if (MyIn == MPI_PROC_NULL) CRNCount = 0;
    if (MyOut == MPI_PROC_NULL) DRNCount = 0;
    if (MyLeftIn == MPI_PROC_NULL) ERNCount = 0;
    if (MyRightIn == MPI_PROC_NULL) FRNCount = 0;
    if (MyLeftOut == MPI_PROC_NULL) GRNCount = 0;
    if (MyRightOut == MPI_PROC_NULL) HRNCount = 0;
    //cout << "ID = " << id << " Neighbors are " << MyLeft << " " << MyRight << " " << MyIn << " " << MyOut << endl;
    //cout << "possible nuclei rank " << id << " no ghost nodes: " << PossibleNuclei_ThisRank << endl;
    PossibleNuclei_ThisRank += (ARNCount + BRNCount + CRNCount + DRNCount + ERNCount + FRNCount + GRNCount + HRNCount);
    //cout << "additional nuclei on rank " << id << " : " << ARNCount << " " << BRNCount << " " << CRNCount << " " << DRNCount << " " << ERNCount << " " << FRNCount << " " << GRNCount << " " << HRNCount << " new total " << PossibleNuclei_ThisRank << endl;
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (id == 0) cout << "C" << endl;
    // Remove delay cells not bordering others
    for (int RankZ=1; RankZ<nz-1; RankZ++) {
        for (int RankX=1; RankX<MyXSlices-1; RankX++) {
            for (int RankY=1; RankY<MyYSlices-1; RankY++) {
                int D3D1ConvPosition = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                if (CellType(D3D1ConvPosition) == Delayed) {
                    // Check to see if this cell is at the interface
                    int LCount = 0;
                    // Which neighbors should be iterated over?
                    int ItBounds = FindItBounds(RankX,RankY,MyXSlices,MyYSlices);
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
                        if ((CellType(NeighborD3D1ConvPosition) != Solid)&&(CellType(NeighborD3D1ConvPosition) != Wall)) {
                            LCount++;
                        }
                    }
                    if (LCount == 0) {
                        // This cell is returned to solid type
                        CellType(D3D1ConvPosition) = Solid;
                        //cout << RankX << " " << RankY << " " << RankZ << " is solid" << endl;
                    }
                }
            }
        }
    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (id == 0) cout << "D" << endl;
    // Second rim of wall cells for spot melt scans
//    for (int k=0; k<nz; k++)  {
//        for(int i=0; i<MyXSlices; i++) {
//            for(int j=0; j<MyYSlices; j++) {
//                int GlobalX = i + MyXOffset;
//                int GlobalY = j + MyYOffset;
//                if ((GlobalX == 1)||(GlobalX == nx-2)||(GlobalY == 1)||(GlobalY == ny-2)) {
//                    int CAGridLocation = k*MyXSlices*MyYSlices + i*MyYSlices + j;
//                    CellType(CAGridLocation) = Wall;
//                    GrainID(CAGridLocation) = 0;
//                }
//            }
//        }
//    }
    
}


// After initializing grain structure and filling ghost nodes, the known potential nucleation sites are placed into the nucleation data structures
// Each nucleation event is assigned a time step, beyond which if the associated cell is not solid or actve, the event occurs
// This data is synced across MPI ranks, for nucleation events that occur in the ghost nodes
void NucleiInit(int DecompositionStrategy, int MyXSlices, int MyYSlices, int nz, int id, double dTN, double dTsigma, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, int &PossibleNuclei_ThisRank, ViewI::HostMirror NucleiLocation, ViewI::HostMirror NucleationTimes, int* GrainOrientation, ViewI::HostMirror CellType, ViewI::HostMirror GrainID, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange) {

    // Counts and buffers for sending/recieving nucleation data from ghost nodes
    int ACount = 0;
    int BCount = 0;
    int CCount = 0;
    int DCount = 0;
    int ECount = 0;
    int FCount = 0;
    int GCount = 0;
    int HCount = 0;
    vector <int> ANuc, BNuc, CNuc, DNuc, ENuc, FNuc, GNuc, HNuc;

    // Gaussian distribution of nucleation undercooling
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(dTN,dTsigma);
    for (int i=0; i<120*id; i++) {
        double Discard = distribution(generator);
    }
    
    // Collect data for ghost nodes' nucleation events
    int NEvent = 0;
    for (int i=0; i<MyXSlices*MyYSlices*nz; i++)  {
        if (CellType(i) == LiqSol) {
            NucleiLocation(NEvent) = i;
            // Undercooling for this nucleation event
            double LocNucUnd = distribution(generator);
            // Time steps to reach this undercooling after cell goes below the liquidus
            int TimeToNucUnd = CritTimeStep(i) + round(LocNucUnd/UndercoolingChange(i));
            NucleationTimes(NEvent) = max(CritTimeStep(i),TimeToNucUnd);
            //cout << "Nuc ID " << id << " Cell " << i << " Time " <<  CritTimeStep(i) << " Adj " << round(LocNucUnd/UndercoolingChange(i)) << " or " << CritTimeStep(i) << endl;
            // Determine if other MPI ranks need information about this potential nucleation event
            // If so, store the location (X,Y,Z), GrainID, and nucleation time step value to be sent
            int RankZ = floor(i/(MyXSlices*MyYSlices));
            int Rem = i % (MyXSlices*MyYSlices);
            int RankX = floor(Rem/MyYSlices);
            int RankY = Rem % MyYSlices;
            if (DecompositionStrategy == 1) {
                if (RankY == 1) {
                    ACount++;
                    ANuc.push_back(RankX);
                    ANuc.push_back(RankZ);
                    ANuc.push_back(TimeToNucUnd);
                    ANuc.push_back(GrainID(i));
                }
                else if (RankY == MyYSlices-2) {
                    BCount++;
                    BNuc.push_back(RankX);
                    BNuc.push_back(RankZ);
                    BNuc.push_back(TimeToNucUnd);
                    BNuc.push_back(GrainID(i));
                }
            }
            else {
                if (RankY == 1) {
                    // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                    if (RankX == MyXSlices-2) {
                        
                        ECount++;
                        ENuc.push_back(RankZ);
                        ENuc.push_back(TimeToNucUnd);
                        ENuc.push_back(GrainID(i));
                        
                        CCount++;
                        CNuc.push_back(RankY);
                        CNuc.push_back(RankZ);
                        CNuc.push_back(TimeToNucUnd);
                        CNuc.push_back(GrainID(i));
                        
                        ACount++;
                        ANuc.push_back(RankX);
                        ANuc.push_back(RankZ);
                        ANuc.push_back(TimeToNucUnd);
                        ANuc.push_back(GrainID(i));
                    }
                    else if (RankX == 1) {
                        
                        GCount++;
                        GNuc.push_back(RankZ);
                        GNuc.push_back(TimeToNucUnd);
                        GNuc.push_back(GrainID(i));
                        
                        DCount++;
                        DNuc.push_back(RankY);
                        DNuc.push_back(RankZ);
                        DNuc.push_back(TimeToNucUnd);
                        DNuc.push_back(GrainID(i));
                        
                        ACount++;
                        ANuc.push_back(RankX);
                        ANuc.push_back(RankZ);
                        ANuc.push_back(TimeToNucUnd);
                        ANuc.push_back(GrainID(i));
                    }
                    else if ((RankX > 1)&&(RankX < MyXSlices-2)) {
                        // This is being sent to MyLeft
                        ACount++;
                        ANuc.push_back(RankX);
                        ANuc.push_back(RankZ);
                        ANuc.push_back(TimeToNucUnd);
                        ANuc.push_back(GrainID(i));
                    }
                }
                else if (RankY == MyYSlices-2) {
                    // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                    if (RankX == MyXSlices-2) {
                        FCount++;
                        FNuc.push_back(RankZ);
                        FNuc.push_back(TimeToNucUnd);
                        FNuc.push_back(GrainID(i));
                        
                        CCount++;
                        CNuc.push_back(RankY);
                        CNuc.push_back(RankZ);
                        CNuc.push_back(TimeToNucUnd);
                        CNuc.push_back(GrainID(i));
                        
                        BCount++;
                        BNuc.push_back(RankX);
                        BNuc.push_back(RankZ);
                        BNuc.push_back(TimeToNucUnd);
                        BNuc.push_back(GrainID(i));
                    }
                    else if (RankX == 1) {
                        HCount++;
                        HNuc.push_back(RankZ);
                        HNuc.push_back(TimeToNucUnd);
                        HNuc.push_back(GrainID(i));
                        
                        DCount++;
                        DNuc.push_back(RankY);
                        DNuc.push_back(RankZ);
                        DNuc.push_back(TimeToNucUnd);
                        DNuc.push_back(GrainID(i));
                        
                        BCount++;
                        BNuc.push_back(RankX);
                        BNuc.push_back(RankZ);
                        BNuc.push_back(TimeToNucUnd);
                        BNuc.push_back(GrainID(i));
                    }
                    else if ((RankX > 1)&&(RankX < MyXSlices-2)) {
                        BCount++;
                        BNuc.push_back(RankX);
                        BNuc.push_back(RankZ);
                        BNuc.push_back(TimeToNucUnd);
                        BNuc.push_back(GrainID(i));
                    }
                }
                else if ((RankX == 1)&&(RankY > 1)&&(RankY < MyYSlices-2)) {
                    DCount++;
                    DNuc.push_back(RankY);
                    DNuc.push_back(RankZ);
                    DNuc.push_back(TimeToNucUnd);
                    DNuc.push_back(GrainID(i));
                }
                else if ((RankX == MyXSlices-2)&&(RankY > 1)&&(RankY < MyYSlices-2)) {
                    CCount++;
                    CNuc.push_back(RankY);
                    CNuc.push_back(RankZ);
                    CNuc.push_back(TimeToNucUnd);
                    CNuc.push_back(GrainID(i));
                }
            }
            NEvent++;
        }
    }
    //cout << "RANK " << id << " starting placement at " << NEvent << " last position is " << PossibleNuclei_ThisRank-1 << endl;
    // Determine whether or not ghost node information transfer needs to take place
    int ARCount, BRCount, CRCount, DRCount, ERCount, FRCount, GRCount, HRCount;

    // Send BCount, Recieve ARCount (send to the right, recieve on the left)
    MPI_Sendrecv(&BCount,1,MPI_INT,MyRight,0,&ARCount,1,MPI_INT,MyLeft,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    
    // Send ACount, Recieve BRCount (send to the left, recieve on the right)
    MPI_Sendrecv(&ACount,1,MPI_INT,MyLeft,0,&BRCount,1,MPI_INT,MyRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    
    if (DecompositionStrategy != 1) {

        // Send CCount, Recieve DRCount (send into the plane, recieve out of the plane)
        MPI_Sendrecv(&CCount,1,MPI_INT,MyIn,0,&DRCount,1,MPI_INT,MyOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        
        // Send DCount, Recieve CRCount (send out of the plane, recieve into the plane)
        MPI_Sendrecv(&DCount,1,MPI_INT,MyOut,0,&CRCount,1,MPI_INT,MyIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        
        
        // Send HCount, Recieve ERCount
        MPI_Sendrecv(&HCount,1,MPI_INT,MyRightOut,0,&ERCount,1,MPI_INT,MyLeftIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        
        // Send ECount, Recieve HRCount
        MPI_Sendrecv(&ECount,1,MPI_INT,MyLeftIn,0,&HRCount,1,MPI_INT,MyRightOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        
        // Send GCount, Recieve FRCount
        MPI_Sendrecv(&GCount,1,MPI_INT,MyLeftOut,0,&FRCount,1,MPI_INT,MyRightIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        
        // Send FCount, Recieve GRCount
        MPI_Sendrecv(&FCount,1,MPI_INT,MyRightIn,0,&GRCount,1,MPI_INT,MyLeftOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
        
    if (MyLeft == MPI_PROC_NULL) ARCount = 0;
    if (MyRight == MPI_PROC_NULL) BRCount = 0;
    if (MyIn == MPI_PROC_NULL) CRCount = 0;
    if (MyOut == MPI_PROC_NULL) DRCount = 0;
    if (MyLeftIn == MPI_PROC_NULL) ERCount = 0;
    if (MyRightIn == MPI_PROC_NULL) FRCount = 0;
    if (MyLeftOut == MPI_PROC_NULL) GRCount = 0;
    if (MyRightOut == MPI_PROC_NULL) HRCount = 0;
    
    //cout << "ID is " << id << " has " << ARCount << " and " << BRCount << endl;
    
    // Buffers for recieving ghost node data
    int GhostNodesAR[4*ARCount];
    int GhostNodesBR[4*BRCount];
    int GhostNodesCR[4*CRCount];
    int GhostNodesDR[4*DRCount];
    int GhostNodesER[3*ERCount];
    int GhostNodesFR[3*FRCount];
    int GhostNodesGR[3*GRCount];
    int GhostNodesHR[3*HRCount];
    
    //MPI_Barrier(MPI_COMM_WORLD);
    //cout << "ID = " << id << " A to send " << ACount << " B to send " << BCount << " A to recieve " << ARCount << " B to recieve " << BRCount << endl;
    // Collect ghost node data and send to other ranks- left and right
    if (ACount > 0) {
        int GhostNodesA[4*ACount];
        for (int i=0; i<4*ACount; i++) {
            GhostNodesA[i] = ANuc[i];
            //cout << "Value " << i << " from rank " << id << " is " << GhostNodesA[i] << endl;
        }
        if (BRCount == 0) {
            // Sending data to id = id - 1 only
            MPI_Send(&GhostNodesA,ACount*4,MPI_INT,MyLeft,0,MPI_COMM_WORLD);
        }
        else {
            // Sending data to id = id - 1 and recieving data from id = id + 1
            MPI_Sendrecv(&GhostNodesA,ACount*4,MPI_INT,MyLeft,0,&GhostNodesBR,BRCount*4,MPI_INT,MyRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        // cout << "ID = " << id << " ACount = " << ACount << " ABuf = " << BufACount << endl;
    }
    else if (BRCount > 0) {
        // Recieving data from id = id + 1 only
        MPI_Recv(&GhostNodesBR,BRCount*4,MPI_INT,MyRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    //cout << "Collect B Start" << endl;
   // if (id == 0) cout << " BCOUNT = " << BCount << endl;
    if (BCount > 0) {
        int GhostNodesB[4*BCount];
        for (int i=0; i<4*BCount; i++) {
            GhostNodesB[i] = BNuc[i];
        }
        if (ARCount == 0) {
            // Sending data to id = id + 1 only
            MPI_Send(&GhostNodesB,BCount*4,MPI_INT,MyRight,1,MPI_COMM_WORLD);
            //if (id == 0) cout << " Rank 0 sent data starting with " << GhostNodesB[0] << endl;
        }
        else {
            // Sending data to id = id + 1 and recieving data from id = id - 1
            MPI_Sendrecv(&GhostNodesB,BCount*4,MPI_INT,MyRight,1,&GhostNodesAR,ARCount*4,MPI_INT,MyLeft,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            //if (id == 1) cout << " Rank 1 recieved data starting with " << GhostNodesAR[0] << endl;
        }
    }
    else if (ARCount > 0) {
        // Recieving data from id = id - 1 only
        MPI_Recv(&GhostNodesAR,ARCount*4,MPI_INT,MyLeft,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        //if (id == 1) cout << " Rank 1 recieved data starting with " << GhostNodesAR[0] << endl;
    }


    if (DecompositionStrategy != 1) {
        // Collect ghost node data and send to other ranks- in and out
        if (CCount > 0) {
            int GhostNodesC[4*CCount];
            for (int i=0; i<4*CCount; i++) {
                GhostNodesC[i] = CNuc[i];
            }
            
            if (DRCount == 0) {
                // Sending data only
                MPI_Send(&GhostNodesC,CCount*4,MPI_INT,MyIn,0,MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
                MPI_Sendrecv(&GhostNodesC,CCount*4,MPI_INT,MyIn,0,&GhostNodesDR,DRCount*4,MPI_INT,MyOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        else if (DRCount > 0) {
            // Recieving data only
            MPI_Recv(&GhostNodesDR,DRCount*4,MPI_INT,MyOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        
        if (DCount > 0) {
            int GhostNodesD[4*DCount];
            for (int i=0; i<4*DCount; i++) {
                GhostNodesD[i] = DNuc[i];
            }
            if (CRCount == 0) {
                // Sending data only
                MPI_Send(&GhostNodesD,DCount*4,MPI_INT,MyOut,1,MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
                MPI_Sendrecv(&GhostNodesD,DCount*4,MPI_INT,MyOut,1,&GhostNodesCR,CRCount*4,MPI_INT,MyIn,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            //    cout << "ID = " << id << " DCount = " << DCount << " DBuf = " << BufDCount << endl;
        }
        else if (CRCount > 0) {
            // Recieving data only
            MPI_Recv(&GhostNodesCR,CRCount*4,MPI_INT,MyIn,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        
        // Collect ghost node data and send to other ranks- MyLeftIn and MyRightOut
        if (ECount > 0) {
            int GhostNodesE[3*ECount];
            for (int i=0; i<3*ECount; i++) {
                GhostNodesE[i] = ENuc[i];
            }
            if (HRCount == 0) {
                // Sending data only
                MPI_Send(&GhostNodesE,ECount*3,MPI_INT,MyLeftIn,0,MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
                MPI_Sendrecv(&GhostNodesE,ECount*3,MPI_INT,MyLeftIn,0,&GhostNodesHR,HRCount*3,MPI_INT,MyRightOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            //  cout << "ID = " << id << " ECount = " << ECount << " EBuf = " << BufECount << endl;
        }
        else if (HRCount > 0) {
            // Recieving data only
            MPI_Recv(&GhostNodesHR,HRCount*3,MPI_INT,MyRightOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        
        if (HCount > 0) {
            int GhostNodesH[3*HCount];
            for (int i=0; i<3*HCount; i++) {
                GhostNodesH[i] = HNuc[i];
            }
            if (ERCount == 0) {
                // Sending data only
                MPI_Send(&GhostNodesH,HCount*3,MPI_INT,MyRightOut,0,MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
                MPI_Sendrecv(&GhostNodesH,HCount*3,MPI_INT,MyRightOut,0,&GhostNodesER,ERCount*3,MPI_INT,MyLeftIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        else if (ERCount > 0) {
            // Recieving data only
            MPI_Recv(&GhostNodesER,ERCount*3,MPI_INT,MyLeftIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        
        // Collect ghost node data and send to other ranks- MyRightIn and MyLeftOut
        if (FCount > 0) {
            int GhostNodesF[4*FCount];
            for (int i=0; i<4*FCount; i++) {
                GhostNodesF[i] = FNuc[i];
            }
            if (GRCount == 0) {
                // Sending data only
                MPI_Send(&GhostNodesF,FCount*3,MPI_INT,MyRightIn,1,MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
                MPI_Sendrecv(&GhostNodesF,FCount*3,MPI_INT,MyRightIn,1,&GhostNodesGR,GRCount*3,MPI_INT,MyLeftOut,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        else if (GRCount > 0) {
            // Recieving data only
            MPI_Recv(&GhostNodesGR,GRCount*3,MPI_INT,MyLeftOut,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        
        if (GCount > 0) {
            int GhostNodesG[3*GCount];
            for (int i=0; i<3*GCount; i++) {
                GhostNodesG[i] = GNuc[i];
            }
            if (FRCount == 0) {
                // Sending data only
                MPI_Send(&GhostNodesG,GCount*3,MPI_INT,MyLeftOut,1,MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
                MPI_Sendrecv(&GhostNodesG,GCount*3,MPI_INT,MyLeftOut,1,&GhostNodesFR,FRCount*3,MPI_INT,MyRightIn,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        else if (FRCount > 0) {
            // Recieving data only
            MPI_Recv(&GhostNodesFR,FRCount*3,MPI_INT,MyRightIn,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    }
   // MPI_Barrier(MPI_COMM_WORLD);
   // cout << "Collect Done" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    // Place ghost node data recieved from the left (if needed)
    if (ARCount > 0) {
        for (int i=0; i<ARCount; i++) {
            int RankX = GhostNodesAR[4*i];
            int RankY = 0;
            int RankZ = GhostNodesAR[4*i+1];
            int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
            NucleiLocation(NEvent) = CellLocation;
            NucleationTimes(NEvent) = GhostNodesAR[4*i+2];
            CellType(CellLocation) = LiqSol;
            GrainID(CellLocation) = GhostNodesAR[4*i+3];
            //            cout << "GN Nuc ID " << id << " Cell " << CellLocation << " Time " << GhostNodesAR[4*i+2] << " GID " << GhostNodesAR[4*i+3] << endl;
            NEvent++;
        }
    }

    // Place ghost node data recieved from the right (if needed)
    if (BRCount > 0) {
        for (int i=0; i<BRCount; i++) {
            int RankX = GhostNodesBR[4*i];
            int RankY = MyYSlices-1;
            int RankZ = GhostNodesBR[4*i+1];
            int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
            NucleiLocation(NEvent) = CellLocation;
            NucleationTimes(NEvent) = GhostNodesBR[4*i+2];
            CellType(CellLocation) = LiqSol;
            GrainID(CellLocation) = GhostNodesBR[4*i+3];
            //            cout << "GN Nuc ID " << id << " Cell " << CellLocation << " Time " << GhostNodesAR[4*i+2] << " GID " << GhostNodesAR[4*i+3] << endl;
            NEvent++;
        }
    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    cout << "Data placed" << endl;
    if (DecompositionStrategy != 1) {
        // Place ghost node data recieved from in plane (if needed)
        if (CRCount > 0) {
            for (int i=0; i<CRCount; i++) {
                int RankX = MyXSlices-1;
                int RankY = GhostNodesCR[4*i];
                int RankZ = GhostNodesCR[4*i+1];
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                NucleiLocation(NEvent) = CellLocation;
                NucleationTimes(NEvent) = GhostNodesCR[4*i+2];
                CellType(CellLocation) = LiqSol;
                GrainID(CellLocation) = GhostNodesCR[4*i+3];
                //            cout << "GN Nuc ID " << id << " Cell " << CellLocation << " Time " << GhostNodesAR[4*i+2] << " GID " << GhostNodesAR[4*i+3] << endl;
                NEvent++;
            }
        }

        // Place ghost node data recieved from out of plane (if needed)
        if (DRCount > 0) {
            for (int i=0; i<DRCount; i++) {
                int RankX = 0;
                int RankY = GhostNodesDR[4*i];
                int RankZ = GhostNodesDR[4*i+1];
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                NucleiLocation(NEvent) = CellLocation;
                NucleationTimes(NEvent) = GhostNodesDR[4*i+2];
                CellType(CellLocation) = LiqSol;
                GrainID(CellLocation) = GhostNodesDR[4*i+3];
                      //      cout << "GN Nuc ID " << id << " Cell " << CellLocation << " Time " << GhostNodesAR[4*i+2] << " GID " << GhostNodesAR[4*i+3] << endl;
                NEvent++;
            }
        }

        // Place ghost node data recieved from left and out of plane (if needed)
        if (ERCount > 0) {
            for (int i=0; i<ERCount; i++) {
                int RankX = MyXSlices-1;
                int RankY = 0;
                int RankZ = GhostNodesER[3*i];
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                NucleiLocation(NEvent) = CellLocation;
                NucleationTimes(NEvent) = GhostNodesER[3*i+1];
                CellType(CellLocation) = LiqSol;
                GrainID(CellLocation) = GhostNodesER[3*i+2];
                   //         cout << "GN Nuc ID " << id << " Cell " << CellLocation << " Time " << GhostNodesAR[4*i+2] << " GID " << GhostNodesAR[4*i+3] << endl;
                NEvent++;
            }
        }

        // Place ghost node data recieved from right and out of plane (if needed)
        if (FRCount > 0) {
            for (int i=0; i<FRCount; i++) {
                int RankX = MyXSlices-1;
                int RankY = MyYSlices-1;
                int RankZ = GhostNodesFR[3*i];
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                NucleiLocation(NEvent) = CellLocation;
                NucleationTimes(NEvent) = GhostNodesFR[3*i+1];
                CellType(CellLocation) = LiqSol;
                GrainID(CellLocation) = GhostNodesFR[3*i+2];
                  //          cout << "GN Nuc ID " << id << " Cell " << CellLocation << " Time " << GhostNodesAR[4*i+2] << " GID " << GhostNodesAR[4*i+3] << endl;
                NEvent++;
            }
        }

        // Place ghost node data recieved from left and into plane (if needed)
        if (GRCount > 0) {
            for (int i=0; i<GRCount; i++) {
                int RankX = 0;
                int RankY = 0;
                int RankZ = GhostNodesGR[3*i];
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                NucleiLocation(NEvent) = CellLocation;
                NucleationTimes(NEvent) = GhostNodesGR[3*i+1];
                CellType(CellLocation) = LiqSol;
                GrainID(CellLocation) = GhostNodesGR[3*i+2];
                    //        cout << "GN Nuc ID " << id << " Cell " << CellLocation << " Time " << GhostNodesAR[4*i+2] << " GID " << GhostNodesAR[4*i+3] << endl;
                NEvent++;
            }
        }

        if (HRCount > 0) {
            for (int i=0; i<HRCount; i++) {
                int RankX = 0;
                int RankY = MyYSlices-1;
                int RankZ = GhostNodesHR[3*i];
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                NucleiLocation(NEvent) = CellLocation;
                NucleationTimes(NEvent) = GhostNodesHR[3*i+1];
                CellType(CellLocation) = LiqSol;
                GrainID(CellLocation) = GhostNodesHR[3*i+2];
                    //        cout << "GN Nuc ID " << id << " Cell " << CellLocation << " Time " << GhostNodesAR[4*i+2] << " GID " << GhostNodesAR[4*i+3] << endl;
                NEvent++;
            }
        }
    }
}


void LayerSetup(string SubstrateFileName, int layernumber, int LayerHeight, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int nz, int LocalDomainSize, int id, int np, ViewF::HostMirror DiagonalLength, ViewI::HostMirror CellType, ViewI::HostMirror TriangleIndex, ViewI::HostMirror GrainID, int* GrainID_Stored, ViewF::HostMirror CritDiagonalLength, ViewF::HostMirror DOCenter, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange, ViewF::HostMirror UndercoolingCurrent, bool* Melted, bool* Melted_Stored, bool LayerwiseTemeperature, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferC, Buffer2D BufferD, Buffer2D BufferE, Buffer2D BufferF, Buffer2D BufferG, Buffer2D BufferH, Buffer2D BufferAR, Buffer2D BufferBR, Buffer2D BufferCR, Buffer2D BufferDR, Buffer2D BufferER, Buffer2D BufferFR, Buffer2D BufferGR, Buffer2D BufferHR, int BufSizeX, int BufSizeY, int BufSizeZ) {
    
    int LayerOffset = layernumber*LayerHeight*MyXSlices*MyYSlices;
    for (int k=1; k<=LayerHeight; k++)  {
        for(int i=0; i<MyXSlices; i++) {
            for(int j=0; j<MyYSlices; j++) {
                GrainID_Stored[LayerOffset+(k-1)*MyXSlices*MyYSlices+i*MyYSlices+j] = GrainID(k*MyXSlices*MyYSlices+i*MyYSlices+j);
                Melted_Stored[LayerOffset+(k-1)*MyXSlices*MyYSlices+i*MyYSlices+j] = Melted[k*MyXSlices*MyYSlices+i*MyYSlices+j];
            }
        }
    }
    
    int* GrainID_Temp = new int[nz*MyXSlices*MyYSlices];
    bool* Melted_Temp = new bool[nz*MyXSlices*MyYSlices];
    for (int i=0; i<nz*MyXSlices*MyYSlices; i++) {
        GrainID_Temp[i] = GrainID(i);
        Melted_Temp[i] = Melted[i];
    }
    
    // Shift substrate grain IDs for new layer
    for (int k=1; k<=nz-LayerHeight-2; k++)  {
        for(int i=0; i<MyXSlices; i++) {
            for(int j=0; j<MyYSlices; j++) {
                // Shift substrate in unmelted region
                GrainID(k*MyXSlices*MyYSlices+i*MyYSlices+j) = GrainID_Temp[(k+LayerHeight)*MyXSlices*MyYSlices+i*MyYSlices+j];
                Melted[k*MyXSlices*MyYSlices+i*MyYSlices+j] = Melted_Temp[(k+LayerHeight)*MyXSlices*MyYSlices+i*MyYSlices+j];
            }
        }
    }
    delete [] GrainID_Temp;
    delete [] Melted_Temp;

    if (LayerwiseTemeperature) {
        // Clear undercooling/critical time step values from past layer
        for (int i=0; i<LocalDomainSize; i++) {
            UndercoolingChange(i) = 0.0;
            UndercoolingCurrent(i) = 0.0;
            CritTimeStep(i) = 0;
        }
    }
    for (int i=0; i<LocalDomainSize; i++) {
        CellType(i) = 0.0;
        DiagonalLength(i) = 0.0;
    }
    for (int i=0; i<3*LocalDomainSize; i++) {
        DOCenter(i) = 0.0;
    }
    for (int i=0; i<26*LocalDomainSize; i++) {
        CritDiagonalLength(i) = 0.0;
    }
    for (int i=0; i<72*LocalDomainSize; i++) {
        TriangleIndex(i) = 0.0;
    }
    
    // Reset buffers
    Kokkos::parallel_for ("XZBufReset",BufSizeX*BufSizeZ, KOKKOS_LAMBDA (const int& i) {
        for (int j=0; j<5; j++) {
            BufferA(i,j) = 0.0;
            BufferAR(i,j) = 0.0;
            BufferB(i,j) = 0.0;
            BufferBR(i,j) = 0.0;
        }
    });
    Kokkos::parallel_for ("YZBufReset",BufSizeY*BufSizeZ, KOKKOS_LAMBDA (const int& i) {
        for (int j=0; j<5; j++) {
            BufferC(i,j) = 0.0;
            BufferCR(i,j) = 0.0;
            BufferD(i,j) = 0.0;
            BufferDR(i,j) = 0.0;
        }
    });
    Kokkos::parallel_for ("ZBufReset",BufSizeZ, KOKKOS_LAMBDA (const int& i) {
        for (int j=0; j<5; j++) {
            BufferE(i,j) = 0.0;
            BufferER(i,j) = 0.0;
            BufferF(i,j) = 0.0;
            BufferFR(i,j) = 0.0;
            BufferG(i,j) = 0.0;
            BufferGR(i,j) = 0.0;
            BufferH(i,j) = 0.0;
            BufferHR(i,j) = 0.0;
        }
    });
    
}
