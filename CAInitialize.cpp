#include "header.h"
using namespace std;
// Initializes input parameters, mesh, temperature field, and grain structures for CA simulations

void MasterInputRead(int &DecompositionStrategy, double &deltax, double &AConst, double &BConst, double &CConst, double &NMax, double &dTN, double &dTsigma, string &BaseFileName, string &GrainOrientationFile, string &TemperatureDataType, double &InitialGrainWidth) {
    
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
    // Initial grain width
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    InitialGrainWidth = atof(str2.c_str());
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
    // Output base file name
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    BaseFileName = ValueRead.substr(found+1,string::npos);
    // File of grain orientations
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    GrainOrientationFile = ValueRead.substr(found+1,string::npos);
    InputData.close();

}

void AInputRead(string &AFile1, string &AFile2, string &AFile3, int &NRatio, int &NumberOfLayers, int &LayerHeight) {

    ifstream InputData;
    string Colon = ":";
    InputData.open("AInputs.txt");
    bool SkippingLines = true;
    while(SkippingLines) {
        string dummyline;
        getline(InputData,dummyline);
        if (dummyline == "*****")
            SkippingLines = false;
    }
    string ValueRead;
    // Data file containing liquidus times
    getline(InputData,ValueRead);
    std::size_t found = ValueRead.find(Colon);
    AFile1 = ValueRead.substr(found+1,string::npos);
    // Data file containing cooling rates
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    AFile2 = ValueRead.substr(found+1,string::npos);
    // Data file containing other heat transport data
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    AFile3 = ValueRead.substr(found+1,string::npos);
    // Ratio of deltax to deltat
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    string str2 = ValueRead.substr(found+1,string::npos);
    NRatio = stoi(str2,nullptr,10);
    // Number of layers
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    NumberOfLayers = stoi(str2,nullptr,10);
    // Layer height
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    LayerHeight = stoi(str2,nullptr,10);
    InputData.close();
}

void CInputRead(double &G, double &R, int &nx, int &ny, int &nz, int &NRatio, int &NumberOfLayers) {
    
    ifstream InputData;
    string Colon = ":";
    InputData.open("CInputs.txt");
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
    // Domain size in y
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    ny = stoi(str2,nullptr,10);
    // Domain size in z
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    nz = stoi(str2,nullptr,10);
    nz = nz+2;
    // Ratio of deltax to deltat
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    NRatio = stoi(str2,nullptr,10);
    InputData.close();
    InputData.close();
    NumberOfLayers = 1;
}

void UInputRead(double &G, double &R, double &DomUndercooling, int &nx, int &ny, int &nz, int &NumberOfLayers) {
    
    ifstream InputData;
    string Colon = ":";
    InputData.open("UInputs.txt");
    bool SkippingLines = true;
    while(SkippingLines) {
        string dummyline;
        getline(InputData,dummyline);
        if (dummyline == "*****")
            SkippingLines = false;
    }
    string ValueRead;
    // Domain Undercooling (constant)
    getline(InputData,ValueRead);
    std::size_t found = ValueRead.find(Colon);
    std::string str2 = ValueRead.substr(found+1,string::npos);
    DomUndercooling = atof(str2.c_str());
    // Domain size in x
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    nx = stoi(str2,nullptr,10);
    // Domain size in y
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    ny = stoi(str2,nullptr,10);
    // Domain size in z
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
    str2 = ValueRead.substr(found+1,string::npos);
    nz = stoi(str2,nullptr,10);
    nz = nz+2;
    InputData.close();
    G = 0;
    R = 0;
    NumberOfLayers = 1;
}

void RInputRead(string &tempfile, double &HT_deltax, double &deltat, int &NumberOfLayers, int &LayerOffset, string &BaseFileName) {
    
    ifstream InputData;
    string Colon = ":";
    InputData.open("RInputs.txt");
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
    std::size_t found = ValueRead.find(Colon);
    tempfile = ValueRead.substr(found+1,string::npos);
    // Heat transport mesh size
    getline(InputData,ValueRead);
    found = ValueRead.find(Colon);
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
    InputData.close();
    
}

void ParallelMeshInit(double &G, double &R, int DecompositionStrategy, int (&NeighborX)[26], int (&NeighborY)[26], int (&NeighborZ)[26], int (&ItList)[9][26], string TemperatureDataType, int ierr, int id, int np, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset,int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn, int &MyLeftOut, int &MyRightIn, int &MyRightOut, double &deltax, double HT_deltax, double &deltat, int &nx, int &ny, int &nz, int &ProcessorsInXDirection, int &ProcessorsInYDirection, string AFile3, string tempfile, float &XMin, float &XMax, float &YMin, float &YMax, float &ZMin, float &ZMax, double DomUndercooling, double NRatio) {
        
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
    if (TemperatureDataType == "A") {
        // Open file containing analytical solution run data
        ifstream TFieldData;
        TFieldData.open(AFile3);
        string Colon = ":";
        string ValueRead;
        
        // Effective input power
        getline(TFieldData,ValueRead);
        std::size_t found = ValueRead.find(Colon);
        std::string str2 = ValueRead.substr(found+1,string::npos);
        double InPower = atof(str2.c_str());
        // Beam velocity (m/s)
        getline(TFieldData,ValueRead);
        found = ValueRead.find(Colon);
        str2 = ValueRead.substr(found+1,string::npos);
        double BeamVelocity =  atof(str2.c_str());
        // Substrate preheat temperature
        getline(TFieldData,ValueRead);
        found = ValueRead.find(Colon);
        str2 = ValueRead.substr(found+1,string::npos);
        double TPreheat =  atof(str2.c_str());
        // Thermal conductivity of substrate
        getline(TFieldData,ValueRead);
        found = ValueRead.find(Colon);
        str2 = ValueRead.substr(found+1,string::npos);
        double ThermalCond = atof(str2.c_str());
        // Thermal diffusivity of substrate
        getline(TFieldData,ValueRead);
        found = ValueRead.find(Colon);
        str2 = ValueRead.substr(found+1,string::npos);
        double ThermalDiff = atof(str2.c_str());
        // Cell size (in m)
        getline(TFieldData,ValueRead);
        found = ValueRead.find(Colon);
        str2 = ValueRead.substr(found+1,string::npos);
        deltax = atof(str2.c_str())*pow(10,-6);
        // Domain size in x
        getline(TFieldData,ValueRead);
        found = ValueRead.find(Colon);
        str2 = ValueRead.substr(found+1,string::npos);
        nx = stoi(str2,nullptr,10);
        // Domain size in y
        getline(TFieldData,ValueRead);
        found = ValueRead.find(Colon);
        str2 = ValueRead.substr(found+1,string::npos);
        ny = stoi(str2,nullptr,10);
        // Domain size in z
        getline(TFieldData,ValueRead);
        found = ValueRead.find(Colon);
        str2 = ValueRead.substr(found+1,string::npos);
        nz = stoi(str2,nullptr,10);
        nz = nz+2;
        // Close file
        TFieldData.close();
        
        // deltat based on the known deltax and beam velocity
        deltat = deltax/(BeamVelocity*NRatio);
        
        if (id == 0) {
            cout << "CA model of analytical heat transport problem solidification" << endl;
            cout << "Effective input power (W/m2): " << InPower << endl;
            cout << "Beam velocity (m/s): " << BeamVelocity << endl;
            cout << "Preheat temperature (K): " << TPreheat << endl;
            cout << "Thermal conductivity: " << ThermalCond << endl;
            cout << "Thermal diffusivity: " << ThermalDiff << endl;
            cout << "Cell size (m): " << deltax << endl;
            cout << "Time step (s): " << deltat << endl;
            cout << "Domain size: " << nx << " by " << ny << " by " << nz << endl;
            cout << "================================================================" << endl;
        }
    }
    else if (TemperatureDataType == "C") {
        
        // deltat based on the known deltax, G, and R
        deltat = deltax/((R/G)*NRatio);
        
        if (id == 0) {
            cout << "CA model of constrained solidification problem" << endl;
            cout << "Thermal gradient (K/m): " << G << endl;
            cout << "Cooling rate (K/s): " << R << endl;
            cout << "Cell size (m): " << deltax << endl;
            cout << "Time step (s): " << deltat << endl;
            cout << "Domain size: " << nx << " by " << ny << " by " << nz << endl;
            cout << "================================================================" << endl;
        }
    }
    
    else if (TemperatureDataType == "U") {
        
        deltat = deltax;
        
        if (id == 0) {
            cout << "CA model of equiaxed grain growth" << endl;
            cout << "Domain Undercooling (K): " << DomUndercooling << endl;
            cout << "Cell size (m): " << deltax << endl;
            cout << "Domain size: " << nx << " by " << ny << " by " << nz << endl;
            cout << "================================================================" << endl;
        }
        
    }
    
    else {

        // Determine mesh size needed based on OpenFOAM/Truchas data
        // Read geometry data from OpenFOAM/Truchas
        int nx_HT, ny_HT, nz_HT; // OpenFOAM/Truchas mesh limits

        XMin = 10000.0;
        YMin = 10000.0;
        ZMin = 10000.0;
        XMax = -10000.0;
        YMax = -10000.0;
        ZMax = -10000.0;
        ifstream Geom;
        Geom.open(tempfile);
        while (!Geom.eof()) {
            string s;
            getline(Geom,s);
            bool ReadingLine = true;
            int i = 0;
            int j = 0;
            while (ReadingLine) {

                if (s.empty()) break;
                char C = s.at(i);
                // If this character is not a space, convert from string
                if (!isblank(C)) {
                    int FirstChar = i;
                    int LastChar;
                    bool ReadingValue = true;
                    while (ReadingValue) {
                        i++;
                        char C2 = s.at(i);
                        if (isblank(C2)) {
                            LastChar = i;
                            ReadingValue = false;
                        }
                    }
                    string NewDataS = s.substr(FirstChar,LastChar-FirstChar);
                    //cout << NewDataS << endl;
                    float MeshData = atof(NewDataS.c_str());
                    //if (id == 0) cout << MeshData << endl;
                    if (j == 0) {
                        if (XMin > MeshData) XMin = MeshData;
                        if (XMax < MeshData) XMax = MeshData;
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
                // Otherwise advance to the next character
                else {
                    i++;
                }
            }
        }

        Geom.close();

        // Is the input in m (OpenFOAM) or mm (Truchas)?
        if (TemperatureDataType == "T") {
            XMin = XMin/1000.0;
            XMax = XMax/1000.0;
            YMin = YMin/1000.0;
            YMax = YMax/1000.0;
            ZMin = ZMin/1000.0;
            ZMax = ZMax/1000.0;
        }

        nx_HT = round((XMax-XMin)/HT_deltax);
        ny_HT = round((YMax-YMin)/HT_deltax);
        nz_HT = round((ZMax-ZMin)/HT_deltax);
        
        nx = round((double)(nx_HT)*(double)(HT_deltax)/(double)(deltax) + 4.0);
        ny = round((double)(ny_HT)*(double)(HT_deltax)/(double)(deltax) + 4.0);
        nz = round((double)(nz_HT)*(double)(HT_deltax)/(double)(deltax) + 3.0);
            
        if (id == 0) {
            if (TemperatureDataType == "O") cout << "CA model of OpenFOAM heat transport problem solidification" << endl;
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

    InitialDecomposition(DecompositionStrategy, nx, ny, ProcessorsInXDirection, ProcessorsInYDirection, id, np, MyXSlices, MyYSlices, MyXOffset, MyYOffset, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyLeftOut, MyRightIn, MyRightOut);

    MyXOffset = XOffsetCalc(id,nx,ProcessorsInXDirection,ProcessorsInYDirection,np,DecompositionStrategy);
    MyXSlices = XMPSlicesCalc(id,nx,ProcessorsInXDirection,ProcessorsInYDirection,np,DecompositionStrategy);

    MyYOffset = YOffsetCalc(id,ny,ProcessorsInYDirection,np,DecompositionStrategy);
    MyYSlices = YMPSlicesCalc(id,ny,ProcessorsInYDirection,np,DecompositionStrategy);

}

void TempInit(double G, double R, int DecompositionStrategy, int (&NeighborX)[26], int (&NeighborY)[26], int (&NeighborZ)[26], int (&ItList)[9][26], string TemperatureDataType, int ierr, int id, int np, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset,int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn, int &MyLeftOut, int &MyRightIn, int &MyRightOut, double deltax, double HT_deltax, double deltat, int &nx, int &ny, int &nz, int &ProcessorsInXDirection, int &ProcessorsInYDirection, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange, ViewF::HostMirror UndercoolingCurrent, double DomUndercooling, string tempfile, string FileA1, string FileA2, float &XMin, float &XMax, float &YMin, float &YMax, float &ZMin, float &ZMax) {

   if ((TemperatureDataType == "O")||(TemperatureDataType == "T")) {
            
        // Temperature data read
        int HTtoCAratio = HT_deltax/deltax; // OpenFOAM/CA cell size ratio

        vector <vector <vector <double> > > CritTS, CritTL;
        for (int k=0; k<nz; k++) {
            vector <vector <double> > TemperatureXX;
            for (int i=0; i<nx; i++) {
                vector <double> TemperatureX;
                for (int j=0; j<ny; j++) {
                    TemperatureX.push_back(-1.0);
                }
                TemperatureXX.push_back(TemperatureX);
            }
            CritTS.push_back(TemperatureXX);
            CritTL.push_back(TemperatureXX);
        }

        double UnitConversion;
        if (TemperatureDataType == "O") UnitConversion = 1;
        else UnitConversion = 1000;
       
        ifstream Geom;
        Geom.open(tempfile);
        while (!Geom.eof()) {
            int XInt, YInt, ZInt;
            double MyX, MyY, MyZ, FirstTime, LastTime;
            string s;
            getline(Geom,s);
            bool ReadingLine = true;
            int i = 0;
            int j = 0;
            while (ReadingLine) {
                
                if (s.empty()) break;
                char C = s.at(i);
                // If this character is not a space, convert from string
                if (!isblank(C)) {
                    int FirstChar = i;
                    int LastChar;
                    bool ReadingValue = true;
                    while (ReadingValue) {
                        i++;
                        //cout << i << endl;
                        char C2 = s.at(i);
                        if ((isblank(C2))||(i == s.length()-1)) {
                            LastChar = i;
                            ReadingValue = false;
                        }
                        //cout << i << endl;
                    }
                    //cout << FirstChar << " " << LastChar << endl;
                    string NewDataS = s.substr(FirstChar,LastChar-FirstChar+1);
                    float MeshData = atof(NewDataS.c_str());
                    //if (id == 0) cout << MeshData << endl;
                    if (j == 0) {
                        MyX = atof(NewDataS.c_str())/UnitConversion; // Only divide by 1000 if from Truchas
                        XInt = round((MyX-XMin)/deltax) + 2; // + HTtoCAratio;
                        // cout << XInt << endl;
                    }
                    else if (j == 1) {
                        //if (id == 0) cout << "Y " << MeshData << endl;
                        MyY = atof(NewDataS.c_str())/UnitConversion; // Only divide by 1000 if from Truchas
                        YInt = round((MyY-YMin)/deltax) + 2;// + HTtoCAratio;
                        // cout << MyY << endl;
                    }
                    else if (j == 2) {
                        //if (id == 0) cout << "Z " << MeshData << endl;
                        MyZ = atof(NewDataS.c_str())/UnitConversion; // Only divide by 1000 if from Truchas
                        ZInt = round(abs(MyZ-ZMin)/deltax)+ 2;
                        // cout << ZInt << endl;
                    }
                    else if (j == 3) {
                        FirstTime = stod(NewDataS.c_str());
                        //if (id == 0) cout << ZInt << endl;
                        CritTL[ZInt][XInt][YInt] = FirstTime/UnitConversion; // Truchas in ms
                        //cout << CritTL[ZInt][XInt][YInt] << endl;
                    }
                    else if (j == 4) {
                        LastTime = stod(NewDataS.c_str());
                        CritTS[ZInt][XInt][YInt] = LastTime/UnitConversion; // Truchas in ms
                        // if (id == 0) cout << XInt << " " << YInt << " " << ZInt << " " << FirstTime/UnitConversion << " " << LastTime/UnitConversion << endl;
                        if (LastTime < FirstTime) {
                            cout << "Negative cooling rate! " << endl;
                            double Temp = FirstTime;
                            FirstTime = LastTime;
                            LastTime = Temp;
                            CritTS[ZInt][XInt][YInt] = LastTime/UnitConversion; // Truchas in ms
                            CritTL[ZInt][XInt][YInt] = FirstTime/UnitConversion; // Truchas in ms
                        }
                        ReadingLine = false;
                    }
                    j++;
                }
                // Otherwise advance to the next character
                else {
                    i++;
                }
            }
        }
        Geom.close();

        // Data interpolation between heat transport and CA grids, if necessary
        if (HTtoCAratio != 1) {
            for (int k=2; k<nz-1; k++) {
                int LowZ = k - ((k-2) % HTtoCAratio);
                int HighZ = LowZ + HTtoCAratio;
                double FHighZ = (double)(k - LowZ)/(double)(HTtoCAratio);
                double FLowZ = 1.0 - FHighZ;
                for (int i=2; i<nx-2; i++) {
                    int LowX =  i - ((i-2) % HTtoCAratio);
                    int HighX = LowX + HTtoCAratio;
                    double FHighX = (double)(i - LowX)/(double)(HTtoCAratio);
                    double FLowX = 1.0 - FHighX;
                    for (int j=2; j<ny-2; j++) {
                        int LowY = j - ((j-2) % HTtoCAratio);
                        int HighY = LowY + HTtoCAratio;
                        double FHighY = (float)(j - LowY)/(float)(HTtoCAratio);
                        double FLowY = 1.0 - FHighY;
                        //cout << LowZ << " " << HighZ << " " << LowX << " " << HighX << " " << LowY << " " << HighY << endl;
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
                        //cout << LowZ << " " << HighZ << " " << LowX << " " << HighX << " " << LowY << " " << HighY << endl;
                    }
                }
            }
        }

        // Convert CritTL, CritTS matrices into CritTimeStep and UndercoolingChange (change in undercooling with time step)
        int MCTS = 0;
        int MinCTS = 10000000;
        for (int k=0; k<nz; k++) {
            for (int ii=0; ii<nx; ii++) {
                if ((ii >= MyXOffset)&&(ii < MyXOffset+MyXSlices)) {
                    for (int jj=0; jj<ny; jj++) {
                        if ((jj >= MyYOffset)&&(jj < MyYOffset+MyYSlices)) {
                            int Adj_i = ii - MyXOffset;
                            int Adj_j = jj - MyYOffset;
                            double CTLiq = CritTL[k][ii][jj];
                            int Coord3D1D = k*MyXSlices*MyYSlices + Adj_i*MyYSlices + Adj_j;
                            if (CTLiq > 0)  {
                                CritTimeStep(Coord3D1D) = round(CTLiq/deltat);
                                UndercoolingChange(Coord3D1D) = (1610.0-1420.0)*deltat/(CritTS[k][ii][jj] - CTLiq);

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
        }
        int GMax, GMin;
        MPI_Reduce(&MCTS, &GMax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&MinCTS, &GMin, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        if (id == 0) cout << "Min Crit time step is = " << GMin << endl;
        if (id == 0) cout << "Max Crit time step is = " << GMax << endl;
    }
    else if (TemperatureDataType == "A") {
        // Using analytical solution for temperature field to drive simulation
        int LowerX = MyXOffset;
        int UpperX = MyXOffset+MyXSlices;
        int LowerY = MyYOffset;
        int UpperY = MyYOffset+MyYSlices;
       

        if (LowerX < 0) LowerX = 0;
        if (LowerY < 0) LowerY = 0;
        if (UpperX > nx) UpperX = nx;
        if (UpperY > ny) UpperY = ny;
        //cout << "ID = " << id << " LowerX = " << LowerX << " UpperX = " << UpperX << " LowerY = " << LowerY << " UpperY " << UpperY << endl;
        // Open file and read the relevant lines
        ifstream Rosenthal1, Rosenthal2;
        Rosenthal1.open(FileA1);
        Rosenthal2.open(FileA2);
       
        //long long int line = 0;
        int ActSizeZ = nz-2; // Actual size of domain in z direction without walls
        int MCTS = 0;
        int MinCTS = 10000000;
        // Ignore data at x coordinates outside of subdomain
        for (int i=0; i<LowerX*ny*ActSizeZ; i++) {
            // Ignore lines
            // This data is not needed by this rank
            string dummyline;
            getline(Rosenthal1,dummyline);
            getline(Rosenthal2,dummyline);
         }
         for (int i=LowerX; i<UpperX; i++) {
            // Ignore data at y coordinates below the region of interest
            for (int j=0; j<LowerY*ActSizeZ; j++) {
               // Ignore lines
               // This data is not needed by this rank
               string dummyline;
               getline(Rosenthal1,dummyline);
               getline(Rosenthal2,dummyline);
            }
            for (int j=LowerY*ActSizeZ; j<UpperY*ActSizeZ; j++) {
                // This data is needed by this rank
                int LocalXPosition = i - MyXOffset;
                int LocalYPosition = j/ActSizeZ - MyYOffset;
                int LocalZPosition = j % ActSizeZ + 1;
                string PointOfInterest1,PointOfInterest2;
                getline(Rosenthal1,PointOfInterest1);
                getline(Rosenthal2,PointOfInterest2);
                double PointOfInterest1_D = atof(PointOfInterest1.c_str());
                double PointOfInterest2_D = atof(PointOfInterest2.c_str());
                CritTimeStep(LocalZPosition*MyXSlices*MyYSlices + MyYSlices*LocalXPosition + LocalYPosition) = round(PointOfInterest1_D/deltat);
                if (CritTimeStep(LocalZPosition*MyXSlices*MyYSlices + MyYSlices*LocalXPosition + LocalYPosition) > MCTS) MCTS = CritTimeStep(LocalZPosition*MyXSlices*MyYSlices + MyYSlices*LocalXPosition + LocalYPosition);
                if (CritTimeStep(LocalZPosition*MyXSlices*MyYSlices + MyYSlices*LocalXPosition + LocalYPosition) < MinCTS) MinCTS = CritTimeStep(LocalZPosition*MyXSlices*MyYSlices + MyYSlices*LocalXPosition + LocalYPosition);
                if (PointOfInterest2_D > 0) {
                    UndercoolingChange(LocalZPosition*MyXSlices*MyYSlices + MyYSlices*LocalXPosition + LocalYPosition) = (1620.0-1410.0)*deltat/PointOfInterest2_D;
                }
                else {
                    UndercoolingChange(LocalZPosition*MyXSlices*MyYSlices + MyYSlices*LocalXPosition + LocalYPosition) = 0.0;
                }
            }
            for (int j=UpperY*ActSizeZ; j<ny*ActSizeZ; j++) {
                // This data is not needed by this rank
                string dummyline;
                getline(Rosenthal1,dummyline);
                getline(Rosenthal2,dummyline);
            }
       }
       int GMax, GMin;
       MPI_Reduce(&MCTS, &GMax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
       MPI_Reduce(&MinCTS, &GMin, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
       if (id == 0) cout << "Min Crit time step is = " << GMin << endl;
       if (id == 0) cout << "Max Crit time step is = " << GMax << endl;
       Rosenthal1.close();
       Rosenthal2.close();
       MPI_Barrier(MPI_COMM_WORLD);
   }
   else if (TemperatureDataType == "U") {
       // unconstrained solidification
       // Initialize constant undercooling field, no cooling rate
       for (int k=0; k<nz; k++) {
           for (int i=0; i<MyXSlices; i++) {
               for (int j=0; j<MyYSlices; j++) {
                   UndercoolingCurrent(k*MyXSlices*MyYSlices + i*MyYSlices + j) = DomUndercooling;
                   UndercoolingChange(k*MyXSlices*MyYSlices + i*MyYSlices + j) = 0;
                   CritTimeStep(k*MyXSlices*MyYSlices + i*MyYSlices + j) = 0;
               }
           }
       }
   }
   else if (TemperatureDataType == "C") {
       // Contrained solidification
       // Initialize temperature field in Z direction with thermal gradient G set in input file
       for (int k=0; k<nz; k++) {
           for (int i=0; i<MyXSlices; i++) {
               for (int j=0; j<MyYSlices; j++) {
                   UndercoolingCurrent(k*MyXSlices*MyYSlices + i*MyYSlices + j) = -(k-1)*G*deltax;
                   UndercoolingChange(k*MyXSlices*MyYSlices + i*MyYSlices + j) = R*deltat;
                   CritTimeStep(k*MyXSlices*MyYSlices + i*MyYSlices + j) = (int)((-(k-1)*G*deltax)/(R*deltat));
               }
           }
       }
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
//
///*
// Initializes cell types for the case of a surface with random orientation square grains in a melt
// with a thermal gradient in the direction opposite the surface
//*/

void ConstrainedGrains(string TemperatureDataType, int InitialGrainWidth, int NGrainOrientations, int DecompositionStrategy, int ProcessorsInXDirection, int ProcessorsInYDirection, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int id, int np, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, int ItList[9][26], int NeighborX[26], int NeighborY[26], int NeighborZ[26], vector <int> &NucLocI, vector <int> &NucLocJ, vector <int> &NucLocK, vector <int> &NucleationTimes, vector <float> &NucleationUndercooling, int* GrainOrientation, float* GrainUnitVector, ViewF::HostMirror DiagonalLength, ViewC::HostMirror CellType, ViewI::HostMirror TriangleIndex, ViewI::HostMirror GrainID, ViewF::HostMirror CritDiagonalLength, ViewF::HostMirror DOCenter, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange, double deltax, double NMax, double dTN, double dTsigma, int &NextLayer_FirstSubstrateGrainID, int &NextLayer_FirstNucleatedGrainID, int &ACount, int &BCount, int &CCount, int &DCount, int &ECount, int &FCount, int &GCount, int &HCount) {
    
    std::seed_seq seed = {765,4,1111};
    std::array<unsigned,5> sequence;
    seed.generate(sequence.begin(),sequence.end());
    mt19937_64 gen(seed);
    
    uniform_real_distribution<double> dis(0.0, 1.0);
    srand((unsigned)id);
    // Convert initial grain spacing to a grain density
    double NWall = 1.0/(InitialGrainWidth*InitialGrainWidth*InitialGrainWidth*pow(10,-18));
    double WallProb = NWall*deltax*deltax*deltax;
    double BulkProb = NMax*deltax*deltax*deltax;
    if (id == 0) cout << "Fraction of wall sites active: " << WallProb << endl;
    if (id == 0) cout << "Fraction of melt pool sites to potentially be activated: " << BulkProb << endl;
    int SubstrateGrains = 0;
    int NucleatedGrains = 0;

    if (TemperatureDataType == "U") {
        
        // Initialize cell types - randomly placed nuclei
        for (int k=0; k<nz; k++)  {
            for(int i=0; i<MyXSlices; i++) {
                for(int j=0; j<MyYSlices; j++) {
                    int GlobalX = i + MyXOffset;
                    int GlobalY = j + MyYOffset;
                    if ((GlobalX > -1)&&(GlobalX < nx)&&(GlobalY > -1)&&(GlobalY < ny)&&(k > 0)&&(k < nz-1)) {
                        if ((i > 0)&&(i < MyXSlices-1)&&(j > 0)&&(j < MyYSlices-1)) {
                            double R = dis(gen);
                            if (R < WallProb) {
                                SubstrateGrains++;
                                CellType(k*MyXSlices*MyYSlices + i*MyYSlices + j) = 'A';
                            }
                            else {
                                GrainID(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 0;
                                CellType(k*MyXSlices*MyYSlices + i*MyYSlices + j) = 'L';
                            }
                        }
                        else {
                            GrainID(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 0;
                            CellType(k*MyXSlices*MyYSlices + i*MyYSlices + j) = 'L';
                        }
                    }
                    else {
                        CellType(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 'W';
                    }
                }
            }
        }
    }
    else if (TemperatureDataType == "C") {
        
        // Initialize cell types - epitaxial grains only on the bottom surface
        for (int k=0; k<nz; k++) {
            for(int i=0; i<MyXSlices; i++) {
                for(int j=0; j<MyYSlices; j++) {
                    int GlobalX = i + MyXOffset;
                    int GlobalY = j + MyYOffset;
                    if ((GlobalX > -1)&&(GlobalX < nx)&&(GlobalY > -1)&&(GlobalY < ny)&&(k > 0)&&(k < nz-1)) {
                        if ((i > 0)&&(i < MyXSlices-1)&&(j > 0)&&(j < MyYSlices-1)&&(k == 1)) {
                        double R = dis(gen);
                        if (R < WallProb) {
                            SubstrateGrains++;
                            CellType(k*MyXSlices*MyYSlices + i*MyYSlices + j) = 'A';
                        }
                        else {
                            GrainID(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 0;
                            CellType(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 'D';
                        }
                    }
                    else {
                        GrainID(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 0;
                        CellType(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 'D';
                    }
                }
                else {
                    CellType(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 'W';
                }
             }
          }
       }
        
        // Potential nucleated grains
        for (int k=2; k<nz; k++) {
            for(int i=0; i<MyXSlices; i++) {
                for(int j=0; j<MyYSlices; j++) {
                    if ((i > 0)&&(i < MyXSlices-1)&&(j > 0)&&(j < MyYSlices-1)) {
                        double R = dis(gen);
                        if (R < BulkProb) {
                            NucleatedGrains++;
                            CellType(k*MyXSlices*MyYSlices + i*MyYSlices + j) = 'N';
                        }
                    }
                }
            }
        }
        
    }
    else {

        // Initialize cell types - epitaxial sites can be anywhere at the interface, bulk sites in melt pool area
        for (int k=0; k<nz; k++)  {
            for(int i=0; i<MyXSlices; i++) {
                for(int j=0; j<MyYSlices; j++) {
                    int GlobalX = i + MyXOffset;
                    int GlobalY = j + MyYOffset;
                    // Non-wall cells
                    if ((GlobalX > -1)&&(GlobalX < nx)&&(GlobalY > -1)&&(GlobalY < ny)&&(k > 0)&&(k < nz-1)) {
                        
                        if ((i > 0)&&(i < MyXSlices-1)&&(j > 0)&&(j < MyYSlices-1)) {
                            // Cells that are either 'D' or 'N'
                            if (CritTimeStep(k*MyXSlices*MyYSlices + i*MyYSlices + j) > 0) {
                                double R = dis(gen);
                                if (R < BulkProb) {
                                    NucleatedGrains++;
                                    CellType(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 'N';
                                }
                                else {
                                    CellType(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 'D';
                                }
                            }
                            else {
                                // This could be an active substrate grain
                                double R = dis(gen);
                                if (R < WallProb) {
                                
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
                                        CellType(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 'S';
                                        GrainID(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 0;
                                    }
                                    else {
                                        // At the interface
                                        CellType(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 'A';
                                        CritTimeStep(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 0;
                                        UndercoolingChange(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 0.1;
                                        SubstrateGrains++;
                                    }
                                }
                            }
                        }
                        else {
                            if (CritTimeStep(k*MyXSlices*MyYSlices + i*MyYSlices + j) > 0) {
                                CellType(k*MyXSlices*MyYSlices + i*MyYSlices + j) = 'D';
                                GrainID(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 0;
                            }
                            else {
                                CellType(k*MyXSlices*MyYSlices + i*MyYSlices + j) = 'S';
                            }
                        }
                    }
                    else {

                        CellType(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 'W';
                        //if (id == 0) cout << " W Cell " << i << " " << j << " " << k << endl;
                        GrainID(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 0;
                    }
                }
            }
        }
    }

    int TotalSubstrateGrains, TotalNucleatedGrains;
    MPI_Reduce(&SubstrateGrains,&TotalSubstrateGrains,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&NucleatedGrains,&TotalNucleatedGrains,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    if (id == 0) cout << "Number of substrate grains: " << TotalSubstrateGrains << endl;
    if (id == 0) cout << "Number of potential nucleated grains: " << TotalNucleatedGrains << endl;

    // Assign GrainIDs for substrate grains (positive values) and nucleated grains (negative values)
    // Grains for epitaxial growth
    int MyFirstSGrainID;
    if (id == 0) {
        int SBuf = SubstrateGrains+1;
        MPI_Send(&SBuf,1,MPI_INT,1,0,MPI_COMM_WORLD);
        MyFirstSGrainID = 1;
    }
    else if (id == np-1) {
        int RBuf;
        MPI_Recv(&RBuf,1,MPI_INT,np-2,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MyFirstSGrainID = RBuf;
        NextLayer_FirstSubstrateGrainID = SubstrateGrains + MyFirstSGrainID;
    }
    else {
        int RBuf;
        MPI_Recv(&RBuf,1,MPI_INT,id-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MyFirstSGrainID = RBuf;
        int SBuf = SubstrateGrains + MyFirstSGrainID;
        MPI_Send(&SBuf,1,MPI_INT,id+1,0,MPI_COMM_WORLD);
    }
    MPI_Bcast(&NextLayer_FirstSubstrateGrainID, 1, MPI_INT, np-1, MPI_COMM_WORLD);
    
    // Grains for nucleated growth
    int MyFirstNGrainID;
    if (id == 0) {
        int SBuf = -NucleatedGrains-1;
        MPI_Send(&SBuf,1,MPI_INT,1,0,MPI_COMM_WORLD);
        MyFirstNGrainID = -1;
    }
    else if (id == np-1) {
        int RBuf;
        MPI_Recv(&RBuf,1,MPI_INT,np-2,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MyFirstNGrainID = RBuf;
        NextLayer_FirstNucleatedGrainID = MyFirstNGrainID - NucleatedGrains;
    }
    else {
        int RBuf;
        MPI_Recv(&RBuf,1,MPI_INT,id-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MyFirstNGrainID = RBuf;
        int SBuf = MyFirstNGrainID - NucleatedGrains;
        MPI_Send(&SBuf,1,MPI_INT,id+1,0,MPI_COMM_WORLD);
    }
    MPI_Bcast(&NextLayer_FirstNucleatedGrainID, 1, MPI_INT, np-1, MPI_COMM_WORLD);

    ACount = 0;
    BCount = 0;
    CCount = 0;
    DCount = 0;
    ECount = 0;
    FCount = 0;
    GCount = 0;
    HCount = 0;

    // Assign Grain IDs to cells - epitaxial and nucleated grains
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(dTN,dTsigma);
    
    int GCounter = MyFirstSGrainID;
    int NCounter = MyFirstNGrainID;
    
    //if (id == 1) {
    for (int RankZ=1; RankZ<nz-1; RankZ++) {
        for (int RankX=1; RankX<MyXSlices-1; RankX++) {
            for (int RankY=1; RankY<MyYSlices-1; RankY++) {
                int D3D1ConvPosition = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                if (CellType(D3D1ConvPosition) == 'A') {
                    GrainID(D3D1ConvPosition) = GCounter;
                    GCounter++;
                    int GlobalX = RankX + MyXOffset;
                    int GlobalY = RankY + MyYOffset;
                    int MyGrainID = GrainID(D3D1ConvPosition);
                    NewGrain(MyXSlices,MyYSlices,nz,RankX,RankY,RankZ,MyGrainID,GlobalX,GlobalY,CellType,GrainID,DiagonalLength,DOCenter);
                    // The orientation for the new grain will depend on its Grain ID
                    int MyOrientation = GrainOrientation[((abs(MyGrainID) - 1) % NGrainOrientations)];
                    // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
                    // (xp,yp,zp) is the new cell's center on the global grid
                    double xp = GlobalX + 0.5;
                    double yp = GlobalY + 0.5;
                    double zp = RankZ + 0.5;
                    CritDiagLengthCalc(xp,yp,zp,MyOrientation,RankX,RankY,RankZ,D3D1ConvPosition,DOCenter(3*D3D1ConvPosition),DOCenter(3*D3D1ConvPosition+1),DOCenter(3*D3D1ConvPosition+2),NeighborX,NeighborY,NeighborZ, GrainUnitVector,TriangleIndex,CritDiagonalLength);

                    if (DecompositionStrategy == 1) {
                        if (RankY == 1) {
                            ACount++;
                            CellType(D3D1ConvPosition) = '1';
                        }
                        else if (RankY == MyYSlices-2) {
                            BCount++;
                            CellType(D3D1ConvPosition) = '1';
                        }
                    }
                    else {
                        if (RankY == 1) {
                            // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                            if (RankX == MyXSlices-2) {
                                ECount++;
                                CCount++;
                                ACount++;
                                CellType(D3D1ConvPosition) = '3';
                            }
                            else if (RankX == 1) {
                                GCount++;
                                DCount++;
                                ACount++;
                                CellType(D3D1ConvPosition) = '3';
                            }
                            else if ((RankX > 1)&&(RankX < MyXSlices-2)) {
                                // This is being sent to MyLeft
                                ACount++;
                                CellType(D3D1ConvPosition) = '1';
                            }
                        }
                        else if (RankY == MyYSlices-2) {
                            // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                            if (RankX == MyXSlices-2) {
                                FCount++;
                                CCount++;
                                BCount++;
                                CellType(D3D1ConvPosition) = '3';
                            }
                            else if (RankX == 1) {
                                HCount++;
                                DCount++;
                                BCount++;
                                CellType(D3D1ConvPosition) = '3';
                            }
                            else if ((RankX > 1)&&(RankX < MyXSlices-2)) {
                                BCount++;
                                CellType(D3D1ConvPosition) = '1';
                            }
                        }
                        else if ((RankX == 1)&&(RankY > 1)&&(RankY < MyYSlices-2)) {
                            DCount++;
                            CellType(D3D1ConvPosition) = '1';
                            //if (id == 0) cout << "RANK 0 LISTED " << MyNeighborX << " " << MyNeighborY << " " << MyNeighborZ << endl;
                        }
                        else if ((RankX == MyXSlices-2)&&(RankY > 1)&&(RankY < MyYSlices-2)) {
                            CCount++;
                            CellType(D3D1ConvPosition) = '1';
                        }
                    }
                }
                else if (CellType(D3D1ConvPosition) == 'N') {
                    // Each rank places the "possible" nuclei in random locations
                    NucLocI.push_back(RankX);
                    NucLocJ.push_back(RankY);
                    NucLocK.push_back(RankZ);
                    GrainID(RankZ*MyXSlices*MyYSlices+RankX*MyYSlices+RankY) = NCounter;
                    NCounter--;
                    // This site will nucleate at a certain time step if unclaimed by another grain
                    // Select local nucleation undercooling from Gaussian distribution
                    double LocNucUnd = distribution(generator);
                    NucleationUndercooling.push_back(LocNucUnd);
                    int CritNTimeStep = CritTimeStep(RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY) + round((LocNucUnd/UndercoolingChange(RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY)));
                    NucleationTimes.push_back(CritNTimeStep);
                }
            }
        }
    }
    
    // Remove delay cells not bordering others
    for (int RankZ=1; RankZ<nz-1; RankZ++) {
        for (int RankX=1; RankX<MyXSlices-1; RankX++) {
            for (int RankY=1; RankY<MyYSlices-1; RankY++) {
                int D3D1ConvPosition = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                if (CellType(D3D1ConvPosition) == 'D') {
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
                        if ((CellType(NeighborD3D1ConvPosition) == 'D')||(CellType(NeighborD3D1ConvPosition) == 'A')||(CellType(NeighborD3D1ConvPosition) == 'N')) {
                            LCount++;
                        }
                    }
                    if (LCount == 0) {
                        // This cell is returned to solid type
                        CellType(D3D1ConvPosition) = 'S';
                        GrainID(D3D1ConvPosition) = 0;
                    }
                }
            }
        }
    }

}

// After initializing grain structure and filling ghost nodes, calculate ghost nodes' grain data for nucleation events
void ConstrainedGrainsUpdate(int DecompositionStrategy, int MyXSlices, int MyYSlices, int id, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, vector <int> &NucLocI, vector <int> &NucLocJ, vector <int> &NucLocK, vector <int> &NucleationTimes, vector <float> &NucleationUndercooling, int* GrainOrientation, ViewC::HostMirror CellType, ViewI::HostMirror GrainID) {

    int ACount = 0;
    int BCount = 0;
    int CCount = 0;
    int DCount = 0;
    int ECount = 0;
    int FCount = 0;
    int GCount = 0;
    int HCount = 0;
    vector <float> ANuc, BNuc, CNuc, DNuc, ENuc, FNuc, GNuc, HNuc;
    

    for (int nnumber=0; nnumber<NucLocI.size(); nnumber++) {
        float LocGO = (float)(GrainID[NucLocK[nnumber]*MyXSlices*MyYSlices + NucLocI[nnumber]*MyYSlices + NucLocJ[nnumber]]);
        
        // Nucleation time, undercooling, location (X or Y and Z), orientation
        if (DecompositionStrategy == 1) {
            if (NucLocJ[nnumber] == 1) {
                ACount++;
                ANuc.push_back(NucleationTimes[nnumber]);
                ANuc.push_back(NucleationUndercooling[nnumber]);
                ANuc.push_back(NucLocI[nnumber]);
                ANuc.push_back(NucLocK[nnumber]);
                ANuc.push_back(LocGO);
            }
            else if (NucLocJ[nnumber] == MyYSlices-2) {
                BCount++;
                BNuc.push_back(NucleationTimes[nnumber]);
                BNuc.push_back(NucleationUndercooling[nnumber]);
                BNuc.push_back(NucLocI[nnumber]);
                BNuc.push_back(NucLocK[nnumber]);
                BNuc.push_back(LocGO);
            }
        }
        else {
            if (NucLocJ[nnumber] == 1) {
                // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                if (NucLocI[nnumber] == MyXSlices-2) {
                    ECount++;
                    ENuc.push_back(NucleationTimes[nnumber]);
                    ENuc.push_back(NucleationUndercooling[nnumber]);
                    ENuc.push_back(NucLocK[nnumber]);
                    ENuc.push_back(LocGO);
        
                    CCount++;
                    CNuc.push_back(NucleationTimes[nnumber]);
                    CNuc.push_back(NucleationUndercooling[nnumber]);
                    CNuc.push_back(NucLocJ[nnumber]);
                    CNuc.push_back(NucLocK[nnumber]);
                    CNuc.push_back(LocGO);
                    
                    ACount++;
                    ANuc.push_back(NucleationTimes[nnumber]);
                    ANuc.push_back(NucleationUndercooling[nnumber]);
                    ANuc.push_back(NucLocI[nnumber]);
                    ANuc.push_back(NucLocK[nnumber]);
                    ANuc.push_back(LocGO);
                }
                else if (NucLocI[nnumber] == 1) {
                    GCount++;
                    GNuc.push_back(NucleationTimes[nnumber]);
                    GNuc.push_back(NucleationUndercooling[nnumber]);
                    GNuc.push_back(NucLocK[nnumber]);
                    GNuc.push_back(LocGO);
                    
                    DCount++;
                    DNuc.push_back(NucleationTimes[nnumber]);
                    DNuc.push_back(NucleationUndercooling[nnumber]);
                    DNuc.push_back(NucLocJ[nnumber]);
                    DNuc.push_back(NucLocK[nnumber]);
                    DNuc.push_back(LocGO);
                    
                    ACount++;
                    ANuc.push_back(NucleationTimes[nnumber]);
                    ANuc.push_back(NucleationUndercooling[nnumber]);
                    ANuc.push_back(NucLocI[nnumber]);
                    ANuc.push_back(NucLocK[nnumber]);
                    ANuc.push_back(LocGO);
                }
                else if ((NucLocI[nnumber] > 1)&&(NucLocI[nnumber] < MyXSlices-2)) {
                    // This is being sent to MyLeft
                    ACount++;
                    ANuc.push_back(NucleationTimes[nnumber]);
                    ANuc.push_back(NucleationUndercooling[nnumber]);
                    ANuc.push_back(NucLocI[nnumber]);
                    ANuc.push_back(NucLocK[nnumber]);
                    ANuc.push_back(LocGO);
                }
            }
            else if (NucLocJ[nnumber] == MyYSlices-2) {
                // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                if (NucLocI[nnumber] == MyXSlices-2) {
                    FCount++;
                    FNuc.push_back(NucleationTimes[nnumber]);
                    FNuc.push_back(NucleationUndercooling[nnumber]);
                    FNuc.push_back(NucLocK[nnumber]);
                    FNuc.push_back(LocGO);
                    
                    CCount++;
                    CNuc.push_back(NucleationTimes[nnumber]);
                    CNuc.push_back(NucleationUndercooling[nnumber]);
                    CNuc.push_back(NucLocJ[nnumber]);
                    CNuc.push_back(NucLocK[nnumber]);
                    CNuc.push_back(LocGO);
                    
                    BCount++;
                    BNuc.push_back(NucleationTimes[nnumber]);
                    BNuc.push_back(NucleationUndercooling[nnumber]);
                    BNuc.push_back(NucLocI[nnumber]);
                    BNuc.push_back(NucLocK[nnumber]);
                    BNuc.push_back(LocGO);
                }
                else if (NucLocI[nnumber] == 1) {
                    HCount++;
                    HNuc.push_back(NucleationTimes[nnumber]);
                    HNuc.push_back(NucleationUndercooling[nnumber]);
                    HNuc.push_back(NucLocK[nnumber]);
                    HNuc.push_back(LocGO);
                    
                    DCount++;
                    DNuc.push_back(NucleationTimes[nnumber]);
                    DNuc.push_back(NucleationUndercooling[nnumber]);
                    DNuc.push_back(NucLocJ[nnumber]);
                    DNuc.push_back(NucLocK[nnumber]);
                    DNuc.push_back(LocGO);
                    BCount++;
                    BNuc.push_back(NucleationTimes[nnumber]);
                    BNuc.push_back(NucleationUndercooling[nnumber]);
                    BNuc.push_back(NucLocI[nnumber]);
                    BNuc.push_back(NucLocK[nnumber]);
                    BNuc.push_back(LocGO);
                }
                else if ((NucLocI[nnumber] > 1)&&(NucLocI[nnumber] < MyXSlices-2)) {
                    BCount++;
                    BNuc.push_back(NucleationTimes[nnumber]);
                    BNuc.push_back(NucleationUndercooling[nnumber]);
                    BNuc.push_back(NucLocI[nnumber]);
                    BNuc.push_back(NucLocK[nnumber]);
                    BNuc.push_back(LocGO);
                }
            }
            else if ((NucLocI[nnumber] == 1)&&(NucLocJ[nnumber] > 1)&&(NucLocJ[nnumber] < MyYSlices-2)) {
                DCount++;
                DNuc.push_back(NucleationTimes[nnumber]);
                DNuc.push_back(NucleationUndercooling[nnumber]);
                DNuc.push_back(NucLocJ[nnumber]);
                DNuc.push_back(NucLocK[nnumber]);
                DNuc.push_back(LocGO);
            }
            else if ((NucLocI[nnumber] == MyXSlices-2)&&(NucLocJ[nnumber] > 1)&&(NucLocJ[nnumber] < MyYSlices-2)) {
                CCount++;
                CNuc.push_back(NucleationTimes[nnumber]);
                CNuc.push_back(NucleationUndercooling[nnumber]);
                CNuc.push_back(NucLocJ[nnumber]);
                CNuc.push_back(NucLocK[nnumber]);
                CNuc.push_back(LocGO);
            }
        }
    }
    
   // cout << "ER" << endl;
        
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
    
    
    // Buffers for recieving ghost node data
    float GhostNodesAR[5*ARCount];
    float GhostNodesBR[5*BRCount];
    float GhostNodesCR[5*CRCount];
    float GhostNodesDR[5*DRCount];
    float GhostNodesER[4*ERCount];
    float GhostNodesFR[4*FRCount];
    float GhostNodesGR[4*GRCount];
    float GhostNodesHR[4*HRCount];
    
    
    // Collect ghost node data and send to other ranks- left and right
    if (ACount > 0) {
        float GhostNodesA[5*ACount];
        for (int i=0; i<5*ACount; i++) {
            GhostNodesA[i] = ANuc[i];
        }
        if (BRCount == 0) {
            // Sending data to id = id - 1 only
            MPI_Send(&GhostNodesA,ACount*5,MPI_FLOAT,MyLeft,0,MPI_COMM_WORLD);
        }
        else {
            // Sending data to id = id - 1 and recieving data from id = id + 1
            MPI_Sendrecv(&GhostNodesA,ACount*5,MPI_FLOAT,MyLeft,0,&GhostNodesBR,BRCount*7,MPI_FLOAT,MyRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        // cout << "ID = " << id << " ACount = " << ACount << " ABuf = " << BufACount << endl;
    }
    else if (BRCount > 0) {
        // Recieving data from id = id + 1 only
        MPI_Recv(&GhostNodesBR,BRCount*5,MPI_FLOAT,MyRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    
    if (BCount > 0) {
        float GhostNodesB[5*BCount];
        for (int i=0; i<5*BCount; i++) {
            GhostNodesB[i] = BNuc[i];
        }
        if (ARCount == 0) {
            // Sending data to id = id + 1 only
            MPI_Send(&GhostNodesB,BCount*5,MPI_FLOAT,MyRight,1,MPI_COMM_WORLD);
        }
        else {
            // Sending data to id = id + 1 and recieving data from id = id - 1
            MPI_Sendrecv(&GhostNodesB,BCount*5,MPI_FLOAT,MyRight,1,&GhostNodesAR,ARCount*7,MPI_FLOAT,MyLeft,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    }
    else if (ARCount > 0) {
        // Recieving data from id = id - 1 only
        MPI_Recv(&GhostNodesAR,ARCount*5,MPI_FLOAT,MyLeft,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    
    if (DecompositionStrategy != 1) {
        // Collect ghost node data and send to other ranks- in and out
        if (CCount > 0) {
            float GhostNodesC[5*CCount];
            for (int i=0; i<5*CCount; i++) {
                GhostNodesC[i] = CNuc[i];
            }
            
            if (DRCount == 0) {
                // Sending data only
                MPI_Send(&GhostNodesC,CCount*5,MPI_FLOAT,MyIn,0,MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
                MPI_Sendrecv(&GhostNodesC,CCount*5,MPI_FLOAT,MyIn,0,&GhostNodesDR,DRCount*7,MPI_FLOAT,MyOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        else if (DRCount > 0) {
            // Recieving data only
            MPI_Recv(&GhostNodesDR,DRCount*5,MPI_FLOAT,MyOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        
        if (DCount > 0) {
            float GhostNodesD[5*DCount];
            for (int i=0; i<5*DCount; i++) {
                GhostNodesD[i] = DNuc[i];
            }
            if (CRCount == 0) {
                // Sending data only
                MPI_Send(&GhostNodesD,DCount*5,MPI_FLOAT,MyOut,1,MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
                MPI_Sendrecv(&GhostNodesD,DCount*5,MPI_FLOAT,MyOut,1,&GhostNodesCR,CRCount*7,MPI_FLOAT,MyIn,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            //    cout << "ID = " << id << " DCount = " << DCount << " DBuf = " << BufDCount << endl;
        }
        else if (CRCount > 0) {
            // Recieving data only
            MPI_Recv(&GhostNodesCR,CRCount*5,MPI_FLOAT,MyIn,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        
        // Collect ghost node data and send to other ranks- MyLeftIn and MyRightOut
        if (ECount > 0) {
            float GhostNodesE[4*ECount];
            for (int i=0; i<4*ECount; i++) {
                GhostNodesE[i] = ENuc[i];
            }
            if (HRCount == 0) {
                // Sending data only
                MPI_Send(&GhostNodesE,ECount*4,MPI_FLOAT,MyLeftIn,0,MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
                MPI_Sendrecv(&GhostNodesE,ECount*4,MPI_FLOAT,MyLeftIn,0,&GhostNodesHR,HRCount*6,MPI_FLOAT,MyRightOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            //  cout << "ID = " << id << " ECount = " << ECount << " EBuf = " << BufECount << endl;
        }
        else if (HRCount > 0) {
            // Recieving data only
            MPI_Recv(&GhostNodesHR,HRCount*4,MPI_FLOAT,MyRightOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        
        if (HCount > 0) {
            float GhostNodesH[4*HCount];
            for (int i=0; i<4*HCount; i++) {
                GhostNodesH[i] = HNuc[i];
            }
            if (ERCount == 0) {
                // Sending data only
                MPI_Send(&GhostNodesH,HCount*4,MPI_FLOAT,MyRightOut,0,MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
                MPI_Sendrecv(&GhostNodesH,HCount*4,MPI_FLOAT,MyRightOut,0,&GhostNodesER,ERCount*6,MPI_FLOAT,MyLeftIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        else if (ERCount > 0) {
            // Recieving data only
            MPI_Recv(&GhostNodesER,ERCount*4,MPI_FLOAT,MyLeftIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        
        // Collect ghost node data and send to other ranks- MyRightIn and MyLeftOut
        if (FCount > 0) {
            float GhostNodesF[4*FCount];
            for (int i=0; i<4*FCount; i++) {
                GhostNodesF[i] = FNuc[i];
            }
            if (GRCount == 0) {
                // Sending data only
                MPI_Send(&GhostNodesF,FCount*4,MPI_FLOAT,MyRightIn,1,MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
                MPI_Sendrecv(&GhostNodesF,FCount*4,MPI_FLOAT,MyRightIn,1,&GhostNodesGR,GRCount*6,MPI_FLOAT,MyLeftOut,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        else if (GRCount > 0) {
            // Recieving data only
            MPI_Recv(&GhostNodesGR,GRCount*4,MPI_FLOAT,MyLeftOut,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        
        if (GCount > 0) {
            float GhostNodesG[4*GCount];
            for (int i=0; i<4*GCount; i++) {
                GhostNodesG[i] = GNuc[i];
            }
            if (FRCount == 0) {
                // Sending data only
                MPI_Send(&GhostNodesG,GCount*4,MPI_FLOAT,MyLeftOut,1,MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
                MPI_Sendrecv(&GhostNodesG,GCount*4,MPI_FLOAT,MyLeftOut,1,&GhostNodesFR,FRCount*6,MPI_FLOAT,MyRightIn,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        else if (FRCount > 0) {
            // Recieving data only
            MPI_Recv(&GhostNodesFR,FRCount*4,MPI_FLOAT,MyRightIn,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    }
    

    // Place ghost node data recieved from the left (if needed)
    if (ARCount > 0) {
        for (int i=0; i<ARCount; i++) {
            NucleationTimes.push_back(GhostNodesAR[5*i]);
            NucleationUndercooling.push_back(GhostNodesAR[5*i+1]);
            int RankX = (int)(GhostNodesAR[5*i+2]);
            int RankY = 0;
            int RankZ = (int)(GhostNodesAR[5*i+3]);
            NucLocI.push_back(RankX);
            NucLocJ.push_back(RankY);
            NucLocK.push_back(RankZ);
            int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
            CellType(CellLocation) = 'N';
            GrainID(CellLocation) = (int)(GhostNodesAR[5*i+4]);
        }
    }

    // Place ghost node data recieved from the right (if needed)
    if (BRCount > 0) {
        for (int i=0; i<BRCount; i++) {
            NucleationTimes.push_back(GhostNodesBR[5*i]);
            NucleationUndercooling.push_back(GhostNodesBR[5*i+1]);
            int RankX = (int)(GhostNodesBR[5*i+2]);
            int RankY = MyYSlices-1;
            int RankZ = (int)(GhostNodesBR[5*i+3]);
            NucLocI.push_back(RankX);
            NucLocJ.push_back(RankY);
            NucLocK.push_back(RankZ);
            int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
            CellType(CellLocation) = 'N';
            GrainID(CellLocation) = (int)(GhostNodesBR[5*i+4]);
        }
    }

    if (DecompositionStrategy != 1) {
        // Place ghost node data recieved from in plane (if needed)
        if (CRCount > 0) {
            for (int i=0; i<CRCount; i++) {
                NucleationTimes.push_back(GhostNodesCR[5*i]);
                NucleationUndercooling.push_back(GhostNodesCR[5*i+1]);
                int RankX = MyXSlices-1;
                int RankY = (int)(GhostNodesCR[5*i+2]);
                int RankZ = (int)(GhostNodesCR[5*i+3]);
                NucLocI.push_back(RankX);
                NucLocJ.push_back(RankY);
                NucLocK.push_back(RankZ);
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                CellType(CellLocation) = 'N';
                GrainID(CellLocation) = (int)(GhostNodesCR[5*i+4]);
            }
        }

        // Place ghost node data recieved from out of plane (if needed)
        if (DRCount > 0) {
            for (int i=0; i<DRCount; i++) {
                NucleationTimes.push_back(GhostNodesDR[5*i]);
                NucleationUndercooling.push_back(GhostNodesDR[5*i+1]);
                int RankX = 0;
                int RankY = (int)(GhostNodesDR[5*i+2]);
                int RankZ = (int)(GhostNodesDR[5*i+3]);
                NucLocI.push_back(RankX);
                NucLocJ.push_back(RankY);
                NucLocK.push_back(RankZ);
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                CellType(CellLocation) = 'N';
                GrainID(CellLocation) = (int)(GhostNodesDR[5*i+4]);
            }
        }

        // Place ghost node data recieved from left and out of plane (if needed)
        if (ERCount > 0) {
            for (int i=0; i<ERCount; i++) {
                NucleationTimes.push_back(GhostNodesER[4*i]);
                NucleationUndercooling.push_back(GhostNodesER[4*i+1]);
                int RankX = MyXSlices-1;
                int RankY = 0;
                int RankZ = (int)(GhostNodesER[4*i+2]);
                NucLocI.push_back(RankX);
                NucLocJ.push_back(RankY);
                NucLocK.push_back(RankZ);
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                CellType(CellLocation) = 'N';
                GrainID(CellLocation) = (int)(GhostNodesER[4*i+3]);
            }
        }

        // Place ghost node data recieved from right and out of plane (if needed)
        if (FRCount > 0) {
            for (int i=0; i<FRCount; i++) {
                NucleationTimes.push_back(GhostNodesFR[4*i]);
                NucleationUndercooling.push_back(GhostNodesFR[4*i+1]);
                int RankX = MyXSlices-1;
                int RankY = MyYSlices-1;
                int RankZ = (int)(GhostNodesFR[4*i+2]);
                NucLocI.push_back(RankX);
                NucLocJ.push_back(RankY);
                NucLocK.push_back(RankZ);
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                CellType(CellLocation) = 'N';
                GrainID(CellLocation) = (int)(GhostNodesFR[4*i+3]);
            }
        }

        // Place ghost node data recieved from left and into plane (if needed)
        if (GRCount > 0) {
            for (int i=0; i<GRCount; i++) {
                NucleationTimes.push_back(GhostNodesGR[4*i]);
                NucleationUndercooling.push_back(GhostNodesGR[4*i+1]);
                int RankX = 0;
                int RankY = 0;
                int RankZ = (int)(GhostNodesGR[4*i+2]);
                NucLocI.push_back(RankX);
                NucLocJ.push_back(RankY);
                NucLocK.push_back(RankZ);
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                CellType(CellLocation) = 'N';
                GrainID(CellLocation) = (int)(GhostNodesGR[4*i+3]);
            }
        }

        if (HRCount > 0) {
            for (int i=0; i<HRCount; i++) {
                NucleationTimes.push_back(GhostNodesHR[4*i]);
                NucleationUndercooling.push_back(GhostNodesHR[4*i+1]);
                int RankX = 0;
                int RankY = MyYSlices-1;
                int RankZ = (int)(GhostNodesHR[4*i+2]);
                NucLocI.push_back(RankX);
                NucLocJ.push_back(RankY);
                NucLocK.push_back(RankZ);
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                CellType(CellLocation) = 'N';
                GrainID(CellLocation) = (int)(GhostNodesHR[4*i+3]);
            }
        }
    }

}






void LayerSetup(int layernumber, int LayerHeight, int InitialGrainWidth, int NGrainOrientations, int DecompositionStrategy, int ProcessorsInXDirection, int ProcessorsInYDirection, int nx, int ny, int nz, string TemperatureDataType, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int id, int np, int ierr, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyLeftOut, int MyRightIn, int MyRightOut, int ItList[9][26], int NeighborX[26], int NeighborY[26], int NeighborZ[26], vector <int> &NucLocI, vector <int> &NucLocJ, vector <int> &NucLocK, vector <int> &NucleationTimes, vector <float> &NucleationUndercooling, int* GrainOrientation, float* GrainUnitVector, ViewF::HostMirror DiagonalLength, ViewC::HostMirror CellType, ViewI::HostMirror TriangleIndex, ViewI::HostMirror GrainID, int* GrainID_Stored, ViewF::HostMirror CritDiagonalLength, ViewF::HostMirror DOCenter, ViewI::HostMirror CritTimeStep, ViewF::HostMirror UndercoolingChange, ViewF::HostMirror UndercoolingCurrent, string tempfile, double deltax, double deltat, double NMax, double dTN, double dTsigma, string FileA1, string FileA2, int XMin, int XMax, int YMin, int YMax, int ZMin, int ZMax, int &NextLayer_FirstSubstrateGrainID, int &NextLayer_FirstNucleatedGrainID, int &ACount, int &BCount, int &CCount, int &DCount, int &ECount, int &FCount, int &GCount, int &HCount) {
    
    std::seed_seq seed = {765,4,1111};
    std::array<unsigned,5> sequence;
    seed.generate(sequence.begin(),sequence.end());
    mt19937_64 gen(seed);
    
    uniform_real_distribution<double> dis(0.0, 1.0);
    srand((unsigned)id);
    // Convert initial grain spacing to a grain density
    double NWall = 1.0/(InitialGrainWidth*InitialGrainWidth*InitialGrainWidth*pow(10,-18));
    double WallProb = NWall*deltax*deltax*deltax;
    double BulkProb = NMax*deltax*deltax*deltax;
    int SubstrateGrains = 0;
    int NucleatedGrains = 0;
    
    int LayerOffset = layernumber*LayerHeight*MyXSlices*MyYSlices;
    
    for (int k=1; k<=LayerHeight; k++)  {
        for(int i=0; i<MyXSlices; i++) {
            for(int j=0; j<MyYSlices; j++) {
                GrainID_Stored[LayerOffset+(k-1)*MyXSlices*MyYSlices+i*MyYSlices+j] = GrainID(k*MyXSlices*MyYSlices+i*MyYSlices+j);
            }
        }
    }
    
    int* GrainID_Temp = new int[nz*MyXSlices*MyYSlices];
    for (int i=0; i<nz*MyXSlices*MyYSlices; i++) {
        GrainID_Temp[i] = GrainID(i);
    }

    // Shift substrate grain IDs for new layer
    for (int k=1; k<=nz-LayerHeight-2; k++)  {
        for(int i=0; i<MyXSlices; i++) {
            for(int j=0; j<MyYSlices; j++) {
                GrainID(k*MyXSlices*MyYSlices+i*MyYSlices+j) = GrainID_Temp[(k+LayerHeight)*MyXSlices*MyYSlices+i*MyYSlices+j];
            }
        }
    }
    // No new grains in the deposited material - only epitaxial growth from previous layer
    for (int k=nz-LayerHeight-1; k<nz-1; k++)  {
        for(int i=0; i<MyXSlices; i++) {
            for(int j=0; j<MyYSlices; j++) {
                GrainID(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 0;
                UndercoolingCurrent(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 0;
            }
        }
    }

    delete [] GrainID_Temp;
    
    // Initialize cell types - epitaxial sites can be anywhere at the interface, bulk sites in melt pool area
    for (int k=0; k<nz; k++)  {
        for(int i=0; i<MyXSlices; i++) {
            for(int j=0; j<MyYSlices; j++) {
                int GlobalX = i + MyXOffset;
                int GlobalY = j + MyYOffset;
                // Non-wall cells
                if ((GlobalX > -1)&&(GlobalX < nx)&&(GlobalY > -1)&&(GlobalY < ny)&&(k > 0)&&(k < nz-1)) {
                    
                    if ((i > 0)&&(i < MyXSlices-1)&&(j > 0)&&(j < MyYSlices-1)) {
                        // Cells that are either 'D' or 'N'
                        if (CritTimeStep(k*MyXSlices*MyYSlices + i*MyYSlices + j) > 0) {
                            double R = dis(gen);
                            if (R < BulkProb) {
                                NucleatedGrains++;
                                CellType(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 'N';
                            }
                            else {
                                CellType(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 'D';
                            }
                        }
                        else {
                            // This could be an active substrate grain
                            double R = dis(gen);
                            if ((R < BulkProb)||(k < nz-LayerHeight)) {
                                
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
                                    CellType(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 'S';
                                }
                                else {
                                    if (k >= nz-LayerHeight) {
                                        // At the interface
                                        CellType(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 'A';
                                        CritTimeStep(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 0;
                                        UndercoolingChange(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 0.1;
                                    }
                                    else {
                                        double R2 = dis(gen);
                                        if (R2 < WallProb) {
                                            SubstrateGrains++;
                                            CellType(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 'A';
                                            CritTimeStep(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 0;
                                            UndercoolingChange(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 0.1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    else {
                        if (CritTimeStep(k*MyXSlices*MyYSlices + i*MyYSlices + j) > 0) {
                            CellType(k*MyXSlices*MyYSlices + i*MyYSlices + j) = 'D';
                            GrainID(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 0;
                        }
                        else {
                            CellType(k*MyXSlices*MyYSlices + i*MyYSlices + j) = 'S';
                        }
                    }
                }
                else {
                    CellType(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 'W';
                    GrainID(k*MyXSlices*MyYSlices+i*MyYSlices+j) = 0;
                }
            }
        }
    }
    
    int ThisLayer_FirstSubstrateGrainID = NextLayer_FirstSubstrateGrainID;
    int ThisLayer_FirstNucleatedGrainID = NextLayer_FirstNucleatedGrainID;
    
    int TotalSubstrateGrains, TotalNucleatedGrains;
    MPI_Reduce(&SubstrateGrains,&TotalSubstrateGrains,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&NucleatedGrains,&TotalNucleatedGrains,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    if (id == 0) cout << "Number of substrate grains: " << TotalSubstrateGrains << endl;
    if (id == 0) cout << "Number of potential nucleated grains: " << TotalNucleatedGrains << endl;
    
    // Assign GrainIDs for substrate grains (positive values) and nucleated grains (negative values)
    // Grains for epitaxial growth
    int MyFirstSGrainID;
    if (id == 0) {
        int SBuf = SubstrateGrains+ThisLayer_FirstSubstrateGrainID;
        MPI_Send(&SBuf,1,MPI_INT,1,0,MPI_COMM_WORLD);
        MyFirstSGrainID = ThisLayer_FirstSubstrateGrainID;
    }
    else if (id == np-1) {
        int RBuf;
        MPI_Recv(&RBuf,1,MPI_INT,np-2,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MyFirstSGrainID = RBuf;
        NextLayer_FirstSubstrateGrainID = SubstrateGrains + MyFirstSGrainID;
    }
    else {
        int RBuf;
        MPI_Recv(&RBuf,1,MPI_INT,id-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MyFirstSGrainID = RBuf;
        int SBuf = SubstrateGrains + MyFirstSGrainID;
        MPI_Send(&SBuf,1,MPI_INT,id+1,0,MPI_COMM_WORLD);
    }
    MPI_Bcast(&NextLayer_FirstSubstrateGrainID, 1, MPI_INT, np-1, MPI_COMM_WORLD);
    
    // Grains for nucleated growth
    int MyFirstNGrainID;
    if (id == 0) {
        int SBuf = -NucleatedGrains-ThisLayer_FirstNucleatedGrainID;
        MPI_Send(&SBuf,1,MPI_INT,1,0,MPI_COMM_WORLD);
        MyFirstNGrainID = ThisLayer_FirstNucleatedGrainID;
    }
    else if (id == np-1) {
        int RBuf;
        MPI_Recv(&RBuf,1,MPI_INT,np-2,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MyFirstNGrainID = RBuf;
        NextLayer_FirstNucleatedGrainID = MyFirstNGrainID - NucleatedGrains;
    }
    else {
        int RBuf;
        MPI_Recv(&RBuf,1,MPI_INT,id-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MyFirstNGrainID = RBuf;
        int SBuf = MyFirstNGrainID - NucleatedGrains;
        MPI_Send(&SBuf,1,MPI_INT,id+1,0,MPI_COMM_WORLD);
    }
    MPI_Bcast(&NextLayer_FirstNucleatedGrainID, 1, MPI_INT, np-1, MPI_COMM_WORLD);

    ACount = 0;
    BCount = 0;
    CCount = 0;
    DCount = 0;
    ECount = 0;
    FCount = 0;
    GCount = 0;
    HCount = 0;
    
    // Assign Grain IDs to cells - epitaxial and nucleated grains
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(dTN,dTsigma);
    
    int GCounter = MyFirstSGrainID;
    int NCounter = MyFirstNGrainID;
    
    for (int RankZ=1; RankZ<nz-1; RankZ++) {
        for (int RankX=1; RankX<MyXSlices-1; RankX++) {
            for (int RankY=1; RankY<MyYSlices-1; RankY++) {
                int D3D1ConvPosition = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                if (CellType(D3D1ConvPosition) == 'A') {
                    // Grains in new deposited substrate
                    if (RankZ >= nz-LayerHeight) {
                           GrainID(D3D1ConvPosition) = GCounter;
                           GCounter++;
                    }
                    int GlobalX = RankX + MyXOffset;
                    int GlobalY = RankY + MyYOffset;
                    int MyGrainID = GrainID(D3D1ConvPosition);
                    NewGrain(MyXSlices,MyYSlices,nz,RankX,RankY,RankZ,MyGrainID,GlobalX,GlobalY,CellType,GrainID,DiagonalLength,DOCenter);
                    // The orientation for the new grain will depend on its Grain ID
                    int MyOrientation = GrainOrientation[((abs(MyGrainID) - 1) % NGrainOrientations)];
                    // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
                    // (xp,yp,zp) is the new cell's center on the global grid
                    double xp = GlobalX + 0.5;
                    double yp = GlobalY + 0.5;
                    double zp = RankZ + 0.5;
                    CritDiagLengthCalc(xp,yp,zp,MyOrientation,RankX,RankY,RankZ,D3D1ConvPosition, DOCenter(3*D3D1ConvPosition),DOCenter(3*D3D1ConvPosition+1),DOCenter(3*D3D1ConvPosition+2),NeighborX,NeighborY,NeighborZ, GrainUnitVector,TriangleIndex,CritDiagonalLength);
                    if (DecompositionStrategy == 1) {
                        if (RankY == 1) {
                            ACount++;
                            CellType(D3D1ConvPosition) = '1';
                        }
                        else if (RankY == MyYSlices-2) {
                            BCount++;
                            CellType(D3D1ConvPosition) = '1';
                        }
                    }
                    else {
                        if (RankY == 1) {
                            // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                            if (RankX == MyXSlices-2) {
                                ECount++;
                                CCount++;
                                ACount++;
                                CellType(D3D1ConvPosition) = '3';
                            }
                            else if (RankX == 1) {
                                GCount++;
                                DCount++;
                                ACount++;
                                CellType(D3D1ConvPosition) = '3';
                            }
                            else if ((RankX > 1)&&(RankX < MyXSlices-2)) {
                                // This is being sent to MyLeft
                                ACount++;
                                CellType(D3D1ConvPosition) = '1';
                            }
                        }
                        else if (RankY == MyYSlices-2) {
                            // This is also potentially being sent to MyLeftIn/MyLeftOut/MyIn/MyOut
                            if (RankX == MyXSlices-2) {
                                FCount++;
                                CCount++;
                                BCount++;
                                CellType(D3D1ConvPosition) = '3';
                            }
                            else if (RankX == 1) {
                                HCount++;
                                DCount++;
                                BCount++;
                                CellType(D3D1ConvPosition) = '3';
                            }
                            else if ((RankX > 1)&&(RankX < MyXSlices-2)) {
                                BCount++;
                                CellType(D3D1ConvPosition) = '1';
                            }
                        }
                        else if ((RankX == 1)&&(RankY > 1)&&(RankY < MyYSlices-2)) {
                            DCount++;
                            CellType(D3D1ConvPosition) = '1';
                            //if (id == 0) cout << "RANK 0 LISTED " << MyNeighborX << " " << MyNeighborY << " " << MyNeighborZ << endl;
                        }
                        else if ((RankX == MyXSlices-2)&&(RankY > 1)&&(RankY < MyYSlices-2)) {
                            CCount++;
                            CellType(D3D1ConvPosition) = '1';
                        }
                    }
                }
                else if (CellType[D3D1ConvPosition] == 'N') {
                    // Each rank places the "possible" nuclei in random locations
                    NucLocI.push_back(RankX);
                    NucLocJ.push_back(RankY);
                    NucLocK.push_back(RankZ);
                    GrainID(RankZ*MyXSlices*MyYSlices+RankX*MyYSlices+RankY) = NCounter;
                    NCounter--;
                    // This site will nucleate at a certain time step if unclaimed by another grain
                    // Select local nucleation undercooling from Gaussian distribution
                    double LocNucUnd = distribution(generator);
                    NucleationUndercooling.push_back(LocNucUnd);
                    int CritNTimeStep = CritTimeStep(RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY) + round((LocNucUnd/UndercoolingChange(RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY)));
                    NucleationTimes.push_back(CritNTimeStep);
                }
            }
        }
    }

    // Remove delay cells not bordering others
    for (int RankZ=1; RankZ<nz-1; RankZ++) {
        for (int RankX=1; RankX<MyXSlices-1; RankX++) {
            for (int RankY=1; RankY<MyYSlices-1; RankY++) {
                int D3D1ConvPosition = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                if (CellType(D3D1ConvPosition) == 'D') {
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
                        if ((CellType(NeighborD3D1ConvPosition) == 'D')||(CellType(NeighborD3D1ConvPosition) == 'A')||(CellType(NeighborD3D1ConvPosition) == 'N')) {
                            LCount++;
                        }
                    }
                    if (LCount == 0) {
                        // This cell is returned to solid type
                        CellType(D3D1ConvPosition) = 'S';
                        GrainID(D3D1ConvPosition) = 0;
                    }
                }
            }
        }
    }

                                
}



