#include "header.h"

#include <iostream>
#include <string>
#include <regex>

using namespace std;
// Initializes input parameters, mesh, temperature field, and grain structures for CA simulations

// Skip initial lines in input files.
void skipLines( std::ifstream &stream )
{
    std::string line;
    while(getline(stream, line))
    {
        if (line == "*****")
            break;
    }
}

// Verify the required input was included with the correct format.
std::string parseInput( std::ifstream &stream, std::string key )
{
    std::string line;
    std::getline(stream, line);
    std::size_t colon = line.find(":");
    std::string actual_key = line.substr(0, colon);

    // Check for keyword
    if ( actual_key.find( key ) == std::string::npos )
    {
        string error = "Required input \"" + key + "\" not found in input file.";
        throw std::runtime_error( error );
    }
    // Check for colon seperator
    if ( colon == std::string::npos )
    {
        string error = "Input \"" + key + "\" must be separated from value by \":\"." ;
        throw std::runtime_error( error );
    }

    // Remove whitespace
    std::string val = line.substr(colon+1,string::npos);
    std::regex r("\\s+");
    val = std::regex_replace(val, r, "");

    return val;
}

// Verify the required boolean input was included with the correct format.
bool parseInputBool( std::ifstream &stream, std::string key )
{
    std::string val = parseInput( stream, key );
    if ( val == "N" )
    {
        return false;
    }
    else if ( val == "Y" )
    {
        return true;
    }
    else
    {
        string error = "Input \"" + key + "\" must be \"Y\" or \"N\".";
        throw std::runtime_error( error );
    }
}

// Read ExaCA input file.
void InputReadFromFile(int id, string InputFile, string &SimulationType, int &DecompositionStrategy, double &AConst, double &BConst, double &CConst, double &DConst, double& FreezingRange, double &deltax, double &NMax, double &dTN, double &dTsigma, string &OutputFile, string &GrainOrientationFile, string &tempfile, int &TempFilesInSeries, bool& TruchasMultilayer, string &ExtraWalls, double &HT_deltax, string &TemperatureDataSource, double &deltat, int &NumberOfLayers, int &LayerHeight, string &SubstrateFileName, double &G, double &R, int &nx, int &ny, int &nz, double &FractSurfaceSitesActive, string &PathToOutput, int &NumberOfTruchasRanks, bool (&FilesToPrint)[6], bool &PrintFilesYN) {

    size_t backslash = InputFile.find_last_of("/");
    string FilePath = InputFile.substr(0, backslash);

    ifstream InputData, MaterialData;
    string Colon = ":";
    string Quote = "'";
    InputData.open(InputFile);
    if (id == 0) cout << "Input file " << InputFile << " opened" << endl;
    skipLines(InputData);

    std::string val;

    // Simulation Type
    SimulationType = parseInput(InputData, "Problem type");

    // Decomposition strategy
    val = parseInput(InputData, "Decomposition strategy");
    DecompositionStrategy = stoi(val,nullptr,10);
    
    // Material (opening a separate file to obtain values for A, B, C, and D for the interfacial reponse function)
    std::string MaterialFile = parseInput(InputData, "Material");
    MaterialData.open(FilePath + "/Materials/" + MaterialFile);
    skipLines(MaterialData);

    // Interfacial response function A
    val = parseInput(MaterialData, "A");
    AConst = atof(val.c_str());
    // Interfacial response function B
    val = parseInput(MaterialData, "B");
    BConst = atof(val.c_str());
    // Interfacial response function C
    val = parseInput(MaterialData, "C");
    CConst = atof(val.c_str());
    // Interfacial response function D
    val = parseInput(MaterialData, "D");
    DConst = atof(val.c_str());
    // Alloy freezing range
    val = parseInput(MaterialData, "Alloy freezing range");
    FreezingRange = atof(val.c_str());

    MaterialData.close();

    // CA cell size
    val = parseInput(InputData, "Cell size");
    deltax = atof(val.c_str())*pow(10,-6);
    
    // Nucleation density
    val = parseInput(InputData, "Heterogeneous nucleation density");
    double NRead = atof(val.c_str());
    NMax = NRead*pow(10,12);
    
    // Mean nucleation undercooling
    val = parseInput(InputData, "Mean nucleation undercooling");
    dTN = atof(val.c_str());
    
    // Standard deviation of nucleation undercooling
    val = parseInput(InputData, "Standard deviation of nucleation undercooling");
    dTsigma = atof(val.c_str());
    
    // Path to output
    PathToOutput = parseInput(InputData, "Path to output");

    // Output base file name
    OutputFile = parseInput(InputData, "Output file base name");
    
    // File of grain orientations
    GrainOrientationFile = parseInput(InputData, "File of grain orientations");
    GrainOrientationFile = FilePath + "/Substrate/" + GrainOrientationFile;

    if (SimulationType == "R") {
        // Read input arguments for a reduced temperature data format solidification problem
        if (id == 0) cout << "CA Simulation using reduced temperature data from file(s)" << endl;
        // Heat transport mesh size
        val = parseInput(InputData, "Heat transport data mesh size");
        HT_deltax = atof(val.c_str())*pow(10,-6);
        if (id == 0) cout << "Mesh size of the temperature data is " << HT_deltax << " microns" << endl;
        
        // Time step (s)
        val = parseInput(InputData, "Time step");
        deltat = atof(val.c_str())*pow(10,-6);
        if (id == 0) cout << "The time step is " << val << " microseconds" << endl;
        
        // Name of substrate file
        SubstrateFileName = parseInput(InputData, "Substrate file name");
        SubstrateFileName = FilePath + "/Substrate/" + SubstrateFileName;
        if (id == 0) cout << "The substrate file used is " << SubstrateFileName << endl;
        
        // Usage of burst buffer/parallel file system for temperature data - if 'Y', values for LayerHeight, NumberOfLayers are unused, temperature filename is overwritten
        TruchasMultilayer = parseInputBool(InputData, "Burst buffer temperature data");
        if (TruchasMultilayer) {
            // Location of temperature data and name of files given in script
            char* InPath;
            char* InFileName;
            InPath = getenv ("PATH_TO_INPUT");
            InFileName = getenv ("TRUCHAS_INPUT_FILENAME");
            string InPathS(InPath);
            string InFileNameS(InFileName);
            tempfile.clear();
            tempfile = InPathS;
            tempfile += InFileNameS;
            if (id == 0) cout << "Using script parameters for temperature data transfer, overriding CA input file parameters based on input data" << endl;
            if (id == 0) cout << "Reading temperature files at location/starting with " << tempfile << endl;
            // Uses Truchas
            TemperatureDataSource = "T";

            // File format using the direct coupling involves a different file for each layer and Truchas
            // NumberOfLayers/NumberOfTruchasRanks require checking files, LayerHeight unknown until files are read
            bool LayerCheck = true;
            NumberOfLayers = 0;
            while (LayerCheck) {
                string LayerNumberString;
                if (NumberOfLayers <= 9) LayerNumberString = "0000" + to_string(NumberOfLayers);
                else if ((NumberOfLayers >= 10)&&(NumberOfLayers <= 99)) LayerNumberString = "000" + to_string(NumberOfLayers);
                else if ((NumberOfLayers >= 100)&&(NumberOfLayers <= 999)) LayerNumberString = "00" + to_string(NumberOfLayers);
                else if ((NumberOfLayers >= 1000)&&(NumberOfLayers <= 9999)) LayerNumberString = "0" + to_string(NumberOfLayers);
                else LayerNumberString = to_string(NumberOfLayers);
                ifstream TestOpen;
                string TestFile = tempfile + "." + LayerNumberString + ".00000";
                TestOpen.open(TestFile);
                if (TestOpen) {
                    // File exists, try the next layer
                    NumberOfLayers++;
                }
                else {
                    // The previous layer was the last one - exit the loop
                    LayerCheck = false;
                }
            }
            // Find number of heat transport model MPI ranks with data
            bool RankCheck = true;
            NumberOfTruchasRanks = 0;
            while (RankCheck) {
                string RankNumberString;
                if (NumberOfTruchasRanks <= 9) RankNumberString = "0000" + to_string(NumberOfTruchasRanks);
                else if ((NumberOfTruchasRanks >= 10)&&(NumberOfTruchasRanks <= 99)) RankNumberString = "000" + to_string(NumberOfTruchasRanks);
                else if ((NumberOfTruchasRanks >= 100)&&(NumberOfTruchasRanks <= 999)) RankNumberString = "00" + to_string(NumberOfTruchasRanks);
                else if ((NumberOfTruchasRanks >= 1000)&&(NumberOfTruchasRanks <= 9999)) RankNumberString = "0" + to_string(NumberOfTruchasRanks);
                else RankNumberString = to_string(NumberOfTruchasRanks);
                ifstream TestOpen;
                string TestFile = tempfile + ".00000." + RankNumberString;
                TestOpen.open(TestFile);
                if (TestOpen) {
                    // File exists, try the next layer
                    NumberOfTruchasRanks++;
                }
                else {
                    // The previous layer was the last one - exit the loop
                    RankCheck = false;
                }
            }
            if (id == 0) cout << "Reading in " << NumberOfLayers << " layers of Truchas data as printed from " << NumberOfTruchasRanks << " ranks" << endl;
        }
        else {
            // File containing temperature data
            tempfile = parseInput(InputData, "Temperature filename");
            tempfile = FilePath + "/Temperatures/" + tempfile;

            // Temperature files in series
            val = parseInput(InputData, "Number of temperature files");
            TempFilesInSeries = stoi(val,nullptr,10);
            
            if (id == 0) {
                cout << "Temperature data file(s) is/are " << tempfile << " , and there are " << TempFilesInSeries << " in the series" << endl;
            }
            
            // Usage of second set of wall cells around temperature field (for spot melt problems, where the melt pool boundaries are right at the walls)
            ExtraWalls = parseInputBool(InputData, "Extra set of wall cells around temperature field");
            
            // OpenFOAM or Truchas as temperature data source (for non-scipt based coupling)?
            TemperatureDataSource = parseInput(InputData, "Source of input length unit");
            
            if (id == 0) {
                if (TemperatureDataSource == "O")
                    cout << "Units of temperature coordinates are meters (OpenFOAM)" << endl;
                else if (TemperatureDataSource == "T")
                    cout << "Units of temperature coordinates are millimeters (Truchas)" << endl;
                else
                    throw std::runtime_error("Input \"Source of input length unit\" must be \"T\" or \"O\".");
            }
            
            // Number of layers (for non script-based coupling)
            val = parseInput(InputData, "Number of layers");
            NumberOfLayers = stoi(val,nullptr,10);
            
            // Layer height (for non script-based coupling)
            val = parseInput(InputData, "Offset between layers");
            LayerHeight = stoi(val,nullptr,10);
            if (id == 0) cout << "A total of " << NumberOfLayers << " of solidification offset by " << LayerHeight << " CA cells will be simulated" << endl;
            
        }
    }
    else if (SimulationType == "C") {
        // Read input arguments for a constrained growth solidification problem
        NumberOfLayers = 1;
        LayerHeight = nz;
        
        // Thermal gradient
        val = parseInput(InputData, "Thermal gradient");
        G = atof(val.c_str());
        
        // Cooling rate
        val = parseInput(InputData, "Cooling rate");
        R = atof(val.c_str());
        if (id == 0) cout << "CA Simulation using a fixed thermal gradient of " << G << " K/m and a cooling rate of " << R << " K/s" << endl;
        
        // Domain size in x
        val = parseInput(InputData, "Domain size in x");
        nx = stoi(val,nullptr,10);
        nx = nx+2;
        
        // Domain size in y
        val = parseInput(InputData, "Domain size in y");
        ny = stoi(val,nullptr,10);
        ny = ny+2;
        
        // Domain size in z
        val = parseInput(InputData, "Domain size in z");
        nz = stoi(val,nullptr,10);
        nz = nz+2;
        if (id == 0) cout << "The domain size is " << nx << " by " << ny << " by " << nz << " cells" << endl;
        
        // delta t (using ratio between cooling rate R, thermal gradient G, and cell size delta x
        val = parseInput(InputData, "\"N\"");
        int NRatio = stoi(val,nullptr,10);
        deltat = deltax/(NRatio*(R/G));
        if (id == 0) cout << "The time step is " << deltat*pow(10,6) << " microseconds" << endl;
        
        // Fraction of bottom surface sites active
        val = parseInput(InputData, "Fraction of surface sites active");
        FractSurfaceSitesActive = atof(val.c_str());
        if (id == 0) cout << "The fraction of CA cells at the bottom surface that are active is " << FractSurfaceSitesActive << endl;
    }
    // Which files should be printed?
    string FileYN;
    // Orientations file
    FilesToPrint[0] = parseInputBool(InputData, "Print file of grain orientations");
    // csv file of grain ids
    FilesToPrint[1] = parseInputBool(InputData, "Print csv file of grain id values");
    // csv file of x,y,z,grainid for ExaConstit
    FilesToPrint[2] = parseInputBool(InputData, "Print csv file of ExaConstit-formatted grain id values");
    // vtk file of grain misorientations
    FilesToPrint[3] = parseInputBool(InputData, "Print Paraview vtk file of grain misorientation values");
    // file of grain areas
    FilesToPrint[4] = parseInputBool(InputData, "Print file of grain area values");
    // file of weighted grain areas
    FilesToPrint[5] = parseInputBool(InputData, "Print file of weighted grain area value");
    if ((!(FilesToPrint[0]))&&(!(FilesToPrint[1]))&&(!(FilesToPrint[2]))&&(!(FilesToPrint[3]))&&(!(FilesToPrint[4]))&&(!(FilesToPrint[5]))) PrintFilesYN = false;
        else PrintFilesYN = true;
    InputData.close();
    
    if (id == 0) {
        cout << "Decomposition Strategy is " << DecompositionStrategy << endl;
        cout << "Material simulated is " << MaterialFile << ", interfacial response function constants are A = " << AConst << ", B = " << BConst << ", C = " << CConst << ", and D = " << DConst << endl;
        cout << "CA cell size is " << deltax*pow(10,6) << " microns" << endl;
        cout << "Nucleation density is " << NRead << " x 10^12 per m^3" << endl;
        cout << "Mean nucleation undercooling is " << dTN << " K, standard deviation of distribution is " << dTsigma << "K" << endl;
    }
    
}

void ParallelMeshInit(int DecompositionStrategy, ViewI_H NeighborX, ViewI_H NeighborY, ViewI_H NeighborZ, ViewI2D_H ItList, string SimulationType, int id, int np, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset,int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn, int &MyLeftOut, int &MyRightIn, int &MyRightOut, double &deltax, int &nx, int &ny, int &nz, int &ProcessorsInXDirection, int &ProcessorsInYDirection, string tempfile, float &XMin, float &XMax, float &YMin, float &YMax, float &ZMin, float &ZMax, string TemperatureDataSource, int &LayerHeight, int NumberOfLayers, int TempFilesInSeries, float* ZMinLayer, float* ZMaxLayer, int* FirstValue, vector <float> &RawData, bool TruchasMultilayer, int NumberOfTruchasRanks) {
        
    // Assignment of neighbors around a cell "X" is as follows (in order of closest to furthest from cell "X")
    // Neighbors 0 through 8 are in the -Y direction
    // Neighbors 9 through 16 are in the XY plane with cell X
    // Neighbors 17 through 25 are in the +Y direction

    NeighborX(0) = 0; NeighborY(0) = -1; NeighborZ(0) = 0;
    NeighborX(1) = 1; NeighborY(1) = -1; NeighborZ(1) = 0;
    NeighborX(2) = -1; NeighborY(2) = -1; NeighborZ(2) = 0;
    NeighborX(3) = 0; NeighborY(3) = -1; NeighborZ(3) = 1;
    NeighborX(4) = 0; NeighborY(4) = -1; NeighborZ(4) = -1;
    NeighborX(5) = -1; NeighborY(5) = -1; NeighborZ(5) = 1;
    NeighborX(6) = 1; NeighborY(6) = -1; NeighborZ(6) = 1;
    NeighborX(7) = -1; NeighborY(7) = -1; NeighborZ(7) = -1;
    NeighborX(8) = 1; NeighborY(8) = -1; NeighborZ(8) = -1;
    
    NeighborX(9) = 0; NeighborY(9) = 0; NeighborZ(9) = 1;
    NeighborX(10) = 0; NeighborY(10) = 0; NeighborZ(10) = -1;
    NeighborX(11) = 1; NeighborY(11) = 0; NeighborZ(11) = 1;
    NeighborX(12) = -1; NeighborY(12) = 0; NeighborZ(12) = 1;
    NeighborX(13) = 1; NeighborY(13) = 0; NeighborZ(13) = -1;
    NeighborX(14) = -1; NeighborY(14) = 0; NeighborZ(14) = -1;
    NeighborX(15) = 1; NeighborY(15) = 0; NeighborZ(15) = 0;
    NeighborX(16) = -1; NeighborY(16) = 0; NeighborZ(16) = 0;
    
    NeighborX(17) = 0; NeighborY(17) = 1; NeighborZ(17) = 0;
    NeighborX(18) = 1; NeighborY(18) = 1; NeighborZ(18) = 0;
    NeighborX(19) = -1; NeighborY(19) = 1; NeighborZ(19) = 0;
    NeighborX(20) = 0; NeighborY(20) = 1; NeighborZ(20) = 1;
    NeighborX(21) = 0; NeighborY(21) = 1; NeighborZ(21) = -1;
    NeighborX(22) = 1; NeighborY(22) = 1; NeighborZ(22) = 1;
    NeighborX(23) = -1; NeighborY(23) = 1; NeighborZ(23) = 1;
    NeighborX(24) = 1; NeighborY(24) = 1; NeighborZ(24) = -1;
    NeighborX(25) = -1; NeighborY(25) = 1; NeighborZ(25) = -1;
    
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
        ItList(0,i) = i;
    }
    // If Y coordinate is on lower edge, Case 1: iteration over only neighbors 9-25 possible
    int Pos = 0;
    for (int i=0; i<=25; i++) {
        if (NeighborY(i) != -1) {
            ItList(1,Pos) = i;
            Pos++;
        }
    }
    // If Y coordinate is on upper edge, Case 2: iteration over only neighbors 0-16 possible
    Pos = 0;
    for (int i=0; i<=25; i++) {
        if (NeighborY(i) != 1) {
            ItList(2,Pos) = i;
            Pos++;
        }
    }
    // If X coordinate is on lower edge, Case 3: iteration over only neighbors 0,1,3,4,6,8,9,10,11,13,15,17,18,20,21,22,24
    Pos = 0;
    for (int i=0; i<=25; i++) {
        if (NeighborX(i) != -1) {
            ItList(3,Pos) = i;
            Pos++;
        }
    }
    // If X coordinate is on upper edge, Case 4: iteration over only neighbors 0,2,3,4,5,7,9,10,12,14,16,17,19,20,21,23,25
    Pos = 0;
    for (int i=0; i<=25; i++) {
        if (NeighborX(i) != 1) {
            ItList(4,Pos) = i;
            Pos++;
        }
    }
    // If X/Y coordinates are on lower edge, Case 5: iteration over only neighbors 9,10,11,13,15,17,18,20,21,22,24
    Pos = 0;
    for (int i=0; i<=25; i++) {
        if ((NeighborX(i) != -1)&&(NeighborY(i) != -1)) {
            ItList(5,Pos) = i;
            Pos++;
        }
    }
    // If X coordinate is on upper edge/Y on lower edge, Case 6:
    Pos = 0;
    for (int i=0; i<=25; i++) {
        if ((NeighborX(i) != 1)&&(NeighborY(i) != -1)) {
            ItList(6,Pos) = i;
            Pos++;
        }
    }
    // If X coordinate is on lower edge/Y on upper edge, Case 7:
    Pos = 0;
    for (int i=0; i<=25; i++) {
        if ((NeighborX(i) != -1)&&(NeighborY(i) != 1)) {
            ItList(7,Pos) = i;
            Pos++;
        }
    }
    // If X/Y coordinates are on upper edge, Case 8:
    Pos = 0;
    for (int i=0; i<=25; i++) {
        if ((NeighborX(i) != 1)&&(NeighborY(i) != 1)) {
            ItList(8,Pos) = i;
            Pos++;
        }
    }

    if (SimulationType == "R") {
        // nx, ny, and nz are unknown prior to reading the relevant temperature file(s)
        // Determine mesh size needed based on OpenFOAM/Truchas/Analytical model data
        // Read geometry data from OpenFOAM/Truchas
        // int nx_HT, ny_HT, nz_HT; // OpenFOAM/Truchas mesh limits

        XMin = 1000000.0;
        YMin = 1000000.0;
        ZMin = 1000000.0;
        XMax = -1000000.0;
        YMax = -1000000.0;
        ZMax = -1000000.0;

        // Read and find bounds of temperature data
        // Differences between reading data from the burst buffer/parallel file systems vs as a standard series:
        int LayersToRead, HTRanksWithData;
        if (TruchasMultilayer) {
            LayersToRead = NumberOfLayers;
            HTRanksWithData = NumberOfTruchasRanks;
        }
        else {
            LayersToRead = min(NumberOfLayers,TempFilesInSeries); // was given in input file
            HTRanksWithData = 1; // Data was previously consolidated
        }
        
	std::size_t NumberOfHTDataPoints = 0; // Counter variable for pieces of data read
        
        // Is the input in m and s (OpenFOAM) or mm and ms (Truchas)?
        double UnitConversion;
        if (TemperatureDataSource == "O") UnitConversion = 1;
        else UnitConversion = 1000;
        
        for (int LayerReadCount=1; LayerReadCount<=LayersToRead; LayerReadCount++) {
        
            ZMinLayer[LayerReadCount-1] = 1000000.0;
            ZMaxLayer[LayerReadCount-1] = -1000000.0;
            FirstValue[LayerReadCount-1] = NumberOfHTDataPoints;
            if (id == 0) cout << "Start value layer " << LayerReadCount-1 << " is " << NumberOfHTDataPoints << endl;
            //float SmallestTime_ThisLayer = 10000000;
            //float LargestTime_ThisLayer =  10000000;
            for (int n=0; n<HTRanksWithData; n++) {
                
                string tempfile_thislayer;
                
                if (!(TruchasMultilayer)) {
                    if (TempFilesInSeries > 1) {
                        string NextLayerFileS = to_string(LayerReadCount);
                        int NextLayerFile = LayerReadCount % TempFilesInSeries;
                        if (NextLayerFile == 0) NextLayerFile = TempFilesInSeries;
                        NextLayerFileS = to_string(NextLayerFile);
                        std::size_t found = tempfile.find_last_of("/");
                        string FPath = tempfile.substr(0,found+1);
                        string FName = tempfile.substr(found+1);
                        tempfile_thislayer = FPath + NextLayerFileS + FName;
                    }
                    else {
                        tempfile_thislayer = tempfile;
                    }
                }
                else {
                   // Read through all files that have data for this layer, with knowledge of where the files are located
                    string LayerNumberString;
                    if (LayerReadCount-1 <= 9) LayerNumberString = "0000" + to_string(LayerReadCount-1);
                    else if ((LayerReadCount-1 >= 10)&&(LayerReadCount-1 <= 99)) LayerNumberString = "000" + to_string(LayerReadCount-1);
                    else if ((LayerReadCount-1 >= 100)&&(LayerReadCount-1 <= 999)) LayerNumberString = "00" + to_string(LayerReadCount-1);
                    else if ((LayerReadCount-1 >= 1000)&&(LayerReadCount-1 <= 9999)) LayerNumberString = "0" + to_string(LayerReadCount-1);
                    else LayerNumberString = to_string(LayerReadCount-1);
        
                    string RankNumberString;
                    if (n <= 9) RankNumberString = "0000" + to_string(n);
                    else if ((n >= 10)&&(n <= 99)) RankNumberString = "000" + to_string(n);
                    else if ((n >= 100)&&(n <= 999)) RankNumberString = "00" + to_string(n);
                    else if ((n >= 1000)&&(n <= 9999)) RankNumberString = "0" + to_string(n);
                    else RankNumberString = to_string(n);
                    
                    tempfile_thislayer = tempfile + "." + LayerNumberString + "." + RankNumberString;
                }

                //if (id == 0) cout << "Reading file " << tempfile_thislayer << endl;
                ifstream Geom;

                Geom.open(tempfile_thislayer);
                if (TruchasMultilayer) {
                    // 3 header lines
                    string s;
                    getline(Geom,s);
                    getline(Geom,s);
                    getline(Geom,s);
                }
                while (!Geom.eof()) {
                    string s;
                    getline(Geom,s);
                    if (s.empty()) break;
                    string XVal, YVal, ZVal, TLVal, TSVal;
                    int Subdivisions[4];
                    int SubdivisionCount = 0;
                    for (std::size_t i=1; i<s.length(); i++) {
                        char ThisChar = s.at(i);
                        char PrevChar = s.at(i-1);
                        // If this character is blank and the previous one was not, this is a spot to subdivide the string
                        if ((isblank(ThisChar))&&(!isblank(PrevChar))) {
                            Subdivisions[SubdivisionCount] = i;
                            SubdivisionCount++;
                            if (SubdivisionCount == 4) {
                                // Last division was made
                                break;
                            }
                        }
                    }
                    //if (id == 0) cout << "Subdivisions at " << Subdivisions[0] << " " << Subdivisions[1] << " " << Subdivisions[2] << " " << Subdivisions[3] << endl;
                    for (int i=0; i<5; i++) {
                        int stringstart;
                        int stringlength;
                        if (i == 0) stringstart = 0;
                        else stringstart = Subdivisions[i-1];
                        if (i == 4) stringlength = s.length();
                        else stringlength = Subdivisions[i]-stringstart;
                        string MeshDataS = s.substr(stringstart,stringlength);
                        //if (id == 0) cout << s << " SUBVAL = " << MeshDataS << endl;
                        float MeshData = atof(MeshDataS.c_str());
                        RawData[NumberOfHTDataPoints] = MeshData/UnitConversion;
                        if (i == 0) {
                            if (XMin > RawData[NumberOfHTDataPoints]) XMin = RawData[NumberOfHTDataPoints];
                            if (XMax < RawData[NumberOfHTDataPoints]) XMax = RawData[NumberOfHTDataPoints];
                        }
                        else if (i == 1) {
                            //if (id == 0) cout << "Y " << MeshData << endl;
                            if (YMin > RawData[NumberOfHTDataPoints]) YMin = RawData[NumberOfHTDataPoints];
                            if (YMax < RawData[NumberOfHTDataPoints]) YMax = RawData[NumberOfHTDataPoints];
                        }
                        else if (i == 2) {
                            // If not reading from the parallel file system/burst buffer directly from Truchas, Z values need to be adjusted based on which layer is being simulated (layers are offset by "LayerHeight"
                            // If reading from the parallel file system/burst buffer, Z values are not adjusted
                            float AdjustedZ;
                            if (!(TruchasMultilayer)&&(TempFilesInSeries > 1)) {
                                AdjustedZ = RawData[NumberOfHTDataPoints] + deltax*LayerHeight*(LayerReadCount-1);
                            }
                            else {
                                AdjustedZ = RawData[NumberOfHTDataPoints];
                            }
                            //if (id == 0) cout << LayerReadCount-1 << " " << s << " Value is " << AdjustedZ << endl;
                            if (ZMin > AdjustedZ) ZMin = AdjustedZ;
                            if (ZMax < AdjustedZ) ZMax = AdjustedZ;
                            if (ZMinLayer[LayerReadCount-1] > AdjustedZ) ZMinLayer[LayerReadCount-1] = AdjustedZ;
                            if (ZMaxLayer[LayerReadCount-1] < AdjustedZ) ZMaxLayer[LayerReadCount-1] = AdjustedZ;
                        }
                        else if (i == 3) {
                            // Liquidus time
                            //if (SmallestTime_ThisLayer > RawData[NumberOfHTDataPoints]) SmallestTime_ThisLayer = RawData[NumberOfHTDataPoints];
                        }
                        else if (i == 4) {
                            //if ((id == 0)&&(LayerReadCount == 2)) cout << "Truchas Rank " << n << " " << " " << RawData[NumberOfHTDataPoints-4] << " " << RawData[NumberOfHTDataPoints-3] << " " << RawData[NumberOfHTDataPoints-2] << " " << RawData[NumberOfHTDataPoints-1] << " " << RawData[NumberOfHTDataPoints] << endl;
                            // Solidus time
                            //if (LargestTime_ThisLayer < RawData[NumberOfHTDataPoints]) LargestTime_ThisLayer = RawData[NumberOfHTDataPoints];
                        }
                        if (NumberOfHTDataPoints >= RawData.size()-5) {
                            int OldSize = RawData.size();
                            RawData.resize(OldSize+1000000);
                        }
                        NumberOfHTDataPoints++;
                    }
                }
                Geom.close();
            } // End loop over all files read for a particular layer
            
            // Sync layer Z bounds across processors
            float MyRankMinZ = ZMinLayer[LayerReadCount-1];
            float MyRankMaxZ = ZMaxLayer[LayerReadCount-1];
            float GlobalMinZ, GlobalMaxZ;
            MPI_Allreduce(&MyRankMinZ,&GlobalMinZ,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
            MPI_Allreduce(&MyRankMaxZ,&GlobalMaxZ,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
            ZMinLayer[LayerReadCount-1] = GlobalMinZ;
            ZMaxLayer[LayerReadCount-1] = GlobalMaxZ;
            if (id == 0) cout << "Layer = " << LayerReadCount << " Z Bounds are " << ZMinLayer[LayerReadCount-1] << " " << ZMaxLayer[LayerReadCount-1]  << endl;
            
        } // End loop over all files read for all layers
        RawData.resize(NumberOfHTDataPoints);

        if (!(TruchasMultilayer)) {
            // Extend domain in Z (build) direction if the number of layers are simulated is greater than the number of temperature files read
            if (NumberOfLayers > TempFilesInSeries) {
                for (int LayerReadCount=TempFilesInSeries; LayerReadCount<NumberOfLayers; LayerReadCount++) {
                    if (TempFilesInSeries == 1) {
                        // Only one temperature file was read, so the upper Z bound should account for an additional "NumberOfLayers-1" worth of data
                        // Since all layers have the same temperature data, each layer's "ZMinLayer" is just translated from that of the first layer
                        ZMinLayer[LayerReadCount] = ZMinLayer[LayerReadCount-1] + deltax*LayerHeight;
                        ZMaxLayer[LayerReadCount] = ZMaxLayer[LayerReadCount-1] + deltax*LayerHeight;
                        FirstValue[LayerReadCount] = FirstValue[LayerReadCount-1];
                        ZMax += deltax*LayerHeight;
                    }
                    else {
                        // "TempFilesInSeries" temperature files was read, so the upper Z bound should account for an additional "NumberOfLayers-TempFilesInSeries" worth of data
                        int RepeatedFile = (LayerReadCount) % TempFilesInSeries;
                        int RepeatUnit = LayerReadCount/TempFilesInSeries;
                        ZMinLayer[LayerReadCount] = ZMinLayer[RepeatedFile] + RepeatUnit*TempFilesInSeries*deltax*LayerHeight;
                        ZMaxLayer[LayerReadCount] = ZMaxLayer[RepeatedFile] + RepeatUnit*TempFilesInSeries*deltax*LayerHeight;
                        FirstValue[LayerReadCount] = FirstValue[RepeatedFile];
                        ZMax += deltax*LayerHeight;
                    }
                }
            }
        }
        else {
            // Determine offset between layers based on data read (assuming this is the same for each layer read)
            if (NumberOfLayers > 1) LayerHeight = round((ZMaxLayer[1] - ZMaxLayer[0])/deltax);
        }

        // Ratio of temperature-CA grid sizes
        //int HTratio = round(HT_deltax/deltax);
        
        // CA nodes in each direction (+2 for wall cells at the boundaries) (+2 for solid cells at X/Y boundaries, +1 for solid cells at lower Z boundary)
        nx = round((XMax-XMin)/deltax) + 1 + 4;
        ny = round((YMax-YMin)/deltax) + 1 + 4;
        nz = round((ZMax-ZMin)/deltax) + 1 + 3;

        if (id == 0) {
            cout << "Domain size: " << nx << " by " << ny << " by " << nz << endl;
            cout << "X Limits of domain: " << XMin << " and " << XMax << endl;
            cout << "Y Limits of domain: " << YMin << " and " << YMax << endl;
            cout << "Z Limits of domain: " << ZMin << " and " << ZMax << endl;
            cout << "================================================================" << endl;
        }
    }

    InitialDecomposition(DecompositionStrategy, nx, ny, ProcessorsInXDirection, ProcessorsInYDirection, id, np, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyLeftOut, MyRightIn, MyRightOut);

    MyXOffset = XOffsetCalc(id,nx,ProcessorsInXDirection,ProcessorsInYDirection,DecompositionStrategy);
    MyXSlices = XMPSlicesCalc(id,nx,ProcessorsInXDirection,ProcessorsInYDirection,DecompositionStrategy);

    MyYOffset = YOffsetCalc(id,ny,ProcessorsInYDirection,np,DecompositionStrategy);
    MyYSlices = YMPSlicesCalc(id,ny,ProcessorsInYDirection,np,DecompositionStrategy);
   // cout << "ID = " << id << " X RANGE = " << MyXOffset << " TO = " << MyXOffset+MyXSlices-1 << " ; YRANGE = " << MyYOffset << " TO = " << MyYOffset+MyYSlices-1 << endl;
}

void TempInit(int layernumber, int TempFilesInSeries, double G, double R, string SimulationType, int, int id, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset, double deltax, double HT_deltax, double deltat, int &nx, int &ny, int &nz, ViewI_H CritTimeStep, ViewF_H UndercoolingChange, ViewF_H UndercoolingCurrent, float XMin, float YMin, float ZMin, bool* Melted, float* ZMinLayer, float* ZMaxLayer, int LayerHeight, int NumberOfLayers, int &nzActive, int &ZBound_Low, int &ZBound_High, int* FinishTimeStep, double FreezingRange, ViewI_H LayerID, int* FirstValue, vector <float> RawData, bool TruchasMultilayer) {

    if (SimulationType == "C") {
        
        // Contrained solidification test problem
        ZBound_Low = 0;
        ZBound_High = nz-1;
        nzActive = nz;
        
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
            
        // Initialize temperature views to 0
        for (int k=0; k<nz; k++) {
           for (int i=0; i<MyXSlices; i++) {
               for (int j=0; j<MyYSlices; j++) {
                   int Coord3D1D = k*MyXSlices*MyYSlices + i*MyYSlices + j;
                   CritTimeStep(Coord3D1D) = 0;
                   UndercoolingChange(Coord3D1D) = 0.0;
                   UndercoolingCurrent(Coord3D1D) = 0.0;
               }
           }
        }
       
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

       if (layernumber == -1) {
           // No sites have melted yet
           for (int i=0; i<MyXSlices*MyYSlices*nz; i++) {
               Melted[i] = false;
               LayerID(i) = -1;
           }
       }

       float LayerwiseTSOffset = 0;
       if (id == 0) cout << "Raw data has " << RawData.size() << " values" << endl;
       
       for (int LayerCounter=0; LayerCounter<NumberOfLayers; LayerCounter++) {
           
           double SmallestTime = 1000000000;
           double SmallestTime_Global = 1000000000;
           double LargestTime = 0;
           double LargestTime_Global = 0;
           
           // How many CA cells in the vertical direction are needed to hold this layer's temperature data?
           int nzTempValuesThisLayer = round((ZMaxLayer[LayerCounter] - ZMinLayer[LayerCounter])/deltax) + 1; // (note this doesn't include the 2 rows of wall/active cells at the bottom surface)
           if (id == 0) cout << "Initializing temporary temperature data structures with " << nzTempValuesThisLayer << " cells in z direction" << endl;
           if (id == 0) cout << "Layer " << LayerCounter << " rank " << id << " ZMin this layer is " << ZMinLayer[LayerCounter] << endl;
           vector <vector <vector <double> > > CritTS, CritTL;
           for (int k=0; k<nzTempValuesThisLayer; k++) {
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

           // Data was already read into the "RawData" temporary data structure
           // Determine which section of "RawData" is relevant for this layer of the overall domain
           int StartRange = FirstValue[LayerCounter];
           int EndRange;
           if (!(TruchasMultilayer)) {
               int RepeatedFile = LayerCounter % TempFilesInSeries;
               if (RepeatedFile != TempFilesInSeries-1) EndRange = FirstValue[RepeatedFile+1];
               else {
                   EndRange = RawData.size();
               }
           }
           else {
               if (LayerCounter < NumberOfLayers-1) EndRange = FirstValue[LayerCounter+1];
               else EndRange = RawData.size();
           }
           int XInt, YInt, ZInt;
           if (id == 0) cout << "Range for layer " << LayerCounter << " is " << StartRange << " to " << EndRange << endl;
           MPI_Barrier(MPI_COMM_WORLD);
           for (int i=StartRange; i<EndRange; i++) {

               int Pos = i % 5;
               if (Pos == 0) {
                   XInt = round((RawData[i]-XMin)/deltax) + 2;
               }
               else if (Pos == 1) {
                   YInt = round((RawData[i]-YMin)/deltax) + 2;
               }
               else if (Pos == 2) {
                   if (!(TruchasMultilayer))  ZInt = round((RawData[i] + deltax*LayerHeight*LayerCounter - ZMinLayer[LayerCounter])/deltax);
                   else ZInt = round((RawData[i] - ZMinLayer[LayerCounter])/deltax);
               }
               else if (Pos == 3) {
                   // Determine if this X and Y fall in this rank's range of interest
                   if ((XInt >= LowerXBound)&&(XInt <= UpperXBound)&&(YInt >= LowerYBound)&&(YInt <= UpperYBound)) {
                       if (RawData[i] <= 0) RawData[i] = 1;
                       CritTL[ZInt][XInt-LowerXBound][YInt-LowerYBound] = RawData[i];
                       //if ((id == 0)&&(LayerCounter < 2)) cout << "ID = " << id << " i = " << XInt << " j = " << YInt << " k = " << ZInt << endl;
                       if (RawData[i] < SmallestTime) SmallestTime = RawData[i];
                   }
               }
               else if (Pos == 4) {
                   if ((XInt >= LowerXBound)&&(XInt <= UpperXBound)&&(YInt >= LowerYBound)&&(YInt <= UpperYBound)) {
                       if (RawData[i] <= 0) RawData[i] = 1;
                       CritTS[ZInt][XInt-LowerXBound][YInt-LowerYBound] = RawData[i];
                       if (RawData[i] > LargestTime) LargestTime = RawData[i];
                   }
               }
           }
           // If reading data from files without a script, time values start at 0 for each layer
           // If reading data with input from a script time values each layer are continuous, are should be renormalized to 0 for each layer
           MPI_Reduce(&LargestTime, &LargestTime_Global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
           MPI_Bcast(&LargestTime_Global,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
           MPI_Reduce(&SmallestTime, &SmallestTime_Global, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
           MPI_Bcast(&SmallestTime_Global,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
           if (TruchasMultilayer) {
               if (id == 0) cout << "Smallest time for layer " << LayerCounter << " is " << SmallestTime_Global << endl;
               LayerwiseTSOffset = SmallestTime_Global;
           }
           if (id == 0) cout << "Largest time globally for layer " << LayerCounter << " is " << LargestTime_Global << endl;
           FinishTimeStep[LayerCounter] = round((LargestTime_Global-LayerwiseTSOffset)/deltat);
           if (id == 0) cout << " Layer " << LayerCounter << " FINISH TIME STEP IS " << FinishTimeStep[LayerCounter] << endl;
           //TimeOffset = LargestTime_Global;
           if (id == 0) cout << "Layer " << LayerCounter << " temperatures read" << endl;
           // Data interpolation between heat transport and CA grids, if necessary
           if (HTtoCAratio != 1) {
               for (int k=0; k<nzTempValuesThisLayer; k++) {
                   int LowZ = k - (k % HTtoCAratio);
                   int HighZ = LowZ + HTtoCAratio;
                   double FHighZ = (double)(k - LowZ)/(double)(HTtoCAratio);
                   double FLowZ = 1.0 - FHighZ;
                   if (HighZ > nzTempValuesThisLayer-1) HighZ = LowZ;
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
                           //if ((id == 0)&&(k == 16)) cout << LowZ << " " << HighZ << " nzTempvalsthislayer = " << nzTempValuesThisLayer << endl;
                       }
                   }
               }
           }
           MPI_Barrier(MPI_COMM_WORLD);
           if (id == 0) cout << "Interpolation done" << endl;
           
           // Convert CritTL, CritTS matrices into CritTimeStep and UndercoolingChange (change in undercooling with time step)
           // "ZMin" is the global Z coordinate that corresponds to cells at Z = 2 (Z = 0 is the domain's bottom wall, Z = 1 are the active cells just outside of the melt pool)
           // "ZMax" is the global Z coordinate that corresponds to cells at Z = nz-2 (Z = nz-1 is the domain's top wall)
           if (LayerCounter == 0) {
               ZBound_Low = 0;
               ZBound_High = round((ZMinLayer[LayerCounter]-ZMin)/deltax) + 2 + nzTempValuesThisLayer-1;
               nzActive = ZBound_High - ZBound_Low + 1;
           }
           if (id == 0) cout << "Layer " << LayerCounter << " data belongs to global z coordinates of " << round((ZMinLayer[LayerCounter]-ZMin)/deltax)+2 << " through " << round((ZMinLayer[LayerCounter]-ZMin)/deltax) + 2 + nzTempValuesThisLayer-1 << endl;
           
           for (int k=0; k<nzTempValuesThisLayer; k++) {
               for (int ii=LowerXBound; ii<=UpperXBound; ii++) {
                   for (int jj=LowerYBound; jj<=UpperYBound; jj++) {
                       if ((ii >= MyXOffset)&&(ii < MyXOffset+MyXSlices)&&(jj >= MyYOffset)&&(jj < MyYOffset+MyYSlices)) {
                           int Adj_i = ii - MyXOffset;
                           int Adj_j = jj - MyYOffset;
                           //if ((id == 16)&&(CritTL[k][ii-LowerXBound][jj-LowerYBound] < 4.5)&&(LayerCounter == 19)) cout << " i = " << Adj_i << " j = " << Adj_j << " k = " << k << " Val = " << CritTL[k][ii-LowerXBound][jj-LowerYBound] << " LTS = " << LayerwiseTSOffset << endl;
                           double CTLiq = CritTL[k][ii-LowerXBound][jj-LowerYBound] - LayerwiseTSOffset;
                           double CTSol = CritTS[k][ii-LowerXBound][jj-LowerYBound] - LayerwiseTSOffset;
                           //if (CTSol < CTLiq) cout << "ID = " << id << " CTLiq = " << CTLiq << " CTSol = " << CTSol << endl;
                           if (CTLiq > 0)  {
                               // Where does this layer's temperature data belong on the global (including all layers) grid?
                               // Adjust Z coordinate by ZMin
                               int ZOffset = round((ZMinLayer[LayerCounter]-ZMin)/deltax) + k + 2;
                               int Coord3D1D = ZOffset*MyXSlices*MyYSlices + Adj_i*MyYSlices + Adj_j;
                               Melted[Coord3D1D] = true;
                               CritTimeStep(Coord3D1D) = round(CTLiq/deltat);
                               LayerID(Coord3D1D) = LayerCounter;
                               UndercoolingChange(Coord3D1D) = FreezingRange*deltat/(CTSol - CTLiq);
                               //if ((id == 1)&&(Adj_i == 68)&&(Adj_j == 5)&&(LayerCounter == 1)) cout << " Z = " << ZOffset << " , CTS = " << CritTimeStep(Coord3D1D) << " UC = " << UndercoolingChange(Coord3D1D) << " LiqT/SolT = " << CTLiq << " " << CTSol << " ts = " << deltat << " FR = " << FreezingRange << endl;
//                               if (UndercoolingChange(Coord3D1D) > 0.25) {
//                                   UndercoolingChange(Coord3D1D) = 0.25;
//                               }
                           }
                       }
                   }
               }
           }
        } // End read over all temperature files and placement of data

       if (id == 0) cout << "First layer Z bounds are " << ZBound_Low << " and " << ZBound_High << endl;
    }


}
////*****************************************************************************/

// Initialize grain orientations and unit vectors

void OrientationInit(int id, int NGrainOrientations, ViewI_H GrainOrientation, ViewF_H GrainUnitVector, string GrainOrientationFile) {
    
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
            
            GrainUnitVector(9*i + 3*UVNumber + Comp) = ReadGO;
            
//            GrainUnitVector[18*i + 3*UVNumber + Comp] = ReadGO;
//            GrainUnitVector[18*i + 3*(UVNumber+1) + Comp] = -ReadGO;
            Comp++;
            if (Comp > 2) {
                Comp = 0;
                UVNumber++;
                //UVNumber = UVNumber + 2;
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
        GrainOrientation(h) = GrainOrientation_master[h];
    }
    
}


////*****************************************************************************/
////*****************************************************************************/
//
///*
// Initializes cell types where the substrate comes from a file
//*/

void GrainInit(int layernumber, string SimulationType, string SubstrateFileName, double FractSurfaceSitesActive, int NGrainOrientations, int DecompositionStrategy, int nx, int ny, int nz, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int id, int np, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, ViewI2D_H ItList, ViewI_H NeighborX, ViewI_H NeighborY, ViewI_H NeighborZ, ViewI_H GrainOrientation, ViewF_H GrainUnitVector, ViewF_H DiagonalLength, ViewI_H CellType, ViewI_H GrainID, ViewF_H CritDiagonalLength, ViewF_H DOCenter, ViewI_H CritTimeStep, ViewF_H UndercoolingChange, bool* Melted, double deltax, double NMax, int &NextLayer_FirstNucleatedGrainID, int &PossibleNuclei_ThisRank, int ZBound_High, int ZBound_Low, string ExtraWalls) {
    
    mt19937_64 gen(id);//2*id);//*234) ; //123);//234); //(id*(layernumber+1));
    uniform_real_distribution<double> dis(0.0, 1.0);
    
    // Convert initial grain spacing to a grain density
    double BulkProb = NMax*deltax*deltax*deltax;
    if (id == 0) cout << "Fraction of heterogenous nucleation sites to potentially be activated: " << BulkProb << endl;
    
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
                int CAGridLocation = k*MyXSlices*MyYSlices + i*MyYSlices + j;
                if ((GlobalX == -1)||(GlobalX == nx)||(GlobalY == -1)||(GlobalY == ny)||(k == 0)||(k == nz-1)) {
                    CellType(CAGridLocation) = Wall;
                    GrainID(CAGridLocation) = 0;
                }
                else {
                    CellType(CAGridLocation) = Solid;
                }
            }
        }
    }
    
    if (SimulationType == "C") {
        
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
        if (id == 0) cout << "Opened substrate file " << SubstrateFileName << endl;
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
        if ((id == 0)&&(nzS < nz)) cout << "Warning: only " << nzS << " layers of substrate data for a simulation of " << nz << " total layers" << endl;

        // Assign GrainID values to cells that are part of the substrate
        // Cells that border the melted region are type active, others are type solid
        for (int k=0; k<nzS; k++) {
            if (k == nz) break;
            for (int j=0; j<nyS; j++) {
                for (int i=0; i<nxS; i++) {
                    string GIDVal;
                    getline(Substrate,GIDVal);
                    if ((i >= Substrate_LowX)&&(i < Substrate_HighX)&&(j >= Substrate_LowY)&&(j < Substrate_HighY)) {
                        int CAGridLocation;
                        CAGridLocation = k*MyXSlices*MyYSlices + (i-MyXOffset)*MyYSlices + (j-MyYOffset);
                        if (CritTimeStep(CAGridLocation) == 0) {
                            GrainID(CAGridLocation) = stoi(GIDVal,nullptr,10);
                        }
                        else {
                            GrainID(CAGridLocation) = 0;
                        }
                    }
                }
            }
        }
        Substrate.close();
//        if (nz > nzS) {
//            for (int k=nzS; k<nz; k++) {
//                for (int j=0; j<nyS; j++) {
//                    for (int i=0; i<nxS; i++) {
//                        if ((i >= Substrate_LowX)&&(i < Substrate_HighX)&&(j >= Substrate_LowY)&&(j < Substrate_HighY)) {
//                            int CAGridLocation;
//                            CAGridLocation = k*MyXSlices*MyYSlices + (i-MyXOffset)*MyYSlices + (j-MyYOffset);
//                            if (CritTimeStep(CAGridLocation) == 0) {
//                                GrainID(CAGridLocation) = floorf(NGrainOrientations*(float) rand()/RAND_MAX);
//                            }
//                            else {
//                                GrainID(CAGridLocation) = 0;
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        if (id == 0) cout << "Substrate file read complete" << endl;

//         Single crystal substrate
//        int Substrate_LowX = MyXOffset;
//        int Substrate_HighX = MyXOffset+MyXSlices;
//        int Substrate_LowY = MyYOffset;
//        int Substrate_HighY = MyYOffset+MyYSlices;
//        for (int k=0; k<nz; k++) {
//            for (int j=0; j<ny; j++) {
//                for (int i=0; i<nx; i++) {
//                    if ((i >= Substrate_LowX)&&(i < Substrate_HighX)&&(j >= Substrate_LowY)&&(j < Substrate_HighY)) {
//                        int CAGridLocation;
//                        CAGridLocation = k*MyXSlices*MyYSlices + (i-MyXOffset)*MyYSlices + (j-MyYOffset);
//                        if (CritTimeStep(CAGridLocation) == 0) {
//                            GrainID(CAGridLocation) = 9589; //stoi(GIDVal,nullptr,10);
//                        }
//                        else {
//                            GrainID(CAGridLocation) = 0;
//                        }
//                    }
//                }
//            }
//        }
        
        if (ExtraWalls == "Y") {
            if (id == 0) cout << "Extra wall cells around domain" << endl;
            // Extra set of wall cells around edges for spot melt problem
            for (int k=0; k<nz; k++)  {
                for(int i=0; i<MyXSlices; i++) {
                    for(int j=0; j<MyYSlices; j++) {
                        int GlobalX = i + MyXOffset;
                        int GlobalY = j + MyYOffset;
                        if ((GlobalX == 0)||(GlobalX == nx-1)||(GlobalY == 0)||(GlobalY == ny-1)||(GlobalX == 1)||(GlobalX == nx-2)||(GlobalY == 1)||(GlobalY == ny-2)) {
                            int CAGridLocation = k*MyXSlices*MyYSlices + i*MyYSlices + j;
                            CellType(CAGridLocation) = Wall;
                            GrainID(CAGridLocation) = 0;
                        }
                    }
                }
            }
        }
    
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
                                int l = ItList(ItBounds,ll);
                                // Local coordinates of adjacent cell center
                                int MyNeighborX = i + NeighborX(l);
                                int MyNeighborY = j + NeighborY(l);
                                int MyNeighborZ = k + NeighborZ(l);
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

    for (int GlobalZ=1; GlobalZ<nz-1; GlobalZ++) {
        //cout << RankZ << endl;
        for (int RankX=0; RankX<MyXSlices; RankX++) {
            for (int RankY=0; RankY<MyYSlices; RankY++) {
                long int D3D1ConvPositionGlobal = GlobalZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                if (CellType(D3D1ConvPositionGlobal) == Active) {
                    //cout << "Coords" << RankX << " " << RankY << " " << RankZ << endl;
                    
                    // If part of the active domain, calculate Critical diagonal lengths
                    if (GlobalZ <= ZBound_High) {
                        int GlobalX = RankX + MyXOffset;
                        int GlobalY = RankY + MyYOffset;
                        int RankZ = GlobalZ - ZBound_Low;
                        int D3D1ConvPosition = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                        int MyGrainID = GrainID(D3D1ConvPositionGlobal);
                        DiagonalLength(D3D1ConvPosition) = 0.01;
                        DOCenter(3*D3D1ConvPosition) = GlobalX + 0.5;
                        DOCenter(3*D3D1ConvPosition+1) = GlobalY + 0.5;
                        DOCenter(3*D3D1ConvPosition+2) = GlobalZ + 0.5;

                        // The orientation for the new grain will depend on its Grain ID
                        int MyOrientation = GrainOrientation(((abs(MyGrainID) - 1) % NGrainOrientations));
                        // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
                        // (xp,yp,zp) is the new cell's center on the global grid
                        double xp = GlobalX + 0.5;
                        double yp = GlobalY + 0.5;
                        double zp = GlobalZ + 0.5;
                        
                        float cx = DOCenter((long int)(3*D3D1ConvPosition));
                        float cy = DOCenter((long int)(3*D3D1ConvPosition+1));
                        float cz = DOCenter((long int)(3*D3D1ConvPosition+2));
                    
                        // Calculate critical diagonal lengths for the new active cell located at (xp,yp,zp) on the local grid
                        // For each neighbor (l=0 to 25), calculate which octahedron face leads to cell capture
                        // Calculate critical octahedron diagonal length to activate each nearest neighbor, as well as the coordinates of the triangle vertices on the capturing face
                         for (int n=0; n<26; n++)  {

                             // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the image of the center of a neighbor cell
                             double x0 = xp + NeighborX(n) - cx;
                             double y0 = yp + NeighborY(n) - cy;
                             double z0 = zp + NeighborZ(n) - cz;
                             // mag0 is the magnitude of (x0,y0,z0)
                             double mag0 = pow(pow(x0,2.0) + pow(y0,2.0) + pow(z0,2.0),0.5);
                             
                             // Calculate unit vectors for the octahedron that intersect the new cell center
                             double Diag1X, Diag1Y, Diag1Z, Diag2X, Diag2Y, Diag2Z, Diag3X, Diag3Y, Diag3Z;
                             double Angle1 = (GrainUnitVector(9*MyOrientation)*x0 + GrainUnitVector(9*MyOrientation + 1)*y0 + GrainUnitVector(9*MyOrientation + 2)*z0)/mag0;
                             if (Angle1 < 0) {
                                 Diag1X = GrainUnitVector(9*MyOrientation);
                                 Diag1Y = GrainUnitVector(9*MyOrientation + 1);
                                 Diag1Z = GrainUnitVector(9*MyOrientation + 2);
                             }
                             else {
                                 Diag1X = -GrainUnitVector(9*MyOrientation);
                                 Diag1Y = -GrainUnitVector(9*MyOrientation + 1);
                                 Diag1Z = -GrainUnitVector(9*MyOrientation + 2);
                             }
                             double Angle2 = (GrainUnitVector(9*MyOrientation + 3)*x0 + GrainUnitVector(9*MyOrientation + 4)*y0 + GrainUnitVector(9*MyOrientation + 5)*z0)/mag0;
                             if (Angle2 < 0) {
                                 Diag2X = GrainUnitVector(9*MyOrientation + 3);
                                 Diag2Y = GrainUnitVector(9*MyOrientation + 4);
                                 Diag2Z = GrainUnitVector(9*MyOrientation + 5);
                             }
                             else {
                                 Diag2X = -GrainUnitVector(9*MyOrientation + 3);
                                 Diag2Y = -GrainUnitVector(9*MyOrientation + 4);
                                 Diag2Z = -GrainUnitVector(9*MyOrientation + 5);
                             }

                             double Angle3 = (GrainUnitVector(9*MyOrientation + 6)*x0 + GrainUnitVector(9*MyOrientation + 7)*y0 + GrainUnitVector(9*MyOrientation + 8)*z0)/mag0;
                             if (Angle3 < 0) {
                                 Diag3X = GrainUnitVector(9*MyOrientation + 6);
                                 Diag3Y = GrainUnitVector(9*MyOrientation + 7);
                                 Diag3Z = GrainUnitVector(9*MyOrientation + 8);
                             }
                             else {
                                 Diag3X = -GrainUnitVector(9*MyOrientation + 6);
                                 Diag3Y = -GrainUnitVector(9*MyOrientation + 7);
                                 Diag3Z = -GrainUnitVector(9*MyOrientation + 8);
                             }

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
                             float CDLVal = pow(pow(ParaT*Diag1X,2.0) + pow(ParaT*Diag1Y,2.0) + pow(ParaT*Diag1Z,2.0),0.5);
                             //                                if ((normx*Diag1X+normy*Diag1Y+normz*Diag1Z) == 0.0) {
                             //                                    printf("Captured cell : %d %d %d %f %d %d %d %f %f %f",MyNeighborX,MyNeighborY,MyNeighborZ,mag0,index1,index2,index3,normx,normy,normz);
                             //                                }
                             CritDiagonalLength((long int)(26)*D3D1ConvPosition+(long int)(n)) = CDLVal;
                        }
                    }
                }
                else if (CellType(D3D1ConvPositionGlobal) == LiqSol) {
                    // Mark and count the number of nucleation events to be sent to other ranks
                    GrainID(D3D1ConvPositionGlobal) = NCounter;
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
                        int l = ItList(ItBounds,ll);
                        // Local coordinates of adjacent cell center
                        int MyNeighborX = RankX + NeighborX(l);
                        int MyNeighborY = RankY + NeighborY(l);
                        int MyNeighborZ = RankZ + NeighborZ(l);
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
                //if (RankZ == nz-2) cout << " ID = " << id << " Melted = " << Melted[RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY] << endl;
            }
        }
    }
    
}


// After initializing grain structure and filling ghost nodes, the known potential nucleation sites are placed into the nucleation data structures
// Each nucleation event is assigned a time step, beyond which if the associated cell is not solid or actve, the event occurs
// This data is synced across MPI ranks, for nucleation events that occur in the ghost nodes
void NucleiInit(int DecompositionStrategy, int MyXSlices, int MyYSlices, int nz, int id, double dTN, double dTsigma, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, ViewI_H NucleiLocation, ViewI_H NucleationTimes, ViewI_H CellType, ViewI_H GrainID, ViewI_H CritTimeStep, ViewF_H UndercoolingChange) {

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
        distribution(generator);
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
            int RankZ = i/(MyXSlices*MyYSlices);
            int Rem = i % (MyXSlices*MyYSlices);
            int RankX = Rem/MyYSlices;
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
    ViewI_H GhostNodesAR("bufferAR", 4*ARCount);
    ViewI_H GhostNodesBR("bufferBR", 4*BRCount);
    ViewI_H GhostNodesCR("bufferCR", 4*CRCount);
    ViewI_H GhostNodesDR("bufferDR", 4*DRCount);
    ViewI_H GhostNodesER("bufferER", 3*ERCount);
    ViewI_H GhostNodesFR("bufferFR", 3*FRCount);
    ViewI_H GhostNodesGR("bufferGR", 3*GRCount);
    ViewI_H GhostNodesHR("bufferHR", 3*HRCount);
    
    //MPI_Barrier(MPI_COMM_WORLD);
    //cout << "ID = " << id << " A to send " << ACount << " B to send " << BCount << " A to recieve " << ARCount << " B to recieve " << BRCount << endl;
    // Collect ghost node data and send to other ranks- left and right
    if (ACount > 0) {
        ViewI_H GhostNodesA("bufferA", 4*ACount);
        for (int i=0; i<4*ACount; i++) {
            GhostNodesA(i) = ANuc[i];
            //cout << "Value " << i << " from rank " << id << " is " << GhostNodesA[i] << endl;
        }
        if (BRCount == 0) {
            // Sending data to id = id - 1 only
	    MPI_Send(GhostNodesA.data(),ACount*4,MPI_INT,MyLeft,0,MPI_COMM_WORLD);
        }
        else {
            // Sending data to id = id - 1 and recieving data from id = id + 1
	    MPI_Sendrecv(GhostNodesA.data(),ACount*4,MPI_INT,MyLeft,0,GhostNodesBR.data(),BRCount*4,MPI_INT,MyRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        // cout << "ID = " << id << " ACount = " << ACount << " ABuf = " << BufACount << endl;
    }
    else if (BRCount > 0) {
        // Recieving data from id = id + 1 only
        MPI_Recv(GhostNodesBR.data(),BRCount*4,MPI_INT,MyRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    //cout << "Collect B Start" << endl;
   // if (id == 0) cout << " BCOUNT = " << BCount << endl;
    if (BCount > 0) {
        ViewI_H GhostNodesB("bufferB", 4*BCount);
        for (int i=0; i<4*BCount; i++) {
            GhostNodesB(i) = BNuc[i];
        }
        if (ARCount == 0) {
            // Sending data to id = id + 1 only
            MPI_Send(GhostNodesB.data(),BCount*4,MPI_INT,MyRight,1,MPI_COMM_WORLD);
            //if (id == 0) cout << " Rank 0 sent data starting with " << GhostNodesB[0] << endl;
        }
        else {
            // Sending data to id = id + 1 and recieving data from id = id - 1
	    MPI_Sendrecv(GhostNodesB.data(),BCount*4,MPI_INT,MyRight,1,GhostNodesAR.data(),ARCount*4,MPI_INT,MyLeft,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            //if (id == 1) cout << " Rank 1 recieved data starting with " << GhostNodesAR[0] << endl;
        }
    }
    else if (ARCount > 0) {
        // Recieving data from id = id - 1 only
        MPI_Recv(GhostNodesAR.data(),ARCount*4,MPI_INT,MyLeft,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        //if (id == 1) cout << " Rank 1 recieved data starting with " << GhostNodesAR[0] << endl;
    }


    if (DecompositionStrategy != 1) {
        // Collect ghost node data and send to other ranks- in and out
        if (CCount > 0) {
	    ViewI_H GhostNodesC("bufferC", 4*CCount);
            for (int i=0; i<4*CCount; i++) {
                GhostNodesC(i) = CNuc[i];
            }
            
            if (DRCount == 0) {
                // Sending data only
	        MPI_Send(GhostNodesC.data(),CCount*4,MPI_INT,MyIn,0,MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
                MPI_Sendrecv(GhostNodesC.data(),CCount*4,MPI_INT,MyIn,0,GhostNodesDR.data(),DRCount*4,MPI_INT,MyOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        else if (DRCount > 0) {
            // Recieving data only
	    MPI_Recv(GhostNodesDR.data(),DRCount*4,MPI_INT,MyOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        
        if (DCount > 0) {
	    ViewI_H GhostNodesD("bufferD", 4*DCount);
            for (int i=0; i<4*DCount; i++) {
                GhostNodesD(i) = DNuc[i];
            }
            if (CRCount == 0) {
                // Sending data only
                MPI_Send(GhostNodesD.data(),DCount*4,MPI_INT,MyOut,1,MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
	        MPI_Sendrecv(GhostNodesD.data(),DCount*4,MPI_INT,MyOut,1,GhostNodesCR.data(),CRCount*4,MPI_INT,MyIn,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            //    cout << "ID = " << id << " DCount = " << DCount << " DBuf = " << BufDCount << endl;
        }
        else if (CRCount > 0) {
            // Recieving data only
	    MPI_Recv(GhostNodesCR.data(),CRCount*4,MPI_INT,MyIn,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        
        // Collect ghost node data and send to other ranks- MyLeftIn and MyRightOut
        if (ECount > 0) {
	    ViewI_H GhostNodesE("bufferE", 3*ECount);
            for (int i=0; i<3*ECount; i++) {
                GhostNodesE(i) = ENuc[i];
            }
            if (HRCount == 0) {
                // Sending data only
                MPI_Send(GhostNodesE.data(),ECount*3,MPI_INT,MyLeftIn,0,MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
	        MPI_Sendrecv(GhostNodesE.data(),ECount*3,MPI_INT,MyLeftIn,0,GhostNodesHR.data(),HRCount*3,MPI_INT,MyRightOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            //  cout << "ID = " << id << " ECount = " << ECount << " EBuf = " << BufECount << endl;
        }
        else if (HRCount > 0) {
            // Recieving data only
	    MPI_Recv(GhostNodesHR.data(),HRCount*3,MPI_INT,MyRightOut,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        
        if (HCount > 0) {
	    ViewI_H GhostNodesH("bufferH", 3*HCount);
            for (int i=0; i<3*HCount; i++) {
                GhostNodesH(i) = HNuc[i];
            }
            if (ERCount == 0) {
                // Sending data only
                MPI_Send(GhostNodesH.data(),HCount*3,MPI_INT,MyRightOut,0,MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
	        MPI_Sendrecv(GhostNodesH.data(),HCount*3,MPI_INT,MyRightOut,0,GhostNodesER.data(),ERCount*3,MPI_INT,MyLeftIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        else if (ERCount > 0) {
            // Recieving data only
	    MPI_Recv(GhostNodesER.data(),ERCount*3,MPI_INT,MyLeftIn,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        
        // Collect ghost node data and send to other ranks- MyRightIn and MyLeftOut
        if (FCount > 0) {
	    ViewI_H GhostNodesF("bufferF", 4*FCount);
            for (int i=0; i<4*FCount; i++) {
                GhostNodesF(i) = FNuc[i];
            }
            if (GRCount == 0) {
                // Sending data only
                MPI_Send(GhostNodesF.data(),FCount*3,MPI_INT,MyRightIn,1,MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
	        MPI_Sendrecv(GhostNodesF.data(),FCount*3,MPI_INT,MyRightIn,1,GhostNodesGR.data(),GRCount*3,MPI_INT,MyLeftOut,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        else if (GRCount > 0) {
            // Recieving data only
	    MPI_Recv(GhostNodesGR.data(),GRCount*3,MPI_INT,MyLeftOut,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        
        if (GCount > 0) {
	    ViewI_H GhostNodesG("bufferG", 3*GCount);
            for (int i=0; i<3*GCount; i++) {
                GhostNodesG(i) = GNuc[i];
            }
            if (FRCount == 0) {
                // Sending data only
                MPI_Send(GhostNodesG.data(),GCount*3,MPI_INT,MyLeftOut,1,MPI_COMM_WORLD);
            }
            else {
                // Sending data and recieving data
	        MPI_Sendrecv(GhostNodesG.data(),GCount*3,MPI_INT,MyLeftOut,1,GhostNodesFR.data(),FRCount*3,MPI_INT,MyRightIn,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        else if (FRCount > 0) {
            // Recieving data only
	    MPI_Recv(GhostNodesFR.data(),FRCount*3,MPI_INT,MyRightIn,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    }
   // MPI_Barrier(MPI_COMM_WORLD);
   // cout << "Collect Done" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    // Place ghost node data recieved from the left (if needed)
    if (ARCount > 0) {
        for (int i=0; i<ARCount; i++) {
            int RankX = GhostNodesAR(4*i);
            int RankY = 0;
            int RankZ = GhostNodesAR(4*i+1);
            int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
            NucleiLocation(NEvent) = CellLocation;
            NucleationTimes(NEvent) = GhostNodesAR(4*i+2);
            CellType(CellLocation) = LiqSol;
            GrainID(CellLocation) = GhostNodesAR(4*i+3);
            //            cout << "GN Nuc ID " << id << " Cell " << CellLocation << " Time " << GhostNodesAR[4*i+2] << " GID " << GhostNodesAR[4*i+3] << endl;
            NEvent++;
        }
    }

    // Place ghost node data recieved from the right (if needed)
    if (BRCount > 0) {
        for (int i=0; i<BRCount; i++) {
            int RankX = GhostNodesBR(4*i);
            int RankY = MyYSlices-1;
            int RankZ = GhostNodesBR(4*i+1);
            int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
            NucleiLocation(NEvent) = CellLocation;
            NucleationTimes(NEvent) = GhostNodesBR(4*i+2);
            CellType(CellLocation) = LiqSol;
            GrainID(CellLocation) = GhostNodesBR(4*i+3);
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
                int RankY = GhostNodesCR(4*i);
                int RankZ = GhostNodesCR(4*i+1);
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                NucleiLocation(NEvent) = CellLocation;
                NucleationTimes(NEvent) = GhostNodesCR(4*i+2);
                CellType(CellLocation) = LiqSol;
                GrainID(CellLocation) = GhostNodesCR(4*i+3);
                //            cout << "GN Nuc ID " << id << " Cell " << CellLocation << " Time " << GhostNodesAR[4*i+2] << " GID " << GhostNodesAR[4*i+3] << endl;
                NEvent++;
            }
        }

        // Place ghost node data recieved from out of plane (if needed)
        if (DRCount > 0) {
            for (int i=0; i<DRCount; i++) {
                int RankX = 0;
                int RankY = GhostNodesDR(4*i);
                int RankZ = GhostNodesDR(4*i+1);
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                NucleiLocation(NEvent) = CellLocation;
                NucleationTimes(NEvent) = GhostNodesDR(4*i+2);
                CellType(CellLocation) = LiqSol;
                GrainID(CellLocation) = GhostNodesDR(4*i+3);
                      //      cout << "GN Nuc ID " << id << " Cell " << CellLocation << " Time " << GhostNodesAR[4*i+2] << " GID " << GhostNodesAR[4*i+3] << endl;
                NEvent++;
            }
        }

        // Place ghost node data recieved from left and out of plane (if needed)
        if (ERCount > 0) {
            for (int i=0; i<ERCount; i++) {
                int RankX = MyXSlices-1;
                int RankY = 0;
                int RankZ = GhostNodesER(3*i);
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                NucleiLocation(NEvent) = CellLocation;
                NucleationTimes(NEvent) = GhostNodesER(3*i+1);
                CellType(CellLocation) = LiqSol;
                GrainID(CellLocation) = GhostNodesER(3*i+2);
                   //         cout << "GN Nuc ID " << id << " Cell " << CellLocation << " Time " << GhostNodesAR[4*i+2] << " GID " << GhostNodesAR[4*i+3] << endl;
                NEvent++;
            }
        }

        // Place ghost node data recieved from right and out of plane (if needed)
        if (FRCount > 0) {
            for (int i=0; i<FRCount; i++) {
                int RankX = MyXSlices-1;
                int RankY = MyYSlices-1;
                int RankZ = GhostNodesFR(3*i);
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                NucleiLocation(NEvent) = CellLocation;
                NucleationTimes(NEvent) = GhostNodesFR(3*i+1);
                CellType(CellLocation) = LiqSol;
                GrainID(CellLocation) = GhostNodesFR(3*i+2);
                  //          cout << "GN Nuc ID " << id << " Cell " << CellLocation << " Time " << GhostNodesAR[4*i+2] << " GID " << GhostNodesAR[4*i+3] << endl;
                NEvent++;
            }
        }

        // Place ghost node data recieved from left and into plane (if needed)
        if (GRCount > 0) {
            for (int i=0; i<GRCount; i++) {
                int RankX = 0;
                int RankY = 0;
                int RankZ = GhostNodesGR(3*i);
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                NucleiLocation(NEvent) = CellLocation;
                NucleationTimes(NEvent) = GhostNodesGR(3*i+1);
                CellType(CellLocation) = LiqSol;
                GrainID(CellLocation) = GhostNodesGR(3*i+2);
                    //        cout << "GN Nuc ID " << id << " Cell " << CellLocation << " Time " << GhostNodesAR[4*i+2] << " GID " << GhostNodesAR[4*i+3] << endl;
                NEvent++;
            }
        }

        if (HRCount > 0) {
            for (int i=0; i<HRCount; i++) {
                int RankX = 0;
                int RankY = MyYSlices-1;
                int RankZ = GhostNodesHR(3*i);
                int CellLocation = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                NucleiLocation(NEvent) = CellLocation;
                NucleationTimes(NEvent) = GhostNodesHR(3*i+1);
                CellType(CellLocation) = LiqSol;
                GrainID(CellLocation) = GhostNodesHR(3*i+2);
                    //        cout << "GN Nuc ID " << id << " Cell " << CellLocation << " Time " << GhostNodesAR[4*i+2] << " GID " << GhostNodesAR[4*i+3] << endl;
                NEvent++;
            }
        }
    }
    cout << "(" << id << ": " << NEvent << ") " << flush;
}


void DomainShiftAndResize(int id, int MyXSlices, int MyYSlices, int &ZShift, int &ZBound_Low, int &ZBound_High, int &nzActive, int LocalDomainSize, int &LocalActiveDomainSize, int &BufSizeZ, int LayerHeight, ViewI CellType, int layernumber, ViewI LayerID) {
    
    int ZBound_LowOld = ZBound_Low;
    
    // The top "top" of the active domain is a shift of "LayerHeight" from the previous domain top
    ZBound_High += LayerHeight;
    
    // The new "bottom" of the active domain is located just below the lowest active cells remaining in the domain
    int NewMin;
    Kokkos::parallel_reduce("MinReduce", LocalDomainSize, KOKKOS_LAMBDA (const int& D3D1ConvPosition, int& lmin) {
        if (CellType(D3D1ConvPosition) == Active) {
            // Check Z position of this active cell
            int RankZ = D3D1ConvPosition/(MyXSlices*MyYSlices);
            if (RankZ < lmin) lmin = RankZ;
        }
    }, Kokkos::Min<int>(NewMin));
    NewMin--;
    MPI_Allreduce(&NewMin,&ZBound_Low,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    //if (layernumber == 2) cout << "On rank " << id << " the new domain bottom should be at " << NewMin+1 << endl;
    // Shift in +Z direction for the bottom of the active region
    ZShift = ZBound_Low - ZBound_LowOld;
    
    if (id == 0) cout << "New domain bottom at Z = " << ZBound_Low << ", a shift of " << ZShift << " cells" << endl;
    
    // The new "top" of the active domain is located at the highest location with cells solidifying during the next layer
    int NewMax;
    Kokkos::parallel_reduce("MaxReduce", LocalDomainSize, KOKKOS_LAMBDA (const int& D3D1ConvPosition, int& lmax) {
        if (LayerID(D3D1ConvPosition) == layernumber+1) {
            // Check Z position of this active cell
            int RankZ = D3D1ConvPosition/(MyXSlices*MyYSlices);
            if (RankZ > lmax) lmax = RankZ;
        }
    }, Kokkos::Max<int>(NewMax));
    MPI_Allreduce(&NewMax,&ZBound_High,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    
    // Change in active region data structures' sizes
    nzActive = ZBound_High - ZBound_Low + 1;
    LocalActiveDomainSize = MyXSlices*MyYSlices*nzActive;
    
    // Change in height of buffers
    BufSizeZ = nzActive;
    

    if (id == 0) cout << "New active domain height is " << nzActive << endl;
    
}

void LayerSetup(int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int LocalActiveDomainSize, ViewI GrainOrientation, int NGrainOrientations, ViewF GrainUnitVector, ViewI NeighborX, ViewI NeighborY, ViewI NeighborZ, ViewF DiagonalLength, ViewI CellType, ViewI GrainID, ViewF CritDiagonalLength, ViewF DOCenter, Buffer2D BufferA, Buffer2D BufferB, Buffer2D BufferC, Buffer2D BufferD, Buffer2D BufferE, Buffer2D BufferF, Buffer2D BufferG, Buffer2D BufferH, Buffer2D BufferAR, Buffer2D BufferBR, Buffer2D BufferCR, Buffer2D BufferDR, Buffer2D BufferER, Buffer2D BufferFR, Buffer2D BufferGR, Buffer2D BufferHR, int BufSizeX, int BufSizeY, int BufSizeZ, int &ZBound_Low, ViewI Locks) {

    // Reset active cell data structures
    Kokkos::parallel_for("DelDLData",LocalActiveDomainSize, KOKKOS_LAMBDA (const int& i) {
       DiagonalLength(i) = 0;
    });

    Kokkos::parallel_for("DelDOData",3*LocalActiveDomainSize, KOKKOS_LAMBDA (const int& i) {
        DOCenter(i) = 0;
    });

    Kokkos::parallel_for("DelCDLData",26*LocalActiveDomainSize, KOKKOS_LAMBDA (const int& i) {
        CritDiagonalLength(i) = 0;
    });
//    Kokkos::fence();
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (id == 0) cout << "X Act structures remade" << endl;
    
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
//    Kokkos::fence();
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (id == 0) cout << "X Buffers remade" << endl;
    
    Kokkos::parallel_for("NewActiveCellInit",LocalActiveDomainSize, KOKKOS_LAMBDA (const int& D3D1ConvPosition) {
        // Initialize active cell data structures for those that are now part of the active domain
        int RankZ = D3D1ConvPosition/(MyXSlices*MyYSlices);
        int Rem = D3D1ConvPosition % (MyXSlices*MyYSlices);
        int RankX = Rem/MyYSlices;
        int RankY = Rem % MyYSlices;
        int GlobalZ = RankZ + ZBound_Low;
        int GlobalD3D1ConvPosition = GlobalZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
        if (CellType(GlobalD3D1ConvPosition) == Active) {
            
            int GlobalX = RankX + MyXOffset;
            int GlobalY = RankY + MyYOffset;

            int MyGrainID = GrainID(GlobalD3D1ConvPosition);
            DiagonalLength(D3D1ConvPosition) = 0.01;
            DOCenter(3*D3D1ConvPosition) = GlobalX + 0.5;
            DOCenter(3*D3D1ConvPosition+1) = GlobalY + 0.5;
            DOCenter(3*D3D1ConvPosition+2) = GlobalZ + 0.5;

            // The orientation for the new grain will depend on its Grain ID
            int MyOrientation = GrainOrientation(((abs(MyGrainID) - 1) % NGrainOrientations));
            // Calculate critical values at which this active cell leads to the activation of a neighboring liquid cell
            // (xp,yp,zp) is the new cell's center on the global grid
            double xp = GlobalX + 0.5;
            double yp = GlobalY + 0.5;
            double zp = GlobalZ + 0.5;
            
            float cx = DOCenter((long int)(3*D3D1ConvPosition));
            float cy = DOCenter((long int)(3*D3D1ConvPosition+1));
            float cz = DOCenter((long int)(3*D3D1ConvPosition+2));

            // Calculate critical diagonal lengths for the new active cell located at (xp,yp,zp) on the local grid
            // For each neighbor (l=0 to 25), calculate which octahedron face leads to cell capture
            // Calculate critical octahedron diagonal length to activate each nearest neighbor, as well as the coordinates of the triangle vertices on the capturing face
            for (int n=0; n<26; n++)  {

                // (x0,y0,z0) is a vector pointing from this decentered octahedron center to the image of the center of a neighbor cell
                double x0 = xp + NeighborX(n) - cx;
                double y0 = yp + NeighborY(n) - cy;
                double z0 = zp + NeighborZ(n) - cz;
                // mag0 is the magnitude of (x0,y0,z0)
                double mag0 = pow(pow(x0,2.0) + pow(y0,2.0) + pow(z0,2.0),0.5);
                
                // Calculate unit vectors for the octahedron that intersect the new cell center
                double Diag1X, Diag1Y, Diag1Z, Diag2X, Diag2Y, Diag2Z, Diag3X, Diag3Y, Diag3Z;
                double Angle1 = (GrainUnitVector(9*MyOrientation)*x0 + GrainUnitVector(9*MyOrientation + 1)*y0 + GrainUnitVector(9*MyOrientation + 2)*z0)/mag0;
                if (Angle1 < 0) {
                    Diag1X = GrainUnitVector(9*MyOrientation);
                    Diag1Y = GrainUnitVector(9*MyOrientation + 1);
                    Diag1Z = GrainUnitVector(9*MyOrientation + 2);
                }
                else {
                    Diag1X = -GrainUnitVector(9*MyOrientation);
                    Diag1Y = -GrainUnitVector(9*MyOrientation + 1);
                    Diag1Z = -GrainUnitVector(9*MyOrientation + 2);
                }
                double Angle2 = (GrainUnitVector(9*MyOrientation + 3)*x0 + GrainUnitVector(9*MyOrientation + 4)*y0 + GrainUnitVector(9*MyOrientation + 5)*z0)/mag0;
                if (Angle2 < 0) {
                    Diag2X = GrainUnitVector(9*MyOrientation + 3);
                    Diag2Y = GrainUnitVector(9*MyOrientation + 4);
                    Diag2Z = GrainUnitVector(9*MyOrientation + 5);
                }
                else {
                    Diag2X = -GrainUnitVector(9*MyOrientation + 3);
                    Diag2Y = -GrainUnitVector(9*MyOrientation + 4);
                    Diag2Z = -GrainUnitVector(9*MyOrientation + 5);
                }

                double Angle3 = (GrainUnitVector(9*MyOrientation + 6)*x0 + GrainUnitVector(9*MyOrientation + 7)*y0 + GrainUnitVector(9*MyOrientation + 8)*z0)/mag0;
                if (Angle3 < 0) {
                    Diag3X = GrainUnitVector(9*MyOrientation + 6);
                    Diag3Y = GrainUnitVector(9*MyOrientation + 7);
                    Diag3Z = GrainUnitVector(9*MyOrientation + 8);
                }
                else {
                    Diag3X = -GrainUnitVector(9*MyOrientation + 6);
                    Diag3Y = -GrainUnitVector(9*MyOrientation + 7);
                    Diag3Z = -GrainUnitVector(9*MyOrientation + 8);
                }

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
                float CDLVal = pow(pow(ParaT*Diag1X,2.0) + pow(ParaT*Diag1Y,2.0) + pow(ParaT*Diag1Z,2.0),0.5);
                //                                if ((normx*Diag1X+normy*Diag1Y+normz*Diag1Z) == 0.0) {
                //                                    printf("Captured cell : %d %d %d %f %d %d %d %f %f %f",MyNeighborX,MyNeighborY,MyNeighborZ,mag0,index1,index2,index3,normx,normy,normz);
                //                                }
                CritDiagonalLength((long int)(26)*D3D1ConvPosition+(long int)(n)) = CDLVal;
                //                                if (CDLVal == 0.0) printf("Zero CDLVal : %d %d %d %d %f %d %d %d %f %f %f",MyNeighborX,MyNeighborY,MyNeighborZ,n,mag0,index1,index2,index3,normx,normy,normz);
            }
        }
    });
//    Kokkos::fence();
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (id == 0) cout << "X new active cells remade" << endl;
    
    // Reset lock values
    Kokkos::parallel_for("LockInit",LocalActiveDomainSize, KOKKOS_LAMBDA (const int& D3D1ConvPosition) {
        int RankZ = D3D1ConvPosition/(MyXSlices*MyYSlices);
        int Rem = D3D1ConvPosition % (MyXSlices*MyYSlices);
        int RankX = Rem/MyYSlices;
        int RankY = Rem % MyYSlices;
        int GlobalZ = ZBound_Low + RankZ;
        int GlobalD3D1ConvPosition = GlobalZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
        if ((CellType(GlobalD3D1ConvPosition) == Delayed)||(CellType(GlobalD3D1ConvPosition) == LiqSol)||(CellType(GlobalD3D1ConvPosition) == Liquid)) Locks(D3D1ConvPosition) = 1;
        else Locks(D3D1ConvPosition) = 0;
    });
//    Kokkos::fence();
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (id == 0) cout << "X locks remade" << endl;
    
}

//void EraseTop(int NextLayerNumber, int TempFilesInSeries, int DecompositionStrategy, string SimulationType, int id, int np, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset, double deltax, double HT_deltax, double deltat, int &nx, int &ny, int &nz, int &ProcessorsInXDirection, int &ProcessorsInYDirection, ViewI_H GrainID, string tempfile, float XMin, float XMax, float YMin, float YMax, float ZMin, float ZMax, bool* Melted, string TemperatureDataSource, int LayersSimulatedAtOnce, int LayerHeight, int NumberOfLayers, ViewI_H CellType, ViewI_H CritTimeStep, ViewF_H UndercoolingChange) {
//    
//    // Temperature data read
//    int HTtoCAratio = HT_deltax/deltax; // OpenFOAM/CA cell size ratio
//    
//    // With row/col 0 being wall cells and row/col 1 being solid cells outside of the melted area, the domain starts at row/col 2
//    // As the wall cells are not part of the physical domain (solid cells at row/col 1 are defined as X = Y = 0,
//    // the melted region domain data starts at X = Y = deltax, with data points at X or Y = deltax + N*HT_deltax through X or Y = nx-3 or ny-3
//    
//    // The X and Y bounds are the region (for this MPI rank) of the physical domain that needs to be read
//    // Extends past the actual spatial extent of the local domain for purposes of interpolating from HT_deltax to deltax
//    int LowerXBound, LowerYBound, UpperXBound, UpperYBound;
//    if (MyXOffset <= 2) LowerXBound = 2;
//    else LowerXBound = MyXOffset - ((MyXOffset-2) % HTtoCAratio);
//    if (MyYOffset <= 2) LowerYBound = 2;
//    else LowerYBound = MyYOffset - ((MyYOffset-2) % HTtoCAratio);
//    
//    if (MyXOffset+MyXSlices-1 >= nx-3) UpperXBound = nx-3;
//    else UpperXBound = MyXOffset + MyXSlices - 1 + HTtoCAratio - ((MyXOffset + (MyXSlices-1) - 2) % HTtoCAratio);
//    if (MyYOffset+MyYSlices-1 >= ny-3) UpperYBound = ny-3;
//    else UpperYBound = MyYOffset + MyYSlices - 1 + HTtoCAratio - ((MyYOffset + (MyYSlices-1) - 2) % HTtoCAratio);
//    
//    vector <vector <vector <int> > > Data;
//    for (int k=0; k<nz; k++) {
//        vector <vector <int> > TemperatureXX;
//        for (int i=LowerXBound; i<=UpperXBound; i++) {
//            vector <int> TemperatureX;
//            for (int j=LowerYBound; j<=UpperYBound; j++) {
//                TemperatureX.push_back(0);
//            }
//            TemperatureXX.push_back(TemperatureX);
//        }
//        Data.push_back(TemperatureXX);
//    }
//    
//    double UnitConversion;
//    if (TemperatureDataSource == "O") UnitConversion = 1;
//    else UnitConversion = 1000;
//    
//    ifstream Geom;
//    string thislayertempfile;
//    
//    if (LayersSimulatedAtOnce > 1) {
//        
//        // "layernumber" is the layer that would be depositied next
//        int NextLayerFile = NumberOfLayers+1;
//        string NextLayerFileS;
//        NextLayerFile = NextLayerFile % TempFilesInSeries;
//        if (NextLayerFile == 0) NextLayerFile = TempFilesInSeries;
//        NextLayerFileS = to_string(NextLayerFile);
//        thislayertempfile = NextLayerFileS + tempfile;
//        ifstream TemperatureFileRead;
//        TemperatureFileRead.open(thislayertempfile);
//        if (!TemperatureFileRead.good()) {
//            if (id == 0) cout << "Temperature file " << thislayertempfile << " not found" << endl;
//            TemperatureFileRead.close();
//        }
//    }
//    else thislayertempfile = tempfile;
//    
//    if (id == 0) cout << "Layer " << NumberOfLayers+1 << " temperature file is: " << thislayertempfile << endl;
//    double PseudoLayerOffsetZ;
//    if (LayersSimulatedAtOnce > 1) {
//        PseudoLayerOffsetZ = NumberOfLayers*LayerHeight*deltax;
//    }
//    else {
//        PseudoLayerOffsetZ = LayerHeight*deltax;
//    }
//
//    Geom.open(thislayertempfile);
//    if (id == 0) cout << "Layer " << NumberOfLayers << " LayerHeight = " << LayerHeight << " Z Off = " << PseudoLayerOffsetZ/deltax << " cells" << endl;
//    while (!Geom.eof()) {
//        string s;
//        getline(Geom,s);
//        bool ReadingLine = true;
//        int i = 0;
//        int j = 0;
//        int XInt, YInt, ZInt;
//        double MyX, MyY, MyZ;
//        
//        int FirstChar = 0;
//        int LastChar;
//        while (ReadingLine) {
//            if (s.empty()) break;
//            char C = s.at(i);
//            // If this character is not a space, convert from string
//            if ((isblank(C))||(i == s.length()-1)) {
//                string NewDataS;
//                if (j < 4) {
//                    LastChar = i;
//                    NewDataS = s.substr(FirstChar,LastChar-FirstChar);
//                    if (j < 3) {
//                        bool Searching = true;
//                        while (Searching) {
//                            //cout << "N searching: Char " << i << " is " << s.at(i) << " of line " << s << endl;
//                            i++;
//                            C = s.at(i);
//                            if (!(isblank(C))) {
//                                FirstChar = i;
//                                Searching = false;
//                            }
//                        }
//                    }
//                }
//                else {
//                    i = s.length()-1;
//                    bool Searching = true;
//                    while (Searching) {
//                        //cout << "Space searching: Char " << i << " is " << s.at(i) << " of line " << s << endl;
//                        C = s.at(i);
//                        if (isblank(C)) {
//                            FirstChar = i;
//                            Searching = false;
//                        }
//                        i--;
//                    }
//                    //cout << "Chars " << FirstChar << " + " << s.length()-FirstChar << endl;
//                    NewDataS = s.substr(FirstChar,s.length()-FirstChar);
//                    //cout << NewDataS << endl;
//                }
//                if (j == 0) {
//                    MyX = atof(NewDataS.c_str())/UnitConversion; // Only divide by 1000 if from Truchas
//                    XInt = round((MyX-XMin)/deltax) + 2; // + HTtoCAratio;
//                    //if (id == 0) cout << XInt << endl;
//                }
//                else if (j == 1) {
//                    //if (id == 0) cout << "Y " << MeshData << endl;
//                    MyY = atof(NewDataS.c_str())/UnitConversion; // Only divide by 1000 if from Truchas
//                    YInt = round((MyY-YMin)/deltax) + 2;// + HTtoCAratio;
//                    //if (id == 0) cout << YInt << endl;
//                }
//                else if (j == 2) {
//                    //if (id == 0) cout << "Z " << MeshData << endl;
//                    MyZ = atof(NewDataS.c_str())/UnitConversion + PseudoLayerOffsetZ; // Only divide by 1000 if from Truchas
//                    ZInt = round((MyZ-ZMin)/deltax) + 2;
//                }
//                else if (j == 3) {
//                    // Determine if this X and Y fall in this rank's range of interest
//                    if ((XInt >= LowerXBound)&&(XInt <= UpperXBound)&&(YInt >= LowerYBound)&&(YInt <= UpperYBound)&&(ZInt <= nz-1)) {
//                        Data[ZInt][XInt-LowerXBound][YInt-LowerYBound] = 1;
//                    }
//                    
//                }
//                else if (j == 4) {
//                    ReadingLine = false;
//                }
//                j++;
//            }
//            // advance to the next character
//            i++;
//        }
//    }
//
//    // Data interpolation between heat transport and CA grids, if necessary
//    if (HTtoCAratio != 1) {
//        for (int k=2; k<=nz-2; k++) {
//            
//            int LowZ = k - ((k-2) % HTtoCAratio);
//            int HighZ = LowZ + HTtoCAratio;
//
//            if (HighZ > nz-2) HighZ = nz-2;
//            
//            for (int i=0; i<=UpperXBound-LowerXBound; i++) {
//                int LowX =  i - (i % HTtoCAratio);
//                int HighX = LowX + HTtoCAratio;
//
//                if (HighX >= UpperXBound-LowerXBound) HighX = UpperXBound-LowerXBound;
//                
//                for (int j=0; j<=UpperYBound-LowerYBound; j++) {
//                    int LowY = j - (j % HTtoCAratio);
//                    int HighY = LowY + HTtoCAratio;
//
//                    if (HighY >= UpperYBound-LowerYBound) HighY = UpperYBound-LowerYBound;
//
//                    double Pt1 = Data[LowZ][LowX][LowY];
//                    double Pt2 = Data[LowZ][HighX][LowY];
//                    double Pt3 = Data[LowZ][LowX][HighY];
//                    double Pt4 = Data[LowZ][HighX][HighY];
//                    double Pt5 = Data[HighZ][LowX][LowY];
//                    double Pt6 = Data[HighZ][HighX][LowY];
//                    double Pt7 = Data[HighZ][LowX][HighY];
//                    double Pt8 = Data[HighZ][HighX][HighY];
//                    if ((Pt1 > 0)&&(Pt2 > 0)&&(Pt3 > 0)&&(Pt4 > 0)&&(Pt5 > 0)&&(Pt6 > 0)&&(Pt7 > 0)&&(Pt8 > 0)) {
//                        Data[k][i][j] = 1;
//                    }
//                    
//                }
//            }
//        }
//    }
//
//    int MaxCTS = 0;
//    int GMax;
//    for (int k=0; k<nz; k++) {
//        for (int ii=LowerXBound; ii<=UpperXBound; ii++) {
//            for (int jj=LowerYBound; jj<=UpperYBound; jj++) {
//                if ((ii >= MyXOffset)&&(ii < MyXOffset+MyXSlices)&&(jj >= MyYOffset)&&(jj < MyYOffset+MyYSlices)) {
//                    int Adj_i = ii - MyXOffset;
//                    int Adj_j = jj - MyYOffset;
//                    double CTLiq = Data[k][ii-LowerXBound][jj-LowerYBound];
//                    int Coord3D1D = k*MyXSlices*MyYSlices + Adj_i*MyYSlices + Adj_j;
//                    if (CTLiq > 0)  {
//                        GrainID(Coord3D1D) = -1;
//                        CellType(Coord3D1D) = Wall;
//                        CritTimeStep(Coord3D1D) = 0;
//                        UndercoolingChange(Coord3D1D) = 0;
//                    }
//                    if (CritTimeStep(Coord3D1D) > MaxCTS) MaxCTS = CritTimeStep(Coord3D1D);
//                }
//            }
//        }
//    }
//    MPI_Reduce(&MaxCTS, &GMax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
//    MPI_Bcast(&GMax,1,MPI_INT,0,MPI_COMM_WORLD);
//    if (id == 0) cout << "Largest critical time step after remelting is " << GMax << endl;
//}
