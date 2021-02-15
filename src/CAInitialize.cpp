#include "header.h"

#include <iostream>
#include <string>
#include <regex>

using namespace std;
// Initializes input parameters, mesh, temperature field, and grain structures for CA simulations

//*****************************************************************************/
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

    std::vector<std::string> deprecated_inputs = {"Burst buffer", "Source of input length unit" };
    for (auto di : deprecated_inputs)
    if ( actual_key.find( di ) != std::string::npos )
    {
        std::cout << "WARNING - this input has been deprecated and has no effect, \""
                  << actual_key << "\"" << std::endl;

        // Ignore this input and get another line
        std::getline(stream, line);
        colon = line.find(":");
        actual_key = line.substr(0, colon);
    }

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

//*****************************************************************************/
// Read ExaCA input file.
void InputReadFromFile(int id, string InputFile, string &SimulationType, int &DecompositionStrategy, double &AConst, double &BConst, double &CConst, double &DConst, double& FreezingRange, double &deltax, double &NMax, double &dTN, double &dTsigma, string &OutputFile, string &GrainOrientationFile, string &tempfile, int &TempFilesInSeries, bool &ExtraWalls, double &HT_deltax, bool &RemeltingYN, double &deltat, int &NumberOfLayers, int &LayerHeight, string &SubstrateFileName, double &G, double &R, int &nx, int &ny, int &nz, double &FractSurfaceSitesActive, string &PathToOutput, bool (&FilesToPrint)[6], bool &PrintFilesYN) {

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

        // Burst buffer/Truchas multilayer simulation input (no longer supported)
        
        // File containing temperature data
        tempfile = parseInput(InputData, "Temperature filename");
        tempfile = FilePath + "/Temperatures/" + tempfile;

        // Temperature files in series
        val = parseInput(InputData, "Number of temperature files");
        TempFilesInSeries = stoi(val,nullptr,10);
        
        // For now, assuming no remelting
        RemeltingYN = false;
        
        if (id == 0) {
            cout << "Temperature data file(s) is/are " << tempfile << " , and there are " << TempFilesInSeries << " in the series" << endl;
        }
        // Check to ensure all temperature files exist - obtain number of temperature data values and data units from each file
        for (int i=0; i<TempFilesInSeries; i++) {
            string CurrentTempFile = tempfile;
            if (TempFilesInSeries > 1) {
                int NextLayerFile = i + 1;
                string NextLayerFileS = to_string(NextLayerFile);
                std::size_t found = tempfile.find_last_of("/");
                string FPath = tempfile.substr(0,found+1);
                string FName = tempfile.substr(found+1);
                CurrentTempFile = FPath + NextLayerFileS + FName;
            }
            ifstream TemperatureFile;
            TemperatureFile.open(CurrentTempFile);
            if (TemperatureFile.is_open()) {
                if (id == 0) cout << "Successfully opened file" << endl;
                TemperatureFile.close();
            }
            else {
                if (id == 0) cout << "Failed to open file " << CurrentTempFile << endl;
                if (id == 0) throw std::runtime_error("Input \"Could not find or open temperature file(s)\"  \".");
            }
        }
    
        // Usage of second set of wall cells around temperature field (for spot melt problems, where the melt pool boundaries are right at the walls)
        ExtraWalls = parseInputBool(InputData, "Extra set of wall cells around temperature field");
        
        // Number of layers (for non script-based coupling)
        val = parseInput(InputData, "Number of layers");
        NumberOfLayers = stoi(val,nullptr,10);
        
        // Layer height (for non script-based coupling)
        val = parseInput(InputData, "Offset between layers");
        LayerHeight = stoi(val,nullptr,10);
        if (id == 0) cout << "A total of " << NumberOfLayers << " of solidification offset by " << LayerHeight << " CA cells will be simulated" << endl;

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

//*****************************************************************************/
void ParallelMeshInit(int DecompositionStrategy, ViewI_H NeighborX, ViewI_H NeighborY, ViewI_H NeighborZ, ViewI2D_H ItList, string SimulationType, int id, int np, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset,int &MyLeft, int &MyRight, int &MyIn, int &MyOut, int &MyLeftIn, int &MyLeftOut, int &MyRightIn, int &MyRightOut, double &deltax, double HT_deltax, int &nx, int &ny, int &nz, int &ProcessorsInXDirection, int &ProcessorsInYDirection, string tempfile, float &XMin, float &XMax, float &YMin, float &YMax, float &ZMin, float &ZMax, float FreezingRange, int &LayerHeight, int NumberOfLayers, int TempFilesInSeries, int &NumberOfTemperatureDataPoints, float* ZMinLayer, float* ZMaxLayer, int* FirstValue, int* LastValue, vector <float> &RawData) {
        
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
        // Two passes through reading temperature data files- the first pass ("Loop 0" only reads the headers to determine units and X/Y/Z bounds of the simulaton domain
        // Using the X/Y/Z bounds of the simulation domain, nx, ny, and nz can be calculated and the domain decomposed among MPI processes
        // The second pass ("Loop 1") reads the actual X/Y/Z/liquidus time/solidus time or X/Y/Z/liquidus time/cooling rate data and each rank stores the data relevant to itself in "RawData"
        // Note that if solidus time data is given, the liquidus/solidus times are used to calculate a cooling rate - solidus time values are not placced into "RawData"

        XMin = 1000000.0;
        YMin = 1000000.0;
        ZMin = 1000000.0;
        XMax = -1000000.0;
        YMax = -1000000.0;
        ZMax = -1000000.0;
        string LengthUnits, TimeUnits;
        bool UseCoolingRate;
        
        for (int Loop=0; Loop<=1; Loop++) {

            int LayersToRead = min(NumberOfLayers,TempFilesInSeries); // was given in input file
            for (int LayerReadCount=1; LayerReadCount<=LayersToRead; LayerReadCount++) {

                string tempfile_thislayer;
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
                ifstream TemperatureFile;
                TemperatureFile.open(tempfile_thislayer);
                
                if (Loop == 0) {
                    // Read the header line data
                    // First line is the number of temperature data points (currently not used by the CA code, but may be in the future)
                    string DummyVar_TemperatureDataPoints = parseInput(TemperatureFile, "temperature");
                    // Second line is the number of remelting events in the file (should be 0 as we currently are assuming no remelting)
                    ZMinLayer[LayerReadCount-1] = 1000000.0;
                    ZMaxLayer[LayerReadCount-1] = -1000000.0;
                    string val = parseInput(TemperatureFile, "remelting");
                    int RemeltingDummyVar = stoi(val,nullptr,10);
                    if ((id == 0)&&(RemeltingDummyVar != 0)) cout << "WARNING: ExaCA does not currently allow for multiple remelting events, remelting events was set to " << RemeltingDummyVar << " by design, errors will likely result" << endl;
                    // Third line contains the units of length (m or mm) and time (s or ms)
                    val = parseInput(TemperatureFile, "Units");
                    std::size_t findcommaseparator = val.find(",");
                    LengthUnits = val.substr(0,findcommaseparator);
                    TimeUnits = val.substr(findcommaseparator+1,val.length()-findcommaseparator);
                    std::regex r("\\s+");
                    LengthUnits = std::regex_replace(LengthUnits, r, "");
                    TimeUnits = std::regex_replace(TimeUnits, r, "");
                    if (((LengthUnits != "mm")&&(LengthUnits != "m"))||(TimeUnits != "ms")&&(TimeUnits != "s")) {
                        string error = "Length units \"" + LengthUnits + "\" or time units \"" + TimeUnits + "\" are not valid- must be mm or m for length and ms or s for time";
                        throw std::runtime_error( error );
                    }
                    if (id == 0) cout << "Units of length are " << LengthUnits << " and units of time are " << TimeUnits << endl;
                    // Fourth line contains XYZ bounds for the global temperature domain contained within this file
                    string GlobalBoundsLine;
                    getline(TemperatureFile,GlobalBoundsLine);
                    string XMinString, XMaxString, YMinString, YMaxString, ZMinString, ZMaxString;
                    size_t XStart = GlobalBoundsLine.find("[");
                    size_t XSeparator = GlobalBoundsLine.find(",");
                    size_t XEnd = GlobalBoundsLine.find(",",XSeparator+1);
                    XMinString = GlobalBoundsLine.substr(XStart+1,XSeparator-XStart-1);
                    XMaxString = GlobalBoundsLine.substr(XSeparator+1,XEnd-XSeparator-1);
                    size_t YSeparator = GlobalBoundsLine.find(",",XEnd+1);
                    size_t YEnd = GlobalBoundsLine.find(",",YSeparator+1);
                    YMinString = GlobalBoundsLine.substr(XEnd+1,YSeparator-XEnd-1);
                    YMaxString = GlobalBoundsLine.substr(YSeparator+1,YEnd-YSeparator-1);
                    size_t ZSeparator = GlobalBoundsLine.find(",",YEnd+1);
                    size_t ZEnd = GlobalBoundsLine.find("]");
                    ZMinString = GlobalBoundsLine.substr(YEnd+1,ZSeparator-YEnd-1);
                    ZMaxString = GlobalBoundsLine.substr(ZSeparator+1,ZEnd-ZSeparator-1);
                    if ((XStart == std::string::npos)||(XSeparator == std::string::npos)||(XEnd == std::string::npos)||(YSeparator == std::string::npos)||(YEnd == std::string::npos)||(ZSeparator == std::string::npos)||(ZEnd == std::string::npos)) {
                        string error = "XYZ data size line improperly formatted in input file: should be [xmin, xmax; ymin, ymax; zmin, zmax]";
                        throw std::runtime_error( error );
                    }
                    if ((Loop == 0)&&(id == 0)) cout << "Layer = " << LayerReadCount << " Z Bounds are " << XMinString << " " << XMaxString  << endl;
                    float XMin_ThisLayer = atof(XMinString.c_str());
                    float XMax_ThisLayer = atof(XMaxString.c_str());
                    float YMin_ThisLayer = atof(YMinString.c_str());
                    float YMax_ThisLayer = atof(YMaxString.c_str());
                    float ZMin_ThisLayer = atof(ZMinString.c_str());
                    float ZMax_ThisLayer = atof(ZMaxString.c_str());
                    if ((Loop == 0)&&(id == 0)) cout << "Layer = " << LayerReadCount << " Z Bounds are " << ZMin_ThisLayer << " " << ZMax_ThisLayer  << endl;
                    if (LengthUnits == "mm") {
                        XMin_ThisLayer = XMin_ThisLayer/1000.0;
                        XMax_ThisLayer = XMax_ThisLayer/1000.0;
                        YMin_ThisLayer = YMin_ThisLayer/1000.0;
                        YMax_ThisLayer = YMax_ThisLayer/1000.0;
                        ZMin_ThisLayer = ZMin_ThisLayer/1000.0;
                        ZMax_ThisLayer = ZMax_ThisLayer/1000.0;
                    }
                    // Make sure fifth line contains all required column names: x, y, z, tl, and either ts or cr
                    // Check mix of lower and uppercase letters just in case
                    string HeaderLine;
                    getline(TemperatureFile,HeaderLine);
                    std::size_t xf = HeaderLine.find("x");
                    if (xf == std::string::npos) xf = HeaderLine.find("X");
                    std::size_t yf = HeaderLine.find("y");
                    if (yf == std::string::npos) yf = HeaderLine.find("Y");
                    std::size_t zf = HeaderLine.find("z");
                    if (zf == std::string::npos) zf = HeaderLine.find("Z");
                    std::size_t col1 = HeaderLine.find("tl");
                    if (col1 == std::string::npos) col1 = HeaderLine.find("Tl");
                    if (col1 == std::string::npos) col1 = HeaderLine.find("TL");
                    if ((xf == std::string::npos)||(yf == std::string::npos)||(zf == std::string::npos)||(col1 == std::string::npos)) {
                        string error = "One of the first four required column headers did not appear in the temperature input file: these are x, y, z, tl";
                        throw std::runtime_error( error );
                    }
                    std::size_t col2at1 = HeaderLine.find("ts");
                    std::size_t col2at2 = HeaderLine.find("Ts");
                    std::size_t col2at3 = HeaderLine.find("TS");
                    std::size_t col2at4 = HeaderLine.find("cr");
                    std::size_t col2at5 = HeaderLine.find("CR");
                    if ((col2at1 == std::string::npos)&&(col2at2 == std::string::npos)&&(col2at3 == std::string::npos)&&(col2at4 == std::string::npos)&&(col2at5 == std::string::npos)) {
                        string error = "The fifth column header must be either ts or cr";
                        throw std::runtime_error( error );
                    }
                    else {
                        if ((col2at1 == std::string::npos)&&(col2at2 == std::string::npos)&&(col2at3 == std::string::npos)) UseCoolingRate = true;
                        else UseCoolingRate = false;
                    }
                    // Based on the input file's layer offset, adjust ZMin/ZMax from the temperature data coordinate system to the multilayer CA coordinate system
                    // Check to see in the XYZ bounds for this layer are also limiting for the entire multilayer CA coordinate system
                    ZMin_ThisLayer += deltax*LayerHeight*(LayerReadCount-1);
                    ZMax_ThisLayer += deltax*LayerHeight*(LayerReadCount-1);
                    if (XMin_ThisLayer < XMin) XMin = XMin_ThisLayer;
                    if (XMax_ThisLayer > XMax) XMax = XMax_ThisLayer;
                    if (YMin_ThisLayer < YMin) YMin = YMin_ThisLayer;
                    if (YMax_ThisLayer > YMax) YMax = YMax_ThisLayer;
                    if (ZMin_ThisLayer < ZMin) ZMin = ZMin_ThisLayer;
                    if (ZMax_ThisLayer > ZMax) ZMax = ZMax_ThisLayer;
                    ZMinLayer[LayerReadCount-1] = ZMin_ThisLayer;
                    ZMaxLayer[LayerReadCount-1] = ZMax_ThisLayer;
                    if (id == 0) cout << "Layer = " << LayerReadCount << " Z Bounds are " << ZMin_ThisLayer << " " << ZMax_ThisLayer  << endl;
                }
                else {
                    // Store raw data relevant to each rank in the vector structure RawData
                    // This additional section of code will be obsolete once the X/Y/Z domain bounds are included in the headers of files
                    
                    // With row/col 0 being wall cells and row/col 1 being solid cells outside of the melted area, the domain starts at row/col 2
                    // As the wall cells are not part of the physical domain (solid cells at row/col 1 are defined as X = Y = 0,
                    // the melted region domain data starts at X = Y = deltax, with data points at X or Y = deltax + N*HT_deltax through X or Y = nx-3 or ny-3
                    
                    // The X and Y bounds are the region (for this MPI rank) of the physical domain that needs to be read
                    // Extends past the actual spatial extent of the local domain for purposes of interpolating from HT_deltax to deltax
                    FirstValue[LayerReadCount-1] = NumberOfTemperatureDataPoints;
                    int HTtoCAratio = HT_deltax/deltax; // OpenFOAM/CA cell size ratio
                    int LowerXBound, LowerYBound, UpperXBound, UpperYBound;
                    if (MyXOffset <= 2) LowerXBound = 2;
                    else LowerXBound = MyXOffset - ((MyXOffset-2) % HTtoCAratio);
                    if (MyYOffset <= 2) LowerYBound = 2;
                    else LowerYBound = MyYOffset - ((MyYOffset-2) % HTtoCAratio);
                    
                    if (MyXOffset+MyXSlices-1 >= nx-3) UpperXBound = nx-3;
                    else UpperXBound = MyXOffset + MyXSlices - 1 + HTtoCAratio - ((MyXOffset + (MyXSlices-1) - 2) % HTtoCAratio);
                    if (MyYOffset+MyYSlices-1 >= ny-3) UpperYBound = ny-3;
                    else UpperYBound = MyYOffset + MyYSlices - 1 + HTtoCAratio - ((MyYOffset + (MyYSlices-1) - 2) % HTtoCAratio);
                    
                    // Second pass through the files - ignore header lines (there are now 5)
                    string DummyLine;
                    for (int line=0; line<5; line++) {
                        getline(TemperatureFile,DummyLine);
                    }
                    // Read data from the remaining lines
                    while (!TemperatureFile.eof()) {
                       string s;
                       getline(TemperatureFile,s);
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
                          int XInt, YInt;
                          float XConverted, YConverted;
                          if (i == 0) stringstart = 0;
                          else stringstart = Subdivisions[i-1];
                          if (i == 4) stringlength = s.length();
                          else stringlength = Subdivisions[i]-stringstart;
                          string MeshDataS = s.substr(stringstart,stringlength);
                          //if (id == 0) cout << s << " SUBVAL = " << MeshDataS << endl;
                          float MeshData = atof(MeshDataS.c_str());
                          if (i == 0) {
                              if (LengthUnits == "mm") XConverted = MeshData/1000.0;
                              else XConverted = MeshData;
                              XInt = round((XConverted-XMin)/deltax) + 2;
                          }
                          else if (i == 1) {
                              if (LengthUnits == "mm") YConverted = MeshData/1000.0;
                              else YConverted = MeshData;
                              YInt = round((YConverted-YMin)/deltax) + 2;
                          }
                          else if (i == 2) {
                              if ((XInt >= LowerXBound)&&(XInt <= UpperXBound)&&(YInt >= LowerYBound)&&(YInt <= UpperYBound)) {
                                  // This data point is inside the bounds of interest for this MPI rank
                                  RawData[NumberOfTemperatureDataPoints] = XConverted;
                                  NumberOfTemperatureDataPoints++;
                                  RawData[NumberOfTemperatureDataPoints] = YConverted;
                                  NumberOfTemperatureDataPoints++;
                                  if (LengthUnits == "mm") RawData[NumberOfTemperatureDataPoints] = MeshData/1000.0;
                                  else RawData[NumberOfTemperatureDataPoints] = MeshData;
                                  NumberOfTemperatureDataPoints++;
                              }
                              else {
                                  // This data point is outside the bounds of interest for this MPI rank - set "i" outside loop bounds to exit
                                  i = 6;
                              }
                          }
                          else if (i == 3) {
                              // Liquidus time
                              if ((XInt >= LowerXBound)&&(XInt <= UpperXBound)&&(YInt >= LowerYBound)&&(YInt <= UpperYBound)) {
                                  if (TimeUnits == "ms") RawData[NumberOfTemperatureDataPoints] = MeshData/1000.0;
                                  else RawData[NumberOfTemperatureDataPoints] = MeshData;
                                  NumberOfTemperatureDataPoints++;
                              }
                          }
                          else if (i == 4) {
                              // Solidus time OR cooling rate (convert all to a positive cooling rate)
                              if ((XInt >= LowerXBound)&&(XInt <= UpperXBound)&&(YInt >= LowerYBound)&&(YInt <= UpperYBound)) {
                                  if (UseCoolingRate) {
                                      if (TimeUnits == "ms") RawData[NumberOfTemperatureDataPoints] = MeshData/1000.0;
                                      else RawData[NumberOfTemperatureDataPoints] = MeshData;
                                  }
                                  else {
                                      float TimeLiq = RawData[NumberOfTemperatureDataPoints-1];
                                      float TimeSol;
                                      if (TimeUnits == "ms") TimeSol = MeshData/1000.0;
                                      else TimeSol = MeshData;
                                      RawData[NumberOfTemperatureDataPoints] = abs(FreezingRange/(TimeSol-TimeLiq));
                                  }
                                  NumberOfTemperatureDataPoints++;
                              }
                          }
                          if (NumberOfTemperatureDataPoints >= RawData.size()-5) {
                              int OldSize = RawData.size();
                              RawData.resize(OldSize+1000000);
                          }
                      }
                   }
                   LastValue[LayerReadCount-1] = NumberOfTemperatureDataPoints;
                }
                TemperatureFile.close();
            } // End loop over all files read for all layers

            if (Loop == 0) {
                // Extend domain in Z (build) direction if the number of layers are simulated is greater than the number of temperature files read
                if (NumberOfLayers > TempFilesInSeries) {
                    for (int LayerReadCount=TempFilesInSeries; LayerReadCount<NumberOfLayers; LayerReadCount++) {
                        if (TempFilesInSeries == 1) {
                            // Only one temperature file was read, so the upper Z bound should account for an additional "NumberOfLayers-1" worth of data
                            // Since all layers have the same temperature data, each layer's "ZMinLayer" is just translated from that of the first layer
                            ZMinLayer[LayerReadCount] = ZMinLayer[LayerReadCount-1] + deltax*LayerHeight;
                            ZMaxLayer[LayerReadCount] = ZMaxLayer[LayerReadCount-1] + deltax*LayerHeight;
                            ZMax += deltax*LayerHeight;
                        }
                        else {
                            // "TempFilesInSeries" temperature files was read, so the upper Z bound should account for an additional "NumberOfLayers-TempFilesInSeries" worth of data
                            int RepeatedFile = (LayerReadCount) % TempFilesInSeries;
                            int RepeatUnit = LayerReadCount/TempFilesInSeries;
                            ZMinLayer[LayerReadCount] = ZMinLayer[RepeatedFile] + RepeatUnit*TempFilesInSeries*deltax*LayerHeight;
                            ZMaxLayer[LayerReadCount] = ZMaxLayer[RepeatedFile] + RepeatUnit*TempFilesInSeries*deltax*LayerHeight;
                            ZMax += deltax*LayerHeight;
                        }
                    }
                }
            
                // Now at the conclusion of "Loop 0", the decomposition can be performed as the domain bounds are known (all header lines from all files have been read)
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
                InitialDecomposition(DecompositionStrategy, nx, ny, ProcessorsInXDirection, ProcessorsInYDirection, id, np, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyLeftOut, MyRightIn, MyRightOut);

                MyXOffset = XOffsetCalc(id,nx,ProcessorsInXDirection,ProcessorsInYDirection,DecompositionStrategy);
                MyXSlices = XMPSlicesCalc(id,nx,ProcessorsInXDirection,ProcessorsInYDirection,DecompositionStrategy);

                MyYOffset = YOffsetCalc(id,ny,ProcessorsInYDirection,np,DecompositionStrategy);
                MyYSlices = YMPSlicesCalc(id,ny,ProcessorsInYDirection,np,DecompositionStrategy);
            }
            else {
                RawData.resize(NumberOfTemperatureDataPoints);
                // Determine start values for each layer's data within "RawData"
                if (NumberOfLayers > TempFilesInSeries) {
                    for (int LayerReadCount=TempFilesInSeries; LayerReadCount<NumberOfLayers; LayerReadCount++) {
                        if (TempFilesInSeries == 1) {
                            // Since all layers have the same temperature data, each layer's "ZMinLayer" is just translated from that of the first layer
                            FirstValue[LayerReadCount] = FirstValue[LayerReadCount-1];
                            LastValue[LayerReadCount] = LastValue[LayerReadCount-1];
                        }
                        else {
                            // All layers have different temperature data but in a repeating pattern
                            int RepeatedFile = (LayerReadCount) % TempFilesInSeries;
                            FirstValue[LayerReadCount] = FirstValue[RepeatedFile];
                            LastValue[LayerReadCount] = LastValue[RepeatedFile];
                        }
                    }
                }
            }
        } // End for loop iterating over file reads twice, Loop 0 for header lines and decomposing the domain, Loop 1 for storing the temperature data
    }

    else {
        // For constrained solidification test problems, nx/ny/nz are already known and the decomposition can be performed immediately
        if (id == 0) {
            cout << "Constrained solidification domain size: " << nx << " by " << ny << " by " << nz << endl;
            cout << "================================================================" << endl;
        }
        InitialDecomposition(DecompositionStrategy, nx, ny, ProcessorsInXDirection, ProcessorsInYDirection, id, np, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyLeftOut, MyRightIn, MyRightOut);

        MyXOffset = XOffsetCalc(id,nx,ProcessorsInXDirection,ProcessorsInYDirection,DecompositionStrategy);
        MyXSlices = XMPSlicesCalc(id,nx,ProcessorsInXDirection,ProcessorsInYDirection,DecompositionStrategy);

        MyYOffset = YOffsetCalc(id,ny,ProcessorsInYDirection,np,DecompositionStrategy);
        MyYSlices = YMPSlicesCalc(id,ny,ProcessorsInYDirection,np,DecompositionStrategy);
    }
}

//*****************************************************************************/
void TempInit(int layernumber, int TempFilesInSeries, double G, double R, string SimulationType, int id, int &MyXSlices, int &MyYSlices, int &MyXOffset, int &MyYOffset, double deltax, double HT_deltax, double deltat, int &nx, int &ny, int &nz, ViewI_H CritTimeStep, ViewF_H UndercoolingChange, ViewF_H UndercoolingCurrent, float XMin, float YMin, float ZMin, bool* Melted, float* ZMinLayer, float* ZMaxLayer, int LayerHeight, int NumberOfLayers, int &nzActive, int &ZBound_Low, int &ZBound_High, int* FinishTimeStep, double FreezingRange, ViewI_H LayerID, int* FirstValue, int* LastValue, vector <float> RawData) {

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
       
       for (int LayerCounter=0; LayerCounter<NumberOfLayers; LayerCounter++) {
           
           double SmallestTime = 1000000000;
           double SmallestTime_Global = 1000000000;
           double LargestTime = 0;
           double LargestTime_Global = 0;
           
           // How many CA cells in the vertical direction are needed to hold this layer's temperature data?
           int nzTempValuesThisLayer = round((ZMaxLayer[LayerCounter] - ZMinLayer[LayerCounter])/deltax) + 1; // (note this doesn't include the 2 rows of wall/active cells at the bottom surface)
           if (id == 0) cout << "Initializing temporary temperature data structures with " << nzTempValuesThisLayer << " cells in z direction" << endl;
           if (id == 0) cout << "Layer " << LayerCounter << " rank " << id << " ZMin this layer is " << ZMinLayer[LayerCounter] << endl;
           vector <vector <vector <double> > > CR, CritTL;
           for (int k=0; k<nzTempValuesThisLayer; k++) {
               vector <vector <double> > TemperatureXX;
               for (int i=LowerXBound; i<=UpperXBound; i++) {
                   vector <double> TemperatureX;
                   for (int j=LowerYBound; j<=UpperYBound; j++) {
                       TemperatureX.push_back(-1.0);
                   }
                   TemperatureXX.push_back(TemperatureX);
               }
               CR.push_back(TemperatureXX);
               CritTL.push_back(TemperatureXX);
           }
           // Data was already read into the "RawData" temporary data structure
           // Determine which section of "RawData" is relevant for this layer of the overall domain
           int StartRange = FirstValue[LayerCounter];
           int EndRange = LastValue[LayerCounter];
           int XInt, YInt, ZInt;
           if (id == 0) cout << "Range for layer " << LayerCounter << " on rank 0 is " << StartRange << " to " << EndRange << endl;
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
                    ZInt = round((RawData[i] + deltax*LayerHeight*LayerCounter - ZMinLayer[LayerCounter])/deltax);
                }
                else if (Pos == 3) {
                    CritTL[ZInt][XInt-LowerXBound][YInt-LowerYBound] = RawData[i];
                    if (RawData[i] < SmallestTime) SmallestTime = RawData[i];
                }
                else if (Pos == 4) {
                    CR[ZInt][XInt-LowerXBound][YInt-LowerYBound] = RawData[i];
                    float SolidusTime = CritTL[ZInt][XInt-LowerXBound][YInt-LowerYBound] + FreezingRange/CR[ZInt][XInt-LowerXBound][YInt-LowerYBound];
                    if (SolidusTime > LargestTime) LargestTime = SolidusTime;
                }
            }
            // If reading data from files without a script, time values start at 0 for each layer
            // If reading data with input from a script time values each layer are continuous, are should be renormalized to 0 for each layer
            MPI_Reduce(&LargestTime, &LargestTime_Global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Bcast(&LargestTime_Global,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Reduce(&SmallestTime, &SmallestTime_Global, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
            MPI_Bcast(&SmallestTime_Global,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

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
                           Pt1 = CR[LowZ][LowX][LowY];
                           Pt2 = CR[LowZ][HighX][LowY];
                           Pt12 = FLowX*Pt1 + FHighX*Pt2;
                           Pt3 = CR[LowZ][LowX][HighY];
                           Pt4 = CR[LowZ][HighX][HighY];
                           Pt34 = FLowX*Pt3 + FHighX*Pt4;
                           Pt1234 = Pt12*FLowY + Pt34*FHighY;
                           Pt5 = CR[HighZ][LowX][LowY];
                           Pt6 = CR[HighZ][HighX][LowY];
                           Pt56 = FLowX*Pt5 + FHighX*Pt6;
                           Pt7 = CR[HighZ][LowX][HighY];
                           Pt8 = CR[HighZ][HighX][HighY];
                           Pt78 = FLowX*Pt7 + FHighX*Pt8;
                           Pt5678 = Pt56*FLowY + Pt78*FHighY;
                           if ((Pt1 > 0)&&(Pt2 > 0)&&(Pt3 > 0)&&(Pt4 > 0)&&(Pt5 > 0)&&(Pt6 > 0)&&(Pt7 > 0)&&(Pt8 > 0)) {
                               CR[k][i][j] = Pt1234*FLowZ + Pt5678*FHighZ;
                           }
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
                           double CTLiq = CritTL[k][ii-LowerXBound][jj-LowerYBound] - LayerwiseTSOffset;
                           if (CTLiq > 0)  {
                               // Where does this layer's temperature data belong on the global (including all layers) grid?
                               // Adjust Z coordinate by ZMin
                               int ZOffset = round((ZMinLayer[LayerCounter]-ZMin)/deltax) + k + 2;
                               int Coord3D1D = ZOffset*MyXSlices*MyYSlices + Adj_i*MyYSlices + Adj_j;
                               Melted[Coord3D1D] = true;
                               CritTimeStep(Coord3D1D) = round(CTLiq/deltat);
                               LayerID(Coord3D1D) = LayerCounter;
                               UndercoolingChange(Coord3D1D) = abs(CR[k][ii-LowerXBound][jj-LowerYBound])*deltat;
                           }
                       }
                   }
               }
           }
        } // End read over all temperature files and placement of data

       if (id == 0) cout << "First layer Z bounds are " << ZBound_Low << " and " << ZBound_High << endl;
    }


}
//*****************************************************************************/
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
            
            Comp++;
            if (Comp > 2) {
                Comp = 0;
                UVNumber++;
            }
        }
    }
    O.close();
    
    // The grain orientations that correspond to each Grain ID must be the same across all ranks
    // Shuffle list of "NGrainOrientation" orientations
    int* GrainOrientation_master = new int[NGrainOrientations];
    if (id == 0) {
        for (int h=0; h<NGrainOrientations; h++) {
            GrainOrientation_master[h] = h;
        }
    }
    MPI_Bcast(&GrainOrientation_master[0], NGrainOrientations, MPI_INT, 0, MPI_COMM_WORLD);
    for (int h=0; h<NGrainOrientations; h++) {
        GrainOrientation(h) = GrainOrientation_master[h];
    }
    
}


//*****************************************************************************/
// Initializes cell types where the substrate comes from a file
void GrainInit(int layernumber, string SimulationType, string SubstrateFileName, double FractSurfaceSitesActive, int NGrainOrientations, int DecompositionStrategy, int nx, int ny, int nz, int LocalActiveDomainSize, int MyXSlices, int MyYSlices, int MyXOffset, int MyYOffset, int id, int np, int MyLeft, int MyRight, int MyIn, int MyOut, int MyLeftIn, int MyRightIn, int MyLeftOut, int MyRightOut, ViewI2D_H ItList, ViewI_H NeighborX, ViewI_H NeighborY, ViewI_H NeighborZ, ViewI_H GrainOrientation, ViewF_H GrainUnitVector, ViewF_H DiagonalLength, ViewI_H CellType, ViewI_H GrainID, ViewF_H CritDiagonalLength, ViewF_H DOCenter, ViewI_H CritTimeStep, ViewF_H UndercoolingChange, bool* Melted, double deltax, double NMax, int &NextLayer_FirstNucleatedGrainID, int &PossibleNuclei_ThisRank, int ZBound_High, int ZBound_Low, bool ExtraWalls) {
    
    mt19937_64 gen(id);
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
                        CellType(CAGridLocation) = Liquid;
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
                int SBuf = FirstEpitaxialGrainID+SubstrateActCells_ThisRank;
                MPI_Send(&SBuf,1,MPI_INT,1,0,MPI_COMM_WORLD);
            }
            else if (id == np-1) {
                int RBuf;
                MPI_Recv(&RBuf,1,MPI_INT,np-2,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                FirstEpitaxialGrainID = RBuf;
            }
            else {
                int RBuf;
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
        if (nz > nzS) {
            for (int k=nzS; k<nz; k++) {
                for (int j=0; j<nyS; j++) {
                    for (int i=0; i<nxS; i++) {
                        if ((i >= Substrate_LowX)&&(i < Substrate_HighX)&&(j >= Substrate_LowY)&&(j < Substrate_HighY)) {
                            int CAGridLocation;
                            CAGridLocation = k*MyXSlices*MyYSlices + (i-MyXOffset)*MyYSlices + (j-MyYOffset);
                            GrainID(CAGridLocation) = 0;
                        }
                    }
                }
            }
        }
        if (id == 0) cout << "Substrate file read complete" << endl;
        
        if (ExtraWalls) {
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
                           // This is a liquid cell or a nuclei, if not in a ghost node and the RNG places one at this site
                           double R = dis(gen);
                           if (R < BulkProb) {
                               if ((i != 0)&&(i != MyXSlices-1)&&(j != 0)&&(j != MyYSlices-1)) {
                                   PossibleNuclei_ThisRank++;
                                   CellType(CAGridLocation) = LiqSol;
                               }
                               else {
                                   CellType(CAGridLocation) = Liquid;
                               }
                           }
                           else {
                               CellType(CAGridLocation) = Liquid;
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
                                if ((CritTimeStep(NeighborD3D1ConvPosition) > 0)&&(CellType(NeighborD3D1ConvPosition) != Wall)) {
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

    // Set up active cell octahedra for growth, mark cell data to be communicated across ranks in ghost nodes
    // Assign GrainIDs to nuclei sites
    int NCounter = MyFirstNGrainID;

    // Nonactive cells should start with diagonal lengths of 0
    for (int i=0; i<LocalActiveDomainSize; i++) {
        DiagonalLength(i) = 0.0;
    }
    for (int GlobalZ=1; GlobalZ<nz-1; GlobalZ++) {
        for (int RankX=0; RankX<MyXSlices; RankX++) {
            for (int RankY=0; RankY<MyYSlices; RankY++) {
                long int D3D1ConvPositionGlobal = GlobalZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                if (CellType(D3D1ConvPositionGlobal) == Active) {

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

    PossibleNuclei_ThisRank += (ARNCount + BRNCount + CRNCount + DRNCount + ERNCount + FRNCount + GRNCount + HRNCount);

    // Remove delay cells not bordering others
    for (int RankZ=1; RankZ<nz-1; RankZ++) {
        for (int RankX=1; RankX<MyXSlices-1; RankX++) {
            for (int RankY=1; RankY<MyYSlices-1; RankY++) {
                int D3D1ConvPosition = RankZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
                if (CellType(D3D1ConvPosition) == Liquid) {
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
                    }
                }
            }
        }
    }
}

//*****************************************************************************/
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
    ViewI_H GhostNodesAR("bufferAR", 4*ARCount);
    ViewI_H GhostNodesBR("bufferBR", 4*BRCount);
    ViewI_H GhostNodesCR("bufferCR", 4*CRCount);
    ViewI_H GhostNodesDR("bufferDR", 4*DRCount);
    ViewI_H GhostNodesER("bufferER", 3*ERCount);
    ViewI_H GhostNodesFR("bufferFR", 3*FRCount);
    ViewI_H GhostNodesGR("bufferGR", 3*GRCount);
    ViewI_H GhostNodesHR("bufferHR", 3*HRCount);
    
    //MPI_Barrier(MPI_COMM_WORLD);

    // Collect ghost node data and send to other ranks- left and right
    if (ACount > 0) {
        ViewI_H GhostNodesA("bufferA", 4*ACount);
        for (int i=0; i<4*ACount; i++) {
            GhostNodesA(i) = ANuc[i];
        }
        if (BRCount == 0) {
            // Sending data to id = id - 1 only
	    MPI_Send(GhostNodesA.data(),ACount*4,MPI_INT,MyLeft,0,MPI_COMM_WORLD);
        }
        else {
            // Sending data to id = id - 1 and recieving data from id = id + 1
	    MPI_Sendrecv(GhostNodesA.data(),ACount*4,MPI_INT,MyLeft,0,GhostNodesBR.data(),BRCount*4,MPI_INT,MyRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    }
    else if (BRCount > 0) {
        // Recieving data from id = id + 1 only
        MPI_Recv(GhostNodesBR.data(),BRCount*4,MPI_INT,MyRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }

    if (BCount > 0) {
        ViewI_H GhostNodesB("bufferB", 4*BCount);
        for (int i=0; i<4*BCount; i++) {
            GhostNodesB(i) = BNuc[i];
        }
        if (ARCount == 0) {
            // Sending data to id = id + 1 only
            MPI_Send(GhostNodesB.data(),BCount*4,MPI_INT,MyRight,1,MPI_COMM_WORLD);
        }
        else {
            // Sending data to id = id + 1 and recieving data from id = id - 1
	    MPI_Sendrecv(GhostNodesB.data(),BCount*4,MPI_INT,MyRight,1,GhostNodesAR.data(),ARCount*4,MPI_INT,MyLeft,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    }
    else if (ARCount > 0) {
        // Recieving data from id = id - 1 only
        MPI_Recv(GhostNodesAR.data(),ARCount*4,MPI_INT,MyLeft,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
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
            NEvent++;
        }
    }

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
                NEvent++;
            }
        }
    }
    cout << "(" << id << ": " << NEvent << ") " << flush;
}


//*****************************************************************************/
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
    if (id == 0) cout << "New domain top at Z = " << ZBound_High << endl;
    // Change in active region data structures' sizes
    nzActive = ZBound_High - ZBound_Low + 1;
    LocalActiveDomainSize = MyXSlices*MyYSlices*nzActive;
    
    // Change in height of buffers
    BufSizeZ = nzActive;
    

    if (id == 0) cout << "New active domain height is " << nzActive << endl;
    
}

//*****************************************************************************/
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
                CritDiagonalLength((long int)(26)*D3D1ConvPosition+(long int)(n)) = CDLVal;
            }
        }
    });
    
    // Reset lock values
    Kokkos::parallel_for("LockInit",LocalActiveDomainSize, KOKKOS_LAMBDA (const int& D3D1ConvPosition) {
        int RankZ = D3D1ConvPosition/(MyXSlices*MyYSlices);
        int Rem = D3D1ConvPosition % (MyXSlices*MyYSlices);
        int RankX = Rem/MyYSlices;
        int RankY = Rem % MyYSlices;
        int GlobalZ = ZBound_Low + RankZ;
        int GlobalD3D1ConvPosition = GlobalZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
        if ((CellType(GlobalD3D1ConvPosition) == LiqSol)||(CellType(GlobalD3D1ConvPosition) == Liquid)) Locks(D3D1ConvPosition) = 1;
        else Locks(D3D1ConvPosition) = 0;
    });
}
