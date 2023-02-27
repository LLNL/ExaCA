#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Command line inputs are the list of temperature files associated with the ensemble of simulations
# See top level README for more details regarding these inputs
# Generate input files for ExaCA using a level 3 tasmanian sparse grid and 3 ExaCA
# input variables: heterogenous nucleation density, mean substrate grain size,
# and mean nucleation undercooling
import Tasmanian
import sys

# Generate Tasmanian grid
grid = Tasmanian.SparseGrid()
grid.makeLocalPolynomialGrid(3, 3, 3, 1, 'localp')
points = grid.getPoints()

# Min, max nucleation density (N0Min * 10^12 /m^3 and N0Max * 10^12 /m^3)
N0Min = 100 # 10^14
N0Max = 10000 # 10^16
N0Mean = 0.5 * (N0Min + N0Max)
N0Dev = N0Mean - N0Min

# Min, max mean nucleation undercooling (dTNMin and dTNMax, in K relative to
# the alloy liquidus temperature)
dTNMin = 5
dTNMax = 65
dTNMean = 0.5 * (dTNMin + dTNMax)
dTNDev = dTNMean - dTNMin

# Min, max substrate grain size (S0Min and S0Max, in microns)
S0Min = 2.5
S0Max = 25
S0Mean = 0.5 * (S0Min + S0Max)
S0Dev = S0Mean - S0Min

# Check number of command line arguments given: should be at least 2
NumCommandLineArgs = len(sys.argv);
if (NumCommandLineArgs == 1):
    print('Error: At least one temperature data file associated with these simulations is required')
    sys.exit(1)

# Get list of temperature files from the command line
TemperatureFiles = "["
for tfile in range(1, NumCommandLineArgs):
    TemperatureFiles += "\"" + str(sys.argv[tfile]) + "\"";
    if (tfile == NumCommandLineArgs-1):
        TemperatureFiles += "]";
    else:
        TemperatureFiles += ",";

# Write ExaCA input files to the examples subdirectory
for filenumber in range(1, 70):
    # Heterogenous nucleation density for this ensemble member
    N0ThisMember = N0Mean + N0Dev * points[filenumber-1,0]
    # Mean nucleation undercooling for this ensemble member
    dTNThisMember = dTNMean + dTNDev * points[filenumber-1,1]
    # Substrate grain spacing for this ensemble member
    S0ThisMember = S0Mean + S0Dev * points[filenumber-1,2]
    
    # Write to example file number "filenumber"
    filename = "examples/Inp_TasmanianTest_" + str(filenumber) + ".json"
    with open(filename, "w") as f:
        OutputData = ["{ \n",
        "   \"SimulationType\": \"RM\",\n",
        "   \"MaterialFileName\": \"Inconel625.json\",\n",
        "   \"GrainOrientationFile\": \"GrainOrientationVectors.csv\",\n",
        "   \"RandomSeed\": 0.0, \n",
        "   \"Domain\": { \n",
        "      \"CellSize\": 2.5, \n",
        "      \"TimeStep\": 0.125, \n",
        "      \"NumberOfLayers\": 56, \n",
        "      \"LayerOffset\": 8 \n",
        "   }, \n",
        "   \"Nucleation\": { \n",
        "      \"Density\": " + str(N0ThisMember) + ",\n",
        "      \"MeanUndercooling\": " + str(dTNThisMember) + ",\n",
        "      \"StDev\": 0.5 \n"
        "   }, \n",
        "   \"TemperatureData\": { \n",
        "      \"TemperatureFiles\": " + TemperatureFiles + "\n"
        "   }, \n",
        "   \"Substrate\": { \n",
        "      \"MeanSize\": " + str(S0ThisMember) + ",\n",
        "      \"PowderDensity\": 0 \n"
        "   }, \n",
        "   \"Printing\": { \n",
        "      \"PathToOutput\": \"./\" ,\n",
        "      \"OutputFile\": \"TasmanianTest_" + str(filenumber) + "\"" ",\n",
        "      \"PrintBinary\": true,\n",
        "      \"PrintFieldsInit\": [],\n",
        "      \"PrintFieldsFinal\": [\"GrainID\", \"LayerID\", \"GrainMisorientation\"] \n",
        "   } \n"]
        f.writelines(OutputData)
