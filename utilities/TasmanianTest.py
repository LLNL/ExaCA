#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Takes two command line inputs: the first is the base temperature file name, and
# the second is the number of temperature files in the series
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

# Check number of command line arguments given: should be 2
if len(sys.argv) == 3:
    BaseTemperatureFilename = str(sys.argv[1])
    TempFilesInSeries = str(sys.argv[2])
else:
    print('Error: exactly 2 command line arguments are required')
    sys.exit(1)

# Temperature filename without extension
StartExt = BaseTemperatureFilename.rfind('.')
BaseTemperatureFilenameNoExt = BaseTemperatureFilename[:StartExt]

# Write ExaCA input files to the examples subdirectory
for filenumber in range(1, 70):
    # Heterogenous nucleation density for this ensemble member
    N0ThisMember = N0Mean + N0Dev * points[filenumber-1,0]
    # Mean nucleation undercooling for this ensemble member
    dTNThisMember = dTNMean + dTNDev * points[filenumber-1,1]
    # Substrate grain spacing for this ensemble member
    S0ThisMember = S0Mean + S0Dev * points[filenumber-1,2]
    
    # Write to example file number "filenumber"
    filename = "examples/Inp_" + BaseTemperatureFilenameNoExt + "EnsembleMember" + str(filenumber) + ".txt"
    with open(filename, "w") as f:
        OutputData = ["ExaCA input file written in python for an ensemble of CA calculations \n",
        "***** \n",
        "Problem type: R \n",
        "Decomposition strategy: 1 \n",
        "Material:Inconel625 \n",
        "Cell size: 2.5 \n",
        "Heterogeneous nucleation density: " + str(N0ThisMember) + "\n",
        "Mean nucleation undercooling: " + str(dTNThisMember) + "\n",
        "Standard deviation of nucleation undercooling: 0.5 \n",
        "Path to output:./ \n",
        "Output file base name: " + BaseTemperatureFilenameNoExt + "_ExaCAEnsMem_" + str(filenumber) + "\n",
        "File of grain orientations:GrainOrientationVectors_Robert.csv \n",
        "Heat transport data mesh size: 2.5 \n",
        "Time step: 0.125 \n",
        "Substrate grain spacing: " + str(S0ThisMember) + "\n",
        "Path to temperature file(s): examples/Temperatures \n",
        "Temperature filename(s): " + BaseTemperatureFilename + "\n",
        "Number of temperature files: " + TempFilesInSeries + "\n",
        "Number of layers: 56 \n",
        "Offset between layers: 8 \n",
        "***Output data printing options: (Y or N) which data should be printed*** \n",
        "Print file of grain misorientation values: Y \n",
        "Print file of all ExaCA data: Y \n",
        "Print default RVE output: Y \n",
        "Debug check (reduced): N \n",
        "Debug check (extensive): N \n"]
        f.writelines(OutputData)
