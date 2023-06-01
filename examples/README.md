# ExaCA problem types and auxiliary files
ExaCA currently can model three types of problems, two of which have the option of whether or not to include multiple melting and solidification events in cells:

* Problem type C is a directional solidification problem, with the bottom surface initialized with some fraction of sites home to epitaxial grains and at the liquidus temperature and a positive thermal gradient in the +Z direction. The domain is then cooled at a constant rate. 
* Problem type S is an array of hemispherical spots, with the number of spots in X, Y, and the number of layers for which the pattern is repeated (offset by a specified number of cells in the positive Z direction) specified. This problem type also uses fixed thermal gradient magnitude and cooling rate for each spot. 
* Problem type R is a custom solidification problem using time-temperature history file(s) (default location is `examples/Temperatures`). The format of these files are as follows:
    * The first line should be the names of the columns: x, y, z, tm, tl, cr
    * Each line following the first should have six comma-separated values corresponding to x, y, z, tm, tl, cr. x, y, and z are cell coordinates, in meters, of a given location in the simulation. The spacing between locations should correpond to a Cartesian grid, with a cell size equivalent to that specified in the input file. For each time that an x,y,z coordinate went above and below the liqiuidus temperature of the alloy during a heat transport simulation, a tm (time at which the point went above the liquidus), tl (time at which the point went below the liquidus), and cr (instantaneous cooling rate at the liquidus) should be recorded. As meters and seconds are the units used, and the cell size and time step tend to be on the order of micrometers and microseconds, it is recommended that this data be given as double precision values to avoid truncation of values
    * If an x,y,z coordinate melted and solidified multiple times, it should appear in the file multiple times on separate lines. The order of the lines do not matter, except that the header line must be before any data.
    * The top surface (the largest Z coordinate in a file) is assumed to be flat. Additionally, if multiple temperature files are being used (for example, a scan pattern consisting of 10 layers of repeating even and odd file data), the Z coordinate corresponding to this flat top surface should be the same for all files.
    * Alternatively, if a time-temperature history file has the extension `.catemp`, it will be parsed as a binary string. The binary form for these files does not contain commas, newlines, nor a header, but consists of sequential x,y,z,tm,tl,cr,x,y,z,tm,tl,cr... data as double precision values (little endian). This is often a significantly smaller file size than the standard format, and will be faster to read during initialization.
    * Problem types SM and RM modify problem types S and R to include multiple melting and solidification events per cell. For problem types S and R all cells that will eventually undergo melting are initialized as liquid, and only the final time that a given cell goes below the liquidus temperature is considered. To obtain the most accurate results, all melting and solidification events should be considered; however, for some problem geometries, the microstructure resulting from only considering the final solidification event in each cell is a reasonable approximation (and faster)

All problem types rely on two files in addition to the main input file. First,
a file containing the interfacial response function data governing
solidification rate as a function of undercooling is required. An example is
`examples/Materials/Inconel625`; if a different interfacial response function
of the same form is desired, a new Materials file can be created using the
Inconel625 file as a template and passed to the main input file. Second, a file
of grain orientations is required. An example is
`examples/Substrate/GrainOrientationVectors.csv`: the first line is the
number of orientations (10000), and each additional line is a list of unit
vectors corresponding to a cubic grain's <001> directions in the form 'x1, y1,
z1, x2, y2, z2, x3, y3, z3', where the coordinate system used is taken as the
ExaCA reference frame. The distribution of orientations is approximately even.
Like the material file, the orientation file could be swapped out with one
consisting of more (or fewer) orientations, following
`GrainOrientationVectors.csv` as a template. Both of these material and
orientation file examples are installed with the executable, making it possible
to simplfy use the file name in the input file. Custom files must either be
added to the ExaCA CMake build, use an absolute file path, or a path relative
to the ExaCA source. It should also be noted that if using a custom file of 
grain orientation vectors in an ExaCA simulations, corresponding files that 
list the orientations in bunge Euler angle form (analogous to 
`examples/Substrate/GrainOrientationEulerAnglesBungeZXZ.csv`) and as RGB values 
corresponding to the orientations' inverse pole figure-mapped RGB values 
(analogous to `examples/Substrate/GrainOrientationRGB_IPF-Z.csv`) will be 
necessary to run the analysis executable; see `analysis/README.md` for more 
details.

If the main input file is a .txt file, problems of type R or RM will rely on a 
third file for temperature input, with the path and name of this file given in 
the master input file. Examples of these temperature field assembly files are 
given in `examples/Temperatures/T_SimpleRaster.txt` and 
`examples/Temperatures/T_AMBenchMultilayer.txt`. These files should always use 
an absolute file path or path relative to the ExaCA source. If the main input 
file is a .json file, this additional temperature input file is not required, 
as the relevant information is stored within the main input file.

# ExaCA input files
The .json files in the examples subdirectory are provided on the command line to let ExaCA know which problem is being simulated. The .txt files can also be used (see the next section for information on the format of these files), but compatibility with files of this older format will be removed in a future release. The json files contain different input sections, with certain inputs only used for certain problem types. Note that some of these inputs have been renamed or redefined from the old format, these changes are noted where appropriate. Each top level section is required for all problem types. For optional inputs, the default value used is noted.

## Top level inputs
|Input                   | Equivalent from old input file format | Details |
|------------------------|---------------------------------------|---------|
| SimulationType         | Problem Type                          |         |
|                        |                                       | `SingleGrain` for solidification of a single grain at the domain center, continuing until it reaches a domain edge
|                        |                                       | C for directional solidification (thermal gradient in build direction, fixed cooling rate)
|                        |                                       | S for spot melt array problem (fixed thermal gradient/constant cooling rate for each hemispherical spot)
|                        |                                       | R for use of temperature data provided in the appropriate format (see README file in examples/Temperatures)
|                        |                                       | M should be appended to problem type if multiple melting and solidifcation events are desired (i.e, SM or RM)
| MaterialFileName       |Material                               | Name of material file in examples/Materials used (see README file in examples/Materials)
| GrainOrientationFile   | File of grain orientations            | File listing rotation matrix components used in assigning orientations to grains (see README file in examples/Substrate)
| RandomSeed                | Random seed for grains and nuclei generation | Value of type double used as the seed to generate baseplate, powder, and nuclei details (default value is 0.0 if not provided)
| Domain                 | N/A | Section for parameters that describe the simulation domain for the given problem type (see below for second level inputs)
| Nucleation             | N/A | Section for parameters that describe nucleation (see below for second level inputs)
| TemperatureData        | N/A | Section for parameters/files governing the temperature field for the given problem type (see below for second level inputs)
| Substrate              | N/A | Section for parameters/files governing the edge boundary conditions (see below for second level inputs)
| Printing               | N/A | Section for parameters/file names for output data (see below for second level inputs)

## Domain inputs
| Input        | Equivalent from      | Relevant problem | Details |
|              | old input file format| type(s)          |         |    
|--------------| ---------------------|------------------| --------|
|CellSize      | deltax               | All              | CA cell size, in microns
|TimeStep      | deltat               | All              | CA time step, in microseconds (note previously for problem type C, this was derived from deltax, G, and R)
|Nx            | Domain size in x     | C                | Domain size in x, in cells
|Ny            | Domain size in y     | C                | Domain size in y, in cells
|Nz            | Domain size in z     | C                | Domain size in z, in cells
|NumberOfLayers| Number of layers     | S, R             | Number of layers for which the temperature pattern will be repeated
|LayerOffset   | Offset between layers| S, R             | If numberOfLayers > 1, the offset (in cells) in the +Z direction for each layer of the temperature pattern
|NSpotsX       | Number of spots in x | S                | Number of spots in the x direction
|NSpotsY       | Number of spots in y | S                | Number of spots in the y direction
|RSpots        | Radii of spots       | S                | Spot radii, in microns
|SpotOffset    | Offset between spot centers | S         | Offset of spot centers along the x and y axes, in microns     

## Nucleation inputs
| Input        | Equivalent from      | Relevant problem | Details |
|              | old input file format| type(s)          |         |    
|--------------| ---------------------|------------------| --------|
|Density           | Heterogeneous nucleation density             | C, S, R | Density of heterogenous nucleation sites in the liquid (evenly distributed among cells that are liquid or undergo melting), normalized by 1 x 10^12 m^-3
|MeanUndercooling  | Mean nucleation undercooling                 | C, S, R | Mean nucleation undercooling (relative to the alloy liquidus temperature) for activation of nucleation sites (Gaussian distribution)
|StDevUndercooling | Standard deviation of nucleation undercooling| C, S, R | Standard deviation of nucleation undercooling (Gaussian distribution), in K

## Temperature inputs
| Input        | Equivalent from      | Relevant problem | Details |
|              | old input file format| type(s)          |         |    
|--------------| ---------------------|------------------| --------|
|G             | Thermal gradient     | SingleGrain, C, S | Thermal gradient in the build (+Z) directions, in K/m
|R             | Cooling rate         | SingleGrain, C, S | Cooling rate (uniform across the domain), in K/s
|InitUndercooling | N/A               | SingleGrain      | Undercooling at the location of the seeded grain 
|HeatTransferCellSize | Heat transport data mesh size| R    | By default, equal to deltax, and cannot be used if remelting is considered. deltax must divide evenly into HTdeltax
|LayerwiseTempRead| Discard temperature data and reread temperature files after each layer | R | If set to Y, the appropriate temperature data will be read during each layer's initialization, stored temporarily, and discarded. If set to N, temperature data for all layers will be read and stored during code initialization, and initialization of each layer will be performed using this stored temperature data. This option is only applicable to simulations with remelting; simulations without remelting (and simulations where this input is not given) default to N. Setting this to Y is only recommended if a large quantity of temperature data is read by ExaCA (for example, a 10 layer simulation where each layer's temperature data comes from a different file).
|TemperatureFiles | N/A | R | List of files corresponding to each layer's temperature data, in the form ["filename1.csv","filename2.csv",...]. If the number of entries is less than numberOfLayers, the list is repeated. Note that if the Z coordinate of the top surface for each data set has the layer offset applied, layerOffset in the "Domain" section of the input file should be set to 0, to avoid offsetting the layers twice.

## Substrate inputs
| Input        | Equivalent from      | Relevant problem | Details |
|              | old input file format| type(s)          |         |    
|--------------| ---------------------|------------------| --------|
|FractionSurfaceSitesActive | Fraction surface sites active| C        | What fraction of cells at the bottom surface of the domain are the source of a grain?
|MeanSize      | Substrate grain spacing      | S, R     | Mean spacing between grain centers in the baseplate/substrate (in microns) (see note (a))
|SubstrateFilename | Substrate filename       | S, R     | Path to and filename for substrate data (see note (a))
|PowderDensity | Density of powder surface sites active | S, R | Density of sites in the powder layer to be assigned as the home of a unique grain, normalized by 1 x 10^12 m^-3 (default value is 1/(CA cell size ^3) (see note (b))
|ExtendSubstrateThroughPower| Extend baseplate through layers | S, R | true/false value: Whether to use the baseplate microstructure as the boundary condition for the entire height of the simulation (defaults to false) (see note (b))
|GrainOrientation | N/A               | SingleGrain      | Which orientation from the orientation's file is assigned to the grain (starts at 0). Default is 0 

(a) One of these inputs must be provided, but not both
(b) This is optional, but if this is given, "extendSubstrateThroughPower" must be set to false

## Printing inputs
| Input        | Equivalent from      | Relevant problem | Details |
|              | old input file format| type(s)          |         |    
|--------------| ---------------------|------------------| --------|
| PathToOutput | Path to output       | All              | File path location for the output files
| OutputFile   | Output file base name| All              | All output files will begin with the string specified on this line
| PrintBinary  | Print vtk data as binary |All           | Whether or not ExaCA vtk output data should be printed as big endian binary data, or as ASCII characters (defaults to false)
| PrintExaConstitSize | Size of the default RVE, in CA cells | R | Length of the cubic representative volume element (RVE) data for ExaConstit, taken from the domain center in X and Y, and at the domain top in Z excluding the final layer's grain structure. If not given (or given a value of 0), the RVE will not be printed 
| PrintFieldsInit     | N/A | All, except with remelting | Which fields to print following initialization, in the form ["Field1","Field2",...]. Currently supported options are GrainID, LayerID, CellType, CritTimeStep, UndercoolingChange, UndercoolingCurrent
| PrintFieldsFinal    | N/A | All | Which fields to print at the end of the run, in the form ["Field1","Field2",...]. Currently supported options are GrainID, LayerID, GrainMisorientation, UndercoolingCurrent
| PrintIntermediateOutput | Print intermediate output frames | All | If this section is given, intermediate code output will be printed, where liquid CA cells are assigned the value -1, and non-liquid cells are assigned a value 0-62 depending on their associated grain's misorientation relative to the build (+Z) direction. This is an optional section that otherwise defaults to false, with second level inputs "Frequency" and "PrintIdleFrames"
| PrintIntermediateOutput:Frequency | Increment to separate frames | All | The number of microseconds defining the ExaCA output intermediate data increment. A value of 0 will mean that no intermediate output is printed
| PrintIntermediateOutput:PrintIdleFrames | Intermediate output even if system is unchanged from previous state | All | Whether or not ExaCA should print intermediate output regardless of whether the simulation has changed from the last frame

# ExaCA input files: Old format
The .txt files in the examples subdirectory are provided on the command line to let ExaCA know which problem is being simulated. Any lines prior to the first string of asterisks are ignored, as are any lines starting with an asterisk.

## All problems
The below lines are required regardless of problem type. They can be in any order, except Problem type must be given first:

|Input                   | Details |
|------------------------|---------|
| Problem type           | C for directional solidification (thermal gradient in build direction, fixed cooling rate)
|                        | S for spot melt array problem (fixed thermal gradient/constant cooling rate for each hemispherical spot)
|                        | R for use of temperature data provided in the appropriate format (see README file in examples/Temperatures)
|                        | M should be appended to problem type if multiple melting and solidifcation events are desired (i.e, SM or RM)
| Material               | Name of material file in examples/Materials used (see README file in examples/Materials)
| Cell size              | CA cell size, in microns
| Heterogeneous nucleation density     | Density of heterogenous nucleation sites in the liquid (evenly distributed among cells that are liquid or undergo melting), normalized by 1 x 10^12 m^-3
| Mean nucleation undercooling       | Mean nucleation undercooling (relative to the alloy liquidus temperature) for activation of nucleation sites (Gaussian distribution)
| Standard deviation of nucleation undercooling| Standard deviation of nucleation undercooling (Gaussian distribution), in K
| Path to output        | File path location for the output files 
| Output file base name | All output files will begin with the string specified on this line
| File of grain orientations | File listing rotation matrix components used in assigning orientations to grains (see README file in examples/Substrate)
| file of final undercooling values | File listing the undercooling at which each cell became solid 
| Print file of grain misorientation values | (Y or N) Print a Paraview vtk file of grain misorientation relative to the build (+Z) direction: epitaxial grains will have values 0-62, while nucleated grains have values 100-162 (the +100 offset is used to differentiate the two types of grains)
| Print file of all ExaCA data | (Y or N) Print Paraview vtk files of all relevant ExaCA cell parameters (GrainID, LayerID, Melted) for use in post-processing

## Problem type specific inputs

### Problem type C
All inputs are required.

|Input                       | Details |
|----------------------------| --------|
| Thermal gradient  | Thermal gradient in the build (+Z) directions, in K/m
| Cooling rate      | Cooling rate (uniform across the domain), in K/s
| Time step ratio (from steady-state Velocity)   | Used to set the time step: time step = cell size/(cooling rate/thermal gradient) x N)
| Domain size in x  | Domain size in x
| Domain size in y  | Domain size in y
| Domain size in z  | Domain size in z
| Fraction surface sites active| What fraction of cells at the bottom surface of the domain are the source of a grain?

### Problem type S or SM
Some additional inputs are optional while others are required

|Input                                   | Required Y/N | Details |
|----------------------------------------|--------------|---------|
| Thermal gradient                       | Y            | Thermal gradient in each hemispherical spot, in K/m
| Cooling rate                           | Y            | Cooling rate (uniform for each spot), in K/s
| Time step ratio                        | Y            | Used to set the time step: time step = cell size/(cooling rate/thermal gradient) x N)
| Number of spots in x                   | Y            | Number of spots in the x direction
| Number of spots in y                   | Y            | Number of spots in the y direction
| Offset between spot centers            | Y            | Offset of spot centers along the x and y axes, in microns
| Radii of spots                         | Y            | Spot radii, in microns
| Number of layers                       | Y            | Number of times this pattern is repeated, offset in the +Z (build) direction. Max value is 32767
| Offset between layers                  | Y            | If Number of layers > 1, the number of CA cells should separate adjacent layers
| Substrate grain spacing                | See note (a) | Mean spacing between grain centers in the baseplate/substrate (in microns)
| Substrate filename                     | See note (a) | Path to and filename for substrate data
| Extend baseplate through layers        | N            | Value should be Y or N: Whether to use the baseplate microstructure as the boundary condition for the entire height of the simulation (default value is N)
| Density of powder surface sites active | See note (b) | Density of sites in the powder layer to be assigned as the home of a unique grain, normalized by 1 x 10^12 m^-3 (default value is 1/(CA cell size ^3)

(a) One of these inputs must be provided, but not both
(b) This is optional, but if this is given, "Extend baseplate through layers" must be set to N

### Problem type R or RM
Some additional inputs are optional while others are required

|Input                                                        | Required Y/N | Details |
|-------------------------------------------------------------| -------------|---------|
| Substrate grain spacing                                     | See note (a) | Mean spacing between grain centers in the baseplate/substrate (in microns)
| Substrate filename                                          | See note (a) | Path to and filename for substrate data
| Time step                                                   | Y            | CA time step, in microseconds
| Path to and name of temperature field assembly instructions | Y            | Additional file containing information on the temperature data
| default RVE output                                          | N            | Whether or not to print representative volume element (RVE) data for ExaConstit, for a 0.5 cubic mm region in the domain center in X and Y, and at the domain top in Z excluding the final layer's grain structure
| Default RVE size, in CA cells                               | N            | Length of the RVE, defaulting to the equivalent of 0.5 by 0.5 by 0.5 mm if not given
| Extend baseplate through layers                             | N            | Value should be Y or N: Whether to use the baseplate microstructure as the boundary condition for the entire height of the simulation (default value is N)
| Density of powder surface sites active                      | See note (b) | Density of sites in the powder layer to be assigned as the home of a unique grain, normalized by 1 x 10^12 m^-3 (default value is 1/(CA cell size ^3)

(a) One of these inputs must be provided, but not both
(b) This is optional, but if this is given, Extend baseplate through layers must be set to N

Additional information regarding the temperature field is provided inside the auxiliary file specified on the "Path to and name of temperature field assembly instructions" line. A temperature field construction file contains the following additional specifications:

|Input              | Required Y/N | Details |
|-------------------| - |---------|
| Number of layers  | Y | Number of times the files specified in the second half of this file will be repeated, with each temperature field offset in the +Z (build) direction
| Offset between layers | Y | If Number of layers > 1, the number of CA cells should separate adjacent layers
| Heat transport data mesh size | N | Resolution of temperature data provided, in microns (if argument not provided, assumed to be equal to CA cell size)
| Discard temperature data and reread temperature files after each layer | N | If set to Y, the appropriate temperature data will be read during each layer's initialization, stored temporarily, and discarded. If set to N, temperature data for all layers will be read and stored during code initialization, and initialization of each layer will be performed using this stored temperature data. This option is only applicable to simulations with remelting; simulations without remelting (and simulations where this input is not given) default to N. Setting this to Y is only recommended if a large quantity of temperature data is read by ExaCA (for example, a 10 layer simulation where each layer's temperature data comes from a different file).

A comment line starting with an asterisk separates the first half of the file, containing the above data, from the bottom half. The bottom half of the file consists of the temperature files (including the paths) used in construction of the temperature field. If there is one file, that temperature field will be repeated, offset by "Offset between Layers" cells in the build direction, for "Number of layers" layers. If there are multiple files, those temperature fields will be repeated in the same manner. For example, if there are two lines below the asterisks, "Even.txt" and "Odd.txt", Offset between layers = 5, and Number of layers = 7, layers 0, 2, 4, and 6 will use "Even.txt" data and layers 1, 3, and 5 will use "Odd.txt" data. ExaCA will offset the Z coordinates of each layer by 5 cells relative to the previous one; as a result, "Odd.txt" should not have a built in offset in the Z direction from "Even.txt", as this would result in the offset being added in twice. Examples temperature construction files are given in `examples/Temperatures/T_SimpleRaster.txt` and `examples/Temperatures/T_AMBenchMultilayer.txt`. 
The deprecated form for temperature field input data, where these 3 input lines exist in the top level input file, alongside inputs "Number of temperature files in series: N" and "Temperature filename(s): Data.txt" (which would indicate reading temperature data from files "1Data.txt", "2Data.txt".... "NData.txt", is still allowed but will be removed in a future release.

## Additional optional inputs for all problem types 
These values govern the printing intermediate data, for debugging or visualization, either following initialization or at specified increments during simulation. Debug check options are currently only available for simulations that do not include multiple melting/soldification events per cell (i.e., only problem types C, S, and R, not SM nor RM)

|Input                       | Details |
|----------------------------|---------|
| Debug check (reduced)  | (Y or N) Print data for main Kokkos views (CellType, LayerID, CritTimeStep) following initialization for debug check
| Debug check (extensive)| (Y or N) Print data for all main data structures (CellType, GrainID, CritTimeStep, UndercoolingCurrent, UndercoolingChange, Melted, LayerID) following initialization for debug check
| Print intermediate output frames | (Y or N) Print intermediate code output, where liquid CA cells are assigned the value -1, and non-liquid cells are assigned a value 0-62 depending on their associated grain's misorientation relative to the build (+Z) direction
| Increment to separate frames | If Print intermediate frames = Y, the number of microseconds defining the ExaCA output intermediate data increment
| Intermediate output even if system is unchanged from previous state | (Y or N) If Print intermediate frames = Y, whether or not ExaCA should print intermediate output regardless of whether the simulation has changed from the last frame (if Print intermediate frames = N, Print intermediate frames strict should also be = N)
| Random seed for grains and nuclei generation | Value of type double used as the seed to generate baseplate, powder, and nuclei details (default value is 0.0 if not provided)
| Print vtk data as binary | Whether or not ExaCA vtk output data should be printed as big endian binary data, or as ASCII characters (default value is false if not provided)
