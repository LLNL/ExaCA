# ExaCA problem types and auxiliary files
ExaCA currently can model three types of problems, two of which have the option of whether or not to include multiple melting and solidification events in cells:

* Problem type C is a directional solidification problem, with the bottom surface initialized with some fraction of sites home to epitaxial grains and at the liquidus temperature and a positive thermal gradient in the +Z direction. The domain is then cooled at a constant rate. 
* Problem type S is an array of hemispherical spots, with the number of spots in X, Y, and the number of layers for which the pattern is repeated (offset by a specified number of cells in the positive Z direction) specified. This problem type also uses fixed thermal gradient magnitude and cooling rate for each spot. 
* Problem type R is a custom solidification problem using time-temperature history file(s) (default location is `examples/Temperatures`). The format of these files are as follows:
    * The first line should be the names of the columns: x, y, z, tm, tl, cr
    * Each line following the first should have six comma-separated values corresponding to x, y, z, tm, tl, cr. x, y, and z are cell coordinates, in meters, of a given location in the simulation. The spacing between locations should correpond to a Cartesian grid, with a cell size equivalent to that specified in the input file. For each time that an x,y,z coordinate went above and below the liqiuidus temperature of the alloy during a heat transport simulation, a tm (time at which the point went above the liquidus), tl (time at which the point went below the liquidus), and cr (instantaneous cooling rate at the liquidus) should be recorded. 
    * If an x,y,z coordinate melted and solidified multiple times, it should appear in the file multiple times on separate lines. The order of the lines do not matter, except that the header line must be before any data.
    * The top surface (the largest Z coordinate in a file) is assumed to be flat. Additionally, if multiple temperature files are being used (for example, a scan pattern consisting of 10 layers of repeating even and odd file data), the Z coordinate corresponding to this flat top surface should be the same for all files.
* Problem types SM and RM modify problem types S and R to include multiple melting and solidification events per cell. For problem types S and R all cells that will eventually undergo melting are initialized as liquid, and only the final time that a given cell goes below the liquidus temperature is considered. To obtain the most accurate results, all melting and solidification events should be considered; however, for some problem geometries, the microstructure resulting from only considering the final solidification event in each cell is a reasonable approximation (and faster))

All problem types rely on two files in addition to the main input file. First,
a file containing the interfacial response function data governing
solidification rate as a function of undercooling is required. An example is
`examples/Materials/Inconel625`; if a different interfacial response function
of the same form is desired, a new Materials file can be created using the
Inconel625 file as a template and passed to the main input file. Second, a file
of grain orientations is required. An example is
`examples/Substrate/GrainOrientationVectors_Robert.csv`: the first line is the
number of orientations (10000), and each additional line is a list of unit
vectors corresponding to a cubic grain's <001> directions in the form 'x1, y1,
z1, x2, y2, z2, x3, y3, z3', where the coordinate system used is taken as the
ExaCA reference frame. The distribution of orientations is approximately even.
Like the material file, the orientation file could be swapped out with one
consisting of more (or fewer) orientations, following
`GrainOrientationVectors_Robert.csv` as a template. Both of these material and
orientation file examples are installed with the executable, making it possible
to simplfy use the file name in the input file. Custom files must either be
added to the ExaCA CMake build, use an absolute file path, or a path relative
to the ExaCA source.

Problems of type R or RM rely on a third file for temperature input, with the path
and name of this file given in the master input file. Examples of these
temperature field assembly files are given in
`examples/Temperatures/T_SimpleRaster.txt` and
`examples/Temperatures/T_AMBenchMultilayer.txt`. These files should always use
an absolute file path or path relative to the ExaCA source.

# ExaCA input files
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
| Extend baseplate through layers                             | N            | Value should be Y or N: Whether to use the baseplate microstructure as the boundary condition for the entire height of the simulation (default value is N)
| Density of powder surface sites active                      | See note (b) | Density of sites in the powder layer to be assigned as the home of a unique grain, normalized by 1 x 10^12 m^-3 (default value is 1/(CA cell size ^3)

(a) One of these inputs must be provided, but not both
(b) This is optional, but if this is given, Extend baseplate through layers must be set to N

The input "Extra set of wall cells", used in previous versions of ExaCA to create a set of wall cells "padding" the temperature data at the domain's X and Y boundaries, is deprecated and no longer has an effect.
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
