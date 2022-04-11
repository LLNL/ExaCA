# ExaCA problem types and auxiliary files
ExaCA currently can model three types of problems:

* Problem type C is a directional solidification problem, with the bottom surface initialized with some fraction of sites home to epitaxial grains and at the liquidus temperature and a positive thermal gradient in the +Z direction. The domain is then cooled at a constant rate. 
* Problem type S is an array of hemispherical spots, with the number of spots in X, Y, and the number of layers for which the pattern is repeated (offset by a specified number of cells in the positive Z direction) specified. This problem type also uses fixed thermal gradient magnitude and cooling rate for each spot. 
* Problem type R is a custom solidification problem using time-temperature history file(s) (default location is `examples/Temperatures`). The time-temperature history can consist of one file, repeated "Number of Layers" times, with each additional layer offset by "Layer height" cells in the +Z direction, or a pattern of temperature files. For example, if a simulation with alternating even and odd layer time-temperature histories is desired, the file names should be prepended with "1" and "2", and the "Number of temperature files in series" line in the input file should be set to 2. The temperature files should be of the following form:
    ** The first line should be the names of the columns: x, y, z, tm, tl, cr
    ** Each line following the first should have six comma-separated values corresponding to x, y, z, tm, tl, cr. x, y, and z are cell coordinates, in meters, of a given location in the simulation. The spacing between locations should correpond to a Cartesian grid, with a cell size equivalent to that specified in the input file. For each time that an x,y,z coordinate went above and below the liqiuidus temperature of the alloy during a heat transport simulation, a tm (time at which the point went above the liquidus), tl (time at which the point went below the liquidus), and cr (instantaneous cooling rate at the liquidus) should be recorded. 
    ** If an x,y,z coordinate melted and solidified multiple times, it should appear in the file multiple times on separate lines. The order of the lines do not matter, except that the header line must be before any data.
    ** The top surface (the largest Z coordinate in a file) is assumed to be flat. Additionally, if multiple temperature files are being used (for example, a scan pattern consisting of 10 layers of repeating even and odd file data), the Z coordinate corresponding to this flat top surface should be the same for all files.
    
All problem types will rely on two additional files. `examples/Materials/Inconel625` is a file containing the interfacial response function data governing solidification rate as a function of undercooling - if a different interfacial response function of the same form is desired, this line in an input file can be changed, and a new Materials file can be created using the Inconel625 file as a template. `examples/Substrate/GrainOrientationVectors_Robert.csv` is a file of grain orientations: the first line is the number of orientations (10000), and each additional line is a list of unit vectors corresponding to a cubic grain's <001> directions in the form 'x1, y1, z1, x2, y2, z2, x3, y3, z3', where the coordinate system used is taken as the ExaCA reference frame. The distribution of orientations is approximately even. Like the material file, the orientation file could be swapped out with one consisting of more (or fewer) orientations, following GrainOrientationVectors_Robert.csv as a template.

# ExaCA input files
The .txt files in the examples subdirectory are provided on the command line to let ExaCA know which problem is being simulated. Any lines prior to the first string of asterisks are ignored, as are any lines starting with an asterisk.

## All problems
The below lines are required regardless of problem type. They can be in any order, except Problem type must be given first:

|Input                   | Details |
|------------------------|---------|
| Problem type           | C for directional solidification (thermal gradient in build direction, fixed cooling rate)
|                        | S for spot melt array problem (fixed thermal gradient/constant cooling rate for each hemispherical spot)
|                        | R for use of temperature data provided in the appropriate format (see README file in examples/Temperatures)
| Decomposition strategy | 1 for a 1D domain decomposition along the Y direction
|                        | 2 for a decomposition along the Y direction, with a single partiton along the X direction
|                        | 3 for a decomposition with roughly equal partitions along the X and Y directions
|                        | If there is an odd number of processors, decomposition pattern defaults to 1
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

### Problem type S
All inputs are required - with the exception of sub grain size and sub filename: one of the two should be provided, but not both.

|Input                       | Details |
|----------------------------|---------|
| Thermal gradient  | Thermal gradient in each hemispherical spot, in K/m
| Cooling rate      | Cooling rate (uniform for each spot), in K/s
| Time step ratio   | Used to set the time step: time step = cell size/(cooling rate/thermal gradient) x N)
| Number of spots in x    | Number of spots in the x direction
| Number of spots in y    | Number of spots in the y direction
| Offset between spot centers | Offset of spot centers along the x and y axes, in microns
| Radii of spots    | Spot radii, in microns
| Number of layers  | Number of times this pattern is repeated, offset in the +Z (build) direction
| Offset between layers | If Number of layers > 1, the number of CA cells should separate adjacent layers
| Substrate grain spacing | Mean spacing between grain centers in the baseplate/substrate (in microns)
| Substrate filename| Filename for substrate data (either this OR Sub grain size should be provided, but not both)


### Problem type R
Some additional inputs are optional while others are required. As was the case for Problem type S, sub grain size and sub filename are both optional inputs, but one of the two must be provided.

|Input              | Required Y/N | Details |
|-------------------| - |---------|
| Substrate grain spacing | Y/N | Mean spacing between grain centers in the baseplate/substrate (in microns) (either this OR Sub filename should be given, but not both)
| Substrate filename      | Y/N | Filename (including path) for substrate data (either this OR Sub grain size should be provided, but not both)
| Time step         | Y | CA time step, in microseconds
| Temperature filename(s) | Y | Name of temperature file(s) (assumed glob before the given name if Number of temperature files > 1)
| Number of temperature files | Y | Number of temperature files used
| Number of layers  | Y | Number of times the pattern of "Temperature filename"/Number of temperature files pattern is repeated, offset in the +Z (build) direction
| Offset between layers | Y | If Number of layers > 1, the number of CA cells should separate adjacent layers
| Extra set of wall cells in lateral domain directions | Y | (Y or N value) should wall cells bound the domain in X and Y, rather than substrate grains (Y or N)?
| Heat transport data mesh size | N | Resolution of temperature data provided, in microns (if argument not provided, assumed to be equal to CA cell size)
| Path to temperature file(s) | N | Location of temperature data (if not provided, assumed to be located in examples/Temperatures)
| default RVE output | N | Whether or not to print representative volume element (RVE) data for ExaConstit, for a 0.5 cubic mm region in the domain center in X and Y, and at the domain top in Z excluding the final layer's grain structure

## Additional optional inputs for all problem types 
These values govern the printing intermediate data, for debugging or visualization, either following initialization or at specified increments during simulation

|Input                       | Details |
|----------------------------|---------|
| Debug check (reduced)  | (Y or N) Print data for main Kokkos views (CellType, LayerID, CritTimeStep) following initialization for debug check
| Debug check (extensive)| (Y or N) Print data for all main data structures (CellType, GrainID, CritTimeStep, UndercoolingCurrent, UndercoolingChange, Melted, LayerID) following initialization for debug check
| Print intermediate output frames | (Y or N) Print intermediate code output, where liquid CA cells are assigned the value -1, and non-liquid cells are assigned a value 0-62 depending on their associated grain's misorientation relative to the build (+Z) direction
| Increment to separate frames | If Print intermediate frames = Y, the number of microseconds defining the ExaCA output intermediate data increment
| Intermediate output even if system is unchanged from previous state | (Y or N) If Print intermediate frames = Y, whether or not ExaCA should print intermediate output regardless of whether the simulation has changed from the last frame (if Print intermediate frames = N, Print intermediate frames strict should also be = N)
| Random seed for grains and nuclei generation | Value of type double used as the seed to generate baseplate, powder, and nuclei details (default value is 0.0 if not provided)
