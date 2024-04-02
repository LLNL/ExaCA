# ExaCA problem types and auxiliary files
ExaCA currently can model three types of problems:

* Problem type C is a directional solidification problem, with the bottom surface initialized with some fraction of sites home to epitaxial grains and at the liquidus temperature and a positive thermal gradient in the +Z direction. The domain is then cooled at a constant rate. 
* Problem type `Spot` is a hemispherical spot, with fixed thermal gradient magnitude and cooling rate as solidification proceeds from the outer edge of the spot towards the center. The ability to simulate multilayer arrays of overlapping spots was deprecated in version 1.2 and removed after version 1.3
* Problem type `SingleGrain` is an initial nuclei at the domain center growing each time step until a domain edge is reached
* Problem type R or RM is a custom solidification problem using time-temperature history file(s) (default location is `examples/Temperatures`). Again note that some example problems require [external data](https://github.com/LLNL/ExaCA-Data). The format of these files are as follows:
    * The first line should be the names of the columns: x, y, z, tm, tl, cr
    * Each line following the first should have six comma-separated values corresponding to x, y, z, tm, tl, cr. x, y, and z are cell coordinates, in meters, of a given location in the simulation. The spacing between locations should correspond to a Cartesian grid, with a cell size equivalent to that specified in the input file. For each time that an x,y,z coordinate went above and below the liqiuidus temperature of the alloy during a heat transport simulation, a tm (time at which the point went above the liquidus), tl (time at which the point went below the liquidus), and cr (instantaneous cooling rate at the liquidus) should be recorded. As meters and seconds are the units used, and the cell size and time step tend to be on the order of micrometers and microseconds, it is recommended that this data be given as double precision values to avoid truncation of values
    * If an x,y,z coordinate melted and solidified multiple times, it should appear in the file multiple times on separate lines. The order of the lines do not matter, except that the header line must be before any data.
    * The top surface (the largest Z coordinate in a file) is assumed to be flat. Additionally, if multiple temperature files are being used (for example, a scan pattern consisting of 10 layers of repeating even and odd file data), the Z coordinate corresponding to this flat top surface should be the same for all files.
    * Alternatively, if a time-temperature history file has the extension `.catemp`, it will be parsed as a binary string. The binary form for these files does not contain commas, newlines, nor a header, but consists of sequential x,y,z,tm,tl,cr,x,y,z,tm,tl,cr... data as double precision values (little endian). This is often a significantly smaller file size than the standard format, and will be faster to read during initialization.
    * The M is no longer required in the problem type, as problem type R is now the same as RM (it now includes multiple melting and solidification events per cell)

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
to simplify use the file name in the input file. Custom files must either be
added to the ExaCA CMake build, use an absolute file path, or a path relative
to the ExaCA source. It should also be noted that if using a custom file of 
grain orientation vectors in an ExaCA simulations, corresponding files that 
list the orientations in bunge Euler angle form (analogous to 
`examples/Substrate/GrainOrientationEulerAnglesBungeZXZ.csv`) and as RGB values 
corresponding to the orientations' inverse pole figure-mapped RGB values 
(analogous to `examples/Substrate/GrainOrientationRGB_IPF-Z.csv`) will be 
necessary to run the analysis executable; see `analysis/README.md` for more 
details.

# ExaCA input files
The .json files in the examples subdirectory are provided on the command line to let ExaCA know which problem is being simulated. The json files contain different input sections, with certain inputs only used for certain problem types. For optional inputs, the default value used is noted.

## Top level inputs
|Input                   | Details |
|------------------------|---------|
| SimulationType         |
|                        | C for directional solidification (thermal gradient in build direction, fixed cooling rate)
|                        | S for spot melt array problem (fixed thermal gradient/constant cooling rate for each hemispherical spot)
|                        | R for use of temperature data provided in the appropriate format ([see below](#temperature-inputs))
|                        | `SingleGrain` for solidification of a single grain at the domain center, continuing until it reaches a domain edge
| MaterialFileName       | Name of material file in examples/Materials used
| GrainOrientationFile   | File listing rotation matrix components used in assigning orientations to grains ([see below](#substrate-inputs))
| RandomSeed             | Value of type double used as the seed to generate baseplate, powder, and nuclei details (default value is 0.0 if not provided)
| Domain                 | Section for parameters that describe the simulation domain for the given problem type ([see below](#domain-inputs))
| Nucleation             | Section for parameters that describe nucleation ([see below](#nucleation-inputs))
| TemperatureData        | Section for parameters/files governing the temperature field for the given problem type ([see below](#temperature-inputs))
| Substrate              | Section for parameters/files governing the edge boundary conditions ([see below](#substrate-inputs))
| Printing               | Section for parameters/file names for output data ([see below](#printing-inputs))

## Domain inputs
| Input        |Relevant problem type(s)| Details |
|--------------|------------------------|---------|
|CellSize      | All                    | CA cell size, in microns
|TimeStep      | All                    | CA time step, in microseconds (note previously for problem type C, this was derived from deltax, G, and R)
|Nx            | C, SingleGrain         | Domain size in x, in cells
|Ny            | C, SingleGrain         | Domain size in y, in cells
|Nz            | C, SingleGrain         | Domain size in z, in cells
|NumberOfLayers| R                      | Number of layers for which the temperature pattern will be repeated
|LayerOffset   | R                      | If numberOfLayers > 1, the offset (in cells) in the +Z direction for each layer of the temperature pattern
|SpotRadius    | Spot                   | Spot radius, in microns

## Nucleation inputs
| Input            |Relevant problem type(s)| Details |
|------------------|------------------------|---------|
|Density           | C, Spot, R             | Density of heterogeneous nucleation sites in the liquid (evenly distributed among cells that are liquid or undergo melting), normalized by 1 x 10^12 m^-3
|MeanUndercooling  | C, Spot, R             | Mean nucleation undercooling (relative to the alloy liquidus temperature) for activation of nucleation sites (Gaussian distribution)
|StDevUndercooling | C, Spot, R             | Standard deviation of nucleation undercooling (Gaussian distribution), in K

## Temperature inputs
| Input        | Relevant problem type(s)| Details |
|--------------|-------------------------|---------|
|G             | C, Spot, SingleGrain    | Thermal gradient in the build (+Z) directions, in K/m
|R             | C, Spot, SingleGrain    | Cooling rate (uniform across the domain), in K/s
|LayerwiseTempRead | R                   | If set to Y, the appropriate temperature data will be read during each layer's initialization, stored temporarily, and discarded. If set to N, temperature data for all layers will be read and stored during code initialization, and initialization of each layer will be performed using this stored temperature data. This option is only applicable to simulations with remelting; simulations without remelting (and simulations where this input is not given) default to N. Setting this to Y is only recommended if a large quantity of temperature data is read by ExaCA (for example, a 10 layer simulation where each layer's temperature data comes from a different file).
|TemperatureFiles | R                    | List of files corresponding to each layer's temperature data, in the form ["filename1.csv","filename2.csv",...]. If the number of entries is less than numberOfLayers, the list is repeated. Note that if the Z coordinate of the top surface for each data set has the layer offset applied, layerOffset in the "Domain" section of the input file should be set to 0, to avoid offsetting the layers twice.
|InitUndercooling | C, SingleGrain       | For SingleGrain, this is the undercooling at the location of the seeded grain. For problem type C, this is an optional argument (defaulting to zero) for the initial undercooling at the domain's bottom surface

## Substrate inputs
| Input        | Relevant problem type(s))| Details |
|--------------| -------------------------|---------|
|SurfaceSiteFraction | C           | What fraction of cells at the bottom surface of the domain are the source of a grain? (see note (b-i))
|SurfaceSiteDensity         | C           | Density, in nuclei/µm^2, of grains at the bottom surface of the domain (see note (b-ii))
|GrainLocationsX | C           | List of grain locations in X on the bottom surface of the domain (see note (b-iii))
|GrainLocationsY | C           | List of grain locations in Y on the bottom surface of the domain (see note (b-iii))
|GrainIDs | C           | GrainID values for each grain in (X,Y) (see note (b3))
|FillBottomSurface | C  | Optionally assign all cells on the bottom surface the grain ID of the closest grain (defaults to false)
|MeanSize      | Spot, R                  | Mean spacing between grain centers in the baseplate/substrate (in microns) (see note (a))
|SubstrateFilename |  Spot, R             | Path to and filename for substrate data (see note (a))
|PowderDensity | Spot, R                  | Density of sites in the powder layer to be assigned as the home of a unique grain, normalized by 1 x 10^12 m^-3 (default value is 1/(CA cell size ^3) (see note (a))
|ExtendSubstrateThroughPowder| R          | true/false value: Whether to use the baseplate microstructure as the boundary condition for the entire height of the simulation (defaults to false) (see note (a))
| BaseplateTopZ   | R                     | The Z coordinate that marks the top of the baseplate/boundary of the baseplate with the powder. If not given, Z = 0 microns will be assumed to be the baseplate top if ExtendSubstrateThroughPowder = false (If ExtendSubstrateThroughPowder = true, the entire domain will be initialized with the baseplate grain structure)
|GrainOrientation | SingleGrain           | Which orientation from the orientation's file is assigned to the grain (starts at 0). Default is 0 
|InitOctahedronSize | All           | Initial size of the octahedra that represent the solid-liquid interface when solidifiation first begins locally. Given as a fraction of a cell size, must be at least 0 and smaller than 1. Default is 0.01

(a) One of these inputs must be provided, but not both
(b) These represent 3 different ways of initializing the interface. Only the inputs corresponding to one mode (i, ii, or iii) should be given. Mode (i) assigns a GrainID value to randomly selected cells at the domain's bottom surface according to the input fixed fraction of sites active. However, this does not guarantee the same grains to be at the same physical (x,y) location if the cell size or domain size are changed. Mode (ii) determines a number of grains based on the input density, and assigns physical (x,y) locations. As only one GrainID per cell is allowed, a large density with a large cell size may lead to underresolution of the desired density - however, this approach does guarantee that the same grains will be initialized at the same physical (x,y) locations in the simulation with changes to cell size or domain size. Mode (iii) requires 3 inputs and allows manual initialization of the substrate with of lists of cell coordinates in X, cell coordinates in Y, and GrainID values.

## Printing inputs
| Input        | Relevant problem type(s))| Details |
|--------------| -------------------------|---------|
| PathToOutput | All                      | File path location for the output files
| OutputFile   | All                      | All output files will begin with the string specified on this line
| PrintBinary  | All                      | Whether or not ExaCA vtk output data should be printed as big endian binary data, or as ASCII characters (defaults to false)
| PrintExaConstitSize | R                 | Length of the cubic representative volume element (RVE) data for ExaConstit, taken from the domain center in X and Y, and at the domain top in Z excluding the final layer's grain structure. If not given (or given a value of 0), the RVE will not be printed 
| Intralayer   | All | Optional section for printing the state of the simulation during a given layer of a multilayer problem/during a single layer problem
| Intralayer: Increment | All | Increment, in time steps, at which intermediate output should be printed. If 0, will only print the state of the system at the start of each layer
| Intralayer: Fields | All | Fields to print during intralayer increments. Currently supported options are "GrainID", "LayerID", "GrainMisorientation", "UndercoolingCurrent", "UndercoolingSolidificationStart", ""MeltTimeStep", "CritTimeStep", "UndercoolingChange", "CellType", "DiagonalLength", "SolidificationEventCounter", "NumberOfSolidificationEvents"
| Intralayer: PrintIdleFrames | All | Whether or not ExaCA should print intermediate output regardless of whether the simulation has changed from the last frame
| Interlayer   | All | List of options for printing the state of the system following a given layer, or at the end of the run
| Interlayer: Layers | All | List of layers (starting at 0 and through "NumberOfLayers-1") following which the state of the simulation should be printed. If not given (or for non-multilayer problems), defaults to printing only after the full simulation has completed
| Interlayer: Increment | All | If "Interlayer: Layers" is not given, this option enables printing of interlayer output starting at layer 0 and repeating at the specified increment. The full simulation results following the final layer will always be printed.
| Interlayer: Fields | All | Fields to print following layers. Currently supported options are "GrainID", "LayerID", "GrainMisorientation", "UndercoolingCurrent", "UndercoolingSolidificationStart", "MeltTimeStep", "CritTimeStep", "UndercoolingChange", "CellType", "DiagonalLength", "SolidificationEventCounter", "NumberOfSolidificationEvents"

Here, GrainMisorientation is not the misorientation of the grain itself, but rather the misorientation of the grain's nearest <100> crystallographic direction with the +Z direction. For cells that are liquid (possible only for intermediate state print, as the final state will only have solid cells), -1 is printed as the misorienatation. Misorientations for grains from the baseplate or powder layer are between 0-62 (degrees, rounded to nearest integer), and cells that are associated with nucleated grains are assigned values between 100-162 to differentiate them. Additionally, 200 is printed as the misorientation for cells in the powder layer that have not been assigned a grain ID.
