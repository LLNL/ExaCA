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

# ExaCA input files
The .json files in the examples subdirectory are provided on the command line to let ExaCA know which problem is being simulated. The json files contain different input sections, with certain inputs only used for certain problem types. For optional inputs, the default value used is noted.

## Top level inputs
|Input                   | Details |
|------------------------|---------|
| SimulationType         |
|                        | C for directional solidification (thermal gradient in build direction, fixed cooling rate)
|                        | S for spot melt array problem (fixed thermal gradient/constant cooling rate for each hemispherical spot)
|                        | R for use of temperature data provided in the appropriate format (see README file in examples/Temperatures)
|                        | M should be appended to problem type if multiple melting and solidifcation events are desired (i.e, SM or RM)
| MaterialFileName       | Name of material file in examples/Materials used (see README file in examples/Materials)
| GrainOrientationFile   | File listing rotation matrix components used in assigning orientations to grains (see README file in examples/Substrate)
| RandomSeed             | Value of type double used as the seed to generate baseplate, powder, and nuclei details (default value is 0.0 if not provided)
| Domain                 | Section for parameters that describe the simulation domain for the given problem type (see below for second level inputs)
| Nucleation             | Section for parameters that describe nucleation (see below for second level inputs)
| TemperatureData        | Section for parameters/files governing the temperature field for the given problem type (see below for second level inputs)
| Substrate              | Section for parameters/files governing the edge boundary conditions (see below for second level inputs)
| Printing               | Section for parameters/file names for output data (see below for second level inputs)

## Domain inputs
| Input        |Relevant problem type(s)| Details |
|--------------|------------------------|---------|
|CellSize      | All                    | CA cell size, in microns
|TimeStep      | All                    | CA time step, in microseconds (note previously for problem type C, this was derived from deltax, G, and R)
|Nx            | C                      | Domain size in x, in cells
|Ny            | C                      | Domain size in y, in cells
|Nz            | C                      | Domain size in z, in cells
|NumberOfLayers| S, R                   | Number of layers for which the temperature pattern will be repeated
|LayerOffset   | S, R                   | If numberOfLayers > 1, the offset (in cells) in the +Z direction for each layer of the temperature pattern
|NSpotsX       | S                      | Number of spots in the x direction
|NSpotsY       | S                      | Number of spots in the y direction
|RSpots        | S                      | Spot radii, in microns
|SpotOffset    | S                      | Offset of spot centers along the x and y axes, in microns     

## Nucleation inputs
| Input            |Relevant problem type(s)| Details |
|------------------|------------------------|---------|
|Density           | All                    | Density of heterogenous nucleation sites in the liquid (evenly distributed among cells that are liquid or undergo melting), normalized by 1 x 10^12 m^-3
|MeanUndercooling  | All                    | Mean nucleation undercooling (relative to the alloy liquidus temperature) for activation of nucleation sites (Gaussian distribution)
|StDevUndercooling | All                    | Standard deviation of nucleation undercooling (Gaussian distribution), in K

## Temperature inputs
| Input        | Relevant problem type(s)| Details |
|--------------|-------------------------|---------|
|G             | C, S                    | Thermal gradient in the build (+Z) directions, in K/m
|R             | C, S                    | Cooling rate (uniform across the domain), in K/s
|HeatTransferCellSize |  R               | By default, equal to deltax, and cannot be used if remelting is considered. deltax must divide evenly into HTdeltax
|LayerwiseTempRead | R                   | If set to Y, the appropriate temperature data will be read during each layer's initialization, stored temporarily, and discarded. If set to N, temperature data for all layers will be read and stored during code initialization, and initialization of each layer will be performed using this stored temperature data. This option is only applicable to simulations with remelting; simulations without remelting (and simulations where this input is not given) default to N. Setting this to Y is only recommended if a large quantity of temperature data is read by ExaCA (for example, a 10 layer simulation where each layer's temperature data comes from a different file).
|TemperatureFiles | R                    | List of files corresponding to each layer's temperature data, in the form ["filename1.csv","filename2.csv",...]. If the number of entries is less than numberOfLayers, the list is repeated. Note that if the Z coordinate of the top surface for each data set has the layer offset applied, layerOffset in the "Domain" section of the input file should be set to 0, to avoid offsetting the layers twice.

## Substrate inputs
| Input        | Relevant problem type(s))| Details |
|--------------| -------------------------|---------|
|FractionSurfaceSitesActive | C           | What fraction of cells at the bottom surface of the domain are the source of a grain?
|MeanSize      | S, R                     | Mean spacing between grain centers in the baseplate/substrate (in microns) (see note (a))
|SubstrateFilename |  S, R                | Path to and filename for substrate data (see note (a))
|PowderDensity | S, R                     | Density of sites in the powder layer to be assigned as the home of a unique grain, normalized by 1 x 10^12 m^-3 (default value is 1/(CA cell size ^3) (see note (b))
|ExtendSubstrateThroughPower| S, R        | true/false value: Whether to use the baseplate microstructure as the boundary condition for the entire height of the simulation (defaults to false) (see note (b))

(a) One of these inputs must be provided, but not both
(b) This is optional, but if this is given, "extendSubstrateThroughPower" must be set to false

## Printing inputs
| Input        | Relevant problem type(s))| Details |
|--------------| -------------------------|---------|
| PathToOutput | All                      | File path location for the output files
| OutputFile   | All                      | All output files will begin with the string specified on this line
| PrintBinary  | All                      | Whether or not ExaCA vtk output data should be printed as big endian binary data, or as ASCII characters (defaults to false)
| PrintExaConstitSize | R                 | Length of the cubic representative volume element (RVE) data for ExaConstit, taken from the domain center in X and Y, and at the domain top in Z excluding the final layer's grain structure. If not given (or given a value of 0), the RVE will not be printed 
| PrintFieldsInit     | All, except with remelting | Which fields to print following initialization, in the form ["Field1","Field2",...]. Currently supported options are GrainID, LayerID, CellType, CritTimeStep, UndercoolingChange, UndercoolingCurrent
| PrintFieldsFinal    | All               | Which fields to print at the end of the run, in the form ["Field1","Field2",...]. Currently supported options are GrainID, LayerID, GrainMisorientation, UndercoolingCurrent
| PrintIntermediateOutput | All           | If this section is given, intermediate code output will be printed, where liquid CA cells are assigned the value -1, and non-liquid cells are assigned a value 0-62 depending on their associated grain's misorientation relative to the build (+Z) direction. This is an optional section that otherwise defaults to false, with second level inputs "Frequency" and "PrintIdleFrames"
| PrintIntermediateOutput:Frequency       | All | The number of microseconds defining the ExaCA output intermediate data increment. A value of 0 will mean that no intermediate output is printed
| PrintIntermediateOutput:PrintIdleFrames | All | Whether or not ExaCA should print intermediate output regardless of whether the simulation has changed from the last frame
