# ExaCA grain analysis post-processing

## Running the grain_analysis executable
Running the `grain_analysis` executable after generation of an ExaCA data set (i.e., a .vtk file that at minimum contains the field `GrainID`, and the associated .json file of log info) allows additional insights into the microstructure generated. When running the analysis script, the path to and name of an associated analysis inputs file, as well as the path to and base name of the microstructure to be analyzed, should be given on the command line. For example:

```
./build/install/bin/ExaCA-GrainAnalysis analysis/examples/AnalyzeDirS.txt TestProblemDirS
```

The analysis files, such as `examples/AnalyzeDirS.json`, allow analysis of multiple regions of the microstructure. Under the top level header `Regions`, the user-specified names of the regions for analysis are given (in the case of `examples/AnalyzeDirS.json`, these are `RepresentativeVolume`, `XYCross`, and `YZCross`). At the third level, analysis options for each individual region are given.

1D, 2D, and 3D regions can be analyzed depending on the bounds provided. The `units` input for each region should be either Meters or Cells, depending on the values used as the bounds. For example, if `xBounds`, `yBounds`, and `zBounds` are all provided (in the form `[lower, upper]`), and the lower and upper bounds for each direction do not match, the region is analyzed as a volume. If one direction has equivalent lower and upper bounds, the region is analyzed as an area, and if two directions have equivalent lower and upper bounds, the region is analyzed as a line. If no bounds are given for a direction, the analysis will default to including the entire domain extent in said direction.

Once the bounds of each region are identified, the analysis options can be specified within JSON arrays under `printAvgStats`, `printPerGrainStats`, and `printPerZCoordinateStats`. 
* Values specified in `printStats` will print average quantities to the console as well as to a file names `[MicrostructureBaseFilename]_[RegionName]_QoIs.txt`
* Values specified in `printPerGrainStats` will print data for each individual grain ID to a csv file of per-grain data named `[MicrostructureBaseFilename]_[RegionName]_grains.csv`
* Values specified in `printPerZCoordinateStats` will be printed in additional separate files.

| Output                | Compatible options            | Details
|=======================|===============================|====================
| GrainTypeFractions    | printAvgStats                    | Prints the fraction of the region consisting of nucleated grains, and the fraction that did not undergo melting
| Misorientation        | printAvgStats/printPerGrainStats | Prints the misorientation of the grain's <001> directions with the cardinal directions
| Size                  | printAvgStats/printPerGrainStats | Prints the grain size (length in microns, area in square microns, or volume in cubic microns, depending on the dimensionality of the region)
| BuildTransAspectRatio | printAvgStats/printPerGrainStats | Prints the grain aspect ratio, calculated as the extent of the grain in the build (Z) direction divided by the average of the extents in the X and Y directions
| XExtent               | printAvgStats/printPerGrainStats | Prints the extent of the grain (in microns) in the X direction
| YExtent               | printAvgStats/printPerGrainStats | Prints the extent of the grain (in microns) in the Y direction
| ZExtent               | printAvgStats/printPerGrainStats | Prints the extent of the grain (in microns) in the Z direction
| IPFZ-RGB              | printPerGrainStats            | Prints the R,G,B values (each between 0 and 1) corresponding to the inverse pole figure mapping of the grain orientation, relative to the build (Z) direction
| MeanGrainArea         | printPerZCoordinateStats            | Prints the mean grain cross-sectional (XY) area for each Z coordinate in the simulation to a file `[MicrostructureBaseFilename]_GrainAreas.csv`
| MeanWeightedGrainArea | printPerZCoordinateStats            | Prints the mean grain cross-sectional (XY) area, weighted by the grain area itself, for each 5th Z coordinate in the simulation to a file `[MicrostructureBaseFilename]_WeightedGrainAreas.csv` (Note: this option will be removed in a future release)

Additional analysis options for certain region types can be specified by setting them to `true` in the analysis input file (if the option does not appear in the input file, it is turned off by default)

| Output                | Compatible options            | Details
|=======================|===============================|====================
| PrintExaConstitYN           | volume                | Prints the grain ID data into an RVE usable by ExaConstit for constitutive properties simulation
| PrintPoleFigureYN           | area/volume           | Prints the grain orientation frequency data to a file which can be further analyzed in Matlab to generate pole figure and inverse pole figure data using the MTEX library
| PrintInversePoleFigureMapYN | area                  | Prints the grain euler angle data as a function of location in the cross-section to a file which can be further analyzed in Matlab to map the orientations to inverse pole figure-colored values (EBSD-like) using the MTEX library

## MATLAB post-processing
The pole figure data and inverse pole figure coloring data generated from Section 2 cross-sections or the volume in Section 3 can be plotted using Matlab and the MTEX toolbox (https://mtex-toolbox.github.io/) using the scripts and colormaps located in the `utilities/MTEX` folder. Those files are as follows:

* `PlotPoleFigure.m` takes an appropriate input file from the analysis post-processing script, a colormap file, and an appropriate upper limit for the colorman (in multiples of uniform distribution) and plots the 100, 110, and 111 pole figures as well as the X, Y, and Z inverse pole figures
* `PlotIPFColoredSection.m` takes an appropriate input file from the analysis post-processing script (note that the format is not the same as the input files for plotting pole figures), smooths the grain structure data, and plots the microstructure colored using the inverse pole figure X, Y, and Z maps (cubic crystal geometry)
* `ColorMap4.m` and `ColorMap8.m` are example colormap files to be used as an input for `PlotPoleFigure.m` (note that the appropriate corresponding upper limits of 4 and 8, respectively, should be used as inputs to `PlotPoleFigure.m` as well)

Three MTEX files resulting from running the ExaCA and analysis executables for the directional solidification example problem can be accessed at https://github.com/LLNL/ExaCA-Data and used as test input in these post-processing scripts.
