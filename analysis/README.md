# ExaCA grain analysis post-processing

Running the `grain_analysis` executable after generation of an ExaCA data set (i.e., the associated .vtk file that results from setting `Print Paraview vtk file: Y` in the ExaCA input file, and the associated .log file) allows additional insights into the microstructure generated. When running the analysis script, the path to and name of an associated analysis inputs file should be given on the command line: for example:

```
./build/install/bin/grain_analysis analysis/examples/AnalyzeDirS.txt
```

The analysis files, such as `examples/AnalyzeDirS.txt`, have multiple distinct sections, starting with lines consisting of 5 asterisks, and ending with lines consisting of 10 asterisks. Everything above the first line with an asterisk is part of a header, and ignored by the program.

## Section 0: Required inputs

This is where the required inputs for analysis are given; these inputs are:

* Path to/name of log file associated with the microstructure of interest (no file extension should be given; .json assumed, use of .log is deprecated)
* Path to/name of microstructure (.vtk) file
* Path to/base (without file extension) name of data files of output resulting from this analysis

Three deprecated inputs may also be given here, though these will be ignored in a future release; these inputs are:
* Path to file of grain orientations (rotation matrix form)
* Path to file of grain orientations (Bunge Euler angle form ZXZ)
* Path to file corresponding to grain orientation IPF-Z colors as fractional RGB values

Values for these inputs should be separated by a colon from the rest of the line. The paths to the 3 grain orientation files will only be be used if they cannot be extracted from the log file. The grain orientation files "examples/Substrate/GrainOrientationVectors.csv", "examples/Substrate/GrainOrientationEulerAnglesBungeZXZ.csv", and "examples/Substrate/GrainOrientationRGB_IPF-Z.csv" are used by the example problems. If a custom file of grain orientation vectors was used by ExaCA to simulate a microstructure (i.e., "path/to/file/GrainOrientationVectors_customname.csv"), this requires that the corresponding 2 orientation files are also provided ("path/to/file/GrainOrientationEulerAnglesBungeZXZ_customname.csv" and "path/to/file/GrainOrientationRGB_IPF-Z_customname.csv") and that they are properly formatted (see `examples/README.md` for more details).
 If the orientation file data is not given in the log file, and "Path to file of grain orientations (Bunge Euler angle form ZXZ)" is not given in the analysis inputs file, the deprecated form of printing ExaCA cross-section data is assumed. If "Path to file corresponding to grain orientation IPF-Z colors as fractional RGB values" is not given, the default file (examples/Substrate/GrainOrientationRGB_IPF-Z.csv) will be used to map orientations to RGB values. In the future these paths will be required and missing input will result in a runtime error.

## Section 1: ExaConstit representative volume element (RVE) printing

Representative volume elements can be carved out from the ExaCA results and used as input to the crystal plasticity FEM code ExaConstit (https://github.com/LLNL/ExaConstit). Lines in this section have the format [size, center x, center y, center z] - multiple RVEs can be printed as long as they are formatted as above and on separate lines.
* Size of the RVE domain is given in cells (must be less than domain dimensions)
* Center x, center y, and center z are the location (in cells) of the center of this RVE (If "D" is given as any of the inputs, the RVE will be taken to be as close to the domain center in x and y as possible, and as close to the top of the domain in z as possible while avoiding the last layer's microstructure)
* If no ExaConstit data should be printed, no lines should exist between the short and long lines of asterisks here

## Section 2: Data for plotting microstructure using MTEX

Cross-sections of the ExaCA microstructure can be specified for printing in a form usable by MTEX, for plotting inverse pole figure-coloring of the grain structure and/or pole figures. Each line should consist of 4 comma-separated values:
* Value 1 is the plane of interest, should be XY, YZ, or XZ
* Value 2 is the out of plane coordinate of the plane of interest (for example, "XY, 100" would be the XY plane located at Z = 100, in CA cells)
  - If "D" is given as the out of plane coordinate, the plane that bisects the microstructure is chosen
  - If "END" is given as the out of plane coordinate, the largest possible out of plane coordinate is used ("XY, END" gives the top surface)
* Value 3 is either Y or N, whether or not pole figure data should be printed for this cross-section
* Value 4 is either Y or N, whether or not EBSD-like inverse pole figure coloring data should be printed for this cross-section

If no cross-section data should be printed, no lines should exist between the short and long lines of asterisks. It should also be noted that a deprecated input form, consisting of square brackets enclosing "Value1" and "Value2" from above, separated by commas, is also accepted though it will be removed in a future release. Use of this deprecated input form will print data as a series of "X coordinate, Y coordinate, GrainID" values: additional post-processing would be needed to convert this into a usable form by MTEX.

## Section 3: Y/N Options for data analysis of representative regions

The first 3 lines in this section allow specification of a volume of the microstructure data for additional analysis. These should take the form:

X Bounds (in cells): [XMin, XMax]
Y Bounds (in cells): [YMin, YMax]
Z Bounds (in cells): [ZMin, ZMax]

The specified bounds XMin, XMax, YMin, YMax, ZMin, ZMax should not exceed the bounds of the data set itself. Putting `D` as the value for XMin, YMin, XMax, or YMax will start the volume 50 cells from that edge of the data set. Putting `D` as ZLow will attempt to start from the lowest Z that doesn't include substrate based on the XY size. Putting `D` as ZHigh will attempt to end at the Z coordinate as close to the top of the domain in z as possible while avoiding the last layer's microstructure. If either of these aren't possible, it will default to using the entirety of the domain in Z)

The remaining lines are yes/no options for calculating various aspects of the microstructure. Regardless of toggled options, code will do some (relatively quick) default analysis on the XY cross-section as close to the top of the domain in z as possible while avoiding the last layer's microstructure, as well as on the representative region (spanning the domain in Z) as a whole. If the final option, "Print volume for pole figure analysis in MTEX", is marked with a `Y`, the pole figure data will either be printed in the form directly usable by MTEX as in Section 2 (if "Path to file of grain orientations (Bunge Euler angle form ZXZ)" was given in Section 0), or it will be printed as a histogram of orientations. The number of lines in the histogram file is the number of orientations given in the specified file "Path to file of grain orientations (rotation matrix form)" from Section 0, and the values are the number of cells associated with that orientation appear in the volume specified.

## MATLAB post-processing
The pole figure data and inverse pole figure coloring data generated from Section 2 cross-sections or the volume in Section 3 can be plotted using Matlab and the MTEX toolbox (https://mtex-toolbox.github.io/) using the scripts and colormaps located in the `utilities/MTEX` folder. Those files are as follows:

* `PlotPoleFigure.m` takes an appropriate input file from the analysis post-processing script, a colormap file, and an appropriate upper limit for the colorman (in multiples of uniform distribution) and plots the 100, 110, and 111 pole figures as well as the X, Y, and Z inverse pole figures
* `PlotIPFColoredSection.m` takes an appropriate input file from the analysis post-processing script (note that the format is not the same as the input files for plotting pole figures), smooths the grain structure data, and plots the microstructure colored using the inverse pole figure X, Y, and Z maps (cubic crystal geometry)
* `ColorMap4.m` and `ColorMap8.m` are example colormap files to be used as an input for `PlotPoleFigure.m` (note that the appropriate corresponding upper limits of 4 and 8, respectively, should be used as inputs to `PlotPoleFigure.m` as well)

Three MTEX files resulting from running the ExaCA and analysis executables for the directional solidification example problem can be accessed at https://github.com/LLNL/ExaCA-Data and used as test input in these post-processing scripts.
