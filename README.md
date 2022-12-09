# ExaCA
## An exascale-capable cellular automaton for nucleation and grain growth
ExaCA is a cellular automata (CA) code for grain growth under additive
manufacturing conditions by ExaAM within the Exascale Computing Project.

## Build
ExaCA-Kokkos uses Kokkos and MPI for parallelism.

### Dependencies

|Dependency | Version | Required | Details|
|---------- | ------- |--------  |------- |
|CMake      | 3.11+   | Yes      | Build system
|Kokkos     | 3.0+   | Yes      | Provides portable on-node parallelism.
|MPI        | GPU Aware if CUDA/HIP Enabled | Yes     | Message Passing Interface
|[nlohmann_json](https://github.com/nlohmann/json) | 3.10+| No       | Input parsing (Note: this will become a required dependency in a future release)
|GoogleTest | 1.10+   | No       | Unit test framework
|CUDA       | 9+      | No       | Programming model for NVIDIA GPUs
|HIP       | 3.5+      | No       | Programming model for AMD GPUs

Kokkos and MPI are available on many systems; if not, obtain the desired
versions:
```
git clone https://github.com/kokkos/kokkos.git --branch 3.4.00
```

### Backend options
Note that ExaCA runs with the default enabled Kokkos backend
 (see https://github.com/kokkos/kokkos/wiki/Initialization).

ExaCA has been tested with Serial, OpenMP, Pthreads, CUDA, and HIP backends.

### Build CPU

First, if Kokkos is not already built on your system, build Kokkos:
```
# Change this path to Kokkos source
cd ./kokkos
mkdir build
cd build
# Check the CPU architecture flag matches the hardware
cmake \
  -D CMAKE_BUILD_TYPE="Release" \
  -D CMAKE_INSTALL_PREFIX=install \
  -D Kokkos_ENABLE_OPENMP=ON \
  -D Kokkos_ARCH_POWER9=ON \
  .. ;
make install
cd ../..
```
Note that there are other host backends available. The Kokkos architecture flag must match the hardware you run on and will improve performance, if used.

Then build ExaCA, including the path to the Kokkos build:
```
# Change this path to Kokkos installation
export KOKKOS_INSTALL_DIR=./kokkos/build/install

# Change this path to ExaCA source
cd ./ExaCA
mkdir build
cd build
cmake \
  -D CMAKE_BUILD_TYPE="Release" \
  -D CMAKE_PREFIX_PATH=$KOKKOS_INSTALL_DIR \
  -D CMAKE_INSTALL_PREFIX=install \
  ..;
make install
cd ../..
```

One additional CMake option can be set to enable JSON input parsing:
`ExaCA_ENABLE_JSON`. The library will be automatically downloaded if requested,
but not found. Note that this option will be removed in a future release when
plain text input files are removed (and JSON is required).

### Build CUDA

If running on NVIDIA GPUs, build Kokkos with additional inputs:
```
# Change this path to Kokkos source
cd ./kokkos
mkdir build
cd build
# Check the GPU architecture flag matches the hardware
cmake \
  -D CMAKE_BUILD_TYPE="Release" \
  -D CMAKE_CXX_COMPILER=../bin/nvcc_wrapper \
  -D CMAKE_INSTALL_PREFIX=install \
  -D Kokkos_ENABLE_CUDA=ON \
  -D Kokkos_ENABLE_CUDA_LAMBDA=ON \
  -D Kokkos_ARCH_VOLTA70=ON \
  .. ;
make install
cd ../..
```
Note the two flags needed for the `Kokkos::Cuda` backend. The Kokkos architecture flag must match the hardware you run on and will improve performance. By default, the host will use `Kokkos::Serial`; other parallel host backends can also be used, e.g. by adding `-D Kokkos_ENABLE_OPENMP`.

Build ExaCA, this time with the Kokkos compiler wrapper:
```
# Change this path to Kokkos installation
export KOKKOS_INSTALL_DIR=./kokkos/build/install

# Change this path to ExaCA source
cd ./ExaCA
mkdir build
cd build
export EXACA_INSTALL_DIR=`pwd`/install
cmake \
  -D CMAKE_BUILD_TYPE="Release" \
  -D CMAKE_CXX_COMPILER=$KOKKOS_INSTALL_DIR/bin/nvcc_wrapper \
  -D CMAKE_PREFIX_PATH=$KOKKOS_INSTALL_DIR \
  -D CMAKE_INSTALL_PREFIX=install \
  ..;
make install
```

### Build HIP
Again, first build Kokkos, this time with the `hipcc` compiler:
```
cd ./kokkos
mkdir build
cd build
cmake \
    -D CMAKE_BUILD_TYPE="Release" \
    -D CMAKE_CXX_COMPILER=hipcc \
    -D CMAKE_INSTALL_PREFIX=install \
    -D Kokkos_ENABLE_HIP=ON \
    -D Kokkos_ARCH_VEGA908=ON \
    .. ;
make install
```
And build ExaCA with the same compiler:
```
# Change this path to Kokkos installation
export KOKKOS_INSTALL_DIR=./kokkos/build/install

cd ./ExaCA
mkdir build
cd build
cmake \
    -D CMAKE_BUILD_TYPE="Release" \
    -D CMAKE_CXX_COMPILER=hipcc \
    -D CMAKE_PREFIX_PATH="$KOKKOS_INSTALL" \
    -D CMAKE_INSTALL_PREFIX=install \
    .. ;
make install
```

## Run

ExaCA-Kokkos runs using an input file, passed on the command line. Example problems are provided in the `examples/` directory - A separate README file located in the `examples/` directory goes into more detail on the problem types, the optional and required arguments needed for each problem type, and additional files used by ExaCA. The example input files present in this repository are:
 * `Inp_DirSolidification.txt`: simulates grain growth from a surface with a fixed thermal gradient and cooling rate
 * `Inp_SmallDirSolidification.txt`: a smaller and simpler version of the previous
 * `Inp_SpotMelt.txt`: simulates overlapping spot melts with fixed a fixed thermal gradient and cooling rate, where cells are only allowed to undergo solidification one time (i.e., overlap regions only solidify when they've cooled below the liquidus for the final time in the simulation)
 * `Inp_SmallSpotMelt.txt`: a smaller and simpler version of the previous
 * `Inp_SpotMelt_RM.txt`: simulates overlapping spot melts with fixed a fixed thermal gradient and cooling rate, where cells in the overlap region are allowed to melt and solidify as many times as needed
 * `Inp_SmallSpotMelt_RM.txt`: a smaller and simpler version of the previous

Example problems only possible with external data (available via https://github.com/LLNL/ExaCA-Data):
 * `Inp_AMBenchMultilayer.txt`: simulates 4 layers of a representative even-odd layer alternating scan pattern for AM builds
 * `Inp_SimpleRaster.txt`: simulates a single layer consisting of four overlapping melt pools
 * `Inp_TwoLineTwoLayer.txt`: simulates two layers consisting of segments of two overlapping melt pools

Run by calling the created executable with an ExaCA input file:
```
mpiexec -n 1 ./build/install/bin/ExaCA-Kokkos examples/Inp_DirSolidification.txt
```
## Automated input file generation using Tasmanian (https://tasmanian.ornl.gov/)
Within the `utilities` directory, an example python script for the generation of an ensemble of input files is available. By running the example script `TasmanianTest.py`, 69 ExaCA input files are generated with a range of heterogenous nucleation density, mean nucleation undercooling, and mean substrate grain size values, based on the ranges in python code (N0Min-N0Max, dTNMin-dTNMax, and S0Min-S0Max), respectively. Running the python script from the ExaCA source directory, via the command
```
python utilities/TasmanianTest.py TemperatureData 2
```
the script will generate an ensemble of input files in the `examples` directory, with "Temperature filename(s): TemperatureData" and "Number of temperature files: 2" based on the command line inputs to the python script. The output file name for ExaCA simulations run using the input file `examples/Inp_TemperatureDataEnsembleMember($N)`, generated as one of N = 1 to 69 input files that resulted from running the example python script, will be `TemperatureData_ExaCAEnsMem_($N)`. Other CA inputs, such as the path to temperature data, or the time step) must be adjusted manually inside of the python script. Separate instances of ExaCA can be run with each ensemble member to probe microstructure dependency on nucleation and substrate.

## Output and post-processing analysis

If the "Print file of grain misorientations" option is turned on within an input file, ExaCA will output a scalar field "Angle_z" as a vtk file ending with "Misorientations.vtk". Angle_z corresponds to the orientation (in degrees) of a given grain relative to the positive Z direction in a simulation (the thermal gradient direction for directional solidification problems, the build/layer offset direction for other problems). Epitaxial grains (from the initial grain structure or powder layer) are assigned values between 0 and 62.7, while nucleated grains (not present in the initial grain structure) are assigned values between 100 and 162.7 (the offset of 100 is simply used to ensure the two types of grains are differentiated, but a nucleated grain with Angle_z = 135 actually has a misorientation of 35 degrees).

If the "Print Paraview vtk file" option is turned on within an input file, post-processing can be performed on the output data set. This functionality is a separate executable from ExaCA, located in the `analysis/` directory and is linked to the ExaCA library for input utilities. 

Specifying debug check options can be done to print various ExaCA data fields to files following simulation initialization. The "reduced" debug check will print "CritTimeStep" (the time step at which each cell goes below the liquidus for the final time), "LayerID" (the layer associated with each cell going below the liquidus for the final time, with layer 0 being the first layer, and -1 for all cells that did not undergo solidification) and "CellType" (integers corrsponding to cell types specified in src/CAtypes.hpp). The "extensive" debug check will, in addition to the "reduced" data fields, also print "UndercoolingChange" (the rate at which a cell cools per time step after reaching its "CritTimeStep" value), "UndercoolingCurrent" (the initial undercooling of each cell), and "Melted" (1 for cells that are part of the melt pool, 0 for cells that are not).

ExaCA can optionally print the system state at intermediate time values as part of a series of vtk files that can be read by Paraview to make animations, if the "Print intermediate output frames" option is turned on. "Increment to separate frames" is the separation between intermediate output files in microseconds - if there is a long time period between solidification events (such as two overlapping melt pools formed via line scan with a long dwell time between them), setting "Intermediate output even if system is unchanged from previous state" to off will skip printing of those files.

Running ExaCA for the test problem `Inp_DirSolidification.txt` yields the output files `TestProblemDirS.vtk` (containing LayerID, GrainID, and Melted data) and `TestProblemDirS.log` (containing information regarding the simulation parameters used, simulation dimensions, and some timing data). To analyze this data, run `grain_analysis` (installed in the same location as `ExaCA-Kokkos`), with one command line argument pointing to the analysis input file. Within the `analysis/examples` directory, there are example analysis input files:

```
./build/install/bin/grain_analysis analysis/examples/AnalyzeDirS.txt
```
Note that the path to the files needed for analysis, e.g. `TestProblemDirS.vtk` and `TestProblemDirS.log`, are configurable inputs within the analysis input file.

## Citing ExaCA

If you use ExaCA in your work, please cite the following paper:
```
@article{rolchigo2022,
  title = {ExaCA: A performance portable exascale cellular automata application for alloy solidification modeling},
  journal = {Computational Materials Science},
  volume = {214},
  pages = {111692},
  year = {2022},
  issn = {0927-0256},
  doi = {https://doi.org/10.1016/j.commatsci.2022.111692},
  url = {https://www.sciencedirect.com/science/article/pii/S0927025622004189},
  author = {Matt Rolchigo and Samuel Temple Reeve and Benjamin Stump and Gerald L. Knapp and John Coleman and Alex Plotkowski and James Belak},
}
```

If you would like to cite the software itself, cite the current release or
version used from [Zenodo](https://zenodo.org/record/7126168#.Yzs_TCHMJOU).

## Contributing

We encourage you to contribute to ExaCA. Please check the
[contribution guidelines](CONTRIBUTING.md).

## License

ExaCA is distributed under an [MIT license](LICENSE).

## Release

LLNL-CODE-821827
