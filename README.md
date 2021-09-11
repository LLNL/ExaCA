# ExaCA
## An exascale-capable cellular automaton for nucleation and grain growth
ExaCA is a cellular automata (CA) code for grain growth under additive
manufacturing conditions by ExaAM within the Exascale Computing Project.

## Build
ExaCA-Kokkos uses Kokkos and MPI for parallelism.

### Dependencies

|Dependency | Version | Required | Details|
|---------- | ------- |--------  |------- |
|CMake      | 3.9+    | Yes       | Build system
|Kokkos     | 3.0+   | Yes      | Provides portable on-node parallelism.
|MPI        | GPU Aware if CUDA/HIP Enabled | Yes     | Message Passing Interface
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

ExaCA-Kokkos runs using an input file, passed on the command line. Example problems are provided in the `examples/` directory:

 * `Inp_DirSolidification.txt`: simulates grain growth from a surface with a fixed thermal gradient and cooling rate
 * `Inp_SmallDirSolidification.txt`: a smaller and simpler version of the previous
 * `Inp_SpotMelt.txt`: simulates overlapping spot melts with fixed a fixed thermal gradient and cooling rate
 * `Inp_SmallSpotMelt.txt`: a smaller and simpler version of the previous

Example problems only possible with external data:
 * `Inp_AMBenchMultilayer.txt`: simulates 4 layers of a representative even-odd layer alternating scan pattern for AM builds
 * `Inp_SimpleRaster.txt`: simulates a single layer consisting of four overlapping melt pools

Run by calling the created executable with an ExaCA input file:
```
mpiexec -n 1 ./build/install/bin/ExaCA-Kokkos examples/Inp_DirSolidification.txt
```
## Post-processing analysis

If the "Print Paraview vtk file of all ExaCA data for post-processing" options is toggled within an input file, post-processing can be performed on the output data set using the code inside of the `analysis/` directory. This code is compiled as a separate executable from, but linked to, ExaCA. Running ExaCA for the test problem  `Inp_DirSolidification.txt` yields the output files `TestProblemDirS.vtk` and `TestProblemDirS.log`. To analyze this data, first compile the source files, then run the analysis program executable, with one command line argument pointing to the analysis file of interest. Within the `analysis/examples` directory, there should be a `.txt` file with analysis options corresponding to the problem of interest. In this case, to analyze the results of the directional solidification test problem:

```
./build/install/bin/grain_analysis analysis/examples/AnalyzeDirS.txt
```
The path to the associated files `TestProblemDirS.vtk` and `TestProblemDirS.log` is an input within `TestProblemDirS.txt`.

## Contributing

We encourage you to contribute to ExaCA! Please check the
[contribution guidelines](CONTRIBUTING.md).

## License

ExaCA is distributed under an [MIT license](LICENSE).

## Release

LLNL-CODE-821827
