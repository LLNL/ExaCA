# ExaCA-Kokkos
ExaCA is a cellular automata (CA) code for grain growth under additive
manufacturing conditions by ExaAM within the Exascale Computing Project.

## Build
ExaCA-Kokkos uses Kokkos and MPI for parallelism.

### Dependencies

|Dependency | Version | Required | Details|
|---------- | ------- |--------  |------- |
|CMake      | 3.9+    | Yes       | Build system
|Kokkos     | 3.0+   | Yes      | Provides portable on-node parallelism.
|MPI        | GPU Aware if CUDA Enabled | Yes     | Message Passing Interface
|CUDA       | 9+      | No       | Programming model for NVIDIA GPUs

Kokkos and MPI are available on many systems; if not, obtain the desired
versions:
```
git clone https://github.com/kokkos/kokkos.git --branch 3.0.00
```

### Backend options
Note that ExaCA runs with the default enabled Kokkos backend
 (see https://github.com/kokkos/kokkos/wiki/Initialization).

ExaCA has been tested with Serial, OpenMP, Pthreads, and Cuda backends.

### Build CPU

First, if Kokkos is not already built on your system, build Kokkos:
```
# Change this path to Kokkos source
cd ./kokkos
mkdir build
cd build
cmake \
  -D CMAKE_BUILD_TYPE="Release" \
  -D CMAKE_INSTALL_PREFIX=install \
  -D Kokkos_ENABLE_OPENMP=ON \
  \
  .. ;
make install
cd ../..
```

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
  \
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
cmake \
  -D CMAKE_BUILD_TYPE="Release" \
  -D CMAKE_CXX_COMPILER=../bin/nvcc_wrapper \
  -D CMAKE_INSTALL_PREFIX=install
  -D Kokkos_ENABLE_CUDA=ON \
  \
  .. ;
make install
cd ../..
```

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
  \
  ..;
make install
cd ../..
```

## Run

ExaCA-Kokkos runs using an input file, passed on the command line. Three examples problems are given in the `examples/` directory:

 * `Inp_AMBenchMultilayer.txt` simulates 4 layers of a representative even-odd layer alternating scan pattern for AM builds
 * `Inp_SimpleRaster.txt` simulates a single layer consisting of four overlapping melt pools
 * `Inp_DirSolidification.txt` does not use a thermal profile for a beam melting problem, but rather simulates grain growth from a surface with a fixed thermal gradient and cooling rate
 * `Inp_SmallDirSolidification.txt` the smallest and simplest example problem

Run by calling the created executable from the ExaCA root directory:
```
mpiexec -n 1 ./build/install/bin/ExaCA-Kokkos examples/Inp_DirSolidification.txt
```
