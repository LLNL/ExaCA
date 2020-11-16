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

### Build Kokkos

First, if Kokkos is not already built on your system, build Kokkos:
```
# Change this path to Kokkos source
cd ./kokkos
mkdir build
cd build
export KOKKOS_INSTALL_DIR=`pwd`/install
cmake \
  -D CMAKE_BUILD_TYPE="Release" \
  -D CMAKE_INSTALL_PREFIX=$KOKKOS_INSTALL_DIR \
  -D Kokkos_ENABLE_OPENMP=ON \
  \
  .. ;
make install
cd ../..
```
Note that ExaCA runs with the default enabled Kokkos backend
 (https://github.com/kokkos/kokkos/wiki/Initialization).

### Build ExaCA

Then build ExaCA, including the path to the Kokkos build:
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
  -D Kokkos_DIR=$KOKKOS_INSTALL_DIR \
  -D CMAKE_PREFIX_PATH=$EXACA_INSTALL_DIR \
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

Run by simply calling the created executable:
```
mpiexec -n 1 ./build/install/bin/ExaCA-Kokkos examples/Inp_DirSolidification.txt
```
