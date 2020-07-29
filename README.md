# ExaCA-Kokkos
ExaCA is a cellular automata (CA) code for grain growth under additive
manufacturing (AM) conditions by ExaAM within the Exascale Computing Project.

## Build
ExaCA-Kokkos uses GNU make and has minimal dependencies.

|Dependency | Version | Required | Details|
|---------- | ------- |--------  |------- |
|Kokkos     | 2.7.0+  | Yes      | Provides portable on-node parallelism
|MPI        | GPU Aware if CUDA Enabled | Yes     | Message Passing Interface
|CUDA       | 9+      | No       | Programming model for NVIDIA GPUs

Kokkos and MPI are available on many systems; if not, obtain the desired
versions:
```
git clone https://github.com/kokkos/kokkos.git --branch 2.9.00
```

Ensure `Makefile` settings match desired build, e.g. for CUDA Volta GPUs:
```
KOKKOS_PATH = ${HOME}/kokkos
KOKKOS_DEVICES = "Cuda"
KOKKOS_ARCH = "Volta70"
```
and run `make` to build ExaCA-Kokkos and Kokkos together.


## Run
ExaCA-Kokkos runs using an input file, the name of which is set by the
environmental variable `CAINPUT` either from the terminal or inside a job
script. Three examples problems are given in the `Examples` directory:

* `Inp_AMBenchMultilayer.txt` simulates 4 layers of a representative even-odd
   layer alternating scan pattern for AM builds
* `Inp_SimpleRaster.txt` simulates a single layer consisting of four
   overlapping melt pools
* `Inp_DirSolidification.txt` does not use a thermal profile for a beam
   melting problem, but rather simulates grain growth from a surface with a
   fixed thermal gradient and cooling rate

Run by choosing an input file and calling the executable:
```
export CAINPUT=Inp_DirSolidification.txt
mpiexec -n 1 ./ExaCA-Kokkos.cuda
```
