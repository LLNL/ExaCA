# ExaCA-Kokkos
ExaCA is a cellular automata (CA) code for grain growth under additive
manufacturing conditions by ExaAM within the Exascale Computing Project.

## Build
ExaCA-Kokkos uses GNU make and has minimal dependencies.

|Dependency | Version | Required | Details|
|---------- | ------- |--------  |------- |
|Kokkos     | 2.9.0   | Yes      | Provides portable on-node parallelism.
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
By default, ExaCA-Kokkos runs with settings in `examples/MasterInputs.txt`.

Run by simply calling the created executable:
```
mpiexec -n 1 ./ExaCA-Kokkos.cuda
```
