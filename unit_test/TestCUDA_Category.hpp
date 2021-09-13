#ifndef EXACA_TEST_CUDA_CATEGORY_HPP
#define EXACA_TEST_CUDA_CATEGORY_HPP

#define TEST_CATEGORY cuda
#define TEST_EXECSPACE Kokkos::Cuda
#define TEST_MEMSPACE Kokkos::CudaSpace
#define TEST_DEVICE Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>

#endif
