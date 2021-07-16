#ifndef EXACA_TEST_CUDAUVM_CATEGORY_HPP
#define EXACA_TEST_CUDAUVM_CATEGORY_HPP

#include <Kokkos_Cuda.hpp>

#include <gtest/gtest.h>

#define TEST_CATEGORY cuda_uvm
#define TEST_EXECSPACE Kokkos::Cuda
#define TEST_MEMSPACE Kokkos::CudaUVMSpace
#define TEST_DEVICE Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>

#endif
