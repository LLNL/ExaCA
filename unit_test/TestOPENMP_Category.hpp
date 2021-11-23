#ifndef EXACA_TEST_OPENMP_CATEGORY_HPP
#define EXACA_TEST_OPENMP_CATEGORY_HPP

#define TEST_CATEGORY openmp
#define TEST_EXECSPACE Kokkos::OpenMP
#define TEST_MEMSPACE Kokkos::HostSpace
#define TEST_DEVICE Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>

#endif
