#ifndef EXACA_TEST_THREADS_CATEGORY_HPP
#define EXACA_TEST_THREADS_CATEGORY_HPP

#define TEST_CATEGORY threads
#define TEST_EXECSPACE Kokkos::Threads
#define TEST_MEMSPACE Kokkos::HostSpace
#define TEST_DEVICE Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>

#endif
