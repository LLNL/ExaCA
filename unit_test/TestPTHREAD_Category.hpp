#ifndef EXACA_TEST_PTHREAD_CATEGORY_HPP
#define EXACA_TEST_PTHREAD_CATEGORY_HPP

#define TEST_CATEGORY pthread
#define TEST_EXECSPACE Kokkos::Threads
#define TEST_MEMSPACE Kokkos::HostSpace
#define TEST_DEVICE Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>

#endif
