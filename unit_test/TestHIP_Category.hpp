#ifndef EXACA_TEST_HIP_CATEGORY_HPP
#define EXACA_TEST_HIP_CATEGORY_HPP

#define TEST_CATEGORY hip
#define TEST_EXECSPACE Kokkos::Experimental::HIP
#define TEST_MEMSPACE Kokkos::Experimental::HIPSpace
#define TEST_DEVICE Kokkos::Device<Kokkos::Experimental::HIP, Kokkos::Experimental::HIPSpace>

#endif
