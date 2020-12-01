cmake_minimum_required(VERSION 3.9)
include(GNUInstallDirs)

#
# Find Kokkos and set compiler to Kokkos' nvcc_wrapper
#
# NOTE: I do this BEFORE calling 'project' so that you don't have to specify
# this compiler on the command line
#
list(APPEND CMAKE_PREFIX_PATH
    "${Kokkos_DIR}"
    "${Kokkos_DIR}/lib64/cmake/Kokkos"
)
set(Kokkos_ROOT ${Kokkos_DIR})
find_package(Kokkos 3 REQUIRED)
if (${Kokkos_ENABLE_CUDA})
    set(CMAKE_CXX_COMPILER "${Kokkos_ROOT}/bin/nvcc_wrapper")
endif()
set(ENV{KOKKOS_CMAKE_CXX_STANDARD} ${CMAKE_CXX_STANDARD})

project(ExaCA CXX)

find_package(MPI REQUIRED)

add_subdirectory(src)

include(CTest)
#add_subdirectory(unit_test)
if(${BUILD_TESTING})
    add_test(NAME ExaCA-Test COMMAND ExaCA-Kokkos)
endif()
