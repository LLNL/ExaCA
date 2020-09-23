cmake_minimum_required(VERSION 3.9)

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

set(EXACA_SOURCES
    CAfunctions.cpp
    CAghostnodes.cpp
    CAInitialize.cpp
    CAPrint.cpp
    CAupdate.cpp
    Run_Reduced.cpp
)

set(EXACA_HEADERS header.h)

add_library(ExaCA STATIC ${EXACA_SOURCES} ${EXACA_HEADERS})
target_link_libraries(ExaCA PUBLIC
    MPI::MPI_CXX
    Kokkos::kokkos
)

include(CTest)

if(${BUILD_TESTING})
    add_executable(ExaCA-Kokkos main.cpp)
    target_link_libraries(ExaCA-Kokkos ExaCA)
    add_test(NAME ExaCA-Test COMMAND ExaCA-Kokkos)
endif(${BUILD_TESTING})
