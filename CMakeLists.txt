cmake_minimum_required(VERSION 3.9)
include(GNUInstallDirs)

project(ExaCA CXX)

find_package(Kokkos 3 REQUIRED)
find_package(MPI REQUIRED)

add_subdirectory(src)

include(CTest)
#add_subdirectory(unit_test)
if(${BUILD_TESTING})
    add_test(NAME ExaCA-Test COMMAND ExaCA-Kokkos)
endif()
