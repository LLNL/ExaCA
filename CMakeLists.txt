cmake_minimum_required(VERSION 3.9)

project(ExaCA CXX)

include(GNUInstallDirs)

find_package(Kokkos 3 REQUIRED)
find_package(MPI REQUIRED)

set(EXACA_SUPPORTED_DEVICES SERIAL PTHREAD OPENMP CUDA)

add_subdirectory(src)

if(ExaCA_ENABLE_TESTING)
  add_subdirectory(unit_test)
  enable_testing()
endif()

# clang-format
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package( CLANG_FORMAT 10)
if(CLANG_FORMAT_FOUND)
  file(GLOB_RECURSE FORMAT_SOURCES src/*.cpp src/*.hpp unit_test/*.cpp unit_test/*.hpp)
  add_custom_target(format
    COMMAND ${CLANG_FORMAT_EXECUTABLE} -i -style=file ${FORMAT_SOURCES}
    DEPENDS ${FORMAT_SOURCES})
endif()
