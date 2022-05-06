cmake_minimum_required(VERSION 3.9)

project(ExaCA CXX)

include(GNUInstallDirs)

find_package(Kokkos 3 REQUIRED)
find_package(MPI REQUIRED)

# NOTE: this order is specifically used in the unit tests.
set(EXACA_SUPPORTED_DEVICES SERIAL PTHREAD OPENMP CUDA HIP)

add_subdirectory(src)
add_subdirectory(analysis)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/ExaCAconfig.cmakein
  ${CMAKE_CURRENT_BINARY_DIR}/ExaCAconfig.cmake @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/ExaCAconfig.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ExaCA)

# Install data files
file(GLOB MATERIALS examples/Materials/Inconel625)
file(GLOB SUBSTRATE examples/Substrate/GrainOrientationVectors*.csv)
install(FILES ${MATERIALS}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/ExaCA/Materials)
install(FILES ${SUBSTRATE}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/ExaCA/Substrate)

option(ExaCA_ENABLE_TESTING "Build unit tests" OFF)

if(ExaCA_ENABLE_TESTING)
  find_package(GTest 1.10 REQUIRED)
  # Workaround for FindGTest module in CMake older than 3.20
  if(TARGET GTest::gtest)
    set(gtest_target GTest::gtest)
  elseif(TARGET GTest::GTest)
    set(gtest_target GTest::GTest)
  else()
    message(FATAL_ERROR "bug in GTest find module workaround")
  endif()
  enable_testing()
  add_subdirectory(unit_test)
endif()

# clang-format
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package( CLANG_FORMAT 10)
if(CLANG_FORMAT_FOUND)
  file(GLOB_RECURSE FORMAT_SOURCES src/*.[c,h]pp analysis/*.[c,h]pp unit_test/*.[c,h]pp)
  add_custom_target(format
    COMMAND ${CLANG_FORMAT_EXECUTABLE} -i -style=file ${FORMAT_SOURCES}
    DEPENDS ${FORMAT_SOURCES})
endif()
