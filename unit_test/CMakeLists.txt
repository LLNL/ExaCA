##---------------------------------------------------------------------------##
# Download and unpack googletest
##---------------------------------------------------------------------------##
set(GTEST_URL "https://github.com/google/googletest/archive/release-1.10.0.tar.gz"
    CACHE STRING "URL for GTest tarball")
mark_as_advanced(GTEST_URL)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL ${GTEST_URL}
  URL_MD5         ecd1fa65e7de707cd5c00bdac56022cd
  )

FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
  FetchContent_Populate(googletest)
  add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif()

# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

set(gtest_args --gtest_color=yes)

##--------------------------------------------------------------------------##
## Add backends
##--------------------------------------------------------------------------##
set(EXACA_TEST_DEVICES)
foreach(_device ${EXACA_SUPPORTED_DEVICES})
  if(Kokkos_ENABLE_${_device})
    list(APPEND EXACA_TEST_DEVICES ${_device})
  endif()
endforeach()

##--------------------------------------------------------------------------##
## Create tests without Kokkos
##--------------------------------------------------------------------------##
macro(ExaCA_add_tests_nobackend)
  cmake_parse_arguments(EXACA_UNIT_TEST "" "" "NAMES" ${ARGN})
  set(EXACA_UNIT_TEST_MPIEXEC_NUMPROCS 1)
  foreach( _np 2 4 )
    if(MPIEXEC_MAX_NUMPROCS GREATER_EQUAL ${_np})
      list(APPEND EXACA_UNIT_TEST_MPIEXEC_NUMPROCS ${_np})
    endif()
  endforeach()
  if(MPIEXEC_MAX_NUMPROCS GREATER 4)
    list(APPEND EXACA_UNIT_TEST_MPIEXEC_NUMPROCS ${MPIEXEC_MAX_NUMPROCS})
  endif()
  set(EXACA_UNIT_TEST_MAIN unit_test_main_nobackend.cpp)

  set(_dir ${CMAKE_CURRENT_SOURCE_DIR})
  foreach(_test ${EXACA_UNIT_TEST_NAMES})
    set(_target ExaCA_${_test}_test)
    add_executable(${_target} tst${_test}.cpp ${EXACA_UNIT_TEST_MAIN})
    target_link_libraries(${_target} ExaCA gtest)

    foreach(_np ${EXACA_UNIT_TEST_MPIEXEC_NUMPROCS})
      add_test(NAME ${_target}_np_${_np} COMMAND
        ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_np} ${MPIEXEC_PREFLAGS}
        ${_target} ${MPIEXEC_POSTFLAGS} ${gtest_args})
    endforeach()
  endforeach()
endmacro()

##--------------------------------------------------------------------------##
## Create main tests
##--------------------------------------------------------------------------##
macro(ExaCA_add_tests)
  cmake_parse_arguments(EXACA_UNIT_TEST "" "" "NAMES" ${ARGN})
  set(EXACA_UNIT_TEST_MPIEXEC_NUMPROCS 1)
  foreach( _np 2 4 )
    if(MPIEXEC_MAX_NUMPROCS GREATER_EQUAL ${_np})
      list(APPEND EXACA_UNIT_TEST_MPIEXEC_NUMPROCS ${_np})
    endif()
  endforeach()
  if(MPIEXEC_MAX_NUMPROCS GREATER 4)
    list(APPEND EXACA_UNIT_TEST_MPIEXEC_NUMPROCS ${MPIEXEC_MAX_NUMPROCS})
  endif()
  set(EXACA_UNIT_TEST_NUMTHREADS 1)
  foreach( _nt 2 4 )
    if(MPIEXEC_MAX_NUMPROCS GREATER_EQUAL ${_nt})
      list(APPEND EXACA_UNIT_TEST_NUMTHREADS ${_nt})
    endif()
  endforeach()
  set(EXACA_UNIT_TEST_MAIN unit_test_main.cpp)
  foreach(_device ${EXACA_TEST_DEVICES})
    set(_dir ${CMAKE_CURRENT_BINARY_DIR}/${_device})
    file(MAKE_DIRECTORY ${_dir})
    foreach(_test ${EXACA_UNIT_TEST_NAMES})
      set(_file ${_dir}/tst${_test}_${_device}.cpp)
      file(WRITE ${_file}
        "#include <Test${_device}_Category.hpp>\n"
        "#include <tst${_test}.hpp>\n"
      )
      set(_target ExaCA_${_test}_test_${_device})
      add_executable(${_target} ${_file} ${EXACA_UNIT_TEST_MAIN})
      target_include_directories(${_target} PRIVATE ${_dir}
        ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
      target_link_libraries(${_target} ExaCA gtest)

      foreach(_np ${EXACA_UNIT_TEST_MPIEXEC_NUMPROCS})
        if(_device STREQUAL PTHREAD OR _device STREQUAL OPENMP)
          foreach(_thread ${EXACA_UNIT_TEST_NUMTHREADS})
            add_test(NAME ${_target}_np_${_np}_nt_${_thread} COMMAND
              ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_np} ${MPIEXEC_PREFLAGS}
              ${_target} ${MPIEXEC_POSTFLAGS} ${gtest_args} --kokkos-threads=${_thread})
          endforeach()
        else()
	    add_test(NAME ${_target}_np_${_np} COMMAND
              ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_np} ${MPIEXEC_PREFLAGS}
              ${_target} ${MPIEXEC_POSTFLAGS} ${gtest_args})
        endif()
      endforeach()
    endforeach()
  endforeach()
endmacro()

ExaCA_add_tests_nobackend(NAMES Init)

ExaCA_add_tests(NAMES SubstrateInitFGS)
