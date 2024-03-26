configure_file(CAconfig.hpp.cmakein CAconfig.hpp)

file(GLOB EXACA_HEADERS GLOB *.hpp)
file(GLOB EXACA_SOURCES GLOB *.cpp)

install(FILES ${EXACA_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CAconfig.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(ExaCA-Core ${EXACA_SOURCES})
add_library(ExaCA::Core ALIAS ExaCA-Core)

target_link_libraries(ExaCA-Core PUBLIC
    MPI::MPI_CXX
    Kokkos::kokkos
    nlohmann_json::nlohmann_json
)
if(ExaCA_ENABLE_FINCH)
  target_link_libraries(ExaCA-Core PUBLIC Finch::Core Finch::ScanPaths)
endif()

target_include_directories(ExaCA-Core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS ExaCA-Core
  EXPORT ExaCA_Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(EXPORT ExaCA_Targets
  FILE ExaCA_Targets.cmake
  NAMESPACE ExaCA::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ExaCA)
