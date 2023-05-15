configure_file(CAconfig.hpp.cmakein CAconfig.hpp)

set(EXACA_HEADERS
    CAfunctions.hpp
    CAghostnodes.hpp
    CAinitialize.hpp
    CAinterfacialresponse.hpp
    CAparsefiles.hpp
    CAprint.hpp
    CAtypes.hpp
    CAupdate.hpp
    ExaCA.hpp
    runCA.hpp
)
set(EXACA_SOURCES
    CAfunctions.cpp
    CAghostnodes.cpp
    CAinitialize.cpp
    CAparsefiles.cpp
    CAprint.cpp
    CAupdate.cpp
    runCA.cpp
)

install(FILES ${EXACA_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CAconfig.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(ExaCA STATIC ${EXACA_SOURCES} ${EXACA_HEADERS})
target_link_libraries(ExaCA PUBLIC
    MPI::MPI_CXX
    Kokkos::kokkos
    nlohmann_json::nlohmann_json
)

target_include_directories(ExaCA PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS ExaCA
  EXPORT ExaCA_Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_library(ExaCA::ExaCA ALIAS ExaCA)

add_executable(ExaCA-Kokkos main.cpp)
target_link_libraries(ExaCA-Kokkos ExaCA)
install(TARGETS ExaCA-Kokkos DESTINATION ${CMAKE_INSTALL_BINDIR})
