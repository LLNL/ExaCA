set(EXAGA_HEADERS
    GAprint.hpp
    GAutils.hpp
    GArepresentativeregion.hpp
)
set(EXAGA_SOURCES
    GAprint.cpp
    GAutils.cpp
)

install(FILES ${EXAGA_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(ExaGA STATIC ${EXAGA_SOURCES} ${EXAGA_HEADERS})
target_link_libraries(ExaGA PUBLIC ExaCA)

target_include_directories(ExaGA PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS ExaGA
  EXPORT ExaGA_Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_library(ExaCA::ExaGA ALIAS ExaGA)

add_executable(grain_analysis runGA.cpp)
target_link_libraries(grain_analysis ExaGA)
install(TARGETS grain_analysis DESTINATION ${CMAKE_INSTALL_BINDIR})