set(EXAGA_HEADERS
    header.hpp
)
set(EXAGA_SOURCES
    runGA.cpp
    GAprint.cpp
    GAutils.cpp
)

add_executable(grain_analysis ${EXAGA_SOURCES})
target_link_libraries(grain_analysis ExaCA)
install(TARGETS grain_analysis DESTINATION ${CMAKE_INSTALL_BINDIR})
