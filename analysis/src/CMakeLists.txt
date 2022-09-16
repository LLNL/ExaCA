set(EXAGA_HEADERS
    runGA.hpp
    GAprint.hpp
    GAutils.hpp
)
set(EXAGA_SOURCES
    runGA.cpp
    GAprint.cpp
    GAutils.cpp
)
set(EXAGA_AMB_HEADERS
    runGA_AMB.hpp
    GAprint_AMB.hpp
    GAutils.hpp
)
set(EXAGA_AMB_SOURCES
    runGA_AMB.cpp
    GAprint_AMB.cpp
    GAutils.cpp
)

add_executable(grain_analysis ${EXAGA_SOURCES})
target_link_libraries(grain_analysis ExaCA)
install(TARGETS grain_analysis DESTINATION ${CMAKE_INSTALL_BINDIR})
add_executable(amb_grain_analysis ${EXAGA_AMB_SOURCES})
target_link_libraries(amb_grain_analysis ExaCA)
install(TARGETS amb_grain_analysis DESTINATION ${CMAKE_INSTALL_BINDIR})
