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
    GAprint.hpp
    GAutils.hpp
)
set(EXAGA_AMB_SOURCES
    runGA_AMB.cpp
    GAprint.cpp
    GAutils.cpp
)

add_executable(grain_analysis ${EXAGA_SOURCES})
target_link_libraries(grain_analysis ExaCA)
install(TARGETS grain_analysis DESTINATION ${CMAKE_INSTALL_BINDIR})
add_executable(grain_analysis_amb ${EXAGA_AMB_SOURCES})
target_link_libraries(grain_analysis_amb ExaCA)
install(TARGETS grain_analysis_amb DESTINATION ${CMAKE_INSTALL_BINDIR})
