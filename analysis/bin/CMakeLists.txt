add_executable(grain_analysis runGA.cpp)
target_link_libraries(grain_analysis ExaGA)
install(TARGETS grain_analysis DESTINATION ${CMAKE_INSTALL_BINDIR})
