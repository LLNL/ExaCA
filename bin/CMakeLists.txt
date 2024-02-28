add_executable(ExaCA run.cpp)
target_link_libraries(ExaCA ExaCA-Core)
install(TARGETS ExaCA DESTINATION ${CMAKE_INSTALL_BINDIR})
