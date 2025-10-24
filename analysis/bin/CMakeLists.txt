add_executable(ExaCA-GrainAnalysis runGA.cpp)
target_link_libraries(ExaCA-GrainAnalysis ExaCA-Analysis)
install(TARGETS ExaCA-GrainAnalysis DESTINATION ${CMAKE_INSTALL_BINDIR})
