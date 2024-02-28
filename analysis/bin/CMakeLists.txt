add_executable(ExaCA-GrainAnalysis runGA.cpp)
target_link_libraries(ExaCA-GrainAnalysis ExaCA-Analysis)
install(TARGETS ExaCA-GrainAnalysis DESTINATION ${CMAKE_INSTALL_BINDIR})

# Symlink for backwards compatibility
add_custom_target(grain_analysis ALL
                  COMMAND ${CMAKE_COMMAND} -E create_symlink ExaCA-GrainAnalysis grain_analysis)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/grain_analysis
        DESTINATION ${CMAKE_INSTALL_BINDIR})
