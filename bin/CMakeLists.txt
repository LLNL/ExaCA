add_executable(ExaCA main.cpp)
target_link_libraries(ExaCA ExaCA-Core)
install(TARGETS ExaCA DESTINATION ${CMAKE_INSTALL_BINDIR})

# Symlink for backwards compatibility
add_custom_target(ExaCA-Kokkos ALL
                  COMMAND ${CMAKE_COMMAND} -E create_symlink ExaCA ExaCA-Kokkos)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ExaCA-Kokkos
        DESTINATION ${CMAKE_INSTALL_BINDIR})
