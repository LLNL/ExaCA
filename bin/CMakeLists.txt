add_executable(ExaCA run.cpp)
target_link_libraries(ExaCA ExaCA-Core)
install(TARGETS ExaCA DESTINATION ${CMAKE_INSTALL_BINDIR})

if(ExaCA_ENABLE_FINCH)
  add_executable(Finch-ExaCA runCoupled.cpp)
  target_link_libraries(Finch-ExaCA ExaCA-Core)
  install(TARGETS Finch-ExaCA DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
