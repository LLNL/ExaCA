set(EXAGA_HEADERS
    header.h
)
set(EXAGA_SOURCES
    DataAnalysis_PrintRoutines.cpp
    FileReadFunctions.cpp
)

install(FILES ${EXAGA_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/analysis)

add_library(ExaGA STATIC ${EXAGA_SOURCES} ${EXAGA_HEADERS})
target_link_libraries(ExaGA PUBLIC ExaCA)
target_include_directories(ExaGA PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS ExaGA
  EXPORT ExaCA_Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/analysis
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/analysis)

add_library(ExaCA::ExaGA ALIAS ExaGA)

add_executable(runExaGA DataAnalysis_main.cpp)
target_link_libraries(runExaGA ExaGA)
install(TARGETS runExaGA DESTINATION ${CMAKE_INSTALL_BINDIR})
