set(EXACA_HEADERS
    CAfunctions.hpp
    CAghostnodes.hpp
    CAinitialize.hpp
    CAprint.hpp
    CAupdate.hpp
    runCA.hpp
    runCAremelt.hpp
)
set(EXACA_SOURCES
    CAfunctions.cpp
    CAghostnodes.cpp
    CAinitialize.cpp
    CAprint.cpp
    CAupdate.cpp
    runCA.cpp
    runCAremelt.cpp
)

install(FILES ${EXACA_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(ExaCA STATIC ${EXACA_SOURCES} ${EXACA_HEADERS})
target_link_libraries(ExaCA PUBLIC
    MPI::MPI_CXX
    Kokkos::kokkos
)
target_include_directories(ExaCA PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ExaCA DESTINATION lib)

add_executable(ExaCA-Kokkos main.cpp)
target_link_libraries(ExaCA-Kokkos ExaCA)
install(TARGETS ExaCA-Kokkos DESTINATION bin)
