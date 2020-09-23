set(EXACA_SOURCES
    CAfunctions.cpp
    CAghostnodes.cpp
    CAInitialize.cpp
    CAPrint.cpp
    CAupdate.cpp
    Run_Reduced.cpp
)

set(EXACA_HEADERS header.h)

add_library(ExaCA STATIC ${EXACA_SOURCES} ${EXACA_HEADERS})
target_link_libraries(ExaCA PUBLIC
    MPI::MPI_CXX
    Kokkos::kokkos
)
install(TARGETS ExaCA DESTINATION lib)

add_executable(ExaCA-Kokkos main.cpp)
target_link_libraries(ExaCA-Kokkos ExaCA)
install(TARGETS ExaCA-Kokkos DESTINATION bin)
