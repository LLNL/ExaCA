add_subdirectory(src)

if(ExaCA_ENABLE_TESTING)
  add_subdirectory(unit_test)
endif()

add_subdirectory(bin)
