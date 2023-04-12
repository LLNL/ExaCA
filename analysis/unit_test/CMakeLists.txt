include(${TEST_HARNESS_DIR}/test_harness.cmake)

ExaCA_add_tests(PACKAGE ExaGA NAMES Utils KokkosUtils)
if(ExaCA_ENABLE_JSON)
   ExaCA_add_tests(PACKAGE ExaGA NAMES RepresentativeRegion)
endif()