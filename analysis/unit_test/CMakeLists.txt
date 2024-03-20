include(${TEST_HARNESS_DIR}/test_harness.cmake)

ExaCA_add_tests(PACKAGE ExaCA-Analysis NAMES Utils KokkosUtils)
ExaCA_add_tests(PACKAGE ExaCA-Analysis NAMES RepresentativeRegion)
