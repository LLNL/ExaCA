include(${TEST_HARNESS_DIR}/test_harness.cmake)

exaca_add_tests(PACKAGE ExaCA-Analysis NAMES Utils KokkosUtils)
exaca_add_tests(PACKAGE ExaCA-Analysis NAMES RepresentativeRegion)
