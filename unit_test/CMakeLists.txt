include(${TEST_HARNESS_DIR}/test_harness.cmake)

# --------------------------------------------------------------------------##
# Install example files used in tests
# --------------------------------------------------------------------------##
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/../examples/Inp_SmallDirSolidification.json
  ${CMAKE_CURRENT_BINARY_DIR}/Inp_SmallDirSolidification.json COPYONLY)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/../examples/Inp_DirSolidification.json
  ${CMAKE_CURRENT_BINARY_DIR}/Inp_DirSolidification.json COPYONLY)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/../examples/Inp_TwoGrainDirSolidification.json
  ${CMAKE_CURRENT_BINARY_DIR}/Inp_TwoGrainDirSolidification.json COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../examples/Inp_SpotMelt.json
               ${CMAKE_CURRENT_BINARY_DIR}/Inp_SpotMelt.json COPYONLY)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/../examples/Inp_SpotMelt.json
  ${CMAKE_CURRENT_BINARY_DIR}/Inp_SmallSpotMelt_Steel2phase.json COPYONLY)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/../examples/Inp_SmallEquiaxedGrain.json
  ${CMAKE_CURRENT_BINARY_DIR}/Inp_SmallEquiaxedGrain.json COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../examples/Materials/Inconel625.json
               ${CMAKE_CURRENT_BINARY_DIR}/Inconel625.json COPYONLY)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/../examples/Materials/Inconel625_Quadratic.json
  ${CMAKE_CURRENT_BINARY_DIR}/Inconel625_Quadratic.json COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../examples/Materials/SS316.json
               ${CMAKE_CURRENT_BINARY_DIR}/SS316.json COPYONLY)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/../examples/Materials/SS316L_AusteniteFerrite.json
  ${CMAKE_CURRENT_BINARY_DIR}/SS316L_AusteniteFerrite.json COPYONLY)

exaca_add_tests(
  PACKAGE
  ExaCA-Core
  NAMES
  InterfacialResponse
  Inputs
  Orientation
  Parse
  Print)

exaca_add_tests(
  MPI
  PACKAGE
  ExaCA-Core
  NAMES
  CellData
  Grid
  Interface
  Nucleation
  Temperature
  Update)
