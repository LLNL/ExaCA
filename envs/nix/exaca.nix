{
  src, version,

  lib, stdenv,

  cmake,

  kokkos, openmpi, nlohmann_json, gtest
}:

stdenv.mkDerivation {
  pname = "exaca";
  inherit version src;

  nativeBuildInputs = [
    cmake
  ];

  buildInputs = [
    kokkos
    nlohmann_json
    gtest
  ];

  propagatedBuildInputs = [
    openmpi
  ];

  cmakeFlags = [
    (lib.cmakeBool "BUILD_SHARED_LIBS" true)
    (lib.cmakeBool "ExaCA_REQUIRE_EXTERNAL_JSON" true)
    (lib.cmakeBool "ExaCA_ENABLE_TESTING" false)
  ];

  doCheck = true;

  checkPhase = ''
    mkdir -p test
    cd test
    cp $src/examples/Inp_SmallDirSolidification.json ./
    substituteInPlace Inp_SmallDirSolidification.json --replace-fail "Inconel625.json" "$src/examples/Materials/Inconel625.json"
    substituteInPlace Inp_SmallDirSolidification.json --replace-fail "GrainOrientationVectors.csv" "$src/examples/Substrate/GrainOrientationVectors.csv"
    ../bin/ExaCA Inp_SmallDirSolidification.json
    test -e TestProblemSmallDirS.json
    test -e TestProblemSmallDirS_Misorientations.vtk
    test -e TestProblemSmallDirS.vtk
    cd ..
    \rm -rf test
  '';

}
