{
  cmake,
  kokkos,
  openmpi,
  nlohmann_json,
  stdenv,
  version,
  src,
}:
stdenv.mkDerivation {

  pname = "exaca";
  inherit version;
  inherit src;

  CMAKE_TLS_VERIFY = 0;

  nativeBuildInputs = [
    cmake
  ];

  buildInputs = [
    kokkos
    openmpi
    nlohmann_json
  ];

  propagatedBuildInputs = [
    openmpi
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
  '';
  
}
