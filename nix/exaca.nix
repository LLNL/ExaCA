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
}
