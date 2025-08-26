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

}
