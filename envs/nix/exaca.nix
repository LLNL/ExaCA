{
  cmake,
  kokkos,
  openmpi,
  stdenv,
  version,
  src
}:
stdenv.mkDerivation {

  pname = "exaca";
  inherit version;
  inherit src;

  CMAKE_TLS_VERIFY=0;
      
  nativeBuildInputs = [
    cmake
  ];

  buildInputs = [
    kokkos
    openmpi
  ];

  propagatedBuildInputs = [
    openmpi
  ];
  
}
