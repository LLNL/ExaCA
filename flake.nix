## See NIX.md for help getting started with Nix

{
  description = "An exascale-capable cellular automaton for nucleation and grain growth";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-25.05";
    utils.url   = "github:numtide/flake-utils";
  };

  outputs = inputs @ { self, utils, ... }: utils.lib.eachDefaultSystem (system: rec {
    config = rec {
      pkgs = import inputs.nixpkgs {
        inherit system;        
        config.cudaSupport = true;
        config.allowUnfree = true;
      };
    };

    exaca = with config; pkgs.stdenv.mkDerivation rec {
      pname = "exaca";
      version = "latest";

      src = self;

      CMAKE_TLS_VERIFY=0;
      
      nativeBuildInputs = [
        pkgs.cmake
      ];

      buildInputs = [
        pkgs.kokkos
        pkgs.openmpi
      ];

      propagatedBuildInputs = [
        pkgs.openmpi
      ];

    };
      
    packages = rec {
      default = exaca;
    };

    devShells = with config; rec {
      default = exacaDev;

      exacaDev = pkgs.mkShell rec {
        name = "exaca-dev";

        packages = with pkgs; [
          git
          clang-tools
        ] ++ self.outputs.packages.${system}.default.buildInputs
          ++ self.outputs.packages.${system}.default.nativeBuildInputs
          ++ self.outputs.packages.${system}.default.propagatedBuildInputs
          ++ pkgs.lib.optionals (pkgs.stdenv.hostPlatform.isLinux) [gdb cntr];

        # Ensure the locales point at the correct archive location.
        LOCALE_ARCHIVE = pkgs.lib.optional (pkgs.stdenv.hostPlatform.isLinux) (
          "${pkgs.glibcLocales}/lib/locale/locale-archive"
        );
        
      };
    };
    
  });

}
