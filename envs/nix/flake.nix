## See NIX.md for help getting started with Nix

{
  description = "An exascale-capable cellular automaton for nucleation and grain growth";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/25.05";
    parts.url = "github:hercules-ci/flake-parts";
  };

  outputs = inputs @ { self, parts, ... }: (
    parts.lib.mkFlake { inherit inputs; } {
      systems = [
        "x86_64-linux"
        "aarch64-linux"
      ];

      perSystem = { pkgs, ... }: {
        packages = rec {
          default = exaca;

          exaca = pkgs.callPackage ./exaca.nix {
            src = ../..;
            version = "master";
          };

        };

        imports = [
          ./dev.nix
        ];
      };
    }
  );
}
