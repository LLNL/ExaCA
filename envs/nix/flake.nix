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

    lib = with config; {
      callPackage = set: pkgs.lib.callPackageWith (pkgs // set);
    };

    derivations = with config; rec {

      callPackage = lib.callPackage {};

      exaca = {

        devel = callPackage ./exaca.nix {
          src = self;
          version = self.shortRev or self.dirtyShortRev;
        };

        stable = callPackage ./exaca.nix (rec {
          src = pkgs.fetchFromGitHub {
            owner = "LLNL";
            repo  = "ExaCA";
            rev   = "${version}";
            hash  = "sha256-X21yP+sqxR/iM/4N/qKucB2hMmBLf00bVtlCS0QGVQw=";
          };
          version = "2.0.2";
        });

      };
      
    };
    
    packages = with config; rec {
      default = exaca.stable;

      inherit (derivations) exaca;
    };

    devShells = with config; rec {
      default = exacaDev;

      exacaDev = pkgs.mkShell rec {
        name = "exaca-dev";

        packages = with pkgs; [
          git
          clang-tools
        ] ++ pkgs.lib.optionals (pkgs.stdenv.hostPlatform.isLinux) [gdb cntr];

        inputsFrom = [
          self.outputs.packages.${system}.default
        ];

        # Ensure the locales point at the correct archive location.
        LOCALE_ARCHIVE = pkgs.lib.optional (pkgs.stdenv.hostPlatform.isLinux) (
          "${pkgs.glibcLocales}/lib/locale/locale-archive"
        );
        
      };
    };
    
  });

}
