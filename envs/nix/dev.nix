{ self', pkgs, lib, ... }:

{
  devShells = rec {
    default = exacaDev;

    exacaDev = pkgs.mkShell {
      name = "exaca-dev";

      packages = with pkgs; [
        git
        clang-tools
      ] ++ lib.optionals (pkgs.stdenv.hostPlatform.isLinux) [
        gdb
      ];

      inputsFrom = [
        self'.packages.default
      ];

      # Ensure the locales point at the correct archive location.
      LOCALE_ARCHIVE = lib.optional (pkgs.stdenv.hostPlatform.isLinux) (
        "${pkgs.glibcLocales}/lib/locale/locale-archive"
      );
    };
  };
}
