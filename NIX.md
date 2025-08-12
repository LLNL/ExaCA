## Nix

First install the [Nix package manager][NIX] and then enable [Flakes][Flakes].
Alternatively, check out the [Determinate Systems Installer][Determinate] for
an out of the box experience. See [nix.dev][nix.dev] for more help with Nix.

To get a shell with ExaCA temporarily installed, run:

    $ nix shell github:LLNL/ExaCA
    # ExaCA now available
    $ ExaCA --help

To install this permanently, run:

    $ nix profile install github:LLNL/ExaCA

To build from a working copy use `nix develop` and run CMake manually:

    $ nix develop
    $ cmake -B build
    $ cmake --build build

[NIX]: https://nixos.org/download.html
[Flakes]: https://nixos.wiki/wiki/Flakes
[nix.dev]: https://nix.dev
[Determinate]: https://github.com/DeterminateSystems/nix-installer
