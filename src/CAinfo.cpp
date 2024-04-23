// Copyright Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAinfo.hpp"

#include <string>

// Functions for printing for ExaCA/Kokkos version
std::string version() { return ExaCA_VERSION; }

std::string gitCommitHash() { return ExaCA_GIT_COMMIT_HASH; }

std::string kokkosVersion() { return ExaCA_Kokkos_VERSION_STRING; }
