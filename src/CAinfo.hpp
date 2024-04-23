// Copyright Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_INFO_HPP
#define EXACA_INFO_HPP

#include "CAconfig.hpp"

#include <string>

// Functions for printing for ExaCA/Kokkos version
std::string version();
std::string gitCommitHash();
std::string kokkosVersion();

#endif
