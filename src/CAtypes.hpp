// Copyright Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_TYPES_HPP
#define EXACA_TYPES_HPP

#include <Kokkos_Core.hpp>

enum TypeNames {
    Wall = 0,
    Solid = 1,
    Active = 2,
    TemporaryUpdate = 3,
    TemporaryInit = 4,
    Liquid = 5,
    TempSolid = 6,
    FutureActive = 7,
    ActiveFailedBufferLoad = 8,
    FutureLiquid = 9,
    LiquidFailedBufferLoad = 10
};

#endif
