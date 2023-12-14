// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
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

typedef Kokkos::View<float *> ViewF;
using exe_space = Kokkos::DefaultExecutionSpace::execution_space;
typedef typename exe_space::array_layout layout;
typedef Kokkos::View<float *, layout, Kokkos::HostSpace> ViewF_H;
typedef Kokkos::View<short ***, layout, Kokkos::HostSpace> ViewS3D_H;
typedef Kokkos::View<int *, layout, Kokkos::HostSpace> ViewI_H;
typedef Kokkos::View<int ***, layout, Kokkos::HostSpace> ViewI3D_H;

#endif
