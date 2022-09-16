// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
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
    FutureActive = 7
};

// Use Kokkos::DefaultExecutionSpace
typedef Kokkos::View<float *> ViewF;
typedef Kokkos::View<int *> ViewI;
typedef Kokkos::View<int **> ViewI2D;
typedef Kokkos::View<int *, Kokkos::MemoryTraits<Kokkos::Atomic>> View_a;
typedef Kokkos::View<float ***> ViewF3D;

using exe_space = Kokkos::DefaultExecutionSpace::execution_space;
using device_memory_space = Kokkos::DefaultExecutionSpace::memory_space;
typedef typename exe_space::array_layout layout;
typedef Kokkos::View<float *, layout, Kokkos::HostSpace> ViewF_H;
typedef Kokkos::View<float **, layout, Kokkos::HostSpace> ViewF2D_H;
typedef Kokkos::View<float ***, layout, Kokkos::HostSpace> ViewF3D_H;
typedef Kokkos::View<int *, layout, Kokkos::HostSpace> ViewI_H;
typedef Kokkos::View<int **, layout, Kokkos::HostSpace> ViewI2D_H;
typedef Kokkos::View<int ***, layout, Kokkos::HostSpace> ViewI3D_H;

typedef Kokkos::Array<int, 26> NList;

struct CommData {
    int GrainID;
    float DiagonalLength;
    float DOCenterX;
    float DOCenterY;
    float DOCenterZ;
};
typedef Kokkos::View<CommData *> ViewBuffer;
#endif
