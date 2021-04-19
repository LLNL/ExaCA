#ifndef EXACA_TYPES_HPP
#define EXACA_TYPES_HPP

#include <Kokkos_Core.hpp>

enum TypeNames { Wall = 0, Solid = 1, Active = 2, TemporaryUpdate = 3, TemporaryInit = 4, Liquid = 5, Ghost = 6 };

// Use Kokkos::DefaultExecutionSpace
typedef Kokkos::View<float *> ViewF;
typedef Kokkos::View<int *> ViewI;
typedef Kokkos::View<int **> ViewI2D;
typedef Kokkos::View<int *, Kokkos::MemoryTraits<Kokkos::Atomic>> View_a;
typedef Kokkos::View<float **> Buffer2D;
typedef Kokkos::View<int ***> Buffer3D; // Used in ghost node initialization of integer structures CellType and GrainID
typedef Kokkos::View<float *> TestView;

using exe_space = Kokkos::DefaultExecutionSpace::execution_space;
typedef typename exe_space::array_layout layout;
typedef Kokkos::View<float *, layout, Kokkos::HostSpace> ViewF_H;
typedef Kokkos::View<int *, layout, Kokkos::HostSpace> ViewI_H;
typedef Kokkos::View<int **, layout, Kokkos::HostSpace> ViewI2D_H;

#endif
