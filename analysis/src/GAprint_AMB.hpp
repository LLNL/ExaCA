// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef GA_PRINT_HPP
#define GA_PRINT_HPP

#include "CAtypes.hpp"
#include <Kokkos_Core.hpp>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

void WriteEulerAngleData_AMB(std::ofstream &Grainplot, ViewI3D_H &GrainID, int i, int j, int k,
                             int NumberOfOrientations, double deltax, ViewF_H &GrainEulerAngles,
                             int &NucleatedGrainCells);
void PrintXYCrossSectionSeries_AMB(std::string BaseFileName, int nx, int ny, int nz, int NumberOfOrientations,
                                   ViewI3D_H &GrainID, double deltax, ViewF_H &GrainEulerAngles, double ZMin);
void PrintDefaultCrossSectionData_AMB(std::string BaseFileName, int nx, int ny, int nz, int NumberOfOrientations,
                                      ViewI3D_H GrainID, double deltax, ViewF_H GrainEulerAngles,
                                      ViewF_H GrainUnitVector, double ZMin);

#endif
