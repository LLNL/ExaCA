// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_INIT_HPP
#define EXACA_INIT_HPP

#include "CAinputs.hpp"
#include "CAprint.hpp"
#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

#include <nlohmann/json.hpp>

#include <iostream>
#include <string>
#include <vector>

//*****************************************************************************/
// TODO: Turn into an orientation struct used by ExaCA and the analysis routine
// Initialize grain orientations and unit vectors
template <typename ViewTypeFloat>
void OrientationInit(int, int &NGrainOrientations, ViewTypeFloat &GrainOrientationData,
                     std::string GrainOrientationFile, int ValsPerLine = 9) {

    // Read file of grain orientations
    std::ifstream O;
    O.open(GrainOrientationFile);

    // Line 1 is the number of orientation values to read (if not specified already)
    std::string ValueRead;
    getline(O, ValueRead);
    NGrainOrientations = getInputInt(ValueRead);

    // Temporary host view for storing grain orientations read from file
    using view_type_host = typename ViewTypeFloat::HostMirror;
    view_type_host GrainOrientationData_Host(Kokkos::ViewAllocateWithoutInitializing("GrainOrientationData_H"),
                                             ValsPerLine * NGrainOrientations);
    // Populate data structure for grain orientation data
    for (int i = 0; i < NGrainOrientations; i++) {
        std::vector<std::string> ParsedLine(ValsPerLine);
        std::string ReadLine;
        if (!getline(O, ReadLine))
            break;
        splitString(ReadLine, ParsedLine, ValsPerLine);
        // Place the 3 grain orientation angles or 9 rotation matrix components into the orientation data view
        for (int Comp = 0; Comp < ValsPerLine; Comp++) {
            GrainOrientationData_Host(ValsPerLine * i + Comp) = getInputFloat(ParsedLine[Comp]);
        }
    }
    O.close();

    // Resize device view and orientation data to device
    Kokkos::realloc(GrainOrientationData, ValsPerLine * NGrainOrientations);
    using memory_space = typename ViewTypeFloat::memory_space;
    GrainOrientationData = Kokkos::create_mirror_view_and_copy(memory_space(), GrainOrientationData_Host);
}

#endif
