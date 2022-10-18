// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_IRF_HPP
#define EXACA_IRF_HPP

#include "CAparsefiles.hpp"
#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Using for compatibility with device math functions.
using std::max;

// Cubic interfacial repsonse function.
struct InterfacialResponseFunction {

    double FreezingRange;
    double A;
    double B;
    double C;
    double D;

    // Constructor. Read data from file.
    InterfacialResponseFunction(std::string MaterialFile) {
        std::ifstream MaterialData;
        MaterialData.open(MaterialFile);
        skipLines(MaterialData, "*****");
        std::string val;
        // Interfacial response function A, B, C, D, and the solidification range
        // for the alloy The order of these is important: "Alloy freezing range"
        // should be before "A", as a search for "A" in either string will return
        // true
        std::vector<std::string> MaterialInputs = {
            "Alloy freezing range", // Required input 0
            "A",                    // Required input 1
            "B",                    // Required input 2
            "C",                    // Required input 3
            "D",                    // Required input 4
        };
        int NumMaterialInputs = MaterialInputs.size();
        std::vector<std::string> MaterialInputsRead(NumMaterialInputs);
        while (std::getline(MaterialData, val)) {
            // Check if this is one of the expected inputs - otherwise throw an error
            bool FoundInput = parseInputFromList(val, MaterialInputs, MaterialInputsRead, NumMaterialInputs);
            if (!(FoundInput)) {
                std::string error = "Error: Unexpected line " + val + " present in material file " + MaterialFile +
                                    " : file should only contain A, B, C, D, and Alloy freezing range";
                throw std::runtime_error(error);
            }
        }

        FreezingRange = getInputDouble(MaterialInputsRead[0]);
        A = getInputDouble(MaterialInputsRead[1]);
        B = getInputDouble(MaterialInputsRead[2]);
        C = getInputDouble(MaterialInputsRead[3]);
        D = getInputDouble(MaterialInputsRead[4]);

        MaterialData.close();
    }

    // Compute velocity from local undercooling.
    KOKKOS_INLINE_FUNCTION
    double compute(const double LocU) const {
        double V = A * pow(LocU, 3.0) + B * pow(LocU, 2.0) + C * LocU + D;
        return max(0.0, V);
    }

    void normalize(const double deltat, const double deltax) {
        A *= deltat / deltax;
        B *= deltat / deltax;
        C *= deltat / deltax;
        D *= deltat / deltax;
    }

    std::string print() {
        std::stringstream out;
        out << "Interfacial response cubic function parameters used were " << (A) << ", " << (B) << ", " << (C) << ", "
            << (D) << ", and the alloy freezing range was " << (FreezingRange);
        return out.str();
    }
};

#endif
