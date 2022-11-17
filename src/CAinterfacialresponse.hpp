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

#include <nlohmann/json.hpp>

// Using for compatibility with device math functions.
using std::max;

// "quadratic", "cubic", or "exponential" interfacial repsonse function
// Quadratic function has the form V = A*(Undercooling)^2 + B*(Undercooling) + C
// Cubic function has the form V = A*(Undercooling)^3 + B*(Undercooling)^2 + C*(Undercooling) + D
// Expontential function has the form of V = A * (Undercooling)^B + C
struct InterfacialResponseFunction {

    double FreezingRange;
    double A;
    double B;
    double C;
    double D;
    int function;
    enum IRFforms { Cubic = 0, Quadratic = 1, Expontential = 2 };

    // Constructor. Read data from file.
    InterfacialResponseFunction(int id, std::string MaterialFile, const double deltat, const double deltax) {
        std::ifstream MaterialData;
        MaterialData.open(MaterialFile);
        std::string firstline;
        getline(MaterialData, firstline);
        if (firstline.find('{') != std::string::npos) {
            // New json input format for material data detected
            MaterialData.close();
            parseMaterialJson(MaterialFile);
        }
        else {
            // Old input file format for material data (cubic function)
            if (id == 0)
                std::cout << "Warning: Old (non-json) input file format detected for the material data file, "
                             "compatability with this is now deprecated and will be removed in a future release"
                          << std::endl;
            function = 0;
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

        normalize(deltat, deltax);
    }

    // Used for reading material file in new json format
    void parseMaterialJson(std::string MaterialFile) {
        std::ifstream MaterialData(MaterialFile);
        nlohmann::json data = nlohmann::json::parse(MaterialData);
        A = data["fittingParams"]["A"];
        B = data["fittingParams"]["B"];
        C = data["fittingParams"]["C"];
        if (data["function"] == "cubic") {
            D = data["fittingParams"]["D"];
            function = 0;
        }
        else if ((data["function"] == "quadratic") || (data["function"] == "exponential")) {
            // D should not have been given, this functional form only takes 3 input fitting parameters
            if (data["fittingParams"]["D"] != nullptr) {
                std::string error = "Error: functional form of this type takes only A, B, and C as inputs";
                throw std::runtime_error(error);
            }
            if (data["function"] == "quadratic")
                function = 1;
            else if (data["function"] == "exponential")
                function = 2;
        }
        else
            throw std::runtime_error("Error: Unrecognized functional form for interfacial response function, currently "
                                     "supported options are quadratic, cubic, and exponential");
        FreezingRange = data["freezingRange"];
    }

    // Compute velocity from local undercooling.
    KOKKOS_INLINE_FUNCTION
    double compute(const double LocU) const {
        double V;
        switch (function) {
        case 0:
            V = A * pow(LocU, 3.0) + B * pow(LocU, 2.0) + C * LocU + D;
            break;
        case 1:
            V = A * pow(LocU, 2.0) + B * LocU + C;
            break;
        case 2:
            V = A * pow(LocU, B) + C;
            break;
        }
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
