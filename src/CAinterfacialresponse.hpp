// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_IRF_HPP
#define EXACA_IRF_HPP

#include "CAparsefiles.hpp"

#include <Kokkos_Core.hpp>

#include <nlohmann/json.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Interfacial repsonse function with various functional forms.
struct InterfacialResponseFunction {

    float freezing_range;
    float A;
    float B;
    float C;
    float D = 0.0;
    enum IRFtypes {
        cubic = 0,
        quadratic = 1,
        power = 2,
    };
    int function = cubic;

    // Constructor
    InterfacialResponseFunction(const int id, const std::string material_file, const double deltat,
                                const double deltax) {
        parseMaterial(id, material_file);
        normalize(deltat, deltax);
    }

    // Used for reading material file in new json format
    void parseMaterial(const int id, const std::string material_file) {
        if (id == 0)
            std::cout << "Parsing material file using json input format" << std::endl;
        std::ifstream material_data(material_file);
        nlohmann::json data = nlohmann::json::parse(material_data);
        A = data["coefficients"]["A"];
        B = data["coefficients"]["B"];
        C = data["coefficients"]["C"];
        std::string functionform = data["function"];
        if (functionform == "cubic") {
            D = data["coefficients"]["D"];
            function = cubic;
        }
        else if ((functionform == "quadratic") || (functionform == "power")) {
            // D should not have been given, this functional form only takes 3 input fitting parameters
            if (data["coefficients"]["D"] != nullptr) {
                std::string error = "Error: functional form of this type takes only A, B, and C as inputs";
                throw std::runtime_error(error);
            }
            if (functionform == "quadratic")
                function = quadratic;
            else if (functionform == "power")
                function = power;
        }
        else
            throw std::runtime_error("Error: Unrecognized functional form for interfacial response function, currently "
                                     "supported options are quadratic, cubic, and exponential");
        freezing_range = data["freezing_range"];
    }

    void normalize(const double deltat, const double deltax) {
        if (function == cubic) {
            // Normalize all 4 coefficients: V = A*x^3 + B*x^2 + C*x + D
            A *= static_cast<float>(deltat / deltax);
            B *= static_cast<float>(deltat / deltax);
            C *= static_cast<float>(deltat / deltax);
            D *= static_cast<float>(deltat / deltax);
        }
        else if (function == quadratic) {
            // Normalize the 3 relevant coefficients: V = A*x^2 + B*x + C
            A *= static_cast<float>(deltat / deltax);
            B *= static_cast<float>(deltat / deltax);
            C *= static_cast<float>(deltat / deltax);
        }
        else if (function == power) {
            // Normalize only the leading and last coefficient: V = A*x^B + C
            A *= static_cast<float>(deltat / deltax);
            C *= static_cast<float>(deltat / deltax);
        }
    }

    // Compute velocity from local undercooling.
    // functional form is assumed to be cubic if not explicitly given in input file
    KOKKOS_INLINE_FUNCTION
    float compute(const float loc_u) const {
        float V;
        if (function == quadratic)
            V = A * Kokkos::pow(loc_u, 2.0) + B * loc_u + C;
        else if (function == power)
            V = A * Kokkos::pow(loc_u, B) + C;
        else
            V = A * Kokkos::pow(loc_u, 3.0) + B * Kokkos::pow(loc_u, 2.0) + C * loc_u + D;
        return Kokkos::fmax(0.0, V);
    }

    std::string functionName() {
        // Not storing string due to Cuda warnings when constructing on device.
        if (function == cubic)
            return "cubic";
        else if (function == quadratic)
            return "quadratic";
        else if (function == power)
            return "power";

        // Should never make it here
        return "none";
    }
    // json format for interfacial response function printing
    std::string print() {
        std::stringstream out;
        out << "   \"InterfacialResponse\": {" << std::endl;
        out << "       \"Function\": "
            << "\"" << functionName() << "\"," << std::endl;
        out << "       \"A\": " << (A) << "," << std::endl;
        out << "       \"B\": " << (B) << "," << std::endl;
        out << "       \"C\": " << (C) << "," << std::endl;
        if (function == cubic)
            out << "       \"D\": " << (D) << "," << std::endl;
        out << "       \"FreezingRange\": " << (freezing_range) << std::endl;
        out << "   },";
        return out.str();
    }
};

#endif
