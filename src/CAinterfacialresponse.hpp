// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_IRF_HPP
#define EXACA_IRF_HPP

#include "CAconfig.hpp"
#include "CAparsefiles.hpp"
#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifdef ExaCA_ENABLE_JSON
#include <nlohmann/json.hpp>
#endif

// Using for compatibility with device math functions.
using std::max;

// Type tags for different functional forms.
struct Quadratic {};
struct Cubic {};
struct Exponential {};

// "quadratic", "cubic", or "exponential" interfacial repsonse function
// Quadratic function has the form V = A*(Undercooling)^2 + B*(Undercooling) + C
// Cubic function has the form V = A*(Undercooling)^3 + B*(Undercooling)^2 + C*(Undercooling) + D
// Exponential function has the form of V = A * (Undercooling)^B + C
struct InterfacialResponseFunctionBase {

    double FreezingRange;
    double A;
    double B;
    double C;
    double D = 0;

    // Constructor for 4 parameter forms.
    InterfacialResponseFunctionBase(const double _A, const double _B, const double _C, const double _D, const double FR,
                                    const double deltat, const double deltax)
        : FreezingRange(FR)
        , A(_A)
        , B(_B)
        , C(_C)
        , D(_D) {
        normalize(deltat, deltax);
    }
    // Constructor for 3 parameter forms.
    InterfacialResponseFunctionBase(const double _A, const double _B, const double _C, const double FR,
                                    const double deltat, const double deltax)
        : FreezingRange(FR)
        , A(_A)
        , B(_B)
        , C(_C) {
        normalize(deltat, deltax);
    }

    virtual ~InterfacialResponseFunctionBase() = default;

    KOKKOS_INLINE_FUNCTION virtual double compute(const double) const = 0;

    double scale(const double V) { return max(0.0, V); }

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

template <typename FunctionType>
struct InterfacialResponseFunction;

template <>
struct InterfacialResponseFunction<Cubic> : InterfacialResponseFunctionBase {
    using base_type = InterfacialResponseFunctionBase;

    using base_type::A;
    using base_type::B;
    using base_type::C;
    using base_type::D;
    using base_type::FreezingRange;

    using base_type::base_type;

    // Compute velocity from local undercooling.
    KOKKOS_INLINE_FUNCTION double compute(const double LocU) const {
        return A * pow(LocU, 3.0) + B * pow(LocU, 2.0) + C * LocU + D;
    }
};

template <>
struct InterfacialResponseFunction<Quadratic> : InterfacialResponseFunctionBase {
    using base_type = InterfacialResponseFunctionBase;

    using base_type::A;
    using base_type::B;
    using base_type::C;
    using base_type::D;
    using base_type::FreezingRange;

    using base_type::base_type;

    // Compute velocity from local undercooling.
    KOKKOS_INLINE_FUNCTION double compute(const double LocU) const { return A * pow(LocU, 2.0) + B * LocU + C; }
};

template <>
struct InterfacialResponseFunction<Exponential> : InterfacialResponseFunctionBase {
    using base_type = InterfacialResponseFunctionBase;

    using base_type::A;
    using base_type::B;
    using base_type::C;
    using base_type::D;
    using base_type::FreezingRange;

    using base_type::base_type;

    // Compute velocity from local undercooling.
    KOKKOS_INLINE_FUNCTION double compute(const double LocU) const { return A * pow(LocU, B) + C; }
};

#ifdef ExaCA_ENABLE_JSON
// Used for reading material file in new json format
inline std::shared_ptr<InterfacialResponseFunctionBase> parseMaterialJson(std::string MaterialFile, const double deltat,
                                                                          const double deltax) {
    std::ifstream MaterialData(MaterialFile);
    nlohmann::json data = nlohmann::json::parse(MaterialData);
    double FR = data["freezingRange"];
    double A = data["fittingParams"]["A"];
    double B = data["fittingParams"]["B"];
    double C = data["fittingParams"]["C"];
    double D = -1;
    if (data["function"] == "cubic") {
        D = data["fittingParams"]["D"];
        return std::make_shared<InterfacialResponseFunction<Cubic>>(A, B, C, D, FR, deltat, deltax);
    }
    else if ((data["function"] == "quadratic") || (data["function"] == "exponential")) {
        // D should not have been given, this functional form only takes 3 input fitting parameters
        if (data["fittingParams"]["D"] != nullptr) {
            std::string error = "Error: functional form of this type takes only A, B, and C as inputs";
            throw std::runtime_error(error);
        }
        if (data["function"] == "quadratic")
            return std::make_shared<InterfacialResponseFunction<Quadratic>>(A, B, C, FR, deltat, deltax);
        else if (data["function"] == "exponential")
            return std::make_shared<InterfacialResponseFunction<Exponential>>(A, B, C, FR, deltat, deltax);
    }

    throw std::runtime_error("Error: Unrecognized functional form for interfacial response function, currently "
                             "supported options are quadratic, cubic, and exponential");
}
#endif

inline std::shared_ptr<InterfacialResponseFunctionBase> createIRF(int id, std::string MaterialFile, const double deltat,
                                                                  const double deltax) {
    std::ifstream MaterialData;
    MaterialData.open(MaterialFile);
    std::string firstline;
    getline(MaterialData, firstline);
    if (firstline.find('{') != std::string::npos) {
#ifdef ExaCA_ENABLE_JSON
        // New json input format for material data detected
        MaterialData.close();
        return parseMaterialJson(MaterialFile, deltat, deltax);
#else
        throw std::runtime_error("Cannot use JSON input file without ExaCA_ENABLE_JSON=ON");
#endif
    }
    else {
        // Old input file format for material data (cubic function)
        if (id == 0)
            std::cout << "Warning: Old (non-json) input file format detected for the material data file, compatability "
                         "with this is now deprecated and will be removed in a future release. Note this format only "
                         "allows cubic functions."
                      << std::endl;
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

        double FreezingRange = getInputDouble(MaterialInputsRead[0]);
        double A = getInputDouble(MaterialInputsRead[1]);
        double B = getInputDouble(MaterialInputsRead[2]);
        double C = getInputDouble(MaterialInputsRead[3]);
        double D = getInputDouble(MaterialInputsRead[4]);

        MaterialData.close();

        return std::make_shared<InterfacialResponseFunction<Cubic>>(A, B, C, D, FreezingRange, deltat, deltax);
    }
}

#endif
