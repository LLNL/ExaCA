// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_IRF_HPP
#define EXACA_IRF_HPP

#include "CAinputs.hpp"

#include <Kokkos_Core.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Interfacial response function with various functional forms.
struct InterfacialResponseFunction {

    InterfacialResponseInputs _inputs;

    // Constructor
    InterfacialResponseFunction(const double deltat, const double deltax, const InterfacialResponseInputs inputs)
        : _inputs(inputs) {
        normalize(deltat, deltax);
    }

    void normalize(const double deltat, const double deltax) {
        if (_inputs.function == _inputs.cubic) {
            // Normalize all 4 coefficients: V = A*x^3 + B*x^2 + C*x + D
            _inputs.A *= static_cast<float>(deltat / deltax);
            _inputs.B *= static_cast<float>(deltat / deltax);
            _inputs.C *= static_cast<float>(deltat / deltax);
            _inputs.D *= static_cast<float>(deltat / deltax);
        }
        else if (_inputs.function == _inputs.quadratic) {
            // Normalize the 3 relevant coefficients: V = A*x^2 + B*x + C
            _inputs.A *= static_cast<float>(deltat / deltax);
            _inputs.B *= static_cast<float>(deltat / deltax);
            _inputs.C *= static_cast<float>(deltat / deltax);
        }
        else if (_inputs.function == _inputs.power) {
            // Normalize only the leading and last coefficient: V = A*x^B + C
            _inputs.A *= static_cast<float>(deltat / deltax);
            _inputs.C *= static_cast<float>(deltat / deltax);
        }
    }

    // Compute velocity from local undercooling.
    // functional form is assumed to be cubic if not explicitly given in input file
    KOKKOS_INLINE_FUNCTION
    float compute(const float loc_u) const {
        float V;
        if (_inputs.function == _inputs.quadratic)
            V = _inputs.A * Kokkos::pow(loc_u, 2.0) + _inputs.B * loc_u + _inputs.C;
        else if (_inputs.function == _inputs.power)
            V = _inputs.A * Kokkos::pow(loc_u, _inputs.B) + _inputs.C;
        else
            V = _inputs.A * Kokkos::pow(loc_u, 3.0) + _inputs.B * Kokkos::pow(loc_u, 2.0) + _inputs.C * loc_u +
                _inputs.D;
        return Kokkos::fmax(0.0, V);
    }

    std::string functionName() {
        // Not storing string due to Cuda warnings when constructing on device.
        if (_inputs.function == _inputs.cubic)
            return "cubic";
        else if (_inputs.function == _inputs.quadratic)
            return "quadratic";
        else if (_inputs.function == _inputs.power)
            return "power";

        // Should never make it here
        return "none";
    }

    auto A() { return _inputs.A; }
    auto B() { return _inputs.B; }
    auto C() { return _inputs.C; }
    auto D() { return _inputs.D; }
    auto freezingRange() { return _inputs.freezing_range; }
};

#endif
