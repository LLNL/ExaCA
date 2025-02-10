// Copyright Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_IRF_HPP
#define EXACA_IRF_HPP

#include "CAinputdata.hpp"

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
        for (int phase_num = 0; phase_num < _inputs.num_phases; phase_num++)
            normalize(deltat, deltax, phase_num);
    }

    void normalize(const double deltat, const double deltax, const int phase_num = 0) {
        if (_inputs.function[phase_num] == _inputs.cubic) {
            // Normalize all 4 coefficients: V = A*x^3 + B*x^2 + C*x + D
            _inputs.A[phase_num] *= static_cast<float>(deltat / deltax);
            _inputs.B[phase_num] *= static_cast<float>(deltat / deltax);
            _inputs.C[phase_num] *= static_cast<float>(deltat / deltax);
            _inputs.D[phase_num] *= static_cast<float>(deltat / deltax);
        }
        else if (_inputs.function[phase_num] == _inputs.quadratic) {
            // Normalize the 3 relevant coefficients: V = A*x^2 + B*x + C
            _inputs.A[phase_num] *= static_cast<float>(deltat / deltax);
            _inputs.B[phase_num] *= static_cast<float>(deltat / deltax);
            _inputs.C[phase_num] *= static_cast<float>(deltat / deltax);
        }
        else if (_inputs.function[phase_num] == _inputs.power) {
            // Normalize only the leading and last coefficient: V = A*x^B + C
            _inputs.A[phase_num] *= static_cast<float>(deltat / deltax);
            _inputs.C[phase_num] *= static_cast<float>(deltat / deltax);
        }
    }

    // Compute velocity from local undercooling.
    // functional form is assumed to be cubic if not explicitly given in input file
    KOKKOS_INLINE_FUNCTION
    float compute(const float loc_u, const int phase_num = 0) const {
        float V;
        assert(phase_num < _inputs.num_phases);
        if (_inputs.function[phase_num] == _inputs.quadratic)
            V = _inputs.A[phase_num] * Kokkos::pow(loc_u, 2.0) + _inputs.B[phase_num] * loc_u + _inputs.C[phase_num];
        else if (_inputs.function[phase_num] == _inputs.power)
            V = _inputs.A[phase_num] * Kokkos::pow(loc_u, _inputs.B[phase_num]) + _inputs.C[phase_num];
        else
            V = _inputs.A[phase_num] * Kokkos::pow(loc_u, 3.0) + _inputs.B[phase_num] * Kokkos::pow(loc_u, 2.0) +
                _inputs.C[phase_num] * loc_u + _inputs.D[phase_num];
        return Kokkos::fmax(0.0, V);
    }

    std::string functionName(const int phase_num = 0) {
        // Not storing string due to Cuda warnings when constructing on device.
        if (_inputs.function[phase_num] == _inputs.cubic)
            return "cubic";
        else if (_inputs.function[phase_num] == _inputs.quadratic)
            return "quadratic";
        else if (_inputs.function[phase_num] == _inputs.power)
            return "power";

        // Should never make it here
        return "none";
    }

    // Get the solid phase based on the faster growing phase at the input undercooling loc_u. If this is a single phase
    // problem, return 0 for the default phase
    KOKKOS_INLINE_FUNCTION
    int getPreferredPhase_Nucleation(const float loc_u) const {
        if (_inputs.num_phases == 1)
            return 0;
        else {
            const double V_p0 = compute(loc_u, 0);
            const double V_p1 = compute(loc_u, 1);
            if (V_p0 > V_p1)
                return 0;
            else
                return 1;
        }
    }

    // Get the solid phase based on the number of phases in the material and the input transformation rule
    KOKKOS_INLINE_FUNCTION
    int getPreferredPhase_Activation(const int local_phase_id) const {
        if (_inputs.transformation == _inputs.solidification)
            return 0;
        else
            return local_phase_id;
    }

    auto A(const int phase_num = 0) { return _inputs.A[phase_num]; }
    auto B(const int phase_num = 0) { return _inputs.B[phase_num]; }
    auto C(const int phase_num = 0) { return _inputs.C[phase_num]; }
    auto D(const int phase_num = 0) { return _inputs.D[phase_num]; }
    auto freezingRange(const int phase_num = 0) { return _inputs.freezing_range[phase_num]; }
    auto num_phases() { return _inputs.num_phases; }
};

#endif
