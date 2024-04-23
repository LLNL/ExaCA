// Copyright 2021-2024 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_ORIENTATION_HPP
#define EXACA_ORIENTATION_HPP

#include "CAparsefiles.hpp"

#include "mpi.h"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>

// Data representing the crystallographic orientations corresponding to GrainID values
template <typename MemorySpace>
struct Orientation {

    using memory_space = MemorySpace;
    using view_type_float = Kokkos::View<float *, memory_space>;
    using view_type_float_host = typename view_type_float::HostMirror;

    // Using the default exec space for this memory space.
    using execution_space = typename memory_space::execution_space;

    // Grain orientation data views - either rotation matrices (9 vals per orientation), Euler angles (3 vals per
    // orientation), or RGB values (3 vals per orientation). Unit vectors are stored on the device and no host copy is
    // maintained in the struct, Euler and RGB representations are only used on the host and no device copy is
    // maintained
    int n_grain_orientations;
    view_type_float grain_unit_vector;
    view_type_float_host grain_bunge_euler_host, grain_rgb_ipfz_host;

    Orientation(const int id, const std::string grain_unit_vector_file, const bool init_euler_rgb_vals)
        : grain_unit_vector(view_type_float(Kokkos::ViewAllocateWithoutInitializing("grain_unit_vector"), 1))
        , grain_bunge_euler_host(
              view_type_float_host(Kokkos::ViewAllocateWithoutInitializing("grain_bunge_euler_host"), 1))
        , grain_rgb_ipfz_host(view_type_float_host(Kokkos::ViewAllocateWithoutInitializing("grain_rgb_ipfz_host"), 1)) {

        // Get unit vectors from the grain orientations file (9 vals per line) and store in temporary host view
        view_type_float_host grain_unit_vector_host_ = getOrientationsFromFile(grain_unit_vector_file, 9, false);
        // Copy unit vector data from temporary host view
        grain_unit_vector = Kokkos::create_mirror_view_and_copy(memory_space(), grain_unit_vector_host_);

        // If needed, repeat this process for Euler and RGB representation of the orientations
        if (init_euler_rgb_vals) {
            // Other orientation files should start with the strings GrainOrientationEulerAnglesBungeZXZ and
            // GrainOrientationRGB_IPF-Z, but be proceeded with a unique string (such as _1e6.csv) that matches the
            // corresponding file of grain unit vectors. In the case of using the default file
            // GrainOrientationVectors.csv, this will be an empty string
            getEulerIPFZFilenames(grain_unit_vector_file);
        }
        if (id == 0)
            std::cout << "Done with orientation initialization " << std::endl;
    }

    // Get grain orientations from the specified file and the number of grain orientations and return a host view
    // containing the data
    view_type_float_host getOrientationsFromFile(const std::string orientation_file, const int vals_per_line,
                                                 const bool check_n_orientations) {

        // Read file of grain orientations
        std::ifstream orientation_input_stream;
        orientation_input_stream.open(orientation_file);

        // Line 1 is the number of orientation values to read (if check_n_orientations is false, store this value as
        // n_orientations in the struct. Otherwise, read this value and check that it matches n_orientations)
        std::string n_grain_orientations_read;
        getline(orientation_input_stream, n_grain_orientations_read);
        if (check_n_orientations) {
            int n_grain_orientations_check = getInputInt(n_grain_orientations_read);
            if (n_grain_orientations != n_grain_orientations_check)
                throw std::runtime_error(
                    "Error: The number of orientations given on the first line of all orientation files must match");
        }
        else
            n_grain_orientations = getInputInt(n_grain_orientations_read);

        // Resize view for storing grain orientations read from file based on the known value for n_grain_orientations
        view_type_float_host orientation_data_host(Kokkos::ViewAllocateWithoutInitializing("orientation_data_host"),
                                                   vals_per_line * n_grain_orientations);

        // Populate data structure for grain orientation data and return the view
        for (int i = 0; i < n_grain_orientations; i++) {
            std::vector<std::string> parsed_line(vals_per_line);
            std::string read_line;
            if (!getline(orientation_input_stream, read_line))
                break;
            splitString(read_line, parsed_line, vals_per_line);
            // Place the 3 grain orientation angles or 9 rotation matrix components into the orientation data view
            for (int comp = 0; comp < vals_per_line; comp++) {
                orientation_data_host(vals_per_line * i + comp) = getInputFloat(parsed_line[comp]);
            }
        }
        orientation_input_stream.close();
        return orientation_data_host;
    }

    // Get names of and read other orientation files
    void getEulerIPFZFilenames(const std::string grain_unit_vector_file) {

        std::string euler_angles_filename, rgb_filename;
        if (grain_unit_vector_file.find("GrainOrientationVectors.csv") != std::string::npos) {
            // Default files for euler angles and RGB mapping
            euler_angles_filename = "GrainOrientationEulerAnglesBungeZXZ.csv";
            rgb_filename = "GrainOrientationRGB_IPF-Z.csv";
        }
        else {
            // Custom files for euler angles and RGB mapping based on rotation filename
            std::size_t baseorientation_startpos = grain_unit_vector_file.find_last_of("/");
            std::size_t endpos = grain_unit_vector_file.find_last_of(".");
            std::string pathtoorientations = grain_unit_vector_file.substr(0, baseorientation_startpos + 1);
            std::string baseorientationname =
                grain_unit_vector_file.substr(baseorientation_startpos + 1, endpos - baseorientation_startpos - 1);
            std::size_t customname_startpos = baseorientationname.find_last_of("_");
            std::string customname =
                baseorientationname.substr(customname_startpos + 1, endpos - customname_startpos - 1);
            euler_angles_filename = pathtoorientations + "GrainOrientationEulerAnglesBungeZXZ_" + customname + ".csv";
            rgb_filename = pathtoorientations + "GrainOrientationRGB_IPF-Z_" + customname + ".csv";
        }
        // Full path to install location was already given for the rotations file, need to check install location for
        // other files and ensure they are not empty
        euler_angles_filename = checkFileInstalled(euler_angles_filename, 0);
        rgb_filename = checkFileInstalled(rgb_filename, 0);
        checkFileNotEmpty(euler_angles_filename);
        checkFileNotEmpty(rgb_filename);

        // Read and store additional orientation data in the appropriate views
        grain_bunge_euler_host = getOrientationsFromFile(euler_angles_filename, 3, true);
        grain_rgb_ipfz_host = getOrientationsFromFile(rgb_filename, 3, true);
    }

    // Create a view of size "NumberOfOrientations" of the misorientation of each possible grain orientation with the X,
    // Y, or Z directions (dir = 0, 1, or 2, respectively)
    view_type_float_host misorientationCalc(const int dir) {

        view_type_float_host grain_misorientation(Kokkos::ViewAllocateWithoutInitializing("GrainMisorientation"),
                                                  n_grain_orientations);
        view_type_float_host grain_unit_vector_host =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), grain_unit_vector);

        // Find the smallest possible misorientation between the specified direction, and this grain orientations' 6
        // possible 001 directions (where 54.7356 degrees is the largest possible misorientation between a 001 and a
        // given cardinal direction for the cubic crystal system)
        for (int n = 0; n < n_grain_orientations; n++) {
            float misorientation_angle_min = 54.7356;
            for (int ll = 0; ll < 3; ll++) {
                float misorientation =
                    Kokkos::abs((180 / M_PI) * Kokkos::acos(Kokkos::abs(grain_unit_vector_host(9 * n + 3 * ll + dir))));
                if (misorientation < misorientation_angle_min) {
                    misorientation_angle_min = misorientation;
                }
            }
            grain_misorientation(n) = misorientation_angle_min;
        }
        return grain_misorientation;
    }
};

// Inline functions
// Get the grain ID from the repeat number of a grain from a given grain ID and the number of possible orientations
KOKKOS_INLINE_FUNCTION int getGrainID(const int my_grain_orientation, const int my_grain_number,
                                      const int n_grain_orientations) {
    int my_grain_id = n_grain_orientations * (Kokkos::abs(my_grain_number) - 1) + my_grain_orientation;
    if (my_grain_number < 0)
        my_grain_id = -my_grain_id;
    return my_grain_id;
}

// Get the orientation of a grain from a given grain ID and the number of possible orientations
// By default, start indexing at 0 (GrainID of 1 has Orientation number 0), optionally starting at 1
// GrainID of 0 is a special case - has orientation 0 no matter what (only used when reconstructing grain ID in
// getGrainID)
KOKKOS_INLINE_FUNCTION int getGrainOrientation(const int my_grain_id, const int n_grain_orientations,
                                               bool start_at_zero = true) {
    int my_orientation;
    if (my_grain_id == 0)
        my_orientation = 0;
    else {
        my_orientation = (Kokkos::abs(my_grain_id) - 1) % n_grain_orientations;
        if (!(start_at_zero))
            my_orientation++;
    }
    return my_orientation;
}

// Get the repeat number for the orientation of a grain with a given grain ID
// 1, 2, 3... or -1, -2, -3...
KOKKOS_INLINE_FUNCTION int getGrainNumber(int my_grain_id, int n_grain_orientations) {
    int my_grain_number = (Kokkos::abs(my_grain_id) - 1) / n_grain_orientations + 1;
    if (my_grain_id < 0)
        my_grain_number = -my_grain_number;
    return my_grain_number;
}

#endif
