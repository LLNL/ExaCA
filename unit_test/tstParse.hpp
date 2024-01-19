// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <Kokkos_Core.hpp>

#include "CAconfig.hpp"
#include "CAparsefiles.hpp"
#include "CAprint.hpp"

#include <gtest/gtest.h>

#include "mpi.h"

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
// file_read_tests
//---------------------------------------------------------------------------//
void testReadWrite(bool print_read_binary) {

    // Make lists of some int and float data
    int int_data[5] = {-2, 0, 2, 4, 6};
    float float_data[5] = {-1.0, 0.0, 1.0, 2.0, 3.0};

    // Write data as binary to be used as input
    std::ofstream test_int_data;
    std::ofstream test_float_data;
    if (print_read_binary) {
        test_int_data.open("TestIntData.txt", std::ios::out | std::ios::binary);
        test_float_data.open("TestFloatData.txt", std::ios::out | std::ios::binary);
    }
    else {
        test_int_data.open("TestIntData.txt");
        test_float_data.open("TestFloatData.txt");
    }
    for (int n = 0; n < 5; n++) {
        // Write to files
        writeData(test_int_data, int_data[n], print_read_binary, true);
        writeData(test_float_data, float_data[n], print_read_binary, true);
    }
    test_int_data.close();
    test_float_data.close();

    // Read data and convert back to ints and floats, compare to original values
    std::ifstream test_int_data_read;
    test_int_data_read.open("TestIntData.txt");
    std::ifstream test_float_data_read;
    test_float_data_read.open("TestFloatData.txt");
    // For reading ASCII data, obtain the lines from the files first, then parse the string stream at the spaces
    if (print_read_binary) {
        for (int n = 0; n < 5; n++) {
            int int_to_compare = readBinaryData<int>(test_int_data_read, true);
            float float_to_compare = readBinaryData<float>(test_float_data_read, true);
            // Compare to expected values
            EXPECT_EQ(int_to_compare, int_data[n]);
            EXPECT_FLOAT_EQ(float_to_compare, float_data[n]);
        }
    }
    else {
        std::string intline, floatline;
        getline(test_int_data_read, intline);
        getline(test_float_data_read, floatline);
        std::istringstream intss(intline);
        std::istringstream floatss(floatline);
        for (int n = 0; n < 5; n++) {
            // Get values from string stream
            int int_to_compare = parseASCIIData<float>(intss);
            float float_to_compare = parseASCIIData<float>(floatss);
            // Compare to expected values
            EXPECT_EQ(int_to_compare, int_data[n]);
            EXPECT_FLOAT_EQ(float_to_compare, float_data[n]);
        }
    }
}
//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, parse) {
    // test functions for reading and writing data as binary (true) and ASCII (false)
    testReadWrite(true);
    testReadWrite(false);
}
} // end namespace Test
