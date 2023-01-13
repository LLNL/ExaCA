// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "GAutils.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>

// Given an input vector of integer Grain ID values, return an output vector consisting of the unique Grain ID values,
// sorted from lowest to highest. Store and print the number of grains to the console and to the output file stream
std::vector<int> getUniqueGrains(std::ofstream &QoIs, const std::vector<int> GrainIDVector, int &NumberOfGrains) {
    std::vector<int> UniqueGrainIDVector = GrainIDVector;
    std::sort(UniqueGrainIDVector.begin(), UniqueGrainIDVector.end());
    std::vector<int>::iterator it;
    it = std::unique(UniqueGrainIDVector.begin(), UniqueGrainIDVector.end());
    UniqueGrainIDVector.resize(std::distance(UniqueGrainIDVector.begin(), it));
    std::cout << "The number of grains is " << NumberOfGrains << std::endl;
    QoIs << "The number of grains is " << NumberOfGrains << std::endl;
    return UniqueGrainIDVector;
}

// Given an input vector of integer Grain ID values "GrainIDVector" (of size "regionSize"), and an input vector of the
// unique Grain ID values "UniqueGrainIDVector" (of size "numberOfGrains", which is less than or equal to "regionSize"),
// return a third vector "GrainSizes" listing the size of each of the "numberOfGrains" grains
std::vector<int> getGrainSizes(const std::vector<int> GrainIDVector, const std::vector<int> UniqueGrainIDVector,
                               const int numberOfGrains) {

    std::vector<int> GrainSizeVector(numberOfGrains) = {0};
    for (int nn = 0; nn < numberOfGrains; nn++) {
        GrainSizeVector[nn] = std::count(GrainIDVector.begin(), GrainIDVector.end(), UniqueGrainIDVector[nn]);
    }
    return GrainSizeVector;
}

// Given the 3D grain structure "GrainID", determine the extent in the direction specified of each of the
// "numberOfGrains" unique grains from the volume bounded by [XLow, XHigh], [YLow, YHigh], [ZLow, ZHigh]
std::vector<int> getGrainLengths(ViewI3D_H GrainID, const std::vector<int> UniqueGrainIDVector,
                                 std::vector<int> GrainSizeVector, const int numberOfGrains, const int XLow,
                                 const int XHigh, const int YLow, const int YHigh, const int ZLow, const int ZHigh,
                                 std::string Direction) {

    std::vector<int> GrainLength(numberOfGrains);
    for (int n = 0; n < numberOfGrains; n++) {
        int grainSize = GrainSizeVector[n];
        std::vector<int> GrainCoordinate(grainSize);
        int count = 0;
        for (int k = ZLow; k <= ZHigh; k++) {
            for (int i = XLow; i <= XHigh; i++) {
                for (int j = YLow; j <= YHigh; j++) {
                    if (GrainID(k, i, j) == UniqueGrainIDVector[n]) {
                        if (Direction == "X")
                            GrainCoordinate[count] = i;
                        else if (Direction == "Y")
                            GrainCoordinate[count] = j;
                        else if (Direction == "Z")
                            GrainCoordinate[count] = k;
                        count++;
                    }
                }
            }
        }
    }
    return GrainLength;
}
