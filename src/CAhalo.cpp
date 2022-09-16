// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "CAhalo.hpp"
#include <Kokkos_Core.hpp>

//-constructor: buffers for ghost node data should be initialized with zeros
Halo::Halo(const int BufSizeX, const int BufSizeZ) {

    // Buffers for sending data in the +/- Y directions (used for 1D decomposition)
    BufferSouthSend = ViewBuffer("BufferSouthSend", BufSizeX * BufSizeZ);
    BufferSouthRecv = ViewBuffer("BufferSouthRecv", BufSizeX * BufSizeZ);
    BufferNorthSend = ViewBuffer("BufferNorthSend", BufSizeX * BufSizeZ);
    BufferNorthRecv = ViewBuffer("BufferNorthRecv", BufSizeX * BufSizeZ);
}

// Reset elements in the buffer views to zeros
void Halo::reset(const int Buffersize) {

    Kokkos::parallel_for(
        "ViewReset", Buffersize, KOKKOS_LAMBDA(const int &i) {
            BufferSouthSend(i).GrainID = 0;
            BufferSouthRecv(i).GrainID = 0;
            BufferNorthSend(i).GrainID = 0;
            BufferNorthRecv(i).GrainID = 0;
            BufferSouthSend(i).DOCenterX = 0;
            BufferSouthRecv(i).DOCenterX = 0;
            BufferNorthSend(i).DOCenterX = 0;
            BufferNorthRecv(i).DOCenterX = 0;
            BufferSouthSend(i).DOCenterY = 0;
            BufferSouthRecv(i).DOCenterY = 0;
            BufferNorthSend(i).DOCenterY = 0;
            BufferNorthRecv(i).DOCenterY = 0;
            BufferSouthSend(i).DOCenterZ = 0;
            BufferSouthRecv(i).DOCenterZ = 0;
            BufferNorthSend(i).DOCenterZ = 0;
            BufferNorthRecv(i).DOCenterZ = 0;
            BufferSouthSend(i).DiagonalLength = 0;
            BufferSouthRecv(i).DiagonalLength = 0;
            BufferNorthSend(i).DiagonalLength = 0;
            BufferNorthRecv(i).DiagonalLength = 0;
        });
}

Halo::~Halo() {}
