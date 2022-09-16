// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_HALO_HPP
#define EXACA_HALO_HPP

#include "CAtypes.hpp"
#include <Kokkos_Core.hpp>

#include "mpi.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Struct specifying ExaCA buffers for halo exchange operations

struct Halo {

  public:
    ViewBuffer BufferSouthSend, BufferSouthRecv, BufferNorthSend, BufferNorthRecv;
    // Constructor
    Halo(const int BufSizeX, const int BufSizeZ);
    // Load data (GrainID, DOCenter, DiagonalLength) into ghost nodes if the given RankY is associated with a 1D halo
    // region
    KOKKOS_INLINE_FUNCTION void loadghostnodes(const int GhostGID, const float GhostDOCX, const float GhostDOCY,
                                               const float GhostDOCZ, const float GhostDL, const int BufSizeX,
                                               const int MyYSlices, const int RankX, const int RankY, const int RankZ,
                                               const bool AtNorthBoundary, const bool AtSouthBoundary) const {

        if ((RankY == 1) && (!(AtSouthBoundary))) {
            int GNPosition = RankZ * BufSizeX + RankX;
            BufferSouthSend(GNPosition).GrainID = GhostGID;
            BufferSouthSend(GNPosition).DOCenterX = GhostDOCX;
            BufferSouthSend(GNPosition).DOCenterY = GhostDOCY;
            BufferSouthSend(GNPosition).DOCenterZ = GhostDOCZ;
            BufferSouthSend(GNPosition).DiagonalLength = GhostDL;
        }
        else if ((RankY == MyYSlices - 2) && (!(AtNorthBoundary))) {
            int GNPosition = RankZ * BufSizeX + RankX;
            BufferNorthSend(GNPosition).GrainID = GhostGID;
            BufferNorthSend(GNPosition).DOCenterX = GhostDOCX;
            BufferNorthSend(GNPosition).DOCenterY = GhostDOCY;
            BufferNorthSend(GNPosition).DOCenterZ = GhostDOCZ;
            BufferNorthSend(GNPosition).DiagonalLength = GhostDL;
        }
    }
    // Reset to zeros
    void reset(const int Buffersize);
    // Destructor
    ~Halo();
};

#endif
