// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_GHOST_HPP
#define EXACA_GHOST_HPP

#include "CAcelldata.hpp"
#include "CAfunctions.hpp"
#include "CAtypes.hpp"

#include <Kokkos_Core.hpp>

// Load data (GrainID, DOCenter, DiagonalLength) into ghost nodes if the given RankY is associated with a 1D halo region
// Uses check to ensure that the buffer position does not reach the buffer size - if it does, return false (otherwise
// return true) but keep incrementing the send size counters for use resizing the buffers in the future
KOKKOS_INLINE_FUNCTION bool loadghostnodes(const int GhostGID, const float GhostDOCX, const float GhostDOCY,
                                           const float GhostDOCZ, const float GhostDL, ViewI SendSizeNorth,
                                           ViewI SendSizeSouth, const int ny_local, const int coord_x,
                                           const int coord_y, const int coord_z, const bool AtNorthBoundary,
                                           const bool AtSouthBoundary, Buffer2D BufferSouthSend,
                                           Buffer2D BufferNorthSend, int NGrainOrientations, int BufSize) {
    bool DataFitsInBuffer = true;
    if ((coord_y == 1) && (!(AtSouthBoundary))) {
        int GNPositionSouth = Kokkos::atomic_fetch_add(&SendSizeSouth(0), 1);
        if (GNPositionSouth >= BufSize)
            DataFitsInBuffer = false;
        else {
            BufferSouthSend(GNPositionSouth, 0) = static_cast<float>(coord_x);
            BufferSouthSend(GNPositionSouth, 1) = static_cast<float>(coord_z);
            BufferSouthSend(GNPositionSouth, 2) =
                static_cast<float>(getGrainOrientation(GhostGID, NGrainOrientations, false));
            BufferSouthSend(GNPositionSouth, 3) = static_cast<float>(getGrainNumber(GhostGID, NGrainOrientations));
            BufferSouthSend(GNPositionSouth, 4) = GhostDOCX;
            BufferSouthSend(GNPositionSouth, 5) = GhostDOCY;
            BufferSouthSend(GNPositionSouth, 6) = GhostDOCZ;
            BufferSouthSend(GNPositionSouth, 7) = GhostDL;
        }
    }
    else if ((coord_y == ny_local - 2) && (!(AtNorthBoundary))) {
        int GNPositionNorth = Kokkos::atomic_fetch_add(&SendSizeNorth(0), 1);
        if (GNPositionNorth >= BufSize)
            DataFitsInBuffer = false;
        else {
            BufferNorthSend(GNPositionNorth, 0) = static_cast<float>(coord_x);
            BufferNorthSend(GNPositionNorth, 1) = static_cast<float>(coord_z);
            BufferNorthSend(GNPositionNorth, 2) =
                static_cast<float>(getGrainOrientation(GhostGID, NGrainOrientations, false));
            BufferNorthSend(GNPositionNorth, 3) = static_cast<float>(getGrainNumber(GhostGID, NGrainOrientations));
            BufferNorthSend(GNPositionNorth, 4) = GhostDOCX;
            BufferNorthSend(GNPositionNorth, 5) = GhostDOCY;
            BufferNorthSend(GNPositionNorth, 6) = GhostDOCZ;
            BufferNorthSend(GNPositionNorth, 7) = GhostDL;
        }
    }
    return DataFitsInBuffer;
}
void ResetSendBuffers(int BufSize, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, ViewI SendSizeNorth,
                      ViewI SendSizeSouth);
int ResizeBuffers(Buffer2D &BufferNorthSend, Buffer2D &BufferSouthSend, Buffer2D &BufferNorthRecv,
                  Buffer2D &BufferSouthRecv, ViewI SendSizeNorth, ViewI SendSizeSouth, ViewI_H SendSizeNorth_Host,
                  ViewI_H SendSizeSouth_Host, int OldBufSize, int NumCellsBufferPadding = 25);
void ResetBufferCapacity(Buffer2D &BufferNorthSend, Buffer2D &BufferSouthSend, Buffer2D &BufferNorthRecv,
                         Buffer2D &BufferSouthRecv, int NewBufSize);
void RefillBuffers(int nx, int nz_layer, int ny_local, CellData<device_memory_space> &cellData,
                   Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, ViewI SendSizeNorth, ViewI SendSizeSouth,
                   bool AtNorthBoundary, bool AtSouthBoundary, ViewF DOCenter, ViewF DiagonalLength,
                   int NGrainOrientations, int BufSize);
void GhostNodes1D(int, int, int NeighborRank_North, int NeighborRank_South, int nx, int ny_local, int y_offset,
                  NList NeighborX, NList NeighborY, NList NeighborZ, CellData<device_memory_space> &cellData,
                  ViewF DOCenter, ViewF GrainUnitVector, ViewF DiagonalLength, ViewF CritDiagonalLength,
                  int NGrainOrientations, Buffer2D BufferNorthSend, Buffer2D BufferSouthSend, Buffer2D BufferNorthRecv,
                  Buffer2D BufferSouthRecv, int BufSize, ViewI SendSizeNorth, ViewI SendSizeSouth);

#endif
