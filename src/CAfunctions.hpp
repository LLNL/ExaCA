// Copyright 2021 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_FUNCTIONS_HPP
#define EXACA_FUNCTIONS_HPP

double CrossP1(double TestVec1[3], double TestVec2[3]);
double CrossP2(double TestVec1[3], double TestVec2[3]);
double CrossP3(double TestVec1[3], double TestVec2[3]);
int FindItBounds(int RankX, int RankY, int MyXSlices, int MyYSlices);
int MaxIndex(double TestVec3[6]);
int XMPSlicesCalc(int p, int nx, int ProcessorsInXDirection, int ProcessorsInYDirection, int DecompositionStrategy);
int XOffsetCalc(int p, int nx, int ProcessorsInXDirection, int ProcessorsInYDirection, int DecompositionStrategy);
int YMPSlicesCalc(int p, int ny, int ProcessorsInYDirection, int np, int DecompositionStrategy);
int YOffsetCalc(int p, int ny, int ProcessorsInYDirection, int np, int DecompositionStrategy);
void AddGhostNodes(int DecompositionStrategy, int NeighborRank_West, int NeighborRank_East, int NeighborRank_North,
                   int NeighborRank_South, int &XRemoteMPSlices, int &RemoteXOffset, int &YRemoteMPSlices,
                   int &RemoteYOffset);
double MaxVal(double TestVec3[6], int NVals);
void InitialDecomposition(int &DecompositionStrategy, int nx, int ny, int &ProcessorsInXDirection,
                          int &ProcessorsInYDirection, int id, int np, int &NeighborRank_North, int &NeighborRank_South,
                          int &NeighborRank_East, int &NeighborRank_West, int &NeighborRank_NorthEast,
                          int &NeighborRank_NorthWest, int &NeighborRank_SouthEast, int &NeighborRank_SouthWest,
                          bool &AtNorthBoundary, bool &AtSouthBoundary, bool &AtEastBoundary, bool &AtWestBoundary);
void XYLimitCalc(int &LLX, int &LLY, int &ULX, int &ULY, int MyXSlices, int MyYSlices, int NeighborRank_South,
                 int NeighborRank_North, int NeighborRank_East, int NeighborRank_West);

#endif
