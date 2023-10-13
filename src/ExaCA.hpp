// Copyright 2021-2023 Lawrence Livermore National Security, LLC and other ExaCA Project Developers.
// See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef EXACA_HPP
#define EXACA_HPP

#include "CAfunctions.hpp"
#include "CAghostnodes.hpp"
#include "CAinfo.hpp"
#include "CAinitialize.hpp"
#include "CAinputs.hpp"
#include "CAinterfacialresponse.hpp"
#include "CAnucleation.hpp"
#include "CAparsefiles.hpp"
#include "CAprint.hpp"
#include "CAtemperature.hpp"
#include "CAtypes.hpp"
#include "CAupdate.hpp"

#include <string>

void RunProgram_Reduced(int id, int np, std::string InputFile);

#endif
