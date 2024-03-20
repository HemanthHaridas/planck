#pragma once

#include "../molecule/molecule.h"

/*-----------------------------------------------------------------------------
 * Planck
 * Copyright (C) 2024 Hemanth Haridas, University of Utah
 * Contact: hemanthhari23@gmail.com
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or a later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 ----------------------------------------------------------------------------*/

void gaussianProducts(cxx_Primitives *primtiveGTO_a, cxx_Primitives *primtiveGTO_b, cxx_gptResults *gptResults);
void overlapPrimitives(cxx_Primitives *primitiveGTO_a, cxx_Primitives *primitiveGTO_b, cxx_gptResults *gptResults, std::uint64_t indexA, std::uint64_t indexB, cxx_Integral *integralResult);
// void overlapCartesians(cxx_Calculator *scfCalculator);
void overlap(cxx_Calculator *scfCalculator);