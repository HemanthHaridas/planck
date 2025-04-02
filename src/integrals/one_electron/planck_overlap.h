#pragma once

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

#include <system_error>

#include "../../base/planck_base.h"

std::double_t computeOverlap1(cxx_Contracted *contractedGaussianA, cxx_Contracted *contractedGaussianB, std::error_code *errorFlag, std::string *errorMessage);
void computeOverlap(cxx_Calculator *planckCalculator, std::error_code *errorFlag, std::string *errorMessage);
void computePrimitive_Huzinaga(cxx_Primitive *primitiveA, std::double_t *locA, std::int64_t *shellA, cxx_Primitive *primitiveB, std::double_t *locB, std::int64_t *shellB, cxx_Gaussians *productGaussian, std::double_t *primitiveOverlaps);
void computePrimitive_ObaraSaika(cxx_Primitive *primitiveA, std::double_t *locA, std::int64_t *shellA, cxx_Primitive *primitiveB, std::double_t *locB, std::int64_t *shellB, cxx_Gaussians *productGaussian, std::double_t *primitiveOverlaps);