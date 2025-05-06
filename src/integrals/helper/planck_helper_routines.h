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
#include "../../math/planck_math.h"
#include "../../lookup/planck_lookup.h"
#include "../huzinaga/planck_huzinaga.h"

void computeGaussianProduct(cxx_Contracted *contractedGaussianA, cxx_Contracted *contractedGaussianB, std::vector<cxx_Gaussians> *productGaussians);
cxx_Gaussians computeGaussianProduct(cxx_Primitive primitiveA, const std::double_t xA, const std::double_t yA, const std::double_t zA, cxx_Primitive primitiveB, const std::double_t xB, const std::double_t yB, const std::double_t zB);
std::double_t boysFunction(std::uint64_t boysIndex, std::double_t boysParam);
std::vector<eriKet> schwatrzSceening(cxx_Calculator *planckCalculator, Eigen::Tensor<std::double_t, 4> &electronicMatrix);
