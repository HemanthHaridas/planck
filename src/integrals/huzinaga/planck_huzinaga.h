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

#include <numeric>
#include <system_error>
#include "../../base/planck_base.h"

namespace Huzinaga
{
    std::double_t computeNuclear(std::double_t *atomCoords, std::uint64_t *atomicNumbers, std::uint64_t nAtoms, cxx_Contracted *contractedGaussianA, cxx_Contracted *contractedGaussianB, std::error_code *errorFlag, std::string *errorMessage);
    std::double_t computeKinetic(cxx_Contracted *contractedGaussianA, cxx_Contracted *contractedGaussianB, std::error_code *errorFlag, std::string *errorMessage);
    std::double_t computeOverlap(cxx_Contracted *contractedGaussianA, cxx_Contracted *contractedGaussianB, std::error_code *errorFlag, std::string *errorMessage);

    void computePrimitive(cxx_Primitive *primitiveA, std::double_t xA, std::double_t yA, std::double_t zA, std::int64_t lxA, std::int64_t lyA, std::int64_t lzA, cxx_Primitive *primitiveB, std::double_t xB, std::double_t yB, std::double_t zB, std::int64_t lxB, std::int64_t lyB, std::int64_t lzB, cxx_Gaussians *productGaussian, std::double_t *primitiveOverlaps);
    void computeNuclear(std::double_t *atomCoords, std::uint64_t *atomicNumbers, std::uint64_t nAtoms, cxx_Primitive *primitiveA, std::double_t xA, std::double_t yA, std::double_t zA, std::int64_t lxA, std::int64_t lyA, std::int64_t lzA, cxx_Primitive *primitiveB, std::double_t xB, std::double_t yB, std::double_t zB, std::int64_t lxB, std::int64_t lyB, std::int64_t lzB, cxx_Gaussians *productGaussian, std::double_t *primitiveOverlaps);
    
    void computeNuclearPrimitive(std::double_t atomCoord, std::double_t coordA, std::int64_t shellA, std::double_t coordB, std::int64_t shellB, std::double_t gammaExp, std::double_t gaussCoord, std::vector <std::double_t> *primitiveResults);
    void computeNuclearPrimitive1(std::int64_t expandIndex, std::int64_t shellA, std::double_t coordA, std::int64_t shellB, std::double_t coordB, std::double_t gaussCoord, std::double_t *result);
    void computeNuclearPrimitive2(std::int64_t indexA, std::int64_t indexB, std::int64_t indexC, std::int64_t shellA, std::double_t coordA, std::int64_t shellB, std::double_t coordB, std::double_t atomCoord, std::double_t gaussCoord, std::double_t gammaExp, std::double_t *result);
}