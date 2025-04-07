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
#include <cstring>
#include <system_error>
#include "../../base/planck_base.h"

namespace Huzinaga
{
    // main functions
    std::double_t computeNuclear(std::double_t *atomCoords, std::uint64_t *atomCharges, std::uint64_t nAtoms, cxx_Contracted *contractedGaussianA, cxx_Contracted *contractedGaussianB, std::error_code *errorFlag, std::string *errorMessage);
    std::double_t computeKinetic(cxx_Contracted *contractedGaussianA, cxx_Contracted *contractedGaussianB, std::error_code *errorFlag, std::string *errorMessage);
    std::double_t computeOverlap(cxx_Contracted *contractedGaussianA, cxx_Contracted *contractedGaussianB, std::error_code *errorFlag, std::string *errorMessage);

    // only for electron nuclear integrals
    // std::double_t computeNuclearPrimitive(std::double_t *atomCoords, std::uint64_t *atomicNumbers, std::uint64_t nAtoms, cxx_Primitive *primitiveA, std::double_t xA, std::double_t yA, std::double_t zA, std::int64_t lxA, std::int64_t lyA, std::int64_t lzA, cxx_Primitive *primitiveB, std::double_t xB, std::double_t yB, std::double_t zB, std::int64_t lxB, std::int64_t lyB, std::int64_t lzB, cxx_Gaussians *productGaussian);
    // void computeNuclearPrimitive0(std::double_t atomCoord, std::double_t coordA, std::int64_t shellA, std::double_t coordB, std::int64_t shellB, std::double_t gamma, std::double_t epsilon, std::double_t gaussCoord, std::vector <nuclearInt> *primitiveResults);
    // std::double_t computeNuclearPrimitive1(const std::int64_t expandIndex, const std::int64_t shellA, const std::double_t coordA, const std::int64_t shellB, const std::double_t coordB, const std::double_t gaussCoord);
    // std::double_t computeNuclearPrimitive2(const std::int64_t indexA, const std::int64_t indexB, const std::int64_t indexC, const std::int64_t shellA, const std::double_t coordA, const std::int64_t shellB, const std::double_t coordB, const std::double_t atomCoord, const std::double_t gaussCoord, const std::double_t gamma, const std::double_t epsilon);

    // for overlap and kinetic energy integrals
    void computePrimitive(cxx_Primitive *primitiveA, std::double_t xA, std::double_t yA, std::double_t zA, std::int64_t lxA, std::int64_t lyA, std::int64_t lzA, cxx_Primitive *primitiveB, std::double_t xB, std::double_t yB, std::double_t zB, std::int64_t lxB, std::int64_t lyB, std::int64_t lzB, cxx_Gaussians *productGaussian, std::double_t *primitiveOverlaps);
    std::double_t expansionCoeff2(const std::int64_t indexA, const std::int64_t indexB, const std::int64_t indexC, const std::int64_t shellA, const std::double_t centerA, const std::int64_t shellB, const std::double_t centerB, const std::double_t atomCenter, const std::double_t gaussCenter, std::double_t gamma);
    std::double_t expansionCoeff1(const std::int64_t expIndex, const std::int64_t shellA, const std::double_t centerA, const std::int64_t shellB, const std::double_t centerB, const std::double_t gaussCenter);
    std::double_t computePrimitive(cxx_Primitive *primitiveA, const std::double_t xA, const std::double_t yA, const std::double_t zA, const std::int64_t lxA, const std::int64_t lyA, const std::int64_t lzA, cxx_Primitive *primitiveB, const std::double_t xB, const std::double_t yB, const std::double_t zB, const std::int64_t lxB, const std::int64_t lyB, const std::int64_t lzB, const std::double_t *atomCoords, const std::uint64_t *atomCharges, const std::uint64_t nAtoms);
    // std::double_t computePrimitive(cxx_Primitive *primitiveA, cxx_Primitive *primitiveB, const std::double_t *atomCoords, const std::uint64_t *atomCharges, const std::uint64_t nAtoms);
}