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

#include "../../base/planck_base.h"
#include "../../math/planck_math.h"
#include "../helper/planck_helper_routines.h"

namespace ObaraSakia
{
    namespace Overlap
    {
        std::double_t computePrimitive1D(std::double_t centerA, std::double_t exponentA, std::int64_t shellA, std::double_t centerB, std::double_t exponentB, std::int64_t shellB, std::double_t gaussianCenter, std::double_t gaussianExponent);
        std::double_t computePrimitive3D(cxx_Primitive primitiveA, std::double_t xA, std::double_t yA, std::double_t zA, std::int64_t lxA, std::int64_t lyA, std::int64_t lzA, cxx_Primitive primitiveB, std::double_t xB, std::double_t yB, std::double_t zB, std::int64_t lxB, std::int64_t lyB, std::int64_t lzB, std::double_t gaussianCenterX, std::double_t gaussianCenterY, std::double_t gaussianCenterZ, std::double_t gaussianExponent);
        std::double_t computeContracted(cxx_Contracted contractedGaussianA, cxx_Contracted contractedGaussianB);
    };
};