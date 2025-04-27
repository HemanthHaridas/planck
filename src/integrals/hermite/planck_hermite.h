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
#include "../helper/planck_helper_routines.h"

namespace Hermite
{
    std::double_t expansionCoeff1(std::int64_t shell1, std::int64_t shell2, std::int64_t nodes, std::double_t centerA, std::double_t centerB, std::double_t expA, std::double_t expB);
};