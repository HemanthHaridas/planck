#pragma once

#include <system_error>

#include "../external/libmsym/install/include/libmsym/msym.h"
#include "../external/libmsym/install/include/libmsym/msym_EXPORTS.h"
#include "../external/libmsym/install/include/libmsym/msym_error.h"
#include "planck_base.h"

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

// void calculateInertia(cxx_Molecule *inputMolecule, std::error_code *errorFlag, std::string *errorMessage);
void detectSymmetry(cxx_Molecule *inputMolecule, std::uint64_t nAtoms, std::error_code *errorFlag, std::string *errorMessage);