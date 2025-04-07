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

#include <cmath>
#include <map>
#include <string>
#include <cstdint>

/*
 * This is a lookup table for atomic masses, atomic numbers and atomic radii for
 * elements from Hydrogen to Einsteinium.
 *
 * Data obtained from PubChem
 */

extern std::map<std::string, double> atomicMass;
extern std::map<std::string, std::uint64_t> atomicNumber;
extern std::map<std::string, std::double_t> atomicRadius;
extern const std::double_t boysTable[][66];