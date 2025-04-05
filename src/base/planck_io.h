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

#include <fstream>
#include <iostream>
#include <sstream>
#include <system_error>
#include <functional>
#include "planck_base.h"

void readInput(std::fstream *filePointer, cxx_Calculator *planckCalculator, cxx_Molecule *inputMolecule, std::error_code *errorFlag, std::string *errorMessage);
void dumpInput(cxx_Calculator *planckCalculator, cxx_Molecule *inputMolecule);
void dumpIntegral(cxx_Calculator *planckCalculator, cxx_Integrals integralType, std::string inputFile);