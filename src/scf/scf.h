#pragma once

#include "../molecule/molecule.h"
#include "overlap.h"

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

void readXML(std::fstream *xmlPointer, cxx_Molecule *inputMolecule, cxx_Calculator *scfCalculator, std::error_code *errorFlag, std::string *errorMessage);
void writeXML_GPT(std::fstream *xmlPointer, const cxx_Calculator *scfCalculator, std::error_code *errorFlag, std::string *errorMessage);