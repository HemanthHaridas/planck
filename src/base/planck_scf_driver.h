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

#include "../integrals/helper/planck_helper_routines.h"
#include "../integrals/huzinaga/planck_huzinaga.h"

// #include "../integrals/planck_integrals.h"

// void scfEngine(cxx_Calculator *planckCalculator, std::error_code *errorFlag, std::string *errorMessage);
struct scfEngine: private cxx_Calculator, private cxx_Molecule
{
public:
    std::uint64_t cycle = 0;
    std::error_code errorFlag;
    std::string errorMessage;
    std::double_t scfEnergy;
    // std::vector <eriShell> eriShells;

public:
    std::int64_t scfCycle();

private:
    std::int64_t computeOverlap();
    std::int64_t computeKinetic();
    std::int64_t computeNuclear();
    std::int64_t computeElectronic();
    std::int64_t schwartzScreening();

    std::int64_t generateCoreHamiltonian();
    std::int64_t generateFockMatrix();
    std::int64_t performDIIS();
};
