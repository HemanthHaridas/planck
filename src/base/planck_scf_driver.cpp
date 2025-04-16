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

#include "planck_scf_driver.h"

// void scfEngine(cxx_Molecule *inputMolecule, cxx_Calculator *planckCalculator, std::error_code *errorFlag, std::string *errorMessage)
// {
//     // first compute overlap and kinetic integrals
//     computeOverlap(planckCalculator, errorFlag, errorMessage);
//     computeKinetic(planckCalculator, errorFlag, errorMessage);

//     // generate core hamiltonian

//     // now compute electron nuclear integrals
//     computeNuclear(inputMolecule->standard_coordinates, inputMolecule->atom_numbers, planckCalculator->total_atoms, planckCalculator, errorFlag, errorMessage);
    
//     // now compute electron electron integrals
//     computeElectronic(planckCalculator, errorFlag, errorMessage);
// }

std::int64_t scfEngine::computeOverlap(cxx_Calculator *planckCalculator)
{
    // increment the cycle
    this->cycle++;

    // reset error buffers
    this->errorFlag.clear();
    this->errorMessage.clear();

    // now check if the pointer to planckCalculator is valid
    if (planckCalculator == nullptr)
    {
        this->errorFlag = std::make_error_code(std::errc::invalid_argument);
        this->errorMessage = "Unable to find the calculator object.";
        return EXIT_FAILURE;
    }
    
    // now check if SCF cycles is within limit
    if (this->cycle >= planckCalculator->max_scf)
    {
        this->errorFlag = std::make_error_code(std::errc::value_too_large);
        this->errorMessage = "Maximum number of SCF cycles reached. Please check the results.";
        return EXIT_FAILURE;
    }
    
    // now calculate overlap and kinetic integrals
    this->computeOverlap(planckCalculator);
    this->computeKinetic(planckCalculator);
}