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

#include <iostream>
#include <iomanip>
#include "planck_integrals.h"

void computeOverlap(cxx_Calculator *planckCalculator, std::error_code *errorFlag, std::string *errorMessage)
{
    // first clear the error buffers
    errorFlag->clear();
    errorMessage->clear();

    // now loop over contracted basis sets
    for (std::uint64_t ii = 0; ii < planckCalculator->total_basis; ii++)
    {
        for (std::uint64_t jj = 0; jj < planckCalculator->total_basis; jj++)
        {
            planckCalculator->overlap[ii * planckCalculator->total_basis + jj] = Huzinaga::computeOverlap(&planckCalculator->calculation_set[ii], &planckCalculator->calculation_set[jj], errorFlag, errorMessage);
            // std::cout << std::setw(10) << std::left << ii << std::setw(10) << std::left << jj << std::setw(5) << std::setprecision(3) << std::right << planckCalculator->overlap[ii * planckCalculator->total_basis + jj] << "\n";
        }
        // std::cout << "\n";
    }
}

void computeKinetic(cxx_Calculator *planckCalculator, std::error_code *errorFlag, std::string *errorMessage)
{
    // first clear the error buffers
    errorFlag->clear();
    errorMessage->clear();

   // now loop over contracted basis sets
    for (std::uint64_t ii = 0; ii < planckCalculator->total_basis; ii++)
    {
        for (std::uint64_t jj = 0; jj < planckCalculator->total_basis; jj++)
        {
            planckCalculator->kinetic[ii * planckCalculator->total_basis + jj] = Huzinaga::computeKinetic(&planckCalculator->calculation_set[ii], &planckCalculator->calculation_set[jj], errorFlag, errorMessage);
            // std::cout << std::setw(10) << std::left << ii << std::setw(10) << std::left << jj << std::setw(5) << std::setprecision(2) << std::right << planckCalculator->kinetic[ii * planckCalculator->total_basis + jj] << "\n";
        }
        // std::cout << "\n";
    }    
}

void computeNuclear(std::double_t *atomCoords, std::uint64_t  *atomNumbers, std::uint64_t nAtoms, cxx_Calculator *planckCalculator, std::error_code *errorFlag, std::string *errorMessage)
{
    // first clear the error buffers
    errorFlag->clear();
    errorMessage->clear();

    // now loop over contracted basis sets
    for (std::uint64_t ii = 0; ii < planckCalculator->total_basis; ii++)
    {
        for (std::uint64_t jj = 0; jj < planckCalculator->total_basis; jj++)
        {
            planckCalculator->nuclear[ii * planckCalculator->total_basis + jj] = Huzinaga::computeNuclear(atomCoords, atomNumbers, nAtoms, &planckCalculator->calculation_set[ii], &planckCalculator->calculation_set[jj], errorFlag, errorMessage);
            std::cout << std::setw(10) << std::left << ii << std::setw(10) << std::left << jj << std::setw(10) << std::setprecision(2) << std::right << planckCalculator->nuclear[ii * planckCalculator->total_basis + jj] << "\n";
        }
        std::cout << "\n";
    }    
}