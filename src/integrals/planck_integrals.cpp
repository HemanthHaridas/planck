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
#include "planck_integrals.h"

void IntegralEngine::computeOverlap(cxx_Calculator &planckCalculator, Eigen::MatrixXd &overlapMatrix)
{
    std::uint64_t nBasis = planckCalculator.total_basis;

    for (std::uint64_t row = 0; row < nBasis; row++)
    {
        for (std::uint64_t col = 0; col < nBasis; col++)
        {
            std::double_t value = Hermite::Overlap::computeContracted(planckCalculator.calculation_set[row], planckCalculator.calculation_set[col]);
            overlapMatrix(row, col) = value;
            overlapMatrix(col, row) = value;
        }
    }
}

// void IntegralEngine::computeKinetic(cxx_Calculator &planckCalculator, Eigen::MatrixXd &kineticMatrix)
// {
//     std::uint64_t nBasis = planckCalculator.total_basis;

//     for (std::uint64_t row = 0; row < nBasis; row++)
//     {
//         for (std::uint64_t col = 0; col <= row; col++)
//         {
//             std::double_t value = Huzinaga::Kinetic::computeContracted(planckCalculator.calculation_set[row], planckCalculator.calculation_set[col]);
//             kineticMatrix(row, col) = value;
//             kineticMatrix(col, row) = value;
//         }
//     }
// }