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

#include "planck_overlap.h"

void computeOverlap(cxx_Calculator *planckCalculator, std::error_code *errorFlag, std::string *errorMessage)
{
    // first clear the error buffers
    errorFlag->clear();
    errorMessage->clear();

    // now allocate buffers for overlap
    planckCalculator->overlap = (std::double_t *)malloc(sizeof(std::double_t) * planckCalculator->total_basis * planckCalculator->total_basis);

    // now loop over contracted basis sets
    for (std::uint64_t ii = 0; ii < planckCalculator->total_basis; ii++)
    {
        std::uint64_t nPrimA = planckCalculator->calculation_set[ii].contracted_GTO.size();
        std::double_t xA = planckCalculator->calculation_set[ii].location_x;
        std::double_t yA = planckCalculator->calculation_set[ii].location_y;
        std::double_t zA = planckCalculator->calculation_set[ii].location_z;
        for (std::uint64_t jj = 0; jj < planckCalculator->total_basis; jj++)
        {
            std::uint64_t nPrimB = planckCalculator->calculation_set[jj].contracted_GTO.size();
            std::double_t xB = planckCalculator->calculation_set[jj].location_x;
            std::double_t yB = planckCalculator->calculation_set[jj].location_y;
            std::double_t zB = planckCalculator->calculation_set[jj].location_z;
            // now loop over primitives
            for (std::uint64_t ij = 0; ij < nPrimA; ij++)
            {
                for (std::uint64_t ji = 0; ji < nPrimB; ji++)
                {
                }
            }
        }
    }
}