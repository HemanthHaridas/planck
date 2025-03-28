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
#include "helper/planck_helper_routines.h"

void computeOverlap(cxx_Calculator *planckCalculator, std::error_code *errorFlag, std::string *errorMessage)
{
    // first clear the error buffers
    errorFlag->clear();
    errorMessage->clear();

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
                std::double_t exponentA = planckCalculator->calculation_set[ii].contracted_GTO[ij].primitive_exp;
                for (std::uint64_t ji = 0; ji < nPrimB; ji++)
                {
                    // need six variables to hold intermediates
                    std::double_t *gaussianCenterX;
                    std::double_t *gaussianCenterY;
                    std::double_t *gaussianCenterZ;
                    std::double_t *gaussianIntegralX;
                    std::double_t *gaussianIntegralY;
                    std::double_t *gaussianIntegralZ;

                    std::double_t exponentB = planckCalculator->calculation_set[jj].contracted_GTO[ji].primitive_exp;
                    computeGaussianProduct(gaussianCenterX, gaussianCenterY, gaussianCenterZ, gaussianIntegralX, gaussianIntegralY, gaussianIntegralZ, errorFlag, errorMessage);
                }
            }
        }
    }
}