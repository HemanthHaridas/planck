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

#include "planck_helper_routines.h"

void computeGaussianProduct(cxx_Calculator *planckCalculator, std::error_code *errorFlag, std::string *errorMessage)
{
    // first clear the error buffers
    errorFlag->clear();
    errorMessage->clear();

    // now set the size of the buffers
    planckCalculator->gaussian_centers = (std::double_t *)malloc(sizeof(std::double_t) * planckCalculator->total_primitives * planckCalculator->total_primitives);
    planckCalculator->gaussian_exps = (std::double_t *)malloc(sizeof(std::double_t) * planckCalculator->total_primitives * planckCalculator->total_primitives);

    // now compute the gaussian centers and products
    for (std::uint64_t ii = 0; ii < planckCalculator->total_basis; ii++)
    {
        for (std::uint64_t jj = 0; jj < planckCalculator->total_basis; jj++)
        {
            /* code */
        }
    }
    
}