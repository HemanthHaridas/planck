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
#include "planck_helper_routines.h"

void computeGaussianProduct(cxx_Contracted *contractedGaussianA, cxx_Contracted *contractedGaussianB, std::vector<cxx_Gaussians> *productGaussians, std::error_code *errorFlag, std::string *errorMessage)
{
    // first clear the error buffers
    errorFlag->clear();
    errorMessage->clear();

    std::double_t xA = contractedGaussianA->location_x;
    std::double_t yA = contractedGaussianA->location_y;
    std::double_t zA = contractedGaussianA->location_z;

    std::double_t xB = contractedGaussianB->location_x;
    std::double_t yB = contractedGaussianB->location_y;
    std::double_t zB = contractedGaussianB->location_z;

    std::double_t aux = dotproduct(xA, yA, zA, xB, yB, zB);

    // first clear the product gaussians
    productGaussians->clear();

    // now loop over primitives
    for (std::uint64_t ii = 0; ii < contractedGaussianA->contracted_GTO.size(); ii++)
    {
        std::double_t expA = contractedGaussianA->contracted_GTO[ii].primitive_exp;
        for (std::uint64_t jj = 0; jj < contractedGaussianB->contracted_GTO.size(); jj++)
        {
            std::double_t expB = contractedGaussianA->contracted_GTO[jj].primitive_exp;
            cxx_Gaussians gaussian;
            
            gaussian.gaussian_center[0] = ((expA * xA) + (expB * xB)) / (expA + expB);
            gaussian.gaussian_center[1] = ((expA * yA) + (expB * yB)) / (expA + expB);
            gaussian.gaussian_center[2] = ((expA * zA) + (expB * zB)) / (expA + expB);

            gaussian.gaussian_exponent = (expA * expB) / (expA + expB);
            gaussian.gaussian_integral = exp(-1 * gaussian.gaussian_exponent * aux);

            productGaussians->push_back(gaussian);
            // std::cout << xA << " " << yA << " " << zA << " " << expA << " " << xB << " " << yB << " " << zB << " " << expB << " " << gaussian.gaussian_integral << "\n";
        }
    }
}
