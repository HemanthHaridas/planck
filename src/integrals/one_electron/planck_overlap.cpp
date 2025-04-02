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
#include "planck_overlap.h"
#include "../helper/planck_helper_routines.h"

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
            planckCalculator->overlap[ii * planckCalculator->total_basis + jj] = computeOverlap1(&planckCalculator->calculation_set[ii], &planckCalculator->calculation_set[jj], errorFlag, errorMessage);
        }
    }
}

std::double_t computeOverlap1(cxx_Contracted *contractedGaussianA, cxx_Contracted *contractedGaussianB, std::error_code *errorFlag, std::string *errorMessage)
{
    std::vector<cxx_Gaussians> productGaussians;
    std::double_t xA = contractedGaussianA->location_x;
    std::double_t yA = contractedGaussianA->location_y;
    std::double_t zA = contractedGaussianA->location_z;

    std::int64_t lxA = contractedGaussianA->shell_x;
    std::int64_t lyA = contractedGaussianA->shell_y;
    std::int64_t lzA = contractedGaussianA->shell_z;

    std::double_t xB = contractedGaussianB->location_x;
    std::double_t yB = contractedGaussianB->location_y;
    std::double_t zB = contractedGaussianB->location_z;

    std::int64_t lxB = contractedGaussianB->shell_x;
    std::int64_t lyB = contractedGaussianB->shell_y;
    std::int64_t lzB = contractedGaussianB->shell_z;

    // containers
    std::double_t locA[] = {xA, yA, zA};
    std::double_t locB[] = {xB, yB, zB};
    std::int64_t shellA[] = {lxA, lyA, lzA};
    std::int64_t shellB[] = {lxB, lyB, lzB};
    std::double_t primtiveOverlaps[4]; // primitive overlap in 3 directions and the contracted value

    std::double_t primtiveOverlap = 0.0;
    computeGaussianProduct(contractedGaussianA, contractedGaussianB, &productGaussians, errorFlag, errorMessage);

    for (std::uint64_t ii = 0; ii < contractedGaussianA->contracted_GTO.size(); ii++)
    {
        for (std::uint64_t jj = 0; jj < contractedGaussianB->contracted_GTO.size(); jj++)
        {
            // compute the integrals over primitives
            computePrimitive_Huzinaga(
                &contractedGaussianA->contracted_GTO[ii], locA, shellA, 
                &contractedGaussianB->contracted_GTO[jj], locB, shellB, 
                &productGaussians[ii * contractedGaussianB->contracted_GTO.size() + jj], primtiveOverlaps
                );
            
            // combine the primtive overlaps
            primtiveOverlap = primtiveOverlap + primtiveOverlaps[3];
            // std::cout << primtiveOverlaps[0] << " " << primtiveOverlaps[1] << " " << primtiveOverlaps[2] << " " << primtiveOverlaps[3] << "\n";
        }
    }
    return primtiveOverlap;
}

// based on the reference implementation in dx.doi.org/doi:10.3888/tmj.14-3
// Evaluation of Gaussian Molecular Integrals I. Overlap Integrals
// original reference in 
//
void computePrimitive_Huzinaga(cxx_Primitive *primitiveA, std::double_t *locA, std::int64_t *shellA, cxx_Primitive *primitiveB, std::double_t *locB, std::int64_t *shellB, cxx_Gaussians *productGaussian, std::double_t *primitiveOverlaps)
{
    std::double_t aux = 0.0;
    memset(primitiveOverlaps, 0, 4 * sizeof(std::double_t));

    // first do x-diection
    for (std::uint64_t ii = 0; ii <= shellA[0]; ii++)
    {
        for (std::uint64_t jj = 0; jj <= shellB[0]; jj++)
        {
            if ((ii + jj) % 2 != 0)
            {
                continue;
            }

            aux = combination(shellA[0], ii) * combination(shellB[0], jj) * doublefactorial(ii + jj - 1);
            aux = aux * pow(productGaussian->gaussian_center[0] - locA[0], shellA[0] - ii) * pow(productGaussian->gaussian_center[0] - locB[0], shellB[0] - jj);
            aux = aux / pow(2 * (primitiveA->primitive_exp + primitiveB->primitive_exp), 0.5 * (ii + jj));
            primitiveOverlaps[0] = primitiveOverlaps[0] + aux;
            // std::cout << ii << " " << jj << " " << primitiveOverlaps[0] << "\n";
        }
    }
    primitiveOverlaps[0] = primitiveOverlaps[0] * pow(M_PI / (primitiveA->primitive_exp + primitiveB->primitive_exp), 0.5);

    // then do y-diection
    for (std::uint64_t ii = 0; ii <= shellA[1]; ii++)
    {
        for (std::uint64_t jj = 0; jj <= shellB[1]; jj++)
        {
            if ((ii + jj) % 2 != 0)
            {
                continue;
            }

            aux = combination(shellA[1], ii) * combination(shellB[1], jj) * doublefactorial(ii + jj - 1);
            aux = aux * pow(productGaussian->gaussian_center[1] - locA[1], shellA[1] - ii) * pow(productGaussian->gaussian_center[1] - locB[1], shellB[0] - jj);
            aux = aux / pow(2 * (primitiveA->primitive_exp + primitiveB->primitive_exp), 0.5 * (ii + jj));
            primitiveOverlaps[1] = primitiveOverlaps[1] + aux;
            // std::cout << ii << " " << jj << " " << primitiveOverlaps[1] << "\n";
        }
    }
    primitiveOverlaps[1] = primitiveOverlaps[1] * pow(M_PI / (primitiveA->primitive_exp + primitiveB->primitive_exp), 0.5);

    // then do z-diection
    for (std::uint64_t ii = 0; ii <= shellA[2]; ii++)
    {
        for (std::uint64_t jj = 0; jj <= shellB[2]; jj++)
        {
            if ((ii + jj) % 2 != 0)
            {
                continue;
            }

            aux = combination(shellA[2], ii) * combination(shellB[2], jj) * doublefactorial(ii + jj - 1);
            aux = aux * pow(productGaussian->gaussian_center[2] - locA[2], shellA[2] - ii) * pow(productGaussian->gaussian_center[2] - locB[2], shellB[2] - jj);
            aux = aux / pow(2 * (primitiveA->primitive_exp + primitiveB->primitive_exp), 0.5 * (ii + jj));
            primitiveOverlaps[2] = primitiveOverlaps[2] + aux;
            // std::cout << ii << " " << jj << " " << primitiveOverlaps[2] << "\n";
        }
    }
    primitiveOverlaps[2] = primitiveOverlaps[2] * pow(M_PI / (primitiveA->primitive_exp + primitiveB->primitive_exp), 0.5);

    // now contract the integral
    primitiveOverlaps[3] = primitiveOverlaps[0] * primitiveOverlaps[1] * primitiveOverlaps[2] * productGaussian->gaussian_integral[4];
    primitiveOverlaps[3] = primitiveOverlaps[3] * primitiveA->orbital_coeff * primitiveA->orbital_norm;
    primitiveOverlaps[3] = primitiveOverlaps[3] * primitiveB->orbital_coeff * primitiveB->orbital_norm;
}

void computePrimitive_ObaraSakia(cxx_Primitive *primitiveA, std::double_t *locA, std::int64_t *shellA, cxx_Primitive *primitiveB, std::double_t *locB, std::int64_t *shellB, cxx_Gaussians *productGaussian, std::double_t *primitiveOverlaps)
{
    // precompute eta and distances
    std::double_t eta = primitiveA->primitive_exp + primitiveB->primitive_exp;
    std::double_t xGA = productGaussian->gaussian_center[0] - locA[0];
    std::double_t yGA = productGaussian->gaussian_center[1] - locA[1];
    std::double_t zGA = productGaussian->gaussian_center[2] - locA[2];

    std::double_t xAB = locB[0] - locA[0];
    std::double_t yAB = locB[1] - locA[1];
    std::double_t zAB = locB[2] - locA[2];
}