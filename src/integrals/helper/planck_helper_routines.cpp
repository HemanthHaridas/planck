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
#include <boost/math/special_functions/hypergeometric_1F1.hpp>

#include "planck_helper_routines.h"

void computeGaussianProduct(cxx_Contracted *contractedGaussianA, cxx_Contracted *contractedGaussianB, std::vector<cxx_Gaussians> *productGaussians)
{
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
            std::double_t expB = contractedGaussianB->contracted_GTO[jj].primitive_exp;
            cxx_Gaussians gaussian;

            gaussian.gaussian_center[0] = ((expA * xA) + (expB * xB)) / (expA + expB);
            gaussian.gaussian_center[1] = ((expA * yA) + (expB * yB)) / (expA + expB);
            gaussian.gaussian_center[2] = ((expA * zA) + (expB * zB)) / (expA + expB);

            gaussian.gaussian_exponent = (expA * expB) / (expA + expB);
            gaussian.gaussian_integral[0] = exp(-1 * gaussian.gaussian_exponent * (xA - xB) * (xA - xB));
            gaussian.gaussian_integral[1] = exp(-1 * gaussian.gaussian_exponent * (yA - yB) * (yA - yB));
            gaussian.gaussian_integral[2] = exp(-1 * gaussian.gaussian_exponent * (zA - zB) * (zA - zB));
            gaussian.gaussian_integral[3] = exp(-1 * gaussian.gaussian_exponent * aux);

            productGaussians->push_back(gaussian);
            // std::cout << xA << " " << yA << " " << zA << " " << expA << " " << xB << " " << yB << " " << zB << " " << expB << " " << gaussian.gaussian_integral << "\n";
        }
    }
}

cxx_Gaussians computeGaussianProduct(cxx_Primitive *primitiveA, const std::double_t xA, const std::double_t yA, const std::double_t zA, cxx_Primitive *primitiveB, const std::double_t xB, const std::double_t yB, const std::double_t zB)
{
    cxx_Gaussians gaussian;
    std::double_t aux = dotproduct(xA, yA, zA, xB, yB, zB);
    gaussian.gaussian_center[0] = ((primitiveA->primitive_exp * xA) + (primitiveB->primitive_exp * xB)) / (primitiveA->primitive_exp + primitiveB->primitive_exp);
    gaussian.gaussian_center[1] = ((primitiveA->primitive_exp * yA) + (primitiveB->primitive_exp * yB)) / (primitiveA->primitive_exp + primitiveB->primitive_exp);
    gaussian.gaussian_center[2] = ((primitiveA->primitive_exp * zA) + (primitiveB->primitive_exp * zB)) / (primitiveA->primitive_exp + primitiveB->primitive_exp);

    gaussian.gaussian_exponent = (primitiveA->primitive_exp * primitiveB->primitive_exp) / (primitiveA->primitive_exp + primitiveB->primitive_exp);

    gaussian.gaussian_integral[0] = exp(-1 * gaussian.gaussian_exponent * (xA - xB) * (xA - xB));
    gaussian.gaussian_integral[1] = exp(-1 * gaussian.gaussian_exponent * (yA - yB) * (yA - yB));
    gaussian.gaussian_integral[2] = exp(-1 * gaussian.gaussian_exponent * (zA - zB) * (zA - zB));
    gaussian.gaussian_integral[3] = exp(-1 * gaussian.gaussian_exponent * aux);

    return gaussian;
}

std::double_t boysFunction(std::uint64_t boysIndex, std::double_t boysParam)
{
    // for some reason the indices in the boys lookup table are swapped
    std::uint64_t yIndex = boysIndex;
    std::uint64_t xIndex = boysParam / 0.1;

    // if boysParam in out of bounds
    // compute it using Kummer's hyperconfluent function
    if (xIndex >= MAXM)
    {
        return boost::math::hypergeometric_1F1(0.5 + boysIndex, 1.5 + boysIndex, -1 * boysParam) / (2 * boysIndex + 1);
    }

    // now compute the integral for boys function as a six term taylor series
    std::double_t boysValue = boysTable[xIndex][yIndex];
    std::double_t delta = boysParam - (0.1 * xIndex);

    for (std::uint64_t ii = 1; ii <= 6; ii++)
    {
        boysValue = boysValue + (pow(-1, ii) * (boysTable[xIndex][yIndex + ii] / factorial(ii)) * pow(delta, ii));
    }
    return boysValue;
}

std::vector<eriKet> schwatrzSceening(cxx_Calculator *planckCalculator, Eigen::Tensor<std::double_t, 4> &electronicMatrix)
{
    std::int64_t nBasis = planckCalculator->total_basis;
    // electronicMatrix.resize(nBasis, nBasis, nBasis, nBasis);

    // calculate (ab|ab) and store
    for (std::int64_t row = 0; row < nBasis; row++)
    {
        for (std::int64_t col = row; col < nBasis; col++)
        {
            // has 8-fold permutation symmetry
            std::double_t integral = Huzinaga::computeElectronic(
                &planckCalculator->calculation_set[row], &planckCalculator->calculation_set[col],
                &planckCalculator->calculation_set[row], &planckCalculator->calculation_set[col]);

            // Store integral and apply 8-fold symmetry
            electronicMatrix(row, col, row, col) = integral;
            electronicMatrix(row, row, col, col) = integral;
            electronicMatrix(row, col, col, row) = integral;
            electronicMatrix(row, col, row, col) = integral;

            electronicMatrix(col, row, row, col) = integral;
            electronicMatrix(col, col, row, row) = integral;
            electronicMatrix(col, row, col, row) = integral;
            electronicMatrix(col, row, row, col) = integral;
        }
    }

    // find max of (ab|ab) and store it
    // now multiply the integrals with the maximum value
    Eigen::Tensor<std::double_t, 0> schwartzMax = electronicMatrix.maximum();
    Eigen::Tensor<std::double_t, 0> schwartzMin = electronicMatrix.minimum();
    std::double_t schwartzValue = abs(schwartzMax(0)) > abs(schwartzMin(0)) ? abs(schwartzMax(0)) : abs(schwartzMin(0));

    std::vector<eriKet> eriScreened;

    // now screen the integrals based on cutoff
    for (std::int64_t row = 0; row < nBasis; row++)
    {
        for (std::int64_t col = row; col < nBasis; col++)
        {
            eriKet ket(row, col); // create a tuple of indices
            std::double_t value = abs(electronicMatrix(row, col, row, col));
            std::double_t eval = sqrt(schwartzValue) * sqrt(value);

            if (eval >= planckCalculator->tol_eri)
            {
                eriScreened.push_back(ket);
            }
        }
    }
    return eriScreened;
}