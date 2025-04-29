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
#include "planck_obarasakia.h"

std::double_t ObaraSakia::Overlap::computePrimitive1D(std::double_t centerA, std::double_t exponentA, std::int64_t shellA, std::double_t centerB, std::double_t exponentB, std::int64_t shellB, std::double_t gaussianCenter, std::double_t gaussianExponent)
{
    if ((shellA < 0) || (shellB < 0))
    {
        return 0.0;
    }

    else if ((shellA == 0) && (shellB == 0))
    {
        // std::cout << shellA << " " << shellB << " " << exponentA << " " << exponentB << " " << gaussianExponent << " " << exp(-1 * gaussianExponent * (centerA - centerB) * (centerA - centerB)) << "\n";
        return exp(-1 * gaussianExponent * (centerA - centerB) * (centerA - centerB));
    }

    else if (shellB == 0 && shellA > 0)
    {
        return (
            (gaussianCenter - centerA) * ObaraSakia::Overlap::computePrimitive1D(centerA, exponentA, shellA, centerB, exponentB, shellB - 1, gaussianCenter, gaussianExponent) +
            (1 / (2 * (exponentA + exponentB))) * ((shellA - 1) * ObaraSakia::Overlap::computePrimitive1D(centerA, exponentA, shellA - 2, centerB, exponentB, shellB, gaussianCenter, gaussianExponent) +
                                                   (shellB + 0) * ObaraSakia::Overlap::computePrimitive1D(centerA, exponentA, shellA - 1, centerB, exponentB, shellB - 1, gaussianCenter, gaussianExponent)));
    }

    else
    {
        return (
            (gaussianCenter - centerB) * ObaraSakia::Overlap::computePrimitive1D(centerA, exponentA, shellA - 1, centerB, exponentB, shellB, gaussianCenter, gaussianExponent) +
            (1 / (2 * (exponentA + exponentB))) * ((shellA + 0) * ObaraSakia::Overlap::computePrimitive1D(centerA, exponentA, shellA - 1, centerB, exponentB, shellB - 1, gaussianCenter, gaussianExponent) +
                                                   (shellB - 1) * ObaraSakia::Overlap::computePrimitive1D(centerA, exponentA, shellA, centerB, exponentB, shellB - 2, gaussianCenter, gaussianExponent)));
    }
}

std::double_t ObaraSakia::Overlap::computePrimitive3D(
    cxx_Primitive primitiveA, std::double_t xA, std::double_t yA, std::double_t zA, std::int64_t lxA, std::int64_t lyA, std::int64_t lzA,
    cxx_Primitive primitiveB, std::double_t xB, std::double_t yB, std::double_t zB, std::int64_t lxB, std::int64_t lyB, std::int64_t lzB,
    std::double_t gaussianCenterX, std::double_t gaussianCenterY, std::double_t gaussianCenterZ, std::double_t gaussianExponent)
{
    std::double_t xDir = ObaraSakia::Overlap::computePrimitive1D(xA, primitiveA.primitive_exp, lxA, xB, primitiveB.primitive_exp, lxB, gaussianCenterX, gaussianExponent);
    std::double_t yDir = ObaraSakia::Overlap::computePrimitive1D(yA, primitiveA.primitive_exp, lyA, yB, primitiveB.primitive_exp, lyB, gaussianCenterY, gaussianExponent);
    std::double_t zDir = ObaraSakia::Overlap::computePrimitive1D(zA, primitiveA.primitive_exp, lzA, zB, primitiveB.primitive_exp, lzB, gaussianCenterZ, gaussianExponent);

    xDir = xDir * 1;
    yDir = yDir * 1;
    zDir = zDir * primitiveA.orbital_coeff * primitiveA.orbital_norm * primitiveB.orbital_coeff * primitiveB.orbital_norm;
    // std::cout << lxA << " " << lxB << " " << xDir << " " << lyA << " " << lyB << " " << " " << yDir << " " << zDir << "\n";
    return xDir * yDir * zDir * pow(M_PI / (primitiveA.primitive_exp + primitiveB.primitive_exp), 1.5);
}

std::double_t ObaraSakia::Overlap::computeContracted(cxx_Contracted contractedGaussianA, cxx_Contracted contractedGaussianB)
{
    // data for first shell
    std::double_t xA = contractedGaussianA.location_x;
    std::double_t yA = contractedGaussianA.location_y;
    std::double_t zA = contractedGaussianA.location_z;

    std::int64_t lxA = contractedGaussianA.shell_x;
    std::int64_t lyA = contractedGaussianA.shell_y;
    std::int64_t lzA = contractedGaussianA.shell_z;

    std::uint64_t nA = contractedGaussianA.contracted_GTO.size();

    // data for second shell
    std::double_t xB = contractedGaussianB.location_x;
    std::double_t yB = contractedGaussianB.location_y;
    std::double_t zB = contractedGaussianB.location_z;

    std::int64_t lxB = contractedGaussianB.shell_x;
    std::int64_t lyB = contractedGaussianB.shell_y;
    std::int64_t lzB = contractedGaussianB.shell_z;

    std::uint64_t nB = contractedGaussianB.contracted_GTO.size();

    std::double_t integral = 0.0;

    for (std::uint64_t ii = 0; ii < nA; ii++)
    {
        for (std::uint64_t jj = 0; jj < nB; jj++)
        {
            cxx_Gaussians productGaussianAB = computeGaussianProduct(contractedGaussianA.contracted_GTO[ii], xA, yA, zA, contractedGaussianB.contracted_GTO[jj], xB, yB, zB);

            std::double_t value = ObaraSakia::Overlap::computePrimitive3D(
                contractedGaussianA.contracted_GTO[ii], xA, yA, zA, lxA, lyA, lzA,
                contractedGaussianB.contracted_GTO[jj], xB, yB, zB, lxB, lyB, lzB,
                productGaussianAB.gaussian_center[0], productGaussianAB.gaussian_center[1], productGaussianAB.gaussian_center[2], productGaussianAB.gaussian_exponent);

            // contract the integrals
            // value = value * contractedGaussianA.contracted_GTO[ii].orbital_coeff * contractedGaussianA.contracted_GTO[ii].orbital_norm;
            // value = value * contractedGaussianB.contracted_GTO[jj].orbital_coeff * contractedGaussianB.contracted_GTO[jj].orbital_norm;
            // std::cout << contractedGaussianA.contracted_GTO[ii].orbital_coeff * contractedGaussianA.contracted_GTO[ii].orbital_norm << " " << contractedGaussianB.contracted_GTO[jj].orbital_coeff * contractedGaussianB.contracted_GTO[jj].orbital_norm << "\n";
            // collect the values;
            integral = integral + (value * 1);
        }
        // std::cout << "\n";
    }
    return integral;
}