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
#include "planck_hermite.h"

std::double_t Hermite::Overlap::computePrimitive1D(std::double_t exponentA, std::double_t centerA, std::int64_t shellA, std::double_t exponentB, std::double_t centerB, std::int64_t shellB, std::int64_t hermiteNodes)
{
    // std::cout << exponentA << " " << centerA << " " << shellA << " " << exponentB << " " << centerB << " " << shellB << " " << hermiteNodes << "\n";
    std::double_t combExp = exponentA + exponentB;
    std::double_t gaussExp = (exponentA * exponentB) / combExp;

    if ((hermiteNodes < 0) || (hermiteNodes > (shellA + shellB)))
    {
        return 0;
    }

    else if ((shellA == 0) && (shellB == 0) && (hermiteNodes == 0))
    {
        // return 1;
        return exp(-1 * gaussExp * (centerA - centerB) * (centerA - centerB));
    }

    else if (shellB == 0)
    {
        return (
            (1 / (2 * (exponentA + exponentB))) * Hermite::Overlap::computePrimitive1D(exponentA, centerA, shellA - 1, exponentB, centerB, shellB, hermiteNodes - 1) -
            (exponentB * (centerA - centerB) / (exponentA + exponentB)) * Hermite::Overlap::computePrimitive1D(exponentA, centerA, shellA - 1, exponentB, centerB, shellB, hermiteNodes) +
            (hermiteNodes + 1) * Hermite::Overlap::computePrimitive1D(exponentA, centerA, shellA - 1, exponentB, centerB, shellB, hermiteNodes + 1));
    }

    else
    {
        return (
            (1 / (2 * (exponentA + exponentB))) * Hermite::Overlap::computePrimitive1D(exponentA, centerA, shellA, exponentB, centerB, shellB - 1, hermiteNodes - 1) +
            (exponentA * (centerA - centerB) / (exponentA + exponentB)) * Hermite::Overlap::computePrimitive1D(exponentA, centerA, shellA, exponentB, centerB, shellB - 1, hermiteNodes) +
            (hermiteNodes + 1) * Hermite::Overlap::computePrimitive1D(exponentA, centerA, shellA, exponentB, centerB, shellB - 1, hermiteNodes + 1));
    }
}

std::double_t Hermite::Overlap::computePrimitive3D(
    cxx_Primitive primitiveA, std::double_t xA, std::double_t yA, std::double_t zA, std::int64_t lxA, std::int64_t lyA, std::int64_t lzA,
    cxx_Primitive primitiveB, std::double_t xB, std::double_t yB, std::double_t zB, std::int64_t lxB, std::int64_t lyB, std::int64_t lzB)
{
    std::double_t xDir = Hermite::Overlap::computePrimitive1D(primitiveA.primitive_exp, xA, lxA, primitiveB.primitive_exp, xB, lxB, 0);
    std::double_t yDir = Hermite::Overlap::computePrimitive1D(primitiveA.primitive_exp, yA, lyA, primitiveB.primitive_exp, yB, lyB, 0);
    std::double_t zDir = Hermite::Overlap::computePrimitive1D(primitiveA.primitive_exp, zA, lzA, primitiveB.primitive_exp, zB, lzB, 0);

    xDir = xDir * 1;
    yDir = yDir * 1;
    zDir = zDir * primitiveA.orbital_coeff * primitiveA.orbital_norm * primitiveB.orbital_coeff * primitiveB.orbital_norm;
    return xDir * yDir * zDir * pow(M_PI / (primitiveA.primitive_exp + primitiveB.primitive_exp), 1.5);
}

std::double_t Hermite::Overlap::computeContracted(cxx_Contracted contractedGaussianA, cxx_Contracted contractedGaussianB)
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

            std::double_t value = Hermite::Overlap::computePrimitive3D(
                contractedGaussianA.contracted_GTO[ii], xA, yA, zA, lxA, lyA, lzA,
                contractedGaussianB.contracted_GTO[jj], xB, yB, zB, lxB, lyB, lzB);

            // contract the integrals
            // value = value * contractedGaussianA.contracted_GTO[ii].orbital_coeff * contractedGaussianA.contracted_GTO[ii].orbital_norm;
            // value = value * contractedGaussianB.contracted_GTO[jj].orbital_coeff * contractedGaussianB.contracted_GTO[jj].orbital_norm;
            // std::cout << contractedGaussianA.contracted_GTO[ii].orbital_coeff * contractedGaussianA.contracted_GTO[ii].orbital_norm << " " << contractedGaussianB.contracted_GTO[jj].orbital_coeff * contractedGaussianB.contracted_GTO[jj].orbital_norm << "\n";
            // collect the values;
            integral = integral + (value * 1);
        }
        // std::cout << "\n";
    }
    std::cout << integral << "\n";
    return integral;
}