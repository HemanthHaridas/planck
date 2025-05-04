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

#include "planck_hermite.h"

///
/// @brief Computes the Hermite overlap integral in one dimension.
/// This function recursively calculates the Hermite overlap integral
/// between two Gaussian primitives along a single dimension using
/// the Hermite polynomial approach.
///
/// @param exponentA Exponent of the first Gaussian function.
/// @param centerA Center coordinate of the first Gaussian function.
/// @param shellA Angular momentum quantum number of the first Gaussian function.
/// @param exponentB Exponent of the second Gaussian function.
/// @param centerB Center coordinate of the second Gaussian function.
/// @param shellB Angular momentum quantum number of the second Gaussian function.
/// @param hermiteNodes Index of the Hermite polynomial term in the recursion.
/// @return The computed Hermite overlap integral.
///

std::double_t Hermite::Overlap::computePrimitive1D(std::double_t exponentA, std::double_t centerA, std::int64_t shellA, std::double_t exponentB, std::double_t centerB, std::int64_t shellB, std::int64_t hermiteNodes)
{
    std::double_t combExp = exponentA + exponentB;
    std::double_t gaussExp = (exponentA * exponentB) / combExp;

    if ((hermiteNodes < 0) || (hermiteNodes > (shellA + shellB)))
    {
        return 0;
    }

    else if ((shellA == 0) && (shellB == 0) && (hermiteNodes == 0))
    {
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

///
/// @brief Computes the three-dimensional Hermite overlap integral.
/// This function calculates the overlap integral between two Gaussian
/// primitives in three dimensions by computing the integrals separately
/// in each coordinate direction.
///
/// @param primitiveA First Gaussian primitive.
/// @param xA X-coordinate of the first primitive.
/// @param yA Y-coordinate of the first primitive.
/// @param zA Z-coordinate of the first primitive.
/// @param lxA X-component of angular momentum for the first primitive.
/// @param lyA Y-component of angular momentum for the first primitive.
/// @param lzA Z-component of angular momentum for the first primitive.
/// @param primitiveB Second Gaussian primitive.
/// @param xB X-coordinate of the second primitive.
/// @param yB Y-coordinate of the second primitive.
/// @param zB Z-coordinate of the second primitive.
/// @param lxB X-component of angular momentum for the second primitive.
/// @param lyB Y-component of angular momentum for the second primitive.
/// @param lzB Z-component of angular momentum for the second primitive.
/// @return The computed three-dimensional overlap integral.
///

std::double_t Hermite::Overlap::computePrimitive3D(
    cxx_Primitive primitiveA, std::double_t xA, std::double_t yA, std::double_t zA, std::int64_t lxA, std::int64_t lyA, std::int64_t lzA,
    cxx_Primitive primitiveB, std::double_t xB, std::double_t yB, std::double_t zB, std::int64_t lxB, std::int64_t lyB, std::int64_t lzB)
{
    std::double_t xDir = Hermite::Overlap::computePrimitive1D(primitiveA.primitive_exp, xA, lxA, primitiveB.primitive_exp, xB, lxB, 0);
    std::double_t yDir = Hermite::Overlap::computePrimitive1D(primitiveA.primitive_exp, yA, lyA, primitiveB.primitive_exp, yB, lyB, 0);
    std::double_t zDir = Hermite::Overlap::computePrimitive1D(primitiveA.primitive_exp, zA, lzA, primitiveB.primitive_exp, zB, lzB, 0);

    std::double_t integral = xDir * yDir * zDir * pow(M_PI / (primitiveA.primitive_exp + primitiveB.primitive_exp), 1.5);

    return xDir * yDir * zDir * pow(M_PI / (primitiveA.primitive_exp + primitiveB.primitive_exp), 1.5);
}

///
/// @brief Computes the contracted overlap integral between two contracted Gaussian functions.
/// This function iterates over the primitive Gaussians in each contracted function
/// and accumulates the contributions to the overlap integral.
///
/// @param contractedGaussianA First contracted Gaussian function.
/// @param contractedGaussianB Second contracted Gaussian function.
/// @return The computed contracted overlap integral.
///

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
            value = value * contractedGaussianA.contracted_GTO[ii].orbital_coeff * contractedGaussianA.contracted_GTO[ii].orbital_norm;
            value = value * contractedGaussianB.contracted_GTO[jj].orbital_coeff * contractedGaussianB.contracted_GTO[jj].orbital_norm;
            integral = integral + value;
        }
    }
    return integral;
}

std::double_t Hermite::Kinetic::computePrimitive3D(
    cxx_Primitive primitiveA, std::double_t xA, std::double_t yA, std::double_t zA, std::int64_t lxA, std::int64_t lyA, std::int64_t lzA,
    cxx_Primitive primitiveB, std::double_t xB, std::double_t yB, std::double_t zB, std::int64_t lxB, std::int64_t lyB, std::int64_t lzB)
{
    std::double_t integral = Hermite::Overlap::computePrimitive3D(primitiveA, xA, yA, zA, lxA, lyA, lzA, primitiveB, xB, yB, zB, lxB, lyB, lzB);
    integral = integral * (primitiveB.primitive_exp * (2 * (lxB + lyB + lzB) + 3));

    integral = integral - ((2 * pow(primitiveB.primitive_exp, 2)) * (Hermite::Overlap::computePrimitive3D(primitiveA, xA, yA, zA, lxA, lyA, lzA, primitiveB, xB, yB, zB, lxB + 2, lyB, lzB)));
    integral = integral - ((2 * pow(primitiveB.primitive_exp, 2)) * (Hermite::Overlap::computePrimitive3D(primitiveA, xA, yA, zA, lxA, lyA, lzA, primitiveB, xB, yB, zB, lxB, lyB + 2, lzB)));
    integral = integral - ((2 * pow(primitiveB.primitive_exp, 2)) * (Hermite::Overlap::computePrimitive3D(primitiveA, xA, yA, zA, lxA, lyA, lzA, primitiveB, xB, yB, zB, lxB, lyB, lzB + 2)));

    integral = integral - ((0.5 * lxB * (lxB - 1)) * (Hermite::Overlap::computePrimitive3D(primitiveA, xA, yA, zA, lxA, lyA, lzA, primitiveB, xB, yB, zB, lxB - 2, lyB, lzB)));
    integral = integral - ((0.5 * lyB * (lyB - 1)) * (Hermite::Overlap::computePrimitive3D(primitiveA, xA, yA, zA, lxA, lyA, lzA, primitiveB, xB, yB, zB, lxB, lyB - 2, lzB)));
    integral = integral - ((0.5 * lzB * (lzB - 1)) * (Hermite::Overlap::computePrimitive3D(primitiveA, xA, yA, zA, lxA, lyA, lzA, primitiveB, xB, yB, zB, lxB, lyB, lzB - 2)));

    std::double_t exponentA = primitiveA.primitive_exp;
    std::double_t exponentB = primitiveB.primitive_exp;
    // std::cout << std::setw(10) << std::setprecision(3) << std::right << xA << std::setw(10) << std::setprecision(3) << std::right << yA << std::setw(10) << std::setprecision(3) << std::right << zA << std::setw(10) << std::setprecision(3) << std::right << xB << std::setw(10) << std::setprecision(3) << std::right << yB << std::setw(10) << std::setprecision(3) << std::right << zB << std::setw(20) << std::setprecision(3) << std::right << integral << std::setw(20) << std::right << lxA << std::setw(20) << std::right << lyA << std::setw(20) << std::right << lzA << std::setw(20) << std::right << lxB << std::setw(20) << std::right << lyB << std::setw(20) << std::right << lzB << "\n";
    // std::cout << std::setw(20) << std::setprecision(3) << std::right << exponentA << std::setw(20) << std::setprecision(3) << std::right << exponentB << std::setw(20) << std::setprecision(3) << std::right << integral << std::setw(20) << std::right << lxA << std::setw(20) << std::right << lxB << std::setw(20) << std::right << lyA << std::setw(20) << std::right << lyB << std::setw(20) << std::right << lzA << std::setw(20) << std::right << lzB << "\n";
    return integral;
}

std::double_t Hermite::Kinetic::computeContracted(cxx_Contracted contractedGaussianA, cxx_Contracted contractedGaussianB)
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

            std::double_t value = Hermite::Kinetic::computePrimitive3D(
                contractedGaussianA.contracted_GTO[ii], xA, yA, zA, lxA, lyA, lzA,
                contractedGaussianB.contracted_GTO[jj], xB, yB, zB, lxB, lyB, lzB);

            // contract the integrals
            value = value * contractedGaussianA.contracted_GTO[ii].orbital_coeff * contractedGaussianA.contracted_GTO[ii].orbital_norm;
            value = value * contractedGaussianB.contracted_GTO[jj].orbital_coeff * contractedGaussianB.contracted_GTO[jj].orbital_norm;
            value = value * 1;

            // collect the values;
            integral = integral + value;
        }
    }
    return integral;
}
