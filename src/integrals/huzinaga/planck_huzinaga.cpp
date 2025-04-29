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

#include "planck_huzinaga.h"

// @brief Computes the one-dimensional primitive integral between two Gaussians.
//
// This function calculates the 1D integral contribution of two Gaussians
// along a specific axis, using their exponents, centers, shell information,
// and the Gaussian product center.
//
// @param exponentA The exponent of the first Gaussian.
// @param centerA The center coordinate of the first Gaussian.
// @param shellA The angular momentum quantum number of the first Gaussian.
// @param exponentB The exponent of the second Gaussian.
// @param centerB The center coordinate of the second Gaussian.
// @param shellB The angular momentum quantum number of the second Gaussian.
// @param gaussianCenter The center of the Gaussian product.
// @return The computed one-dimensional primitive integral value.

std::double_t Huzinaga::Overlap::computePrimitive1D(std::double_t exponentA, std::double_t centerA, std::int64_t shellA, std::double_t exponentB, std::double_t centerB, std::int64_t shellB, std::double_t gaussianCenter)
{
    std::double_t integral = 0.0;
    for (std::int64_t ii = 0; ii < shellA + 1; ii++)
    {
        for (std::int64_t jj = 0; jj < shellB + 1; jj++)
        {
            if ((ii + jj) % 2 == 0)
            {
                std::double_t value = combination(shellA, ii);
                value = value * combination(shellB, jj);
                value = value * doublefactorial(ii + jj - 1);
                value = value * pow(gaussianCenter - centerA, shellA - ii);
                value = value * pow(gaussianCenter - centerB, shellB - jj);
                value = value / pow(2 * (exponentA + exponentB), 0.5 * (ii + jj));
                integral = integral + value;
            }
        }
    }
    return integral;
}

// @brief Computes the primitive integral between two contracted Gaussian-type orbitals.
//
// This function calculates the integral of two contracted Gaussian-type orbitals (GTOs)
// using their primitive components, location, shell information, and normalization factors.
//
// @param contractedGaussianA The first contracted Gaussian, containing its location,
//        shell information, and a list of primitive Gaussians.
// @param contractedGaussianB The second contracted Gaussian, containing its location,
//        shell information, and a list of primitive Gaussians.
// @return The computed value of the integral.

std::double_t Huzinaga::Overlap::computeContracted(cxx_Contracted contractedGaussianA, cxx_Contracted contractedGaussianB)
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
            std::double_t value = Huzinaga::Overlap::computePrimitive3D(
                contractedGaussianA.contracted_GTO[ii], xA, yA, zA, lxA, lyA, lzA,
                contractedGaussianB.contracted_GTO[jj], xB, yB, zB, lxB, lyB, lzB);

            // collect the values;
            integral = integral + value;
        }
    }
    return integral;
}

// @brief Computes the three-dimensional primitive integral between two Gaussian primitives.
//
// This function calculates the 3D integral of two Gaussian primitives, considering their exponents,
// centers, angular momentum quantum numbers, and normalization coefficients.
//
// @param primitiveA The first Gaussian primitive containing its exponent, orbital coefficient, and normalization factor.
// @param xA The x-coordinate of the center of the first primitive.
// @param yA The y-coordinate of the center of the first primitive.
// @param zA The z-coordinate of the center of the first primitive.
// @param lxA The angular momentum quantum number along the x-axis for the first primitive.
// @param lyA The angular momentum quantum number along the y-axis for the first primitive.
// @param lzA The angular momentum quantum number along the z-axis for the first primitive.
// @param primitiveB The second Gaussian primitive containing its exponent, orbital coefficient, and normalization factor.
// @param xB The x-coordinate of the center of the second primitive.
// @param yB The y-coordinate of the center of the second primitive.
// @param zB The z-coordinate of the center of the second primitive.
// @param lxB The angular momentum quantum number along the x-axis for the second primitive.
// @param lyB The angular momentum quantum number along the y-axis for the second primitive.
// @param lzB The angular momentum quantum number along the z-axis for the second primitive.
// @return The computed three-dimensional primitive integral value.

std::double_t Huzinaga::Overlap::computePrimitive3D(
    cxx_Primitive primitiveA, std::double_t xA, std::double_t yA, std::double_t zA, std::int64_t lxA, std::int64_t lyA, std::int64_t lzA,
    cxx_Primitive primitiveB, std::double_t xB, std::double_t yB, std::double_t zB, std::int64_t lxB, std::int64_t lyB, std::int64_t lzB)
{
    cxx_Gaussians productGaussianAB = computeGaussianProduct(primitiveA, xA, yA, zA, primitiveB, xB, yB, zB);

    std::double_t exponentA = primitiveA.primitive_exp;
    std::double_t exponentB = primitiveB.primitive_exp;

    std::double_t xValue = Huzinaga::Overlap::computePrimitive1D(exponentA, xA, lxA, exponentB, xB, lxB, productGaussianAB.gaussian_center[0]);
    std::double_t yValue = Huzinaga::Overlap::computePrimitive1D(exponentA, yA, lyA, exponentB, yB, lyB, productGaussianAB.gaussian_center[1]);
    std::double_t zValue = Huzinaga::Overlap::computePrimitive1D(exponentA, zA, lzA, exponentB, zB, lzB, productGaussianAB.gaussian_center[2]);

    std::double_t integral3D = xValue * yValue * zValue * productGaussianAB.gaussian_integral[3] * pow(M_PI / (exponentA + exponentB), 1.5);
    integral3D = integral3D * primitiveA.orbital_coeff * primitiveA.orbital_norm;
    integral3D = integral3D * primitiveB.orbital_coeff * primitiveB.orbital_norm;

    return integral3D;
}

std::double_t Huzinaga::Kinetic::computePrimitive3D(
    cxx_Primitive primitiveA, std::double_t xA, std::double_t yA, std::double_t zA, std::int64_t lxA, std::int64_t lyA, std::int64_t lzA,
    cxx_Primitive primitiveB, std::double_t xB, std::double_t yB, std::double_t zB, std::int64_t lxB, std::int64_t lyB, std::int64_t lzB)
{
    cxx_Gaussians productGaussianAB = computeGaussianProduct(primitiveA, xA, yA, zA, primitiveB, xB, yB, zB);

    std::double_t exponentA = primitiveA.primitive_exp;
    std::double_t exponentB = primitiveB.primitive_exp;

    // first compute regular overlaps
    std::double_t xOriginal = Huzinaga::Overlap::computePrimitive1D(exponentA, xA, lxA, exponentB, xB, lxB, productGaussianAB.gaussian_center[0]);
    std::double_t yOriginal = Huzinaga::Overlap::computePrimitive1D(exponentA, yA, lyA, exponentB, yB, lyB, productGaussianAB.gaussian_center[1]);
    std::double_t zOriginal = Huzinaga::Overlap::computePrimitive1D(exponentA, zA, lzA, exponentB, zB, lzB, productGaussianAB.gaussian_center[2]);

    // now compute the (-,-,-,-) overlaps
    std::double_t xNN = Huzinaga::Overlap::computePrimitive1D(exponentA, xA, lxA - 1, exponentB, xB, lxB - 1, productGaussianAB.gaussian_center[0]);
    std::double_t yNN = Huzinaga::Overlap::computePrimitive1D(exponentA, yA, lyA - 1, exponentB, yB, lyB - 1, productGaussianAB.gaussian_center[1]);
    std::double_t zNN = Huzinaga::Overlap::computePrimitive1D(exponentA, zA, lzA - 1, exponentB, zB, lzB - 1, productGaussianAB.gaussian_center[2]);

    // now compute the (+,+,+,+) overlaps
    std::double_t xPP = Huzinaga::Overlap::computePrimitive1D(exponentA, xA, lxA + 1, exponentB, xB, lxB + 1, productGaussianAB.gaussian_center[0]);
    std::double_t yPP = Huzinaga::Overlap::computePrimitive1D(exponentA, yA, lyA + 1, exponentB, yB, lyB + 1, productGaussianAB.gaussian_center[1]);
    std::double_t zPP = Huzinaga::Overlap::computePrimitive1D(exponentA, zA, lzA + 1, exponentB, zB, lzB + 1, productGaussianAB.gaussian_center[2]);

    // now compute the (-,+,-,+) overlaps
    std::double_t xNP = Huzinaga::Overlap::computePrimitive1D(exponentA, xA, lxA - 1, exponentB, xB, lxB + 1, productGaussianAB.gaussian_center[0]);
    std::double_t yNP = Huzinaga::Overlap::computePrimitive1D(exponentA, yA, lyA - 1, exponentB, yB, lyB + 1, productGaussianAB.gaussian_center[1]);
    std::double_t zNP = Huzinaga::Overlap::computePrimitive1D(exponentA, zA, lzA - 1, exponentB, zB, lzB + 1, productGaussianAB.gaussian_center[2]);

    // now compute the (+,-,+,-) overlaps
    std::double_t xPN = Huzinaga::Overlap::computePrimitive1D(exponentA, xA, lxA + 1, exponentB, xB, lxB - 1, productGaussianAB.gaussian_center[0]);
    std::double_t yPN = Huzinaga::Overlap::computePrimitive1D(exponentA, yA, lyA + 1, exponentB, yB, lyB - 1, productGaussianAB.gaussian_center[1]);
    std::double_t zPN = Huzinaga::Overlap::computePrimitive1D(exponentA, zA, lzA + 1, exponentB, zB, lzB - 1, productGaussianAB.gaussian_center[2]);

    std::double_t tX = (lxA * lxB * xNN) + (-2 * exponentA * lxB * xPN) + (-2 * exponentB * lxA * xNP) + (4 * exponentA * exponentB * xPP);
    std::double_t tY = (lyA * lyB * yNN) + (-2 * exponentA * lyB * yPN) + (-2 * exponentB * lyA * yNP) + (4 * exponentA * exponentB * yPP);
    std::double_t tZ = (lzA * lzB * zNN) + (-2 * exponentA * lzB * zPN) + (-2 * exponentB * lzA * zNP) + (4 * exponentA * exponentB * zPP);

    std::double_t kX = 0.5 * tX * yOriginal * zOriginal * productGaussianAB.gaussian_integral[3] * pow(M_PI / (exponentA + exponentB), 1.5);
    std::double_t kY = 0.5 * tY * zOriginal * xOriginal * productGaussianAB.gaussian_integral[3] * pow(M_PI / (exponentA + exponentB), 1.5);
    std::double_t kZ = 0.5 * tZ * xOriginal * yOriginal * productGaussianAB.gaussian_integral[3] * pow(M_PI / (exponentA + exponentB), 1.5);

    kX = kX * primitiveA.orbital_coeff * primitiveA.orbital_norm;
    kY = kY * primitiveA.orbital_coeff * primitiveA.orbital_norm;
    kZ = kZ * primitiveA.orbital_coeff * primitiveA.orbital_norm;
}