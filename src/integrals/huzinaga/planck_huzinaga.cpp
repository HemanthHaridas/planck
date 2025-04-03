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
#include "../helper/planck_helper_routines.h"

std::double_t Huzinaga::computeOverlap(cxx_Contracted *contractedGaussianA, cxx_Contracted *contractedGaussianB, std::error_code *errorFlag, std::string *errorMessage)
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

    std::double_t primtiveOverlaps[4]; // primitive overlap in 3 directions and the contracted value
    std::double_t primtiveOverlap = 0.0;

    // compute gaussian products
    computeGaussianProduct(contractedGaussianA, contractedGaussianB, &productGaussians, errorFlag, errorMessage);

    for (std::uint64_t ii = 0; ii < contractedGaussianA->contracted_GTO.size(); ii++)
    {
        for (std::uint64_t jj = 0; jj < contractedGaussianB->contracted_GTO.size(); jj++)
        {
            // compute the integrals over primitives
            Huzinaga::computePrimitive(
                &contractedGaussianA->contracted_GTO[ii], xA, yA, zA, lxA, lyA, lzA,
                &contractedGaussianB->contracted_GTO[jj], xB, yB, zB, lxB, lyB, lzB,
                &productGaussians[ii * contractedGaussianB->contracted_GTO.size() + jj], primtiveOverlaps);

            // combine the primtive overlaps
            primtiveOverlap = primtiveOverlap + primtiveOverlaps[3];
            // std::cout << primtiveOverlaps[0] << " " << primtiveOverlaps[1] << " " << primtiveOverlaps[2] << " " << primtiveOverlaps[3] << "\n";
        }
    }
    return primtiveOverlap;
}

std::double_t Huzinaga::computeKinetic(cxx_Contracted *contractedGaussianA, cxx_Contracted *contractedGaussianB, std::error_code *errorFlag, std::string *errorMessage)
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

    // auxiliary variables
    std::double_t tx = 0, ty = 0, tz = 0;

    std::double_t primtiveOverlaps[5][4]; // primitive overlap in 3 directions and the contracted value
    std::double_t primitiveKinetic[3];

    memset(primitiveKinetic, 0, 3 * sizeof(std::double_t));

    // compute the gaussian products
    // should merge with the overlap integral computation
    // now the overlap integral is recomputed
    computeGaussianProduct(contractedGaussianA, contractedGaussianB, &productGaussians, errorFlag, errorMessage);

    for (std::uint64_t ii = 0; ii < contractedGaussianA->contracted_GTO.size(); ii++)
    {
        std::double_t coeffA = contractedGaussianA->contracted_GTO[ii].orbital_coeff;
        std::double_t expA = contractedGaussianA->contracted_GTO[ii].primitive_exp;
        std::double_t normA = contractedGaussianA->contracted_GTO[ii].orbital_norm;
        for (std::uint64_t jj = 0; jj < contractedGaussianB->contracted_GTO.size(); jj++)
        {
            cxx_Gaussians gaussianProduct = productGaussians[ii * contractedGaussianB->contracted_GTO.size() + jj];
            std::double_t coeffB = contractedGaussianB->contracted_GTO[jj].orbital_coeff;
            std::double_t expB = contractedGaussianB->contracted_GTO[jj].primitive_exp;
            std::double_t normB = contractedGaussianB->contracted_GTO[jj].orbital_norm;

            // compute the integrals over primitives for the derivates
            Huzinaga::computePrimitive(
                &contractedGaussianA->contracted_GTO[ii], xA, yA, zA, lxA, lyA, lzA,
                &contractedGaussianB->contracted_GTO[jj], xB, yB, zB, lxB, lyB, lzB,
                &productGaussians[ii * contractedGaussianB->contracted_GTO.size() + jj], primtiveOverlaps[0]);

            // first do (-,-,-,-,-,-)
            Huzinaga::computePrimitive(
                &contractedGaussianA->contracted_GTO[ii], xA, yA, zA, lxA - 1, lyA - 1, lzA - 1,
                &contractedGaussianB->contracted_GTO[jj], xB, yB, zB, lxB - 1, lyB - 1, lzB - 1,
                &productGaussians[ii * contractedGaussianB->contracted_GTO.size() + jj], primtiveOverlaps[1]);

            // then do (+,+,+,-,-,-)
            Huzinaga::computePrimitive(
                &contractedGaussianA->contracted_GTO[ii], xA, yA, zA, lxA + 1, lyA + 1, lzA + 1,
                &contractedGaussianB->contracted_GTO[jj], xB, yB, zB, lxB - 1, lyB - 1, lzB - 1,
                &productGaussians[ii * contractedGaussianB->contracted_GTO.size() + jj], primtiveOverlaps[2]);

            // then do (-,-,-,+,+,+)
            Huzinaga::computePrimitive(
                &contractedGaussianA->contracted_GTO[ii], xA, yA, zA, lxA - 1, lyA - 1, lzA - 1,
                &contractedGaussianB->contracted_GTO[jj], xB, yB, zB, lxB + 1, lyB + 1, lzB + 1,
                &productGaussians[ii * contractedGaussianB->contracted_GTO.size() + jj], primtiveOverlaps[3]);

            // finally do (+,+,+,+,+,+)
            Huzinaga::computePrimitive(
                &contractedGaussianA->contracted_GTO[ii], xA, yA, zA, lxA + 1, lyA + 1, lzA + 1,
                &contractedGaussianB->contracted_GTO[jj], xB, yB, zB, lxB + 1, lyB + 1, lzB + 1,
                &productGaussians[ii * contractedGaussianB->contracted_GTO.size() + jj], primtiveOverlaps[4]);

            // now contract the kinetic integrals
            tx = lxA * lxB * primtiveOverlaps[1][0];
            tx = tx + (-2 * expA * lxB  * primtiveOverlaps[2][0]);
            tx = tx + (-2 * expB * lxA  * primtiveOverlaps[3][0]);
            tx = tx + ( 4 * expB * expA * primtiveOverlaps[4][0]);

            ty = lyA * lyB * primtiveOverlaps[1][1];
            ty = ty + (-2 * expA * lyB  * primtiveOverlaps[2][1]);
            ty = ty + (-2 * expB * lyA  * primtiveOverlaps[3][1]);
            ty = ty + ( 4 * expB * expA * primtiveOverlaps[4][1]);

            tz = lzA * lzB * primtiveOverlaps[1][2];
            tz = tz + (-2 * expA * lzB  * primtiveOverlaps[2][2]);
            tz = tz + (-2 * expB * lzA  * primtiveOverlaps[3][2]);
            tz = tz + ( 4 * expB * expA * primtiveOverlaps[4][2]);

            primitiveKinetic[0] = primitiveKinetic[0] + (
                0.5 * tx * primtiveOverlaps[0][1] * primtiveOverlaps[0][2] * gaussianProduct.gaussian_integral[3] * pow(M_PI / (expA + expB), 1.5) * normA * coeffA * normB * coeffB
                );
            primitiveKinetic[1] = primitiveKinetic[1] + (
                0.5 * ty * primtiveOverlaps[0][2] * primtiveOverlaps[0][0] * gaussianProduct.gaussian_integral[3] * pow(M_PI / (expA + expB), 1.5) * normA * coeffA * normB * coeffB
                );
            primitiveKinetic[2] = primitiveKinetic[2] + (
                0.5 * tz * primtiveOverlaps[0][0] * primtiveOverlaps[0][1] * gaussianProduct.gaussian_integral[3] * pow(M_PI / (expA + expB), 1.5) * normA * coeffA * normB * coeffB
                );
        }
    }
    return primitiveKinetic[0] + primitiveKinetic[1] + primitiveKinetic[2];
}

// based on the reference implementation in dx.doi.org/doi:10.3888/tmj.14-3
// Evaluation of Gaussian Molecular Integrals I. Overlap Integrals
// original reference in
//
void Huzinaga::computePrimitive(cxx_Primitive *primitiveA, std::double_t xA, std::double_t yA, std::double_t zA, std::int64_t lxA, std::int64_t lyA, std::int64_t lzA, cxx_Primitive *primitiveB, std::double_t xB, std::double_t yB, std::double_t zB, std::int64_t lxB, std::int64_t lyB, std::int64_t lzB, cxx_Gaussians *productGaussian, std::double_t *primitiveOverlaps)
{
    std::double_t aux = 0.0;
    memset(primitiveOverlaps, 0, 4 * sizeof(std::double_t));

    // first do x-diection
    for (std::int64_t ii = 0; ii <= lxA; ii++)
    {
        for (std::int64_t jj = 0; jj <= lxB; jj++)
        {
            if ((ii + jj) % 2 == 0)
            {
                aux = combination(lxA, ii) * combination(lxB, jj) * doublefactorial(ii + jj - 1);
                aux = aux * pow(productGaussian->gaussian_center[0] - xA, lxA - ii) * pow(productGaussian->gaussian_center[0] - xB, lxB - jj);
                aux = aux / pow(2 * (primitiveA->primitive_exp + primitiveB->primitive_exp), 0.5 * (ii + jj));
                primitiveOverlaps[0] = primitiveOverlaps[0] + aux;
            }
            // std::cout << ii << " " << jj << " " << primitiveOverlaps[0] << "\n";
        }
    }
    // primitiveOverlaps[0] = primitiveOverlaps[0] * pow(M_PI / (primitiveA->primitive_exp + primitiveB->primitive_exp), 0.5);

    // then do y-diection
    for (std::int64_t ii = 0; ii <= lyA; ii++)
    {
        for (std::int64_t jj = 0; jj <= lyB; jj++)
        {
            if ((ii + jj) % 2 == 0)
            {
                aux = combination(lyA, ii) * combination(lyB, jj) * doublefactorial(ii + jj - 1);
                aux = aux * pow(productGaussian->gaussian_center[1] - yA, lyA - ii) * pow(productGaussian->gaussian_center[1] - yB, lyB - jj);
                aux = aux / pow(2 * (primitiveA->primitive_exp + primitiveB->primitive_exp), 0.5 * (ii + jj));
                primitiveOverlaps[1] = primitiveOverlaps[1] + aux;
            }
            // std::cout << ii << " " << jj << " " << primitiveOverlaps[1] << "\n";
        }
    }
    // primitiveOverlaps[1] = primitiveOverlaps[1] * pow(M_PI / (primitiveA->primitive_exp + primitiveB->primitive_exp), 0.5);

    // then do z-diection
    for (std::int64_t ii = 0; ii <= lzA; ii++)
    {
        for (std::int64_t jj = 0; jj <= lzB; jj++)
        {
            if ((ii + jj) % 2 == 0)
            {
                aux = combination(lzA, ii) * combination(lzB, jj) * doublefactorial(ii + jj - 1);
                aux = aux * pow(productGaussian->gaussian_center[2] - zA, lzA - ii) * pow(productGaussian->gaussian_center[2] - zB, lzB - jj);
                aux = aux / pow(2 * (primitiveA->primitive_exp + primitiveB->primitive_exp), 0.5 * (ii + jj));
                primitiveOverlaps[2] = primitiveOverlaps[2] + aux;
            }
            // std::cout << ii << " " << jj << " " << primitiveOverlaps[2] << "\n";
        }
    }
    // primitiveOverlaps[2] = primitiveOverlaps[2] * pow(M_PI / (primitiveA->primitive_exp + primitiveB->primitive_exp), 0.5);

    // now contract the integral
    primitiveOverlaps[3] = primitiveOverlaps[0] * primitiveOverlaps[1] * primitiveOverlaps[2] * productGaussian->gaussian_integral[3];
    primitiveOverlaps[3] = primitiveOverlaps[3] * primitiveA->orbital_coeff * primitiveA->orbital_norm;
    primitiveOverlaps[3] = primitiveOverlaps[3] * primitiveB->orbital_coeff * primitiveB->orbital_norm;
    primitiveOverlaps[3] = primitiveOverlaps[3] * pow(M_PI / (primitiveA->primitive_exp + primitiveB->primitive_exp), 1.5);
}