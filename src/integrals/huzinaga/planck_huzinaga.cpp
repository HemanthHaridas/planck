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

// based on the reference implementation in dx.doi.org/doi:10.3888/tmj.14-3
// Evaluation of Gaussian Molecular Integrals I. Overlap Integrals
// original reference in dx.doi.org/10.1143/JPSJ.21.2313
// Gaussian-Expansion Methods for Molecular Integrals
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
        }
    }
    return primtiveOverlap;
}

// based on the reference implementation in dx.doi.org/doi:10.3888/tmj.14-3
// Evaluation of Gaussian Molecular Integrals I. Overlap Integrals
// original reference in dx.doi.org/10.1143/JPSJ.21.2313
// Gaussian-Expansion Methods for Molecular Integrals
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

    std::double_t primtiveOverlaps[5][4];
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
            tx = tx + (-2 * expA * lxB * primtiveOverlaps[2][0]);
            tx = tx + (-2 * expB * lxA * primtiveOverlaps[3][0]);
            tx = tx + (4 * expB * expA * primtiveOverlaps[4][0]);

            ty = lyA * lyB * primtiveOverlaps[1][1];
            ty = ty + (-2 * expA * lyB * primtiveOverlaps[2][1]);
            ty = ty + (-2 * expB * lyA * primtiveOverlaps[3][1]);
            ty = ty + (4 * expB * expA * primtiveOverlaps[4][1]);

            tz = lzA * lzB * primtiveOverlaps[1][2];
            tz = tz + (-2 * expA * lzB * primtiveOverlaps[2][2]);
            tz = tz + (-2 * expB * lzA * primtiveOverlaps[3][2]);
            tz = tz + (4 * expB * expA * primtiveOverlaps[4][2]);

            primitiveKinetic[0] = primitiveKinetic[0] + (0.5 * tx * primtiveOverlaps[0][1] * primtiveOverlaps[0][2] * gaussianProduct.gaussian_integral[3] * pow(M_PI / (expA + expB), 1.5) * normA * coeffA * normB * coeffB);
            primitiveKinetic[1] = primitiveKinetic[1] + (0.5 * ty * primtiveOverlaps[0][2] * primtiveOverlaps[0][0] * gaussianProduct.gaussian_integral[3] * pow(M_PI / (expA + expB), 1.5) * normA * coeffA * normB * coeffB);
            primitiveKinetic[2] = primitiveKinetic[2] + (0.5 * tz * primtiveOverlaps[0][0] * primtiveOverlaps[0][1] * gaussianProduct.gaussian_integral[3] * pow(M_PI / (expA + expB), 1.5) * normA * coeffA * normB * coeffB);
        }
    }
    return primitiveKinetic[0] + primitiveKinetic[1] + primitiveKinetic[2];
}

std::double_t Huzinaga::computeNuclear(std::double_t *atomCoords, std::uint64_t *atomicNumbers, std::uint64_t nAtoms, cxx_Contracted *contractedGaussianA, cxx_Contracted *contractedGaussianB, std::error_code *errorFlag, std::string *errorMessage)
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

    // // auxiliary variables
    // std::double_t tx = 0, ty = 0, tz = 0;

    std::double_t primitiveNuclear[3];

    memset(primitiveNuclear, 0, 3 * sizeof(std::double_t));

    // compute the gaussian products
    computeGaussianProduct(contractedGaussianA, contractedGaussianB, &productGaussians, errorFlag, errorMessage);

    for (std::uint64_t ii = 0; ii < contractedGaussianA->contracted_GTO.size(); ii++)
    {
        for (std::uint64_t jj = 0; jj < contractedGaussianB->contracted_GTO.size(); jj++)
        {
            /* code */
        }
    }
}

// based on the reference implementation in dx.doi.org/doi:10.3888/tmj.14-3
// Evaluation of Gaussian Molecular Integrals I. Overlap Integrals
// original reference in dx.doi.org/10.1143/JPSJ.21.2313
// Gaussian-Expansion Methods for Molecular Integrals
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
        }
    }

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
        }
    };

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
        }
    }

    // now contract the integral
    primitiveOverlaps[3] = primitiveOverlaps[0] * primitiveOverlaps[1] * primitiveOverlaps[2] * productGaussian->gaussian_integral[3];
    primitiveOverlaps[3] = primitiveOverlaps[3] * primitiveA->orbital_coeff * primitiveA->orbital_norm;
    primitiveOverlaps[3] = primitiveOverlaps[3] * primitiveB->orbital_coeff * primitiveB->orbital_norm;
    primitiveOverlaps[3] = primitiveOverlaps[3] * pow(M_PI / (primitiveA->primitive_exp + primitiveB->primitive_exp), 1.5);
}

void Huzinaga::computeNuclear(std::double_t *atomCoords, std::uint64_t *atomicNumbers, std::uint64_t nAtoms, cxx_Primitive *primitiveA, std::double_t xA, std::double_t yA, std::double_t zA, std::int64_t lxA, std::int64_t lyA, std::int64_t lzA, cxx_Primitive *primitiveB, std::double_t xB, std::double_t yB, std::double_t zB, std::int64_t lxB, std::int64_t lyB, std::int64_t lzB, cxx_Gaussians *productGaussian, std::double_t *primitiveOverlaps)
{
    std::double_t gamma = primitiveA->primitive_exp + primitiveB->primitive_exp;
    std::double_t epsilon = pow(4 * gamma, -1);
    std::double_t *gaussCoord = productGaussian->gaussian_center;

    // calculate the maximum array size required to hold intermediates
    std::int64_t nX = lxA + lxB + 1;
    std::int64_t nY = lyA + lyB + 1;
    std::int64_t nZ = lzA + lzB + 1;

    std::int64_t xSize = (pow(nX, 3) / 16) + (pow(nX, 2) / 4) + nX;
    std::int64_t ySize = (pow(nY, 3) / 16) + (pow(nY, 2) / 4) + nY;
    std::int64_t zSize = (pow(nZ, 3) / 16) + (pow(nZ, 2) / 4) + nZ;

    // now setup buffers
    std::vector<std::double_t> xIntermediates(xSize);
    std::vector<std::double_t> yIntermediates(ySize);
    std::vector<std::double_t> zIntermediates(zSize);

    // now calculate the integral between the atomic charge distribution and the gaussians
    for (std::uint64_t ii = 0; ii < nAtoms; ii++)
    {
        // first extract atom coords
        std::double_t xCoord = atomCoords[ii * 3 + 0];
        std::double_t yCoord = atomCoords[ii * 3 + 1];
        std::double_t zCoord = atomCoords[ii * 3 + 2];

        computeNuclearPrimitive(xCoord, xA, lxA, xB, lxB, epsilon, gaussCoord[0], &xIntermediates);
        computeNuclearPrimitive(yCoord, yA, lyA, yB, lyB, epsilon, gaussCoord[1], &yIntermediates);
        computeNuclearPrimitive(zCoord, zA, lzA, zB, lzB, epsilon, gaussCoord[2], &zIntermediates);

        // now compute the combinations
        for (std::uint64_t ij = 0; ij < xIntermediates.size(); ij++)
        {
            for (std::uint64_t jk = 0; jk < yIntermediates.size(); jk++)
            {
                for (std::uint64_t kl = 0; kl < zIntermediates.size(); kl++)
                {
                    
                }
            }
        }
    }
}

void Huzinaga::computeNuclearPrimitive(std::double_t atomCoord, std::double_t coordA, std::int64_t shellA, std::double_t coordB, std::int64_t shellB, std::double_t gammaExp, std::double_t gaussCoord, std::vector<std::double_t> *primitiveResults)
{
    std::int64_t iiMax = shellA + shellB;
    std::double_t result = 0.0;
    for (std::uint64_t ii = 0; ii <= iiMax; ii++)
    {
        std::int64_t jjMax = ii / 2;
        for (std::uint64_t jj = 0; jj <= jjMax; jj++)
        {
            std::int64_t kkMax = (ii - 2 * (jjMax)) / 2;
            for (std::uint64_t kk = 0; kk <= kkMax; kk++)
            {
                Huzinaga::computeNuclearPrimitive2(ii, jj, kk, shellA, coordA, shellB, coordB, atomCoord, gaussCoord, gammaExp, &result);
                primitiveResults->push_back(result);
            }
        }
    }
}

void Huzinaga::computeNuclearPrimitive2(std::int64_t indexA, std::int64_t indexB, std::int64_t indexC, std::int64_t shellA, std::double_t coordA, std::int64_t shellB, std::double_t coordB, std::double_t atomCoord, std::double_t gaussCoord, std::double_t gammaExp, std::double_t *result)
{
    Huzinaga::computeNuclearPrimitive1(indexA, shellA, coordA, shellB, coordB, gaussCoord, result); // compute expansion coefficient

    std::int64_t  boysIndex = indexA - (2 * indexB) - indexC;
    std::double_t boysParam = pow(atomCoord - gaussCoord, 2);
    
    *(result) = *(result)*factorial(indexA) * pow(-1, indexA + indexC) * pow(gaussCoord - atomCoord, indexA - 2 * (indexB + indexC));
    *(result) = *(result)*pow(gammaExp, indexB + indexC);
    *(result) = *(result) / (factorial(indexB) * factorial(indexC));
    *(result) = *(result) / factorial(indexA - 2 * (indexB + indexC));
    *(result) = *(result) * boysFunction(boysIndex, boysParam);
}

void Huzinaga::computeNuclearPrimitive1(std::int64_t expandIndex, std::int64_t shellA, std::double_t coordA, std::int64_t shellB, std::double_t coordB, std::double_t gaussCoord, std::double_t *result)
{
    std::uint64_t iiMin = std::max(static_cast<std::int64_t>(0), expandIndex - shellB);
    std::uint64_t iiMax = std::min(expandIndex, shellA);

    for (std::uint64_t ii = iiMin; ii <= iiMax; ii++)
    {
        *result = combination(shellA, ii) * combination(shellB, expandIndex - ii);
        *result = *(result) + (*(result)*pow(gaussCoord - coordA, shellA - ii) * pow(gaussCoord - shellB, shellB + ii - expandIndex));
    }
}