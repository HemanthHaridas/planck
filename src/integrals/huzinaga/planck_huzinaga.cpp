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
std::double_t Huzinaga::computeOverlap(cxx_Contracted *contractedGaussianA, cxx_Contracted *contractedGaussianB)
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
    computeGaussianProduct(contractedGaussianA, contractedGaussianB, &productGaussians);

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
std::double_t Huzinaga::computeKinetic(cxx_Contracted *contractedGaussianA, cxx_Contracted *contractedGaussianB)
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
    computeGaussianProduct(contractedGaussianA, contractedGaussianB, &productGaussians);

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

// Electron Nuclear Integrals
std::double_t Huzinaga::expansionCoeff2(const std::int64_t indexA, const std::int64_t indexB, const std::int64_t indexC, const std::int64_t shellA, const std::double_t centerA, const std::int64_t shellB, const std::double_t centerB, const std::double_t atomCenter, const std::double_t gaussCenter, std::double_t gamma)
{
    std::double_t epsilon = 1 / (4 * gamma);
    std::double_t expansionCoeff = Huzinaga::expansionCoeff1(indexA, shellA, centerA, shellB, centerB, gaussCenter);
    expansionCoeff = expansionCoeff * factorial(indexA);
    expansionCoeff = expansionCoeff * pow(-1, indexA + indexC);
    expansionCoeff = expansionCoeff * pow(gaussCenter - atomCenter, indexA - (2 * indexB) - (2 * indexC));
    expansionCoeff = expansionCoeff * pow(epsilon, indexB + indexC);
    expansionCoeff = expansionCoeff / factorial(indexC);
    expansionCoeff = expansionCoeff / factorial(indexB);
    expansionCoeff = expansionCoeff / factorial(indexA - (2 * indexB) - (2 * indexC));
    return expansionCoeff;
}

std::double_t Huzinaga::expansionCoeff1(const std::int64_t expIndex, const std::int64_t shellA, const std::double_t centerA, const std::int64_t shellB, const std::double_t centerB, const std::double_t gaussCenter)
{
    std::double_t expansionCoeff = 0.0;
    std::double_t aux;
    std::int64_t cMin = std::max(static_cast<int64_t>(0), expIndex - shellB);
    std::int64_t cMax = std::min(expIndex, shellA);

    for (std::int64_t ii = cMin; ii <= cMax; ii++)
    {
        aux = combination(shellA, ii);
        aux = aux * combination(shellB, expIndex - ii);
        aux = aux * pow(gaussCenter - centerA, shellA - ii);
        aux = aux * pow(gaussCenter - centerB, shellB + ii - expIndex);
        expansionCoeff = expansionCoeff + aux;
    }
    return expansionCoeff;
}

// Benchmark implementation
// Do not uncomment
// std::double_t Huzinaga::computePrimitive(cxx_Primitive *primitiveA, const std::double_t xA, const std::double_t yA, const std::double_t zA, const std::int64_t lxA, const std::int64_t lyA, const std::int64_t lzA, cxx_Primitive *primitiveB, const std::double_t xB, const std::double_t yB, const std::double_t zB, const std::int64_t lxB, const std::int64_t lyB, const std::int64_t lzB, const std::double_t *atomCoords, const std::uint64_t *atomCharges, const std::uint64_t nAtoms)
// {
//     std::double_t primitiveIntegral = 0.0;
//     std::double_t gamma = primitiveA->primitive_exp + primitiveB->primitive_exp;
//     cxx_Gaussians productGaussian = computeGaussianProduct(primitiveA, xA, yA, zA, primitiveB, xB, yB, zB);

//     for (std::uint64_t ii = 0; ii < nAtoms; ii++)
//     {
//         std::double_t xCoord = atomCoords[ii * 3 + 0];
//         std::double_t yCoord = atomCoords[ii * 3 + 1];
//         std::double_t zCoord = atomCoords[ii * 3 + 2];

//         for (std::int64_t xi = 0; xi <= lxA + lxB; xi++)
//         {
//             for (std::int64_t xj = 0; xj <= (xi / 2); xj++)
//             {
//                 for (std::int64_t xk = 0; xk <= (xi - 2 * xj) / 2; xk++)
//                 {
//                     std::double_t xDir = Huzinaga::expansionCoeff2(xi, xj, xk, lxA, xA, lxB, xB, xCoord, productGaussian.gaussian_center[0], gamma);
//                     std::cout << xi << " " << xj << " " << xk << " " << xDir << "\n";
//                     for (std::int64_t yi = 0; yi <= lyA + lyB; yi++)
//                     {
//                         for (std::int64_t yj = 0; yj <= (yi / 2); yj++)
//                         {
//                             for (std::int64_t yk = 0; yk <= (yi - 2 * yj) / 2; yk++)
//                             {
//                                 std::double_t yDir = Huzinaga::expansionCoeff2(yi, yj, yk, lyA, yA, lyB, yB, yCoord, productGaussian.gaussian_center[1], gamma);
//                                 std::cout << yi << " " << yj << " " << yk << " " << yDir << "\n";
//                                 for (std::int64_t zi = 0; zi <= lzA + lzB; zi++)
//                                 {
//                                     for (std::int64_t zj = 0; zj <= (zi / 2); zj++)
//                                     {
//                                         for (std::int64_t zk = 0; zk <= (zi - 2 * zj) / 2; zk++)
//                                         {
//                                             std::double_t zDir = Huzinaga::expansionCoeff2(zi, zj, zk, lzA, zA, lzB, zB, zCoord, productGaussian.gaussian_center[2], gamma);
//                                             std::cout << zi << " " << zj << " " << zk << " " << zDir << "\n";
//                                             std::double_t boysP = dotproduct(xCoord, yCoord, zCoord, productGaussian.gaussian_center[0], productGaussian.gaussian_center[1], productGaussian.gaussian_center[2]);
//                                             std::uint64_t boysI = (xi + yi + zi) - 2 * (xj + yj + zj) - (xk + yk + zk);
//                                             std::double_t boysF = boysFunction(boysI, boysP);
//                                             primitiveIntegral = primitiveIntegral + (xDir * yDir * zDir * boysF) * atomCharges[ii] * -1;
//                                         }
//                                     }
//                                 }
//                             }
//                         }
//                     }
//                 }
//             }
//         }
//     }
//     primitiveIntegral = primitiveIntegral * primitiveA->orbital_coeff * primitiveA->orbital_norm;
//     primitiveIntegral = primitiveIntegral * primitiveB->orbital_coeff * primitiveB->orbital_norm;
//     return primitiveIntegral * (2 * M_PI / gamma) * productGaussian.gaussian_integral[3];
// }

std::double_t Huzinaga::computeNuclear(std::double_t *atomCoords, std::uint64_t *atomCharges, std::uint64_t nAtoms, cxx_Contracted *contractedGaussianA, cxx_Contracted *contractedGaussianB)
{

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

    std::double_t primitiveIntegral = 0.0;

    for (std::uint64_t ii = 0; ii < contractedGaussianA->contracted_GTO.size(); ii++)
    {
        for (std::uint64_t jj = 0; jj < contractedGaussianB->contracted_GTO.size(); jj++)
        {
            std::double_t primitiveIntegrals = Huzinaga::computePrimitive(
                &contractedGaussianA->contracted_GTO[ii], xA, yA, zA, lxA, lyA, lzA,
                &contractedGaussianB->contracted_GTO[jj], xB, yB, zB, lxB, lyB, lzB,
                atomCoords, atomCharges, nAtoms);

            // contract the integral
            primitiveIntegral = primitiveIntegral + primitiveIntegrals;
        }
    }
    return primitiveIntegral;
}

std::vector<nuclearInt> Huzinaga::Intermediates(const std::int64_t shellA, const std::double_t coordA, const std::int64_t shellB, const std::double_t coordB, const std::double_t atomCoord, const std::double_t gaussCoord, const std::double_t gamma)
{
    std::vector<nuclearInt> results;
    // results.reserve();

    for (std::int64_t ii = 0; ii <= shellA + shellB; ii++)
    {
        for (std::int64_t jj = 0; jj <= (ii / 2); jj++)
        {
            for (std::int64_t kk = 0; kk <= (ii - 2 * jj) / 2; kk++)
            {
                nuclearInt result;
                result.result = Huzinaga::expansionCoeff2(ii, jj, kk, shellA, coordA, shellB, coordB, atomCoord, gaussCoord, gamma);
                result.x = ii;
                result.y = jj;
                result.z = kk;
                results.push_back(result);
            }
        }
    }
    return results;
}

std::double_t Huzinaga::computePrimitive(cxx_Primitive *primitiveA, const std::double_t xA, const std::double_t yA, const std::double_t zA, const std::int64_t lxA, const std::int64_t lyA, const std::int64_t lzA, cxx_Primitive *primitiveB, const std::double_t xB, const std::double_t yB, const std::double_t zB, const std::int64_t lxB, const std::int64_t lyB, const std::int64_t lzB, const std::double_t *atomCoords, const std::uint64_t *atomCharges, const std::uint64_t nAtoms)
{
    std::double_t primitiveIntegral = 0.0;
    std::double_t gamma = primitiveA->primitive_exp + primitiveB->primitive_exp;
    cxx_Gaussians productGaussian = computeGaussianProduct(primitiveA, xA, yA, zA, primitiveB, xB, yB, zB);

    for (std::uint64_t ii = 0; ii < nAtoms; ii++)
    {
        std::double_t xCoord = atomCoords[ii * 3 + 0];
        std::double_t yCoord = atomCoords[ii * 3 + 1];
        std::double_t zCoord = atomCoords[ii * 3 + 2];

        std::vector<nuclearInt> xDir = Huzinaga::Intermediates(lxA, xA, lxB, xB, xCoord, productGaussian.gaussian_center[0], gamma);
        std::vector<nuclearInt> yDir = Huzinaga::Intermediates(lyA, yA, lyB, yB, yCoord, productGaussian.gaussian_center[1], gamma);
        std::vector<nuclearInt> zDir = Huzinaga::Intermediates(lzA, zA, lzB, zB, zCoord, productGaussian.gaussian_center[2], gamma);

        for (auto xVal : xDir)
        {
            for (auto yVal : yDir)
            {
                for (auto zVal : zDir)
                {
                    std::double_t boysP = dotproduct(xCoord, yCoord, zCoord, productGaussian.gaussian_center[0], productGaussian.gaussian_center[1], productGaussian.gaussian_center[2]) * gamma;
                    std::uint64_t boysI = (xVal.x + yVal.x + zVal.x) - 2 * (xVal.y + yVal.y + zVal.y) - (xVal.z + yVal.z + zVal.z);
                    std::double_t boysF = boysFunction(boysI, boysP);
                    primitiveIntegral = primitiveIntegral + (xVal.result * yVal.result * zVal.result * boysF) * atomCharges[ii] * -1;
                }
            }
        }
    }
    primitiveIntegral = primitiveIntegral * primitiveA->orbital_coeff * primitiveA->orbital_norm;
    primitiveIntegral = primitiveIntegral * primitiveB->orbital_coeff * primitiveB->orbital_norm;
    return primitiveIntegral * (2 * M_PI / gamma) * productGaussian.gaussian_integral[3];
}

// Electron Repulsion Integrals
std::double_t Huzinaga::computePrimitive(
    cxx_Primitive *primitiveA, const std::double_t xA, const std::double_t yA, const std::double_t zA, const std::int64_t lxA, const std::double_t lyA, const std::double_t lzA,
    cxx_Primitive *primitiveB, const std::double_t xB, const std::double_t yB, const std::double_t zB, const std::int64_t lxB, const std::double_t lyB, const std::double_t lzB,
    cxx_Primitive *primitiveC, const std::double_t xC, const std::double_t yC, const std::double_t zC, const std::int64_t lxC, const std::double_t lyC, const std::double_t lzC,
    cxx_Primitive *primitiveD, const std::double_t xD, const std::double_t yD, const std::double_t zD, const std::int64_t lxD, const std::double_t lyD, const std::double_t lzD)
{
    std::double_t primitiveIntegral = 0.0;

    std::double_t gamma1 = primitiveA->primitive_exp + primitiveB->primitive_exp;
    std::double_t gamma2 = primitiveC->primitive_exp + primitiveD->primitive_exp;
    std::double_t f_gamma = (1 / gamma1) + (1 / gamma2);

    cxx_Gaussians productGaussianAB = computeGaussianProduct(primitiveA, xA, yA, zA, primitiveB, xB, yB, zB);
    cxx_Gaussians productGaussianCD = computeGaussianProduct(primitiveC, xC, yC, zC, primitiveD, xD, yD, zD);

    std::vector<electronInt> xDir = Huzinaga::Intermediates(lxA, xA, lxB, xB, lxC, xC, lxD, xD, productGaussianAB.gaussian_center[0], productGaussianCD.gaussian_center[0], gamma1, gamma2);
    std::vector<electronInt> yDir = Huzinaga::Intermediates(lyA, yA, lyB, yB, lyC, yC, lyD, yD, productGaussianAB.gaussian_center[1], productGaussianCD.gaussian_center[1], gamma1, gamma2);
    std::vector<electronInt> zDir = Huzinaga::Intermediates(lzA, zA, lzB, zB, lzC, zC, lzD, zD, productGaussianAB.gaussian_center[2], productGaussianCD.gaussian_center[2], gamma1, gamma2);

    for (auto xVal : xDir)
    {
        for (auto yVal : yDir)
        {
            for (auto zVal : zDir)
            {
                std::double_t boysP = dotproduct(productGaussianAB.gaussian_center[0], productGaussianAB.gaussian_center[1], productGaussianAB.gaussian_center[2], productGaussianCD.gaussian_center[0], productGaussianCD.gaussian_center[1], productGaussianCD.gaussian_center[2]) * (1 / f_gamma);
                std::uint64_t boysI = (xVal.i + yVal.i + zVal.i + xVal.k + yVal.k + zVal.k) - 2 * (xVal.j + yVal.j + zVal.j + xVal.l + yVal.l + zVal.l) - (xVal.m + yVal.m + zVal.m);
                std::double_t boysF = boysFunction(boysI, boysP);
                primitiveIntegral = primitiveIntegral + (xVal.result * yVal.result * zVal.result * boysF);
            }
        }
    }
    primitiveIntegral = primitiveIntegral * pow(M_PI, 2) * 2;
    primitiveIntegral = primitiveIntegral * (gamma1 * gamma2);
    primitiveIntegral = primitiveIntegral * sqrt(M_PI / (gamma1 + gamma2));
    primitiveIntegral = primitiveIntegral * productGaussianAB.gaussian_integral[3] * productGaussianCD.gaussian_integral[3];
    return primitiveIntegral;
}

std::double_t Huzinaga::expansionCoeff3(const std::int64_t expIndexA, const std::int64_t expIndexB, const std::int64_t shellA, const std::double_t centerA, const std::int64_t shellB, const std::double_t centerB, const std::double_t gaussCoordAB, const std::double_t gamma)
{
    std::double_t result = Huzinaga::expansionCoeff1(expIndexA, shellA, centerA, shellB, centerB, gaussCoordAB);
    result = result * factorial(expIndexA);
    result = result * pow(gamma, expIndexB - expIndexA);
    result = result / factorial(expIndexB);
    result = result / factorial(expIndexA - (2 * expIndexB));
    return result;
}

std::vector<electronInt> Huzinaga::Intermediates(const std::int64_t shellA, const std::double_t coordA, const std::int64_t shellB, const std::double_t coordB, const std::int64_t shellC, const std::double_t coordC, const std::int64_t shellD, const std::double_t coordD, const std::double_t gaussCoordAB, const std::double_t gaussCoordCD, const std::double_t gammaA, const std::double_t gammaB)
{
    std::vector<electronInt> results;
    std::double_t delta = (1 / (4 * gammaA)) + (1 / (4 * gammaB));

    for (std::int64_t ii = 0; ii <= (shellA + shellB); ii++) // l loop
    {
        for (std::int64_t jj = 0; jj <= (ii / 2); jj++) // r loop
        {
            std::double_t aux = Huzinaga::expansionCoeff3(ii, jj, shellA, coordA, shellB, coordB, gaussCoordAB, gammaA) * pow(-1, ii);
            for (std::int64_t kk = 0; kk <= (shellC + shellD); kk++) // l' loop
            {
                for (std::int64_t ll = 0; ll <= (kk / 2); ll++) // r' loop
                {
                    aux = aux * Huzinaga::expansionCoeff3(kk, ll, shellC, coordC, shellD, coordD, gaussCoordCD, gammaB);
                    for (std::int64_t mm = 0; mm <= ((ii + ll - 2 * (jj + kk)) / 2); mm++) // i loop
                    {
                        electronInt result;
                        // calculate the value of the term
                        result.result = aux * pow(-1, mm) * pow(2 * delta, 2 * (jj + ll));
                        result.result = result.result * factorial(ii + kk - 2 * (jj + ll)) * pow(delta, mm);
                        result.result = result.result * pow(gaussCoordAB - gaussCoordCD, ii + kk - 2 * (jj + ll + mm));
                        result.result = result.result / pow(4 * delta, 2 * (ii + kk));
                        result.result = result.result / factorial(mm);
                        result.result = result.result / factorial(ii + kk - 2 * (jj + ll + mm));

                        // set up the indices
                        result.i = ii;
                        result.j = jj;
                        result.k = kk;
                        result.l = ll;
                        result.m = ii;

                        results.push_back(result);
                    }
                }
            }
        }
    }
    return results;
}

std::double_t Huzinaga::computeElectronic(cxx_Contracted *contractedGaussianA, cxx_Contracted *contractedGaussianB, cxx_Contracted *contractedGaussianC, cxx_Contracted *contractedGaussianD)
{
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

    std::double_t xC = contractedGaussianC->location_x;
    std::double_t yC = contractedGaussianC->location_y;
    std::double_t zC = contractedGaussianC->location_z;

    std::int64_t lxC = contractedGaussianC->shell_x;
    std::int64_t lyC = contractedGaussianC->shell_y;
    std::int64_t lzC = contractedGaussianC->shell_z;

    std::double_t xD = contractedGaussianD->location_x;
    std::double_t yD = contractedGaussianD->location_y;
    std::double_t zD = contractedGaussianD->location_z;

    std::int64_t lxD = contractedGaussianD->shell_x;
    std::int64_t lyD = contractedGaussianD->shell_y;
    std::int64_t lzD = contractedGaussianD->shell_z;

    std::double_t primitiveIntegral = 0.0;

    for (std::uint64_t ii = 0; ii < contractedGaussianA->contracted_GTO.size(); ii++)
    {
        for (std::uint64_t jj = 0; jj < contractedGaussianB->contracted_GTO.size(); jj++)
        {
            for (std::uint64_t kk = 0; kk < contractedGaussianC->contracted_GTO.size(); kk++)
            {
                for (std::uint64_t ll = 0; ll < contractedGaussianD->contracted_GTO.size(); ll++)
                {
                    std::double_t primitiveIntegrals = Huzinaga::computePrimitive(
                        &contractedGaussianA->contracted_GTO[ii], xA, yA, zA, lxA, lyA, lzA,
                        &contractedGaussianB->contracted_GTO[jj], xB, yB, zB, lxB, lyB, lzB,
                        &contractedGaussianC->contracted_GTO[kk], xC, yC, zC, lxC, lyC, lzC,
                        &contractedGaussianD->contracted_GTO[ll], xD, yD, zD, lxD, lyD, lzD);
                    
                    primitiveIntegral = primitiveIntegral + primitiveIntegrals;
                }
            }
        }
    }
    return primitiveIntegral;
}