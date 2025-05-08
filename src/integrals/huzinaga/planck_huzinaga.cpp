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

#include "planck_huzinaga.h"

std::double_t Huzinaga::Overlap::expansionIndex1(std::int64_t expIndex, std::int64_t shellA, std::int64_t shellB, std::double_t distPA, std::double_t distPB)
{
    std::int64_t cMax = expIndex <= shellA ? expIndex : shellA;
    std::int64_t cMin = 0 > (expIndex - shellB) ? 0 : (expIndex - shellB);
    std::double_t expansionCoeff = 0.0;

    for (std::int64_t ii = cMin; ii <= cMax; ii++)
    {
        std::double_t aux = combination(shellA, ii);
        aux = aux * combination(shellB, expIndex - ii);
        aux = aux * pow(distPA, shellA - ii);
        aux = aux * pow(distPB, shellB + ii - expIndex);
        expansionCoeff = expansionCoeff + aux;
    }
    return expansionCoeff;
}

std::double_t Huzinaga::Overlap::computePrimitive1D(std::double_t exponentA, std::double_t centerA, std::int64_t shellA, std::double_t exponentB, std::double_t centerB, std::int64_t shellB, std::double_t gaussianCenter)
{
    std::double_t integral = 0.0;

    std::double_t gaussian = ((exponentA * centerA) + (exponentB * centerB)) / (exponentA + exponentB);

    for (std::int64_t ii = 0; ii <= (shellA + shellB) / 2; ii++)
    {
        std::double_t value = doublefactorial((2 * ii) - 1);
        value = value / pow(2 * (exponentA + exponentB), ii);
        value = value * Huzinaga::Overlap::expansionIndex1(2 * ii, shellA, shellB, gaussian - centerA, gaussian - centerB);
        integral = integral + value;
    }

    return integral;
}

// std::double_t Huzinaga::Overlap::computePrimitive1D(std::double_t exponentA, std::double_t centerA, std::int64_t shellA, std::double_t exponentB, std::double_t centerB, std::int64_t shellB, std::double_t gaussianCenter)
// {
//     std::double_t integral = 0.0;
//     for (std::int64_t ii = 0; ii <= shellA; ii++)
//     {
//         for (std::int64_t jj = 0; jj <= shellB; jj++)
//         {
//             if ((ii + jj) % 2 == 0)
//             {
//                 std::double_t value = combination(shellA, ii);
//                 value = value * combination(shellB, jj);
//                 value = value * doublefactorial((ii + jj) - 1);
//                 value = value * pow(gaussianCenter - centerA, shellA - ii);
//                 value = value * pow(gaussianCenter - centerB, shellB - jj);
//                 value = value / pow(2 * (exponentA + exponentB), 0.5 * (ii + jj));
//                 integral = integral + value;
//             }
//         }
//     }
//     return integral;
// }

std::double_t Huzinaga::Overlap::computePrimitive3D(
    cxx_Primitive primitiveA, std::double_t xA, std::double_t yA, std::double_t zA, std::int64_t lxA, std::int64_t lyA, std::int64_t lzA,
    cxx_Primitive primitiveB, std::double_t xB, std::double_t yB, std::double_t zB, std::int64_t lxB, std::int64_t lyB, std::int64_t lzB,
    std::double_t *gaussianCenter, std::double_t gaussianIntegral)
{
    std::double_t xDir = Huzinaga::Overlap::computePrimitive1D(primitiveA.primitive_exp, xA, lxA, primitiveB.primitive_exp, xB, lxB, gaussianCenter[0]);
    std::double_t yDir = Huzinaga::Overlap::computePrimitive1D(primitiveA.primitive_exp, yA, lyA, primitiveB.primitive_exp, yB, lyB, gaussianCenter[1]);
    std::double_t zDir = Huzinaga::Overlap::computePrimitive1D(primitiveA.primitive_exp, zA, lzA, primitiveB.primitive_exp, zB, lzB, gaussianCenter[2]);

    return xDir * yDir * zDir * pow(M_PI / (primitiveA.primitive_exp + primitiveB.primitive_exp), 1.5) * gaussianIntegral;
}

std::double_t Huzinaga::Overlap::computeContracted(cxx_Contracted contractedGaussianA, cxx_Contracted contractedGaussianB)
{
    std::double_t integral = 0.0;

    for (std::uint64_t ii = 0; ii < contractedGaussianA.contracted_GTO.size(); ii++)
    {
        for (std::uint64_t jj = 0; jj < contractedGaussianB.contracted_GTO.size(); jj++)
        {
            cxx_Gaussians productGaussianAB = computeGaussianProduct(
                contractedGaussianA.contracted_GTO[ii], contractedGaussianA.location_x, contractedGaussianA.location_y, contractedGaussianA.location_z,
                contractedGaussianB.contracted_GTO[jj], contractedGaussianB.location_x, contractedGaussianB.location_y, contractedGaussianB.location_z);

            std::double_t value = Huzinaga::Overlap::computePrimitive3D(
                contractedGaussianA.contracted_GTO[ii], contractedGaussianA.location_x, contractedGaussianA.location_y, contractedGaussianA.location_z, contractedGaussianA.shell_x, contractedGaussianA.shell_y, contractedGaussianA.shell_z,
                contractedGaussianB.contracted_GTO[jj], contractedGaussianB.location_x, contractedGaussianB.location_y, contractedGaussianB.location_z, contractedGaussianB.shell_x, contractedGaussianB.shell_y, contractedGaussianB.shell_z,
                productGaussianAB.gaussian_center, productGaussianAB.gaussian_integral[3]);

            value = value * contractedGaussianA.contracted_GTO[ii].orbital_coeff * contractedGaussianA.contracted_GTO[ii].orbital_norm;
            value = value * contractedGaussianB.contracted_GTO[jj].orbital_coeff * contractedGaussianB.contracted_GTO[jj].orbital_norm;
            value = value * 1;

            integral = integral + value;
        }
    }
    return integral;
}

// std::double_t Huzinaga::Kinetic::computePrimitive3D(
//     cxx_Primitive primitiveA, std::double_t xA, std::double_t yA, std::double_t zA, std::int64_t lxA, std::int64_t lyA, std::int64_t lzA,
//     cxx_Primitive primitiveB, std::double_t xB, std::double_t yB, std::double_t zB, std::int64_t lxB, std::int64_t lyB, std::int64_t lzB,
//     std::double_t *gaussianCenter, std::double_t gaussianIntegral)
// {
//     std::double_t integral = Huzinaga::Overlap::computePrimitive3D(primitiveA, xA, yA, zA, lxA, lyA, lzA, primitiveB, xB, yB, zB, lxB, lyB, lzB, gaussianCenter, gaussianIntegral);
//     integral = integral * (primitiveB.primitive_exp * (2 * (lxB + lyB + lzB) + 3));

//     integral = integral - ((2 * pow(primitiveB.primitive_exp, 2)) * (Huzinaga::Overlap::computePrimitive3D(primitiveA, xA, yA, zA, lxA, lyA, lzA, primitiveB, xB, yB, zB, lxB + 2, lyB, lzB, gaussianCenter, gaussianIntegral)));
//     integral = integral - ((2 * pow(primitiveB.primitive_exp, 2)) * (Huzinaga::Overlap::computePrimitive3D(primitiveA, xA, yA, zA, lxA, lyA, lzA, primitiveB, xB, yB, zB, lxB, lyB + 2, lzB, gaussianCenter, gaussianIntegral)));
//     integral = integral - ((2 * pow(primitiveB.primitive_exp, 2)) * (Huzinaga::Overlap::computePrimitive3D(primitiveA, xA, yA, zA, lxA, lyA, lzA, primitiveB, xB, yB, zB, lxB, lyB, lzB + 2, gaussianCenter, gaussianIntegral)));

//     integral = integral - ((0.5 * lxB * (lxB - 1)) * (Huzinaga::Overlap::computePrimitive3D(primitiveA, xA, yA, zA, lxA, lyA, lzA, primitiveB, xB, yB, zB, lxB - 2, lyB, lzB, gaussianCenter, gaussianIntegral)));
//     integral = integral - ((0.5 * lyB * (lyB - 1)) * (Huzinaga::Overlap::computePrimitive3D(primitiveA, xA, yA, zA, lxA, lyA, lzA, primitiveB, xB, yB, zB, lxB, lyB - 2, lzB, gaussianCenter, gaussianIntegral)));
//     integral = integral - ((0.5 * lzB * (lzB - 1)) * (Huzinaga::Overlap::computePrimitive3D(primitiveA, xA, yA, zA, lxA, lyA, lzA, primitiveB, xB, yB, zB, lxB, lyB, lzB - 2, gaussianCenter, gaussianIntegral)));

//     // std::double_t exponentA = primitiveA.primitive_exp;
//     // std::double_t exponentB = primitiveB.primitive_exp;
//     // std::cout << std::setw(10) << std::setprecision(3) << std::right << xA << std::setw(10) << std::setprecision(3) << std::right << yA << std::setw(10) << std::setprecision(3) << std::right << zA << std::setw(10) << std::setprecision(3) << std::right << xB << std::setw(10) << std::setprecision(3) << std::right << yB << std::setw(10) << std::setprecision(3) << std::right << zB << std::setw(20) << std::setprecision(3) << std::right << integral << std::setw(20) << std::right << lxA << std::setw(20) << std::right << lyA << std::setw(20) << std::right << lzA << std::setw(20) << std::right << lxB << std::setw(20) << std::right << lyB << std::setw(20) << std::right << lzB << "\n";
//     // std::cout << std::setw(20) << std::setprecision(3) << std::right << exponentA << std::setw(20) << std::setprecision(3) << std::right << exponentB << std::setw(20) << std::setprecision(3) << std::right << integral << std::setw(20) << std::right << lxA << std::setw(20) << std::right << lxB << std::setw(20) << std::right << lyA << std::setw(20) << std::right << lyB << std::setw(20) << std::right << lzA << std::setw(20) << std::right << lzB << "\n";
//     return integral;
// }

std::double_t Huzinaga::Kinetic::computePrimitive3D(
    cxx_Primitive primitiveA, std::double_t xA, std::double_t yA, std::double_t zA, std::int64_t lxA, std::int64_t lyA, std::int64_t lzA,
    cxx_Primitive primitiveB, std::double_t xB, std::double_t yB, std::double_t zB, std::int64_t lxB, std::int64_t lyB, std::int64_t lzB,
    std::double_t *gaussianCenter, std::double_t gaussianIntegral)
{
    std::double_t overlapX = Huzinaga::Overlap::computePrimitive1D(primitiveA.primitive_exp, xA, lxA, primitiveB.primitive_exp, xB, lxB, gaussianCenter[0]);
    std::double_t overlapY = Huzinaga::Overlap::computePrimitive1D(primitiveA.primitive_exp, yA, lyA, primitiveB.primitive_exp, yB, lyB, gaussianCenter[1]);
    std::double_t overlapZ = Huzinaga::Overlap::computePrimitive1D(primitiveA.primitive_exp, zA, lzA, primitiveB.primitive_exp, zB, lzB, gaussianCenter[2]);

    // now do (-,+)
    std::double_t overlapXMP = Huzinaga::Overlap::computePrimitive1D(primitiveA.primitive_exp, xA, lxA - 1, primitiveB.primitive_exp, xB, lxB + 1, gaussianCenter[0]);
    std::double_t overlapYMP = Huzinaga::Overlap::computePrimitive1D(primitiveA.primitive_exp, yA, lyA - 1, primitiveB.primitive_exp, yB, lyB + 1, gaussianCenter[1]);
    std::double_t overlapZMP = Huzinaga::Overlap::computePrimitive1D(primitiveA.primitive_exp, zA, lzA - 1, primitiveB.primitive_exp, zB, lzB + 1, gaussianCenter[2]);

    // now do (+,+)
    std::double_t overlapXPP = Huzinaga::Overlap::computePrimitive1D(primitiveA.primitive_exp, xA, lxA + 1, primitiveB.primitive_exp, xB, lxB + 1, gaussianCenter[0]);
    std::double_t overlapYPP = Huzinaga::Overlap::computePrimitive1D(primitiveA.primitive_exp, yA, lyA + 1, primitiveB.primitive_exp, yB, lyB + 1, gaussianCenter[1]);
    std::double_t overlapZPP = Huzinaga::Overlap::computePrimitive1D(primitiveA.primitive_exp, zA, lzA + 1, primitiveB.primitive_exp, zB, lzB + 1, gaussianCenter[2]);

    // now do (-,-)
    std::double_t overlapXMM = Huzinaga::Overlap::computePrimitive1D(primitiveA.primitive_exp, xA, lxA - 1, primitiveB.primitive_exp, xB, lxB - 1, gaussianCenter[0]);
    std::double_t overlapYMM = Huzinaga::Overlap::computePrimitive1D(primitiveA.primitive_exp, yA, lyA - 1, primitiveB.primitive_exp, yB, lyB - 1, gaussianCenter[1]);
    std::double_t overlapZMM = Huzinaga::Overlap::computePrimitive1D(primitiveA.primitive_exp, zA, lzA - 1, primitiveB.primitive_exp, zB, lzB - 1, gaussianCenter[2]);

    // now do (+,-)
    std::double_t overlapXPM = Huzinaga::Overlap::computePrimitive1D(primitiveA.primitive_exp, xA, lxA + 1, primitiveB.primitive_exp, xB, lxB - 1, gaussianCenter[0]);
    std::double_t overlapYPM = Huzinaga::Overlap::computePrimitive1D(primitiveA.primitive_exp, yA, lyA + 1, primitiveB.primitive_exp, yB, lyB - 1, gaussianCenter[1]);
    std::double_t overlapZPM = Huzinaga::Overlap::computePrimitive1D(primitiveA.primitive_exp, zA, lzA + 1, primitiveB.primitive_exp, zB, lzB - 1, gaussianCenter[2]);

    std::double_t kx = lxA * lxB * overlapXMM;
    std::double_t ky = lyA * lyB * overlapYMM;
    std::double_t kz = lzA * lzB * overlapZMM;

    kx = kx - (2 * primitiveA.primitive_exp * lxB * overlapXPM);
    ky = ky - (2 * primitiveA.primitive_exp * lyB * overlapYPM);
    kz = kz - (2 * primitiveA.primitive_exp * lzB * overlapZPM);

    kx = kx - (2 * primitiveB.primitive_exp * lxA * overlapXMP);
    ky = ky - (2 * primitiveB.primitive_exp * lyA * overlapYMP);
    kz = kz - (2 * primitiveB.primitive_exp * lzA * overlapZMP);

    kx = kx + (4 * primitiveB.primitive_exp * primitiveA.primitive_exp * overlapXPP);
    ky = ky + (4 * primitiveB.primitive_exp * primitiveA.primitive_exp * overlapYPP);
    kz = kz + (4 * primitiveB.primitive_exp * primitiveA.primitive_exp * overlapZPP);

    std::double_t Kx = 0.5 * kx * overlapY * overlapZ * pow(M_PI / (primitiveA.primitive_exp + primitiveB.primitive_exp), 1.5);
    std::double_t Ky = 0.5 * ky * overlapZ * overlapX * pow(M_PI / (primitiveA.primitive_exp + primitiveB.primitive_exp), 1.5);
    std::double_t Kz = 0.5 * kz * overlapX * overlapY * pow(M_PI / (primitiveA.primitive_exp + primitiveB.primitive_exp), 1.5);

    return (Kx + Ky + Kz) * gaussianIntegral;
}
std::double_t Huzinaga::Kinetic::computeContracted(cxx_Contracted contractedGaussianA, cxx_Contracted contractedGaussianB)
{
    std::double_t integral = 0.0;

    for (std::uint64_t ii = 0; ii < contractedGaussianA.contracted_GTO.size(); ii++)
    {
        for (std::uint64_t jj = 0; jj < contractedGaussianB.contracted_GTO.size(); jj++)
        {
            cxx_Gaussians productGaussianAB = computeGaussianProduct(
                contractedGaussianA.contracted_GTO[ii], contractedGaussianA.location_x, contractedGaussianA.location_y, contractedGaussianA.location_z,
                contractedGaussianB.contracted_GTO[jj], contractedGaussianB.location_x, contractedGaussianB.location_y, contractedGaussianB.location_z);

            std::double_t value = Huzinaga::Kinetic::computePrimitive3D(
                contractedGaussianA.contracted_GTO[ii], contractedGaussianA.location_x, contractedGaussianA.location_y, contractedGaussianA.location_z, contractedGaussianA.shell_x, contractedGaussianA.shell_y, contractedGaussianA.shell_z,
                contractedGaussianB.contracted_GTO[jj], contractedGaussianB.location_x, contractedGaussianB.location_y, contractedGaussianB.location_z, contractedGaussianB.shell_x, contractedGaussianB.shell_y, contractedGaussianB.shell_z,
                productGaussianAB.gaussian_center, productGaussianAB.gaussian_integral[3]);

            value = value * contractedGaussianA.contracted_GTO[ii].orbital_coeff * contractedGaussianA.contracted_GTO[ii].orbital_norm;
            value = value * contractedGaussianB.contracted_GTO[jj].orbital_coeff * contractedGaussianB.contracted_GTO[jj].orbital_norm;
            value = value * 1;

            integral = integral + value;
        }
    }
    return integral;
}