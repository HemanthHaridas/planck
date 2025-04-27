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

#include "planck_hermite.h"

// std::double_t Hermite::expansionCoeff1(const std::int64_t expIndex, const std::int64_t shellA, const std::double_t centerA, const std::int64_t shellB, const std::double_t centerB, const std::double_t gaussCenter)
// {
//     std::double_t expansionCoeff = 0.0;
//     std::double_t aux;
//     std::int64_t cMin = std::max(static_cast<int64_t>(0), expIndex - shellB);
//     std::int64_t cMax = std::min(expIndex, shellA);

//     for (std::int64_t ii = cMin; ii <= cMax; ii++)
//     {
//         aux = combination(shellA, ii);
//         aux = aux * combination(shellB, expIndex - ii);
//         aux = aux * pow(gaussCenter - centerA, shellA - ii);
//         aux = aux * pow(gaussCenter - centerB, shellB + ii - expIndex);
//         expansionCoeff = expansionCoeff + aux;
//     }
//     return expansionCoeff;
// }

std::double_t Hermite::expansionCoeff1(std::int64_t shellA, std::int64_t shellB, std::int64_t nodes, std::double_t centerA, std::double_t centerB, std::double_t expA, std::double_t expB)
{
    std::double_t gamma = expA + expB;
    std::double_t combExp = (expA * expB) / gamma;

    if ((nodes < 0) || (nodes > (shellA + shellB)))
    {
        return 0;
    }
    else if (nodes == shellA == shellB == 0)
    {
        return exp(-1 * combExp * (centerA - centerB) * (centerA - centerB));
    }
    else if (shellB == 0)
    {
        return (1 / (2 * gamma)) * expansionCoeff1(shellA - 1, shellB, nodes - 1, centerA, centerB, expA, expB) -
               (combExp * (centerA - centerB) / expA) * expansionCoeff1(shellA - 1, shellB, nodes, centerA, centerB, expA, expB) +
               (nodes + 1) * expansionCoeff1(shellA - 1, shellB, nodes + 1, centerA, centerB, expA, expB);
    }
    else
    {
        return (1 / (2 * gamma)) * expansionCoeff1(shellA, shellB - 1, nodes - 1, centerA, centerB, expA, expB) -
               (combExp * (centerA - centerB) / expB) * expansionCoeff1(shellA, shellB - 1, nodes, centerA, centerB, expA, expB) +
               (nodes + 1) * expansionCoeff1(shellA, shellB - 1, nodes + 1, centerA, centerB, expA, expB);
    }
}

std::double_t auxiliaryHermite(std::int64_t xNodes, std::int64_t yNodes, std::int64_t zNodes, std::int64_t boysOrder, std::double_t exponent, std::double_t xDist, std::double_t yDist, std::double_t zDist)
{
    std::double_t auxHermite = 0.0;
    if (xNodes == yNodes == zNodes == 0)
    {
        std::double_t boysParam = exponent * ((xDist * xDist) + (yDist * yDist) + (zDist * zDist));
        auxHermite = auxHermite + pow(-2 * exponent, boysOrder) * boysFunction(boysOrder, boysParam);
    }
    else if (xNodes == yNodes == 0)
    {
        if (zNodes > 1)
        {
            auxHermite = auxHermite + (zNodes - 1) * auxiliaryHermite(xNodes, yNodes, zNodes - 2, boysOrder + 1, exponent, xDist, yDist, zDist);
        }
        auxHermite = auxHermite + (zDist)*auxiliaryHermite(xNodes, yNodes, zNodes - 1, boysOrder + 1, exponent, xDist, yDist, zDist);
    }
    else if (xNodes == 0)
    {
        if (yNodes > 1)
        {
            auxHermite = auxHermite + (yNodes - 1) * auxiliaryHermite(xNodes, yNodes - 2, zNodes, boysOrder + 1, exponent, xDist, yDist, zDist);
        }
        auxHermite = auxHermite + (yDist)*auxiliaryHermite(xNodes, yNodes - 1, zNodes, boysOrder + 1, exponent, xDist, yDist, zDist);
    }
    else
    {
        if (xNodes > 1)
        {
            auxHermite = auxHermite + (xNodes - 1) * auxiliaryHermite(xNodes - 2, yNodes, zNodes, boysOrder + 1, exponent, xDist, yDist, zDist);
        }
        auxHermite = auxHermite + (xDist)*auxiliaryHermite(xNodes - 1, yNodes, zNodes, boysOrder + 1, exponent, xDist, yDist, zDist);
    }
    return auxHermite;
}

