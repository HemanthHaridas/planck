#pragma once

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

#include <cmath>
#include <cstdint>
#include <vector>

std::uint64_t factorial(std::int64_t number);
std::uint64_t combination(std::int64_t number, std::int64_t choice);
std::uint64_t doublefactorial(std::int64_t number);
std::double_t dotproduct(std::double_t xA, std::double_t yA, std::double_t zA, std::double_t xB, std::double_t yB, std::double_t zB);

template <typename T>
T reduce2D(const std::vector <T> &vecA, const std::vector <T> &vecB)
{
    T result = 0;

    // now comute the outer product and reduce
    for (std::uint64_t ii = 0; ii < vecA.size(); ii++)
    {
        for (std::uint64_t jj = 0; jj < vecB.size(); jj++)
        {
            result = result + (vecA[ii] * vecB[jj]);
        }
    }
    return result;
}

template <typename T>
T reduce1D(const std::vector <T> &vecA, const T &value)
{
    T result = 0;
    for (std::uint64_t ii = 0; ii < vecA.size(); ii++)
    {
        result = result + (vecA[ii] * value);
    }
    return result;
}