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

#include "planck_math.h"

const inline std::double_t factorial(std::int64_t number)
{
    // check if number is less than zero
    if (number < 0)
    {
        return 0;
    }

    std::double_t result = 1.0;
    for (std::int64_t ii = 1; ii < number; ii++)
    {
        result = result * ii;
    }
    return result;
}

const inline std::double_t combination(std::int64_t number, std::int64_t choice)
{
    // check if number is less than zero
    if (number < choice)
    {
        return 0;
    }

    std::double_t numerator = factorial(number);
    std::double_t denominator = factorial(choice) * factorial(number - choice);
    return (numerator / denominator);
}

const inline std::double_t doublefactorial(std::int64_t number)
{
    // check if number is less than -1
    if (number < -1)
    {
        return 0;
    }

    std::double_t result = 1.0;
    for (std::int64_t ii = number; ii > 1; ii-=2)
    {
        result = result * ii;
    }
    return result;
}

std::double_t dotproduct(std::double_t xA, std::double_t yA, std::double_t zA, std::double_t xB, std::double_t yB, std::double_t zB)
{
    return (xA - xB) * (xA - xB) + (yA - yB) * (yA - yB) + (zA - zB) * (zA - zB);
}