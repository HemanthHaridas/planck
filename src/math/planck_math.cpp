#include "planck_math.h"
#include <iostream>

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

// This implementation of factorial will return 0 if the input
// is less than 0. This is not mathematically correct, but makes
// sense computationally
std::double_t factorial(const std::int64_t number)
{
    std::double_t result;

    if (number < 0)
    {
        result = 0;
    }

    result = 1;
    for (std::int64_t i = 1; i <= number; ++i)
    {
        result = result * i;
    }
    return result;
}

std::double_t doublefactorial(const std::int64_t number)
{
    std::double_t result;

    if (number < 0)
    {
        result = 0;
    }

    result = 1;
    for (std::int64_t i = number; i >= 1; i -= 2)
    {
        result = result * i;
    }
    return result;
}

std::double_t combination(const std::int64_t number, const std::int64_t choice)
{
    std::double_t result;

    if (choice > number)
    {
        result = 0;
    }

    std::double_t numerator;
    std::double_t denominator;

    numerator = factorial(number);
    denominator = factorial(choice);
    result = numerator / denominator;
    return result;
}
