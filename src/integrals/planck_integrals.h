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

// Base files
#include "../base/planck_base.h"

// Huzinaga
#include "huzinaga/planck_huzinaga.h"

// Hermite
#include "hermite/planck_hermite.h"

// Obara Sakia
#include "obarasakia/planck_obarasakia.h"
namespace IntegralEngine
{
    void computeOverlap(cxx_Calculator &planckCalculator, Eigen::MatrixXd &overlapMatrix);
    void computeKinetic(cxx_Calculator &planckCalculator, Eigen::MatrixXd &kineticMatrix);
};