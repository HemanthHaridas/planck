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

#include "planck_base.h"
#include "../integrals/planck_integrals.h"

#include <Eigen/Dense>

void noDiisRHF(scfData *scfInstance, const Eigen::Tensor<std::double_t, 4> &electronicMatrix, const std::uint64_t nElectrons);
void DiisRHF(scfData *scfInstance, const Eigen::Tensor<std::double_t, 4> &electronicMatrix, const std::uint64_t nElectrons, const Eigen::MatrixXd &overlapMatrix, const std::uint64_t diisDim);
void diisEngine(Eigen::MatrixXd &fockMatrix, const Eigen::MatrixXd &orthoMatrix, const Eigen::MatrixXd &densityMatrix, const Eigen::MatrixXd &overlapMatrix, std::vector<Eigen::MatrixXd> &fockMatrices, std::vector<Eigen::MatrixXd> &errorMatrices, const std::uint64_t diisDim);
