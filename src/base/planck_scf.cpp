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

#include "planck_scf.h"

void computeUHFDensity(const std::uint64_t nElectrons, Eigen::MatrixXd &canonicalCoeffs, Eigen::MatrixXd &densityMatrix)
{
    const std::uint64_t nBasis = canonicalCoeffs.rows();

    // compute density from coefficient matrix
    for (std::uint64_t row = 0; row < nBasis; row++)
    {
        for (std::uint64_t col = 0; col < nBasis; col++)
        {
            for (std::uint64_t elec = 0; elec < nElectrons; elec++)
            {
                densityMatrix(row, col) = canonicalCoeffs(row, elec) * canonicalCoeffs(col, elec) + densityMatrix(row, col);
            }
        }
    }
}

void computeRHFDensity(const std::uint64_t nElectrons, Eigen::MatrixXd &canonicalCoeffs, Eigen::MatrixXd &densityMatrix)
{
    const std::uint64_t nBasis = canonicalCoeffs.rows();

    // compute density from coefficient matrix
    for (std::uint64_t row = 0; row < nBasis; row++)
    {
        for (std::uint64_t col = 0; col < nBasis; col++)
        {
            for (std::uint64_t elec = 0; elec < static_cast<std::int64_t>(nElectrons / 2); elec++)
            {
                densityMatrix(row, col) = (2 * canonicalCoeffs(row, elec) * canonicalCoeffs(col, elec)) + densityMatrix(row, col);
            }
        }
    }
}

scfData scfStep()
{
}