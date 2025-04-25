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
    // need to reset density matrix
    densityMatrix.resize(nBasis, nBasis);
    densityMatrix.setZero();

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

void generateHamiltoninan(const Eigen::MatrixXd &densityMatrix, const Eigen::Tensor<std::double_t, 4> &electronicMatrix, Eigen::MatrixXd &hamiltonianMatrix)
{
    const std::int64_t nBasis = densityMatrix.rows();

    for (std::int64_t ii = 0; ii < nBasis; ii++)
    {
        for (std::int64_t jj = 0; jj < nBasis; jj++)
        {
            double sum = 0.0;
            for (std::int64_t kk = 0; kk < nBasis; kk++)
            {
                for (std::int64_t ll = 0; ll < nBasis; ll++)
                {
                    double factor = electronicMatrix(ii, jj, kk, ll) - 0.5 * electronicMatrix(ii, ll, kk, jj);
                    sum += densityMatrix(kk, ll) * factor;
                }
            }
            hamiltonianMatrix(ii, jj) += sum;
        }
    }
}

void noDiisRHF(scfData *scfInstance, const Eigen::Tensor<std::double_t, 4> &electronicMatrix, const std::uint64_t nElectrons)
{
    // first recompute the hamiltonian
    generateHamiltoninan(scfInstance->densityMatrix, electronicMatrix, scfInstance->hamiltonianMatrix);

    // compute fock matrix and orthogonalize it
    scfInstance->fockMatrix = scfInstance->hamiltonianMatrix + scfInstance->coreMatrix;
    scfInstance->orthoFock = scfInstance->orthoMatrix.conjugate().transpose() * scfInstance->fockMatrix * scfInstance->orthoMatrix;

    // solve the eigen value problem
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(scfInstance->orthoFock);
    scfInstance->orbitalEnegies = eigenSolver.eigenvalues();
    scfInstance->orthogonalMO   = eigenSolver.eigenvectors();
    scfInstance->canonicalMO    = scfInstance->orthoMatrix * scfInstance->orthogonalMO;

    // store old density
    Eigen::MatrixXd oldDensity = scfInstance->densityMatrix;

    // recompute density
    computeRHFDensity(nElectrons, scfInstance->canonicalMO, scfInstance->densityMatrix);

    // std::cout << "Old Density" << oldDensity << "\n";
    // std::cout << "New Density" << scfInstance->densityMatrix << "\n";
}