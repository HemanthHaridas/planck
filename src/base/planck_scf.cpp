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
    hamiltonianMatrix.setZero();

    for (std::int64_t ii = 0; ii < nBasis; ii++)
    {
        for (std::int64_t jj = 0; jj < nBasis; jj++)
        {
            for (std::int64_t kk = 0; kk < nBasis; kk++)
            {
                for (std::int64_t ll = 0; ll < nBasis; ll++)
                {
                    hamiltonianMatrix(ii, jj) = hamiltonianMatrix(ii, jj) + (densityMatrix(kk, ll) * (electronicMatrix(ii, jj, kk, ll) - 0.5 * electronicMatrix(ii, ll, kk, jj)));
                    // sum += densityMatrix(kk, ll) * factor;
                }
            }
            // hamiltonianMatrix(ii, jj) += sum;
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
    scfInstance->orthogonalMO = eigenSolver.eigenvectors();
    // std::cout << eigenSolver.eigenvalues() << "\n";
    scfInstance->canonicalMO = scfInstance->orthoMatrix * scfInstance->orthogonalMO;

    // store old density
    Eigen::MatrixXd oldDensity = scfInstance->densityMatrix;

    // recompute density
    computeRHFDensity(nElectrons, scfInstance->canonicalMO, scfInstance->densityMatrix);

    // std::cout << oldDensity << "\n";
    scfInstance->maxDensity = (scfInstance->densityMatrix - oldDensity).maxCoeff();
    scfInstance->rmsDensity = (scfInstance->densityMatrix - oldDensity).norm();
}

void DiisRHF(scfData *scfInstance, const Eigen::Tensor<std::double_t, 4> &electronicMatrix, const std::uint64_t nElectrons, const Eigen::MatrixXd &overlapMatrix, const std::uint64_t diisDim)
{
    // first recompute the hamiltonian
    generateHamiltoninan(scfInstance->densityMatrix, electronicMatrix, scfInstance->hamiltonianMatrix);

    // compute fock matrix and orthogonalize it
    scfInstance->fockMatrix = scfInstance->hamiltonianMatrix + scfInstance->coreMatrix;

    // perform direct iversion of iterative subspace
    diisEngine(scfInstance->fockMatrix, scfInstance->orthoMatrix, scfInstance->densityMatrix, overlapMatrix, scfInstance->fockMatrices, scfInstance->errorVectors, diisDim);
    // std::cout << scfInstance->fockMatrix << "\n";
    scfInstance->orthoFock = scfInstance->orthoMatrix.conjugate().transpose() * scfInstance->fockMatrix * scfInstance->orthoMatrix;

    // solve the eigen value problem
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(scfInstance->orthoFock);
    scfInstance->orbitalEnegies = eigenSolver.eigenvalues();
    scfInstance->orthogonalMO = eigenSolver.eigenvectors();
    // std::cout << eigenSolver.eigenvalues() << "\n";
    scfInstance->canonicalMO = scfInstance->orthoMatrix * scfInstance->orthogonalMO;

    // store old density
    Eigen::MatrixXd oldDensity = scfInstance->densityMatrix;

    // recompute density
    computeRHFDensity(nElectrons, scfInstance->canonicalMO, scfInstance->densityMatrix);

    // std::cout << oldDensity << "\n";
    scfInstance->maxDensity = (scfInstance->densityMatrix - oldDensity).maxCoeff();
    scfInstance->rmsDensity = (scfInstance->densityMatrix - oldDensity).norm();
}

// void soscfRHF()
// {

// }

void diisEngine(Eigen::MatrixXd &fockMatrix, const Eigen::MatrixXd &orthoMatrix, const Eigen::MatrixXd &densityMatrix, const Eigen::MatrixXd &overlapMatrix, std::vector<Eigen::MatrixXd> &fockMatrices, std::vector<Eigen::MatrixXd> &errorMatrices, const std::uint64_t diisDim)
{
    Eigen::MatrixXd tempVector = (fockMatrix * densityMatrix * overlapMatrix) - (overlapMatrix * densityMatrix * fockMatrix);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(overlapMatrix);
    Eigen::VectorXd eigenvalues = solver.eigenvalues();
    Eigen::MatrixXd eigenvectors = solver.eigenvectors();

    // Compute halfeigmat
    Eigen::MatrixXd halfeigmat = eigenvalues.array().inverse().sqrt().matrix().asDiagonal();

    // Compute rightmat
    Eigen::MatrixXd rightmat = eigenvectors * halfeigmat;

    // Compute orthomat
    Eigen::MatrixXd orthomat = rightmat * eigenvectors.transpose();

    // Compute t_errorvector
    Eigen::MatrixXd t_errorvector = (fockMatrix * (densityMatrix * overlapMatrix)) - (overlapMatrix * (densityMatrix * fockMatrix));

    // Compute leftmat
    Eigen::MatrixXd leftmat = orthomat.transpose() * t_errorvector;

    // Compute errorvector
    Eigen::MatrixXd errorVector = leftmat * orthomat;
    // Eigen::MatrixXd errorVector = orthoMatrix.transpose() * tempVector * orthoMatrix;
    // std::cout << errorVector << "\n";

    fockMatrices.push_back(fockMatrix);
    errorMatrices.push_back(errorVector);

    std::uint64_t fockDim = fockMatrices.size();
    // if the number of elements in fock matrix in greater than max fock
    // remove the oldest element
    if (fockDim >= diisDim)
    {
        fockMatrices.erase(fockMatrices.begin());
        errorMatrices.erase(errorMatrices.begin());
        fockDim--;
    }
    Eigen::MatrixXd bMatrix(fockDim + 1, fockDim + 1);
    Eigen::VectorXd residueVect(fockDim + 1);

    bMatrix.row(fockDim).setConstant(-1); // Last row
    bMatrix.col(fockDim).setConstant(-1); // Last column
    bMatrix(fockDim, fockDim) = 0;
    residueVect(fockDim) = -1;
    // std::cout << residueVect << "\n";

    for (std::uint64_t row = 0; row < fockDim; row++)
    {
        for (std::uint64_t col = 0; col <= row; col++)
        {
            std::double_t value = (errorMatrices[row].transpose() * errorMatrices[col]).trace();
            bMatrix(row, col) = value;
            bMatrix(col, row) = value;
        }
    }
    Eigen::VectorXd diisWeights = bMatrix.colPivHouseholderQr().solve(residueVect);
    // std::cout << diisWeights << " " << diisWeights.sum() << "\n";
    // std::cout << bMatrix.rows() << " " << residueVect.size() << "\n";
    Eigen::VectorXd diisModified = diisWeights;
    diisModified(fockDim) = 0;
    assert((diisModified.sum() - 1) < 1e-6);

    fockMatrix.setZero();
    for (std::uint64_t ii = 0; ii < fockDim; ii++)
    {
        fockMatrix = fockMatrix + (fockMatrices[ii] * diisWeights[ii]);
    }
}