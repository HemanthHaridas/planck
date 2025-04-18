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

void scfStep(std::uint64_t *scfStep, cxx_Calculator *planckCalculator, std::error_code *errorFlag, std::string *errorMessage)
{
    // first reset the flags
    errorFlag->clear();
    errorMessage->clear();

    // first check if scfStep is within limits
    if (*scfStep > planckCalculator->max_scf)
    {
        *errorFlag = std::make_error_code(std::errc::argument_out_of_domain);
        *errorMessage = "Maximum Number of SCF steps reached. SCF failed to converge within MAX SCF cycles. Please check the results carefully.";
        return;
    }

    // generate core hamiltonian
    // currently uses the sum of kinetic and electron nuclear integrals
    // should implement Harris functional later for better results
    planckCalculator->coreMatrix = planckCalculator->kineticMatrix + planckCalculator->nuclearMatrix;
    
    // orthogonalization procedure
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(planckCalculator->coreMatrix);
    Eigen::MatrixXd halfMatrix;
    Eigen::MatrixXd rightMatrix;

    halfMatrix.resize(planckCalculator->total_basis, planckCalculator->total_basis);
    rightMatrix.resize(planckCalculator->total_basis, planckCalculator->total_basis);
    planckCalculator->orthoMatrix.resize(planckCalculator->total_basis, planckCalculator->total_basis);

    // compute eigen values and eigen vectors of overlap matrix
    Eigen::VectorXd eigenValues = eigenSolver.eigenvalues();
    Eigen::VectorXd eigenVectors = eigenSolver.eigenvectors();

    // DIIS is enabled by default
    if (!planckCalculator->use_diis)
    {
        // now do a unitary transformation of the eigen vectors matrix
        halfMatrix  = eigenVectors.array().inverse().sqrt().matrix().asDiagonal();
        rightMatrix = eigenVectors * halfMatrix;

        planckCalculator->orthoMatrix = rightMatrix * eigenVectors.transpose();
    }
}