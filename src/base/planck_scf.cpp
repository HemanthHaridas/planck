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

void scfStep(cxx_scfStep *scfStep, cxx_Calculator *planckCalculator, std::error_code *errorFlag, std::string *errorMessage)
{
    // first reset the flags
    errorFlag->clear();
    errorMessage->clear();

    // first check if scfStep is within limits
    if (scfStep->scfStep > planckCalculator->max_scf)
    {
        *errorFlag = std::make_error_code(std::errc::argument_out_of_domain);
        *errorMessage = "Maximum Number of SCF steps reached. SCF failed to converge within MAX SCF cycles. Please check the results carefully.";
        return;
    }

    // do this only once
    if (scfStep->init_scf)
    {
        // generate core hamiltonian
        // currently uses the sum of kinetic and electron nuclear integrals
        // should implement Harris functional later for better results
        scfStep->coreMatrix = scfStep->kineticMatrix + scfStep->nuclearMatrix;

        // do cholesky decomposition to orthogonalize MOs
        Eigen::LLT<Eigen::MatrixXd> lltSolver(scfStep->overlapMatrix);
        scfStep->orthoMatrix = lltSolver.matrixL();
        scfStep->orthoMatrix = scfStep->orthoMatrix.inverse();

        // verify if the transformation is correct
        Eigen::MatrixXd LHS = scfStep->orthoMatrix.transpose() * scfStep->overlapMatrix * scfStep->orthoMatrix;
        Eigen::MatrixXd RHS = Eigen::MatrixXd::Identity(planckCalculator->total_basis, planckCalculator->total_basis);
        assert((LHS - RHS).norm() < 1e-6); // aseert that the difference must be very close to zero

        // reset the flag
        scfStep->init_scf = false;
    }

    // orthogonalize MOs
    scfStep->orthogonalMO = scfStep->orthoMatrix * scfStep->canonicalMO;
}

Eigen::MatrixXd rhfStep()
{

}

Eigen::MatrixXd ufhStep()
{
    
}