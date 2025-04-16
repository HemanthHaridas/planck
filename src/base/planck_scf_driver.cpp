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

#include "planck_scf_driver.h"

// void scfEngine(cxx_Molecule *inputMolecule, cxx_Calculator *planckCalculator, std::error_code *errorFlag, std::string *errorMessage)
// {
//     // first compute overlap and kinetic integrals
//     computeOverlap(planckCalculator, errorFlag, errorMessage);
//     computeKinetic(planckCalculator, errorFlag, errorMessage);

//     // generate core hamiltonian

//     // now compute electron nuclear integrals
//     computeNuclear(inputMolecule->standard_coordinates, inputMolecule->atom_numbers, planckCalculator->total_atoms, planckCalculator, errorFlag, errorMessage);

//     // now compute electron electron integrals
//     computeElectronic(planckCalculator, errorFlag, errorMessage);
// }

std::int64_t scfEngine::scfCycle()
{
    // increment SCF cycle count
    this->cycle++;

    // first clear the error buffers
    this->errorFlag.clear();
    this->errorMessage.clear();

    if (EXIT_SUCCESS != scfEngine::computeOverlap())
        return this->errorFlag.value();
    if (EXIT_SUCCESS != scfEngine::computeKinetic())
        return this->errorFlag.value();
    if (EXIT_SUCCESS != scfEngine::computeNuclear())
        return this->errorFlag.value();
    if (EXIT_SUCCESS != scfEngine::computeElectronic())
        return this->errorFlag.value();
}

std::int64_t scfEngine::computeOverlap()
{
    std::int64_t basisSize = this->total_basis;

    // now iterate over the contracted basis functions
    for (std::uint64_t ii = 0; ii < basisSize; ii++)
    {
        for (std::uint64_t jj = 0; jj < basisSize; jj++)
        {
            this->overlap[ii * basisSize + jj] = Huzinaga::computeOverlap(&this->calculation_set[ii], &this->calculation_set[jj], &this->errorFlag, &this->errorMessage);
        }
    }

    // check if the error code is zero
    if (errorFlag.value() != 0)
    {
        return errorFlag.value();
    }

    // return success
    return EXIT_SUCCESS;
}

std::int64_t scfEngine::computeKinetic()
{
    std::int64_t basisSize = this->total_basis;

    // now iterate over the contracted basis functions
    for (std::uint64_t ii = 0; ii < basisSize; ii++)
    {
        for (std::uint64_t jj = 0; jj < basisSize; jj++)
        {
            this->kinetic[ii * basisSize + jj] = Huzinaga::computeKinetic(&this->calculation_set[ii], &this->calculation_set[jj], &this->errorFlag, &this->errorMessage);
        }
    }

    // check if the error code is zero
    if (errorFlag.value() != 0)
    {
        return errorFlag.value();
    }

    // return success
    return EXIT_SUCCESS;
}

std::int64_t scfEngine::computeNuclear()
{
    std::int64_t basisSize = this->total_basis;

    // now iterate over the contracted basis functions
    for (std::uint64_t ii = 0; ii < basisSize; ii++)
    {
        for (std::uint64_t jj = 0; jj < basisSize; jj++)
        {
            this->nuclear[ii * basisSize + jj] = Huzinaga::computeNuclear(this->standard_coordinates, this->atom_numbers, this->total_atoms, &this->calculation_set[ii], &this->calculation_set[jj], &this->errorFlag, &this->errorMessage);
        }
    }

    // check if the error code is zero
    if (errorFlag.value() != 0)
    {
        return errorFlag.value();
    }

    // return success
    return EXIT_SUCCESS;
}

std::int64_t scfEngine::schwartzScreening()
{
    std::int64_t basisSize = this->total_basis;

    // now iterate over the contracted basis functions
    for (std::uint64_t ii = 0; ii < basisSize; ii++)
    {
        for (std::uint64_t jj = 0; jj < basisSize; jj++)
        {
            this->electronic[ii * basisSize + jj * basisSize + ii * basisSize + jj] = Huzinaga::computeElectronic(
                &this->calculation_set[ii], &this->calculation_set[jj], &this->calculation_set[ii], &this->calculation_set[jj], &this->errorFlag, &this->errorMessage);
        }
    }
}

std::int64_t scfEngine::computeElectronic()
{
    std::int64_t basisSize = this->total_basis;

    // now iterate over the contracted basis functions
    for (std::uint64_t ii = 0; ii < basisSize; ii++)
    {
        for (std::uint64_t jj = 0; jj < basisSize; jj++)
        {
            std::double_t bra = this->electronic[ii * basisSize + jj * basisSize + ii * basisSize + jj];
            for (std::uint64_t kk = 0; kk < basisSize; kk++)
            {
                for (std::uint64_t ll = 0; ll < basisSize; ll++)
                {
                    std::double_t ket = this->electronic[kk * basisSize + ll * basisSize + kk * basisSize + ll];
                    if (((ii * jj) >= (kk * ll)) && (sqrt(bra * ket) <= this->tol_eri))
                    {
                        this->electronic[ii * basisSize + jj * basisSize + kk * basisSize + ll] = Huzinaga::computeElectronic(
                            &this->calculation_set[ii], &this->calculation_set[jj], &this->calculation_set[ii], &this->calculation_set[jj], &this->errorFlag, &this->errorMessage);
                    }
                }
            }
        }
    }

    // check if the error code is zero
    if (errorFlag.value() != 0)
    {
        return errorFlag.value();
    }

    // return success
    return EXIT_SUCCESS;
}