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
#include <iomanip>

#include "../base/planck_base.h"
#include "planck_integrals.h"
#include "./helper/planck_helper_routines.h"

void computeOverlap(cxx_Calculator *planckCalculator, Eigen::MatrixXd &overlapMatrix)
{
    // get the number of basis functions
    std::uint64_t nBasis = planckCalculator->total_basis;
    // overlapMatrix.resize(nBasis, nBasis);

    // now iterate over the contracted GTOs
    for (std::uint64_t row = 0; row < nBasis; row++)
    {
        for (std::uint64_t col = 0; col <= row; col++)
        {
            overlapMatrix(row, col) = Huzinaga::computeOverlap(&planckCalculator->calculation_set[row], &planckCalculator->calculation_set[col]);
            overlapMatrix(col, row) = overlapMatrix(row, col);
        }
    }
}

void computeKinetic(cxx_Calculator *planckCalculator, Eigen::MatrixXd &kineticMatrix)
{
    // get the number of basis functions
    std::uint64_t nBasis = planckCalculator->total_basis;
    // kineticMatrix.resize(nBasis, nBasis);

    // now iterate over the contracted GTOs
    for (std::uint64_t row = 0; row < nBasis; row++)
    {
        for (std::uint64_t col = 0; col <= row; col++)
        {
            kineticMatrix(row, col) = Huzinaga::computeKinetic(&planckCalculator->calculation_set[row], &planckCalculator->calculation_set[col]);
            kineticMatrix(col, row) = kineticMatrix(row, col);
        }
    }
}

void computeNuclear(std::double_t *atomCoords, std::uint64_t *atomCharges, cxx_Calculator *planckCalculator, Eigen::MatrixXd &nuclearMatrix)
{
    // get the number of basis functions
    std::uint64_t nBasis = planckCalculator->total_basis;
    // nuclearMatrix.resize(nBasis, nBasis);

    // now iterate over the contracted GTOs
    for (std::uint64_t row = 0; row < nBasis; row++)
    {
        for (std::uint64_t col = 0; col <= row; col++)
        {
            nuclearMatrix(row, col) = Huzinaga::computeNuclear(atomCoords, atomCharges, planckCalculator->total_atoms, &planckCalculator->calculation_set[row], &planckCalculator->calculation_set[col]);
            nuclearMatrix(col, row) = nuclearMatrix(row, col);
        }
    }
}

std::vector<eriShell> schwartzScreeing(cxx_Calculator *planckCalculator, Eigen::Tensor<std::double_t, 4> &electronMatrix)
{
    std::vector<eriShell> screenedPairs;

    // first get the screened kets
    std::vector<eriKet> eriScreened = schwatrzSceening(planckCalculator, electronMatrix);
    // std::cout << eriScreened.size() << "\n";

    // now generate all possible pairs of (ab|cd)
    for (auto bra : eriScreened)
    {
        std::uint64_t b = std::get<0>(bra) * std::get<1>(bra);
        for (auto ket : eriScreened)
        {
            std::uint64_t k = std::get<0>(ket) * std::get<1>(ket);
            // check if the lexicographic ordering is maintained
            if ((b <= k) && (std::get<0>(bra) <= std::get<1>(bra)) && (std::get<0>(ket) <= std::get<1>(ket)))
            {
                eriShell computeShells(std::get<0>(bra), std::get<1>(bra), std::get<0>(ket), std::get<1>(ket));
                screenedPairs.push_back(computeShells);
            }
        }
    }
    return screenedPairs;
}

void computeElectronic(cxx_Calculator *planckCalculator, Eigen::Tensor<std::double_t, 4> &electronMatrix)
{
    // first do schwartz screening and get initial values
    std::vector<eriShell> screenedPairs = schwartzScreeing(planckCalculator, electronMatrix);

    // print number of screened ERIs
    std::uint64_t nERI = screenedPairs.size(); // This is the naive number of ERIs to be computed
    std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Screened Number of Electron Repulsion Integrals : " << nERI << "\n";

    // now do the full calculations
    for (const auto &braket : screenedPairs)
    {
        std::int64_t ii = std::get<0>(braket);
        std::int64_t jj = std::get<1>(braket);
        std::int64_t kk = std::get<2>(braket);
        std::int64_t ll = std::get<3>(braket);
        
        std::double_t value = Huzinaga::computeElectronic(&planckCalculator->calculation_set[ii], &planckCalculator->calculation_set[jj], &planckCalculator->calculation_set[kk], &planckCalculator->calculation_set[ll]);

        // Assign the value to all configurations in the 8-fold symmetry
        electronMatrix(ii, jj, kk, ll) = value;
        electronMatrix(jj, ii, kk, ll) = value;
        electronMatrix(ii, jj, ll, kk) = value;
        electronMatrix(jj, ii, ll, kk) = value;

        electronMatrix(kk, ll, ii, jj) = value;
        electronMatrix(ll, kk, ii, jj) = value;
        electronMatrix(kk, ll, jj, ii) = value;
        electronMatrix(ll, kk, jj, ii) = value;
    }
}