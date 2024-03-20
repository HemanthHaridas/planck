#pragma once

#include <boost/algorithm/string/erase.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <cmath>
#include <cstdint>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

#include "../auxiliary/tables.h"

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

struct cxx_Primitives
{
    std::double_t orbitalCoeff;
    std::double_t primitiveExp;
    std::double_t orbitalNorms;
    std::double_t locationX;
    std::double_t locationY;
    std::double_t locationZ;

    std::int64_t angularMomentumX;
    std::int64_t angularMomentumY;
    std::int64_t angularMomentumZ;

    std::uint64_t index;
};

struct cxx_gptResults
{
    std::double_t locationX;
    std::double_t locationY;
    std::double_t locationZ;

    std::double_t integralX;
    std::double_t integralY;
    std::double_t integralZ;

    std::uint64_t indexA;
    std::uint64_t indexB;
};

struct cxx_Integral
{
    std::uint64_t indexA;
    std::uint64_t indexB;
    std::double_t result = 0;
};

// Moved the results to a separate block
// Everything except replusion integrals in a 2D matrix.
struct cxx_Results {
    Eigen::Matrix<cxx_gptResults, Eigen::Dynamic, Eigen::Dynamic> gaussianResults;
    Eigen::Matrix<std::double_t, Eigen::Dynamic, Eigen::Dynamic> overlapIntegrals;
    Eigen::Matrix<std::double_t, Eigen::Dynamic, Eigen::Dynamic> kineticIntegrals;
    Eigen::Matrix<std::double_t, Eigen::Dynamic, Eigen::Dynamic> nuclearIntegrals;
    Eigen::Tensor<std::double_t, 4> repulsionIntegrals;
};

struct cxx_Basis
{
    std::vector<cxx_Primitives> cGTO;
    std::uint64_t index;
};

struct cxx_Molecule
{
    Eigen::Matrix<double, Eigen::Dynamic, 3> atomCoordinates;
    Eigen::VectorXd atomMasses;
    Eigen::VectorXi atomNumbers;
    std::vector<std::uint16_t> atomBasis;
};

struct cxx_Calculator
{
    std::int64_t totalCharge;
    std::string calculationBasis;
    std::string calculationTheory;
    std::string calculationType;
    std::uint16_t molMultiplicity;

    // This block is required only if the theory is UHF
    std::uint64_t alphaElectrons;
    std::uint64_t betaElectrons;

    std::uint64_t nAtoms;
    std::uint64_t nBasis = 0;
    std::uint64_t nPrimitives = 0;
    std::vector<cxx_Basis> basisFunctions;
    cxx_Results resultSCF;
};

void readInput(std::fstream *filePointer, cxx_Molecule *inputMolecule, cxx_Calculator *scfCalculator, std::error_code *errorFlag, std::string *errorMessage);
void readBasis(std::fstream *basisPointer, std::string atomNumber, std::string atomIndex, std::error_code *errorFlag, std::string *errorMessage);