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

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

// #include <numbers>

const std::double_t ANGTOBOHR = 1.8897259886;
const std::uint64_t MAXSCF = 120;
const std::uint64_t MAXITER = 120;
const std::double_t TOLSCF = 1.0E-14;
const std::double_t TOLERI = 1.0E-14;
const std::string DEFAULT_BASIS = "sto-3g";
const std::string DEFAULT_THEORY = "rhf";
const std::string DEFAULT_CALC = "energy";
const std::string DEFAULT_COORD = "ang";
const bool USE_DIIS = true;
const bool USE_SYMM = true;

// this is the maximum supported boysindex
// change this value only if you regenerate
// the lookup table for the boysfunction
const std::uint64_t MAXM = 60;

// const std::double_t PI=355/113;

struct cxx_Primitive
{
    // definition of primitive gaussians
    std::double_t primitive_exp;
    std::double_t orbital_coeff;
    std::double_t orbital_norm;

    // index of primitive gaussian
    std::uint64_t pgto_index;
};

struct cxx_Contracted
{
    // location of contracted gaussian
    std::double_t location_x;
    std::double_t location_y;
    std::double_t location_z;

    // angular momenetum of contracted gaussian
    std::int64_t shell_x;
    std::int64_t shell_y;
    std::int64_t shell_z;

    // index of contracted gaussian
    std::uint64_t cgto_index;

    // list of primitive gaussians forming contracted gaussian
    std::vector<cxx_Primitive> contracted_GTO;
};

struct cxx_Molecule
{
    // buffers to hold molecule geometry and associated data
    std::double_t *input_coordinates;
    std::double_t *standard_coordinates;
    std::uint64_t *atom_numbers;
    std::double_t *atom_masses;

    // molecule point group
    std::string point_group;
    bool is_reoriented;
};

struct cxx_Calculator
{
    // variables to hold molecule information
    std::int64_t molecule_charge;
    std::uint64_t alpha_electrons;
    std::uint64_t beta_electrons;
    std::uint64_t molecule_multiplicity;
    std::uint64_t total_electrons;
    std::uint64_t total_atoms;
    bool is_unrestricted;
    bool use_pgsymmetry = USE_SYMM;

    // variables to hold calculation interface
    std::double_t tol_eri = TOLERI;
    std::double_t tol_scf = TOLSCF;

    std::string calculation_basis = DEFAULT_BASIS;
    std::string calculation_theory = DEFAULT_BASIS;
    std::string calculation_type = DEFAULT_CALC;
    std::string coordinate_type = DEFAULT_COORD;
    std::uint64_t max_iter = MAXITER;
    std::uint64_t max_scf = MAXSCF;

    bool use_diis = USE_DIIS;

    std::string basis_path;

    // array to hold the basis set
    std::vector<cxx_Contracted> calculation_set;
    std::uint64_t total_basis;
    std::uint64_t total_primitives;
};

struct cxx_Integrals
{
    // Eigen matrices to hold integrals;
    Eigen::MatrixXd overlapMatrix;
    Eigen::MatrixXd kineticMatrix;
    Eigen::MatrixXd nuclearMatrix;
    Eigen::Tensor<std::double_t, 4> electronicMatrix;
};

struct scfData
{
    // Eigen matrices for intermediates
    Eigen::MatrixXd fockMatrix;
    Eigen::MatrixXd hamiltonianMatrix;
    Eigen::MatrixXd orthoMatrix;
    Eigen::MatrixXd coreMatrix;

    // Eigen matrices for MOs
    Eigen::MatrixXd canonicalMO;
    Eigen::MatrixXd orthogonalMO;
};

struct cxx_Gaussians
{
    std::double_t gaussian_center[3];
    std::double_t gaussian_exponent;
    std::double_t gaussian_integral[4];
};

struct nuclearInt
{
    std::double_t result;
    std::uint64_t x;
    std::uint64_t y;
    std::uint64_t z;
};

struct electronInt
{
    std::double_t result;
    std::uint64_t i;
    std::uint64_t j;
    std::uint64_t k;
    std::uint64_t l;
    std::uint64_t m;
};

typedef std::tuple<std::uint64_t, std::uint64_t, std::uint64_t, std::uint64_t> eriShell;
typedef std::tuple<std::uint64_t, std::uint64_t> eriKet;