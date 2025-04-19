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

#include <boost/date_time/posix_time/posix_time.hpp>
#include <filesystem>
#include <iomanip>

#include "base/planck_base.h"
#include "base/planck_basis.h"
#include "base/planck_io.h"
#include "base/planck_symmetry.h"
#include "math/planck_math.h"

int main(int argc, char const *argv[])
{
    cxx_Calculator planck_calculator;
    cxx_Molecule input_molecule;

    std::error_code error_flag;
    std::string error_message;

    // start of the program
    boost::posix_time::ptime startTime = boost::posix_time::microsec_clock::local_time();
    std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Program Started On : " << startTime << "\n";
    std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Current Working Directory : " << std::filesystem::current_path().string() << "\n";
    std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
    std::cout << std::setw(20) << std::left << "[Planck] " << "\n";

    // check if an input file is provided. Exit if not input file is provided
    if (argc < 2)
    {
        std::cout << std::setw(20) << std::left << "[Error]    <= " << std::left << " Unable To Find An Input File. Please Run Planck As : planck input " << "\n";
        exit(-1);
    }

    // proceed to parsing to the input file
    std::string input_file = argv[1];
    std::fstream file_pointer(input_file);

    // first check if a planck.defaults file is available
    std::fstream defaultPointer("planck.defaults");
    if (defaultPointer)
    {
        std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Found planck.defaults file" << "\n";
        std::getline(defaultPointer, planck_calculator.basis_path);
    }

    // tokenizeInput(&file_pointer, &planck_calculator, &input_molecule, &error_flag, &error_message);
    readInput(&file_pointer, &planck_calculator, &input_molecule, &error_flag, &error_message);

    // check if input file was parsed correctly
    if (error_flag && error_flag.value() != std::make_error_code(std::errc::protocol_error).value())
    {
        std::cout << std::setw(21) << std::left << "[Error]    <= " << std::left << error_message << "\n";
        exit(error_flag.value());
    }

    // need a special case to handle rhf -> uhf error
    if (error_flag.value() == std::make_error_code(std::errc::protocol_error).value())
    {
        std::cout << std::setw(21) << std::left << "[Warning]  <= " << std::setw(35) << std::left << error_message << "\n";
        std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
        std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
    }

    if (planck_calculator.use_pgsymmetry)
    {
        // now symmetrize the molecule and process the errors
        detectSymmetry(&input_molecule, planck_calculator.total_atoms, &error_flag, &error_message);
        if (error_flag)
        {
            std::cout << std::setw(21) << std::left << "[Error]   <=  " << std::left << error_message << "\n";
            exit(error_flag.value());
        }
        // std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
        // std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
    }
    else
    {
        std::cout << std::setw(21) << std::left << "[Warning]  <= " << std::setw(35) << std::left << "Symmetry detection is turned off by request" << "\n";
        input_molecule.is_reoriented = false;
        input_molecule.standard_coordinates = input_molecule.input_coordinates;
        std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
        std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
    }

    // now read the basis sets
    readBasis(&input_molecule, &planck_calculator, &error_flag, &error_message);

    std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
    std::cout << std::setw(20) << std::left << "[Planck] " << "\n";

    if (error_flag)
    {
        std::cout << std::setw(21) << std::left << "[Error]   <=  " << std::left << error_message << "\n";
        exit(error_flag.value());
    }

    // start dumping the input file
    dumpInput(&planck_calculator, &input_molecule);

    // preallocate buffers
    planck_calculator.overlapMatrix.resize(planck_calculator.total_basis, planck_calculator.total_basis);
    planck_calculator.kineticMatrix.resize(planck_calculator.total_basis, planck_calculator.total_basis);
    planck_calculator.nuclearMatrix.resize(planck_calculator.total_basis, planck_calculator.total_basis);
    planck_calculator.electronicMatrix.resize(static_cast<int64_t>(planck_calculator.total_basis), static_cast<int64_t>(planck_calculator.total_basis), static_cast<int64_t>(planck_calculator.total_basis), static_cast<int64_t>(planck_calculator.total_basis));

    // std::uint64_t totalMemory;

    // preallocate the buffers
    // planck_calculator.overlap = (std::double_t *)malloc(sizeof(std::double_t) * planck_calculator.total_basis * planck_calculator.total_basis);
    // planck_calculator.kinetic = (std::double_t *)malloc(sizeof(std::double_t) * planck_calculator.total_basis * planck_calculator.total_basis);
    // planck_calculator.nuclear = (std::double_t *)malloc(sizeof(std::double_t) * planck_calculator.total_basis * planck_calculator.total_basis);
    // planck_calculator.electronic = (std::double_t *)malloc(sizeof(std::double_t) * planck_calculator.total_basis * planck_calculator.total_basis);

    // memset all the buffers to zero to avoid junk values
    // memset(planck_calculator.overlap,    0, sizeof(std::double_t) * planck_calculator.total_basis * planck_calculator.total_basis);
    // memset(planck_calculator.kinetic,    0, sizeof(std::double_t) * planck_calculator.total_basis * planck_calculator.total_basis);
    // memset(planck_calculator.nuclear,    0, sizeof(std::double_t) * planck_calculator.total_basis * planck_calculator.total_basis);
    // memset(planck_calculator.electronic, 0, sizeof(std::double_t) * planck_calculator.total_basis * planck_calculator.total_basis);

    // if theory is uhf => allocate twice big size for fock matrix
    if (planck_calculator.is_unrestricted)
    {
        std::uint64_t nelem = (4 * planck_calculator.total_basis * planck_calculator.total_basis);
        planck_calculator.fockMatrix.resize(nelem, nelem);
    }
    else
    {
        std::uint64_t nelem = (planck_calculator.total_basis * planck_calculator.total_basis);
        planck_calculator.fockMatrix.resize(nelem, nelem);
    }

    // dump integrals
    // dumpIntegral(planck_calculator.overlap, planck_calculator.total_basis * planck_calculator.total_basis, "overlap", input_file);
    // dumpIntegral(planck_calculator.kinetic, planck_calculator.total_basis * planck_calculator.total_basis, "kinetic", input_file);

    // free the allocated buffers
    // free(planck_calculator.overlap);
    // free(planck_calculator.kinetic);
    // free(planck_calculator.nuclear);
    // free(planck_calculator.electronic);
    // free(planck_calculator.fock);

    free(input_molecule.input_coordinates);
    free(input_molecule.standard_coordinates);
    free(input_molecule.atom_masses);
    free(input_molecule.atom_numbers);

    // end of the program
    boost::posix_time::ptime endTime = boost::posix_time::microsec_clock::local_time();
    std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Program Completed On : " << endTime << "\n";
    return 0;
}
