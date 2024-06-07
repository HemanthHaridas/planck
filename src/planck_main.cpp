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

int main(int argc, char const *argv[])
{
    cxx_Calculator planck_calculator;
    cxx_Molecule input_molecule;

    std::error_code error_flag;
    std::string error_message;

    // start of the program
    boost::posix_time::ptime startTime = boost::posix_time::microsec_clock::local_time();
    std::cout << std::setw(20) << std::left << "[Planck] => " << std::setw(35) << std::left << " Program Started On : " << startTime << "\n";
    std::cout << std::setw(20) << std::left << "[Planck] => " << std::setw(35) << std::left << " Current Working Directory : " << std::filesystem::current_path().string() << "\n";
    std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
    std::cout << std::setw(20) << std::left << "[Planck] " << "\n";

    // check if an input file is provided. Exit if not input file is provided
    if (argc < 2)
    {
        std::cout << std::setw(20) << std::left << "[Error]  <= " << std::left << " Unable To Find An Input File. Please Run Planck As : planck input " << "\n";
        exit(-1);
    }

    // proceed to parsing to the input file
    std::string input_file = argv[1];
    std::fstream file_pointer(input_file);

    readInput(&file_pointer, &planck_calculator, &input_molecule, &error_flag, &error_message);

    // check if input file was parsed correctly
    if (error_flag && error_flag.value() != 71)
    {
        std::cout << std::setw(21) << std::left << "[Error]  <=  " << std::left << error_message << "\n";
        exit(error_flag.value());
    }

    // need a special case to handle rhf -> uhf error
    if (error_flag.value() == 71)
    {
        std::cout << std::setw(21) << std::left << "[Warning] <= " << std::setw(35) << std::left << error_message << "\n";
        std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
        std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
    }

    // now symmetrize the molecule and process the errors
    detectSymmetry(&input_molecule, planck_calculator.total_atoms, &error_flag, &error_message);

    if (error_flag)
    {
        std::cout << std::setw(21) << std::left << "[Error]  <=  " << std::left << error_message << "\n";
        exit(error_flag.value());
    }

    // now read the basis sets
    readBasis(&input_molecule, &planck_calculator, &error_flag, &error_message);

    if (error_flag)
    {
        std::cout << std::setw(21) << std::left << "[Error]  <=  " << std::left << error_message << "\n";
        exit(error_flag.value());
    }
    
    // start dumping the input file
    dumpInput(&planck_calculator, &input_molecule);

    // end of the program
    boost::posix_time::ptime endTime = boost::posix_time::microsec_clock::local_time();
    std::cout << std::setw(20) << std::left << "[Planck] => " << std::setw(35) << std::left << " Program Completed On : " << endTime << "\n";
    return 0;
}
