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

#include <iomanip>

#include "../lookup/planck_lookup.h"
#include "planck_io.h"

void readInput(std::fstream *filePointer, cxx_Calculator *planckCalculator, cxx_Molecule *inputMolecule, std::error_code *errorFlag, std::string *errorMessage)
{
    // first reset the error buffers
    errorFlag->clear();
    errorMessage->clear();

    // check if the filepointer is valid
    if (!filePointer || !filePointer->is_open())
    {
        *errorFlag = std::make_error_code(std::errc::io_error);
        *errorMessage = "Unable To Open The Input File. Please Check The Input File Provided.";
        return;
    }

    // tokeinze the input file and parse the data
    tokenizeInput(filePointer, planckCalculator, inputMolecule, errorFlag, errorMessage);

    // now convert the coordinates
    if (planckCalculator->coordinate_type == "ang")
    {
        for (std::uint64_t ii = 0; ii < planckCalculator->total_atoms; ii++)
        {
            inputMolecule->input_coordinates[ii * 3 + 0] = inputMolecule->input_coordinates[ii * 3 + 0] * ANGTOBOHR;
            inputMolecule->input_coordinates[ii * 3 + 1] = inputMolecule->input_coordinates[ii * 3 + 1] * ANGTOBOHR;
            inputMolecule->input_coordinates[ii * 3 + 2] = inputMolecule->input_coordinates[ii * 3 + 2] * ANGTOBOHR;
        }
    }

    // now check if charge and spin multiplicity match
    bool checkMultiplicity = false;
    std::int64_t unpairedElectrons = planckCalculator->molecule_multiplicity - 1;

    if (unpairedElectrons > static_cast<int64_t>(planckCalculator->total_electrons))
    {
        *errorFlag = std::make_error_code(std::errc::invalid_argument);
        *errorMessage = "You Cannot Have " + std::to_string(unpairedElectrons) + " Unpaired Electrons With A Total Electron Count Of " + std::to_string(planckCalculator->total_electrons);
        return;
    }

    for (std::int64_t ii = planckCalculator->total_electrons; ii >= 0; ii -= 2)
    {
        // std::cout << ii << "\n";
        if (unpairedElectrons == ii)
        {
            checkMultiplicity = true;
            planckCalculator->alpha_electrons = (planckCalculator->total_electrons / 2) + (unpairedElectrons / 2);
            planckCalculator->beta_electrons  = planckCalculator->total_electrons - planckCalculator->alpha_electrons;

            if (planckCalculator->alpha_electrons < planckCalculator->beta_electrons)
            {
                std::uint64_t _ = planckCalculator->alpha_electrons;
                planckCalculator->alpha_electrons = (planckCalculator->beta_electrons);
                planckCalculator->beta_electrons = _; 
            }
            break;
        }
    }

    if (!checkMultiplicity)
    {
        *errorFlag = std::make_error_code(std::errc::invalid_argument);
        *errorMessage = "A Combination Of " + std::to_string(planckCalculator->total_electrons) + " And A Multiplicty Of " + std::to_string(planckCalculator->molecule_multiplicity) + " Is Not Allowed.";
        return;
    }

    // if the number of alpha and beta electrons are not equal, but the requested theory is restricted
    if ((planckCalculator->alpha_electrons != planckCalculator->beta_electrons) && (planckCalculator->calculation_theory[0] == 'r'))
    {
        *errorFlag = std::make_error_code(std::errc::protocol_error);
        *errorMessage = "A Combination Of " + std::to_string(planckCalculator->total_electrons) + " And A Multiplicty Of " + std::to_string(planckCalculator->molecule_multiplicity) + " Is Not Allowed For Restricted Shell Calculations. The Calculations Will Be Converted To Unrestricted Shell.";
        planckCalculator->calculation_theory.replace(0, 1, "u");
        return;
    }

    // if successful, set the erroflag to zero
    errorFlag = nullptr;
}

void dumpInput(cxx_Calculator *planckCalculator, cxx_Molecule *inputMolecule)
{
    std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Requested Calculation Type : " << planckCalculator->calculation_type << "\n";
    std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Requested Basis Set : " << planckCalculator->calculation_basis << "\n";
    std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Requested Level Of Theory : " << planckCalculator->calculation_theory << "\n";
    std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
    std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
    std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Number Of Atoms : " << planckCalculator->total_atoms << "\n";
    std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Input Charge : " << planckCalculator->molecule_charge << "\n";
    std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Input Multiplicity : " << planckCalculator->molecule_multiplicity << "\n";
    std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Number Of Electrons : " << planckCalculator->total_electrons << "\n";
    std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
    std::cout << std::setw(20) << std::left << "[Planck] " << "\n";

    // now print the SCF control variables
    std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Max SCF Steps : " << planckCalculator->max_scf << "\n";
    std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Max GEOM Iter : " << planckCalculator->max_iter << "\n";
    std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " SCF Tolerance : " << planckCalculator->tol_scf << "\n";
    std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " ERI Tolerance : " << planckCalculator->tol_eri << "\n";
    std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
    std::cout << std::setw(20) << std::left << "[Planck] " << "\n";

    // check if the calculation is unrestricted, if yes print the number of alpha and beta electrons
    if (planckCalculator->calculation_theory[0] == 'u')
    {
        std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Number Of Alpha Electrons : " << planckCalculator->alpha_electrons << "\n";
        std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Number Of Beta Electrons : " << planckCalculator->beta_electrons << "\n";
        std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
        std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
    }
    
    if (inputMolecule->is_reoriented)
    {
        // std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Use Point Group Symmetry : " << inputMolecule->use_pgsymmetry << "\n";
        std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Detected Point Group : " << inputMolecule->point_group << "\n";
        std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
        std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
    }

    // Now dump input coordinates
    std::cout << std::setw(21) << std::left << "[Planck]   => " << std::setw(80) << std::left << "Begin Input Coordinates" << "\n";
    for (std::uint64_t atom_index = 0; atom_index < planckCalculator->total_atoms; atom_index++)
    {
        std::cout << std::setw(20) << std::left << "[Planck]"
                  << std::setw(20) << std::fixed << std::right << inputMolecule->atom_numbers[atom_index]
                  << std::setw(20) << std::fixed << std::right << inputMolecule->input_coordinates[atom_index * 3 + 0]
                  << std::setw(20) << std::fixed << std::right << inputMolecule->input_coordinates[atom_index * 3 + 1]
                  << std::setw(20) << std::fixed << std::right << inputMolecule->input_coordinates[atom_index * 3 + 2]
                  << "\n";
    }
    std::cout << std::setw(21) << std::left << "[Planck]   => " << std::setw(80) << std::left << "End Input Coordinates" << "\n";
    std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
    std::cout << std::setw(20) << std::left << "[Planck] " << "\n";

    // if the point group is not C1, print reoriented coordinates also
    if (inputMolecule->is_reoriented)
    {
        std::cout << std::setw(21) << std::left << "[Planck]   => " << std::setw(80) << std::left << "Begin Standard Coordinates" << "\n";
        for (std::uint64_t atom_index = 0; atom_index < planckCalculator->total_atoms; atom_index++)
        {
            std::cout << std::setw(20) << std::left << "[Planck]"
                      << std::setw(20) << std::fixed << std::right << inputMolecule->atom_numbers[atom_index]
                      << std::setw(20) << std::fixed << std::right << inputMolecule->standard_coordinates[atom_index * 3 + 0]
                      << std::setw(20) << std::fixed << std::right << inputMolecule->standard_coordinates[atom_index * 3 + 1]
                      << std::setw(20) << std::fixed << std::right << inputMolecule->standard_coordinates[atom_index * 3 + 2]
                      << "\n";
        }
        std::cout << std::setw(21) << std::left << "[Planck]   => " << std::setw(80) << std::left << "End Standard Coordinates" << "\n";
        std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
        std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
    }

    // now print basis set information
    std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Number Of Basis Functions : " << planckCalculator->total_basis << "\n";
    std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Number Of Primtive Gaussians : " << planckCalculator->total_primitives << "\n";
    std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
    std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
}

// TODO : write a better way to store SCF information on file
// void dumpIntegral(cxx_Calculator *plankCalculator, cxx_Integrals integralType, std::string inputFile)
// {
//     // housekeeping
//     std::string fileName = inputFile + ".chk";
//     std::fstream filePointer(fileName, std::ios::out | std::ios::binary | std::ios::trunc);

//     // check if file is valid
//     if (!filePointer || !filePointer.is_open())
//     {
//         std::cout << "unable to open file";
//     }
    
//     // write the size of the integral
//     filePointer.write(reinterpret_cast<char*>(), sizeof(std::uint64_t));

//     // // write the integrals
//     filePointer.write(reinterpret_cast<char*>(integral), sizeof(std::double_t) * dim);

//     // close the file pointer after writing
//     filePointer.close();
// }