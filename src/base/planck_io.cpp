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

    // first read the header section
    std::string headerLine;
    std::getline(*filePointer, headerLine);
    std::stringstream headerBuffer(headerLine);

    // check for optional symmetry keyword
    headerBuffer >> planckCalculator->calculation_type >> planckCalculator->calculation_theory >> planckCalculator->calculation_basis >> planckCalculator->coordinate_type >> inputMolecule->use_pgsymmetry;

    // always use point group symmetry by default
    if (inputMolecule->use_pgsymmetry != 0)
    {
        inputMolecule->use_pgsymmetry = 1;
    }

    // check if defaults.txt file is present
    // if yes, read in the default basis set path
    // else read path from input file
    std::fstream defaultPointer("planck.defaults");
    if (defaultPointer)
    {
        std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Found planck.defaults file" << "\n";
        std::getline(*filePointer, planckCalculator->basis_path);  // will be immediately overwritten by path read from planck.defaults file
        std::getline(defaultPointer, planckCalculator->basis_path);
    }
    else
    {
        std::getline(*filePointer, planckCalculator->basis_path);
    }

    std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Basis sets read from : " << planckCalculator->basis_path << "\n";

    // now read the number of atoms and set up the buffers
    std::getline(*filePointer, headerLine);
    std::stringstream atomBuffer(headerLine);
    atomBuffer >> planckCalculator->total_atoms;
    // std::cout << planckCalculator->total_atoms << "\n";

    inputMolecule->atom_masses = (std::double_t *)malloc(sizeof(std::double_t) * planckCalculator->total_atoms);
    inputMolecule->atom_numbers = (std::uint64_t *)malloc(sizeof(std::uint64_t) * planckCalculator->total_atoms);
    inputMolecule->input_coordinates = (std::double_t *)malloc(sizeof(std::double_t) * planckCalculator->total_atoms * 3);

    std::getline(*filePointer, headerLine);
    std::stringstream infoBuffer(headerLine);
    infoBuffer >> planckCalculator->molecule_charge >> planckCalculator->molecule_multiplicity;

    // check if multiplicity is a positive number >= 1
    if (planckCalculator->molecule_multiplicity < 1)
    {
        *errorFlag = std::make_error_code(std::errc::invalid_argument);
        *errorMessage = "Multiplicity Cannot Be A Value That Is Less Than 1. Please Check Your Input File.";
        return;
    }

    std::uint64_t atomIndex = 0;
    planckCalculator->total_electrons = 0;
    
    while (std::getline(*filePointer, headerLine) && atomIndex < planckCalculator->total_atoms)
    {
        // buffers to hold the atom information
        std::string atomName;
        std::double_t xCoord;
        std::double_t yCoord;
        std::double_t zCoord;

        std::stringstream coordBuffer(headerLine);
        coordBuffer >> atomName >> xCoord >> yCoord >> zCoord;

        inputMolecule->atom_masses[atomIndex] = atomicMass[atomName];
        inputMolecule->atom_numbers[atomIndex] = atomicNumber[atomName];
        inputMolecule->input_coordinates[atomIndex * 3 + 0] = xCoord;
        inputMolecule->input_coordinates[atomIndex * 3 + 1] = yCoord;
        inputMolecule->input_coordinates[atomIndex * 3 + 2] = zCoord;

        planckCalculator->total_electrons += atomicNumber[atomName];
        atomIndex++;
    }
    planckCalculator->total_electrons = planckCalculator->total_electrons + planckCalculator->molecule_charge;

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

    if (unpairedElectrons > static_cast<int64_t>(planckCalculator->total_electrons / 2 + 1))
    {
        *errorFlag = std::make_error_code(std::errc::invalid_argument);
        *errorMessage = "You Cannot Have " + std::to_string(unpairedElectrons) + " Extra Unpaired Electrons With A Total Electron Count Of " + std::to_string(planckCalculator->total_electrons);
        return;
    }

    for (std::int64_t ii = planckCalculator->total_electrons; ii >= 0; ii -= 2)
    {
        // std::cout << ii << "\n";
        if (unpairedElectrons == ii)
        {
            checkMultiplicity = true;
            planckCalculator->alpha_electrons = (planckCalculator->total_electrons / 2) + unpairedElectrons;
            planckCalculator->beta_electrons = planckCalculator->total_electrons - planckCalculator->alpha_electrons;
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
        std::cout << std::setw(20) << std::left << "[Planck]   => " << std::setw(35) << std::left << " Use Point Group Symmetry : " << inputMolecule->use_pgsymmetry << "\n";
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