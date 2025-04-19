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

#include <algorithm>

#include "planck_tokenizer.h"
#include "../lookup/planck_lookup.h"

std::string toLower(const std::string &parsedString)
{
    std::string lowerString = parsedString; // Create a copy to preserve the original string
    std::transform(lowerString.begin(), lowerString.end(), lowerString.begin(), ::tolower);
    return lowerString; // Return the transformed string
}

bool stringToBool(const std::string &parsedString)
{
    std::string upperStr = parsedString;
    std::transform(upperStr.begin(), upperStr.end(), upperStr.begin(), ::toupper); // Convert to uppercase

    if (upperStr == "ON")
    {
        return true; // "ON" maps to true
    }
    else if (upperStr == "OFF")
    {
        return false; // "OFF" maps to false
    }
    else
    {
        throw std::invalid_argument("Invalid string for boolean conversion."); // Handle unexpected input
    }
}

void tokenizeInput(std::fstream *filePointer, cxx_Calculator *planckCalculator, cxx_Molecule *inputMolecule, std::error_code *errorFlag, std::string *errorMessage)
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

    std::string inputLine;
    std::unordered_map<std::string, std::function<void(std::string)>> handlers_setup = {
        {"CALC_TYPE",   [planckCalculator](std::string value){ planckCalculator->calculation_type   = toLower(value); }},
        {"THEORY",      [planckCalculator](std::string value){ planckCalculator->calculation_theory = toLower(value); }},
        {"BASIS",       [planckCalculator](std::string value){ planckCalculator->calculation_basis  = toLower(value); }},
        {"COOR_TYPE",   [planckCalculator](std::string value){ planckCalculator->coordinate_type    = toLower(value); }},
        {"USE_SYMM",    [planckCalculator](std::string value){ planckCalculator->use_pgsymmetry     = stringToBool(value); }},
        {"USE_DIIS",    [planckCalculator](std::string value){ planckCalculator->use_diis           = stringToBool(value); }},
        {"BASIS_PATH",  [planckCalculator](std::string value){ planckCalculator->basis_path         = toLower(value);}}};

    std::unordered_map<std::string, std::function<void(std::string)>> handlers_control = {
        {"MAXITER", [planckCalculator](std::string value){  planckCalculator->max_iter  = std::stoi(value); }},
        {"MAXSCF",  [planckCalculator](std::string value){  planckCalculator->max_scf   = std::stoi(value); }},
        {"TOLSCF",  [planckCalculator](std::string value){  planckCalculator->tol_scf   = std::stod(value); }},
        {"TOLERI",  [planckCalculator](std::string value){  planckCalculator->tol_eri   = std::stod(value); }}};

    // now read the input file
    while (getline(*filePointer, inputLine))
    {
        if (inputLine.starts_with("[") && inputLine.ends_with("]"))
        {
            std::string headerLine = inputLine.substr(1, inputLine.length() - 2);
            if (headerLine.compare("SETUP") == 0) // Use proper comparison
            {
                std::string dataLine;
                while (std::getline(*filePointer, dataLine)) // Correct usage of getline
                {
                    if (dataLine.find("END") != std::string::npos)
                    {
                        break;
                    }
                    // process the inputs

                    std::string controlVariable, controlValue;
                    std::stringstream scfBuffer(dataLine);
                    scfBuffer >> controlVariable >> controlValue;

                    // Check if the control variable exists in the map
                    auto it = handlers_setup.find(controlVariable);
                    if (it != handlers_setup.end())
                    {
                        // Call the corresponding handler to set the variable
                        it->second(controlValue);
                    }
                    else
                    {
                        *errorMessage = "Unknown control variable: " + controlVariable;
                        *errorFlag = std::make_error_code(std::errc::invalid_argument);
                        return;
                    }
                }
            }
            if (headerLine.compare("GEOM") == 0)
            {
                std::string dataLine;
                std::getline(*filePointer, dataLine);
                std::stringstream atomBuffer(dataLine);
                atomBuffer >> planckCalculator->total_atoms;

                inputMolecule->atom_masses = (std::double_t *)malloc(sizeof(std::double_t) * planckCalculator->total_atoms);
                inputMolecule->atom_numbers = (std::uint64_t *)malloc(sizeof(std::uint64_t) * planckCalculator->total_atoms);
                inputMolecule->input_coordinates = (std::double_t *)malloc(sizeof(std::double_t) * planckCalculator->total_atoms * 3);

                std::getline(*filePointer, dataLine);
                std::stringstream infoBuffer(dataLine);
                infoBuffer >> planckCalculator->molecule_charge >> planckCalculator->molecule_multiplicity;
                planckCalculator->total_electrons = 0;

                for (std::uint64_t ii = 0; ii < planckCalculator->total_atoms; ii++)
                {
                    std::string atomName;
                    std::double_t xCoord;
                    std::double_t yCoord;
                    std::double_t zCoord;
                    std::string coorLine;

                    std::getline(*filePointer, coorLine);
                    std::stringstream coordBuffer(coorLine);
                    coordBuffer >> atomName >> xCoord >> yCoord >> zCoord;

                    inputMolecule->atom_masses[ii] = atomicMass[atomName];
                    inputMolecule->atom_numbers[ii] = atomicNumber[atomName];
                    inputMolecule->input_coordinates[ii * 3 + 0] = xCoord;
                    inputMolecule->input_coordinates[ii * 3 + 1] = yCoord;
                    inputMolecule->input_coordinates[ii * 3 + 2] = zCoord;
                    planckCalculator->total_electrons += atomicNumber[atomName];
                }
                planckCalculator->total_electrons = planckCalculator->total_electrons + planckCalculator->molecule_charge;
            }
            if (headerLine.compare("CONTROL") == 0) // Use proper comparison
            {
                std::string dataLine;
                while (std::getline(*filePointer, dataLine)) // Correct usage of getline
                {
                    if (dataLine.find("END") != std::string::npos)
                    {
                        break;
                    }
                    // process the inputs

                    std::string controlVariable, controlValue;
                    std::stringstream scfBuffer(dataLine);
                    scfBuffer >> controlVariable >> controlValue;

                    // Check if the control variable exists in the map
                    auto it = handlers_control.find(controlVariable);
                    if (it != handlers_control.end())
                    {
                        // Call the corresponding handler to set the variable
                        it->second(controlValue);
                    }
                    else
                    {
                        *errorMessage = "Unknown control variable: " + controlVariable;
                        *errorFlag = std::make_error_code(std::errc::invalid_argument);
                        return;
                    }
                }
            }
            if (headerLine.find("END") != std::string::npos)
            {
                continue;
            }
        }
    }
    errorFlag = nullptr;
}