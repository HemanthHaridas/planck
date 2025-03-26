#include <algorithm>
#include <iostream>
#include <iomanip>
#include <set>
#include <string.h>

#include "planck_symmetry.h"

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

void detectSymmetry(cxx_Molecule *inputMolecule, std::uint64_t nAtoms, std::error_code *errorFlag, std::string *errorMessage)
{
    // first clear the error flag and error message buffers
    errorFlag->clear();
    errorMessage->clear();

    msym_thresholds_t sloppy_thresholds = {
        0.08,   // zero
        0.1,    // geometry
        0.1,    // angle
        0.06,   // equivalence
        1.0e-1, // permutation
        1.0e-3, // eigfact
        0.1     // orthogonalization
    };

    std::cout << std::setw(21) << std::left << "[Planck]   => " << std::setw(80) << std::left << "We use libmsym library to detect point groups" << "\n";
    std::cout << std::setw(20) << std::left << "[Planck] " << "\n";
    std::cout << std::setw(20) << std::left << "[Planck] " << "\n";

    // Need to create a msym_element_arry for using libmsym
    // this section of the code is unfortunately in C
    msym_error_t ret = MSYM_SUCCESS;
    msym_element_t *elements = NULL;
    char point_group[6];
    unsigned int length = nAtoms;

    msym_element_t *geometry = NULL;
    msym_element_t *atom;
    atom = (msym_element_t *)malloc(nAtoms * sizeof(msym_element_t));
    memset(atom, 0, nAtoms * sizeof(msym_element_t));

    for (std::uint64_t ii = 0; ii < nAtoms; ii++)
    {
        atom[ii].m = inputMolecule->atom_masses[ii];  // mass
        atom[ii].n = inputMolecule->atom_numbers[ii]; // nuclear charge

        // now set positions
        atom[ii].v[0] = inputMolecule->input_coordinates[ii * 3 + 0];
        atom[ii].v[1] = inputMolecule->input_coordinates[ii * 3 + 1];
        atom[ii].v[2] = inputMolecule->input_coordinates[ii * 3 + 2];
    }
    geometry = atom;
    msym_context symm_context = msymCreateContext();
    msymSetThresholds(symm_context, &sloppy_thresholds);

    // need to initialize the eigen matrix to hold standardized coordinates
    inputMolecule->standard_coordinates = (std::double_t *)malloc(sizeof(std::double_t) * nAtoms * 3);

    // now check the point group of the molecule
    if (MSYM_SUCCESS != (ret = msymSetElements(symm_context, length, geometry)))
    {
        free(geometry);
        *errorMessage = "Unable to set elements.";
        *errorFlag = std::make_error_code(std::errc::invalid_argument);
        return;
    }

    if (MSYM_SUCCESS != (ret = msymFindSymmetry(symm_context)))
    {
        free(geometry);
        *errorMessage = "Unable to find molecule point group. Continuing with C1";
        inputMolecule->point_group = "C1";
        inputMolecule->standard_coordinates = inputMolecule->standard_coordinates;
        inputMolecule->is_reoriented = false;
        return;
    }

    if (MSYM_SUCCESS != (ret = msymGetPointGroupName(symm_context, sizeof(char[6]), point_group)))
    {
        free(geometry);
        return;
    }

    inputMolecule->point_group = point_group;
    free(geometry);

    // if it has a Cinf axis, replace 0 with inf
    if (point_group[1] == '0')
    {
        inputMolecule->point_group = inputMolecule->point_group.replace(1, 1, "inf");
    }

    // symmetrize the molecule before that
    msym_element_t *new_geometry = NULL;
    double symm_error = 0.0;
    int new_n_atoms = 0;

    // first check if the symmetrization is success
    if (MSYM_SUCCESS != (ret = msymSymmetrizeElements(symm_context, &symm_error)))
    {
        *errorMessage = "Unable to symmetrize the molecule to a higher point group.";
        *errorFlag = std::make_error_code(std::errc::invalid_argument);
        return;
    }

    if (MSYM_SUCCESS != (ret = msymGetElements(symm_context, &new_n_atoms, &new_geometry)))
    {
        *errorMessage = "Unable to get the symmetry elements for " + inputMolecule->point_group + "after symmetrization.";
        *errorFlag = std::make_error_code(std::errc::invalid_argument);
        return;
    }

    // align the axes correctly
    if (MSYM_SUCCESS != (ret = msymAlignAxes(symm_context)))
    {
        *errorMessage = "Unable to align the symmetry axes.";
        *errorFlag = std::make_error_code(std::errc::invalid_argument);
        return;
    }
    // now update the new geometry
    for (std::uint64_t ii = 0; ii < nAtoms; ii++)
    {
        inputMolecule->standard_coordinates[ii * 3 + 0] = (new_geometry[ii].v[0]);
        inputMolecule->standard_coordinates[ii * 3 + 1] = (new_geometry[ii].v[1]);
        inputMolecule->standard_coordinates[ii * 3 + 2] = (new_geometry[ii].v[2]);
    }

    inputMolecule->is_reoriented = true;
}