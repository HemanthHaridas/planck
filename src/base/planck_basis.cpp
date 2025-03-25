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

#include "planck_basis.h"

void readBasis(cxx_Molecule *inputMolecule, cxx_Calculator *planckCalculator, std::error_code *errorFlag, std::string *errorMessage)
{
    // first reset the error flag and error message
    errorFlag->clear();
    errorMessage->erase();

    std::uint64_t cgtoIndex = 0;
    std::uint64_t pgtoIndex = 0;

    for (std::uint64_t atomIndex = 0; atomIndex < planckCalculator->total_atoms; atomIndex++)
    {
        std::string basisFile = planckCalculator->basis_path + "/" + planckCalculator->calculation_basis + "-" + std::to_string(inputMolecule->atom_numbers[atomIndex]) + ".xml";
        std::fstream basisPointer(basisFile);
        // std::cout << basisFile << "\n";

        // now check if the basis file exists
        if (!basisPointer.is_open() || !basisPointer)
        {
            *errorFlag = std::make_error_code(std::errc::io_error);
            *errorMessage = "Unable To open basis file";
            return;
        }

        // now read the basis file
        boost::property_tree::ptree basisObject;
        boost::property_tree::read_xml(basisPointer, basisObject);
        boost::property_tree::ptree basisRoot = basisObject.get_child("BasisSet");

        BOOST_FOREACH (const boost::property_tree::ptree::value_type &basisNode, basisRoot)
        {
            // now read each primitive gaussians
            cxx_Contracted basisShell;
            basisShell.shell_x = basisNode.second.get_optional<std::int64_t>("<xmlattr>.AngularMomentumX").get_value_or(0);
            basisShell.shell_y = basisNode.second.get_optional<std::int64_t>("<xmlattr>.AngularMomentumY").get_value_or(0);
            basisShell.shell_z = basisNode.second.get_optional<std::int64_t>("<xmlattr>.AngularMomentumZ").get_value_or(0);

            // first set the location of the contracted gaussian
            basisShell.location_x = inputMolecule->standard_coordinates[atomIndex * 3 + 0];
            basisShell.location_y = inputMolecule->standard_coordinates[atomIndex * 3 + 1];
            basisShell.location_z = inputMolecule->standard_coordinates[atomIndex * 3 + 2];

            BOOST_FOREACH (const boost::property_tree::ptree::value_type &pgtoNode, basisNode.second)
            {
                cxx_Primitive primitiveGTO;

                // if the node is not PGTO, get the angular momentum information
                if (pgtoNode.first != "PGTO")
                {
                    continue;
                }

                // now read the remaining information
                boost::optional<std::double_t> primitiveExp = pgtoNode.second.get_optional<std::double_t>("<xmlattr>.Exponent");
                boost::optional<std::double_t> orbitalCoeff = pgtoNode.second.get_optional<std::double_t>("<xmlattr>.Coefficient");
                boost::optional<std::double_t> orbitalNorms = pgtoNode.second.get_optional<std::double_t>("<xmlattr>.Normalization");

                primitiveGTO.primitive_exp = primitiveExp.get_value_or(0);
                primitiveGTO.orbital_coeff = orbitalCoeff.get_value_or(0);
                primitiveGTO.orbital_norm = orbitalNorms.get_value_or(0);
                primitiveGTO.pgto_index = pgtoIndex;
                basisShell.contracted_GTO.push_back(primitiveGTO);
                pgtoIndex++;
            }
            basisShell.cgto_index = cgtoIndex;
            planckCalculator->calculation_set.push_back(basisShell);
            cgtoIndex++;
        }
    }
    planckCalculator->total_basis = cgtoIndex;
    planckCalculator->total_primitives = pgtoIndex;
}