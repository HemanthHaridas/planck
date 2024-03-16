#include <boost/algorithm/string/erase.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <cmath>
#include <cstdint>
#include <eigen3/Eigen/Core>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

#include "tables.h"

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
};

struct cxx_Basis
{
    std::vector<cxx_Primitives> cGTO;
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
    std::uint64_t nBasis;

    std::vector<cxx_Basis> basisFunctions;
};

void readXML(std::fstream *xmlPointer, cxx_Molecule *inputMolecule, cxx_Calculator *scfCalculator, std::error_code *errorFlag, std::string *errorMessage)
{
    // Need to reset the error flag
    errorFlag->clear();

    // Check if a JobFile.xml is present
    if (!xmlPointer || !xmlPointer->is_open())
    {
        *errorFlag = std::make_error_code(std::errc::io_error);
        *errorMessage = "Unable to open the JobFile.xml file.";
        return;
    }

    // Create a boost property tree to parse the JobFile.xml
    boost::property_tree::ptree jobFile;
    boost::property_tree::read_xml(*xmlPointer, jobFile);

    // Set the root of the xml tree based on Input_Options and parse the options
    boost::property_tree::ptree inputRoot = jobFile.get_child("Data.Input_Options");
    boost::optional<std::string> calculationType = inputRoot.get_optional<std::string>("<xmlattr>.CalculationType");
    boost::optional<std::string> basisSet = inputRoot.get_optional<std::string>("<xmlattr>.BasisSet");
    boost::optional<std::string> theoryType = inputRoot.get_optional<std::string>("<xmlattr>.Theory");

    scfCalculator->calculationType = calculationType.get_value_or("");
    scfCalculator->calculationBasis = basisSet.get_value_or("");
    scfCalculator->calculationTheory = theoryType.get_value_or("");

    // Now start reading the Molecule_Information and parse the options
    boost::property_tree::ptree moleculeRoot = jobFile.get_child("Data.Molecule_Information");
    boost::optional<std::uint64_t> nAtoms = moleculeRoot.get_optional<std::uint64_t>("<xmlattr>.NumberAtoms");
    boost::optional<std::int64_t> totalCharge = moleculeRoot.get_optional<std::int64_t>("<xmlattr>.TotalCharge");
    boost::optional<std::uint64_t> molMultiplicity = moleculeRoot.get_optional<std::uint64_t>("<xmlattr>.MolMultiplicity");

    scfCalculator->nAtoms = nAtoms.value_or(0);
    scfCalculator->totalCharge = totalCharge.value_or(0);
    scfCalculator->molMultiplicity = molMultiplicity.value_or(0);

    // Now start reading the coordinate section
    inputMolecule->atomCoordinates.resize(scfCalculator->nAtoms, 3);
    inputMolecule->atomMasses.resize(scfCalculator->nAtoms);
    inputMolecule->atomNumbers.resize(scfCalculator->nAtoms);
    BOOST_FOREACH (const boost::property_tree::ptree::value_type &childNode, moleculeRoot)
    {
        boost::optional<std::string> atomName = childNode.second.get_optional<std::string>("<xmlattr>.AtomName");
        boost::optional<std::uint64_t> atomIndex = childNode.second.get_optional<std::uint64_t>("<xmlattr>.AtomIndex");
        boost::optional<std::double_t> atomX = childNode.second.get_optional<std::double_t>("<xmlattr>.AtomX");
        boost::optional<std::double_t> atomY = childNode.second.get_optional<std::double_t>("<xmlattr>.AtomY");
        boost::optional<std::double_t> atomZ = childNode.second.get_optional<std::double_t>("<xmlattr>.AtomZ");
        boost::optional<std::uint64_t> atomicNumber = childNode.second.get_optional<std::uint64_t>("<xmlattr>.AtomicNumber");
        boost::optional<std::double_t> atomicMass = childNode.second.get_optional<std::double_t>("<xmlattr>.AtomicMass");

        if (atomName)
        {
            std::uint64_t indexAtom = atomIndex.value_or(0);
            inputMolecule->atomCoordinates(indexAtom, 0) = atomX.value_or(0);
            inputMolecule->atomCoordinates(indexAtom, 1) = atomY.value_or(0);
            inputMolecule->atomCoordinates(indexAtom, 2) = atomZ.value_or(0);
            inputMolecule->atomMasses(indexAtom) = atomicMass.value_or(0);
            inputMolecule->atomNumbers(indexAtom) = atomicNumber.value_or(0);
        }
    }

    // Now start reading the Basis_Information and parse the options
    BOOST_FOREACH (const boost::property_tree::ptree::value_type &basisNode, jobFile.get_child("Data"))
    {
        // If the node is not a Basis set, step over
        if (basisNode.first != "BasisSet")
        {
            continue;
        }

        boost::optional<std::uint64_t> atomNumber = basisNode.second.get_optional<std::uint64_t>("<xmlattr>.AtomNumber");
        boost::optional<std::uint64_t> atomIndex = basisNode.second.get_optional<std::uint64_t>("<xmlattr>.AtomIndex");

        std::uint64_t indexAtom = atomIndex.value_or(0);

        // Iterate over the CGTO nodes
        BOOST_FOREACH (const boost::property_tree::ptree::value_type &cgtoNode, basisNode.second)
        {
            if (cgtoNode.first != "CGTO")
            {
                continue;
            }

            // Create a Shell object to hold the information
            cxx_Basis basisShell;

            // Now iterate over PGTO nodes
            BOOST_FOREACH (const boost::property_tree::ptree::value_type &pgtoNode, cgtoNode.second)
            {
                cxx_Primitives primitiveGTO;

                // First set the locations of the primitive gaussians
                primitiveGTO.locationX = inputMolecule->atomCoordinates(indexAtom, 0);
                primitiveGTO.locationY = inputMolecule->atomCoordinates(indexAtom, 1);
                primitiveGTO.locationZ = inputMolecule->atomCoordinates(indexAtom, 2);

                // Now read the information about primitives
                boost::optional<std::double_t> primitiveExp = pgtoNode.second.get_optional<std::double_t>("<xmlattr>.Exponent");
                boost::optional<std::double_t> orbitalCoeff = pgtoNode.second.get_optional<std::double_t>("<xmlattr>.Coefficient");
                boost::optional<std::double_t> orbitalNorms = pgtoNode.second.get_optional<std::double_t>("<xmlattr>.Normalization");

                primitiveGTO.orbitalCoeff = orbitalCoeff.value_or(0);
                primitiveGTO.primitiveExp = primitiveExp.value_or(0);
                primitiveGTO.orbitalNorms = orbitalNorms.value_or(0);

                // Append the primitives to form contracted gaussians
                basisShell.cGTO.push_back(primitiveGTO);
            }

            // Append contracted gaussians to form basis
            scfCalculator->basisFunctions.push_back(basisShell);
        }
    }
}

int main(int argc, const char *argv[])
{
    std::error_code errorFlag;
    std::string errorMessage;
    cxx_Molecule inputMolecule;
    cxx_Calculator scfCalculator;

    std::fstream xmlPointer("JobFile.xml");
    readXML(&xmlPointer, &inputMolecule, &scfCalculator, &errorFlag, &errorMessage);

    return 0;
}