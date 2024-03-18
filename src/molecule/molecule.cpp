#include "molecule.h"

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

void readInput(std::fstream *filePointer, cxx_Molecule *inputMolecule, cxx_Calculator *scfCalculator, std::error_code *errorFlag, std::string *errorMessage)
{
    // Need to reset the error flag
    errorFlag->clear();

    // Now try reading the input file
    if (!filePointer || !filePointer->is_open())
    {
        *errorFlag = std::make_error_code(std::errc::io_error);
        *errorMessage = "Unable to open the input file. Please check the input file provided";
        return;
    }

    // Read the header section of the input file
    std::string headerLine;
    std::getline(*filePointer, headerLine);
    std::string calculationType;
    std::string basisSet;
    std::string theoryType;
    std::stringstream headerBuffer(headerLine);
    headerBuffer >> calculationType >> basisSet >> theoryType;

    // set the prining settings for boost
    boost::property_tree::xml_writer_settings<std::string> settings('\t', 1);

    // Write to xml file
    boost::property_tree::ptree pTree_Base;
    boost::property_tree::ptree pTree_Input;
    pTree_Input.put("<xmlattr>.CalculationType", calculationType);
    pTree_Input.put("<xmlattr>.BasisSet", basisSet);
    pTree_Input.put("<xmlattr>.Theory", theoryType);
    pTree_Base.add_child("Data.Input_Options", pTree_Input);
    // boost::property_tree::write_xml("JobFile.xml", pTree_Base, std::locale(), settings);

    // Now read number of atoms and charge and multiplicity
    std::getline(*filePointer, headerLine);
    std::uint32_t nAtoms;
    std::stringstream countBuffer(headerLine);
    countBuffer >> nAtoms;

    std::getline(*filePointer, headerLine);
    std::int16_t totalCharge;
    std::uint32_t molMultiplicity;
    std::stringstream infoBuffer(headerLine);
    infoBuffer >> totalCharge >> molMultiplicity;

    // Set relevant information to SCF calculator
    scfCalculator->calculationBasis = basisSet;
    scfCalculator->nAtoms = nAtoms;
    scfCalculator->molMultiplicity = molMultiplicity;
    scfCalculator->totalCharge = totalCharge;

    // Write to xml file
    boost::property_tree::ptree pTree_Molecule;
    pTree_Molecule.put("<xmlattr>.NumberAtoms", nAtoms);
    pTree_Molecule.put("<xmlattr>.TotalCharge", totalCharge);
    pTree_Molecule.put("<xmlattr>.MolMultiplicity", molMultiplicity);
    pTree_Base.add_child("Data.Molecule_Information", pTree_Molecule);
    // boost::property_tree::write_xml("JobFile.xml", pTree_Input, std::locale(), settings);

    // Now initialize the Eigen Matrix to hold the parameters
    inputMolecule->atomCoordinates.resize(nAtoms, 3);
    inputMolecule->atomMasses.resize(nAtoms);
    inputMolecule->atomNumbers.resize(nAtoms);

    // Now read the coordinates and atom names and write the xml file
    std::string atomLine;
    std::string atomName;
    std::double_t atomXCoord;
    std::double_t atomYCoord;
    std::double_t atomZCoord;
    std::uint16_t totalElectrons = 0;
    std::uint16_t atomIndex = 0;

    while (std::getline(*filePointer, atomLine) && atomIndex < nAtoms)
    {
        boost::property_tree::ptree pTree_Coordinates;
        std::stringstream atomBuffer(atomLine);
        atomBuffer >> atomName >> atomXCoord >> atomYCoord >> atomZCoord;

        // Set the coordinates information to the atomCoordinates
        inputMolecule->atomCoordinates(atomIndex, 0) = atomXCoord;
        inputMolecule->atomCoordinates(atomIndex, 1) = atomYCoord;
        inputMolecule->atomCoordinates(atomIndex, 2) = atomZCoord;

        // Set the atom number and mass information to corresponding arrays
        inputMolecule->atomNumbers(atomIndex) = atomicNumber[atomName];
        inputMolecule->atomMasses(atomIndex) = atomicMass[atomName];

        totalElectrons = totalElectrons + atomicNumber[atomName];
        pTree_Coordinates.put("<xmlattr>.AtomName", atomName);
        pTree_Coordinates.put("<xmlattr>.AtomIndex", atomIndex);
        pTree_Coordinates.put("<xmlattr>.AtomX", atomXCoord);
        pTree_Coordinates.put("<xmlattr>.AtomY", atomYCoord);
        pTree_Coordinates.put("<xmlattr>.AtomZ", atomZCoord);
        pTree_Coordinates.put("<xmlattr>.AtomicNumber", atomicNumber[atomName]);
        pTree_Coordinates.put("<xmlattr>.AtomicMass", atomicMass[atomName]);
        pTree_Base.add_child("Data.Molecule_Information.Geometry", pTree_Coordinates);
        atomIndex++;
    }

    // Need to write total number of electrons to the xml file
    totalElectrons = totalElectrons - totalCharge;
    pTree_Base.get_child("Data.Molecule_Information.<xmlattr>").put("TotalElectrons", totalElectrons);
    // boost::property_tree::write_xml("JobFile.xml", pTree_Base, std::locale(), settings);

    // Need to verify if the unpaied electron counts and multiplicity match
    // Update the xml file with the information of alpha and beta electrons
    if (theoryType == "UHF")
    {
        std::uint16_t unpairedElectrons = molMultiplicity - 1;
        bool checkMultiplicity = false;
        for (std::int16_t i = totalElectrons; i >= 0; i -= 2)
        {
            if (unpairedElectrons == i)
            {
                checkMultiplicity = true;
                std::uint16_t alphaElectrons = (totalElectrons / 2) + unpairedElectrons;
                std::uint16_t betaElectrons = totalElectrons - alphaElectrons;

                scfCalculator->alphaElectrons = alphaElectrons;
                scfCalculator->betaElectrons = betaElectrons;

                pTree_Base.get_child("Data.Molecule_Information.<xmlattr>").put("AlphaElectrons", alphaElectrons);
                pTree_Base.get_child("Data.Molecule_Information.<xmlattr>").put("BetaElectrons", betaElectrons);
                // boost::property_tree::write_xml("JobFile.xml", pTree_Base, std::locale(), settings);
                break;
            }
        }
        if (!checkMultiplicity)
        {
            *errorFlag = std::make_error_code(std::errc::invalid_argument);
            *errorMessage = "A combination of " + std::to_string(totalElectrons) + " electrons and multiplicity of " + std::to_string(molMultiplicity) + " is not allowed.";
            return;
        }
    }

    // If no errors were reported, gracefully exit the program
    boost::property_tree::write_xml("JobFile.xml", pTree_Base, std::locale(), settings);
    *errorFlag = std::error_code();
    *errorMessage = "";
}

void readBasis(std::fstream *basisPointer, std::string atomNumber, std::string atomIndex, std::error_code *errorFlag, std::string *errorMessage)
{
    // Need to reset the error flag
    errorFlag->clear();

    // Now try reading the basis file
    if (!basisPointer || !basisPointer->is_open())
    {
        *errorFlag = std::make_error_code(std::errc::io_error);
        *errorMessage = "Unable to open the basis file. The basis file provided does not exist.";
        return;
    }

    // set the prining settings for boost
    boost::property_tree::xml_writer_settings<std::string> settings('\t', 1);
    boost::property_tree::ptree pTree_Basis;
    boost::property_tree::ptree pTree_Job;

    boost::property_tree::read_xml(*basisPointer, pTree_Basis);
    boost::property_tree::read_xml("JobFile.xml", pTree_Job);

    pTree_Basis.get_child("BasisSet").put("<xmlattr>.AtomNumber", atomNumber);
    pTree_Basis.get_child("BasisSet").put("<xmlattr>.AtomIndex", atomIndex);
    pTree_Job.add_child("Data.BasisSet", pTree_Basis.get_child("BasisSet"));

    // This hack is needed to strip the extra newlines and tabs that appear
    // in the new xml file that gets created automatically.
    // I have absolutely no idea why this is happening, but it seems that
    // the newlines at the end of the original xml file
    std::ostringstream streamBuffer;
    std::ofstream newFile("JobFile.xml");

    boost::property_tree::write_xml(streamBuffer, pTree_Job, settings);
    std::string xmlString = streamBuffer.str();
    boost::algorithm::trim(xmlString, std::locale());
    boost::algorithm::erase_all(xmlString, "\n");
    boost::algorithm::erase_all(xmlString, "\t");

    boost::property_tree::ptree outNode;
    std::istringstream outBuffer(xmlString);
    boost::property_tree::read_xml(outBuffer, outNode);
    boost::property_tree::write_xml("JobFile.xml", outNode, std::locale(), settings);
}