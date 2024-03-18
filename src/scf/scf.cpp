#include "../molecule/molecule.h"
#include "scf.h"

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
        boost::optional<std::uint64_t> atomicNumber_copy = childNode.second.get_optional<std::uint64_t>("<xmlattr>.AtomicNumber");
        boost::optional<std::double_t> atomicMass_copy = childNode.second.get_optional<std::double_t>("<xmlattr>.AtomicMass");

        if (atomName)
        {
            std::uint64_t indexAtom = atomIndex.value_or(0);
            inputMolecule->atomCoordinates(indexAtom, 0) = atomX.value_or(0);
            inputMolecule->atomCoordinates(indexAtom, 1) = atomY.value_or(0);
            inputMolecule->atomCoordinates(indexAtom, 2) = atomZ.value_or(0);
            inputMolecule->atomMasses(indexAtom) = atomicMass_copy.value_or(0);
            inputMolecule->atomNumbers(indexAtom) = atomicNumber_copy.value_or(0);
        }
    }

    // Now start reading the Basis_Information and parse the options
    BOOST_FOREACH (const boost::property_tree::ptree::value_type &basisNode, jobFile.get_child("Data"))
    {
        // If the node is not a Basis set, step over
        if (basisNode.first == "BasisSet")
        {

            boost::optional<std::uint64_t> atomNumber = basisNode.second.get_optional<std::uint64_t>("<xmlattr>.AtomNumber");
            boost::optional<std::uint64_t> atomIndex = basisNode.second.get_optional<std::uint64_t>("<xmlattr>.AtomIndex");

            std::uint64_t indexAtom = atomIndex.value_or(0);

            // Iterate over the CGTO nodes
            BOOST_FOREACH (const boost::property_tree::ptree::value_type &cgtoNode, basisNode.second)
            {
                if (cgtoNode.first == "CGTO")
                {
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

                        // For some reason extra PGTOs were being created
                        // This hack is necessary to remove empty shells
                        if (orbitalCoeff.value_or(0) != 0)
                        {
                            primitiveGTO.orbitalCoeff = orbitalCoeff.value_or(0);
                            primitiveGTO.primitiveExp = primitiveExp.value_or(0);
                            primitiveGTO.orbitalNorms = orbitalNorms.value_or(0);

                            // Append the primitives to form contracted gaussians
                            basisShell.cGTO.push_back(primitiveGTO);
                            scfCalculator->nPrimitives++;
                        }
                    }

                    // Append contracted gaussians to form basis
                    scfCalculator->basisFunctions.push_back(basisShell);
                    scfCalculator->nBasis++;
                }
            }
        }
    }
    scfCalculator->gaussianResults.resize(scfCalculator->nPrimitives, scfCalculator->nPrimitives);
    // scfCalculator->gaussianIntegrals.resize(scfCalculator->nPrimitives * 3, scfCalculator->nPrimitives * 3);
}

// Function to write XML data
void writeXML_GPT(std::fstream *xmlPointer, const cxx_Calculator *scfCalculator, std::error_code *errorFlag, std::string *errorMessage)
{
    // Need to clear the error flag
    errorFlag->clear();

    // Now check if the xml file can be read from
    if (!xmlPointer || !xmlPointer->is_open())
    {
        *errorFlag = std::make_error_code(std::errc::io_error);
        *errorMessage = "Unable to open the JobFile.xml file. Please check if the previous step was successfully completed.";
        return;
    }
    
    // set the prining settings for boost
    boost::property_tree::xml_writer_settings<std::string> settings('\t', 1);
    boost::property_tree::ptree pTree_Job;

    boost::property_tree::read_xml("JobFile.xml", pTree_Job);
    
    for (std::uint64_t ii = 0; ii < scfCalculator->nPrimitives; ++ii)
    {
        boost::property_tree::ptree pTree_Node;
        for (std::uint64_t jj = 0; jj < scfCalculator->nPrimitives; ++jj)
        {
            boost::property_tree::ptree pTree_GaussianProduct;
            pTree_GaussianProduct.put("<xmlattr>.LocationX", scfCalculator->gaussianResults(ii, jj).locationX);
            pTree_GaussianProduct.put("<xmlattr>.LocationY", scfCalculator->gaussianResults(ii, jj).locationY);
            pTree_GaussianProduct.put("<xmlattr>.LocationZ", scfCalculator->gaussianResults(ii, jj).locationZ);
            pTree_GaussianProduct.put("<xmlattr>.IntegralX", scfCalculator->gaussianResults(ii, jj).integralX);
            pTree_GaussianProduct.put("<xmlattr>.IntegralY", scfCalculator->gaussianResults(ii, jj).integralY);
            pTree_GaussianProduct.put("<xmlattr>.IntegralZ", scfCalculator->gaussianResults(ii, jj).integralZ);
            pTree_GaussianProduct.put("<xmlattr>.IndexA", ii);
            pTree_GaussianProduct.put("<xmlattr>.IndexB", jj);
            pTree_Node.add_child("Gaussian", pTree_GaussianProduct);
        }
        pTree_Job.add_child("Data.Output.GPT", pTree_Node);
    }

    // This hack is needed to strip the extra newlines and tabs that appear
    // in the new xml file that gets created automatically.
    // I have absolutely no idea why this is happening, but it seems that
    // the newlines at the end of the original xml file
    std::ostringstream streamBuffer;
    // std::ofstream newFile("JobFile.xml");

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

// Do this once and store the data
void gaussianProducts(cxx_Primitives *primtiveGTO_a, cxx_Primitives *primtiveGTO_b, cxx_gptResults *gptResult)
{
    // Calculate the gaussian center and assign them to gaussian products array
    gptResult->locationX = ((primtiveGTO_a->locationX * primtiveGTO_a->primitiveExp) + (primtiveGTO_b->locationX * primtiveGTO_b->primitiveExp)) / (primtiveGTO_a->primitiveExp + primtiveGTO_b->primitiveExp);
    gptResult->locationY = ((primtiveGTO_a->locationY * primtiveGTO_a->primitiveExp) + (primtiveGTO_b->locationY * primtiveGTO_b->primitiveExp)) / (primtiveGTO_a->primitiveExp + primtiveGTO_b->primitiveExp);
    gptResult->locationZ = ((primtiveGTO_a->locationZ * primtiveGTO_a->primitiveExp) + (primtiveGTO_b->locationZ * primtiveGTO_b->primitiveExp)) / (primtiveGTO_a->primitiveExp + primtiveGTO_b->primitiveExp);

    std::double_t gaussianExponent = (primtiveGTO_a->primitiveExp * primtiveGTO_b->primitiveExp) / (primtiveGTO_a->primitiveExp + primtiveGTO_b->primitiveExp);

    // Calculate the gaussian integrals
    gptResult->integralX = exp(-1 * gaussianExponent * (primtiveGTO_a->locationX - primtiveGTO_b->locationX) * (primtiveGTO_a->locationX - primtiveGTO_b->locationX));
    gptResult->integralY = exp(-1 * gaussianExponent * (primtiveGTO_a->locationY - primtiveGTO_b->locationY) * (primtiveGTO_a->locationY - primtiveGTO_b->locationY));
    gptResult->integralZ = exp(-1 * gaussianExponent * (primtiveGTO_a->locationZ - primtiveGTO_b->locationZ) * (primtiveGTO_a->locationZ - primtiveGTO_b->locationZ));
}

