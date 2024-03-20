#include <omp.h>

#include "molecule/molecule.h"
#include "scf/scf.h"
#include "scf/overlap.h"
#include "auxiliary/license.h"

int main(int argc, char const *argv[])
{
    std::error_code errorFlag;
    std::string errorMessage;
    cxx_Molecule inputMolecule;
    cxx_Calculator scfCalculator;

    // Print the GPL3 license agreement

    printLicenseAgreement(licenseAgreement);

    // Check if the input file has been provided
    if (argc < 2)
    {
        errorFlag = std::make_error_code(std::errc::no_such_file_or_directory);
        errorMessage = "No input file has been provided. Run planck as planck inputfile";
        std::cerr << std::setw(10) << std::left << "[Error] :" << std::setw(50) << std::left << errorMessage << "\n";
        exit(errorFlag.value());
    }

    std::string inputFile = argv[1];
    std::fstream filePointer(inputFile);
    readInput(&filePointer, &inputMolecule, &scfCalculator, &errorFlag, &errorMessage);

    // Check if the input file has been read successfully
    if (errorFlag.value())
    {
        std::cerr << "\n";
        std::cerr << std::setw(10) << std::left << "[Error] :" << std::setw(50) << std::left << errorMessage << "\n";
        exit(errorFlag.value());
    }

    // Now read the basis sets
    for (std::uint16_t atomIndex = 0; atomIndex < scfCalculator.nAtoms; ++atomIndex)
    {
        // Updated the path to indicate that the basis folder will be packaged along with the executables
        std::string basisFile = "./basis/" + scfCalculator.calculationBasis + "-" + std::to_string(inputMolecule.atomNumbers(atomIndex)) + ".xml";
        std::fstream basisPointer(basisFile);
        readBasis(&basisPointer, std::to_string(inputMolecule.atomNumbers(atomIndex)), std::to_string(atomIndex), &errorFlag, &errorMessage);

        // Check if the input file has been read successfully
        if (errorFlag.value())
        {
            std::cerr << std::setw(10) << std::left << "[Error] :" << std::setw(50) << std::left << errorMessage << "\n";
            exit(errorFlag.value());
        }
    }

    // Once the basis sets are written to JobFile.xml, read the basis set to start the caluclations
    std::fstream xmlPointer("JobFile.xml");
    readXML(&xmlPointer, &inputMolecule, &scfCalculator, &errorFlag, &errorMessage);

    // First we need to calculate gaussian products and store them
    std::cout << "\n";
    std::cout << std::setw(10) << std::left << "[Planck] There are " << scfCalculator.nPrimitives << " Primtive Gaussians from " << scfCalculator.nBasis << " Contracted Gaussians"
              << "\n";
    std::cout << std::setw(10) << std::left << "[Planck] Allocated " << sizeof(cxx_gptResults) * scfCalculator.resultSCF.gaussianResults.cols() * scfCalculator.resultSCF.gaussianResults.rows() << " B for Gaussian Products."
              << "\n";
    std::cout << "\n";
    
    // Original version

    // #pragma omp parallel for
    //     for (std::uint64_t ii = 0; ii < scfCalculator.nBasis; ++ii)
    //     {
    //         for (std::uint64_t ij = 0; ij < scfCalculator.basisFunctions[ii].cGTO.size(); ++ij)
    //         {
    //             rowIndex = (ii std::cout << std::left <<  scfCalculator.basisFunctions[ii].cGTO.size()) + ij;
    //             cxx_Primitives primtiveGTO_a = scfCalculator.basisFunctions[ii].cGTO[ij];

    //             for (std::uint64_t jj = 0; jj < scfCalculator.nBasis; ++jj)
    //             {
    //                 for (std::uint64_t ji = 0; ji < scfCalculator.basisFunctions[jj].cGTO.size(); ++ji)
    //                 {
    //                     cxx_Primitives primtiveGTO_b = scfCalculator.basisFunctions[jj].cGTO[ji];
    //                     colIndex = (jj std::cout << std::left <<  scfCalculator.basisFunctions[ii].cGTO.size()) + ji;

    //                     cxx_gptResults gptResult;
    //                     gaussianProducts(&primtiveGTO_a, &primtiveGTO_b, &gptResult);
    //                     scfCalculator.resultSCF.gaussianResults(rowIndex, colIndex) = gptResult;
    //                 }
    //             }
    //         }
    //     }

#pragma omp parallel for collapse(2)
    for (std::uint64_t ii = 0; ii < scfCalculator.nBasis; ++ii)
    {
        for (std::uint64_t jj = 0; jj < scfCalculator.nBasis; ++jj)
        {
            for (std::uint64_t ij = 0; ij < scfCalculator.basisFunctions[ii].cGTO.size(); ++ij)
            {
                cxx_Primitives primtiveGTO_a = scfCalculator.basisFunctions[ii].cGTO[ij];

                for (std::uint64_t ji = 0; ji < scfCalculator.basisFunctions[jj].cGTO.size(); ++ji)
                {
                    cxx_Primitives primtiveGTO_b = scfCalculator.basisFunctions[jj].cGTO[ji];
                    cxx_gptResults gptResult;
                    gaussianProducts(&primtiveGTO_a, &primtiveGTO_b, &gptResult);
                    scfCalculator.resultSCF.gaussianResults(gptResult.indexA, gptResult.indexB) = gptResult;
                }
            }
        }
    }
    writeXML_GPT(&xmlPointer, &scfCalculator, &errorFlag, &errorMessage);
    
    return 0;
}