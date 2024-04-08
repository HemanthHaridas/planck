#include "inputoutput.h"

void printCoordinates(cxx_Molecule *moleculeObject)
{
    for (std::uint64_t ii = 0; ii < moleculeObject->nAtoms; ii++)
    {
        std::cout << std::setw(5) << std::left << moleculeObject->atomNumbers(ii); 
        std::cout << std::setw(20) << std::right << std::fixed << moleculeObject->atomCoordinates(ii, 0);
        std::cout << std::setw(20) << std::right << std::fixed << moleculeObject->atomCoordinates(ii, 1);
        std::cout << std::setw(20) << std::right << std::fixed << moleculeObject->atomCoordinates(ii, 2);
        std::cout << "\n";
    }
}