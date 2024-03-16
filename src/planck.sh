#!/bin/bash

# First run the sub program to generate the input file for remaining sections 
./molecule.o $1

# If the JobFile.xml was successfully generated, run the scf module
if [ $? -eq 0 ]; then
    ./scf.o
else
    echo "Planck failed with an error $?"
fi