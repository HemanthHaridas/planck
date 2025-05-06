#!/bin/bash

for file in ./basissets/*/*.basis
do
    python generate_basis.py $file
done
