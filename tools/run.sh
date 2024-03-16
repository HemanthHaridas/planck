#!/bin/bash

for file in ./basissets/*/*.basis
do
    python gbs_xml_converter.py $file
done