#!/bin/bash

for file in *.basis
do
	python gbs_xml_converter.py $file
done
