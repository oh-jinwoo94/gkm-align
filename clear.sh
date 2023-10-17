#!/bin/sh

# remove binary files
for file in $(ls bin/|grep -v .py)
do
	rm bin/${file}
done
cd src; make clean


cd ..; rm -r data


# clear example data
rm example/HBB_LCR/masker_models.txt; rm -r example/HBB_LCR/output_files
