#!/bin/sh

# remove binary files
for file in $(ls bin/|grep -v .py)
do
	rm bin/${file}
done
cd src; make clean


cd ..; rm -r data


# clear example data
# HBB
rm example/HBB_LCR/masker_models.txt; rm -r example/HBB_LCR/output_files

# FADS
rm example/FADS_cluster/masker_models.txt;  rm -r example/FADS_cluster/output_files

# WG
rm -r example/whole_genome/output_files
for file in $(ls example/whole_genome|grep -v .sh)
do 
	rm example/whole_genome/$file
done

git add .; git commit -m "test"; git push gkm main
