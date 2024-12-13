#!/bin/sh


# Download gkm-SVM software package by running the following command:
# ** The lsgkm package is made by Dongwon Lee  **  Forked so that lsgkm update does not affect this pipeline.
git clone https://github.com/oh-jinwoo94/lsgkm.git

# compile 
cd lsgkm/src
make
cd ../..
