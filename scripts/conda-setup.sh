#! /usr/bin/env bash

set -x

#/opt/miniconda/bin/conda init bash
#/opt/miniconda/bin/conda create --name mantid python=3 anaconda
#conda activate mantid
conda config --add channels conda-forge
conda config --add channels mantid
conda install --file requirements.txt
conda install -c mantid/label/nightly mantid-framework python=3

