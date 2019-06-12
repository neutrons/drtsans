#! /usr/bin/env bash

/opt/anaconda/bin/conda create --name mantid python=3 anaconda
/opt/anaconda/bin/conda activate mantid
/opt/anaconda/bin/conda config --add channels conda-forge
/opt/anaconda/bin/conda config --add channels mantid
/opt/anaconda/bin/conda install --file requirements.txt
/opt/anaconda/bin/conda install -c mantid/label/nightly mantid-framework python=3

