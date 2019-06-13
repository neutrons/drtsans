#! /usr/bin/env bash

set -x

conda config --add channels conda-forge
conda config --add channels mantid
conda install --file requirements.txt
conda install ipython pytest
conda install -c mantid/label/nightly mantid-framework python=3

