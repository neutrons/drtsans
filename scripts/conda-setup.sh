#! /usr/bin/env bash

set -x

conda config --add channels conda-forge
conda config --add channels mantid
conda install -q --file requirements.txt
conda install -q ipython pytest
conda install -q -c mantid/label/nightly mantid-framework python=3

