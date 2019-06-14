#! /usr/bin/env bash

set -x

conda config --add channels conda-forge
conda config --add channels mantid
conda install -q --file requirements.txt
conda install -q pytest

