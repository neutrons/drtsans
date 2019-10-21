#!/bin/bash
set -x
source activate mantid
conda install pylint
python -m pylint --disable=C --disable=fixme --disable=no-name-in-module drtsans tests --ignore drtsans/_version.py
