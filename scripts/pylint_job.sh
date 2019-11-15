#!/bin/bash
set -x
source activate mantid
conda install -y pylint
python -m pylint --disable=C --disable=fixme --disable=no-name-in-module drtsans tests --ignore drtsans/_version.py
