#!/bin/bash
set -x
source activate mantid
conda install pylint
python -m pylint --disable=C --disable=no-name-in-module drtsans tests --extension-pkg-whitelist=mantid,numpy --ignore drtsans/_version.py
