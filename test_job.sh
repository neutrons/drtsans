#!/bin/bash

v=$1  # python version
module load mantid/ORNL_SANS_py${v}
virtualenv -p $(which python2) testenv${v}
source testenv${v}/bin/activate
pip install -r requirements_test.txt
pytest -v
