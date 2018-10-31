#!/bin/bash

v=$1  # python version
test_scope=$2  # 'unit' or 'integration'
module load mantid/ORNL_SANS_py${v}
virtualenv -p python${v} testenv${v}
source testenv${v}/bin/activate
pip install -r requirements_test.txt
pytest -v tests/${test_scope}
