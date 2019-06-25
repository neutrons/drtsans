#!/bin/bash

PYTHON_VERSION=$1  # python version
TEST_SCOPE=$2  # 'unit' or 'integration'

N_SUB=8  # number of python subprocesses

module load mantid/master_py${PYTHON_VERSION}

virtualenv -p python${PYTHON_VERSION} testenv${PYTHON_VERSION}
source testenv${PYTHON_VERSION}/bin/activate

pip install -r requirements_test.txt
if [ ${PYTHON_VERSION} = "2.7" ]; then
  pip install -r requirements_2.x.txt
fi

pytest -v -d --tx ${N_SUB}*popen//python=python${PYTHON_VERSION} tests/${TEST_SCOPE}
