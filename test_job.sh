#!/bin/bash

TEST_SCOPE=$1  # 'unit' or 'integration'

N_SUB=8  # number of python subprocesses

#module load mantid/ORNL_SANS_py${PYTHON_VERSION}

pytest -v -d --tx ${N_SUB} tests/${TEST_SCOPE}
