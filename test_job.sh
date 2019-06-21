#!/bin/bash

set -x

TEST_SCOPE=$1  # 'unit' or 'integration'

N_SUB=8  # number of python subprocesses
#source /etc/profiles.d/modules.sh
#module load mantid/ORNL_SANS_py3.6
source activate mantid
cd /opt/sans-backend
pytest -v /opt/sans-backend/tests/${TEST_SCOPE}