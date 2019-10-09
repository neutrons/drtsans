#!/bin/bash

set -x

TEST_SCOPE=$1  # 'unit' or 'integration'

source activate mantid
cd /opt/sans-backend
echo "Writing tests results to $(pwd)/${TEST_SCOPE}_test_results.xml"
pytest -v /opt/sans-backend/tests/${TEST_SCOPE} -n 7 --junitxml=./${TEST_SCOPE}_test_results.xml
