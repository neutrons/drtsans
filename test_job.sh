#!/bin/bash

set -x

TEST_SCOPE=$1  # 'unit' or 'integration'

source activate mantid
cd /opt/sans-backend
echo "Writing tests results to /tmp/sans-backend/${TEST_SCOPE}_test_results.xml"
pytest -v /opt/sans-backend/tests/${TEST_SCOPE} --junitxml=/tmp/${TEST_SCOPE}_test_results.xml
