#!/bin/bash

set -x

TEST_SCOPE=$1  # 'unit' or 'integration'

source activate mantid
cd /opt/sans-backend
pytest -v /opt/sans-backend/tests/${TEST_SCOPE}