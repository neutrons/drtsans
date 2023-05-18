#!/bin/bash

set -x

# Set default values
MARKERS=""
TEST_SCOPE="unit"

# Process command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        -m|--markers)
        MARKERS="-m $2"
        shift # past argument
        shift # past value
        ;;
        -s|--scope)
        TEST_SCOPE="$2"
        shift # past argument
        shift # past value
        ;;
        *)    # unknown option
        echo "Unknown option: $1"
        shift # past argument
        ;;
    esac
done

# the default environment in bash should already be drtsans-dev
source activate drtsans-dev
cd /opt/sans-backend
echo "Writing tests results to $(pwd)/${TEST_SCOPE}_test_results.xml"

pytest --dist loadscope $MARKERS -v /opt/sans-backend/tests/${TEST_SCOPE} -n 4 --junitxml=./${TEST_SCOPE}_test_results.xml
