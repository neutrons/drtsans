#!/bin/bash

# Set default values
MARKERS="not datarepo"
TEST_SCOPE="unit"

# Function to display the help message
display_help() {
    echo "Usage: test_job.sh [OPTIONS]"
    echo "Optional arguments:"
    echo "  -m/--markers VALUE  An valid expression for a pytest marker"
    echo "  -s/--scope   VALUE  One of 'unit' or 'integration'. Default: $TEST_SCOPE"
    echo "  --help              Show this help message"
}

# Process command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        -m|--markers)
        MARKERS="$2"
        shift # past argument
        shift # past value
        ;;
        -s|--scope)
        TEST_SCOPE="$2"
        shift # past argument
        shift # past value
        ;;
        --help)
            display_help
            exit 0
            ;;
        *)    # unknown option
        echo "Unknown option: $1"
        shift # past argument
        ;;
    esac
done

cd /opt/sans-backend
python setup.py develop # have versioningit write drtsans/_version.py
echo "Writing tests results to $(pwd)/${TEST_SCOPE}_test_results.xml"

# Run tests
ARGS_COMMON="-vv --dist loadscope ./tests/${TEST_SCOPE} -n 2 --junitxml=./${TEST_SCOPE}_test_results.xml"
if [ -n "$MARKERS" ]; then
    pytest -m "$MARKERS" ${ARGS_COMMON}
else
    pytest ${ARGS_COMMON}
fi
