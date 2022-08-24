#!/bin/bash

source activate drtsans-dev
cd /opt/sans-backend
python -m build --wheel --no-isolation
check-wheel-contents dist/drtsans-*.whl
