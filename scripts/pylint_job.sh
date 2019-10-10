#!/bin/bash
set -x
source activate mantid
pylint --disable=C --disable=no-name-in-module drtsans tests
