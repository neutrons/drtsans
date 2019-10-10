#!/bin/bash
set -x
source activate mantid
pylint --disable=C drtsans tests
