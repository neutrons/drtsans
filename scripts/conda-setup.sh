#! /usr/bin/env bash

set -x

conda install -q six==1.11.0 \
                 configparser==3.5.0 \
#                 numpy==1.15.0 \
#                 matplotlib==2.2.2 \
#                 sortedcontainers==1.5.7 \
#                 scikit-image==0.14.0 \
#                 h5py==2.8.0 \
#                 pyyaml==3.13 \
                 pytest

