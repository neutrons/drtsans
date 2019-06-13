#! /bin/bash

set -x

wget -O /tmp/miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x /tmp/miniconda.sh
/tmp/miniconda.sh -b -p /opt/miniconda

