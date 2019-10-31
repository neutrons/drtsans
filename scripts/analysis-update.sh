#! /bin/bash
set -expu

export PATH=/SNS/software/miniconda2/bin:$PATH
source activate sans
conda install -q -y  -c mantid/label/nightly mantid-framework
pip install /opt/sans-backend