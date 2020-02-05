#! /bin/bash
set -expu

CONDA_ENV='sans'     # `master` branch only
CONDA_ENV='sans-qa'  # versioned branch to make release candidates from
CONDA_ENV='sans-dev' # `next` branch

export PATH=/SNS/software/miniconda2/bin:$PATH
source activate ${CONDA_ENV}
conda install -q -y  -c mantid/label/nightly mantid-framework
pip install /opt/sans-backend
