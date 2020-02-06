#! /bin/bash
set -expu

case "$1" in

    master) CONDA_ENV='sans'     # `master` branch only
            ;;
    qa)     CONDA_ENV='sans-qa'  # versioned branch to make release candidates from
            ;;
    next)   CONDA_ENV='sans'     # `next` branch FIXME this should be changed back to sans-dev before the release !!!
            ;;
esac

export PATH=/SNS/software/miniconda2/bin:$PATH
source activate ${CONDA_ENV}
conda install -q -y  -c mantid/label/nightly mantid-framework
pip install /opt/sans-backend
