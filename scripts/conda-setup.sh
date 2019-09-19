#! /usr/bin/env bash
set -x

source activate mantid
conda config --add channels conda-forge
conda config --add channels mantid
conda update -n base -c defaults conda
conda install -q -y --file /opt/sans-backend/requirements_dev.txt
conda install -q -y -c mantid/label/nightly mantid-framework=4
conda clean -afy
python -c "import mantid; print(mantid.__version__)"
