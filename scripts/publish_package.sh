# Build conda library
set -ex
echo GITHUB REF $CI_COMMIT_REF_SLUG

# activate conda environment
source activate drtsans-dev

# go to the code directory in /tmp
cp -R /opt/sans-backend /tmp/
cd /tmp/sans-backend/conda.recipe

# setup and build the conda package
conda render .
conda mambabuild --output-folder . . -c mantid/label/nightly -c conda-forge -c defaults

# show what tarballs were created
ls */*.tar.bz2

# verify
conda-verify ./noarch/drtsans-*.tar.bz2

# Deploy tags to anaconda.org
if [ -n "${CI_COMMIT_TAG}" ]; then
    CONDA_LABEL="main"
    if [ "${CI_COMMIT_TAG}" = "*rc*" ]; then
        CONDA_LABEL="rc"
    fi
    echo pushing $CI_COMMIT_REF_SLUG with label $CONDA_LABEL
    anaconda upload --label $CONDA_LABEL ./noarch/drtsans-*.tar.bz2
fi
