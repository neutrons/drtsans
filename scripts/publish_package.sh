# Build conda library
set -e
echo GITHUB REF $CI_COMMIT_REF_SLUG

# cd into the correct directory
THISFILE=$(readlink -f "$0")
DIREC=$(dirname $THISFILE)   # directory of executable
cd "${DIREC}/../conda.recipe"

# setup and build the conda package
conda install -y anaconda-client conda-build conda-verify
conda build --output-folder . . -c neutrons -c mantid/label/nightly

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
