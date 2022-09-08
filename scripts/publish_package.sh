# Build conda library
set -e

echo GITHUB REF $CI_COMMIT_REF_SLUG
conda install -y anaconda-client conda-build conda-verify
cd conda.recipe
conda build --output-folder . . -c neutrons -c mantid/label/nightly

# Verify
conda-verify ./linux-64/drtsans-*.tar.bz2

# Deploy tags to anaconda.org
if [ -n "${CI_COMMIT_TAG}" ]; then
    CONDA_LABEL="main"
    if [ "${CI_COMMIT_TAG}" = "*rc*" ]; then
        CONDA_LABEL="rc"
    fi
    echo pushing $CI_COMMIT_REF_SLUG with label $CONDA_LABEL
    anaconda upload --label $CONDA_LABEL ./linux-64/drtsans-*.tar.bz2
fi
