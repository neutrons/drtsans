# Build conda library
set -ex
echo "GITHUB REF ${CI_COMMIT_REF_SLUG}"

# activate conda environment
source activate drtsans-dev

# go to the code directory in /tmp
cp -R /opt/sans-backend /tmp/
cd /tmp/sans-backend/conda.recipe

# setup and build the conda package
#conda render .
echo "Building conda package"
VERSION=$(python -m versioningit ../) conda mambabuild --output-folder . . -c mantid/label/nightly -c conda-forge || exit 1

# show what tarballs were created
ls */*.tar.bz2

# verify
echo "Verifying conda package"
conda-verify ./noarch/drtsans-*.tar.bz2

# Deploy tags to anaconda.org
if [ -n "${CI_COMMIT_TAG}" ]; then
    if [ -z "${ANACONDA_TOKEN}" ]; then
	echo "ANACONDA_TOKEN is not set"
	exit -1
    fi
    # determine the label for anaconda
    if echo "${CI_COMMIT_TAG}" | grep -q "v.\+rc.\+" ; then
	CONDA_LABEL="rc"
    elif echo "${CI_COMMIT_TAG}" | grep -q "v.\+" ; then
	CONDA_LABEL="main"
    else
	CONDA_LABEL="dev"
    fi
    # push the package
    echo pushing $CI_COMMIT_REF_SLUG with label $CONDA_LABEL
    anaconda --token ${ANACONDA_TOKEN} upload --label $CONDA_LABEL ./noarch/drtsans-*.tar.bz2
else
    echo "Not publishing package to anaconda.org"
fi
