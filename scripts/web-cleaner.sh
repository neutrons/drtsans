#! /bin/bash

declare CI_USER="${1}"
declare CI_PASS="${2}"
declare -a BRANCHES=( $(git ls-remote --heads https://${CI_USER}:${CI_PASS}@code.ornl.gov/sns-hfir-scse/sans/sans-backend.git | cut -f 3 -d '/') )
declare HTTPD_PATH='/var/www/html/sans-backend'

for DIR in ${HTTPD_PATH}/*; do
  if [[ ! "${DIR}" =~ ${HTTPD_PATH}/* ]]; then
    printf "%s is OK!\n" "${DIR}"
  else
    printf "%s is NOT OK!\n" "${DIR}"
  fi
done
  