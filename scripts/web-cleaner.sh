#! /bin/bash

declare CI_USER="${1}"
declare CI_PASS="${2}"
declare -a BRANCHES=( $(git ls-remote --heads https://${CI_USER}:${CI_PASS}@code.ornl.gov/sns-hfir-scse/sans/sans-backend.git | cut -f 3 -d '/') )
declare HTTPD_PATH='/var/www/html/sans-backend'
rm -f /tmp/web-cleaner.txt
touch /tmp/web-cleaner.txt

for BRANCH in "${BRANCHES[@]}"; do
  printf "%s/%s\n" "${HTTPD_PATH}" "${BRANCH}" >> /tmp/web-cleaner.txt
done

cat  /tmp/web-cleaner.txt

for DIR in ${HTTPD_PATH}/*; do
  if ! $(grep ""${HTTPD_PATH}/${BRANCH}"" /tmp/web-cleaner.txt); then
    rm -rf "${HTTPD_PATH}/${BRANCH}"
  fi
done

rm -f /tmp/web-cleaner.txt