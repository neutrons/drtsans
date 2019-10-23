#! /bin/bash
set -expu

NDAV_USER='jbq'
NDAV_USER_ID='11226'
NDAV_GROUP='sns_ndav_team'
NDAV_GROUP_ID='52321'

groupadd -g "${NDAV_GROUP_ID}" "${NDAV_GROUP}"
useradd -s /bin/bash -u "${NDAV_USER_ID}" -m -g "${NDAV_GROUP}" "${NDAV_USER}"
runuser -l "${NDAV_USER}" -c "bash /opt/sans-backend/scripts/analysis-update.sh"