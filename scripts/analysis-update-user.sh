#! /bin/bash
set -expu

NDAV_USER='jbq'
NDAV_USER_ID='11226'
NDAV_GROUP='sns_ndav_team'
NDAV_GROUP_ID='52321'
PIP_DIR='/opt/sans-backend'

groupadd -g "${NDAV_GROUP_ID}" "${NDAV_GROUP}"

CURRENT_OWNER="$(ls -ld ${PIP_DIR} | awk '{print $3}')"
if [[ ! "${CURRENT_OWNER}" == "${NDAV_USER}" ]]; then
  CURRENT_OWNER_ID="$(ls -lnd ${PIP_DIR} | awk '{print $3}')"
  useradd -s /bin/bash -u "${CURRENT_OWNER_ID}" -m -g "${NDAV_GROUP}" "${CURRENT_OWNER}"
  runuser -l "${CURRENT_OWNER}" -c "chown -R ${PIP_DIR}"
fi

useradd -s /bin/bash -u "${NDAV_USER_ID}" -m -g "${NDAV_GROUP}" "${NDAV_USER}"
runuser -l "${NDAV_USER}" -c "bash /opt/sans-backend/scripts/analysis-update.sh"