#! /usr/bin/env bash
set -puex

date +%F-%T
hostname -A
hostname -I

# functions
add_instrument_mount () {
    SRC_DIR=$1
    DEST_DIR=$2
    INSTRUMENT=$3

    declare COUNT=0
    if ! grep -qs "${INSTRUMENT}" /proc/mounts; then
      until mount -t nfs -o "tcp,ro,noatime,rsize=32768,wsize=32768" ${SRC_DIR}/${INSTRUMENT} ${DEST_DIR}/${INSTRUMENT}; do
        sleep 10
        COUNT=$((COUNT + 1))
        if [[ ${COUNT} == 5 ]]; then
          sudo reboot
        fi
      done
    fi
    ls ${DEST_DIR}/${INSTRUMENT} | head -1
}

add_sns_mount () {
  INSTRUMENT=$1
  add_instrument_mount "snsdata.ornl.gov:/stornext/snfs1/instruments" "/SNS" ${INSTRUMENT}
}

# add mounts
add_sns_mount 'EQSANS'
