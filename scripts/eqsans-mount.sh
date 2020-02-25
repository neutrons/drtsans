#! /usr/bin/env bash
set -puex

date +%F-%T
hostname -A
hostname -I
declare COUNT=0
if ! grep -qs 'SNS' /proc/mounts; then
  until mount -t nfs -o "tcp,ro,noatime,rsize=32768,wsize=32768" snsdata.ornl.gov:/stornext/snfs1/instruments /SNS; do
    sleep 10
    COUNT=$((COUNT + 1))
    if [[ ${COUNT} == 5 ]]; then
      sudo reboot
    fi
  done
fi
ls /SNS

declare COUNT=0
if ! grep -qs 'HFIR' /proc/mounts; then
  until mount -t nfs -o "tcp,ro,noatime,rsize=32768,wsize=32768" hfir-storage.ornl.gov:/volumes/datapool/instruments /SNS; do
    sleep 10
    COUNT=$((COUNT + 1))
    if [[ ${COUNT} == 5 ]]; then
      sudo reboot
    fi
  done
fi
ls /HFIR
