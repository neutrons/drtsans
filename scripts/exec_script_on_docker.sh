#! /usr/bin/env bash
set -e

SCRIPT=$1

sudo rm -rf /tmp/sans-backend || true
sudo mkdir -p /tmp/sans-backend
sudo cp -r . /tmp/sans-backend
sudo chmod 777 /tmp/sans-backend
pushd /tmp/sans-backend
docker login --username=$CI_REGISTRY_USER --password=$CI_REGISTRY_PASSWORD $CI_REGISTRY
time docker pull $CONTAINER_URL
time docker run -v /SNS:/SNS -v /HFIR:/HFIR -v $PWD:/opt/sans-backend -t $CONTAINER_URL bash -c "bash /opt/sans-backend/${SCRIPT}"
