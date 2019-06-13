#! /bin/sh

if [ $(docker 2>/dev/null) ]; then
  if [ $(docker login code.ornl.gov:4567) ]; then
    docker run -it CONTAINER_URL pytest tests/unit/new/ornl/sans/hfir/  
  else
    echo "Login failed. Do you have access to this repository?"
    exit 1
  fi
else
  echo "Docker doesn't seem to be working. Is it installed?"
  exit 1
fi