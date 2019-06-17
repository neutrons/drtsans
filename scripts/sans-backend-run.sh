#! /bin/bash

set -o errexit -o noclobber -o nounset
VAR_APP_NAME="${0}"
VAR_INTERACT='false'
VAR_UPDATE='false'

func_help_message() {
  printf "  %s is a wrapper for the ORNL Small Angle Neutron Scattering mantid module.

  Usage: %s [-iuh] script_1.py script_2.py script_n.py
      -i) launches a python terminal rather than running provided scripts.
      -u) forces an update of the application.
      -h) prints this message.\n" "${VAR_APP_NAME}" "${VAR_APP_NAME}"
}

func_chk_perms() {
  if [ ! "$(id -nG | grep -o docker)" ]; then
    if [ "$(id -u)" -ne 0 ]; then
      printf "Error: You must either be in the docker group or run as root (sudo)."
      exit 1
    fi
  fi
}

func_main() {
  while getopts "i:u:-h" OPT; do
    case "${OPT}" in
        i)
            VAR_INTERACT='true';;
        u)
            VAR_UPDATE='true';;
        h)
            func_help_message
            exit 0;;
        -)
            shift
            break;;
        *)
            func_help_message
            printf "Error: Not implemented: %s" "${1}" >&2
            exit 1;;
    esac
    shift
  done
  docker -v 1>/dev/null 2>/dev/null
  if docker -v 1>/dev/null 2>/dev/null; then
    if [ "${VAR_UPDATE}" = 'true' ]; then
      docker pull CONTAINER_URL
    fi
    if docker login code.ornl.gov:4567 2>/dev/null; then
      if ${VAR_INTERACT}; then
        docker run -v "$PWD":/tmp/input -it CONTAINER_URL bash
      else
        docker run -v "$PWD":/tmp/input -t CONTAINER_URL bash -c "bash scripts/run-inputs.sh ${*}"
      fi
    else
      echo "Error: Login failed. Do you have access to this repository?"
      exit 1
    fi
  else
    echo "Error: Docker doesn't seem to be working. Is it installed?"
    exit 1
  fi
}

func_chk_perms
func_main ${@}