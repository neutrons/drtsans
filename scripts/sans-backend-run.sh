#! /bin/sh

set -o errexit -o noclobber -o nounset
VAR_APP_NAME="${0}"
VAR_INTERACT='false'
VAR_UPDATE='false'

func_help_message() {
  printf "  %s is a wrapper for the ORNL Small Angle Neutron Scattering mantid module.

  Usage: %s [-iuh] -- script_1.py script_2.py script_n.py
      -i/--interactive) launches a python terminal rather than running provided scripts.
      -u/--update)      forces an update of the application.
      -h/--help)        prints this message.\n" "${VAR_APP_NAME}" "${VAR_APP_NAME}"
}

func_main() {
  params="$(getopt -o iuh -l iteractive,update,help --name "$0" -- "$@")"
  eval set -- "$params"
  while true; do
    case "${1}" in
        -i|--interactive)
            VAR_INTERACT='true'
            shift
            ;;
        -u|--update)
            VAR_UPDATE='true'
            shift
            ;;
        -h|--help)
            func_help_message
            exit 0
            ;;
        --)
            shift
            break
            ;;
        *)
            func_help_message
            printf "Not implemented: %s" "${1}" >&2
            exit 1
            ;;
    esac
  done
  docker -v 1>/dev/null 2>/dev/null
  if docker -v 1>/dev/null 2>/dev/null; then
    if ${VAR_UPDATE}; then
      docker rmi CONTAINER_URL
    fi
    if docker login code.ornl.gov:4567 2>/dev/null; then
      VAR_TMP_DIR="/tmp/${VAR_APP_NAME}_work_$(date +%s)"
      mkdir -p "${VAR_TMP_DIR}"
      if "${@}"; then
        cp -r $@ "${VAR_TMP_DIR}"/
      fi
      if ${VAR_INTERACT}; then
        docker run -v "${VAR_TMP_DIR}":/tmp/input -it CONTAINER_URL bash
      else
        docker run -v "${VAR_TMP_DIR}":/tmp/input -t CONTAINER_URL bash -c 'find /tmp/input -iname "*.py" -execdir python {} +'
      fi
      if "${VAR_TMP_DIR}"/*; then
        cp -r "${VAR_TMP_DIR}"/* . 2>/dev/null
      fi
    else
      echo "Login failed. Do you have access to this repository?"
      exit 1
    fi
  else
    echo "Docker doesn't seem to be working. Is it installed?"
    exit 1
  fi
}

func_main ${@}