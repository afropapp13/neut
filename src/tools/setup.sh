#!/bin/bash

### Adapted from https://unix.stackexchange.com/questions/4965/keep-duplicates-out-of-path-on-source
function add_to_PATH () {
  for d; do

    d=$(cd -- "$d" && { pwd -P || pwd; }) 2>/dev/null  # canonicalize symbolic links
    if [ -z "$d" ]; then continue; fi  # skip nonexistent directory

    if [ "$d" == "/usr/bin" ] || [ "$d" == "/usr/bin64" ] || [ "$d" == "/usr/local/bin" ] || [ "$d" == "/usr/local/bin64" ]; then
      case ":$PATH:" in
        *":$d:"*) :;;
        *) export PATH=$PATH:$d;;
      esac
    else
      case ":$PATH:" in
        *":$d:"*) :;;
        *) export PATH=$d:$PATH;;
      esac
    fi
  done
}

function add_to_LD_LIBRARY_PATH () {
  for d; do

    d=$(cd -- "$d" && { pwd -P || pwd; }) 2>/dev/null  # canonicalize symbolic links
    if [ -z "$d" ]; then continue; fi  # skip nonexistent directory

    if [ "$d" == "/usr/lib" ] || [ "$d" == "/usr/lib64" ] || [ "$d" == "/usr/local/lib" ] || [ "$d" == "/usr/local/lib64" ]; then
      case ":$LD_LIBRARY_PATH:" in
        *":$d:"*) :;;
        *) export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$d;;
      esac
    else
      case ":$LD_LIBRARY_PATH:" in
        *":$d:"*) :;;
        *) export LD_LIBRARY_PATH=$d:$LD_LIBRARY_PATH;;
      esac
    fi
  done
}

#if it was sourced as . setup.sh then you can't scrub off the end... assume that
#we are in the correct directory.
if ! echo "${BASH_SOURCE}" | grep "/" --silent; then
  SETUPDIR=$(readlink -f $PWD)
else
  SETUPDIR=$(readlink -f ${BASH_SOURCE%/*})
fi

export NEUT=${SETUPDIR}
export NEUT_CARDS=${SETUPDIR}/share/neut/Cards
export NEUT_CRSDAT=${SETUPDIR}/share/neut/crsdat

add_to_PATH ${NEUT}/bin
add_to_LD_LIBRARY_PATH ${NEUT}/lib

unset SETUPDIR