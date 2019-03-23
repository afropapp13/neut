#!/bin/sh
#
#  $Id: imake.profile,v 1.1 2007-01-20 07:35:39 skrep Exp $
#
#  Set the "standard" environment variables for Imake
#

# Make sure the experiment is defined.
export EXPERIMENT
if [ "${EXPERIMENT}x" = "x" ]; then
    EXPERIMENT=""
fi

# Set the private root.  This must be overridden to not use.
export PRIVATE_ROOT
if [ "${PRIVATE_ROOT}x" = "x" ]; then
    PRIVATE_ROOT=""
fi

export PACKAGE
if [ "${PACKAGE}x" = "x" ]; then
    PACKAGE=${EXPERIMENT}
fi

# What sort of machine is this.
export MACHINE
if [ "${MACHINE}" = "" ]; then
    if [ -x ${HOME}/bin/Machine ]; then
	IMAKE_PATH=${HOME}/bin
    elif [ -x ${EXPERIMENT}/bin/Machine ]; then
	IMAKE_PATH=${EXPERIMENT}/bin
    elif [ -x ${PACKAGE}/bin/Machine ]; then
	IMAKE_PATH=${PACKAGE}/bin
    else 
	IMAKE_PATH=DeStDiR
    fi
    MACHINE=`${IMAKE_PATH}/Machine`
fi

# Define where to find the Imake exec.
if [ "${IMAKE_PATH}" != "" ]; then
    PATH="${PATH}:${IMAKE_PATH}"
fi

# Choose the machine type.
if [ "${MACHINE}" = solaris_sparc ]; then
    export IMAKECPP
    IMAKECPP=/usr/ccs/lib/cpp
fi

# Provide a default place to look.
CVS_SEDED_NAME=SeDnAmE

# Initialize the CVSCOSRC variable.
export CVSCOSRC
if [ "${CVSCOSRC}" = "" ]; then
    if [ -d ${PRIVATE_ROOT}/lib/config ]; then
	CVSCOSRC=${PRIVATE_ROOT}/lib/config
    elif [ -d ${PRIVATE_ROOT}/config ]; then
	CVSCOSRC=${PRIVATE_ROOT}/config
    elif [ -d ${HOME}/lib/config ]; then
	CVSCOSRC=${HOME}/lib/config
    elif [ -d ${PACKAGE}/lib/config ]; then
	CVSCOSRC=${PACKAGE}/lib/config
    elif [ -d ${PACKAGE}/config ]; then
	CVSCOSRC=${PACKAGE}/config
    elif [ -d ${EXPERIMENT}/lib/config ]; then
	CVSCOSRC=${EXPERIMENT}/lib/config
    elif [ -d ${EXPERIMENT}/config ]; then
	CVSCOSRC=${EXPERIMENT}/config
    elif [ -d ${CVS_SEDED_NAME} ]; then
	CVSCOSRC=${CVS_SEDED_NAME}
    else
	echo You must set CVSCOSRC before you can use imake.
    fi
fi

# If the experiment is defined, add the paths
if [ "${EXPERIMENT}" != "" ]; then 
    export EXPERIMENT_LEVEL
    if [ "${EXPERIMENT_LEVEL}" = "" ]; then
	EXPERIMENT_LEVEL=pro
    fi

    export EXPERIMENT_ROOT
    if [ "${EXPERIMENT_ROOT}" = "" ]; then
	EXPERIMENT_ROOT=${EXPERIMENT}/${EXPERIMENT_LEVEL}
    fi
fi

# If the PACKAGE is defined, add the paths
if [ "${PACKAGE}" != "" ]; then 
    export PACKAGE_LEVEL
    if [ "${PACKAGE_LEVEL}" = "" ]; then
	PACKAGE_LEVEL=${EXPERIMENT_LEVEL}
    fi

    export PACKAGE_ROOT
    if [ "${PACKAGE_ROOT}" = "" ]; then
	PACKAGE_ROOT=${PACKAGE}/${PACKAGE_LEVEL}
    fi

    export PATH
    if [ "${MACHINE}" != "" ]; then
	PATH="${PATH}:${PACKAGE_ROOT}/bin/${MACHINE}"
	PATH="${PATH}:${PACKAGE_ROOT}/bin"
    fi

    export MANPATH
    if [ "${MANPATH}" = "" ]; then
	MANPATH="${PACKAGE_ROOT}/man:"
    else
	MANPATH="${PACKAGE_ROOT}/man:${MANPATH}"
    fi
fi

if [ "${EXPERIMENT}" != "" ]; then 
    export PATH
    if [ "${MACHINE}" != "" ]; then
	PATH="${PATH}:${EXPERIMENT_ROOT}/bin/${MACHINE}"
	PATH="${PATH}:${EXPERIMENT_ROOT}/bin"
    fi

    export MANPATH
    if [ "${MANPATH}" = "" ]; then
	MANPATH="${EXPERIMENT_ROOT}/man:"
    else
	MANPATH="${EXPERIMENT_ROOT}/man:${MANPATH}"
    fi

fi
