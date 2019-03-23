#!/bin/csh
#
# $Id: imake.cshrc,v 1.1.1.1 2016/11/04 20:29:56 mdunkman Exp $
#
#  Set the "standard" environment variables for Imake
#

if ( ${?EXPERIMENT} == 0 ) then
    echo "You should set EXPERIMENT before sourcing this file."
endif

if ( ${?PRIVATE_ROOT} == 0 ) then
    setenv PRIVATE_ROOT ""
endif

if ( ${?PACKAGE} == 0 ) then
    setenv PACKAGE ${EXPERIMENT}
endif

if ( ${?MACHINE} == 0 ) then
    if ( -x ${HOME}/bin/Machine ) then
        set imake_path=${HOME}/bin
    else if ( -x ${EXPERIMENT}/bin/Machine ) then
        set imake_path=${EXPERIMENT}/bin
    else if ( -x ${PACKAGE}/bin/Machine ) then
        set imake_path=${PACKAGE}/bin
    else
	set imake_path=DeStDiR
    endif
    setenv MACHINE `${imake_path}/Machine`
endif

if ( ${?imake_path} != 0 ) then
    set path = ( $path ${imake_path} )
endif

if ( ${?MACHINE} == 1 ) then
    if ( ${MACHINE} == solaris_sparc ) then
	setenv IMAKECPP /usr/ccs/lib/cpp
    endif
endif

if ( ${?CVSCOSRC} == 0 ) then
    set CVS_SEDED_NAME=SeDnAmE
    if ( -d ${PRIVATE_ROOT}/lib/config ) then
	setenv CVSCOSRC ${PRIVATE_ROOT}/lib/config
    else if ( -d ${PRIVATE_ROOT}/config ) then
	setenv CVSCOSRC ${PRIVATE_ROOT}/config
    else if ( -d ${HOME}/lib/config ) then
	setenv CVSCOSRC ${HOME}/lib/config
    else if ( -d ${PACKAGE}/lib/config ) then
	setenv CVSCOSRC ${PACKAGE}/lib/config
    else if ( -d ${PACKAGE}/config ) then
	setenv CVSCOSRC ${PACKAGE}/config
    else if ( -d ${EXPERIMENT}/lib/config ) then
	setenv CVSCOSRC ${EXPERIMENT}/lib/config
    else if ( -d ${EXPERIMENT}/config ) then
	setenv CVSCOSRC ${EXPERIMENT}/config
    else if ( -d ${CVS_SEDED_NAME} ) then
	setenv CVSCOSRC ${CVS_SEDED_NAME}
    else
	echo You must set CVSCOSRC before you can use imake.
    endif
endif

if ( ${?EXPERIMENT} != 0 ) then
    if ( ${?EXPERIMENT_LEVEL} == 0) then
	setenv EXPERIMENT_LEVEL pro
    endif

    if ( ${?EXPERIMENT_ROOT} == 0) then
	setenv EXPERIMENT_ROOT ${EXPERIMENT}/${EXPERIMENT_LEVEL}
    endif
endif

if ( ${?PACKAGE} != 0 ) then
    if ( ${?PACKAGE_LEVEL} == 0) then
	setenv PACKAGE_LEVEL ${EXPERIMENT_LEVEL}
    endif

    if ( ${?PACKAGE_ROOT} == 0) then
	setenv PACKAGE_ROOT ${PACKAGE}/${PACKAGE_LEVEL}
    endif

    if ( ${?MACHINE} != 0 ) then
	set path = ( ${path} ${PACKAGE_ROOT}/bin/${MACHINE} )
	set path = ( ${path} ${PACKAGE_ROOT}/bin )
    endif

    if ( ${?MANPATH} == 0 ) then
	setenv MANPATH "${PACKAGE_ROOT}/man:"
    else
	setenv MANPATH "${PACKAGE_ROOT}/man:$MANPATH"
    endif

endif
    
if ( ${?EXPERIMENT} != 0 ) then
    if ( ${?MACHINE} != 0 ) then
	set path = ( ${path} ${EXPERIMENT_ROOT}/bin/${MACHINE} )
	set path = ( ${path} ${EXPERIMENT_ROOT}/bin )
    endif
    
    if ( ${?MANPATH} == 0 ) then
	setenv MANPATH "${EXPERIMENT_ROOT}/man:"
    else
	setenv MANPATH "${EXPERIMENT_ROOT}/man:$MANPATH"
    endif

endif
