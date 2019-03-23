#!/bin/csh

setenv SOMEWHERE `pwd`/..
setenv MACHINE `${SOMEWHERE}/neutsmpl/bin/Machine`

#setenv FC g77
#setenv FC gfortran
if (${?FC} == 0) then
	echo "EnvMakeneutsmpl: set environmental variable FC "
	exit 1
endif

if (($FC != "g77")&&(($FC != "gfortran"))) then
    echo "EnvMakeneutsmpl: only g77 and gfortran are allowed for the environmental variable FC"
    exit 1
else
    echo "FC=" $FC
endif
rm -f ${SOMEWHERE}/neutsmpl/config/site.def
rm -f ${SOMEWHERE}/neutsmpl/config_neut.gmk
if ($FC == "g77") then
    ln -s ${SOMEWHERE}/neutsmpl/config/site_g77.def \
			${SOMEWHERE}/neutsmpl/config/site.def
    ln -s ${SOMEWHERE}/neutsmpl/config_g77.gmk \
			${SOMEWHERE}/neutsmpl/config_neut.gmk
else
    ln -s ${SOMEWHERE}/neutsmpl/config/site_gfortran.def \
			${SOMEWHERE}/neutsmpl/config/site.def
    ln -s ${SOMEWHERE}/neutsmpl/config_gfortran.gmk \
			${SOMEWHERE}/neutsmpl/config_neut.gmk
endif

#setenv CERN /cern
#setenv CERN /data3/T2K/panos/vector/CERNLIB
if (${?CERN} == 0) then
	echo "set environmental variable CERN"
	exit 1
endif

#setenv CERN_LEVEL pro
#setenv CERN_LEVEL 2005
if (${?CERN_LEVEL} == 0) then
	echo "set environmental variable CERN_LEVEL"
	exit 1
endif

#setenv ROOTSYS /usr/local/sklib_g77/root
#setenv ROOTSYS=/data3/T2K/nd280rep/stable/v10r11p17/ROOT/v5r30p02n01/Linux-x86_64
if (${?ROOTSYS} == 0) then
	echo "set environmental variable ROOTSYS"
	exit 1
endif

setenv PATH ${CERN}/${CERN_LEVEL}/bin:${ROOTSYS}/bin:$PATH
echo $PATH
setenv LD_LIBRARY_PATH ${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}
echo $LD_LIBRARY_PATH


setenv EXPERIMENT /
setenv PATH ${SOMEWHERE}/neutsmpl/bin:$PATH
rehash

setenv CVSCOSRC ${SOMEWHERE}/neutsmpl/config
setenv PACKAGE ${SOMEWHERE}/neutsmpl/../
setenv PACKAGE_LEVEL ".."

# Get required data files
# Selecting only data files (.dat) and removing color tags
set DATFILES=(`ls -q ${SOMEWHERE}/crsdat | grep .dat | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g"`)
