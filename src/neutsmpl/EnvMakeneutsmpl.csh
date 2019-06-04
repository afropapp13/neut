#!/bin/csh

setenv NEUT_ROOT `pwd`/../..
setenv SOMEWHERE ${NEUT_ROOT}/src
setenv MACHINE `${SOMEWHERE}/neutsmpl/bin/Machine`


#setenv FC g77
setenv FC gfortran
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

${FC} -no-pie -v |& grep "unrecognized" > /dev/null

if ($status == 1) then
   setenv CONFIG_FILE ${SOMEWHERE}/neutsmpl/config_gfortran_nopie.gmk
else
   setenv CONFIG_FILE ${SOMEWHERE}/neutsmpl/config_gfortran.gmk
endif

if ($FC == "g77") then
    ln -s ${SOMEWHERE}/neutsmpl/config/site_g77.def \
			${SOMEWHERE}/neutsmpl/config/site.def
    ln -s ${SOMEWHERE}/neutsmpl/config_g77.gmk \
			${SOMEWHERE}/neutsmpl/config_neut.gmk
else
    ln -s ${SOMEWHERE}/neutsmpl/config/site_gfortran.def \
			${SOMEWHERE}/neutsmpl/config/site.def
    ln -s ${CONFIG_FILE} \
			${SOMEWHERE}/neutsmpl/config_neut.gmk
endif

#setenv CERN /cern
#setenv CERN /usr/local/sklib_gcc4.8.5/cern
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

setenv CERN_ROOT ${CERN}/${CERN_LEVEL}

setenv ROOTSYS /vols/build/t2k/cvw09/root
#setenv ROOTSYS /usr/local/root_v5.34.36
#setenv ROOTSYS /usr/local/sklib_gcc4.8.5/root_v5.28.00h
if (${?ROOTSYS} == 0) then
	echo "set environmental variable ROOTSYS"
	exit 1
endif

setenv PATH ${CERN}/${CERN_LEVEL}/bin:${ROOTSYS}/bin:$PATH
echo $PATH
#setenv LD_LIBRARY_PATH ${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}
#echo $LD_LIBRARY_PATH


setenv EXPERIMENT /
setenv PATH ${SOMEWHERE}/neutsmpl/bin:$PATH
rehash

setenv CVSCOSRC ${SOMEWHERE}/neutsmpl/config
setenv PACKAGE ${SOMEWHERE}/neutsmpl/../
setenv PACKAGE_LEVEL ".."

# Get required data files
# Selecting only data files (.dat) and removing color tags
set DATFILES=(`ls -q ${SOMEWHERE}/crsdat | grep .dat | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g"`)
