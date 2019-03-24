#!/bin/csh -f

#setenv FC gfortran
#setenv FC g77

unset NEUT_ROOT

#############################################################################
setenv NEUTSMPL ../neutsmpl
source ${NEUTSMPL}/EnvMakeneutsmpl.csh

cd ${SOMEWHERE}/zbsfns
find . -name "*C.h" -exec \rm \{} \;
find . -name "*.o" -exec \rm \{} \;
find . -name "*.so" -exec \rm \{} \;
cd ${SOMEWHERE}/t2kflux_zbs

#############################################################################
echo "--------- COMPILING ZBSFNS LIBS --------- "
foreach i ( zbsfns )
 
  cd ${SOMEWHERE}/$i
  make clean
  make install.include 
  make library || exit
  make install.library

end
#############################################################################
cd ${SOMEWHERE}/t2kflux_zbs
make clean
make all || exit
#############################################################################
