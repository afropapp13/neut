#!/bin/csh -f

#setenv FC gfortran
#setenv FC g77

#unset NEUT_ROOT

#############################################################################
setenv NEUTSMPL ../neutsmpl
source ${NEUTSMPL}/EnvMakeneutsmpl.csh

\rm -f nework.h mcgenpar.h necard.h neutparams.h neutmodel.h nefillver.h
\rm -f nrcard.h
\rm -f vcwork.h vcvrtx.h
\rm -f fsihist.h

find . -name "*C.h" -exec \rm \{} \;
find . -name "*.o" -exec \rm \{} \;
find . -name "*.so" -exec \rm \{} \;

#############################################################################
echo "--------- COMPILING ZBSFNS LIBS --------- "
foreach i ( zbsfns )
 
  cd ${SOMEWHERE}/$i
  make clean
  make all
  make install.library

end
#############################################################################
cd ${SOMEWHERE}/t2kflux_zbs

make clean

make necardbmC.h beamntplC.h skheadC.h  || exit
make t2kneut neutntpl || exit
make t2kneut_ndall || exit
make t2kneut_sk skflux_dump || exit
