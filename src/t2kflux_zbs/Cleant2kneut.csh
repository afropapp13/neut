#!/bin/csh

#############################################################################
source ../neutsmpl/EnvMakeneutsmpl.csh

foreach datfile ($DATFILES)
   \rm -f $datfile
end

\rm -f nework.h mcgenpar.h necard.h neutparams.h neutmodel.h nefillver.h
\rm -f nrcard.h
\rm -f vcwork.h vcvrtx.h
\rm -f fsihist.h

find . -name "*C.h" -exec \rm \{} \;
find . -name "*.o" -exec \rm \{} \;
find . -name "*.so" -exec \rm \{} \;

rm -rf ${MACHINE} Makefile Makefile.bak

cd ${SOMEWHERE}/zbsfns
find . -name "*C.h" -exec \rm \{} \;
find . -name "*.o" -exec \rm \{} \;
find . -name "*.so" -exec \rm \{} \;

foreach i ( zbsfns )
 
  cd ${SOMEWHERE}/$i
  \rm -rf ${MACHINE} Makefile
  imake_boot -DQE111 -DSPI111
  make Makefile
  make clean
  \rm -f *C.h
  \rm -rf $MACHINE
  \rm -f Makefile.bak Makefile

end

#############################################################################

setenv NEUTCORE  ../neutcore

#############################################################################
cd ${SOMEWHERE}/t2kflux_zbs
imake_boot
make Makefile
make clean
make -f Makefile.t2kneut_sk clean
\rm -f Makefile Makefile.bak

#############################################################################
