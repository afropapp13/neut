#!/bin/csh -f

#setenv FC gfortran
#setenv FC g77

unset NEUT_ROOT

#############################################################################
setenv NEUTSMPL ../neutsmpl
source ${NEUTSMPL}/EnvMakeneutsmpl.csh

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
cd ${SOMEWHERE}/t2kflux_zbs

#############################################################################

foreach datfile ($DATFILES)
   ln -s ../crsdat/$datfile
end

#############################################################################
echo "--------- COMPILING ZBSFNS LIBS --------- "
foreach i ( zbsfns )
 
  cd ${SOMEWHERE}/$i
  \rm -rf ${MACHINE} Makefile
  imake_boot 
  make Makefile
  make clean
  make includes
  make all || ( echo "Failed in compiling $i" ; exit 2 )
  make install.include
  make install.lib

end
#############################################################################
cd ${SOMEWHERE}/t2kflux_zbs
imake_boot
make Makefile
make clean
make -f Makefile.t2kneut_sk clean
make -f Makefile.t2kneut_sk necardbmC.h beamntplC.h skheadC.h  || exit
make dump_crs_for_shota_main || exit
make t2kneut neutntpl dumpcrs_main dumptotpau || exit
make split_zbsfile || exit
make t2kneut_ndall || exit
make cmpfiles || exit
#############################################################################

ln -s ../neutcore/nework.h     .
ln -s ../neutcore/mcgenpar.h   .
ln -s ../neutcore/necard.h     .
ln -s ../neutcore/neutparams.h .
ln -s ../neutcore/neutmodel.h  .
ln -s ../neutcore/nefillver.h  .

ln -s ../nuccorspl/nrcard.h    .

ln -s ../skmcsvc/vcwork.h      .
ln -s ../skmcsvc/vcvrtx.h      .

ln -s ../nuceff/fsihist.h      .

make -f Makefile.t2kneut_sk t2kneut_sk skflux_dump || exit
