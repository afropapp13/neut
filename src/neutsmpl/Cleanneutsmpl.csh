#!/bin/csh

#############################################################################
source EnvMakeneutsmpl.csh


foreach datfile ($DATFILES)
   \rm -f $datfile
end

cd ${SOMEWHERE}

cd ../
rm -rf lib inc include

cd ${SOMEWHERE}/neutcore
\rm -f pdf804
\rm -f vcwork.h  vcvrtx.h  efpion.h   efomega.h  efpiinth.h nrnuclparam.h fsihist.h
\rm -f nemecvcp.F nemeclpv.F fnmec_select.F fnmec_nieves.F meccrs.h nemechad.F neranbll2.F
\rm -f *C.h
cd ..

cd ${SOMEWHERE}/nuccorspl
\rm -f vcwork.h vcvrtx.h posinnuc.h nework.h efpiinth.h fsihist.h
\rm -f *C.h
cd ..

cd ${SOMEWHERE}/nuceff
\rm -f nework.h neutparams.h necard.h vcwork.h vcvrtx.h posinnuc.h
\rm -f *C.h
cd ..

foreach i ( neutcore nuccorspl nuceff partnuck skmcsvc radcorr tauola )
 
  cd ${SOMEWHERE}/$i
  \rm -rf ${MACHINE} Makefile
  imake_boot -DQE111 -DSPI111
  make Makefile
  make clean
  \rm -f *C.h
  \rm -rf $MACHINE
  \rm -f Makefile.bak Makefile

end

cd ${SOMEWHERE}/neutclass
make -f Makefile.old clean

cd ${SOMEWHERE}/specfunc
make clean

cd ${SOMEWHERE}/neutsmpl
echo ---------------------------------------------------------
\rm -f necard.h necardev.h efpion.h vcvrtx.h
\rm -f *C.h

\rm -rf ${MACHINE} Makefile
imake_boot
make clean
\rm -f Makefile.bak Makefile

\rm -f neutmodel.h nework.h vcwork.h neutparams.h posinnuc.h 

make -f GNUmakefile.neutroot clean
\rm -f neutvect.root

\rm -f config_neut.gmk
\rm -f config/site.def

\rm -f ${SOMEWHERE}/../lib/Linux_pc/*.a

\rm -rf ${SOMEWHERE}/../inc
\rm -rf ${SOMEWHERE}/../include
