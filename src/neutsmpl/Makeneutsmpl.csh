#!/bin/csh

#############################################################################
source EnvMakeneutsmpl.csh

foreach datfile ($DATFILES)
   \rm -f $datfile
end
cd ${SOMEWHERE}/neutcore
\rm -f vcwork.h  vcvrtx.h  efpion.h   efomega.h  efpiinth.h nrnuclparam.h fsihist.h
cd ${SOMEWHERE}/nuccorspl
\rm -f vcwork.h vcvrtx.h posinnuc.h nework.h efpiinth.h fsihist.h
cd ${SOMEWHERE}/nuceff
\rm -f nework.h neutparams.h necard.h vcwork.h vcvrtx.h posinnuc.h
cd ${SOMEWHERE}/neutsmpl
\rm -f necard.h necardev.h fsihist.h
\rm -f neutmodel.h nework.h vcwork.h neutparams.h posinnuc.h 
rm -rf ${MACHINE} Makefile

foreach i ( neutcore nuccorspl nuceff partnuck skmcsvc )
  cd ${SOMEWHERE}/$i
  rm -rf ${MACHINE} Makefile
  find . -name "*C.h" -exec \rm \{} \;
end

cd ${SOMEWHERE}/specfunc
make clean

cd ${SOMEWHERE}
find . -name "*.o" -exec \rm \{} \;
find . -name "*.so" -exec \rm \{} \;
rm -rf */${MACHINE} 

#############################################################################

cd ${SOMEWHERE}/neutcore
\rm -f pdf804
\rm -f vcwork.h  vcvrtx.h  efpion.h   efomega.h  efpiinth.h nrnuclparam.h fsihist.h
ln -s ../skmcsvc/vcwork.h  .
ln -s ../skmcsvc/vcvrtx.h  .
ln -s ../nuceff/efpion.h   .
ln -s ../nuceff/efomega.h  .
ln -s ../nuceff/efpiinth.h .
ln -s ../nuceff/fsihist.h .
ln -s ../nuccorspl/nrnuclparam.h .
ln -s ${CERN}/${CERN_LEVEL}/src/mclibs/pdf/pdf pdf804
ln -s ../mec/nemecvcp.F .
ln -s ../mec/nemeclpv.F .
ln -s ../mec/fnmec_select.F .
ln -s ../mec/fnmec_nieves.F .
ln -s ../mec/meccrs.h
ln -s ../mec/nemechad.F .
ln -s ../mec/neranbll2.F .
cd ..

cd ${SOMEWHERE}/nuccorspl
\rm -f vcwork.h vcvrtx.h posinnuc.h nework.h efpiinth.h
ln -s ../skmcsvc/vcwork.h    .
ln -s ../skmcsvc/vcvrtx.h    .
ln -s ../neutcore/posinnuc.h .
ln -s ../neutcore/nework.h   .
ln -s ../nuceff/efpiinth.h   .
ln -s ../nuceff/fsihist.h    .
cd ..

cd ${SOMEWHERE}/nuceff
\rm -f nework.h neutparams.h necard.h vcwork.h vcvrtx.h
ln -s ../neutcore/nework.h     .
ln -s ../neutcore/neutparams.h .
ln -s ../neutcore/posinnuc.h   .
ln -s ../neutcore/necard.h     .
ln -s ../skmcsvc/vcwork.h      .
ln -s ../skmcsvc/vcvrtx.h      .
cd ..

rm -f ${SOMEWHERE}/../lib/Linux_pc/*.a

rm -rf ${SOMEWHERE}/../inc
rm -rf ${SOMEWHERE}/../include

#############################################################################
echo "--------- COMPILING NEUT LIBS --------- "
foreach i ( neutcore nuccorspl nuceff partnuck skmcsvc tauola )
 
  cd ${SOMEWHERE}/$i
  \rm -rf ${MACHINE} Makefile
  imake_boot  || exit 1
  make Makefile  || exit 1
  make clean  || exit 1
  make includes  || exit 1
  make all || exit 2
  make install.include
  make install.lib

end

echo "--------- COMPILING SPECFUNC --------- "
cd ${SOMEWHERE}/specfunc
make clean
make -e install || exit 2
make all || exit 1

# note: neutclass is not compiled like this, it no longer has a Makefile
#echo "--------- COMPILING NEUTCLASS --------- "
#cd ${SOMEWHERE}/neutclass
#make clean
#make all || exit 2

echo "--------- CLEANING NEUTSMPL --------- "
cd ${SOMEWHERE}/neutsmpl
\rm -f necard.h necardev.h
ln -s ../neutcore/necard.h
ln -s ../skmcsvc/necardev.h

echo "--------- bootstrap NEUTSMPL Makefile --------- "
rm -rf ${MACHINE} Makefile
imake_boot || exit 1

echo "--------- Regen. Headers --------- "
\rm -f neutmodel.h nework.h vcwork.h neutparams.h posinnuc.h 
ln -s ../neutcore/neutmodel.h .
ln -s ../neutcore/nework.h .
ln -s ../skmcsvc/vcwork.h .
ln -s ../skmcsvc/vcvrtx.h .
ln -s ../neutcore/neutparams.h .
ln -s ../neutcore/posinnuc.h .
ln -s ../nuceff/efpion.h .

echo "--------- Make NEUT executables --------- "
make neut neut_ntpl dumptotpau dumpelspau dumpcohcrs dumpcrs || exit 2

echo "--------- Make NEUTROOT --------- "
make -f GNUmakefile.neutroot clean  || exit 1
make -f GNUmakefile.neutroot neutroot2  || exit 1

setenv NEUTCORE  ../neutcore

echo "--------- Make LINKS TO XSEC Data files  --------- "
foreach datfile ($DATFILES)
   ln -s ../crsdat/$datfile $datfile
end

echo "--------- Make Links To Spectral Function Data Files --------- "
rm -rf qelSfData
ln -s ${SOMEWHERE}/crsdat/qelSfData qelSfData
