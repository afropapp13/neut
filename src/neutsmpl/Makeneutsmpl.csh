#!/bin/csh

#############################################################################
source EnvMakeneutsmpl.csh

source Cleanneutsmpl.csh

foreach i ( neutcore nuccorspl nuceff partnuck skmcsvc tauola neutclass radcorr n1p1h ht2p2h)
  cd ${SOMEWHERE}/$i
  make install.include || exit 2
end

foreach i ( neutcore nuccorspl nuceff partnuck skmcsvc tauola neutclass radcorr  n1p1h ht2p2h )
  cd ${SOMEWHERE}/$i
  make all || exit 3
  make install.library || exit 4
end    

echo "--------- COMPILING SPECFUNC --------- "
cd ${SOMEWHERE}/specfunc
make clean
make library || exit 1
make -e install.library || exit 2
make -e exec || exit 3

echo "--------- Make NEUT executables --------- "
cd ${SOMEWHERE}/neutsmpl
make all || exit 4
