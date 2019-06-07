#!/bin/csh

#############################################################################
source EnvMakeneutsmpl.csh

cd ${SOMEWHERE}

rm -rf ${NEUT_ROOT}/lib

rm -rf ${NEUT_ROOT}/inc
rm -rf ${NEUT_ROOT}/include

#############################################################################
echo "--------- CLEANING NEUT LIBS --------- "
foreach i ( neutcore nuccorspl nuceff partnuck specfunc skmcsvc tauola neutclass radcorr n1p1h ht2p2h neutsmpl)
  cd ${SOMEWHERE}/$i
  make clean  || exit 1
end

