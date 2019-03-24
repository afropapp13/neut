#!/bin/csh

#############################################################################
source EnvMakeneutsmpl.csh

cd ${SOMEWHERE}

rm -rf ${SOMEWHERE}/../lib

rm -rf ${SOMEWHERE}/../inc
rm -rf ${SOMEWHERE}/../include

#############################################################################
echo "--------- COMPILING NEUT LIBS --------- "
foreach i ( neutcore nuccorspl nuceff partnuck specfunc skmcsvc tauola neutclass radcorr n1p1h ht2p2h neutsmpl)
  cd ${SOMEWHERE}/$i
  make clean  || exit 1
end

