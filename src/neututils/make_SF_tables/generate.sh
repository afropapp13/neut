#!/bin/csh
source /usr/local/sklib_g77/cshenv_g77
unsetenv NEUT_ROOT
cd /disk/usr2/hayato/offline/sk4/neut_5.3.1/src/specfunc
date > $4
time ./makeTables.exe $1 $2 $3  >>& $4
date > $4

