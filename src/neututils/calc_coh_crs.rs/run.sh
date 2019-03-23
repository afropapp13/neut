#!/bin/csh -f

./main 12 1
mv cross.dat cross.dat.cc.12
./main 14 1
mv cross.dat cross.dat.cc.14
./main 16 1
mv cross.dat cross.dat.cc.16
./main 12 0
mv cross.dat cross.dat.nc
