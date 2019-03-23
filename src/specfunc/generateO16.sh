#!/bin/sh
./makeTables.exe 14 O16 NC  >& makeTables.+14.O16.NC.log &
./makeTables.exe -14 O16 NC >& makeTables.-14.O16.NC.log &
./makeTables.exe 12 O16 NC  >& makeTables.+12.O16.NC.log &
./makeTables.exe -12 O16 NC >& makeTables.-12.O16.NC.log &
./makeTables.exe 16 O16 NC  >& makeTables.+16.O16.NC.log &
./makeTables.exe -16 O16 NC >& makeTables.-16.O16.NC.log &

./makeTables.exe 14 O16 CC  >& makeTables.+14.O16.CC.log &
./makeTables.exe -14 O16 CC >& makeTables.-14.O16.CC.log &
./makeTables.exe 12 O16 CC  >& makeTables.+12.O16.CC.log &
./makeTables.exe -12 O16 CC >& makeTables.-12.O16.CC.log &
./makeTables.exe 16 O16 CC  >& makeTables.+16.O16.CC.log &
./makeTables.exe -16 O16 CC >& makeTables.-16.O16.CC.log &

