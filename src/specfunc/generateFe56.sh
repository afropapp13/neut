#!/bin/sh
./makeTables.exe 14 Fe56 NC  >& makeTables.+14.Fe56.NC.log &
./makeTables.exe -14 Fe56 NC >& makeTables.-14.Fe56.NC.log &
./makeTables.exe 12 Fe56 NC  >& makeTables.+12.Fe56.NC.log &
./makeTables.exe -12 Fe56 NC >& makeTables.-12.Fe56.NC.log &
./makeTables.exe 16 Fe56 NC  >& makeTables.+16.Fe56.NC.log &
./makeTables.exe -16 Fe56 NC >& makeTables.-16.Fe56.NC.log &

./makeTables.exe 14 Fe56 CC  >& makeTables.+14.Fe56.CC.log &
./makeTables.exe -14 Fe56 CC >& makeTables.-14.Fe56.CC.log &
./makeTables.exe 12 Fe56 CC  >& makeTables.+12.Fe56.CC.log &
./makeTables.exe -12 Fe56 CC >& makeTables.-12.Fe56.CC.log &
./makeTables.exe 16 Fe56 CC  >& makeTables.+16.Fe56.CC.log &
./makeTables.exe -16 Fe56 CC >& makeTables.-16.Fe56.CC.log &

