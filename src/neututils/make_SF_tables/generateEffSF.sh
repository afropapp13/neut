#!/bin/bash

#for targ in He3 He4 C12 O16 Ne20 Al27 Ar40 Fe56 Cu63 Zn64 Pb208
#for targ in Zn64 Pb208
for targ in Fe56 Cu63
do
    
    ## CC
    # ./makeEffSFTables.exe 14  $targ CC &> /dev/null & # makeEffSFTables.14.$targ.CC.log &
    # ./makeEffSFTables.exe -14 $targ CC &> /dev/null & # makeEffSFTables.-14.$targ.CC.log &
    # ./makeEffSFTables.exe 12  $targ CC &> /dev/null & # makeEffSFTables.12.$targ.CC.log &
    # ./makeEffSFTables.exe -12 $targ CC &> /dev/null & # makeEffSFTables.-12.$targ.CC.log &
    # ./makeEffSFTables.exe 16  $targ CC &> /dev/null & # makeEffSFTables.16.$targ.CC.log &
    # ./makeEffSFTables.exe -16 $targ CC &> /dev/null & # makeEffSFTables.-16.$targ.CC.log &
    
    ## NC
    ./makeEffSFTables.exe 14  $targ NC &> /dev/null & # makeEffSFTables.14.$targ.NC.log &
    ./makeEffSFTables.exe -14 $targ NC &> /dev/null & # makeEffSFTables.-14.$targ.NC.log &
    ./makeEffSFTables.exe 12  $targ NC &> /dev/null & # makeEffSFTables.12.$targ.NC.log &
    ./makeEffSFTables.exe -12 $targ NC &> /dev/null & # makeEffSFTables.-12.$targ.NC.log &
    ./makeEffSFTables.exe 16  $targ NC &> /dev/null & # makeEffSFTables.16.$targ.NC.log &
    ./makeEffSFTables.exe -16 $targ NC &> /dev/null & # makeEffSFTables.-16.$targ.NC.log &

    sleep 3600
done