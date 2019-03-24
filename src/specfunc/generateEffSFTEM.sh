#!/bin/bash

for targ in O16 Ne20 Al27 Ar40 Fe56 Cu63 Zn64 Pb208
do
    
    ## CC
    # ./makeEffSFTables.exe 14  $targ CC TEM &> /dev/null & # makeEffSFTables.14.$targ.CC.log &
    # ./makeEffSFTables.exe -14 $targ CC TEM &> /dev/null & # makeEffSFTables.-14.$targ.CC.log &
    # ./makeEffSFTables.exe 12  $targ CC TEM &> /dev/null & # makeEffSFTables.12.$targ.CC.log &
    # ./makeEffSFTables.exe -12 $targ CC TEM &> /dev/null & # makeEffSFTables.-12.$targ.CC.log &
    # ./makeEffSFTables.exe 16  $targ CC TEM &> /dev/null & # makeEffSFTables.16.$targ.CC.log &
    # ./makeEffSFTables.exe -16 $targ CC TEM &> /dev/null & # makeEffSFTables.-16.$targ.CC.log &
    
    ## NC
    ./makeEffSFTables.exe 14  $targ NC TEM &> /dev/null & # makeEffSFTables.14.$targ.NC.log &
    ./makeEffSFTables.exe -14 $targ NC TEM &> /dev/null & # makeEffSFTables.-14.$targ.NC.log &
    ./makeEffSFTables.exe 12  $targ NC TEM &> /dev/null & # makeEffSFTables.12.$targ.NC.log &
    ./makeEffSFTables.exe -12 $targ NC TEM &> /dev/null & # makeEffSFTables.-12.$targ.NC.log &
    ./makeEffSFTables.exe 16  $targ NC TEM &> /dev/null & # makeEffSFTables.16.$targ.NC.log &
    ./makeEffSFTables.exe -16 $targ NC TEM &> /dev/null & # makeEffSFTables.-16.$targ.NC.log &   

    sleep 3600
done