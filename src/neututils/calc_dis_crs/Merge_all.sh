#!/bin/bash

HERE=`pwd`
CardFold=${HERE}/CardFiles
TableFold=${HERE}/Output
FinalFold=${HERE}/Merged

mkdir -p ${FinalFold}

for ifile in $CardFold/*.card
do
    ifile=${ifile#"$CardFold/"}
    ifile=${ifile%".card"}
    echo $ifile

    cd ${TableFold}/$ifile
    cp ${HERE}/merge2.sh .
    chmod +x merge2.sh
    ./merge2.sh

    sed -i '1s/^/C\n/' 2212.dat
    sed -i '1s/^/C\n/' 2x12.dat

    cp 2212.dat ${FinalFold}/${ifile}_p.dat
    cp 2x12.dat ${FinalFold}/${ifile}_pn.dat
    
done
