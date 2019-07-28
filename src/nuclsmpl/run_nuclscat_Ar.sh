#!/bin/bash
#
# Example script to run piscat program across multiple cores
#


#ln -sf ../crsdat/nucdens_int_coulomb.dat

FILESTART=0
FILEEND=9

PID=2212
if [[ $PID == 2212 ]]; then
    PID_NAME=proton
elif [[ $PID == 2112 ]]; then
    PID_NAME=neutron
else
    exit
fi

PREFIX_NAME=def_ar
#PREFIX_NAME=coul_c


NPROCS=0
for (( i=$FILESTART; i<=$FILEEND; i++ ))
do

    ( 
	./nuclscat neut.Ar.card ${PREFIX_NAME}_${PID_NAME}_${i}.hbk $PID > ${PREFIX_NAME}_${PID_NAME}_${i}.log 2>&1

	if [[ $? == 0 ]]; then

            h2root ${PREFIX_NAME}_${PID_NAME}_${i}.hbk > /dev/null 2>&1
	    
	    rm ${PREFIX_NAME}_${PID_NAME}_${i}.hbk

#            ./run_piscatana -i ${PREFIX_NAME}_${PID_NAME}_${i}.root -o ${PREFIX_NAME}_${PID_NAME}_ana_${i}.root -p 1 -n o > /dev/null 2>&1
	    
	fi
	
    ) &

    FILESTOHADD[$i]=${PREFIX_NAME}_${PID_NAME}_ana_${i}.root

    NPROCS=$(( $NPROCS + 1 ))
    if [[ $NPROCS == 10 ]]; then
        NPROCS=0
        wait
    fi

done

wait

#hadd -f ${PREFIX_NAME}_${PID_NAME}_ana.root ${FILESTOHADD[@]}

#rm ${FILESTOHADD[@]}
