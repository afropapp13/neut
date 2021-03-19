#!/bin/bash

GROUP=all

HERE=`pwd`
CardFold=${HERE}/CardFiles
bincode=${HERE}/calcapicrs
oFold=${HERE}/Output

targ=(2112 2212)
nu=(12 14 16 -12 -14 -16)


for ifile in $CardFold/*.card
do
    ifile=${ifile#"$CardFold/"}
    ifile=${ifile%".card"}
    echo $ifile

    mkdir -p ${oFold}/$ifile
    cd ${oFold}/$ifile
    
    for itarg in $(seq 0 1)
    do
#	echo ${targ[${itarg}]}
	for inu in $(seq 0 5)
	do
	    #	    echo ${nu[${inu}]}
	    cat > ${ifile}.${targ[${itarg}]}.${nu[${inu}]}.sh <<EOF
#!/bin/sh
time $bincode ${CardFold}/${ifile}.card ${targ[${itarg}]} ${nu[${inu}]} ${targ[${itarg}]}.${nu[${inu}]}.dat   

EOF
	    chmod +x ${ifile}.${targ[${itarg}]}.${nu[${inu}]}.sh
	    qsub -q $GROUP -eo -o ${ifile}.${targ[${itarg}]}.${nu[${inu}]}.log ${ifile}.${targ[${itarg}]}.${nu[${inu}]}.sh
	done	
    done
    cd ${HERE}
done

