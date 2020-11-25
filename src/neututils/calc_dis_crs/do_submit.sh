#!/bin/csh -f

# General setting
set GROUP = atmpd
set HERE = `pwd`
set list = "mode.list"
set line = 1

if (-e $list) then
    set linenum = `cat $list | wc -l`

    set NUCLPID = `awk '{print$1}' $list`
    set NEUTPID = `awk '{print$2}' $list`

    while($line <= $linenum)
	@ NUCL = $NUCLPID[$line]
	@ NEUT = $NEUTPID[$line]

cat > submit.sh<<EOF
#
# Batch mode using NQS
#

# @\$-q $GROUP
# @\$-o $HERE/log/$NUCL.$NEUT.log
# @\$-e $HERE/log/$NUCL.$NEUT.err
cd $HERE
./Linux_pc/calcapicrs neut.card $NUCL $NEUT $NUCL.$NEUT.dat
EOF

	chmod +x submit.sh
	qsub -q $GROUP submit.sh
#        cat submit.sh
	/bin/rm submit.sh

	@ line ++
    end
else
    echo ++++++++++++++++++++++++++
    echo $list not found
    echo ++++++++++++++++++++++++++
endif
