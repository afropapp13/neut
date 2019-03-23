#!/bin/bash
export NEUT_ROOT=/lustre/ific.uv.es/sw/t2k.org/vector2/NEUT/NIWG/neut/branches/neut_5.2.1

if [[ $NEUT_ROOT == "" ]]; then
    echo "Error: You must set the absolute path for NEUT_ROOT in the NEUT_ROOT/src/neutgeom/setup.sh file"
    return
fi
export NEUTGEOM=$NEUT_ROOT/src/neutgeom
export LD_LIBRARY_PATH=$NEUTGEOM:$LD_LIBRARY_PATH
export RANFILE=random.tbl
export CERNLIB=$CERN/$CERN_LEVEL/lib

#link necessary .dat files from neutcore
DATFILES=(`ls ${NEUT_ROOT}/src/crsdat`)

for i in "${DATFILES[@]}"
do
    if [[ $i == *.dat ]]; then 
	    echo "Linking file " $i; 
	    ln -fs $NEUT_ROOT/src/crsdat/$i .
    fi
    if [[ $i == qelSfData ]]; then
	    echo "Linking file " $i; 
	    ln -fs $NEUT_ROOT/src/crsdat/$i .
    fi
done

cp -f $NEUT_ROOT/src/neutsmpl/random.tbl .
cp -f $NEUT_ROOT/src/neutsmpl/neut.card .
