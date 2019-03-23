#!/bin/csh -f

if ( -f cross.dat ) then
    /bin/rm cross.dat
endif

set ITYPE = $1
    @ TYPE = 10 + $ITYPE * 2

    set NUM = 0
    while ($NUM <= 8)
	@ PF = 195 + $NUM * 10
	echo $TYPE $PF

	./main $TYPE  1 $PF
	mv cross.dat tmp1
	./main $TYPE -1 $PF
	mv cross.dat tmp2
	cat tmp1 tmp2 > cross.dat.$TYPE.$PF

	/bin/rm tmp1
	/bin/rm tmp2

	@ NUM ++
    end
    @ ITYPE ++
