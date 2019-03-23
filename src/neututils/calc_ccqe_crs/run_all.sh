#!/bin/csh
foreach i ( 1 2 3 )
    @ TYPE = 10 + $i * 2      
   foreach j ( "1.03" "1.11" "1.21" "1.31" )
       @ AXMASS = `echo $j | awk '{print $1 * 1000}' `
       cd ${TYPE}_$j
	   ../run.sh $i $AXMASS >& ${TYPE}_$j &
	   cd ..
	   echo $AXMASS
   end
end   	   

