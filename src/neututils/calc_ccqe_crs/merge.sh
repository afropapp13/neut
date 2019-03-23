#!/bin/csh -f

\rm -f tmp.*

# NUE
set NUM = 0
echo "MAVAL RMAVAL" > tmp.format.12
awk '{printf"%5d %9.3f\n",$1,$2/1000.}' cross.dat.12.205 >> tmp.format.12
foreach file(cross.dat.12.*)
    @ NUM ++
	set PF=`echo $file | sed "s/cross.dat.12.//"`
	echo 0.$PF > tmp.$NUM.12
    awk '{print$3}' $file >> tmp.$NUM.12
end
paste tmp.format.12 tmp.?.12 > tmp.crs.12

# NUMU
set NUM = 0
echo "MAVAL RMAVAL 0.205" > tmp.format.14
awk '{printf"%5d %9.3f\n",$1,$2/1000.}' cross.dat.14.205 > tmp.format.14
foreach file(cross.dat.14.*)
    @ NUM ++
    awk '{print$3}' $file >> tmp.$NUM.14
end
paste tmp.format.14 tmp.?.14 > tmp.crs.14

# NUTAU
set NUM = 0
echo "MAVAL RMAVAL 0.205" > tmp.format.16
awk '{printf"%5d %9.3f\n",$1,$2/1000.}' cross.dat.16.205 > tmp.format.16
foreach file(cross.dat.16.*)
    @ NUM ++
    awk '{print$3}' $file >> tmp.$NUM.16
end
paste tmp.format.16 tmp.?.16 > tmp.crs.16

cat tmp.crs.12 tmp.crs.14 tmp.crs.16 > crosssection.txt
\rm tmp.*
