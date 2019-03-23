#!/bin/csh -f
setenv RFLIST  rflist.$$
rm -if $RFLIST
touch $RFLIST
echo '10{' > $RFLIST
while ($2 != '')
cat <<! >> $RFLIST
{"$1",LOCAL,,RED,,,"recl=5670 status=old"}
!
shift
end
echo '}' >> $RFLIST
cat $RFLIST
echo $1
#./${MACHINE}/neutntpl_tq
./neutntpl_tq
mv genvec.nt $1
\/bin/rm $RFLIST
