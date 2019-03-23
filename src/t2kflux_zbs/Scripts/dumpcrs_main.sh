#! /bin/csh -f
setenv RFLIST  rflist.$$
echo "Flux dir=" $2
cat <<! >! $RFLIST
./Linux_pc/dumpcrs_main $1
