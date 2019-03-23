#! /bin/csh -f
setenv RFLIST  rflist.$$
echo "Output file=" $3
echo "Flux dir=" $2
cat <<! >! $RFLIST
10{
{"$2/nu.nd3_horn250ka.55.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nu.nd3_horn250ka.156.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
}
20{{"$3",LOCAL,,WRT,,,"recl=5670 status=new"}}
!
./Linux_pc/t2kneut_ndall $1

