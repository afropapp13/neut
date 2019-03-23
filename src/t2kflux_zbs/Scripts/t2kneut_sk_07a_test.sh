#! /bin/csh -f
setenv RFLIST  rflist.$$
echo "Flux dir=" $2
cat <<! >! $RFLIST
10{
{"$2/07a/antinu.sk.1.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
}
20{{"$3",LOCAL,,WRT,,,"recl=5670 status=new"}}
!
./t2kneut_sk $1
