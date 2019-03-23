#! /bin/csh -f
setenv RFLIST  rflist.$$
echo "Output file=" $3
echo "Flux dir=" $2
cat <<! >! $RFLIST
10{
{"$2/nubeam.40gev.nd5.1.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.2.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.3.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.4.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.5.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.6.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.7.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.8.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.9.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.10.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.11.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.12.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.13.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.14.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.15.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.16.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.17.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.18.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.19.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.20.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.21.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.22.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.23.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.24.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.25.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.26.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.27.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.28.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.29.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.30.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.31.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.32.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.33.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.34.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.35.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.36.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.37.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.38.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.39.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.40.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.41.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.42.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.43.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.44.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.45.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.46.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.47.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.48.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.49.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.nd5.50.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
}
20{{"$3",LOCAL,,WRT,,,"recl=5670 status=new"}}
!
./Linux_pc/t2kneut $1
