#! /bin/csh -f
setenv RFLIST  rflist.$$
echo "Flux dir=" $2
cat <<! >! $RFLIST
10{
{"$2/nubeam.40gev.sk.1.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.2.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.3.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.4.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.5.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.6.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.7.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.8.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.9.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.10.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.11.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.12.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.13.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.14.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.15.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.16.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.17.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.18.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.19.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.20.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.21.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.22.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.23.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.24.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.25.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.26.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.27.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.28.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.29.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.30.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.31.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.32.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.33.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.34.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.35.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.36.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.37.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.38.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.39.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.40.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.41.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.42.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.43.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.44.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.45.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.46.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.47.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.48.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.49.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
{"$2/nubeam.40gev.sk.50.hbk",LOCAL,,RED,,,"recl=1024 status=old"}
}
!
./skflux_test $11
