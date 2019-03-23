#! /bin/csh -f
#
# Info: For generating neutrino interaction vectors using pre-filled 
#       histogram file from jnubeam vectors
#
# Usage: ./t2kneut_sk_nu_histo.sh <Card File> <Histogram File> <Output filename>
#
#    where:
#           - <Card File> = A NEUT card file, e.g. in t2kflux_zbs/Cards/
#           - <Histogram File> = File produced by skflux_dump program
#
set DATFILES=(`ls -q ${NEUT_ROOT}/src/crsdat | grep .dat | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g"`)
foreach datfile ($DATFILES)
    ln -sf $NEUT_ROOT/src/crsdat/$datfile .
end

ln -sf $NEUT_ROOT/src/crsdat/qelSfData .

setenv RFLIST  rflist.$$
cat <<! >! $RFLIST
11{
{"$2",LOCAL,,RED,,,"status=old"}
}
20{{"$3",LOCAL,,WRT,,,"recl=5670 status=new"}}
!
${NEUT_ROOT}/src/t2kflux_zbs/./t2kneut_sk $1

