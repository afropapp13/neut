#! /bin/csh -f

setenv RANFILE  random.tbl
setenv RFLIST  rflist.$$.`hostname`
setenv NEUTCORE  ../neutcore

/bin/rm -f ccqe_xsec_ma1.1.dat
/bin/rm -f ccqe_xsec_ma1.2.dat
/bin/rm -f 94org_p.dat
/bin/rm -f 94org_pn.dat
/bin/rm -f 94mod_p.dat
/bin/rm -f 94mod_pn.dat
/bin/rm -f 98org_p.dat
/bin/rm -f 98org_pn.dat
/bin/rm -f 98mod_p.dat
/bin/rm -f 98mod_pn.dat

ln -s $NEUTCORE/ccqe_xsec_ma1.1.dat ccqe_xsec_ma1.1.dat
ln -s $NEUTCORE/94org_p.dat  94org_p.dat
ln -s $NEUTCORE/94org_pn.dat 94org_pn.dat
ln -s $NEUTCORE/94mod_p.dat  94mod_p.dat
ln -s $NEUTCORE/94mod_pn.dat 94mod_pn.dat
ln -s $NEUTCORE/98org_p.dat  98org_p.dat
ln -s $NEUTCORE/98org_pn.dat 98org_pn.dat
ln -s $NEUTCORE/98mod_p.dat  98mod_p.dat
ln -s $NEUTCORE/98mod_pn.dat 98mod_pn.dat

cat <<! >! $RFLIST
80{{"ccqe_xsec_ma1.1.dat",LOCAL,,RED,,,"form=formatted"}}
81{{"94org_p.dat",LOCAL,,RED,,,"form=formatted"}}
82{{"94org_pn.dat",LOCAL,,RED,,,"form=formatted"}}
83{{"94mod_p.dat",LOCAL,,RED,,,"form=formatted"}}
84{{"94mod_pn.dat",LOCAL,,RED,,,"form=formatted"}}
85{{"98org_p.dat",LOCAL,,RED,,,"form=formatted"}}
86{{"98org_pn.dat",LOCAL,,RED,,,"form=formatted"}}
87{{"98mod_p.dat",LOCAL,,RED,,,"form=formatted"}}
88{{"98mod_pn.dat",LOCAL,,RED,,,"form=formatted"}}
!

./Linux_pc/neut $1

/bin/rm $RFLIST
/bin/rm ccqe_xsec_ma1.1.dat
/bin/rm 94org_p.dat
/bin/rm 94org_pn.dat
/bin/rm 94mod_p.dat
/bin/rm 94mod_pn.dat
/bin/rm 98org_p.dat
/bin/rm 98org_pn.dat
/bin/rm 98mod_p.dat
/bin/rm 98mod_pn.dat
