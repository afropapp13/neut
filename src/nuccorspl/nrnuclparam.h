      real*4 apn,ap
      real*4 nrc,nra,nrcnn,nrann,nrdr,nraf,nrcc2,nrdfact
      real*4 nrrmsrad,NRWPARM
      real*4 NRPFSURF
C     add common block used to scale to 16O
      common /nrtaerget/apn,ap,nrc,nra,nrcnn,nrann
     &      ,nrdr,nraf,nrcc2,nrdfact
     &      ,nrrmsrad,NRWPARM,NRPFSURF
C      DATA APN/16./
C      DATA AP/8./
