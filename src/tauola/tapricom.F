      SUBROUTINE TAPRICOM

#include "taucom.h"

      WRITE(6,*) 'DECPAR'
      WRITE(6,*) GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      WRITE(6,*)

      WRITE(6,*) 'PARMAS'
      WRITE(6,*) AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU,
     &           AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1,
     &           AMK,AMKZ,AMKST,GAMKST
      WRITE(6,*)

      WRITE(6,*) 'TAUBRA'
      DO 10 I=1,NCHAN
         WRITE(6,*) I,GAMPRT(I),JLIST(I)
 10   CONTINUE
      WRITE(6,*)

      WRITE(6,*) 'TAUNPI'
      DO 20 I=1,4
        WRITE(6,*)CBRNPI(I),(AMAS(J,I),J=1,6),
     &            (KPI(J,I),J=1,6),MULT(I)
 20   CONTINUE
      WRITE(6,*)

      WRITE(6,*) 'JAKI'
      WRITE(6,*) JAK1,JAK2,JAKP,JAKM,KTOM
      WRITE(6,*)

      WRITE(6,*) 'IDFC'
      WRITE(6,*) IDFF
      WRITE(6,*)

      WRITE(6,*) 'INOUT'
      WRITE(6,*) INUT,IOUT
      WRITE(6,*)

      WRITE(6,*) 'IDPART'
      WRITE(6,*) IA1
      WRITE(6,*)

      WRITE(6,*) 'QEDPRM'
      WRITE(6,*) ALFINV,ALFPI,XK0
      WRITE(6,*)

      WRITE(6,*) 'TAUKLE'
      WRITE(6,*) BRA1,BRK0,BRK0B,BRKS
      WRITE(6,*)

      WRITE(6,*) 'DECOMP'
      DO 30 I=1,NMODE
         WRITE(6,*) I
         WRITE(6,*) (IDFFIN(J,I),J=1,9)
         WRITE(6,'(I6,2X,31A)') MULPIK(I),NAMES(I)
 30   CONTINUE

      RETURN
      END

