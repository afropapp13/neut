C
      program crssctchk
C
C#include "nesolact.h"

      IMPLICIT NONE

      DOUBLE PRECISION     ENERGY
      INTEGER*4  IPAR

      INTEGER*4  IPRTBL(4)
      DATA IPRTBL/  14, -14,  12, -12 /

      integer i,j
C      EXTERNAL FNELSPAU,FNELSCRN,FNELSPAN,RSSPICRS,
C     $         RSETCRS,RSSKCRS,FNELSCRS,FNAPICRS,cohcrsne
C      DOUBLE PRECISION   FNELSPAU,FNELSCRN,FNELSPAN,RSSPICRS,
C     $         RSETCRS,RSSKCRS,FNELSCRS,fnapicrs,cohcrsne
C      DOUBLE PRECISION   R1314

      EXTERNAL FNTOTPAU
      DOUBLE PRECISION   FNTOTPAU

      DOUBLE PRECISION   TCRS(4)

      CALL NECARD

      open(UNIT=10,FILE='newtotpau.tbl',ERR=9999,STATUS='NEW')

      DO 20 J=0,599
         ENERGY=real(J)*0.05+0.025
         DO 10 I=1,4
            IPAR = IPRTBL(I)
            TCRS(I) = FNTOTPAU(IPAR,ENERGY)
 10      CONTINUE
         write(10,'(F7.3,F10.6,F10.6,F10.6,F10.6)')
     $     energy,tcrs(1),tcrs(2),tcrs(3),tcrs(4)

 20   continue
      goto 2000

 9999 write(*,*) 'Failed to open output file totpau.tbl'

 2000 continue

      END

      SUBROUTINE RMARIN(int i1,int i2,int i3)

      write(*,*) "RMARIN: Dummy routine"
      return
      end
