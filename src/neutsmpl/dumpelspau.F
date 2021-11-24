C
      program crssctchk
C
C#include "nesolact.h"

      IMPLICIT NONE

      REAL*4 COST,PHI
      REAL*4    POS(3),DIR(3)

      real*4     ENERGY
      INTEGER*4  IPAR

      INTEGER*4  IPRTBL(4)
      DATA IPRTBL/  14, -14,  12, -12 /

      integer i,j
C      EXTERNAL FNELSPAU,FNELSCRN,FNELSPAN,RSSPICRS,
C     $         RSETCRS,RSSKCRS,FNELSCRS,FNAPICRS,cohcrsne
C      REAL*4   FNELSPAU,FNELSCRN,FNELSPAN,RSSPICRS,
C     $         RSETCRS,RSSKCRS,FNELSCRS,fnapicrs,cohcrsne
C      REAL*4   R1314

C      EXTERNAL FNELSPAU_INTPOL
C      REAL*4   FNELSPAU_INTPOL
      EXTERNAL FNELSPAN
      REAL*4   FNELSPAN

      REAL*4   QECRS(4)
      
      character*80 FILENAME
      data FILENAME/'elspau.tbl'/

      CALL NECARD
      CALL NECARDEV

      IF (IARGC() .GE. 2) THEN
         CALL GETARG(2, FILENAME)
      ENDIF

      open(UNIT=10,FILE=FILENAME,ERR=9999,STATUS='NEW')

      write(10,'(A7,A10,A10,A10,A10)')
     $     'Energy','numu','numubar','nue','nuebar'

      DO 20 J=0,99999
         ENERGY=real(J)*0.01+0.01
         DO 10 I=1,4
            IPAR = IPRTBL(I)
            QECRS(I) = FNELSPAN(ENERGY,IPAR,2212)
C            QECRS(I) = FNELSPAU_INTPOL(ENERGY,IPAR)
 10      CONTINUE
         write(10,'(F7.3,F10.6,F10.6,F10.6,F10.6)')
     $     energy,QEcrs(1),QEcrs(2),QEcrs(3),QEcrs(4)

 20   continue
      goto 2000

 9999 write(*,*) 'Failed to open output file: ', FILENAME

 2000 continue

      CALL NECHADIR(COST,PHI,DIR,POS)
      END
