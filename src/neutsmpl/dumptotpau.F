C
      program crssctchk
C
C#include "nesolact.h"

      IMPLICIT NONE

      DOUBLE PRECISION COST,PHI
      DOUBLE PRECISION    POS(3),DIR(3)

      DOUBLE PRECISION     ENERGY
      INTEGER*4  IPAR

      INTEGER*4  IPRTBL(4)
      DATA IPRTBL/  14, -14,  12, -12 /

      integer i,j

      EXTERNAL FNTOTPAU
      DOUBLE PRECISION FNTOTPAU

      DOUBLE PRECISION   TCRS(4)

      character*80 FILENAME
      FILENAME='totpau.tbl'

      CALL NECARD
      CALL NECARDEV

      IF (IARGC() .GE. 2) THEN
         CALL GETARG(2, FILENAME)
      ENDIF

      open(UNIT=10,FILE=FILENAME,ERR=9999,STATUS='REPLACE')

      write(10,'(A7,A10,A10,A10,A10)')
     $     'Energy','numu','numubar','nue','nuebar'

      DO 20 J=0,199
         ENERGY=real(J)*0.05+0.025
         DO 10 I=1,4
            IPAR = IPRTBL(I)
            TCRS(I) = FNTOTPAU(IPAR,ENERGY)
 10      CONTINUE
         write(10,'(F7.3,F10.6,F10.6,F10.6,F10.6)')
     $     energy,tcrs(1),tcrs(2),tcrs(3),tcrs(4)

 20   continue
      goto 2000

 9999 write(*,*) 'Failed to open output file: ', FILENAME

 2000 continue

      CALL NECHADIR(COST,PHI,DIR,POS)
      END

C      SUBROUTINE RMARIN(int i1,int i2,int i3)
C
C      write(*,*) "RMARIN: Dummy routine"
C      return
C      end
