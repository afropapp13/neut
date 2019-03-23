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

      INTEGER*4  MODE,IMODE
      INTEGER*4  MODNEU(28), MODNEUB(28)

      DATA MODNEU /  1, 11, 12, 13, 21, 31, 32, 33, 34, 41,
     &              51, 51, 52, 16, 36, 22, 42, 43, 23, 44,
     &              45,  0, 26, 46, 17, 38, 39,  2/
      DATA MODNEUB/ -1,-11,-12,-13,-21,-31,-32,-33,-34,-41,
     &              -1,-51,-51,-52,-16,-36,-22,-42,-43,-23,
     &              -44,-45,-26,-46,-17,-38,-39,-2/


      integer i,j
C      EXTERNAL FNELSPAU,FNELSCRN,FNELSPAN,RSSPICRS,
C     $         RSETCRS,RSSKCRS,FNELSCRS,FNAPICRS,cohcrsne
C      REAL*4   FNELSPAU,FNELSCRN,FNELSPAN,RSSPICRS,
C     $         RSETCRS,RSSKCRS,FNELSCRS,fnapicrs,cohcrsne
C      REAL*4   R1314

      EXTERNAL FNrawcrs,fntotpau
      REAL*4   FNrawcrs,fntotpau

      REAL*4   TCRS(29)
          
      character*80 FILENAME
      FILENAME='dumpcrs.tbl'

      CALL NECARD
      CALL NECARDEV

      IF (IARGC() .GE. 2) THEN
         CALL GETARG(2, FILENAME)
      ENDIF

      open(UNIT=10,FILE=FILENAME,ERR=9999,STATUS='REPLACE')

      write(10,'(A7,A7,A10,A10,A10,A10)')
     $     'mode','Energy','numu','numubar','nue','nuebar'

      do 60 I=1,4

      DO 30 J=0,189

         DO 20 IMODE=1,28
            if (j.le.99) then
               ENERGY=real(J)*0.1+0.05
            else
               energy=real(J-100)+10.5
            endif
            IPAR = IPRTBL(I)
            if (IPAR .gt. 0) then
               MODE = MODNEU(IMODE)
            else
               MODE = MODNEUB(IMODE)
            endif
            write(*,*) IMODE,MODE
            TCRS(IMODE) = FNRAWCRS(IPAR,MODE,ENERGY)
 20      continue

         tcrs(29) = fntotpau(ipar,energy)
         
         write(10,'(I4,I7,F7.3,29F10.6)')
     $        ipar,abs(modneu(1)),energy,
     $        tcrs( 1),tcrs( 2),tcrs( 3),tcrs( 4),tcrs( 5),
     $        tcrs( 6),tcrs( 7),tcrs( 8),tcrs( 9),tcrs(10),
     $        tcrs(11),tcrs(12),tcrs(13),tcrs(14),tcrs(15),
     $        tcrs(16),tcrs(17),tcrs(18),tcrs(19),tcrs(20),
     $        tcrs(21),tcrs(22),tcrs(23),tcrs(24),tcrs(25),
     $        tcrs(26),tcrs(27),tcrs(28),tcrs(29)
         
 30   continue
 60   continue
            
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

      
