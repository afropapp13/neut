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
      DATA IPRTBL/  14, 12, -14, -12 /

      INTEGER*4  MODE,IMODE
      INTEGER*4  MODNEU(30), MODNEUB(30)

      DATA MODNEU /  1, 11, 12, 13, 21, 31, 32, 33, 34, 41,
     &              51, 51, 52, 16, 36, 22, 42, 43, 23, 44,
     &              45,  0, 26, 46, 17, 38, 39,  2, 15, 35/
      DATA MODNEUB/ -1,-11,-12,-13,-21,-31,-32,-33,-34,-41,
     &              -1,-51,-51,-52,-16,-36,-22,-42,-43,-23,
     &              -44,-45,-26,-46,-17,-38,-39,-2,-15,-35/


      integer i,j
C      EXTERNAL FNELSPAU,FNELSCRN,FNELSPAN,RSSPICRS,
C     $         RSETCRS,RSSKCRS,FNELSCRS,FNAPICRS,cohcrsne
C      REAL*4   FNELSPAU,FNELSCRN,FNELSPAN,RSSPICRS,
C     $         RSETCRS,RSSKCRS,FNELSCRS,fnapicrs,cohcrsne
C      REAL*4   R1314

      EXTERNAL FNrawcrs,fntotpau
      REAL*4   FNrawcrs,fntotpau

      REAL*4   TCRS(31)
          
      character*80 FILENAME
      FILENAME='dumpcrs.tbl'

      CALL NECARD
      CALL NECARDEV

      IF (IARGC() .GE. 2) THEN
         CALL GETARG(2, FILENAME)
      ENDIF

      open(UNIT=10,FILE=FILENAME,ERR=9999,STATUS='REPLACE')


      do 60 I=1,4
         if ((i.eq.1).or.(i.eq.3)) then
            write(10,'(A4,A7,$)') '#IP','Energy'
            if ( i.eq.1 ) then
               do 100 IMODE=1,28
                  write(10,'(I10,$)') MODNEU(IMODE)
 100           continue
            endif
            if ( i.eq.3 ) then
               do 110 IMODE=1,28
                  write(10,'(I10,$)') MODNEUB(IMODE)
 110           continue
            endif
            write(10,'(A10)') 'TOTAL'
         endif

         DO 30 J=0,189
            
            DO 20 IMODE=1,30
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
 20         continue
            
            tcrs(31) = fntotpau(ipar,energy)
            
            write(10,'(I4,F7.3,29F10.6)')
     $           ipar,energy,
     $           tcrs( 1),tcrs( 2),tcrs( 3),tcrs( 4),tcrs( 5),
     $           tcrs( 6),tcrs( 7),tcrs( 8),tcrs( 9),tcrs(10),
     $           tcrs(11),tcrs(12),tcrs(13),tcrs(14),tcrs(15),
     $           tcrs(16),tcrs(17),tcrs(18),tcrs(19),tcrs(20),
     $           tcrs(21),tcrs(22),tcrs(23),tcrs(24),tcrs(25),
     $           tcrs(26),tcrs(27),tcrs(28),tcrs(29),tcrs(30),
     $           tcrs(31)
            
 30      continue
 60   continue
      
      goto 2000
      
 9999 write(*,*) 'Failed to open output file: ', FILENAME

 2000 continue

      CALL NECHADIR(COST,PHI,DIR,POS)
      END
