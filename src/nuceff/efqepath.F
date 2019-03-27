*     -----------------------------------------------------
      FUNCTION EFQEPATH(PPI,RIN)
*     -----------------------------------------------------
*
*     ( purpose )
*       GIVE INELASTIC SCATTER MEAN FREE PATH **(-1)
*
*     ( input )
*       PPI   : PION MOMENTUM(MeV/c)
*       r     : radius
*
*     ( output )
*       QEPATH : REACTION MEAN FREE PATH **(-1)  UNIT (/FERMI)
*
*     ( Creation Date and Author )
*       2007.11.05 ; G.Mitsuka - add support other target than 16O
*                                using efconv2roxy function
C

      IMPLICIT NONE

#include "necard.h"
#include "efpion.h"

      REAL*4 EFQEPATH

      REAL*4 PPI
      REAL*4 DIFFP,DIFFR
      REAL*4 PATH1,PATH2, PATH3,RATIO

      integer*4 i,j

      REAL*4 REFF

      REAL*4    PATH(8,13),PPITBL(13)
      REAL*4    OPATH(18)
      DATA PPITBL/  0.,180.,200.,225.,250.,275.,300.,325.,350.,
     $            375.,400.,425.,500/
      real efconv2roxy,Rorg,Rorgcut,RMSRADOXY,RIN
      external efconv2roxy
      parameter (RMSRADOXY = 2.730)

      DATA PATH/0,0,0,0,0,0,0,0,
     $  .129,.123,.101,5.66E-02,1.68E-02,3.64E-03,5.33E-04,5.30E-05,
     $  .231,.221,.181,.101,3.03E-02,5.78E-03,6.85E-04,5.51E-05,
     $  .394,.384,.325,.184,5.54E-02,9.27E-03,9.05E-04,6.21E-05,
     $  .525,.518,.459,.272,8.20E-02,1.25E-02,1.09E-03,6.97E-05,
     $  .599,.595,.542,.329,9.53E-02,1.31E-02,1.05E-03,6.41E-05,
     $  .632,.629,.577,.343,9.05E-02,1.08E-02,7.75E-04,4.54E-05,
     $  .639,.636,.581,.331,7.67E-02,7.81E-03,5.16E-04,2.92E-05,
     $  .632,.628,.569,.308,6.27E-02,5.59E-03,3.49E-04,1.92E-05,
     $  .619,.614,.552,.284,5.14E-02,4.15E-03,2.49E-04,1.35E-05,
     $  .603,.598,.534,.263,4.29E-02,3.21E-03,1.88E-04,1.00E-05,
     $  .588,.583,.517,.245,3.66E-02,2.58E-03,1.49E-04,7.86E-06,
     $  .588,.583,.517,.245,3.66E-02,2.58E-03,1.49E-04,7.86E-06/

      DATA OPATH/.881,.752,.666,.612,.578,.570,.571,.580,.590,
     $           .600,.600,.600,.610,.610,.619,.638,.667,.696/

      EFQEPATH = 0.

      REFF=RIN
      if(NUMATOM.ne.16) REFF=efconv2roxy(RIN)

      IF (PPI.GT.425) goto 100

      Rorgcut = 7.*(RMSRAD/RMSRADOXY)
      IF (RIN.GT.Rorgcut) RETURN

      DO 10 I=1,13
         IF (PPI.LT.PPITBL(I)) GOTO 15
 10   CONTINUE
 15   CONTINUE

      DO 20 J=1,15
         IF (REFF.LT.(J-1)) GOTO 25
 20   CONTINUE
 25   CONTINUE

      DIFFP = (PPI-PPITBL(I-1))/(PPITBL(I)-PPITBL(I-1))
      DIFFR = (REFF-(J-2))
      PATH1 = PATH(J-1,I-1)+(PATH(J-1,I)-PATH(J-1,I-1))*DIFFP
      PATH2 = PATH(J,I-1)+(PATH(J,I)-PATH(J,I-1))*DIFFP

      EFQEPATH = (PATH1+(PATH2-PATH1)*DIFFR)

      RETURN

 100  IF (PPI.GT.850) then
         write(*,*) 'efqepath: Error :P(PI) > 850.MeV)'
         return
      endif

      I=INT((PPI-425.)/25.)+1

      DO 120 J=1,15
         IF (REFF.LT.(J-1)) GOTO 125
 120  CONTINUE
 125  CONTINUE

      DIFFR = (REFF-(J-2))
      PATH1 = PATH(J-1,13)
      PATH2 = PATH(J,13)
      PATH3 = PATH(1,13)

      RATIO = (PATH1+(PATH2-PATH1)*DIFFR)/PATH3 * 0.66

      DIFFP = (PPI-(I-1)*25.-425.)/25.
      PATH1 = OPATH(I)
      PATH2 = OPATH(I+1)

      EFQEPATH = RATIO*(PATH1 + (PATH2-PATH1)*DIFFP)

      RETURN

      END

