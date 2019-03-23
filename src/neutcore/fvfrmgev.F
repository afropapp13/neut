************************************************************************
*     ------------------------
      FUNCTION FVFRMGEV(PFGEV)
*     ------------------------
*
*     (Purpose)
*       probability to nucleon to take fermi momentum pfgev
*
*     (Input)
*       PFGEV  : FERMI MOMENTUM ( GEV )
*     
*     (Output)
*       FVFRMGEV : PROBABILITY TO TAKE PFGEV
*
*     (Creation Date and Author)
*       19??.??.?? ; M.NAKAHATA
*       1995.02.08 ; K. KANEYUKI  FOR S.K.
*                    FRMGEV -> FVFRMGEV
*       2001.04.11 ; call fcpswave, vcppwave
*
*       2003.03.18 ; Change to use flat dist.
************************************************************************
      IMPLICIT NONE

#include "neutparams.h"

      REAL*4 FVFRMGEV
      REAL*4 PFGEV

      IF (PFGEV.LE.PFMAX) THEN
         FVFRMGEV = PFGEV*PFGEV
      ELSE
         FVFRMGEV = 0.
      ENDIF

      RETURN
      END

C*      REAL DISS(10),DISP(12)
C*      DATA DISS/4. , 3.8 , 3.4 , 2.6 , 1.9 , 1.2 , 0.75 , 0.4 , 0.2,
C*     &     0./
C*      DATA DISP/0.8 , 1.6 , 3.4 , 5.9 , 5.5 , 4.9 , 3.8 , 2.5 , 1.5 ,
C*     &     0.8 , 0.4 , 0./
C      DATA SNORM/2.552/,PNORM/1.162E1/
C     
C      FVFRMGEV=0.
C      PF=PFGEV*1.E+3
C*      I=IFIX(PF/25.)+1
C*      FRAC=(PF-25.*FLOAT(I-1))/25.
CC
CC     -- S WAVE
CC      
C      IF (PF.LT.225.) THEN
C*         PSWAVE=DISS(I)+(DISS(I+1)-DISS(I))*FRAC
C*         FVFRMGEV=PSWAVE*PF**2/SNORM*0.25
C         CALL VCPSWAVE(PF,WGHT)
C         FVFRMGEV=WGHT/SNORM*0.25
C      ENDIF
CC     
CC     -- P WAVE
CC      
C      IF (PF.LT.275.) THEN
C*         PPWAVE=DISP(I)+(DISP(I+1)-DISP(I))*FRAC
C*         FVFRMGEV=FVFRMGEV+PPWAVE*PF**2/PNORM*0.75
C         CALL VCPPWAVE(PF,WGHT)
C         FVFRMGEV=FVFRMGEV+WGHT/PNORM*0.75
C      ENDIF
C      
C      RETURN
C      END
C
