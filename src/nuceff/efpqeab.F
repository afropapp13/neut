*     ------------------------------------------------------
      SUBROUTINE EFPQEAB(PPI,DL,DXXY,PESC,PQE,PABS,PCOH)
*     ------------------------------------------------------
*
*     ( purpose )
*       GIVE INTERACTION PROBABILITY OF PION
*
*     ( input )
C        PPI    : PION MOMENTUM
C        DL     : LENGTH FROM CENTER OF NUCLEI AND PION POSITION(FERMI)
C        DXXY   : UNIT CALCULATE LENGTH (FERMI)
*
*     ( output )
C        PESC  : ESCAPE PROBABILITY
C        PQE   : QUASI-ELASTIC SCATTERING PROBABILITY
C        PABS  : ABSORPTION PROBABILITY
C
*     (Update)
*     2010.06    ; P.de Perio - Change to get FACTOR from card (FEFCOHF)
*                             - Implement pion scattering tuned paths

#include<efpion.h>
#include "neutparams.h"
#include "necard.h"
C
      REAL*4 DLFTBQ(7),DLFTBA(7),RRTBL(7)
      DATA DLFTBQ/1.0,0.999,0.972,1.02,1.47,3.08,3.08/
      DATA DLFTBA/1.0,0.977,0.917,0.826,0.818,1.11,0.758/
      DATA RRTBL/0.0,1.0,2.0,3.0,4.0,5.0,6.0/

      REAL*4 PROBSUM


 20   PESC=1.
      PQE=0.
      PABS=0.
      IF(DL.GT.CC2)RETURN

      IF (NENEFMODL.EQ.1) THEN
         PQE=DXXY*EFQEPATHTUN(PPI,DL)
         PABS=DXXY*EFABPATHTUN(PPI,DL)
      ELSE
         PQE=DXXY*EFQEPATH(PPI,DL)
         PABS=DXXY*EFABPATH(PPI,DL)
      END IF

C     Scale parameters
      PCOH=PQE*FEFCOHF*FEFALL
      PQE =PQE*FEFQE*FEFALL
      PABS=PABS*FEFABS*FEFALL


      PROBSUM=PQE+PABS+PCOH
      IF (PROBSUM.gt.1.) then
         PQE =PQE/PROBSUM
         PABS=PABS/PROBSUM
         PCOH=PCOH/PROBSUM
      ENDIF
      PESC=1-PQE-PABS-PCOH
#ifdef DEBUG
      write(*,*) "PPI=",PPI," : DL=",DL," : PQE=",PQE,
     $           " : PABS=",PABS," : PCOH=",PCOH,
     $           " : PESC=",PESC
#endif
      RETURN
      END
