*     ------------------------------------------------------
      SUBROUTINE EFPQEABH_PNSEP(PAIMOM,DL,DXXY,NQPAI,NQN,PESC,
     $                          PQE,PINEL,PCOH,PCX)      
*     ------------------------------------------------------
*     
*     ( purpose )
*       GIVE INTERACTION PROBABILITY OF PION AT HIGH ENERGIES
*
*     ( input ) 
C         PAIMOM: PION MOMENTUM (MEV/C)
C         DL    : LENGTH FROM CENTER OF NUCLEI AND PION POSITION(FERMI)
C         DXXY  : UNIT CALCULATE LENGTH (FERMI)
C         NQPAI : CHARGE OF PION
C         NQN   : CHARGE OF NUCLEON
*               : 
*     ( output ) 
C         PESC  : ESCAPE PROBABILITY
C         PEL   : ELASTIC SCATTERING PROBABILITY
C         PINEL : INELASTIC PROBABILITY
C         PCOH  : COHERENT SCATTERING PROBABILITY
*
*     (Creation Date and Author)
*     2007.11.05 ; G.Mitsuka - add conversion with atomic number
*     2010.05    ; P.de Perio - Change to get FACCOH from card (FEFQEHF)
*                             - Changed factor to multiply PELTOT
*                             - Nucleon knockout for nu-interaction and 
*                               not pi-scattering or photoproduction
C     

C     
#include "nework.h"
#include "necard.h"
#include "neutparams.h"
C      COMMON /NVECT/ IMOD,LVECT,IP(20),P(3,20)
#include<efpion.h>

      DATA PRBINT/0.4/
      REAL DENOMINATOR
      REAL*4 PROBSUM
C     
      PESC=1.
      PEL=0.
      PINEL=0.
      PCOH=0.
      PCX=0.
      IF(DL.GT.CC2)RETURN
C     
C     
CC      XGOOUT=FLOAT(LVECT-4-1)*PRBINT
      XGOOUT=FLOAT(NUMNE-4-1)*PRBINT
C     NUMBER OF PROTONS EXPECTED TO GO OUT OF 16O AFTER NUCLEAR
C     INTERACTIONS


C     If rho(r) is taken into account, there is no big changes
C     100.33 assumes that target is 16O
      DENOMINATOR = 100.33*(float(NUMATOM)/16.)
      RHO0=(float(NUMATOM)-1.-(XGOOUT/2.))/DENOMINATOR
      IF (NUSIM.EQ.0) THEN
         RHO0=float(NUMATOM)/DENOMINATOR
      END IF

C     NUCLEAR DENSITY AT THE CENTER OF 16O
C     
      RHO=EFABRHO(DL)
      CALL EFPNCRSH_PNSEP(PAIMOM,NQPAI,NQN,SIGEL,SIGINEL,SIGCX,IERROR)
      IF(IERROR.EQ.1) RETURN
C     
C     CALCULATE PROBABILITY OF EACH PROCESS
      PQETOT  = RHO*DXXY*(SIGEL/10.*RHO0)
      PINEL   = RHO*DXXY*(SIGINEL/10.*RHO0)
      PCX     = RHO*DXXY*(SIGCX/10.*RHO0)

      PQE = PQETOT*FEFQEHF
      PCOH  = PQETOT*(1.-FEFQEHF)

      PINEL = PINEL + PCX*FEFCXHF
      PCX = PCX - PCX*FEFCXHF

      PESC=1.-PQE-PINEL-PCOH-PCX
C     Re-set DXXY using CC2
C      DXXY=CC2/100.
C 
C     Scale parameters
      PQE  =PQE*FEFQEH
      PINEL=PINEL*FEFINEL
      PCOH =PCOH*FEFCOH
      PCX  =PCX*FEFCXH

      PROBSUM=PQE+PINEL+PCOH+PCX
      IF (PROBSUM.gt.1.) then
         PQE  =PQE/PROBSUM
         PINEL=PINEL/PROBSUM
         PCOH =PCOH/PROBSUM
         PCX =PCX/PROBSUM
      ENDIF
      PESC=1-PQE-PINEL-PCOH-PCX

      RETURN
      END
