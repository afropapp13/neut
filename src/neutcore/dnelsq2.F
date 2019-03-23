************************************************************************
*     ----------------------------------
      real*8 FUNCTION DNELSQ2(E,IPAR,Q2)
*     ----------------------------------
*
*     (Purpose)
*       elastic without pauli principle 
*
*     (Input)
*       E      : NEUTRINO ENERGY ( GEV )      
*       IPAR   : NEUTRINO TYPE
*                  12 : NEU-E
*                 -12 : NEU-E-BAR            
*                  14 : NEU-MU
*                 -14 : NEU-MU-BAR            
*                  16 : NEU-TAU
*                 -16 : NEU-TAU-BAR
*       Q2     : Q**2
*     
*     (Output)
*       DNELSQ2 : Q**2 DISTRIBUTION OF ELASTIC SCATTERING
*
*     (Creation Date and Author)
*       1983.11.20 ; M.Nakahata
*       1985.11.08 ; M.Nakahata  Bug fix
*       1987.05.21 ; N.Sato      Tau
*       1988.06.28 ; M.Takita    unified all species      
*       1995.02.02 ; K. KANEYUKI  FOR S.K.
*                    ELSQ2 -> DNELSQ2
*       2000.08.10 ; K. KANEYUKI  CONST 0.166 -> 0.176 ( update G_F )
*       2002.09.03 ; Y.Hayato     MA=1.01 -> 1.1
*       2002.09.03 ; Y.Hayato     MA=1.01 -> 1.1
*       2009.01.11 ; Y.Hayato     G_A value : 1.23 to 1.267
*       2009.02.25 ; Y.Nakajima   Add F_P
*
************************************************************************
C
C    CONST ; G**2*COS(THETAC)**2*M**2/8/PI
C    PM    ; PROTON MASS
C    EM    ; ELECTRON MASS
C    UM    ; MUON MASS
C    TM    ; TAUON MASS
C    AVMASS; COEEFICIENT OF AXIAL VECTOR FORM FACTOR
C                CHANGED FROM 0.95 TO 1.01 BY T.KAJITA 25-JUN-88
C                1.01+-0.03 IS THE AVERAGED BEST VALUE. SEE REF.
C                   (REF.K.L.MILLER ET.AL, PHYS.REV.D.26(1982)537.)
C            Now input from CARD (2010.10 P. de Perio)

      IMPLICIT REAL*8 (A-H,O-Z)

#include "neutmodel.h"
#include "rscons.h"

      REAL FNQ2MAX,FNQ2MIN,EE
      EXTERNAL FNQ2MAX, FNQ2MIN
C      DATA CONST/1.66E-39/
C      DATA CONST/1.66D-1/     MODIFICATION K.Kaneyuki 10-AUG-2000 
      DATA CONST/1.76D-1/      
C      DATA PM/0.938/
C      DATA EM/0.511E-3/
C      DATA UM/0.1057/
C      DATA TM/1.7771/
C      PARAMETER (PM = XMN)
      PARAMETER (EM = XME)
      PARAMETER (UM = XMMU)
      PARAMETER (TM = XMTAU)

C--   taken from A. Bodek et al, Eur. Phys. J. C 71 (2011) 1726 
      REAL TRCORA, TRCORB
      PARAMETER (TRCORA = 6.0)
      PARAMETER (TRCORB = 0.34)

CCCC  DATA TM/1.7842/     MODIFICATION  K.Kaneyuki 2-FEB-95
C      DATA TM/1.7842/
      DATA PIM/0.13957/

      INTEGER QEMODE

      if (IPAR.GT.0) THEN
         PM = XMNE
      else
         PM = XMP
      endif
      AVMASS = XMAQE
      VVMASS = XMVQE
    
      IF(IABS(IPAR).EQ.12)AM=EM
      IF(IABS(IPAR).EQ.14)AM=UM
      IF(IABS(IPAR).EQ.16)AM=TM
C-------
      EE=E
      Q2T1=FNQ2MIN(EE,IPAR)
      IF(Q2.LT.Q2T1) THEN
         DNELSQ2=0.
         RETURN
      END IF
      Q2T2=FNQ2MAX(EE,IPAR)
      IF(Q2.GT.Q2T2) THEN
         DNELSQ2=0.
         RETURN
      END IF
C-------
      SQ2=-Q2
     
      QEMODE = mod(mod(MDLQE,100),10)

      if (QEMODE.eq.3) THEN
         CALL BBBA07(Q2,GPE,GPM,GNE,GNM)
         GVE=GPE-GNE
         GVM=GPM-GNM
      else if (QEMODE.eq.2) THEN
         CALL BBBA05(Q2,GPE,GPM,GNE,GNM)
         GVE=GPE-GNE
         GVM=GPM-GNM
      else if (QEMODE.eq.1) THEN 
         GVE=1./(1.-SQ2/VVMASS**2)**2
         GVM=(1.+3.71)/(1.-SQ2/VVMASS**2)**2
      else
         write(*,*) 'Invalid CCQE mode specified by MDLQE',MDLQE
         stop
      endif

      if (MDLQE/100.eq.1) THEN
         GVM = GVM * sqrt(1+TRCORA * Q2 * exp( SQ2 / TRCORB ))
      endif

      F1V=(GVE-SQ2*GVM/4./PM**2)/(1.-SQ2/4./PM**2)
      F2VG=(GVM-GVE)/(1.-SQ2/4./PM**2)
C      FA=-1.23/(1.-SQ2/AVMASS**2)**2

      if (MDLQEAF.eq.1) then
         FA=-1.267/(1.-SQ2/AVMASS**2)**2
      else if (MDLQEAF.eq.2) then
         FA=faBBBA07(q2)
      else
         write(*,*) 'Invalid CCQE mode specified by MDLQEAF',MDLQEAF
         stop
      endif
      FP=2.*(PM**2)*FA/(-SQ2+PIM**2)
     
C       AAA=(AM**2-SQ2)/4./PM**2*((4.-SQ2/PM**2)*FA**2-(4.+SQ2/PM**2)*
C      &     F1V**2-SQ2/PM**2*F2VG**2*(1.+SQ2/4./PM**2)
C      &     -4.*SQ2*F1V*F2VG/PM**2
C      &     -AM**2/PM**2*((F1V+F2VG)**2+FA**2))
      AAA=(AM**2-SQ2)/4./PM**2*((4.-SQ2/PM**2)*FA**2-(4.+SQ2/PM**2)*
     &     F1V**2-SQ2/PM**2*F2VG**2*(1.+SQ2/4./PM**2)
     &     -4.*SQ2*F1V*F2VG/PM**2
     &     -AM**2/PM**2*((F1V+F2VG)**2+(FA+2.*FP)**2+(SQ2/PM**2-4.)*FP**2))
      BBB=-SQ2/PM**2*FA*(F1V+F2VG)
      CCC=(FA**2+F1V**2-SQ2/PM**2*F2VG**2/4.)/4.
      SMU=4.*PM*E+SQ2-AM**2
     
      IF (IPAR.GT.0) THEN
         DNELSQ2=CONST/E**2*(AAA-BBB*SMU/PM**2+CCC*SMU**2/PM**4)
      ELSE
         DNELSQ2=CONST/E**2*(AAA+BBB*SMU/PM**2+CCC*SMU**2/PM**4)
      ENDIF
      
      RETURN
      END
