************************************************************************
*     ----------------------------------
      real*8 FUNCTION dnncelsq2(E,IPAR,ITARG,Q2)
*     ----------------------------------
*
*     (Purpose)
*       Neutral Current elastic without pauli principle 
*       nu + n -> nu + n
*       nu + p -> nu + p
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
*       ITARG  : Target type
*                2212 : Proton
*                2112 : Proton
*       Q2     : Q**2
*     
*     (Output)
*       DNNCELSQ2 : Q**2 DISTRIBUTION OF NC ELASTIC SCATTERING
*
*     (Creation Date and Author)
*      2011.02.23 ; T.Mori   ( ncn,ncp )
*      2011.03.11 ; Y.Hayato ( combine two )
************************************************************************
C
C    CONST ; G**2*M**2/8/PI
C    WEINBG ; Weingberg angle
C    AM ; lepton mass(=0.)  
C
      IMPLICIT NONE

#include "neutmodel.h"
#include "rscons.h"

      REAL*8    E,Q2
      INTEGER*4 IPAR,ITARG

      REAL FNNCQ2MAX,FNNCQ2MIN,EE
      EXTERNAL FNNCQ2MAX, FNNCQ2MIN

      REAL faBBBA07
      EXTERNAL faBBBA07

      REAL*8 CONST,PM,PIM,WEINBG
      REAL*8 AM

C--   taken from A. Bodek et al, Eur. Phys. J. C 71 (2011) 1726 
      REAL TRCORA, TRCORB
      PARAMETER (TRCORA = 6.0)
      PARAMETER (TRCORB = 0.34)

c      DATA CONST/1.76D-1/      
      DATA CONST/1.85D-1/      !<------ Remove Cabbibo angle
      DATA PM/0.938/
      DATA PIM/0.13957/
   
      DATA WEINBG/0.23117/

      REAL*8 AVMASS,VVMASS
      REAL*8 Q2T1,Q2T2,SQ2
      REAL*8 GPE,GPM,GNE,GNM,GVE,GVM
      REAL*8 F1V,F2VG,FA,FP
      REAL*8 FN1,FN2
      REAL*8 AAA,BBB,CCC,SMU
      REAL*8 F3A,F3V

      REAL*8 ELAXFF
      EXTERNAL ELAXFF
      
      INTEGER SGN
      LOGICAL NC
      
      INTEGER NCELMODE

      AVMASS = XMANCEL
      VVMASS = XMVQE
     
      AM = 0.

C --  Setup info for FA CALC      
      NC = .true.
C      if (SGN.lt.0) then
         SGN = -1
C      else
C         SGN = 1
C      endif
      
C-------
      EE=E
      Q2T1=FNNCQ2MIN(EE,ITARG)
      IF(Q2.LT.Q2T1) THEN
         dnncelsq2=0.
         RETURN
      END IF
      Q2T2=FNNCQ2MAX(EE,ITARG)
      IF(Q2.GT.Q2T2) THEN
         dnncelsq2=0.
         RETURN
      END IF
C-------
      SQ2=-Q2

      NCELMODE = mod(MDLQE,100)/10

      if (NCELMODE.eq.3) THEN           
         CALL BBBA07(Q2,GPE,GPM,GNE,GNM)
      else if (NCELMODE.eq.2) THEN
         CALL BBBA05(Q2,GPE,GPM,GNE,GNM)
      else if (NCELMODE.eq.1) THEN 
         GPE = 1./(1.- SQ2/VVMASS**2)**2
         GPM = 2.793/(1.- SQ2/VVMASS**2)**2
         
         GNE = 0.
         GNM = -1.913/(1.-SQ2/VVMASS**2)**2
      else
         write(*,*) 'Invalid NC Elastic mode specified by MDLQE',MDLQE
         stop
      endif
         
      GVE=GPE-GNE
      GVM=GPM-GNM

      if (mod(MDLQE, 1000)/100.eq.1) THEN
         GVM = GVM * sqrt(1+TRCORA * Q2 * exp( SQ2 / TRCORB ))
      endif

      F1V=(GVE-SQ2*GVM/4./PM**2)/(1.-SQ2/4./PM**2)
      F2VG=(GVM-GVE)/(1.-SQ2/4./PM**2)
      F3V = SCCFV*(-1.) /(1.-SQ2/VVMASS**2)**2

C      FA=-0.5*1.23/(1.-SQ2/AVMASS**2)**2

C      FA=0.5*1.267/(1.-SQ2/AVMASS**2)**2

      FA=ELAXFF(Q2, NC, SGN)
      
C-- Proton, flip sign
      if (ITARG.eq.2212) then
         FA=-1. * FA
      endif

      FP=FPQE*2.*(PM**2)*FA/(-SQ2+PIM**2)
      F3A = SCCFA * FA 
   
C-- Proton or neutron
      if (ITARG.eq.2112) then
         FN1 = (GNE-SQ2/4./PM**2*GNM)/(1.-SQ2/4./PM**2)
         FN2 = (GNM-GNE)/(1.-SQ2/4./PM**2)
      else if (ITARG.eq.2212) then
         FN1 = (GPE-SQ2/4./PM**2*GPM)/(1.-SQ2/4./PM**2)
         FN2 = (GPM-GPE)/(1.-SQ2/4./PM**2)
      else
         write(*,*) 'dnncelsq2.F : ITARG is not neutron nor proton',
     $        ITARG
         stop
      endif


C-- Proton or neutron
      if (ITARG.eq.2112) then
         F1V  = -0.5*F1V - 2.*WEINBG*FN1
         F2VG = -0.5*F2VG - 2.*WEINBG*FN2
      else if (ITARG.eq.2212) then
         F1V  = 0.5*F1V - 2.*WEINBG*FN1
         F2VG = 0.5*F2VG - 2.*WEINBG*FN2
      else
         write(*,*) 'dnncelsq2.F : ITARG is not neutron nor proton',
     $        ITARG
         stop
      endif

C       AAA=(AM**2-SQ2)/4./PM**2*((4.-SQ2/PM**2)*FA**2-(4.+SQ2/PM**2)*
C      &     F1V**2-SQ2/PM**2*F2VG**2*(1.+SQ2/4./PM**2)
C      &     -4.*SQ2*F1V*F2VG/PM**2
C      &     -AM**2/PM**2*((F1V+F2VG)**2+FA**2))
C      AAA=(AM**2-SQ2)/4./PM**2*((4.-SQ2/PM**2)*FA**2-(4.+SQ2/PM**2)*
C     &     F1V**2-SQ2/PM**2*F2VG**2*(1.+SQ2/4./PM**2)
C     &     -4.*SQ2*F1V*F2VG/PM**2 -AM**2/PM**2*((F1V+F2VG)**2
C     &     +(FA+2.*FP)**2+(SQ2/PM**2-4.)*FP**2))
C      BBB=-SQ2/PM**2*FA*(F1V+F2VG)
C      CCC=(FA**2+F1V**2-SQ2/PM**2*F2VG**2/4.)/4.
      AAA=(AM**2-SQ2)/4./PM**2*((4.-SQ2/PM**2)*FA**2-(4.+SQ2/PM**2)*
     &     F1V**2-SQ2/PM**2*F2VG**2*(1.+SQ2/4./PM**2)
     &     -4.*SQ2*F1V*F2VG/PM**2+SQ2/PM**2*(4-SQ2/PM**2)*F3V*2
     &     -AM**2/PM**2*((F1V+F2VG)**2+(FA+2.*FP)**2
     &     +(SQ2/PM**2-4.)*(F3V**2+FP**2)))
      BBB=-SQ2/PM**2*FA*(F1V+F2VG)-AM**2/PM**2*(F1V+SQ2/4./PM**2*F2VG)
     &     *F3V-(FA+SQ2*FP/2./PM**2)*F3A
      CCC=(FA**2+F1V**2-SQ2/PM**2*F2VG**2/4.-SQ2/PM**2*F3A**2)/4.
      SMU=4.*PM*E+SQ2-AM**2
     
      IF (IPAR.GT.0) THEN
       dnncelsq2=CONST/E**2*(AAA-BBB*SMU/PM**2+CCC*SMU**2/PM**4)
      ELSE
       dnncelsq2=CONST/E**2*(AAA+BBB*SMU/PM**2+CCC*SMU**2/PM**4)
      ENDIF

      RETURN
      END
