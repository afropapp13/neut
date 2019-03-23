************************************************************************
*     ----------------------------------------
      SUBROUTINE NECOHVCT(IPAR,ENEU,IMODE,DIR)
*     ----------------------------------------
*
*     (Purpose)
*       Coherent PI0 Production
*
*     (Input)
*       IPAR  : INCIDENT NEUTRINO ID
*       ENEU  : INCIDENT NEUTRINO ENERGY
*       IMODE : INTERACTION MODE
*                +-16 : C.C.
*                +-36 : N.C.
*       DIR(3): INCIDENT NEUTRINO DIRECTION
*
*     (Output)
*
*     (Creation Date and Author)
*       1995.07.01 ; G.Dion & K.Okumura
*       1996.02.01 ; K.Kaneyuki   COHVCT -> NECOHVCT
*       2006.06.12 ; Y.Hayato to use IBOUND instaed of INO16,XIN
*       2007.01.08 ; G.Mitsuka, add Kartavtsev model and fix some bugs
*       2007.11.05 ; G.Mitsuka, make consistency for some functions
*
************************************************************************
      IMPLICIT NONE
      
      integer*4 IPAR,IMODE
      real*4    ENEU
      REAL*4    DIR(3)

      real*4    dcohnerein,dcohnekart,dcohmxnerein,dcohmxnekart,
     & dcohmxneberger,dcohneberger,rlu
      external  dcohnerein,dcohnekart,dcohmxnerein,dcohmxnekart,
     & dcohmxneberger,dcohneberger,rlu

#include "nework.h"
#include "necard.h"
#include "vcwork.h"
#include "posinnuc.h"
#include "neutcrs.h"

C      COMMON /SAVXIN/INO16,XIN(3)

      INTEGER IORG(5)/ 0, 0, 1, 2, 2/
      INTEGER IFLG(5)/-1,-1, 0, 5, 0/
      INTEGER ICRN(5)/ 0, 0, 1, 0, 1/
C     
C      COMMON /NVECT/ IMOD,LVECT,IP(20),PNE(3,20)
C
C     IMOD      : INTERACTION MODE ( = IMODE )
C     LVECT     : NUMBER OF PARTICLES ( = 5 )
C     IP(20)    : PARTICLE INDENTITY
C                       IP(1)  INCIDENT NEUTRINO
C                       IP(2)  O2
C                       IP(3)  OUTGOING LEPTON
C                       IP(4)  O2
C                       IP(5)  PION
C     PNE(3,20) : PARTICLE MOMENTUM
C
C
      REAL*4 pi
      PARAMETER (     PI      = 3.141592653589793D0 )

      REAL*4 XMAX,XMIN
      PARAMETER (     XMAX = 1.     , XMIN = 0.     )

      REAL*4 ZMAX,ZMIN
      PARAMETER (     ZMAX = 1.     , ZMIN = -1.    )

      REAL*4 PHIMAX,PHIMIN
      PARAMETER (     PHIMAX = 2*PI , PHIMIN = 0.   )

      REAL*4 PSIMAX,PSIMIN
      PARAMETER (     PSIMAX = 2*PI , PSIMIN = 0.   )

      real*4 NuclMass

      real*4    dirabs,dummy
      integer*4 isign
      
      real*4    pimass,om,prmass
      real*4    f,fmax,fvalue,ymin,ymax,x,y,z,phi,psi
      real*4    e2,p3abs,w,sa,sp,cp,st,xyabs,epi,p5abs

      integer*4 i

C     !---> Statements and initialization for Kartavtsev model
      integer num,num1,num2
      real*4 Q2,t,cohkartq2min,s,Wmin,Q2min,Q2max,Mpr,Mnt
      external cohkartq2min
      parameter (Q2max = 0.2)
      real*4 numinbuf,numin,numax,cohkartwmax,Wmax,xi,nu
      parameter (xi = 3.)
      real*4 cohkarttmin,tmin,tmax
      external cohkarttmin
      parameter (tmax = 0.05)
C     <---!

C
C     INPUT DATA
C
C      INO16=0
      IBOUND=0

      MODENE = IMODE
      NUMNE = 5
      IPNE(1) = IPAR
C GEANT O16 : 69
C      IPNE(2) = 1000000
C      IPNE(4) = 1000000
C--Modify on Feb.24, 97 --- Hayato
C      IPNE(2) = 69
C      IPNE(4) = 69
C--Modify on Feb.25, 97 --- Hayato
C      IPNE(2) = 10069
C      IPNE(4) = 10069
C--Modify on Mar.16, 98 --- Hayato
C      IPNE(2) = 100069
C      IPNE(4) = 100069
C--Modify on Oct.05, 09 --- Hayato
      IPNE(2) = 1000000000 + NUMBNDP*10000 + NUMATOM * 10
      IPNE(4) = IPNE(2)

      DIRABS = SQRT (DIR(1)**2 + DIR(2)**2 + DIR(3)**2)
      PNE(1,1) = ENEU*DIR(1)/DIRABS
      PNE(2,1) = ENEU*DIR(2)/DIRABS
      PNE(3,1) = ENEU*DIR(3)/DIRABS
C
C     DETERMINATION OF OUTGOING LEPTON AND PION
C
      IF (ABS(IMODE).EQ.16) THEN
         ISIGN = IPNE(1)/ABS(IPNE(1))
C FIX June 3/ Y.H
C         IPNE(3) = ISIGN * ABS(IPNE(1)-1)
         IPNE(3) = ISIGN * (ABS(IPNE(1))-1)
         IPNE(5) = ISIGN * 211
      ELSE IF (ABS(IMODE).EQ.36) THEN
         IPNE(3) = IPNE(1)
         IPNE(5) = 111
      ELSE
         RETURN
      ENDIF

      CALL MCMASSGV(IPNE(5),PIMASS)
      CALL MCMASSGV(IPNE(3),OM)
      CALL MCMASSGV(2212,Mpr)
      CALL MCMASSGV(2112,Mnt)
C      NuclMass = 8*(Mpr+Mnt)
      NuclMass = REAL(NUMATOM)*(Mpr+Mnt)
C
C     MONTE-CARLO PART
C
      CALL MCMASSGV(2112,PRMASS)
C============================
C     Rein & Sehgal model
C============================
      if(NECOHEPI.eq.0) then
         FMAX = DCOHMXNEREIN(ENEU)
 10      CONTINUE
         YMIN = PIMASS/ENEU
         X = XMIN+(XMAX-XMIN)*RLU(DUMMY)
         YMAX = 1./(1.+PRMASS*X/(2.*ENEU))
         Y = YMIN+(YMAX-YMIN)*RLU(DUMMY)
         Z = ZMIN+(ZMAX-ZMIN)*RLU(DUMMY)
         PHI = PHIMIN+(PHIMAX-PHIMIN)*RLU(DUMMY)
         PSI = PSIMIN+(PSIMAX-PSIMIN)*RLU(DUMMY)
         F = FMAX*RLU(DUMMY)
         FVALUE  = DCOHNEREIN(ENEU,IPAR,IMODE,X,Y,Z,PHI)
C         write(6,*) ' necohvct : dcohne = ',fvalue,FMAX,eneu,x,y,z,phi
         IF (F .GE. FVALUE ) GOTO 10
C      write(6,*) ' necohvct : dcohne = ',fvalue,eneu,x,y,z,phi      
C
C     OUTGING LEPTON MOMENTUM CALCULATION PART
C
         E2 = ENEU*(1.-Y)
         P3ABS = SQRT(ABS(E2**2-OM**2))

         Q2 = 2*PRMASS*Eneu*Y*X
         if(abs(modene).eq.16) then
            if ((2*ENEU*SQRT(ENEU**2*(1-Y)**2-OM**2)).ne.0) then
CC(G.Mitsuka, Jan 2006) 
C                W = 1.-(OM**2+Q2)/(2*ENEU*E2)
CC modified by P.Martins
               W =  (2*ENEU**2*(1-Y)-Q2-OM**2)/
     &              (2*ENEU*SQRT(ENEU**2*(1-Y)**2-OM**2))
            else
               goto 10
            endif
         elseif(abs(modene).eq.36) then ! NC
            W = 1.-PRMASS*X*Y/E2            
         endif
C
         SA = SQRT(ABS(1.-W*W))
         CP = COS(PSI)
         SP = SIN(PSI)
         XYABS = SQRT(PNE(1,1)**2+PNE(2,1)**2)
         IF ( XYABS.EQ.0. ) THEN
            PNE(1,3) = P3ABS*SA*SP
            PNE(2,3) = P3ABS*SA*CP
            PNE(3,3) = P3ABS*W
         ELSE
            PNE(1,3) = P3ABS*(W*PNE(1,1)/ENEU+
     &           SA*(PNE(1,1)*PNE(3,1)*CP/ENEU-PNE(2,1)*SP)/XYABS)
            PNE(2,3) = P3ABS*(W*PNE(2,1)/ENEU+
     &           SA*(PNE(2,1)*PNE(3,1)*CP/ENEU+PNE(1,1)*SP)/XYABS)
            PNE(3,3) = P3ABS*(W*PNE(3,1)/ENEU-SA*XYABS*CP/ENEU)
         ENDIF
C
C     PION MOMENTUM CALCULATION PART
         EPI = ENEU*Y
         P5ABS = SQRT(ABS(EPI**2-PIMASS**2))
         
         ST = SQRT(ABS(1.-Z*Z))
         CP = COS(PHI+PSI)
         SP = SIN(PHI+PSI)
         IF ( XYABS.EQ.0. ) THEN
            PNE(1,5) = P5ABS*ST*SP
            PNE(2,5) = P5ABS*ST*CP
            PNE(3,5) = P5ABS*Z
         ELSE
            PNE(1,5) = P5ABS*(Z*PNE(1,1)/ENEU+
     &           ST*(PNE(1,1)*PNE(3,1)*CP/ENEU-PNE(2,1)*SP)/XYABS)
            PNE(2,5) = P5ABS*(Z*PNE(2,1)/ENEU+
     &           ST*(PNE(2,1)*PNE(3,1)*CP/ENEU+PNE(1,1)*SP)/XYABS)
            PNE(3,5) = P5ABS*(Z*PNE(3,1)/ENEU-ST*XYABS*CP/ENEU)
         ENDIF
C============================
C     Berger & Sehgal model
C============================
      elseif(NECOHEPI.eq.2) then
         FMAX = DCOHMXNEBERGER(ENEU)
 12      CONTINUE
         YMIN = PIMASS/ENEU
         X = XMIN+(XMAX-XMIN)*RLU(DUMMY)
         YMAX = 1./(1.+PRMASS*X/(2.*ENEU))
         Y = YMIN+(YMAX-YMIN)*RLU(DUMMY)
         Z = ZMIN+(ZMAX-ZMIN)*RLU(DUMMY)
         PHI = PHIMIN+(PHIMAX-PHIMIN)*RLU(DUMMY)
         PSI = PSIMIN+(PSIMAX-PSIMIN)*RLU(DUMMY)
         F = FMAX*RLU(DUMMY)
         FVALUE  = DCOHNEBERGER(ENEU,IPAR,IMODE,X,Y,Z,PHI)
C         write(6,*) ' necohvct : dcohne = ',fvalue,FMAX,eneu,x,y,z,phi
         IF (F .GE. FVALUE ) GOTO 12
C      write(6,*) ' necohvct : dcohne = ',fvalue,eneu,x,y,z,phi      
C
C     OUTGING LEPTON MOMENTUM CALCULATION PART
C
         E2 = ENEU*(1.-Y)
         P3ABS = SQRT(ABS(E2**2-OM**2))

         Q2 = 2*PRMASS*Eneu*Y*X
         if(abs(modene).eq.16) then
            if ((2*ENEU*SQRT(ENEU**2*(1-Y)**2-OM**2)).ne.0) then
CC(G.Mitsuka, Jan 2006) 
C                W = 1.-(OM**2+Q2)/(2*ENEU*E2)
CC modified by P.Martins
               W =  (2*ENEU**2*(1-Y)-Q2-OM**2)/
     &              (2*ENEU*SQRT(ENEU**2*(1-Y)**2-OM**2))
            else
               goto 12
            endif
         elseif(abs(modene).eq.36) then ! NC
            W = 1.-PRMASS*X*Y/E2            
         endif
C
         SA = SQRT(ABS(1.-W*W))
         CP = COS(PSI)
         SP = SIN(PSI)
         XYABS = SQRT(PNE(1,1)**2+PNE(2,1)**2)
         IF ( XYABS.EQ.0. ) THEN
            PNE(1,3) = P3ABS*SA*SP
            PNE(2,3) = P3ABS*SA*CP
            PNE(3,3) = P3ABS*W
         ELSE
            PNE(1,3) = P3ABS*(W*PNE(1,1)/ENEU+
     &           SA*(PNE(1,1)*PNE(3,1)*CP/ENEU-PNE(2,1)*SP)/XYABS)
            PNE(2,3) = P3ABS*(W*PNE(2,1)/ENEU+
     &           SA*(PNE(2,1)*PNE(3,1)*CP/ENEU+PNE(1,1)*SP)/XYABS)
            PNE(3,3) = P3ABS*(W*PNE(3,1)/ENEU-SA*XYABS*CP/ENEU)
         ENDIF
C
C     PION MOMENTUM CALCULATION PART
         EPI = ENEU*Y
         P5ABS = SQRT(ABS(EPI**2-PIMASS**2))
         
         ST = SQRT(ABS(1.-Z*Z))
         CP = COS(PHI)
         SP = SIN(PHI)
         IF ( XYABS.EQ.0. ) THEN
            PNE(1,5) = P5ABS*ST*SP
            PNE(2,5) = P5ABS*ST*CP
            PNE(3,5) = P5ABS*Z
         ELSE
            PNE(1,5) = P5ABS*(Z*PNE(1,1)/ENEU+
     &           ST*(PNE(1,1)*PNE(3,1)*CP/ENEU-PNE(2,1)*SP)/XYABS)
            PNE(2,5) = P5ABS*(Z*PNE(2,1)/ENEU+
     &           ST*(PNE(2,1)*PNE(3,1)*CP/ENEU+PNE(1,1)*SP)/XYABS)
            PNE(3,5) = P5ABS*(Z*PNE(3,1)/ENEU-ST*XYABS*CP/ENEU)
         ENDIF

C================================
C     Kartavtsev et al. model
C================================
      else
         FMAX = DCOHMXNEKART(ENEU)

C     Integration range over Q2,nu,t are limited
         s = NuclMass**2+2*NuclMass*ENEU
         Wmin = (NuclMass+PIMASS)**2
         Q2min = cohkartq2min(s,Wmin)
         if(Q2min.gt.Q2max.or.Q2min.eq.1.) then
            print*, 'Starange Q2min value in Coherent-pi mode'
            stop
         elseif(Q2min.gt.-1..and.Q2min.lt.0.) then
            Q2min = 0.
         endif
         
C-------------------------------------||
C     Start of kinematics calculation \/
C     Make random parameter inside the allowed range
 11      continue
         X = XMIN+(XMAX-XMIN)*RLU(DUMMY)
         YMIN = Q2min/(2*NuclMass*ENEU*X)
         YMAX = Q2max/(2*NuclMass*ENEU*X)
         Y = YMIN+(YMAX-YMIN)*RLU(DUMMY)
         Z = ZMIN+(ZMAX-ZMIN)*RLU(DUMMY)
         PHI = PHIMIN+(PHIMAX-PHIMIN)*RLU(DUMMY)
         PSI = PSIMIN+(PSIMAX-PSIMIN)*RLU(DUMMY)

         Q2 = 2*NuclMass*Eneu*Y*X ! Q2

         num = num + 1
         if((num-(num/1000000)*1000000).eq.1)
     &        write(*,'("Q2 ",7E10.3)') Q2min,Q2,Q2max,ENEU,x,y
         if((num-(num/10000)*10000).eq.1) then
            num = 1
            goto 11
         endif

C     If Q2 is out of range, go to start point
         if(Q2.lt.Q2min.or.Q2.gt.Q2max) goto 11

         nu = Eneu*y ! nu
         numinbuf = (Wmin+Q2-NuclMass**2)/(2*NuclMass)
         numin=amax1(xi*sqrt(Q2),numinbuf)
         Wmax = cohkartwmax(Q2,s)
         numax=(Wmax+Q2-NuclMass**2)/(2*NuclMass)
C     If nu is out of range, go to start point
         if(Wmax.lt.0..or.numax.lt.0..or.
     &        nu.lt.numin.or.nu.gt.numax) goto 11

         tmin = abs(cohkarttmin(Q2,nu)) 
C     If tmin is out of range, go to start point
        if(tmin.lt.0..or.tmin.gt.tmax) goto 11

C
C     OUTGING LEPTON MOMENTUM CALCULATION PART
C
         E2 = ENEU*(1.-Y)
         P3ABS = SQRT(ABS(E2**2-OM**2))
            W = (2.*ENEU*E2-OM**2-Q2)/2./ENEU/P3ABS
         SA = SQRT(ABS(1.-W*W))
         CP = COS(PSI)
         SP = SIN(PSI)
         XYABS = SQRT(PNE(1,1)**2+PNE(2,1)**2)
         IF ( XYABS.EQ.0. ) THEN
            PNE(1,3) = P3ABS*SA*SP
            PNE(2,3) = P3ABS*SA*CP
            PNE(3,3) = P3ABS*W
         ELSE
            PNE(1,3) = P3ABS*(W*PNE(1,1)/ENEU+
     &           SA*(PNE(1,1)*PNE(3,1)*CP/ENEU-PNE(2,1)*SP)/XYABS)
            PNE(2,3) = P3ABS*(W*PNE(2,1)/ENEU+
     &           SA*(PNE(2,1)*PNE(3,1)*CP/ENEU+PNE(1,1)*SP)/XYABS)
            PNE(3,3) = P3ABS*(W*PNE(3,1)/ENEU-SA*XYABS*CP/ENEU)
         ENDIF
C
C     PION MOMENTUM CALCULATION PART
         EPI = ENEU*Y
         P5ABS = SQRT(ABS(EPI**2-PIMASS**2))
         
         ST = SQRT(ABS(1.-Z*Z))
         CP = COS(PHI+PSI)
         SP = SIN(PHI+PSI)
         IF ( XYABS.EQ.0. ) THEN
            PNE(1,5) = P5ABS*ST*SP
            PNE(2,5) = P5ABS*ST*CP
            PNE(3,5) = P5ABS*Z
         ELSE
            PNE(1,5) = P5ABS*(Z*PNE(1,1)/ENEU+
     &           ST*(PNE(1,1)*PNE(3,1)*CP/ENEU-PNE(2,1)*SP)/XYABS)
            PNE(2,5) = P5ABS*(Z*PNE(2,1)/ENEU+
     &           ST*(PNE(2,1)*PNE(3,1)*CP/ENEU+PNE(1,1)*SP)/XYABS)
            PNE(3,5) = P5ABS*(Z*PNE(3,1)/ENEU-ST*XYABS*CP/ENEU)
         ENDIF
         
         t = abs(-(PNE(1,1)-PNE(1,3)-PNE(1,5))**2 ! t
     &        -(PNE(2,1)-PNE(2,3)-PNE(2,5))**2
     &        -(PNE(3,1)-PNE(3,3)-PNE(3,5))**2)
         num1 = num1 + 1
         num2 = num2 + 1
         if((num1-(num1/100)*100).eq.1) then
            print*, 't ',tmin,t,tmax
            num1 = 1
            goto 11
         elseif((num2-(num2/1000)*1000).eq.1) then
            print*, 't ',tmin,t,tmax
            num2 = 1
            goto 11
         endif

C     If t is out of range, go to start point
         if(t.lt.tmin.or.t.gt.tmax) goto 11

         num = 0
         num1 = 0
         num2 = 0

         F = FMAX*RLU(DUMMY)
C     Calculate the differential cross section with given parameters
         FVALUE  = DCOHNEKART(ENEU,Q2,Epi,t,x,y,z)
         IF (F .GE. FVALUE ) GOTO 11

C     End of kinematics calculation /\
C-----------------------------------||
      endif ! End of model selection

C     Fill common block with kinematic variables
      CRSX = X
      CRSY = Y
      CRSZ = Z
      CRSPHI = PHI
C
C --  FILL FLAGS
C
      DO 100 I = 1,5
         IORGNE(I)=IORG(I)
         IFLGNE(I)=IFLG(I)
         ICRNNE(I)=ICRN(I)
 100  CONTINUE
      
      RETURN
      END
