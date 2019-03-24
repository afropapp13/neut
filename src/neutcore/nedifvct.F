************************************************************************
*     ----------------------------------------
      SUBROUTINE NEDIFVCT(IPAR,ENEU,IMODE,DIR)
*     ----------------------------------------
*
*     (Purpose)
*       Diffractive Pion Production on free protons
*
*     (Input)
*       IPAR  : INCIDENT NEUTRINO ID
*       ENEU  : INCIDENT NEUTRINO ENERGY
*       IMODE : INTERACTION MODE
*                +-15 : C.C.
*                +-35 : N.C.
*       DIR(3): INCIDENT NEUTRINO DIRECTION
*
*     (Output)
*
*     (Creation Date and Author)
*       2016.01.18 ; R. Terri (based on necohvct)
************************************************************************
      IMPLICIT NONE
      
      integer*4 IPAR,IMODE
      real*4    ENEU
      REAL*4    DIR(3)

      real*4    rlu
      external  rlu
      real*8    ddifnerein,ddifmx
      external  ddifnerein,ddifmx

      real*4 tmpnu(3)
      real*4 tmplep(3)
      real*4 tmppi(3)

      REAL    COSTH,SINTH,COSPH,SINPH,DX,DY,DZ

#include "nework.h"
#include "necard.h"
#include "vcwork.h"
#include "posinnuc.h"
#include "neutcrs.h"

C      COMMON /SAVXIN/INO16,XIN(3)

      INTEGER IORG(5)/ 0, 0, 1, 2, 2/
      INTEGER IFLG(5)/-1,-1, 0, 0, 0/
      INTEGER ICRN(5)/ 0, 0, 1, 1, 1/
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
      PARAMETER (     ZMAX = 0.9999999     , ZMIN = -0.9999999    )

      REAL*4 PHIMAX,PHIMIN
      PARAMETER (     PHIMAX = 2*PI , PHIMIN = 0.   )

      REAL*4 PSIMAX,PSIMIN
      PARAMETER (     PSIMAX = 2*PI , PSIMIN = 0.   )

      real*4    dirabs,dummy
      integer*4 isign
      
      real*4    pimass,om,prmass,tmax,tmin,nu,cosnul
      real*4    f,fmax,fvalue,ymin,ymax,x,y,z,phi,psi
      real*4    e2,p3abs,w,sa,sp,cp,st,xyabs,epi,p5abs
      real*4    cospiq,cosnuq,tr

      integer*4 i,count

      real*4 Q2,t,A1min,A2min,A3min,A1max,A2max,A3max
C     <---!

C
C     INPUT DATA
C
C      INO16=0

C      write(*,*) 'starting nedifct'

      count = 0

      IBOUND=0

      MODENE = IMODE
      NUMNE = 5
      IPNE(1) = IPAR
C Target & outgoing nucleon ALWAYS free proton
      IPNE(2) = 2212
      IPNE(4) = IPNE(2)

      DIRABS = SQRT (DIR(1)**2 + DIR(2)**2 + DIR(3)**2)
      PNE(1,1) = ENEU*DIR(1)/DIRABS
      PNE(2,1) = ENEU*DIR(2)/DIRABS
      PNE(3,1) = ENEU*DIR(3)/DIRABS

C Free proton is assumed to be at rest
      PNE(1,2) = 0
      PNE(2,2) = 0
      PNE(3,2) = 0


      tmpnu(1) = 0
      tmpnu(2) = 0
      tmpnu(3) = ENEU

      DZ=DIR(3)
      DX=DIR(1)
      DY=DIR(2)
      COSTH=DZ
      if (abs(dz).eq.1.) then
         sinth = 0.
      else
         SINTH=SQRT(1.-DZ**2)
      endif
      IF (ABS(DX).GT.0.000001) THEN
        COSPH=COS(ATAN2(DY,DX))
        SINPH=SIN(ATAN2(DY,DX)) 
      ELSE
        COSPH=0.
        IF (ABS(DY).GT.0.00001) THEN
          SINPH=1.*(ABS(DY)/DY)
        ELSE
          SINPH=1.
        ENDIF
      ENDIF 

C     DETERMINATION OF OUTGOING LEPTON AND PION
C
      IF (ABS(IMODE).EQ.15) THEN
         ISIGN = IPNE(1)/ABS(IPNE(1))
C FIX June 3/ Y.H
C         IPNE(3) = ISIGN * ABS(IPNE(1)-1)
         IPNE(3) = ISIGN * (ABS(IPNE(1))-1)
         IPNE(5) = ISIGN * 211
      ELSE IF (ABS(IMODE).EQ.35) THEN
         IPNE(3) = IPNE(1)
         IPNE(5) = 111
      ELSE
         RETURN
      ENDIF

      CALL MCMASSGV(IPNE(5),PIMASS)
      CALL MCMASSGV(IPNE(3),OM)
C
C     MONTE-CARLO PART
C
      CALL MCMASSGV(2212,PRMASS)
C============================
C     Rein & Sehgal model
C============================
      if(NEDIFPI.eq.0) then
         FMAX = DDIFMX(ENEU)
 10      CONTINUE
C         YMIN = PIMASS/ENEU
         count = count+1
         if (count .gt. 600000) then 
C            write(*,*) 'nedifvct: 600000 tries, screw it'
            goto 120
         endif 
         X = XMIN+(XMAX-XMIN)*RLU(DUMMY)
         YMAX = 1./(1.+PRMASS*X/(2.*ENEU))
         YMIN = (4.-PRMASS**2)/(2*PRMASS*ENEU*(1-X))
         if (YMIN .gt. YMAX) goto 10
         if (YMIN .lt. PIMASS/ENEU) YMIN = PIMASS/ENEU
         Y = YMIN+(YMAX-YMIN)*RLU(DUMMY)
C First check to see if W>2 GeV
C        W = SQRT(PRMASS**2 + 2*PRMASS*ENEU*Y*(1-X))
C        if (W .le. 2) goto 10
C     have x & y, need to get t! 
         NU = ENEU*Y
         Q2 = 2*ENEU*PRMASS*X*Y
         E2 = ENEU*(1.-Y)
         if (E2 .LE. OM) goto 10
         P3ABS = SQRT(ABS(E2**2-OM**2))

C     sanity check: is cos theta_nulep in proper range?
         cosnul = (E2 - (Q2+OM**2)/(2*ENEU))/P3ABS

         if (ABS(cosnul) .gt. 1) goto 10

         cosnuq = (ENEU - P3ABS*cosnul)/
     $        sqrt(ENEU*ENEU+P3ABS*P3ABS-2*ENEU*P3ABS*cosnul)

         if (ABS(cosnuq) .GE. 1) goto 10


C     NOW we calculate T (should be negative!)
         A1min = 1+(2*NU/PRMASS)+NU*NU/(PRMASS*PRMASS)
     $        -(NU*NU+Q2)/(PRMASS*PRMASS)
         A2min = (1+(NU/PRMASS))*(PIMASS*PIMASS-Q2-2*NU*NU)
     $        +2*NU*(NU*NU+Q2)/PRMASS
         A3min = (PIMASS*PIMASS-Q2-2*NU*NU)*(PIMASS*PIMASS-Q2-2*NU*NU)
     $        -4*(NU*NU+Q2)*(NU*NU-PIMASS*PIMASS)

         A1max = 1+(2*NU/PRMASS)+NU*NU/(PRMASS*PRMASS)
         A2max = (1+(NU/PRMASS))*(PIMASS*PIMASS-Q2-2*NU*NU)
         A3max = (PIMASS*PIMASS-Q2-2*NU*NU)*(PIMASS*PIMASS-Q2-2*NU*NU)

         if (A2min*A2min-A1min*A3min .lt. 0) GOTO 10
         if (A2max*A2max-A1max*A3max .lt. 0) GOTO 10
         tmin = (A2min+sqrt(A2min*A2min-A1min*A3min))/A1min
         tmax = (A2max+sqrt(A2max*A2max-A1max*A3max))/A1max
         if (tmax>tmin) GOTO 10

         T = TMIN+(TMAX-TMIN)*RLU(DUMMY)

         EPI = ENEU*Y + T/(2*PRMASS)
         if (EPI .LE. PIMASS) goto 10
         P5ABS = SQRT(ABS(EPI**2-PIMASS**2))

         if (P5ABS.eq.0) THEN
            P5ABS = 0.0000001
         endif

C check if angle between pion & q is in allowed range
         if ((ENEU**2 + P3ABS**2 - 
     $        2*ENEU*P3ABS*cosnul) .lt. 0 .or.
     $        P5ABS .le. 0) then
            goto 10
         ENDIF

C check if angle between pion & q is in allowed range
         cospiq = (T-PIMASS**2+Q2+2*EPI*ENEU*Y)/(2*P5ABS*
     $        SQRT(ENEU**2 + P3ABS**2 - 
     $        2*ENEU*P3ABS*cosnul))
         if (ABS(cospiq) .ge. 1) GOTO 10

         F = FMAX*RLU(DUMMY)
         FVALUE  = DDIFNEREIN(ENEU,IPAR,IMODE,X,Y,T)
C         write(6,*) ' nedifvct : ddifne = ',fvalue,FMAX,eneu,x,y,z,phi
         IF (F .GE. FVALUE ) GOTO 10
C      write(6,*) ' nedifvct : ddifne = ',fvalue,eneu,x,y,z,phi      
C
C     OUTGING LEPTON MOMENTUM CALCULATION PART
C

C
C        START GETTING KINEMATICS & directions
         Z = ZMIN+(ZMAX-ZMIN)*RLU(DUMMY)
         if (Z.eq.1) then
            z = 0.999999
         else if (z.eq.-1) then
            z = -0.999999
         endif 
C         PHI = PHIMIN+(PHIMAX-PHIMIN)*RLU(DUMMY)
         PSI = PSIMIN+(PSIMAX-PSIMIN)*RLU(DUMMY)

         SA = SQRT(1-cosnul**2)
         CP = COS(PSI)
         SP = SIN(PSI)
         tmplep(1) = P3ABS*SA*CP
         tmplep(2) = P3ABS*SA*SP
         tmplep(3) = P3ABS*cosnul
C
C     PION MOMENTUM CALCULATION PART

         if (ABS((cospiq-Z*cosnuq)/(SQRT(ABS(1-Z*Z))*SQRT(1-cosnuq**2))) 
     $        .gt. 1) GOTO 10


C To allow all possible phis:
         PHI = ACOS((cospiq-Z*cosnuq)/(SQRT(ABS(1-Z*Z))
     $        *SQRT(ABS(1-cosnuq**2))))

         if (RLU(DUMMY) .GT. 0.5) then 
            PHI = PHI+PI
         endif

         ST = SQRT(ABS(1.-Z*Z))
         CP = COS(PHI+PSI)
         SP = SIN(PHI+PSI)
         tmppi(1) = P5ABS*ST*CP
         tmppi(2) = P5ABS*ST*SP
         tmppi(3) = P5ABS*Z
      else if (NEDIFPI .eq. 1) then 
         FMAX = DDIFMX(ENEU)
C         write(*,*) 'nedifcvt: Enu = ',ENEU,' and FMAX=',FMAX
 20      CONTINUE
         count = count+1
         if (count .gt. 600000) then 
C            write(*,*) 'nedifvct: 600000 tries, screw it'
            goto 120
         endif 
         YMIN = PIMASS/ENEU
         X = XMIN+(XMAX-XMIN)*RLU(DUMMY)
         YMAX = 1./(1.+PRMASS*X/(2.*ENEU))
         if (YMIN .gt. YMAX) then
C            if (MOD(count,2).eq.0) write(*,*) 'nedifvct: ymin>ymax'
            goto 20
         endif
         Y = YMIN+(YMAX-YMIN)*RLU(DUMMY)

C     have x & y, need to get t! 
         NU = ENEU*Y
         Q2 = 2*ENEU*PRMASS*X*Y
         E2 = ENEU*(1.-Y)
         if (E2 .LE. OM) goto 20
         P3ABS = SQRT(ABS(E2**2-OM**2))

C     sanity check: is cos theta_nulep in proper range?
         cosnul = (E2 - (Q2+OM**2)/(2*ENEU))/P3ABS

         if (ABS(cosnul) .gt. 1) then
C            if (MOD(count,2).eq.0) write(*,*) 'nedifvct: |cosnul|>1'
            goto 20
         endif
         cosnuq = (ENEU - P3ABS*cosnul)/
     $        sqrt(ENEU*ENEU+P3ABS*P3ABS-2*ENEU*P3ABS*cosnul)

         if (ABS(cosnuq) .ge. 1) then
C            if (MOD(count,2).eq.0) write(*,*) 'nedifvct: |cosnuq|>1'
            goto 20
         endif

C     NOW we calculate T (should be negative!)
         A1min = 1+(2*NU/PRMASS)+NU*NU/(PRMASS*PRMASS)
     $        -(NU*NU+Q2)/(PRMASS*PRMASS)
         A2min = (1+(NU/PRMASS))*(PIMASS*PIMASS-Q2-2*NU*NU)
     $        +2*NU*(NU*NU+Q2)/PRMASS
         A3min = (PIMASS*PIMASS-Q2-2*NU*NU)*(PIMASS*PIMASS-Q2-2*NU*NU)
     $        -4*(NU*NU+Q2)*(NU*NU-PIMASS*PIMASS)

         A1max = 1+(2*NU/PRMASS)+NU*NU/(PRMASS*PRMASS)
         A2max = (1+(NU/PRMASS))*(PIMASS*PIMASS-Q2-2*NU*NU)
         A3max = (PIMASS*PIMASS-Q2-2*NU*NU)*(PIMASS*PIMASS-Q2-2*NU*NU)

         if (A2min*A2min-A1min*A3min .lt. 0) then 
C            if (MOD(count,2).eq.0) write(*,*) 'nedifvct: tmin sqrt'
            GOTO 20
         endif
         if (A2max*A2max-A1max*A3max .lt. 0) then
C            if (MOD(count,2).eq.0) write(*,*) 'nedifvct: tmax sqrt'
            GOTO 20
         endif
         tmin = (A2min+sqrt(A2min*A2min-A1min*A3min))/A1min
         tmax = (A2max+sqrt(A2max*A2max-A1max*A3max))/A1max
         if (tmax>tmin) then
C            if (MOD(count,2).eq.0) write(*,*) 'nedifvct: tmin>tmax'
            GOTO 20
         endif
C         write(*,*) 'nedifvct: now calculating t!'

         T = TMIN+(TMAX-TMIN)*RLU(DUMMY)

         EPI = ENEU*Y + T/(2*PRMASS)
         if (EPI .LE. PIMASS) then
C            if (MOD(count,2).eq.0) write(*,*) 'nedifvct: Epi<mpi'
            goto 20
         endif
         P5ABS = SQRT(ABS(EPI**2-PIMASS**2))

         if (P5ABS.eq.0) THEN
            P5ABS = 0.0000001
         endif

C check if angle between pion & q is in allowed range
         if ((ENEU**2 + P3ABS**2 - 
     $        2*ENEU*P3ABS*cosnul) .lt. 0 .or.
     $        P5ABS .le. 0) then
            goto 20
         ENDIF

         cospiq = (T-PIMASS**2+Q2+2*EPI*ENEU*Y)/(2*P5ABS*
     $        SQRT(ENEU**2 + P3ABS**2 - 
     $        2*ENEU*P3ABS*cosnul))
         if (ABS(cospiq) .ge. 1) then
C            if (MOD(count,2).eq.0) write(*,*) 'nedifvct: |cosqpi|>1'
            GOTO 20
         endif
C         write(*,*) 'nedifct: now we can see if we have a good xsec'

         F = FMAX*RLU(DUMMY)
         FVALUE  = DDIFNEREIN(ENEU,IPAR,IMODE,X,Y,T)
         IF (F .GE. FVALUE ) then
C            if (MOD(count,2).eq.0) write(*,*) 'nedifvct: f>=fvalue'
            GOTO 20
         endif
C         write(6,*) ' nedifvct : ddif=',fvalue,FMAX,F,x,y,t
C      write(6,*) ' nedifvct : ddifne = ',fvalue,eneu,x,y,z,phi      
C
C     OUTGING LEPTON MOMENTUM CALCULATION PART
C

C
C        START GETTING KINEMATICS & directions
         Z = ZMIN+(ZMAX-ZMIN)*RLU(DUMMY)
         if (Z.eq.1) then
            z = 0.9999999
         else if (z.eq.-1) then
            z = -0.999999
         endif 
C         PHI = PHIMIN+(PHIMAX-PHIMIN)*RLU(DUMMY)
         PSI = PSIMIN+(PSIMAX-PSIMIN)*RLU(DUMMY)

         SA = SQRT(1-cosnul**2)
         CP = COS(PSI)
         SP = SIN(PSI)
         tmplep(1) = P3ABS*SA*SP
         tmplep(2) = P3ABS*SA*CP
         tmplep(3) = P3ABS*cosnul
C
C     PION MOMENTUM CALCULATION PART

         if (ABS((cospiq-Z*cosnuq)/(SQRT(ABS(1-Z*Z))*SQRT(1-cosnuq**2))) 
     $        .gt. 1) then
C            if (MOD(count,2).eq.0) write(*,*) 'nedifvct: pimom calc'
            GOTO 20
         endif

C To allow all possible phis:
         PHI = ACOS((cospiq-Z*cosnuq)/(SQRT(ABS(1-Z*Z))
     $        *SQRT(ABS(1-cosnuq**2))))

         if (RLU(DUMMY) .GT. 0.5) then 
            PHI = PHI+PI
         endif

         ST = SQRT(ABS(1.-Z*Z))
         CP = COS(PHI+PSI)
         SP = SIN(PHI+PSI)
         tmppi(1) = P5ABS*ST*SP
         tmppi(2) = P5ABS*ST*CP
         tmppi(3) = P5ABS*Z
      else
C IF any other model, get out
      CONTINUE
      write(*,*) 'nedifvct: how did you get in here?'
C-----------------------------------||
      endif ! End of model selection

C Now rotate everything back into where it should be
      CALL RSROTVEC(tmplep,PNE(1,3),COSTH,SINTH,COSPH,SINPH)
      CALL RSROTVEC(tmppi ,PNE(1,5),COSTH,SINTH,COSPH,SINPH)


C Outgoing proton kinematics!
      PNE(1,4) = PNE(1,1)-PNE(1,3)-PNE(1,5)
      PNE(2,4) = PNE(2,1)-PNE(2,3)-PNE(2,5)
      PNE(3,4) = PNE(3,1)-PNE(3,3)-PNE(3,5)

C     Fill common block with kinematic variables
      CRSX = X
      CRSY = Y
      CRST = T

C     ---- Sanity check
      TR = SQRT(PNE(1,4)*PNE(1,4)+PNE(2,4)*PNE(2,4)
     $     +PNE(3,4)*PNE(3,4)+PRMASS*PRMASS)
      if ((TR-(PRMASS-T/(2.*PRMASS)))/TR .gt. 0.000001) then
C         write(*,*) 'nedifvct: CofE issue'
         if (NEDIFPI .eq. 0) then
            goto 10
         else if (NEDIFPI .eq. 1) then
            goto 20
         endif
      endif
      
C     
C --  FILL FLAGS
C
      DO 100 I = 1,5
         IORGNE(I)=IORG(I)
         IFLGNE(I)=IFLG(I)
         ICRNNE(I)=ICRN(I)
 100  CONTINUE
 120  CONTINUE

      RETURN
      END
