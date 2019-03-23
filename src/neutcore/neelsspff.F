************************************************************************
*     --------------------------------------------
      SUBROUTINE NEELSSPFF(E,IPAR,PFSEL,ICALL,IERR)
*     --------------------------------------------
*
*     (Purpose)
*       Set fermi momentu of nucleon for elastic scattering
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
*       ICALL  : FLAG TO MAKE MAP
*                   0 : MAKE MAP
*     (Output)
*       PFSEL(3) : FERMI MOMENTUM ( GEV )
*
*     (Creation Date and Author)
*       19??.??.?? ; M.NAKAHATA
*       1987.08.14 ; N.SATO     
*                    ELSPF2 -> ELSPF ( SSL2 VR )
*     1995.02.06 ; K. KANEYUKI  FOR S.K.
*                    ELSSPF -> NEELSSPF
*
*     1998.02.24 ; K. Ishihara  sets PFERMI 217 MEV
*     
************************************************************************
      IMPLICIT NONE

#include "neutparams.h"
      
      real e,PFSEL(3)
      integer*4 ipar,icall,ierr

      integer*4 i,j

      real*4 ami

      REAL*4 XMAP(7,6),PF(3),pfabs,pneu(3)
      real*4 beta,gm,ev(3)
      real*4 cost,sint
      real*4 enuc,enstop

      integer*4 np,nc

      real*4 crs,ferm,dpf,dco,dat1,dat2,distri
      real*4 xmax,hit,pfa,coa

      real*4 dum
      real*4 fnelscrs,fvfrmgev,rlu
      external fnelscrs,fvfrmgev,rlu
C     
      IERR=0

      IF (IPAR.gt.0) then
         CALL MCMASS(2112,AMI)
      ELSE
         CALL MCMASS(2212,AMI)
      ENDIF         

      AMI = AMI * 1.e-3

      IF (ICALL.EQ.0) THEN
C     
C     -- MAKE FERMI MOMENTUM -- DIRECTION MAP
C     
         DO 110 J=1,6
            XMAP(1,J)=0.
            XMAP(7,J)=0.
 110     CONTINUE

         DO 120 I=1,5
            PFABS=FLOAT(I)*(PFMAX/5.0)
            DO 130 J=1,6
               COST=-1.+0.4*FLOAT(J-1)
C               IF (abs(COST).eq.1.) then
C     - Bug fix for gfortran (Y. Hayato, July 12, 2011)
               IF ((J.eq.1).or.(J.eq.6)) THEN
                  SINT = 0.
               ELSE
                  SINT=SQRT(1.-COST**2)
               endif
               PF(1)=PFABS*COST
               PF(2)=PFABS*SINT
               PF(3)=0.

C-- calc. enstop
               PNEU(1)=E
               PNEU(2)=0.
               PNEU(3)=0.

               ENUC=SQRT(PFABS**2+AMI**2)
               BETA=PFABS/ENUC
               GM=1./SQRT(1.-BETA**2)
               IF (PFABS.NE.0.) THEN
                  EV(1)=-PF(1)/PFABS
                  EV(2)=-PF(2)/PFABS
                  EV(3)=-PF(3)/PFABS
               ELSE
                  EV(1)=1.
                  EV(2)=0.
                  EV(3)=0.
               ENDIF
               CALL MCVECBST(PNEU,0.,EV,GM)
               ENSTOP=SQRT(PNEU(1)**2+PNEU(2)**2+PNEU(3)**2)

C               CRS=FNELSPF(E,PF,IPAR)
               CRS=FNELSCRS(ENSTOP,IPAR)
               FERM=FVFRMGEV(PFABS)
               XMAP(I+1,J)=CRS*FERM
 130        CONTINUE
 120     CONTINUE
C     
C --  SET FERMI MOMENTUM ON THIS MAP
C     
         XMAX=0.
         DO 300 I=2,6
            DO 310 J=1,6
               IF(XMAX.GT.XMAP(I,J))GO TO 310
               XMAX=XMAP(I,J)
 310        CONTINUE
 300     CONTINUE
         IF(XMAX.EQ.0.) THEN
C     
C --        ERROR RETURN
C     
            WRITE(6,900)
 900        FORMAT(' *** ERROR AT ELSSPF (MAX OF MAP=0)')
            IERR=1
            RETURN
         ENDIF
      ENDIF
      
 320  PFABS=RLU(DUM)*PFMAX
      COST=-1.+RLU(DUM)*2.
      HIT=XMAX*RLU(DUM)
      PFA=PFABS/(PFMAX/5.0)
      COA=(COST+1.)/0.4
      NP=IFIX(PFA)+1
      NC=IFIX(COA)+1
      DPF=PFA+1.-FLOAT(NP)
      DCO=COA+1.-FLOAT(NC)
      DAT1=(1.-DPF)*XMAP(NP,NC)+DPF*XMAP(NP+1,NC)
      if (COST.gt.0.99999999) then
         DAT2=(1.-DPF)*XMAP(NP,NC+1)+DPF*XMAP(NP+1,NC+1)
         DISTRI=(1.-DCO)*DAT1+DCO*DAT2
      else
         DISTRI=DAT1
      endif
      IF(HIT.GT.DISTRI)GO TO 320

      if (abs(cost).eq.1.) then
         sint=0.
      else
         SINT=SQRT(1.-COST**2)
      endif
      PFSEL(1)=PFABS*COST
      PFSEL(2)=PFABS*SINT
      PFSEL(3)=0.

      RETURN
      END
