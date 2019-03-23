************************************************************************
*     --------------------------------------------
      SUBROUTINE NEELSSPF(E,IPAR,PFSEL,ICALL,IERR)
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

      real*4    E,PFSEL(3)
      integer*4 IPAR,ICALL,IERR

      real*4   fnelspf,fvfrmgev,rlu
      external fnelspf,fvfrmgev,rlu


      REAL XMAP(7,6),PF(3)
C     
      integer*4 i,j,np,nc
      real*4    pfabs,cost,sint,crs,ferm,xmax,dum,hit
      real*4    pfa,coa,dpf,dco,dat1,dat2,distri

      integer*4 getpid,kill
#ifndef gFortran
      external  kill
      external  getpid
#endif
      
      integer*4 ipid,iret
      
      ipid = getpid()
C--- use SIGNAL 3 : SIGQUIT to make core-dump file
      iret = kill(ipid,3)
      stop
      RETURN

      IERR=0
      IF (ICALL.EQ.0) THEN
C     
C     -- MAKE FERMI MOMENTUM -- DIRECTION MAP
C     
         DO 110 J=1,6
            XMAP(1,J)=0.
            XMAP(7,J)=0.
 110     CONTINUE

         DO 120 I=1,5
C K.Ishihara
C            PFABS=FLOAT(I)*0.05
C            PFABS=FLOAT(I)*(PFERMI/5.0)
            PFABS=FLOAT(I)*(PFMAX/5.0)
            DO 130 J=1,6
               COST=-1.+0.4*FLOAT(J-1)
               if ((abs(cost).ge.1.).or.(J.eq.1).or.(J.eq.6)) then
                  sint = 0.
               else
                  SINT=SQRT(1.-COST**2)
               endif
               PF(1)=PFABS*COST
               PF(2)=PFABS*SINT
               PF(3)=0.
               CRS=FNELSPF(E,PF,IPAR)
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
      
C 98/02/23 K.Ishihara
C 320  PFABS=RLU(DUM)*0.25
C 320  PFABS=RLU(DUM)*PFERMI
 320  PFABS=RLU(DUM)*PFMAX
      COST=-1.+RLU(DUM)*2.
      HIT=XMAX*RLU(DUM)
      PFA=PFABS/0.05
      COA=(COST+1.)/0.4
      NP=IFIX(PFA)+1
      NC=IFIX(COA)+1
      DPF=PFA+1.-FLOAT(NP)
      DCO=COA+1.-FLOAT(NC)
      DAT1=(1.-DPF)*XMAP(NP,NC)+DPF*XMAP(NP+1,NC)
      DAT2=(1.-DPF)*XMAP(NP,NC+1)+DPF*XMAP(NP+1,NC+1)
      DISTRI=(1.-DCO)*DAT1+DCO*DAT2
      IF(HIT.GT.DISTRI)GO TO 320

      IF (abs(cost).ge.1) then
         sint = 0.
      else
         SINT=SQRT(1.-COST**2)
      endif
      PFSEL(1)=PFABS*COST
      PFSEL(2)=PFABS*SINT
      PFSEL(3)=0.

      RETURN
      END
