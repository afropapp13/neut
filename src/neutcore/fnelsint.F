************************************************************************
*     -------------------------
      FUNCTION FNELSINT(E,IPAR)
*     -------------------------
*
*     (Purpose)
*       GIVE TOTAL CROSS SECTION OF ELASTIC SCATTERING
*       INTEGRATED BY FERMI MOMENTUM
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
*     
*     (Output)
*       FNELSINT : TOTAL CROSS SECTION ( 10^-30cm^2 )
*
*     (Creation Date and Author)
*       1984.01.08 ; M.NAKAHATA
*       1987.06.20 ; N.SATO  FOR TAU
*       1995.02.10 ; K. KANEYUKI  FOR S.K.
*                    ELSINT -> FNELSINT
*     
************************************************************************
      IMPLICIT NONE

#include "neutparams.h"

      real*4    FNELSINT
      REAL*4    E
      INTEGER*4 IPAR

      real*4 fnelspf,fvfrmgev
      external fnelspf,fvfrmgev

      REAL*4 PF(3)

      INTEGER*4 IMAX,JMAX
      PARAMETER (IMAX=10)
      PARAMETER (JMAX=11)
      REAL*4 XMAP(IMAX,JMAX)

      real*4 DPF,dcost,xnorm,cost,sint,crs,ferm,pfabs

      integer*4 i,j
      
C      DPF=PFERMI/FLOAT(IMAX)
      DPF=PFMAX/FLOAT(IMAX)
      DCOST=2./FLOAT(JMAX-1)
      XNORM=0.
      FNELSINT=0.
C     
C     ++ INTEGRATE BY FERMI MOMENTUM
C     
C      write(*,'(A,F10.8,A)') "Energy = ",E,'GeV'

      DO 100 I=1,IMAX
         PFABS=FLOAT(I)*DPF
C         write(*,'(I3$)') int(PFABS*1000.)
         DO 200 J=1,JMAX
            COST=-1.+DCOST*FLOAT(J-1)
            IF ((J.eq.1).or.(j.eq.JMAX)) THEN
               SINT=0.
            ELSE
               SINT=SQRT(1.-COST**2)
            ENDIF
            PF(1)=PFABS*COST
            PF(2)=PFABS*SINT
            PF(3)=0.
            CRS=FNELSPF(E,PF,IPAR)
C**** CRS=ELSPF3(E,PF,IPAR)
C**** CRS=ELSPF4(E,PF,IPAR)
            FERM=FVFRMGEV(PFABS)
C            write(6,*) i,j,crs,ferm,pf,e
C            write(*,'(F7.4$)') CRS*10.
            XMAP(I,J)=CRS*FERM
            FNELSINT=FNELSINT+XMAP(I,J)
            XNORM=XNORM+FERM
  200    CONTINUE
C         write(*,*) ' '
  100 CONTINUE
      FNELSINT=FNELSINT/XNORM
C      write(6,*) fnelsint,xnorm
C     
      RETURN
      END
