***********************************************************************
*     ---------------------------------------------------
      SUBROUTINE RSDIST(IMODE,INLEP,ENU,XINPPR,RETPMU,RETPPI,RETPPR,
     $                  IERR)
*     ---------------------------------------------------
*
*     ( purpose )
*       calculate kinematics of final state particles
*
*     ( input )
*       imode       : interaction mode
*       inlep       : kind of incoming neutrino
*       enu         : energy of incoming neutrino
*       inppr(3)    : momentum of incoming nucleon
*     
*     ( output )
*       retpmu(3)    : momentum of outgoing lepton
*       retppi(3)    : momentum of outgoing pion
*       retppr(3)    : momentum of outgoing nucleon
*
*     ( creation date and author )
*       1993.Dec.  ; for Super-Kamioka by Y.Hayato
*                   (Based on Rein and Sehgal's code)
*       1995.Jul.14; Bug FIX(treatment of IFREE)
*
*     ( comment )
*       Reference:Rein and Sehgal,Ann. of Phys. 133(1981),79-153
*                 Rein           ,Z.Phys.C 35(1987),43-64
*                 Feynman et al. ,Phys.Rev.D3(1971),2706-
***********************************************************************


      IMPLICIT NONE

      INTEGER*4 IMODE,INLEP
      REAL*4    ENU
      REAL*4    XINPPR(3),RETPMU(3),RETPPI(3),RETPPR(3)
      INTEGER*4 IERR

#include "neutmodel.h"
#include "rscons.h"

      REAL*4    XPNU
      REAL*4    DINPPR(3),XINPNU(3)
      REAL*4    TMPPMU(3),TMPPPI(3),TMPPPR(3)
      REAL*4    XMLEP
      REAL*4    GAMMA,EPR,PPR,COSTH,SINTH,COSPHI,SINPHI

      INTEGER*4 I

C---initialize valiables
      IERR = 0

      XMLEP = 0.00000
      if ((IMODE.lt.4).or.((IMODE.gt.10).and.(IMODE.lt.14)))
     $     then
         IF (ABS(INLEP).EQ.12) XMLEP=XME
         IF (ABS(INLEP).EQ.14) XMLEP=XMMU
         IF (ABS(INLEP).EQ.16) XMLEP=XMTAU
      endif
      
      GAMMA = 0.00000

      XINPNU(1) = 0.00000
      XINPNU(2) = 0.00000
      XINPNU(3) = ENU

      EPR = SQRT(XINPPR(1)**2 + XINPPR(2)**2 + XINPPR(3)**2 + XMN2)
      PPR = SQRT(XINPPR(1)**2 + XINPPR(2)**2 + XINPPR(3)**2)
      IF ((PPR.lt.0.00001)) THEN
        DINPPR(1)=1.0
        DINPPR(2)=0.0
        DINPPR(3)=0.0
        XINPPR(1)=0.0
        XINPPR(2)=0.0
        XINPPR(3)=0.0
      ELSE 
        DO 10 I=1,3
          XINPPR(I) = XINPPR(I)
          DINPPR(I) = -1.*XINPPR(I)/PPR
          GAMMA = GAMMA + (XINPPR(I)/EPR)**2
 10     CONTINUE
      ENDIF
      GAMMA = 1./SQRT(1.- GAMMA)

C     Boost to nucleon frame
      CALL MCVECBST(XINPNU,0.000,DINPPR,GAMMA)
      XPNU  = SQRT(XINPNU(1)**2 + XINPNU(2)**2 + XINPNU(3)**2)
C

      COSTH= XINPNU(3)/XPNU
      if (abs(costh).eq.1) then
         sinth = 0.
      else
         SINTH= SQRT(1.-COSTH**2)
      endif
      IF (SINTH.gt.0.00001) THEN
        COSPHI= XINPNU(1)/XPNU/SINTH
        SINPHI= XINPNU(2)/XPNU/SINTH
      ELSE
        COSPHI=1.00000
        SINPHI=0.00000
      ENDIF

C     Get the outgoing kinematics from Rein-Sehgal implementation
      IF (MDLSPI .EQ. 1) THEN
        CALL RSDT33(IMODE,INLEP,XPNU,TMPPMU,TMPPPI,TMPPPR,IERR)
        if (ierr.ne.0) then
           return
        endif
C     Alternatively call Minoo's model (MDLSPI)
      ELSE IF (MDLSPI .EQ. 3) THEN
        CALL MKDT33(IMODE,INLEP,XPNU,TMPPMU,TMPPPI,TMPPPR)
      ELSE
        WRITE(*,*) 'WRONG MDLSPI GIVEN'
        WRITE(*,*) 'YOU GAVE MDLSPI = ', MDLSPI
        WRITE(*,*) '1=REIN-SEHGAL : Global Fermi Gas'
        WRITE(*,*) '2=REIN-SEHGAL : Local Fermi Gas (not implemented)'
        WRITE(*,*) '3=MINOO KABIRNEZHAD'
        STOP
      ENDIF

C     Rotate around
      CALL RS3DRT(TMPPMU,SINTH,COSTH,SINPHI,COSPHI)
      CALL RS3DRT(TMPPPI,SINTH,COSTH,SINPHI,COSPHI)
      CALL RS3DRT(TMPPPR,SINTH,COSTH,SINPHI,COSPHI)

      DO 20 I=1,3
         DINPPR(I) = -1.*DINPPR(I)
         RETPMU(I) = TMPPMU(I)
         RETPPI(I) = TMPPPI(I)
         RETPPR(I) = TMPPPR(I)
 20   CONTINUE

C     Boost out of nucleon frame
      CALL MCVECBST(RETPMU,XMLEP,DINPPR,GAMMA)
      CALL MCVECBST(RETPPI,XMPI,DINPPR,GAMMA)
      CALL MCVECBST(RETPPR,XMN,DINPPR,GAMMA)

      RETURN
      END
      
