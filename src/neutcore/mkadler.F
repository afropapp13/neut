      SUBROUTINE MKAdler(cosThAd, phiAd)
C     Boost particles to get the Adler angles needed as input for
C     Minoo's cross-section calculations. Needed for the reweight
C     routine
C     Use the common block data to calculate the Adler angles
      implicit none

C     Get the filled common blocks
#include "nework.h"
#include "vcwork.h"
#include "posinnuc.h"
#include "rscons.h"

C     The returns
      real cosThAd, phiAd
C     Double versions for internal processing
      double precision cosThAdPr, phiAdPr

C     Nucleon, pion and resonant 4-momentum vectors
      double precision Nucl(4),Pion(4),Res(4)
C     Pion in resonance frame, rotated in Adler resonance frame 
      double precision PiOut(4),PiAdl(4)

C     Outgoing lepton and neutrino 4-momentum vectors
C     N.B. "Out" refers to the 4-momentum in the resonance frame
      double precision Lep(4),Neu(4)

C     The boosted lepton and neutrino
      double precision LepOut(4),NeuOut(4)
C     The 3-vectors of neutrino and lepton after boost
      double precision NeuVect(3),LepVect(3)
C     Looper
      integer I
C     The x axis, y axis and z axis in the Adler frame
      double precision xAxis(3),yAxis(3),zAxis(3)

C     4 particle is nucleon, 5th particle is pion ALWAYS
      DO 10 I=1,3
        Neu(I) = dble(PNE(I, 1))
        Lep(I) = dble(PNE(I, 3))
        Nucl(I) = dble(PNE(I, 4))
        Pion(I) = dble(PNE(I, 5))
        Res(I) = Nucl(I) + Pion(I)
  10  CONTINUE

C     Neutrino momentum
      Neu(4) = sqrt(Neu(1)**2+Neu(2)**2+Neu(3)**2)

C     Set mass of nucleon depending on PDG
      IF (IPNE(4) .EQ. 2112) THEN
        Nucl(4) = dble(XMNE)
      ELSE
        Nucl(4) = dble(XMP)
      ENDIF

      Nucl(4)=XMN

C     Set fourth component to be energy
      Nucl(4) = sqrt(Nucl(1)**2+Nucl(2)**2+
     &               Nucl(3)**2+Nucl(4)**2)

C     Mass of pion (charged)
C     IF (ABS(IPNE(5)).EQ.211) THEN
C       Pion(4) = dble(0.13957)
C     Mass of pion (neutral)
C     ELSE
C       Pion(4) = dble(0.13497)
C     ENDIF
      Pion(4) = XMPI
      Pion(4) = sqrt(Pion(1)**2+Pion(2)**2+
     &               Pion(3)**2+Pion(4)**2)

C     Set lepton mass for lepton energy by looking at PID
      IF (abs(IPNE(3)).EQ.11) THEN
         Lep(4) = dble(XME)
      ELSE IF (abs(IPNE(3)).EQ.13) THEN
         Lep(4) = dble(XMMU)
      ELSE IF (abs(IPNE(3)).EQ.15) THEN
         Lep(4) = dble(XMTAU)
      ELSE 
         Lep(4) = dble(0.)
      endIF

C     Outgoing lepton energy
      Lep(4) = sqrt(Lep(1)**2+Lep(2)**2+
     &              Lep(3)**2+Lep(4)**2)

C     Energy of resonance is sum of pion and nucleon energy which it
C     decayed into
      Res(4) = Pion(4)+Nucl(4)

C     Do Lorentz boost using CERNLIB
C     First argument is 4-vector (P,E), second argument is 4-vector
C     (p,e), third argument is return 4-vector (p',e') --- i.e. boost
C     second argument to first argument frame

C     Boost the neutrino
      call LOREN4_DBL(RES, Neu, NeuOut)

C     Boost the muon
      call LOREN4_DBL(RES, Lep, LepOut)

C     Boost the pion
      call LOREN4_DBL(RES, Pion, PiOut)

C     Now we have the pion, outgoing lepton and incoming neutrino boosted in to the resonance rest frame
C     Make the zaxis defined as the momentum transfer in the resonance
C     rest frame, z = k = k1 - k2 (k1 = neutrino, k2 = muon)
      DO 30 I=1,3
        zAxis(I) = (NeuOut(I) - LepOut(I))
  30  CONTINUE

C     Normalise the z-axis vector to unit
      call normalise(zAxis)

C     Need to reduce the 4-vector to a three vector to use CERNLIB cross
C     product
      DO 50 I=1,3
        NeuVect(I) = NeuOut(I)
        LepVect(I) = LepOut(I)
  50  CONTINUE

C     Have zaxis defined above, now define y axis as k1 x k2
      call CROSS_DBL(NeuVect, LepVect, yAxis)
      call normalise(yAxis)

C     And finally the x-axis as (k1 x k2) x k
      call CROSS_DBL(yAxis, zAxis, xAxis)
      call normalise(xAxis)

C     Now we project the pion onto the new axes and get the adler angles
      DO 100 I=1,3
        PiAdl(I) = 0.
 100  CONTINUE
C     Do the dot product of the boosted pion with the axes
      DO 101 I=1,3
C       The x projection
        PiAdl(1)=PiAdl(1)+PiOut(I)*xAxis(I)
C       The y projection
        PiAdl(2)=PiAdl(2)+PiOut(I)*yAxis(I)
C       The z projection
        PiAdl(3)=PiAdl(3)+PiOut(I)*zAxis(I)
 101  CONTINUE

C     write(*,*) "Calc Pion in Adler: ", PiAdl
C     write(*,*) "Energy:             ", PiOut(4)

C     Now we have the pion and nucleon in the Adler frame defined as
C     PiAdl and NuclAdl

C     Cos theta Adler is THEN defined as the angle between the pion in
C     the resonant rest frame and the z, costheta = z/r
      cosThAdPr=PiAdl(3)/sqrt(PiAdl(1)**2+PiAdl(2)**2+PiAdl(3)**2)

C     Tan phi defined as y/x
C     Need to use atan2 here because inverse tan
      phiAdPr=atan2(PiAdl(2),PiAdl(1))

C     Check if negative -> force positive
C     We define phi from 0 to 2pi, atan2 returns -pi to +pi
      IF (phiAdPr .LT. 0) THEN
        phiAdPr = 2.*PI + phiAdPr
      END IF

C     Finally convert to single precision
      phiAd=phiAdPr
      cosThAd=cosThAdPr

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE NORMALISE(INPUT)
C     helper function to normalise a three vector, e.g. momentum
      implicit none
      double precision input(3), mag
      integer i

C     Get the magnitude of the vector
      mag = sqrt(input(1)**2+input(2)**2+input(3)**2)
C     Normalise the vector by magnitude
      DO 10 I=1,3
        input(I) = input(I)/mag
  10  CONTINUE

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE NORMALISE4D(INPUT)
C     helper function to normalise a three vector, e.g. momentum
      implicit none
      double precision input(4), mag
      integer i

C     Get the magnitude of the vector
      mag = sqrt(input(1)**2+input(2)**2+input(3)**2
     &          +input(4)**2)
C     Normalise the vector by magnitude
      DO 10 I=1,4
        input(I) = input(I)/mag
  10  CONTINUE

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CROSS_DBL(X,Y,Z)
C     helper function to cross-product in double precision
C     Just a copy of CERNLIB CROSS function but better precision
      implicit none
      double precision X(3), Y(3), Z(3)
      double precision Z1, Z2

      Z1=X(2)*Y(3)-X(3)*Y(2)
      Z2=X(3)*Y(1)-X(1)*Y(3)
      Z(3)=X(1)*Y(2)-X(2)*Y(1)
      Z(1)=Z1
      Z(2)=Z2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE LOREN4_DBL(DIR,P4IN,P4OUT)
C     helper function to Lorentz boost in double precision
C     Just a copy of CERNLIB LOREN4function but better precision
      IMPLICIT NONE
      DOUBLE PRECISION PCM2, ONMCM, EPBETA, PROD
      DOUBLE PRECISION DIR(4),P4IN(4),P4OUT(4)
      integer i

      PCM2=DIR(1)*DIR(1)+DIR(2)*DIR(2)+DIR(3)*DIR(3)
      ONMCM=1.D0/ SQRT (DIR(4)*DIR(4)-PCM2)

      EPBETA=P4IN(1)*DIR(1)+P4IN(2)*DIR(2)+P4IN(3)*DIR(3)
      PROD=EPBETA*(DIR(4)*ONMCM-1.D0)/PCM2-P4IN(4)*ONMCM
      P4OUT(4)=ONMCM*(P4IN(4)*DIR(4)-EPBETA)
      DO 50 I=1,3
        P4OUT(I)=P4IN(I)+DIR(I)*PROD
  50  CONTINUE

      RETURN
      END
