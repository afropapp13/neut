***********************************************************************
*     ---------------------------------------------------
      SUBROUTINE RSLZBT(E,Q2,W,XMLEP,THETA,PHI,RTPLEPL,RTPPIL,RTPPRL)
*     ---------------------------------------------------
*
*     ( purpose )
*       calculate kinematics(lorentz boost to LAB. frame)
*
*     ( input )
*       E           : Energy of incident neutrino
*       q2          : momentum transfer
*       W           : invariant mass
*       xmlep       : mass of outgoing lepton
*       theta       : angle of pion in ADLER frame
*       phi         : angle of pion in ADLER frame
*     
*     ( output )
*       rtplepl(3)  : outgoing lepton  momentum
*       rtppil(3)   : outgoing pion    momentum
*       rtpprl(3)   : outgoing nucleon momentum
*
*     ( creation date and author )
*       1993.Dec. ; for Super-Kamioka by Y.Hayato
*                   (Based on Rein and Sehgal's code)
*
*     ( comment )
*       Reference:Rein and Sehgal,Ann. of Phys. 133(1981),79-153
*                 Rein           ,Z.Phys.C 35(1987),43-64
*                 Feynman et al. ,Phys.Rev.D3(1971),2706-
***********************************************************************

      IMPLICIT NONE

      REAL E,Q2,W,THETA,PHI
      REAL XMLEP
      REAL RTPLEPL(3),RTPPIL(3),RTPPRL(3)

#include"rscons.h"

      INTEGER I

      REAL W2,XMLEP2
      REAL ERESL,APRESL,ELEPL,APLEPL
      REAL CSNULP ,SNNULP ,CSNURS,SNNURS
      REAL COSZQA,SINZQA

      REAL GAMMA,PA
      REAL ZAXIS(3)
      REAL RTPNUL(3) ,PNCL(3) , PRESL(3) , QA(3)
      REAL PNUA(3)   ,PNCA(3) , PRESA(3) , PLEPA(3)      , PPIA(3) 
      REAL PPRA(3) , PPRAR(3),PPIAR(3)
      REAL DPRESL(3)                
      REAL AQA
      REAL AQA2

#include "necard.h"

      W2=W*W

      XMLEP2 = XMLEP*XMLEP
      
C     Calculate the various event kinematics in the resonance rest frame
C     All of these can be derived by using four-vector conservation

C     Energy of resonance
      ERESL=(W2+XMN2-Q2)/(2*XMN)
C     Momentum of resonance
      APRESL=SQRT(ERESL**2-W2)
C     Lepton energy
      ELEPL=(2*E*XMN+XMN2+Q2-W2)/(2*XMN)
C     Momentum of lepton
      APLEPL=SQRT(ELEPL**2-XMLEP2)

C     Cosine between neutrino and lepton (get from Q2)
      CSNULP=(Q2+2*E*ELEPL-XMLEP2)/(2*E*APLEPL)
C     Check that the cosine is between -1 and 1
      if ((1.0-csnulp**2) .LE. 0) then
         snnulp = 0.0
         csnulp = 1.0
C     Sine between neutrino and lepton
      else
         SNNULP=SQRT(1.0-CSNULP**2)
      endif

C     Cosine between neutrino and resonance
      CSNURS=(XMN2+XMLEP2-W2+2*E*ERESL-2*ELEPL*XMN)/(2*E*APRESL)
      if ((1.0-csnurs**2) .LE. 0.) then
         snnurs = 0.0
         csnurs = 1.0
C     Sine between neutrino and resonance
      else
         SNNURS=-1.*SQRT(1-CSNURS**2)
      endif

C     Define z-axis as being the neutrino incoming direction
      RTPNUL (1)=0.000
      RTPNUL (2)=0.000
      RTPNUL (3)=E

C     Define z-axis as being the neutrino incoming direction
      PNUA (1)=0.000
      PNUA (2)=0.000
      PNUA (3)=E

      PNCL (1)=0.000
      PNCL (2)=0.000
      PNCL (3)=0.000

      PNCA (1)=0.000
      PNCA (2)=0.000
      PNCA (3)=0.000

C     Get the x-axis projection of the lepton
      RTPLEPL(1)=SNNULP*APLEPL
C     Define the y-axis component to be zero?
      RTPLEPL(2)=0.
C     Get the z-axis projection of the lepton
      RTPLEPL(3)=CSNULP*APLEPL

C     Set the lepton momentum to the above derived quantity
C     N.B. don't know why it's repeated
      PLEPA(1)=RTPLEPL(1) 
      PLEPA(2)=RTPLEPL(2)
      PLEPA(3)=RTPLEPL(3)

C     Do the x-axis projection of the resonance momentum
      PRESL(1)=SNNURS*APRESL
      PRESL(2)=0.
C     Do the z-axis projection of the resonance momentum
      PRESL(3)=CSNURS*APRESL

C     Set the resonance momentum to the derived quantities
      PRESA(1)=PRESL(1) 
      PRESA(2)=PRESL(2)
      PRESA(3)=PRESL(3)

C     Make the unit vector of the resonance 3 vector
      DPRESL(1) =-1.* SNNURS
      DPRESL(2) = 0.
      DPRESL(3) =-1.* CSNURS

C     Get the relativistic gamma factor for the resonance rest frame
C     GAMMA = 1./SQRT(1-(APRESL/SQRT(APRESL**2+W**2))**2)
C     Avoid some maths operators by rearranging the above
      GAMMA = SQRT(1. + (APRESL/W)**2)

C     Boost outgoing from resonance rest frame to lab frame
C     PLEPA, PNUA, PRESA are the returned 3 vectors
      CALL MCVECBST(PLEPA, XMLEP, DPRESL, GAMMA)
      CALL MCVECBST(PNUA, 0.000, DPRESL, GAMMA)
      CALL MCVECBST(PNCA, XMN, DPRESL, GAMMA)
      CALL MCVECBST(PRESA, W, DPRESL, GAMMA)

C     This is the |q3|^2
      AQA2=0.
      DO 10 I=1,3
C     Calculate the q3 vector
         QA(I)=PNUA(I)-PLEPA(I)
C     Calculate the |q3|^2
         AQA2=AQA2+QA(I)**2
 10   CONTINUE

C     The |q3| quantity (q3 is the vector part of q)
      AQA=SQRT(AQA2)
C     The z-axis should be along the q3, make a unit vector of q3
      DO 20 I=1,3
         ZAXIS(I)=QA(I)/AQA
 20   CONTINUE

C     The magnitude of the pion momentum in the resonance rest frame
      PA=SQRT
     $   ((W2**2+XMPI2**2+XMN2**2-2*W2*XMPI2-2*W2*XMN2-2*XMN2*XMPI2)
     $   /(4*W2))

C     Project pion momentum using spherical polar coordinates
C     X component
      PPIAR(1)=PA*SIN(THETA)*COS(PHI)
C     Y component
      PPIAR(2)=PA*SIN(THETA)*SIN(PHI)
C     Z component
      PPIAR(3)=PA*COS(THETA)

C     Make the nucleon in the negative direction of pion
      DO 50 I=1,3
C     CWRET save the pion
         ppirs(I)=ppiar(I)
         PPRAR(I)=-1.*PPIAR(I)
 50   CONTINUE
      ppirs(4) = sqrt(XMPI**2+ppirs(1)**2+ppirs(2)**2+ppirs(3)**2)

      COSZQA=ZAXIS(3)
      SINZQA=ZAXIS(1)

C     Pion momentum projected onto the boosted non-Adler axis
      PPIA(3) = PPIAR(3)*COSZQA-PPIAR(1)*SINZQA
      PPIA(1) = PPIAR(3)*SINZQA+PPIAR(1)*COSZQA
      PPIA(2) = PPIAR(2)

C     Nucleon momentum projected onto the boosted non-Adler axis
      PPRA(3) = PPRAR(3)*COSZQA-PPRAR(1)*SINZQA
      PPRA(1) = PPRAR(3)*SINZQA+PPRAR(1)*COSZQA
      PPRA(2) = PPRAR(2)

      DO 70 I=1,3
         RTPPIL(I)=PPIA(I)
         RTPPRL(I)=PPRA(I)
         DPRESL(I)=-1.*DPRESL(I)
 70   CONTINUE
         
C     Boost the nucleon and pion back into rest frame 
      CALL MCVECBST(RTPPIL,XMPI,DPRESL,GAMMA)
      CALL MCVECBST(RTPPRL,XMN,DPRESL,GAMMA)

      RETURN
      END
      
