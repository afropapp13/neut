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


      W2=W*W

      XMLEP2 = XMLEP*XMLEP
      
      ERESL=(W2+XMN2-Q2)/(2*XMN)
      APRESL=SQRT(ERESL**2-W2)
      ELEPL=(2*E*XMN+XMN2+Q2-W2)/(2*XMN)
      APLEPL=SQRT(ELEPL**2-XMLEP2)


      CSNULP=(Q2+2*E*ELEPL-XMLEP2)/(2*E*APLEPL)
      if ((1.0-csnulp**2).le.0.0) then
         snnulp = 0.0
         csnulp = 1.0
      else
         SNNULP=SQRT(1.0-CSNULP**2)
      endif


      CSNURS=(XMN2+XMLEP2-W2+2*E*ERESL-2*ELEPL*XMN)/(2*E*APRESL)
      if ((1.0-csnurs**2).le.0.0) then
         snnurs = 0.0
         csnurs = 1.0
      else
         SNNURS=-1.*SQRT(1-CSNURS**2)
      endif


      RTPNUL (1)=0.000
      RTPNUL (2)=0.000
      RTPNUL (3)=E

      PNUA (1)=0.000
      PNUA (2)=0.000
      PNUA (3)=E

      PNCL (1)=0.000
      PNCL (2)=0.000
      PNCL (3)=0.000

      PNCA (1)=0.000
      PNCA (2)=0.000
      PNCA (3)=0.000

      RTPLEPL(1)=SNNULP*APLEPL
      RTPLEPL(2)=0.
      RTPLEPL(3)=CSNULP*APLEPL

      PLEPA(1)=RTPLEPL(1) 
      PLEPA(2)=RTPLEPL(2)
      PLEPA(3)=RTPLEPL(3)

      PRESL(1)=SNNURS*APRESL
      PRESL(2)=0.
      PRESL(3)=CSNURS*APRESL

      PRESA(1)=PRESL(1) 
      PRESA(2)=PRESL(2)
      PRESA(3)=PRESL(3)

      DPRESL(1) =-1.* SNNURS
      DPRESL(2) = 0.
      DPRESL(3) =-1.* CSNURS

      GAMMA = 1./SQRT(1-(APRESL/SQRT(APRESL**2+W**2))**2)


      CALL MCVECBST(PLEPA,XMLEP,DPRESL,GAMMA)
      CALL MCVECBST(PNUA,0.000,DPRESL,GAMMA)
      CALL MCVECBST(PNCA,XMN,DPRESL,GAMMA)
      CALL MCVECBST(PRESA,W,DPRESL,GAMMA)


      AQA2=0.
      DO 10 I=1,3
         QA(I)=PNUA(I)-PLEPA(I)
         AQA2=AQA2+QA(I)**2
 10   CONTINUE
      AQA=SQRT(AQA2)
      DO 20 I=1,3
         ZAXIS(I)=QA(I)/AQA
 20   CONTINUE


      PA=SQRT
     $   ((W2**2+XMPI2**2+XMN2**2-2*W2*XMPI2-2*W2*XMN2-2*XMN2*XMPI2)
     $   /(4*W2))


      PPIAR(1)=PA*SIN(THETA)*COS(PHI)
      PPIAR(2)=PA*SIN(THETA)*SIN(PHI)
      PPIAR(3)=PA*COS(THETA)

      DO 50 I=1,3
         PPRAR(I)=-1.*PPIAR(I)
 50   CONTINUE

      COSZQA=ZAXIS(3)
      SINZQA=ZAXIS(1)


C      PpiA(3) = PpiAR(3)*coszQA-PpiAR(1)*sinzQA
C      PpiA(1) = PpiAR(3)*sinzQA+PpiAR(1)*coszQA

C      PprA(3) = PprAR(3)*coszQA-PprAR(1)*sinzQA
C      PprA(1) = PprAR(3)*sinzQA+PprAR(1)*coszQA

      PPIA(3) = PPIAR(3)*COSZQA-PPIAR(1)*SINZQA
      PPIA(1) = PPIAR(3)*SINZQA+PPIAR(1)*COSZQA
      PPIA(2) = PPIAR(2)

      PPRA(3) = PPRAR(3)*COSZQA-PPRAR(1)*SINZQA
      PPRA(1) = PPRAR(3)*SINZQA+PPRAR(1)*COSZQA
      PPRA(2) = PPRAR(2)

      
      DO 70 I=1,3
         RTPPIL(I)=PPIA(I)
         RTPPRL(I)=PPRA(I)
         DPRESL(I)=-1.*DPRESL(I)
 70   CONTINUE
         
      
      CALL MCVECBST(RTPPIL,XMPI,DPRESL,GAMMA)
      CALL MCVECBST(RTPPRL,XMN,DPRESL,GAMMA)


      RETURN
      END
      
