#ifndef EFOMEGA_INCLUDED
      INTEGER*4 NOMEGATL,IPINTOMEGA,IPIORIOMEGA,LOOPOMEGA
      REAL*4    PINTOMEGA,XINTOMEGA
      COMMON /OMEGAOUT/ NOMEGATL,IPINTOMEGA(50),PINTOMEGA(3,50),
     &     XINTOMEGA(3,50),IPIORIOMEGA(50),LOOPOMEGA
#define EFOMEGA_INCLUDED
#endif
