#ifndef EFPIINTH_INCLUDED
      integer*4 NPARTL,IPINT(50),IPIORI(50),LOOPPI
      real*4    PINT(3,50),XINT(3,50),POSPIPROD(3,50)
      COMMON /PIINTH/ NPARTL,IPINT,PINT,XINT,IPIORI,
     &                LOOPPI
#define EFPIINTH_INCLUDED
#endif
