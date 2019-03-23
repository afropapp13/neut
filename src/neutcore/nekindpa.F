**********************************************************************
*     ------------------------
      SUBROUTINE NEKINDPA(MUL)
*     ------------------------
*     
*     ( purpose )
*       SET KIND OF PARTICLES OF HADRON SIDE
*
*     ( input ) 
*       MUL     : MULTIPLICITY
*
*     ( output ) 
*       COMMON NEWORK
*
*     ( creation date and author )
*       198?.??.?? ; First version by M.Nakahata
*       1995.03.03 ; K.Kaneyuki for S.K.
*                    KINDPA -> NEKINDPA
*       2005.03.30 ; M.Ishitsuka
*                    Change the ratio of pi+/-/0 to be uniform
*
**********************************************************************
#include "nework.h"
C     
      ICTOT=NECHARGE(IPNE(2))-NECHARGE(IPNE(3))
C     
   10 IPNE(4)=2212
      RAN=RLU(DUM)
      IF(RAN.GT.0.5)IPNE(4)=2112
      IORGNE(4)=2
      ICRNNE(4)=1
      ICHAD=NECHARGE(IPNE(4))
      DO 100 I=1,MUL-1
         RAN=RLU(DUM)
c         IRAN=IFIX(RAN*1000.)
c         ITYPE=MOD(IRAN,3)
c         IF(ITYPE.EQ.0)IPP=111
c         IF(ITYPE.EQ.1)IPP=211
c         IF(ITYPE.EQ.2)IPP=-211
         IF(RAN.LT.0.333*1.05) THEN
            IPP=211
         ELSE IF(RAN.LT.(0.333*1.05+0.333*1.05)) THEN
            IPP=-211
         ELSE
            IPP=111
         ENDIF
         IPNE(I+4)=IPP
         IORGNE(I+4)=2
         ICRNNE(I+4)=1
         ICHAD=ICHAD+NECHARGE(IPP)
  100 CONTINUE
      IF(ICHAD.NE.ICTOT)GO TO 10
      RETURN
      END
