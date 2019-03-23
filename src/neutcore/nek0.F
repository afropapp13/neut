************************************************************************
*     ------------------------------------------------
      SUBROUTINE NEK0
*     ------------------------------------------------
*
*     (Purpose)
*      Translate K0 ,K0BAR  -> K LONG (50%) , K SHORT (50%)
*
*
*     (Input)
*       COMMON /VCWORK/
*       COMMON /VCVRTX/
*
*     (Output)
*       COMMON /VCWORK/
*       COMMON /VCVRTX/
*
*       K0 : 311 ,K0BAR : -311  KL : 130  KS : 310   
*
*     (Creation Date and Author)
*       1998.03.05    J.Kameda 
*       2006.06.15    Y.Hayato
*                         check icrnvc
*                         copy  interaction positions
*       2011.05.25    Y.Hayato
*                         RNKIND was defined as integer -> KS only
*
*************************************************************************

      implicit none

#include "vcvrtx.h"
#include "vcwork.h"
      
#include "posinnuc.h"

C      INTEGER  RNKIND,I,NVC2,J
      INTEGER  I,NVC2,J
      REAL     DUMMY,RNKIND
  
      integer*4 KTYP

      real*4    RLU
      external  RLU


      NVC2 = NVC

      DO 10 I=1,NVC

         IF ((ICRNVC(I).ne.0).and.(NVC2.lt.MAXVC)) then
            IF (ABS(IPVC(I)).eq.311) then
               KTYP = 310
               RNKIND = RLU(DUMMY)
               IF  (RNKIND.ge.0.5) KTYP = 130
               
               IFLGVC(I)=4
               ICRNVC(I)=0
               CALL UCOPY(POSVC,POSFVC(1,I),3)
               IVTFVC(I)=1
               
               NVC2=NVC2+1
               IPVC(NVC2)=KTYP
               CALL UCOPY(PVC(1,I),PVC(1,NVC2),3)
               IORGVC(NVC2)=I
               IFLGVC(NVC2)=0
               ICRNVC(NVC2)=1
               CALL UCOPY(POSVC,POSIVC(1,NVC2),3)
               IVTIVC(NVC2)=1
               CALL MCMASS(IPVC(NVC2),AMASVC(NVC2))
               
C copy generated position
               do 100 J=1,3
                  POSNUC(J,NVC2) = POSNUC(J,I)
 100           continue

            endif
         ENDIF

 10   continue

      NVC = NVC2
      RETURN
      END
