************************************************************************
*     ------------------------------------------------
      SUBROUTINE NEOMEGA
*     ------------------------------------------------
*
*     (Purpose)
*      DECAY OMEGA MESON  -> PI+ PI- PI0 (88.8%) 
*                            PI0 GAMMA   (8.5%)
*                            PI+ PI-     (2.2%)
*
*     (Input)
*       COMMON /VCWORK/
*       COMMON /VCVRTX/
*
*     (Output)
*       COMMON /VCWORK/
*       COMMON /VCVRTX/
*
*     (Creation Date and Author)
*       2002.08.02    M.Ishitsuka
*
*************************************************************************

#include "vcvrtx.h"
#include "vcwork.h"
#include "posinnuc.h"
      
      INTEGER  I,K,LL,NVC2
      INTEGER MULDCY
      REAL PBUF(3,10)
      INTEGER IPBUF(10)
  
      NVC2 = NVC

      DO I=1,NVC
         IF (IPVC(I).eq.223 .AND. ICRNVC(I).EQ.1) then
            CALL EFOMEGADCY(PVC(1,I),MULDCY,IPBUF,PBUF)

            IFLGVC(I)=1
            ICRNVC(I)=0
            CALL UCOPY(POSVC,POSFVC(1,I),3)
            IVTFVC(I)=1

            DO K=1,MULDCY
               NVC2=NVC2+1
               IPVC(NVC2)=IPBUF(K)
               CALL UCOPY(PBUF(1,K),PVC(1,NVC2),3)
               IORGVC(NVC2)=I
               IFLGVC(NVC2)=0
               ICRNVC(NVC2)=1
               CALL UCOPY(POSVC,POSIVC(1,NVC2),3)
               IVTIVC(NVC2)=1
               CALL MCMASS(IPVC(NVC2),AMASVC(NVC2))
               DO LL=1,3
                  POSNUC(LL,NVC2)=9999.
               ENDDO
            ENDDO
         ENDIF
      ENDDO

      NVC = NVC2

      RETURN
      END
