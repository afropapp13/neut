**********************************************************************
*     ------------------
      SUBROUTINE TAUVECT
*     ------------------
*     
*     ( purpose )
*       execute tau decay and make vector
*
*     ( input ) 
*       COMMON VCWORK, VCVRTX
*
*     ( output ) 
*       COMMON VCWORK, VCVRTX
*
*     ( creation date and author )
*       1995.10.26 ; K.Kaneyuki
*
**********************************************************************
#include "vcwork.h"
#include "vcvrtx.h"
#include "taucom.h"
#include "tahepevt.h"

      REAL DIR0(3),DIR1(3),DIR2(3),P(3)
      IF(NVC.LE.0)RETURN

      NALL = NVC
      DO 10 I=1,NVC
C
C --     TAU ?
C
         IF (ABS(IPVC(I)).EQ.15 .AND. IFLGVC(I).EQ.0) THEN
C
C --        EXECUTE TAU DECAY ( TAUOLA Vr.2.6 )
C
C           ( DIRECTION OF TAU = Z FRAME )
C            
            CALL TAUDECAY(IPVC(I),PVC(1,I))

C
C --        MAKE TAU BANK
C
            CALL NEMKTAUBK

C
C --        FILL COMMON OF TAU
C            
            IF (NHEP.LE.1) RETURN
            IFLGVC(I) = 1
            ICRNVC(I) = 0
            CALL UCOPY(POSIVC(1,I),POSFVC(1,I),3)
            IVTFVC(I) = IVTIVC(I)
C
C --        DETERMINE TAU DIRECTION
C            
            CALL UCOPY(PVC(1,I),DIR0,3)
            DIRABS=SQRT(DIR0(1)**2+DIR0(2)**2+DIR0(3)**2)
            IF (DIRABS.LT.1.0E-5) THEN
               WRITE(6,*) ' TAUVECT : ABS(DIR)=0 '
               CALL RNDIR(DIR0)
            ELSE
               DIR0(1)=DIR0(1)/DIRABS
               DIR0(2)=DIR0(2)/DIRABS
               DIR0(3)=DIR0(3)/DIRABS
            ENDIF
            COST = DIR0(3)
            SINT2= 1.0-COST**2
            IF (SINT2.LT.0.0) SINT2=0.0
            SINT = SQRT(SINT2)
            XYABS=SQRT(DIR0(1)**2+DIR0(2)**2)
            IF (XYABS.LT.1.0E-5) THEN
               COSP=1.0
               SINP=0.0
            ELSE
               COSP=DIR0(1)/XYABS
               SINP=DIR0(2)/XYABS
            ENDIF
C
C --        FILL EACH PARTICLE
C
            NALL2=NALL

C     Used for debugging
C            write (*,*) "NHEP: ", nhep

            DO 20 J=2,NHEP
               if (NALL.ge.99) goto 20
               NALL=NALL+1
               IPVC(NALL)=IDHEP(J)
c               CALL MCMASS(IPVC(J),AMASVC(NALL))
               CALL MCMASS(IPVC(NALL),AMASVC(NALL))
C
C --           ROTATE EACH PARTICLE
C               
               CALL VZERO(DIR1,3)
               CALL VZERO(DIR2,3)
               P(1)=PHEP(1,J)*1.0E+3
               P(2)=PHEP(2,J)*1.0E+3
               P(3)=PHEP(3,J)*1.0E+3
               PABS=SQRT(P(1)**2+P(2)**2+P(3)**2)
               IF (PABS.GT.1.0E-5) THEN
                  DIR1(1)=P(1)/PABS
                  DIR1(2)=P(2)/PABS
                  DIR1(3)=P(3)/PABS
                  DIR2(1)=DIR1(1)*COST+DIR1(3)*SINT
                  DIR2(2)=DIR1(2)
                  DIR2(3)=-DIR1(1)*SINT+DIR1(3)*COST
                  CALL UCOPY(DIR2,DIR1,3)
                  DIR2(1)=DIR1(1)*COSP+DIR1(2)*SINP
                  DIR2(2)=DIR1(1)*SINP-DIR1(2)*COSP
                  DIR2(3)=DIR1(3)
               ENDIF
               PVC(1,NALL)=PABS*DIR2(1)
               PVC(2,NALL)=PABS*DIR2(2)
               PVC(3,NALL)=PABS*DIR2(3)
C
C --           FILL PARAMETERS
C               

               IF (JMOHEP(1,J).EQ.1) IORGVC(NALL)=I
c -- tau decay products: 
c               IF (JMOHEP(1,J).NE.1) IORGVC(NALL)=NALL2+JMOHEP(1,J)-1
               IF (JMOHEP(1,J).NE.1) IORGVC(NALL) = 15

               IFLGVC(NALL)=0
               IF (JDAHEP(1,J).NE.0) IFLGVC(NALL)=1
               IF (ABS(IPVC(NALL)).EQ.12) IFLGVC(NALL)=2
               IF (ABS(IPVC(NALL)).EQ.14) IFLGVC(NALL)=2
               IF (ABS(IPVC(NALL)).EQ.16) IFLGVC(NALL)=2
c -- Fix JDAHEP bug
               IF (ABS(IPVC(NALL)).EQ.11)  IFLGVC(NALL)=0  ! e
               IF (ABS(IPVC(NALL)).EQ.13)  IFLGVC(NALL)=0  ! mu
               IF (ABS(IPVC(NALL)).EQ.211) IFLGVC(NALL)=0  ! pi+/-
               IF (ABS(IPVC(NALL)).EQ.111) IFLGVC(NALL)=0  ! pi0
               IF (ABS(IPVC(NALL)).EQ.321) IFLGVC(NALL)=0  ! K+/-

               ICRNVC(NALL)=0
               IF (IFLGVC(NALL).EQ.0) ICRNVC(NALL)=1
               CALL UCOPY(POSFVC(1,I),POSIVC(1,NALL),3)
               IVTIVC(NALL)=IVTFVC(I)
               IF (IFLGVC(NALL).EQ.1) THEN
                  CALL UCOPY(POSFVC(1,I),POSFVC(1,NALL),3)
                  IVTFVC(NALL)=IVTFVC(I)
               ENDIF

C     Use this to write out the decay particles for testing
C
C
C               write (*,*) "IPVC: ", nall, " Particle: ", ipvc(nall),
C     &         " Flag: ", iflgvc(nall), " Track: ", icrnvc(nall),
C     &         " Iorg: ", iorgvc(nall) 


 20         CONTINUE
         ENDIF
 10   CONTINUE
      NVC=NALL

      RETURN
      END
