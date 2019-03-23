*     -----------------------------------------------------
      SUBROUTINE NERANBLL(AA,PIN,OPOS,XPOS)
*     -----------------------------------------------------
*     
*     ( purpose )
*       Return random position in nuclei
*       (consider formation zone)
*
*       If (PIN == 0)  : XPOS is selected randomly without formation zone
*
*       If (PIN > 0)   : XPOS is modified considering formation zone.
*                        OPOS is used for the original position      
*
*     ( input ) 
*       AA     : number of nucleons (dummy)
*       PIN(3) : momentum of the particle in GeV
*       OPOS(3): initial position
*     ( output ) 
*       XPOS(3) : position (FERMI)
C
C    --- Taken from efranbll
C
C
C     
      IMPLICIT NONE

      REAL*4 AA, pin(3), opos(3)
      real*4 xpos(3)

      real*4 pos(3)
      
#include<neutparams.h>
#include<efpion.h>

C FOR FORMATION ZONES
      REAL FORMLENGTH, RNFORMLENGTH
      REAL MU2, FM2GEV
C      DATA MU2/0.08/, FM2GEV/5.07/
      DATA FM2GEV/5.07/

C--- initialization
      INTEGER*4 ICALL
      DATA ICALL/0/

      real*4     RLU
      external   RLU
      real*4     efrho
      external   efrho

      real*4     RMAX,dismax,R,HIT,DEN,pabs
      integer*4  i,idum,j

C     
      ICALL=ICALL+1
      IF(ICALL.NE.1)GO TO 200
C     
      RMAX=2.5*C
      DISMAX=0.
      DO 100 I=1,50
         R=RMAX/50.*FLOAT(I)
         DEN=EFRHO(R)
         IF(DEN.LT.DISMAX)GO TO 100
         DISMAX=DEN
 100  CONTINUE

 200  continue

      PABS=SQRT(PIN(1)**2+PIN(2)**2+PIN(3)**2)

C---  formation zone correction
      MU2 = FZMU2
C      WRITE(6,*) 'NERANBLL : FZMU2 = ',FZMU2 
      IF(MU2.le.0) THEN
         WRITE(6,*) 'WARNING NERANBLL, not applying formation zone:',
     $   FZMU2
         DO J=1,3
            XPOS(J)=POS(J)
         END DO
         RETURN
      END IF

      IF (PABS.ne.0.) then

C-- PABS != 0 -> use original position from caller function.

         POS(1)=OPOS(1)
         POS(2)=OPOS(2)
         POS(3)=OPOS(3)

C         FORMLENGTH=PABS/1000./MU2/FM2GEV  ! Bug in version<=5.1.4.2
         FORMLENGTH=PABS/MU2/FM2GEV
         RNFORMLENGTH=-1.*FORMLENGTH*LOG(RLU(IDUM))

         DO J=1,3
            IF (IFORMLEN.ge.1) THEN
C-- formation zone
               XPOS(J)=POS(J)+(PIN(J)/PABS)*RNFORMLENGTH
            ELSE
C-- Not use formation zone
               XPOS(J)=POS(J)
            ENDIF
         END DO
          
         RETURN

      else

C-- PABS == 0 -> use original position from caller function.

 300     continue
            R=RLU(IDUM)*RMAX
            HIT=RLU(IDUM)*DISMAX
            DEN=EFRHO(R)
         IF(HIT.GT.DEN)GO TO 300
C     
         CALL RNDIR(POS)
         XPOS(1)=POS(1)*R
         XPOS(2)=POS(2)*R
         XPOS(3)=POS(3)*R

         return

      endif

      write(*,*) 'neranbll:What?'

      END
