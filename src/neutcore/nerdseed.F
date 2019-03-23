      SUBROUTINE NERDSEED
*
*        READ RANDOM-SEED from the file
*
*
      IMPLICIT NONE

      real*4     ranf
      external   ranf

      INTEGER*4  ISEED
      INTEGER*4  ISEED1,ISEED2,istau1,istau2,istau3
      CHARACTER*80 RANFILE

      real*4     DUMMY

      
      call getenv('RANFILE',ranfile)
      write(*,*) "Now opening random seed file(",ranfile,")"
      open(unit=99,file=ranfile,status='unknown')
      READ(99,*) ISEED1,ISEED2,istau1,istau2,istau3

      ISEED = ISEED1

C Now ranf is calling RLU (2008/11/22, hayato)
CC--- use same seed for RANLUX : 
C      call seedranf(iseed)
C      dummy = ranf()

C--- Changed from RNDM to RANLUX
      CALL RLUXGO(3,ISEED1,0,0)

C--- Changed from RAMBO to RM48
      CALL RM48IN(ISEED2,0,0)

      CALL RMARIN(istau1,istau2,istau3)


      RETURN
      END
      
      
